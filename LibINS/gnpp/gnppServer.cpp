/**
*
* @verbatim
History
-1.0  2020-06-11 lvhb created (Network Precise Position,NPP)
@endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file		    gnppServer.cpp
* @brief		main control class for ppp/rtk estimate in network Server
*
* @author       lvhb,BoWang Wuhan University
* @version		1.0.0
				1.1.0 lvhb finished NRTK testing (Comprehensive augmentation) in 2020-10-15
				1.1.1 lvhb added and finished URTK testing (Comprehensive augmentation) in 2020-10-24
				1.1.2 lvhb tested in real-time TCP/IP solution about NRTK/URTK in 2020-11-10
* @date		    2020-06-11
*
*/

#include "gnpp/gnppServer.h"
#include "gutils/gtimesync.h"
#include "gcoders/bncobs.h"
#include "gcoders/bnccorr.h"

using namespace gnut;

namespace great
{
	t_gnppServer::t_gnppServer(t_gsetbase* set, t_gallproc* data, t_glog* log) :
		_set(set),
		_alldata(data),
		_log(log),
		_isBase(false)
	{
		if (dynamic_cast<t_gsetgen*>(_set)->list_base().size()) _isBase = true;
		_isCompAug = true;
		_nppmodel = dynamic_cast<t_gsetproc*>(_set)->npp_model();

		shared_ptr<t_gbiasmodel> precise_bias = make_shared<t_gprecisebias>(_alldata, _log, _set);
		_base_model = make_shared<t_gcombALL>(_set, _log, precise_bias, _alldata);
	}

	t_gnppServer::~t_gnppServer()
	{
		
	}

	bool t_gnppServer::_InitProcAll(const t_gtime& begT, const t_gtime& endT)
	{
		//Step1.0:set site list from gset
		set<string> sites = dynamic_cast<t_gsetgen*>(_set)->recs();
		if (sites.size() == 0) return false;

		_list_site.clear(); _list_site_base.clear(); _mapProc.clear();

		//Step1.1: get site list or base site list
		if (_isBase) {
			_list_site = dynamic_cast<t_gsetgen*>(_set)->list_rover();
			_list_site_base = dynamic_cast<t_gsetgen*>(_set)->list_base();
			if (_list_site.size() != _list_site_base.size()) return false;
		}
		else {
			for (set<string>::iterator it = sites.begin(); it != sites.end(); ++it) {
				_list_site.push_back(*it);
			}
		}
		_nSite = _list_site.size();
		if (_nSite == 0) return false;

		//Setp2:initialize posproc class
		for (unsigned int i = 0; i < _nSite; i++)
		{
			string site, site_base = "";
			site = _list_site[i];
			if (_isBase) site_base = _list_site_base[i];

			_mapProc[i] = make_shared<t_gposproc>(site, site_base, _set, _alldata, _log);
			_mapProc[i]->InitProc(begT, endT);
		}//end site

		//real-time process
		string port = dynamic_cast<t_gsetout*>(_set)->output_port("aug");
		if (port.size() > 0) {
			_augtcp = make_shared<t_gtcp>(str2int(port));
		}
		// initialize augementation file header
		if (!dynamic_cast<t_gsetproc*>(_set)->realtime())
		{
			for (set<string>::iterator it = sites.begin(); it != sites.end(); ++it) {
				_prtAugFileHeader(begT, endT, *it);
			}
		}
		return true;
	}

	bool t_gnppServer::_selMajorSite(const t_gtime& nowT, vector<int>& ind_state)
	{
		_strMajorSite = ""; _Allsducpar.clear();
		map<string, pair<int, int>> tmp_index;
		//select baseline with related major_site
		for (int i = 0; i < _nSite; i++)
		{
			if (ind_state[i] != 1) continue;// fail or float
			//pick major site for DD->UD
			_Allsducpar[i] = _mapProc[i]->GetParam(nowT);
			auto site = _list_site[i];
			auto site_base = _list_site_base[i];

			auto sdamb_nsat = _Allsducpar[i].amb_prns().size();

			if (_strMajorSite == "") {
				_strMajorSite = site_base;
				tmp_index[site_base].first = 1;
				tmp_index[site_base].second = sdamb_nsat;

				tmp_index[site].first = 1;
				tmp_index[site].second = sdamb_nsat;
			}
			else {
				string tmpsite;
				for (auto j = 0; j < 2; j++) {
					if (j) tmpsite = site;		//rover
					else tmpsite = site_base;	//base

					if (tmp_index.find(tmpsite) == tmp_index.end()) {
						tmp_index[tmpsite].first = 1;
						tmp_index[tmpsite].second = sdamb_nsat;
					}
					else {
						tmp_index[tmpsite].first++;
						if (tmp_index[tmpsite].second < sdamb_nsat) tmp_index[tmpsite].second = sdamb_nsat;

						if (tmp_index[tmpsite].first > tmp_index[_strMajorSite].first) {
							_strMajorSite = tmpsite;
						}
						else if (tmp_index[tmpsite].first == tmp_index[_strMajorSite].first)
						{
							if (tmp_index[tmpsite].second > tmp_index[_strMajorSite].second) {
								_strMajorSite = tmpsite;
							}
						}
					}
				}//end j
			}
		}

		int idxmaj = 0;
		for (int i = 0; i < _nSite; i++)
		{
			if (ind_state[i] != 1) continue;// fail or float
			else if (_list_site[i] == _strMajorSite || _list_site_base[i] == _strMajorSite)
			{
				idxmaj++;
			}
		}
		if (idxmaj < 2) return false;

		/* using same base site for setting reference ambiguity */
		if (_strMajorSite != _strMajorSite_pre && _nSite <= 3 ){
			if (_strMajorSite_pre.empty()){
				_strMajorSite_pre = _strMajorSite;
			}
			else if (tmp_index.find(_strMajorSite_pre) != tmp_index.end()) {
				_strMajorSite = _strMajorSite_pre;
			}
		}

		return true;
	}

	bool t_gnppServer::_getOneBLuducParam(const t_gtime& nowT, const int& i, bool bmajor_site)
	{
		gtrace("rtk sduc->uduc for one baseline");

		auto site = _list_site[i];
		auto site_base = _list_site_base[i];

		bool bindex = (site == _strMajorSite || site_base == _strMajorSite);
		if (!bmajor_site) bindex = (site != _strMajorSite && site_base != _strMajorSite);

		if (bindex)
		{
			auto gobs = dynamic_cast<t_gallobs*>((*_alldata)[t_gdata::ALLOBS]);
			double sign = 1.0;
			if (bmajor_site && (site == _strMajorSite))	sign = -1.0;

			auto data = gobs->obs(site, nowT);
			auto data_base = gobs->obs(site_base, nowT);
			if (data.size() == 0 || data_base.size() == 0) return false;
			auto gbias = dynamic_cast<t_gallbias*>((*_alldata)[t_gdata::ALLBIAS]);
			//get flt UC param
			t_gallpar sdparam = _Allsducpar[i];
			set<string> sdamb_prns = sdparam.amb_prns();
			if (sdamb_prns.size() < 4) return false;
			for (auto isite = 0; isite < 2; isite++)
			{
				string tmpsite = "";
				if (isite) tmpsite = site_base;	//BASE
				else tmpsite = site;			//ROVER

				if (_Alluducpar.find(tmpsite) == _Alluducpar.end()) {
					//clk
					double clk = 0.0; int ipar = 0;
					if (!_mapProc[i]->SppProcOneEpoch(nowT, tmpsite, NULL, &clk)) return false;
					t_gpar par(tmpsite, t_gpar::CLK, ++ipar, ""); par.value(clk);
					_Alluducpar[tmpsite].addParam(par);
					//trop
					if (!_isCompAug) {
						t_gpar troppar(tmpsite, t_gpar::TRP, ++ipar, "");
						int idx = sdparam.getParam(tmpsite, t_gpar::TRP, "");
						if (idx != -1) {
							troppar.value(sdparam[idx].value());
							troppar.apriori(sdparam[idx].apriori());
						}
						//else {
						//	troppar.value(0);
						//	troppar.apriori(1E-9);
						//}
						_Alluducpar[tmpsite].addParam(troppar);
					}
					//_Alluducpar[tmpsite].reIndex();
					_list_rec_bl.push_back(tmpsite);
				}
				else if (!bmajor_site) {
					sign = -1.0;
					if (isite) sign = 1.0;
				}
			}//end isite

			for (auto iter = data.begin(); iter != data.end(); ++iter) {
				string sat = iter->sat();
				if (sdamb_prns.find(sat) == sdamb_prns.end()) continue;
				t_gsatdata obsdata_other;
				for (int j = 0; j < data_base.size(); j++) {
					if (data_base[j].sat() == sat) {
						obsdata_other = (data_base)[j]; break;
					}
				}
				if (obsdata_other.sat() == "") continue;
				//dcb 
				if (gbias) {
					iter->apply_bias(*gbias);
					obsdata_other.apply_bias(*gbias);
				}

				vector < GOBSBAND > band = dynamic_cast<t_gsetgnss*>(_set)->band(iter->gsys());
				for (FREQ_SEQ f = FREQ_1; f <= band.size(); f = (FREQ_SEQ)(f + 1))
				{
					GOBSBAND b = band[f - 1];
					t_gpar::t_type amb_type;
					switch (f)
					{
						/*case gnut::FREQ_X:	amb = t_gpar::AMB_IF; break;*/
					case gnut::FREQ_1:	amb_type = t_gpar::AMB_L1; break;
					case gnut::FREQ_2:	amb_type = t_gpar::AMB_L2; break;
					case gnut::FREQ_3:	amb_type = t_gpar::AMB_L3; break;
					case gnut::FREQ_4:	amb_type = t_gpar::AMB_L4; break;
					case gnut::FREQ_5:	amb_type = t_gpar::AMB_L5; break;
					default:	break;
					}
					if (f > FREQ_5) continue;

					int idx = sdparam.getParam(site, amb_type, sat);
					if (idx < 0) continue;
					t_gsatdata dbase, drover;
					if (sign > 0) {
						drover = *iter; dbase = obsdata_other;
					}
					else {
						drover = obsdata_other; dbase = *iter;
					}
					//base amb
					if (_Alluducamb[dbase.site()].find(sat) ==
						_Alluducamb[dbase.site()].end() ||
						_Alluducamb[dbase.site()][sat].find(b) ==
						_Alluducamb[dbase.site()][sat].end()) {
						if (bmajor_site)
						{
							t_gobs gobsL = dbase.select_phase(b);
							t_gobs gobsP = dbase.select_range(b);
							if (gobsL.gobs() == X || gobsP.gobs() == X) continue;
							double obsL = dbase.obs_L(gobsL);
							double obsP = dbase.obs_C(gobsP);//dcb
							//double lam_sat = dbase.wavelength(b);
							_Alluducamb[dbase.site()][sat][b] = 1E-9;// obsL - obsP;// round((obsL - obsP) / lam_sat)* lam_sat;
							if (double_eq(_Alluducamb[dbase.site()][sat][b], 0.0)) _Alluducamb[dbase.site()][sat][b] += 1e-9;
						}
						else {
							if (_Alluducamb[drover.site()].find(sat) ==
								_Alluducamb[drover.site()].end() ||
								_Alluducamb[drover.site()][sat].find(b) ==
								_Alluducamb[drover.site()][sat].end()) {
								continue;//to use same major site
							}
							else {
								_Alluducamb[dbase.site()][sat][b] = _Alluducamb[drover.site()][sat][b]
									- sign * sdparam[idx].value();
							}
						}
					}

					//rover amb
					if (_Alluducamb[drover.site()].find(sat) ==
						_Alluducamb[drover.site()].end() ||
						_Alluducamb[drover.site()][sat].find(b) ==
						_Alluducamb[drover.site()][sat].end()) {
						_Alluducamb[drover.site()][sat][b] = _Alluducamb[dbase.site()][sat][b]
							+ sign * sdparam[idx].value();
					}

					for (auto isite = 0; isite < 2; isite++)
					{
						string tmpsite = "";
						if (isite) tmpsite = dbase.site();	//BASE
						else tmpsite = drover.site();		//ROVER
						if (_Alluducpar[tmpsite].getParam(tmpsite, amb_type, sat) < 0)
						{
							t_gpar newamb(tmpsite, amb_type, _Alluducpar[tmpsite].parNumber() + 1, sat);
							newamb.value(_Alluducamb[tmpsite][sat][b]);
							_Alluducpar[tmpsite].addParam(newamb);
						}
					}//end isite
				}//end f
			}//end data
		}
		else return false;

		return true;
	}

	bool t_gnppServer::_getAllBLuducParam(const t_gtime& nowT, vector<int>& ind_state)
	{
		gtrace("rtk sduc->uduc for all baseline");

		_Alluducpar.clear(); _Alluducamb.clear(); _list_rec_bl.clear();

		if (!_selMajorSite(nowT, ind_state)) return false;
		//_strMajorSite = "HKOH";

		//select baseline with related major_site
		for (int i = 0; i < _nSite; i++)
		{
			if (ind_state[i] != 1) continue;// fail or float
			_getOneBLuducParam(nowT, i, true);
		}

		if (_list_rec_bl.size() < 3 || _nSite > 3){
			//select baseline without related major_site
			for (int i = 0; i < _nSite; i++)
			{
				if (ind_state[i] != 1) continue;// fail or float
				_getOneBLuducParam(nowT, i, false);
			}
		}

		if (_list_rec_bl.size() < 3) return false;

		return true;
	}

	bool t_gnppServer::_getAugData(const t_gtime& nowT, vector<int>& ind_state)
	{
		_AllCorr.clear();
		if (!_isCompAug) { _AlltropAug.clear(); _AllIonAug.clear(); }

		int nSite = _nSite;
		if (_isBase) {
			if (!_getAllBLuducParam(nowT, ind_state)) return false;
			nSite = _list_rec_bl.size();
		}

		auto gobs = dynamic_cast<t_gallobs*>((*_alldata)[t_gdata::ALLOBS]);
		for (int i = 0; i < nSite; i++)
		{
			string site = "";
			if (_isBase) {
				site = _list_rec_bl[i];
			}
			else {
				if (ind_state[i] != 1) continue;// fail or float
				site = _list_site[i];
			}
			if (site == "") continue;
			auto data = gobs->obs(site, nowT);
			if (data.size() == 0) continue;

			//get flt/lsq UC param
			t_gallpar param;
			if (_isBase)	param = _Alluducpar[site];
			else param = _mapProc[i]->GetParam(nowT);

			set<string> amb_prns = param.amb_prns();
			if (amb_prns.size() < 4)	continue;

			// set TRP and ION value 0 and it's precise
			for (unsigned int ipar = 0; ipar < param.parNumber(); ipar++)
			{
				if (param[ipar].str_type().find("TRP") != string::npos ||
					param[ipar].str_type().find("ION") != string::npos)
				{
					param[ipar].value(0);
					param[ipar].apriori(1E-9);
				}
			}//end param

			ostringstream os;
			os << "> " << left << setw(35) << nowT.str("%Y %m %d %H %M %p")
				<< setw(5) << amb_prns.size() << endl;

			t_gtime obsEpo = data.begin()->epoch();

			// use combmodel
			for (auto iter = data.begin(); iter != data.end(); ++iter) {
				string sat = iter->sat();
				if (amb_prns.find(sat) == amb_prns.end())	continue;
				// save omc L,P
				auto band_index = dynamic_cast<t_gcombmodel*>(&(*_base_model))->get_band_index();
				map<FREQ_SEQ, GOBSBAND> crt_bands = band_index[iter->gsys()];
				int nf = crt_bands.size();
				if (nf >= 2)	iter->tb12(true);
				if (nf >= 3)	iter->tb13(true);
				if (nf >= 4)	iter->tb14(true);
				if (nf >= 5)	iter->tb15(true);
				//cout << "CompAu********************" << endl;
				t_glsqEquationMatrix crt_equ;
				bool irc = _base_model->cmb_equ(obsEpo, param, *iter, crt_equ);
				if (irc == false) continue;

				ostringstream rtos;
				rtos << "> AUG " << left  << setw(10) << obsEpo.gwk() << setw(15) << fixed << setprecision(5) << (obsEpo.sow() + obsEpo.dsec())
					<< setw(10) << site;

				os << right << setw(3) << iter->sat();
				rtos << right << setw(3) << iter->sat();

				if (!_isCompAug /*&& _isBase*/) {
					if (_isBase)	param = _Alluducpar[site];
					if (!_getUDIonAug(i, *iter, crt_bands, param)) continue;
					if (!_getUDTropAug(nowT, i, *iter, param)) continue;
				}
				for (const auto& itband : crt_bands)
				{
					GOBSBAND band = itband.second;
					t_gobs  obsP = iter->select_range(band);
					t_gobs  obsL = iter->select_phase(band);

					int indexP = crt_equ.find_equ(site, sat, t_gobscombtype(obsP, OBSCOMBIN::RAW_ALL));
					int indexL = crt_equ.find_equ(site, sat, t_gobscombtype(obsL, OBSCOMBIN::RAW_ALL));
					double atmosP = 0.0, atmosL = 0.0;
					if (indexP >= 0 && indexL >= 0) {
						if (!_isCompAug /*& _isBase*/) {
							atmosP = _AllIonAug[i][sat][band] + _AlltropAug[i][sat];
							atmosL = -_AllIonAug[i][sat][band] + _AlltropAug[i][sat];
						}
						_AllCorr[i][sat][band].first = crt_equ.l[indexP] - atmosP;//range
						_AllCorr[i][sat][band].second = crt_equ.l[indexL] - atmosL;//phase

						os << right
							<< setw(12) << fixed << setprecision(4) << _AllCorr[i][sat][band].first
							<< setw(12) << fixed << setprecision(4) << _AllCorr[i][sat][band].second;
						rtos << right
							<< setw(12) << fixed << setprecision(4) << _AllCorr[i][sat][band].first
							<< setw(12) << fixed << setprecision(4) << _AllCorr[i][sat][band].second;
					}
					else {
						os << right << setw(12) << "" << setw(12) << "";
						rtos << right << setw(12) << "" << setw(12) << "";
					}
				}//end band loop
				if (!_isCompAug/*&&_isBase*/)
				{
					//output ion
					for (const auto& itband : crt_bands) {
						GOBSBAND band = itband.second;
						if (_AllIonAug[i][sat].find(band) != _AllIonAug[i][sat].end()) {
							os << right
								<< setw(12) << fixed << setprecision(4) << _AllIonAug[i][sat][band];
						}
						else {
							os << right << setw(12) << "" << setw(12) << "";
						}
					}//end band
					//output trop
					if (_AlltropAug[i].find(sat) != _AlltropAug[i].end()) {
						os << right << setw(12) << fixed << setprecision(4) << _AlltropAug[i][sat];
					}
					else {
						os << right << setw(12) << "" << setw(12) << "";
					}
					
				}
				
				if (0)//output ambiguity,tested by lvhb 202007
				{
					auto freq_index = dynamic_cast<t_gcombmodel*>(&(*_base_model))->get_freq_index();
					auto crt_freqs = freq_index[iter->gsys()];
					for (const auto& itband : crt_bands)
					{
						GOBSBAND band = itband.second;
						t_gobsgnss gnss(sat);
						if (_isBase)
						{
							os << right << setw(12) << fixed << setprecision(4) << (_Alluducamb[site][sat][band] / gnss.wavelength(band));
						}
						else{
							FREQ_SEQ f = crt_freqs[band];
							t_gpar::t_type amb_type;
							switch (f)
							{
								/*case gnut::FREQ_X:	amb = t_gpar::AMB_IF; break;*/
								case gnut::FREQ_1:	amb_type = t_gpar::AMB_L1; break;
								case gnut::FREQ_2:	amb_type = t_gpar::AMB_L2; break;
								case gnut::FREQ_3:	amb_type = t_gpar::AMB_L3; break;
								case gnut::FREQ_4:	amb_type = t_gpar::AMB_L4; break;
								case gnut::FREQ_5:	amb_type = t_gpar::AMB_L5; break;
								default:	break;
							}
							int ipar = param.getParam(site, amb_type, sat);
							if(ipar<0) continue;
							os << right << setw(12) << fixed << setprecision(4) << (param[ipar].value() / gnss.wavelength(band));
						}

						
					}
				}
				os << endl;
				rtos << endl;
				if (_augtcp) {
					_augtcp->run_send((rtos.str()));
				}
			}//end data loop

			// Print aug data
			
			if (_augfile[site])
			{
				_augfile[site]->write(os.str().c_str(), os.str().size());
				_augfile[site]->flush();
			}
		}//end _nSite in one epoch
		return true;
	}

	bool t_gnppServer::_getUDIonAug(const int& i, t_gsatdata& data,map<FREQ_SEQ, GOBSBAND>& crt_bands, t_gallpar& param)
	{
		GOBSBAND b1 = crt_bands[FREQ_1];
		GOBSBAND b2 = crt_bands[FREQ_2];

		GOBS gobs_L1 = data.select_phase(b1);
		double obs_L1 = data.obs_L(gobs_L1);
		if (gobs_L1 == X || double_eq(obs_L1, 0.0)) return false;

		GOBS gobs_L2 = data.select_phase(b2);
		double obs_L2 = data.obs_L(gobs_L2);
		if (gobs_L2 == X || double_eq(obs_L2, 0.0)) return false;

		double amb_L1, amb_L2;
		int ipar = param.getParam(data.site(), t_gpar::AMB_L1, data.sat());
		if (ipar != -1) {
			double lam = 1.0;// data.wavelength(b1);
			amb_L1 = lam * param[ipar].value();
		}
		else return false;

		ipar = param.getParam(data.site(), t_gpar::AMB_L2, data.sat());
		if (ipar != -1) {
			double lam = 1.0;// data.wavelength(b2);
			amb_L2 = lam * param[ipar].value();
		}
		else return false;

		auto& ion = _AllIonAug[i][data.sat()];
		ion[b1] = (obs_L1 - obs_L2) - (amb_L1 - amb_L2);

		double f1 = data.frequency(b1);
		double f2 = data.frequency(b2);
		double coef = SQR(f2) / (SQR(f1) - SQR(f2));
		ion[b1] *= coef;

		for (const auto& itband : crt_bands)
		{
			GOBSBAND band = itband.second;
			if (band != b1) {
				f2 = data.frequency(band);
				coef = SQR(f1) / SQR(f2);
				ion[band] = coef * ion[b1];
			}
		}

		return true;
	}

	bool t_gnppServer::_getUDTropAug(const t_gtime& now, const int& i, t_gsatdata& data, t_gallpar& param)
	{
		////saas
		//auto trpModStr = dynamic_cast<t_gsetproc*>(_set)->tropo_model();
		//shared_ptr<t_gtropo> tropoModel;
		//if (trpModStr == SAASTAMOINEN)		 tropoModel = make_shared<t_saast>();
		//else if (trpModStr == DAVIS)         tropoModel = make_shared<t_davis>();
		//else if (trpModStr == HOPFIELD)      tropoModel = make_shared<t_hopf>();
		////else if (trpModStr == MOPS)          tropoModel = make_shared<t_blindmops>();
		//else return false;

		string site = data.site();
		auto grec = dynamic_cast<t_gallobj*>((*_alldata)[t_gdata::ALLOBJ])->obj(site);
		// Cartesian coordinates to ellipsodial coordinates
		t_gtriple xyz, ell;
		if (_isBase) {
			if (param.getCrdParam(site, xyz) < 0){
				xyz = grec->crd_arp(now);
			}
		}
		else {
			if (param.getCrdParam(site, xyz) < 0){
				xyz = grec->crd(now);
			}
			xyz += grec->eccxyz(now);
		}

		xyz2ell(xyz, ell, false);

		string sat = data.sat();
		double ele = data.ele();
		double zwd = 0.0, zhd = 0.0;
		int ipar = param.getParam(site, t_gpar::TRP, "");
		if (ipar >= 0)
		{
			zwd = param[ipar].value();
			zhd = param[ipar].zhd;
		}
		else return false;

		auto tropo_mf = dynamic_cast<t_gsetproc*>(_set)->tropo_mf();
		auto& delay = _AlltropAug[i][sat];

		double mfh, mfw, dmfh, dmfw;
		mfh = mfw = dmfh = dmfw = 0.0;
		if (tropo_mf == GMF)
		{
			t_gmf mf;
			mf.gmf(now.mjd(), ell[0], ell[1], ell[2], G_PI / 2.0 - ele,
				mfh, mfw, dmfh, dmfw);
		}
		else if (tropo_mf == COSZ)
		{
			mfh = mfw = 1 / sin(ele);
			dmfh = dmfw = -(cos(ele)) / (sin(ele) * sin(ele));
		}
		delay = mfh * zhd + mfw * zwd;

		int j = param.getParam(site, t_gpar::GRD_N, "");
		int k = param.getParam(site, t_gpar::GRD_E, "");

		if (j >= 0 && k >= 0) {
			double grdN, grdE;
			grdN = grdE = 0.0;
			double azi = data.azi();
			auto grad_mf = dynamic_cast<t_gsetproc*>(_set)->grad_mf();

			if (grad_mf == TILTING) {
				grdN = param[j].value() * dmfw * cos(azi);
				grdE = param[k].value() * dmfw * sin(azi);
			}
			else if (grad_mf == CHEN_HERRING) {
				double mfg = 1.0 / (sin(ele) * tan(ele) + 0.0032);
				grdN = param[j].value() * 1000.0 * mfg * cos(azi);  grdN /= 1000.0;
				grdE = param[k].value() * 1000.0 * mfg * sin(azi);  grdE /= 1000.0;
			}
			else if (grad_mf == BAR_SEVER) {
				double mfg = mfw * (1 / tan(ele));
				grdN = param[j].value() * mfg * cos(azi);
				grdE = param[k].value() * mfg * sin(azi);
			}

			delay += grdN + grdE;
		}

		return true;
	}

	bool t_gnppServer::ProcessBatch(const t_gtime& begT, const t_gtime& endT)
	{
		gtrace("t_gnppServer::processBatch");

		if (!_alldata || begT > endT) {
			write_log_info(_log, 0, "ERROR : beg > end or no proc data");
			return false;
		}

		//Step1: initiating NPP parameter
		_InitProcAll(begT, endT);

		//Step2: pre process post-time(delete by lvhb)

		//Step3: time synchronization
		double sampling = dynamic_cast<t_gsetgen*>(_set)->sampling();
		if (sampling == 0) sampling = dynamic_cast<t_gsetgen*>(_set)->sampling_default();
		int scale = dynamic_cast<t_gsetgen*>(_set)->sampling_scalefc();
		int sign = 1;
		double subint = 0.1;
		if (scale > 0) subint = 1.0 / scale;
		if (sampling > 1) subint = pow(10, floor(log10(sampling)));
		t_gtime now(begT);
		bool time_loop = true;

		while (time_loop)
		{
			// Check whether current time is after begin time and before end time
			if (now < endT || now == endT) time_loop = true;
			else if (now > endT) { time_loop = false; break; }

			//Step3 : time synchronization
			if (now != endT) {
				if (!time_sync(now, sampling, scale, _log)) { // now.sod()%(int)_rtk_set->sampling != 0 ){      
					now.add_dsec(sign * subint / 100);

					continue;
				}
				if (sampling > 1) now.reset_dsec();
			}

			double percent = now.diff(begT) / endT.diff(begT) * 100.0;
			std::cerr << "\r" << now.str_ymdhms() << fixed << setprecision(1) << setw(6) << percent << "%";

			//Step4: process one epoch positon for all station or baseline
			vector<int> ind_state(_nSite);
			int sta_fix = 0;
//#pragma omp parallel for
			for (int i = 0; i < _nSite; i++)
			{

				//3.0 judge data's avaliabily in real-time mode
				if (dynamic_cast<t_gsetproc*>(_set)->realtime()){
					while (!_mapProc[i]->avaliable(now))
						t_gtime::gmsleep(1);
				}
				
				//3.1 preprocess for real time
				_mapProc[i]->Preprocess(now, now, sampling,NULL,NULL);// test result= postpreprocess by lvhb

				//3.2 process for one epoch
				ind_state[i] = _mapProc[i]->ProcessOneEpoch(now, NULL, NULL);
				if (ind_state[i] == 1) {
					sta_fix++;//fix number
				}
			}

			double nfix = 2;
			//cerr << ind_state[0] << " " << ind_state[1] << " " << ind_state[2] << endl;
			if (sta_fix < nfix) {
				if (sampling > 1)  now.add_secs(int(sign * sampling));  // =<1Hz data
				else               now.add_dsec(sign * sampling);
				continue;
			}

			//Step4 calculated or output Augmentation information
			_getAugData(now, ind_state);

			if (sampling > 1)  now.add_secs(int(sign * sampling));  // =<1Hz data
			else               now.add_dsec(sign * sampling);
		}//end time loop

		if (!dynamic_cast<t_gsetproc*>(_set)->realtime()) {
			for (int i = 0; i < _nSite; i++)
			{
				_mapProc[i]->_ExtractResfile();
			}
		}

		return true;
	}

	void t_gnppServer::_prtAugFileHeader(const t_gtime& begT, const t_gtime& endT, string site)
	{
		string aug_tmp = dynamic_cast<t_gsetout*>(_set)->outputs("aug");
		if (!aug_tmp.empty()) {
			substitute(aug_tmp, "$(rec)", site, false);
			if (_augfile[site]) {
				if (_augfile[site]->is_open()) _augfile[site]->close();
			}
			_augfile[site] = make_shared<t_giof>(aug_tmp);
			_augfile[site]->tsys(t_gtime::GPS);
			_augfile[site]->mask(aug_tmp);
			_augfile[site]->append(dynamic_cast<t_gsetout*>(_set)->append());
		}
		else return;

		// add header of aug file
		ostringstream os;
		double sampling = dynamic_cast<t_gsetgen*>(_set)->sampling();
		set<string> _sys_list = dynamic_cast<t_gsetgen*>(_set)->sys();

		map< GSYS, map<FREQ_SEQ, GOBSBAND> > band_index;
		band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(_set)->band_index(gnut::GPS);
		band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(_set)->band_index(gnut::GAL);
		band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(_set)->band_index(gnut::GLO);
		band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(_set)->band_index(gnut::BDS);
		band_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(_set)->band_index(gnut::QZS);

		// receiver crd
		auto gobj = dynamic_cast<t_gallobj*>((*_alldata)[t_gdata::ALLOBJ]);
		t_gtriple snx_crd = gobj->obj(site)->crd(begT);
		t_gtriple xml_crd = dynamic_cast<t_gsetrec*>(_set)->get_crd_xyz(site);
		t_gtriple final_crd(0, 0, 0);
		if (!snx_crd.zero())
		{
			final_crd = snx_crd;
		}
		else if (!xml_crd.zero() && snx_crd.zero())
		{
			final_crd = xml_crd;
		}

		os << left << fixed
			<< setw(60) << "Great v1.0   GREAT"
			<< "PGM / RUN BY"
			<< endl
			<< setw(60) << site
			<< "MARKER NAME"
			<< endl
			<< setw(20) << setprecision(4) << final_crd[0]
			<< setw(20) << setprecision(4) << final_crd[1]
			<< setw(20) << setprecision(4) << final_crd[2]
			<< "APPROX POSITION XYZ"
			<< endl
			<< right << setprecision(3) << setw(10) << sampling << setw(50) << ""
			<< "INTERVAL"
			<< endl
			<< left << setw(60) << begT.str("%Y %m %d %H %M %S")
			<< "TIME OF FIRST OBS"
			<< endl
			<< setw(60) << endT.str("%Y %m %d %H %M %S")
			<< "TIME OF LAST OBS"
			<< endl;


		for (auto iter = _sys_list.begin(); iter != _sys_list.end(); iter++)
		{
			GSYS sys = t_gsys::str2gsys(*iter);
			map<FREQ_SEQ, GOBSBAND> crt_bands = band_index[sys];
			int band_num = crt_bands.size();
			int outnum = (/*_isBase &&*/ !_isCompAug) ? band_num * 3 + 1 : band_num * 2;

			os << left
				<< setw(3) << t_gsys::gsys2char(sys) << setw(4) << outnum;
			for (const auto& itband : crt_bands)
			{
				os << setw(3) << "P" + to_string(itband.second)
					<< setw(3) << "L" + to_string(itband.second);
			}
			if (/*_isBase&&*/!_isCompAug)
			{
				for (const auto& itband : crt_bands)
				{
					os << setw(5) << "ION" + to_string(itband.second);
				}
				os << setw(5) << "TRP1";
			}

			int blank = (/*_isBase &&*/ !_isCompAug) ? 60 - 11 * band_num - 2 : 60 - 6 * band_num - 7;
			os << setw(blank) << "" << "SYS / # / AUG TYPES" << endl;
		}
		os << setw(60) << "" << "END OF HEADER" << endl;

		if (_augfile[site])
		{
			_augfile[site]->write(os.str().c_str(), os.str().size());
			_augfile[site]->flush();
		}
	}

}