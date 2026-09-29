/**
*
* @verbatim
History
-1.0  2020-06-11  Lvhb created (Network Precise Position,NPP)
@endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file		    gppprtkServer.cpp
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

#include "gnpp/gnppClient.h"
#include "gutils/gtimesync.h"

using namespace gnut;

namespace great
{
	t_gnppClient::t_gnppClient(t_gsetbase* set, t_gallproc* data, t_glog* log) :
		_set(set),
		_alldata(data),
		_log(log),
		_isBase(false),
		_isCompAug(true),
		_isWLProc(false)
	{
		_gInterpol = make_shared<t_gInterpol>(set, data, log);
		if (dynamic_cast<t_gsetgen*>(_set)->list_base().size()) _isBase = true;
   		_nppmodel = dynamic_cast<t_gsetproc*>(_set)->npp_model();
		//_obs_mode = dynamic_cast<t_gsetproc*>(_set)->obs_combin();
		//_frequency = dynamic_cast<t_gsetproc*>(_set)->frequency();
		_allaug = dynamic_cast<t_gaug*>((*_alldata)[t_gdata::AUG]);
		_allobj = dynamic_cast<t_gallobj*>((*_alldata)[t_gdata::ALLOBJ]);
		_nppmodel = dynamic_cast<t_gsetproc*>(_set)->npp_model();
		auto obs_mode = dynamic_cast<t_gsetproc*>(_set)->obs_combin();
		auto frequency = dynamic_cast<t_gsetproc*>(_set)->frequency();
		_site = dynamic_cast<t_gsetgen*>(_set)->list_rover()[0];
		t_gtime begt = dynamic_cast<t_gsetgen*>(_set)->beg();
		t_gtime endt = dynamic_cast<t_gsetgen*> (_set)->end();

		if (!_isBase && obs_mode == OBSCOMBIN::IONO_FREE && frequency == 2)_isWLProc = true;

		//Step1:set interpolation site list from gset
		auto sites = dynamic_cast<t_gsetgen*>(_set)->recs();
		_list_ipsite.clear();
		int i = 0;
		for (auto it = sites.begin(); it != sites.end(); it++) {
			if (*it == _site) continue;
			_list_ipsite.push_back(*it);
			//auto crd = dynamic_cast<t_gsetrec*>(_set)->get_crd_xyz(*it);
		}
		_nipSite = _list_ipsite.size();

		_site_vrs = "";
		if (_isBase) {
			_site_vrs = dynamic_cast<t_gsetgen*>(_set)->list_base()[0];
		}

		if (_isWLProc)
		{
			_Proc_WL = make_shared<t_gposproc>(_site, _site_vrs, _set, _alldata, _log);
			_Proc_WL->setObsMode(WL_COMBIN);
			_Proc_WL->resetLsqUpdate(_set, _alldata, _log);
			_Proc_WL->setClient(true);
		}
		_Proc = make_shared<t_gposproc>(_site, _site_vrs, _set, _alldata,_log);
		_Proc->setClient(true);

		if(!dynamic_cast<t_gsetproc*>(_set)->realtime()) {
			_prtAugFileHeader(begt, endt, _site);
		}
	}

	t_gnppClient::~t_gnppClient()
	{
		if (_augfile) { if (_augfile->is_open()) { _augfile->close(); }; delete _augfile; }
	}

	bool t_gnppClient::ProcessBatch(const t_gtime& begT, const t_gtime& endT)
	{
		gtrace("t_gnppClient::processBatch");

		if (!_alldata || begT > endT){
			write_log_info(_log, 0, "ERROR : beg > end or no proc data");
			return false;
		}

		//Step1: initiating NPP parameter
		_Proc->InitProc(begT, endT);
		if (_isWLProc)_Proc_WL->InitProc(begT, endT);

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

			////Step2 : time synchronization
			if (now != endT) {
				if (!time_sync(now, sampling, scale, _log)) { // now.sod()%(int)_rtk_set->sampling != 0 ){      
					now.add_dsec(sign * subint / 100);

					continue;
				}
				if (sampling > 1) now.reset_dsec();
			}
				
			double percent = now.diff(begT) / endT.diff(begT) * 100.0;	
			std::cerr << "\r" << now.str_ymdhms() << fixed << setprecision(1) << setw(6) << percent << "%";

			//Step3.0: judge data's avaliabily in real-time mode 
			if (dynamic_cast<t_gsetproc*>(_set)->realtime()) {

				//observation data
				while (!_Proc->avaliable(now))
					t_gtime::gmsleep(1);

				//augmentation data
				for (int j = 0; j < _nipSite; j++){
					while (now > _allaug->endTime_aug(_list_ipsite[j])) {
						t_gtime::gmsleep(1);
					}
				}
			}

			// Step 3.1 Outliers / Doppler_check 
			if (dynamic_cast<t_gsetproc*>(_set)->slip_model() == SLIPMODEL::TURBO_EDIT)
			{
				_Proc->Preprocess(now, now, sampling, NULL, NULL, true);
				if (_isWLProc)_Proc_WL->Preprocess(now, now, sampling, NULL, NULL, true);
			}
			
			//Step3.2:内插改正数
			if (!_gInterpol->interpolAug(now,_data_rover,_data_base)) {
				if (sampling > 1)  now.add_secs(int(sign * sampling));  // =<1Hz data
				else               now.add_dsec(sign * sampling);
				continue;
			}

			//Step4:preprocess
			_Proc->Preprocess(now, now, sampling, &_data_rover, &_data_base);
			if (_isWLProc)_Proc_WL->Preprocess(now, now, sampling, &_data_rover, &_data_base);
			
			//Step5: process one epoch position for rover station or baseline
			if (_isWLProc){
				int wl_valid = _Proc_WL->ProcessOneEpoch(now, &_data_rover, &_data_base);
				if (wl_valid == 1){
					auto wlparam = _Proc_WL->GetParam(now);
					_Proc->setWLParams(wlparam);
					_Proc->ProcessOneEpoch(now, &_data_rover, &_data_base);
				}
			}
			else{
				_Proc->ProcessOneEpoch(now, &_data_rover, &_data_base);
			}

			if (sampling > 1)  now.add_secs(int(sign * sampling));  // =<1Hz data
			else               now.add_dsec(sign * sampling);
		}//end time loop

		if (!dynamic_cast<t_gsetproc*>(_set)->realtime()) {
			if (!_Proc->_ExtractResfile())return false;
		}
		return true;
	}

	void t_gnppClient::_prtAugFileHeader(const t_gtime& begT, const t_gtime& endT, string site)
	{
		string aug_tmp = dynamic_cast<t_gsetout*>(_set)->outputs("aug");
		if (!aug_tmp.empty()) {
			substitute(aug_tmp, "$(rec)", site, false);
			if (_augfile) {
				if (_augfile->is_open()) _augfile->close();
			}
			_augfile = new t_giof(aug_tmp);
			_augfile->tsys(t_gtime::GPS);
			_augfile->mask(aug_tmp);
			_augfile->append(dynamic_cast<t_gsetout*>(_set)->append());
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
		t_gtriple snx_crd = _allobj->obj(site)->crd(begT);
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

		if (_augfile)
		{
			_augfile->write(os.str().c_str(), os.str().size());
			_augfile->flush();
			
			_gInterpol->printAugFile(_augfile);
		}
	}

}