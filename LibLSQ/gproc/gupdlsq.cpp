/**
*
* @verbatim
History
-1.0 BoWong xjhan glfeng 2019-02-22  creat the file.
-1.1 glfeng  2019-03-02  Adjusted some functions and variables.
-1.2 BoWong  2019-04-09  Adding Doxygen Style Code Remarks.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gupdlsq.cpp
* @brief	compute Wide-Lane, Ex-Wide-Lane and Narrow-Lane upd.
*
* @author   glfeng BoWong xjhan , Wuhan University
* @version	1.0.0
* @date		2019-04-09
*
*/

#include <algorithm>
#include "gproc/gupdlsq.h"
#include "gio/gfile.h"
#include "gcoders/rinexo.h"

#define DEBUG 1

namespace great
{
	t_gupdlsq::t_gupdlsq(UPDTYPE updmode, bool carrier_range, t_gsetbase* set, t_gallproc* data, t_glog* log) :
		_obs_type(OBSCOMBIN::IONO_FREE),
		_bds_codebias_cor(set),
		_outliers_proc(set, log)
	{
		// set settings && log file
		_set = set;
		_glog = log;
		_data = data;

		// set upd estimate mode
		_updmode = updmode;

		// whether write carrier_range
		_carrier_range = carrier_range;

		// IF/UDUC only the two type
		if (dynamic_cast<t_gsetproc*>(_set)->obs_combin() == OBSCOMBIN::RAW_ALL)
		{
			_obs_type = OBSCOMBIN::RAW_ALL;
		}

		_band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::GPS);
		_band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::GAL);
		_band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::GLO);
		_band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::BDS);
		_band_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::QZS);
		
		if (_cal_UPD)
		{
			delete _cal_UPD;
			_cal_UPD = nullptr;
		}
		_cal_UPD = new t_gupd();
		_cal_UPD->set_est_updtype(updmode);

		// set Out
		_setOut();
	}


	t_gupdlsq::~t_gupdlsq()
	{
		if (_cal_UPD) { delete _cal_UPD; _cal_UPD = nullptr; }
		if (_fambres) { if (_fambres->is_open()) { _fambres->close(); }; delete _fambres; _fambres = nullptr; }
	}

	bool t_gupdlsq::GenerateProduct()
	{
		// Encoder IFCB/UPD files
		string path = dynamic_cast<t_gsetout*>(_set)->outputs("upd");

		t_gcoder* gencode = new t_upd(_set, "", 4096);
		gencode->clear();
		gencode->glog(_glog);
		gencode->path(path);
		gencode->add_data("IDX", _cal_UPD);

		t_gfile* gout = new t_gfile;
		gout->glog(_glog);
		gout->path(path);
		gout->coder(gencode);
		write_log_info(_glog, 1, " t_gupdlsq::GenerateProduct :  Writing IFCB/UPD file " + path);
		gout->run_write();
		delete gout; gout = NULL;
		delete gencode; gencode = NULL;

		return true;
	}

	bool t_gupdlsq::_initUPDProcData(t_gallproc * data)
	{
		//gtrace("t_gupdlsq::_initUPDProcData");
		set<string> sys = dynamic_cast<t_gsetgen*>(_set)->sys();
		if (sys.size() > 1) {
			write_log_info(_glog, 1, " ERROR[t_gupdlsq::_initUPDProcData] : Only One SYS can be process every time!");
			return false;
		}
		_sys = *sys.begin();

		// prepare data
		if (_updmode == UPDTYPE::NL)
		{
			_allambupd = dynamic_cast<t_gallambupd*>((*data)[t_gdata::AMBUPD]);
			_allWLUPD = dynamic_cast<t_gupd*>((*data)[t_gdata::UPD]);
			if (_carrier_range) {
				_allambflag = dynamic_cast<t_gallambflag*>((*data)[t_gdata::AMBFLAG]);
				_allambflag13 = dynamic_cast<t_gallambflag*>((*data)[t_gdata::AMBFLAG13]);
			}
			if (!_allambupd || !_allWLUPD)
			{
				write_log_info(_glog, 1, " ERROR[t_gupdlsq::_initUPDProcData] : can not init NL process data.");
				return false;
			}
			_site_ambupd_end = _allambupd->getSitesEndRecord();
		}
		else
		{
			_gifcb = dynamic_cast<t_gifcb*>((*data)[t_gdata::IFCB]);
			_allobs = dynamic_cast<t_gallobs*>((*data)[t_gdata::ALLOBS]);
			_allbias = dynamic_cast<t_gallbias*>((*data)[t_gdata::ALLBIAS]);
			_allambupd = dynamic_cast<t_gallambupd*>((*data)[t_gdata::AMBUPD]);
			_allambflag = dynamic_cast<t_gallambflag*>((*data)[t_gdata::AMBFLAG]);
			_allambflag13 = dynamic_cast<t_gallambflag*>((*data)[t_gdata::AMBFLAG13]);
			_allambflag14 = dynamic_cast<t_gallambflag*>((*data)[t_gdata::AMBFLAG14]);
			_allambflag15 = dynamic_cast<t_gallambflag*>((*data)[t_gdata::AMBFLAG15]);
			if (_obs_type == OBSCOMBIN::IONO_FREE && (!_allobs || !_allbias || !_allambflag)
				|| _obs_type == OBSCOMBIN::RAW_ALL && !_allambupd)
			{
				write_log_info(_glog, 1, " ERROR[t_gupdlsq::_initUPDProcData] : can not init EWL/WL/IFCB process data.");
				return false;
			}
			if (_obs_type == OBSCOMBIN::IONO_FREE && (_updmode == UPDTYPE::EWL || _updmode == UPDTYPE::EWL_epoch) && sys.find("GPS") != sys.end() && !_gifcb
				|| _obs_type == OBSCOMBIN::IONO_FREE && (_updmode == UPDTYPE::EWL || _updmode == UPDTYPE::EWL_epoch) && !_allambflag13)
			{

				write_log_info(_glog, 1, " ERROR[t_gupdlsq::_initUPDProcData] : can not init EWL/IFCB process data.");
				return false;
			}
			if ((_obs_type == OBSCOMBIN::IONO_FREE && _updmode == UPDTYPE::EWL24 && !_allambflag14)
				|| (_obs_type == OBSCOMBIN::IONO_FREE && _updmode == UPDTYPE::EWL25 && !_allambflag15))
			{

				write_log_info(_glog, 1, " ERROR[t_gupdlsq::_initUPDProcData] : can not init EWL24/EWL25 process data.");
				return false;
			}
		}

		// GLO freq ID / BDS codebias
		_allnav = dynamic_cast<t_gallnav*>((*data)[t_gdata::GRP_EPHEM]);
		if (_sys == "GLO")
		{
			_glo_freq_num = _allnav->glo_freq_num();
		}

		write_log_info(_glog, 1, " INIT [t_gupdlsq::_initUPDProcData] : process data in UPD Estimation is prepared.");
		
		return true;
	}

	bool t_gupdlsq::_initSatSite()
	{
		try
		{
			set<string> sat_rm = dynamic_cast<t_gsetgen*>(_set)->sat_rm();
			if (_updmode == UPDTYPE::NL || _obs_type == OBSCOMBIN::RAW_ALL)
			{
				// set parlist<first sats, then sites> + sites + sats according ambupd
				// get sat list
				set<string> sats1, sats2, sats3;
				sats1 = _allambupd->get_sats();
				sats2 = GNSS_SATS()[t_gsys::str2gsys(_sys)];

				// 交集 sats1+sats2并集存于sats3
				set_intersection(sats1.begin(), sats1.end(), sats2.begin(), sats2.end(), inserter(sats3, sats3.begin()));
				// 差集(sats3中属于sat_rm的元素去掉)
				set_difference(sats3.begin(), sats3.end(), sat_rm.begin(), sat_rm.end(), inserter(_sats, _sats.begin()));
				// push sats into _parlist				
				for (auto it = _sats.begin(); it != _sats.end(); it++)
				{
					if (_updmode == UPDTYPE::WL) _cal_UPD->re_init_upd(UPDTYPE::WL, WL_IDENTIFY, *it);   
					if (_updmode == UPDTYPE::EWL) _cal_UPD->re_init_upd(UPDTYPE::EWL, EWL_IDENTIFY, *it); 
					if (_updmode == UPDTYPE::EWL24) _cal_UPD->re_init_upd(UPDTYPE::EWL24, EWL24_IDENTIFY, *it);
					if (_updmode == UPDTYPE::EWL25) _cal_UPD->re_init_upd(UPDTYPE::EWL25, EWL25_IDENTIFY, *it);
					
					_parlist.push_back(*it);
				}
				// get site list && init site wl upd && push sites into _parlist
				_sites = _allambupd->get_sites();
				for (set<string>::iterator iter_site = _sites.begin(); iter_site != _sites.end(); iter_site++)
				{
					
					if (_updmode == UPDTYPE::WL )  _cal_UPD->re_init_upd(UPDTYPE::WL, WL_IDENTIFY, *iter_site);  
					if (_updmode == UPDTYPE::EWL ) _cal_UPD->re_init_upd(UPDTYPE::EWL, EWL_IDENTIFY, *iter_site); 
					if (_updmode == UPDTYPE::EWL24 ) _cal_UPD->re_init_upd(UPDTYPE::EWL24, EWL24_IDENTIFY, *iter_site);
					if (_updmode == UPDTYPE::EWL25 ) _cal_UPD->re_init_upd(UPDTYPE::EWL25, EWL25_IDENTIFY, *iter_site);
					if (_updmode == UPDTYPE::NL) {
						_allWLUPD->re_init_upd(UPDTYPE::WL, WL_IDENTIFY, *iter_site);
						_allWLUPD->re_init_upd(UPDTYPE::EWL, EWL_IDENTIFY, *iter_site);
					}
					_parlist.push_back(*iter_site);
				}
			}
			else
			{
				// set parlist<first sats, then sites> + sites + sats according to log/log13
				set<string> sats1, sats2, sats3, tmp1, tmp2;
				tmp1 = _allambflag->getSatList();
				if (_updmode == UPDTYPE::EWL || _updmode == UPDTYPE::EWL_epoch)
				{
					tmp2 = _allambflag13->getSatList();
					set_intersection(tmp1.begin(), tmp1.end(), tmp2.begin(), tmp2.end(), inserter(sats1, sats1.begin()));
				}
				else if (_updmode == UPDTYPE::EWL24)
				{
					tmp2 = _allambflag14->getSatList();
					set_intersection(tmp1.begin(), tmp1.end(), tmp2.begin(), tmp2.end(), inserter(sats1, sats1.begin()));
				}
				else if (_updmode == UPDTYPE::EWL25)
				{
					tmp2 = _allambflag15->getSatList();
					set_intersection(tmp1.begin(), tmp1.end(), tmp2.begin(), tmp2.end(), inserter(sats1, sats1.begin()));
				}
				else
				{
					sats1 = tmp1;
				}
				sats2 = GNSS_SATS()[t_gsys::str2gsys(_sys)];
				// 交集 sats1+sats2并集存于sats3
				set_intersection(sats1.begin(), sats1.end(), sats2.begin(), sats2.end(), inserter(sats3, sats3.begin()));
				// 差集(sats3中属于sat_rm的元素去掉)
				set_difference(sats3.begin(), sats3.end(), sat_rm.begin(), sat_rm.end(), inserter(_sats, _sats.begin()));
				// push sats into _parlist && init sat wl upd				
				for (auto it = _sats.begin(); it != _sats.end(); it++)
				{
					if (_updmode == UPDTYPE::WL)  _cal_UPD->re_init_upd(UPDTYPE::WL, WL_IDENTIFY, *it);  
					if (_updmode == UPDTYPE::EWL) _cal_UPD->re_init_upd(UPDTYPE::EWL, EWL_IDENTIFY, *it); 
					if (_updmode == UPDTYPE::EWL24) _cal_UPD->re_init_upd(UPDTYPE::EWL24, EWL24_IDENTIFY, *it);
					if (_updmode == UPDTYPE::EWL25) _cal_UPD->re_init_upd(UPDTYPE::EWL25, EWL25_IDENTIFY, *it);
					_parlist.push_back(*it);
				}

				// get site list && init site wl upd && push sites into _parlist
				set<string> sites1, sites2, sites3, sites4;
				sites1 = dynamic_cast<t_gsetgen*>(_set)->recs();
				sites2 = _allambflag->getSiteList();
				sites4 = _allobs->getSiteList();
				set_intersection(sites1.begin(), sites1.end(), sites2.begin(), sites2.end(), inserter(sites3, sites3.begin()));
				set_intersection(sites3.begin(), sites3.end(), sites4.begin(), sites4.end(), inserter(_sites, _sites.begin()));

				//_sites = dynamic_cast<t_gsetgen*>(_set)->recs();
				for (set<string>::iterator iter_site = _sites.begin(); iter_site != _sites.end(); iter_site++)
				{
					if (_updmode == UPDTYPE::WL)  _cal_UPD->re_init_upd(UPDTYPE::WL,  WL_IDENTIFY, *iter_site); 
					if (_updmode == UPDTYPE::EWL) _cal_UPD->re_init_upd(UPDTYPE::EWL, EWL_IDENTIFY, *iter_site); 
					if (_updmode == UPDTYPE::EWL24) _cal_UPD->re_init_upd(UPDTYPE::EWL24, EWL24_IDENTIFY, *iter_site);
					if (_updmode == UPDTYPE::EWL25) _cal_UPD->re_init_upd(UPDTYPE::EWL25, EWL25_IDENTIFY, *iter_site);
	
					_parlist.push_back(*iter_site);
				}
			}
			return true;
		}

		catch (...)
		{
			write_log_info(_glog, 1, " ERROR[t_gupdlsq::_initSatSite] : Unknown Mistake");
			return false;
		}
	}

	void t_gupdlsq::_setOut()
	{
		gtrace("t_gupdlsq::_setOut");
		// set ambres file
		_fambres = new t_giof;
		_fambres->tsys(t_gtime::GPS);
		set<string> sys = dynamic_cast<t_gsetgen*>(_set)->sys();
		GOBSBAND b1 = _band_index[t_gsys::str2gsys(*sys.begin())][FREQ_1];
		GOBSBAND b2 = _band_index[t_gsys::str2gsys(*sys.begin())][FREQ_2];
		printf("t_gupdlsq::_setOut : %s_%s%s-res\n", updmode2str(_updmode).c_str(), sys.begin()->c_str(), gobsband2str(b1).c_str(), gobsband2str(b2).c_str());
		_fambres->mask(updmode2str(_updmode) + gobsband2str(b1) + gobsband2str(b2) 
			+ "_" + *sys.begin() + "-res");
 		_fambres->append(dynamic_cast<t_gsetout*>(_set)->append());
	}

	bool t_gupdlsq::processBatch(const t_gtime& beg, const t_gtime& end, const double &intv)
	{
		/** Prepare data */
		if (!_initUPDProcData(_data)) return false;

		/** Init sat/site info */
		if (!_initSatSite()) return false;

		/** Init wavelen */
		if (!_setWavelength()) return false;

		if (_updmode == UPDTYPE::WL || _updmode == UPDTYPE::EWL || _updmode == UPDTYPE::EWL_epoch
			|| _updmode == UPDTYPE::EWL24 || _updmode == UPDTYPE::EWL25)
		{
			if (!_processWL(beg, end, intv)) return false;
		}
		else if (_updmode == UPDTYPE::NL)
		{
			if (!_processNL(beg, end, intv)) return false;
		}
		else if (_updmode == UPDTYPE::IFCB)
		{
			if (!_processIFCB(beg, end, intv)) return false;
		}
		return true;
	}

	bool t_gupdlsq::_processWL(const t_gtime& beg, const t_gtime& end, const double &intv)
	{
		set<string>::iterator it_site = _sites.begin();
		map <string, t_epo_sat_amb> all_site_amb;  // ewl_epoch
		READ_OFILE_MODE rom = dynamic_cast<t_gsetproc*>(_set)->read_ofile_mode();
		/** Loop sites */
		for (it_site; it_site != _sites.end(); it_site++)
		{
			t_gallobs  obs;  obs.glog(_glog); obs.gset(_set);
			if (rom == READ_OFILE_MODE::REALTIME) {
				_readOfile(it_site,obs);//read current site's obs
				_allobs = &obs;
			}
			

#ifdef DEBUG
			std::cout << "" << endl;
			std::cout << "----- Processing " + *it_site << endl;
			std::cout << "" << endl;
#endif
			/** Init sat_site pair ambiguities */
			auto it_sat = _sats.begin();
			for (; it_sat != _sats.end(); it_sat++)
			{
				_sat_site_amb[*it_site][*it_sat] = 9.90;  /// == xx(ambiguity of satellite-station pair)
			}

			t_all_sats_OW OW;             /// == OW
			t_epo_sat_amb one_site_amb;   /// == wl/ewl

			/** Form WL/EWL observations */
			if (_obs_type == OBSCOMBIN::IONO_FREE)
			{
				if (!t_gupdlsq::_combineObsIF(*it_site, beg, end, intv, one_site_amb, OW)) return false;
			}
			else
			{
				if (!t_gupdlsq::_combineObsALL(*it_site, beg, end, intv, one_site_amb, OW)) return false;
			}
#ifdef DEBUG
			auto itsat = OW.begin();
			for (itsat; itsat != OW.end(); itsat++)
			{
				auto itepo = itsat->second.begin();
				for (itepo; itepo != itsat->second.end(); itepo++) {
					std::cout << updmode2str(_updmode) << " " << itsat->first << setw(16) <<
						(*itepo)->beg_epo << setw(16) << (*itepo)->end_epo << endl;
				}
			}
#endif
			if (_updmode == UPDTYPE::WL || _updmode == UPDTYPE::EWL
				|| _updmode == UPDTYPE::EWL24 || _updmode == UPDTYPE::EWL25)
			{
				/** Calculate EWL/WL ambiguities for UPD estimation */
				if (!t_gupdlsq::_calculateAmbWL(one_site_amb, OW)) return false;

				/** Connect EWL/WL ambiguities of a satellite-station pair */
				if (!t_gupdlsq::_connectAmb(*it_site, OW)) return false;
			}
			else if (_updmode == UPDTYPE::EWL_epoch)
			{
				all_site_amb[*it_site] = one_site_amb;
			}
		}

		t_gtime fixed_time = WL_IDENTIFY;
		if (_updmode == UPDTYPE::EWL) fixed_time = EWL_IDENTIFY;
		if (_updmode == UPDTYPE::EWL24) fixed_time = EWL24_IDENTIFY;
		if (_updmode == UPDTYPE::EWL25) fixed_time = EWL25_IDENTIFY;

		if (_updmode == UPDTYPE::WL || _updmode == UPDTYPE::EWL || _updmode == UPDTYPE::EWL24 || _updmode == UPDTYPE::EWL25)
		{
			/** Search reference station */
			if (!t_gupdlsq::_setRefStation(fixed_time)) return false;

			/** Get initial values of UPDs from UD-ambiguities */
			if (!t_gupdlsq::_initUPD(fixed_time)) return false;

			/** estimate upd of each satellite and station **/
			if (!t_gupdlsq::_updEstimate(fixed_time)) return false;

			//// temp, for BD, C06 (IGSO) as reference			
			//double ref = _cal_UPD->get_upd_value(WL_IDENTIFY, *_sats.begin());
			//double value;
			//auto itsat = _sats.begin();
			//for (itsat; itsat != _sats.end(); itsat++)
			//{
			//	value = _cal_UPD->get_upd_value(WL_IDENTIFY, *itsat) - ref;
			//	_cal_UPD->reset_upd_value(WL_IDENTIFY, *itsat, value);
			//}
			/** Delete site upd result */
			for (set<string>::iterator it_site = _sites.begin(); it_site != _sites.end(); ++it_site)
			{
				string site = *it_site;
				_cal_UPD->delete_upd(_updmode, fixed_time, site);
			}
		}
		else if (_updmode == UPDTYPE::EWL_epoch)
		{
			int iepo = 0;
			t_gtime epoch = beg;
			t_gtime pre_epoch;
			bool is_first = true;
			_cal_UPD->set_valid_beg(_updmode, epoch);
			while (epoch <= end)
			{
				iepo = iepo + 1;
				it_site = _sites.begin();
				pre_epoch.from_mjd(epoch.mjd(), epoch.sod() - intv, epoch.dsec());
				for (it_site; it_site != _sites.end(); it_site++)
				{
					auto it_sat = _sats.begin();
					t_all_sats_OW OW;
					shared_ptr<t_oneway_ambiguity> temp_single_ow;
					_cal_UPD->copy_upd(_updmode, pre_epoch, *it_site, epoch, is_first, true);

					t_epo_sat_amb one_site_amb_tmp;
					for (; it_sat != _sats.end(); it_sat++)
					{
						if (it_site == _sites.begin()) _cal_UPD->copy_upd(_updmode, pre_epoch, *it_sat, epoch, is_first, false);
						one_site_amb_tmp[iepo][*it_sat] = all_site_amb[*it_site][iepo][*it_sat];
						temp_single_ow = make_shared<t_oneway_ambiguity>();
						(*temp_single_ow).beg_epo = iepo;
						(*temp_single_ow).end_epo = iepo;
						(*temp_single_ow).sat = *it_sat;
						OW[*it_sat].push_back(temp_single_ow);
					}
					if (!t_gupdlsq::_calculateAmbWL(one_site_amb_tmp, OW)) return false;
					if (!t_gupdlsq::_connectAmb(*it_site, OW)) return false;
				}
				if (!t_gupdlsq::_setRefStation(epoch)) return false;
				if (!t_gupdlsq::_initUPD(epoch)) return false;
				if (!t_gupdlsq::_updEstimate(epoch)) return false;

				//// temp, for BD, C06 (IGSO) as reference			
				//double ref = _cal_UPD->get_upd_value(epoch, *_sats.begin());
				//double value;
				//auto itsat = _sats.begin();
				//for (itsat; itsat != _sats.end(); itsat++)
				//{
				//	value = _cal_UPD->get_upd_value(epoch, *itsat) - ref;
				//	_cal_UPD->reset_upd_value(epoch, *itsat, value);
				//}
				epoch.from_mjd(epoch.mjd(), epoch.sod() + intv, epoch.dsec());
				is_first = false;
			}
			/** Delete site upd result */
			for (set<string>::iterator it_site = _sites.begin(); it_site != _sites.end(); ++it_site)
			{
				string site = *it_site;
				_cal_UPD->delete_upd(_updmode, epoch - intv, site);
			}
			return true;
		}

#ifdef DEBUG
		/** Output WL/EWL Sats UPD */
		for (set<string>::iterator it_sat = _sats.begin(); it_sat != _sats.end(); ++it_sat)
		{
			string sat = *it_sat;
			std::cout << " " << sat << " " << right << setprecision(3) << setw(10) <<
				_cal_UPD->get_upd_value(_updmode, fixed_time, sat) << setw(10) << _cal_UPD->get_upd_sigma(_updmode, fixed_time, sat)
				<< setw(5) << _cal_UPD->get_upd(_updmode)[fixed_time][sat]->npoint << endl;
		}
#endif
		return true;
	}

	bool t_gupdlsq::_processNL(const t_gtime& beg, const t_gtime& end, const double &intv)
	{
		t_gtime epoch = beg;
		t_gtime pre_epoch;
		t_epoch_amb AMB;
		bool is_first = true;
		int nepo = 0;
		double factor0, factor1, a, N1, N2, N3;
		string sys_idx;
		_allfixedamb = new t_gallambupd();
		t_ambupd_data ambupd_data;
		for (const auto site : _sites) _allfixedamb->addAmbUpd(site, t_gambupd());

		///** Remove the impact of  TAI-tsys difference */
		//epoch = epoch + epoch.tai_tsys(t_gtime::GPS);
		// loop epo
		while (epoch <= end)
		{
#ifdef DEBUG
			std::cout << fixed << " Epoch Time " << setw(6) << epoch.mjd()
				<< setprecision(2) << fixed << setw(10) << epoch.sod() + 0.0 << ":" << "\r";
#endif
			// cerr << epoch.str_ymdhms("", false) << endl;
			pre_epoch.from_mjd(epoch.mjd(), epoch.sod() - intv, epoch.dsec());
			// loop site
			auto itsite = _sites.begin();
			for (itsite; itsite != _sites.end(); itsite++)
			{
				/** use last epoch results as current epoch init */
				_cal_UPD->copy_upd(_updmode, pre_epoch, *itsite, epoch, is_first, true);
				// loop sat
				auto itsat = _sats.begin();
				for (itsat; itsat != _sats.end(); itsat++)
				{
					/** Init sat_site pair ambiguities */
					_sat_site_amb[*itsite][*itsat] = 9.90;
					/** use last epoch results as current epoch init */
					if (itsite == _sites.begin()) _cal_UPD->copy_upd(_updmode, pre_epoch, *itsat, epoch, is_first, false);
					/** Assign value to AMB from _gambupd */
					if (!_allambupd->isAmbUpdValid(epoch, *itsite, *itsat))continue;
					if (_obs_type == OBSCOMBIN::IONO_FREE)
					{
						AMB[*itsite][*itsat].bc = _allambupd->getOneSiteOneEpoOneSatAmbUpd(epoch, *itsite, *itsat)->amb_if;
						AMB[*itsite][*itsat].bwi = _allambupd->getOneSiteOneEpoOneSatAmbUpd(epoch, *itsite, *itsat)->amb_wl;
						AMB[*itsite][*itsat].bw0 = 0;
						AMB[*itsite][*itsat].bw = _allambupd->getOneSiteOneEpoOneSatAmbUpd(epoch, *itsite, *itsat)->amb_wl;
						AMB[*itsite][*itsat].sbw = _allambupd->getOneSiteOneEpoOneSatAmbUpd(epoch, *itsite, *itsat)->amb_wl_sigma;
						AMB[*itsite][*itsat].epo = epoch;
					}
					else
					{
						sys_idx = *itsat;
						if (_sys != "GLO") sys_idx = (*itsat).substr(0, 1);
						a = _sys_wavelen[sys_idx]["L1"] / _sys_wavelen[sys_idx]["L2"];
						factor0 = 1 / (1 - pow(a, 2));
						factor1 = 1 / (pow(1 / a, 2) - 1);
						N1 = _allambupd->getOneSiteOneEpoOneSatAmbUpd(epoch, *itsite, *itsat)->amb_if;
						N2 = _allambupd->getOneSiteOneEpoOneSatAmbUpd(epoch, *itsite, *itsat)->amb_wl;
						N3 = _allambupd->getOneSiteOneEpoOneSatAmbUpd(epoch, *itsite, *itsat)->amb_l3;
						AMB[*itsite][*itsat].bc = factor0 * N1 - factor1 * N2;
						AMB[*itsite][*itsat].bwi = N1 / _sys_wavelen[sys_idx]["L1"] - N2 / _sys_wavelen[sys_idx]["L2"];
						AMB[*itsite][*itsat].bw0 = 0;
						AMB[*itsite][*itsat].bw = N1 / _sys_wavelen[sys_idx]["L1"] - N2 / _sys_wavelen[sys_idx]["L2"];
						AMB[*itsite][*itsat].sbw = _allambupd->getOneSiteOneEpoOneSatAmbUpd(epoch, *itsite, *itsat)->amb_wl_sigma;
						if (!double_eq(N3, 0)) 
							AMB[*itsite][*itsat].bewi = N2 / _sys_wavelen[sys_idx]["L2"] - N3 / _sys_wavelen[sys_idx]["L3"];
						AMB[*itsite][*itsat].epo = epoch;
					}
				}
				/** Calculate NL ambiguities. */
				if (!_calculateAmbNL(*itsite, epoch, AMB)) return false;
			}
			nepo++;
			if (nepo >= 1)   // temp for different test
			{
				if (_jref_site.empty())
				{
					/** Search reference station */
					if (!t_gupdlsq::_setRefStation(epoch)) return false;
				}
				if (!_jref_site.empty() || !_iref_site.empty())
				{
					/** Set NL UPD valid begin time */
					_cal_UPD->set_valid_beg(_updmode, epoch);
					/** Get initial values of UPDs from UD-ambiguities */
					if (!t_gupdlsq::_initUPD(epoch)) return false;
					/** estimate upd of each satellite and station **/
					if (!t_gupdlsq::_updEstimate(epoch)) return false;

					for (const auto &site : _sites) {
						for (const auto &sat : _sats) {
							string str_idx;
							str_idx = sat.substr(0, 1) == "R" ? sat : sat.substr(0, 1);
							if (_AMB_RES[site][sat].wl_fixed && _AMB_RES[site][sat].nl_fixed) {
								double nlupd_s = _cal_UPD->get_upd_value(_updmode, epoch, sat);
								double nlupd_r = _cal_UPD->get_upd_value(_updmode, epoch, site);
								double wlupd_s = _allWLUPD->get_upd_value(UPDTYPE::WL, WL_IDENTIFY, sat);
								double wlupd_r = _allWLUPD->get_upd_value(UPDTYPE::WL, WL_IDENTIFY, site);
								double upd1 = -1.0 * (nlupd_s + nlupd_r) - _sys_wavelen[str_idx]["L1"] / (_sys_wavelen[str_idx]["L2"] - _sys_wavelen[str_idx]["L1"]) * (wlupd_s + wlupd_r);
								double upd2 = -1.0 * (nlupd_s + nlupd_r) - _sys_wavelen[str_idx]["L2"] / (_sys_wavelen[str_idx]["L2"] - _sys_wavelen[str_idx]["L1"]) * (wlupd_s + wlupd_r);
								double amb_l1 = -1.0 * (AMB[site][sat].bn - nlupd_s - nlupd_r);
								ambupd_data.amb_if = round(amb_l1);
								ambupd_data.amb_wl = ambupd_data.amb_if - AMB[site][sat].bwi;
								if (_AMB_RES[site][sat].ewl_fixed) {
									double ewlupd_s = _allWLUPD->get_upd_value(UPDTYPE::EWL, EWL_IDENTIFY, sat);
									double ewlupd_r = _allWLUPD->get_upd_value(UPDTYPE::EWL, EWL_IDENTIFY, site);
									double upd3 = upd2 - ewlupd_s - ewlupd_r;
									ambupd_data.amb_l3 = ambupd_data.amb_wl - AMB[site][sat].bewi;
									//ambupd_data.amb_l3 += upd3;
								}
								//ambupd_data.amb_if += upd1;
								//ambupd_data.amb_wl += upd2;
								_allfixedamb->addAmbUpd(site, epoch, sat, ambupd_data);
							}
						}
					}
#ifdef DEBUG
					/** Output Sats UPD */
					stringstream os1, os2;
					os1 << "SAT_UN X" << fixed << right << setprecision(1) << setw(8) << epoch.sod() + 0.0;
					os2 << "SAT_UN S" << fixed << right << setprecision(1) << setw(8) << epoch.sod() + 0.0;
					for (set<string>::iterator it_sat = _sats.begin(); it_sat != _sats.end(); ++it_sat)
					{
						string sat = *it_sat;
						os1 << fixed << right << setprecision(1) << setw(10) << _cal_UPD->get_upd_value(_updmode, epoch, sat);
						os2 << fixed << right << setprecision(1) << setw(10) << _cal_UPD->get_upd_sigma(_updmode, epoch, sat);
					}
					cout << os1.str() << endl;
					cout << os2.str() << endl;
					for (set<string>::iterator it_sat = _sats.begin(); it_sat != _sats.end(); ++it_sat)
					{
						string sat = *it_sat;
						std::cout << " " << sat << " " << right << setprecision(3) << setw(10) <<
							_cal_UPD->get_upd_value(_updmode, epoch, sat) << setw(10) << _cal_UPD->get_upd_sigma(_updmode, epoch, sat)
							<< setw(5) << _cal_UPD->get_upd(_updmode)[epoch][sat]->npoint << endl;
					}
#endif
					/** Write amb_res file */
					_writeAmbResidual(epoch);//_AMB_RES
				}
			}
			/** Increase time */
			epoch.from_mjd(epoch.mjd(), epoch.sod() + intv, epoch.dsec());
			is_first = false;
		}
		if (_carrier_range && _allfixedamb && _allambflag) {
			_writeCarrierRange(beg, intv);
		}	
		delete _allfixedamb; _allfixedamb = nullptr;
		/** Delete site upd result */
		for (set<string>::iterator it_site = _sites.begin(); it_site != _sites.end(); ++it_site)
		{
			string site = *it_site;
			_cal_UPD->delete_upd(_updmode, epoch - intv, site);
		}
		return true;
	}

	bool t_gupdlsq::_processIFCB(const t_gtime& beg, const t_gtime& end, const double &intv)
	{
		//_sites= dynamic_cast<t_gsetgen*>(_set)->recs();
		auto it_site = _sites.begin();
		map <string, t_epo_sat_amb> all_site_GFIF;
		READ_OFILE_MODE rom = dynamic_cast<t_gsetproc*>(_set)->read_ofile_mode();
		for (it_site; it_site != _sites.end(); it_site++)
		{
			t_gallobs  obs;  obs.glog(_glog); obs.gset(_set);
			if (rom == READ_OFILE_MODE::REALTIME) {
					_readOfile(it_site, obs);//read current site's obs
					_allobs = &obs;
			}
			
			auto it_sat = _sats.begin();
			for (; it_sat != _sats.end(); it_sat++)
			{
				_sat_site_amb[*it_site][*it_sat] = 9.90;
			}

			t_all_sats_OW OW;
			t_epo_sat_amb one_site_GFIF;

			if (!t_gupdlsq::_combineObsGFIF(*it_site, beg, end, intv, one_site_GFIF, OW)) return false;

			auto itsat = OW.begin();
			for (itsat; itsat != OW.end(); itsat++)
			{
				auto itepo = itsat->second.begin();
				for (itepo; itepo != itsat->second.end(); itepo++)
				{
					int iepo = (*itepo)->beg_epo;
					all_site_GFIF[*it_site][iepo][itsat->first] = 0.0;
					for (; iepo < (*itepo)->end_epo; iepo++)
					{
						if (double_eq(one_site_GFIF[iepo][itsat->first], 0.0) || double_eq(one_site_GFIF[iepo + 1][itsat->first], 0.0))
							all_site_GFIF[*it_site][iepo + 1][itsat->first] = 0.0;
						else
							all_site_GFIF[*it_site][iepo + 1][itsat->first] = one_site_GFIF[iepo + 1][itsat->first] - one_site_GFIF[iepo][itsat->first];
					}
				}
			}
		}

		int iepo = 0;
		t_gtime epoch = beg;
		t_gtime pre_epoch;
		_cal_UPD->set_valid_beg(_updmode, epoch);
		while (epoch <= end)
		{
			iepo = iepo + 1;
			it_site = _sites.begin();
			pre_epoch.from_mjd(epoch.mjd(), epoch.sod() - intv, epoch.dsec());
			for (; it_site != _sites.end(); it_site++)
			{
				auto it_sat = _sats.begin();
				t_all_sats_OW OW;
				shared_ptr<t_oneway_ambiguity> temp_single_ow;

				t_epo_sat_amb one_site_GFIF_tmp;
				for (; it_sat != _sats.end(); it_sat++)
				{
					_sat_site_amb[*it_site][*it_sat] = 9.90;
					_cal_UPD->copy_upd(_updmode, pre_epoch, *it_sat, epoch, true, false);

					one_site_GFIF_tmp[iepo][*it_sat] = all_site_GFIF[*it_site][iepo][*it_sat];
					temp_single_ow = make_shared<t_oneway_ambiguity>();
					(*temp_single_ow).beg_epo = iepo;
					(*temp_single_ow).end_epo = iepo;
					(*temp_single_ow).sat = *it_sat;
					OW[*it_sat].push_back(temp_single_ow);
				}
				if (!t_gupdlsq::_calculateAmbWL(one_site_GFIF_tmp, OW)) return false;
				if (!t_gupdlsq::_connectAmb(*it_site, OW)) return false;
			}

			auto it_sat = _sats.begin();
			for (it_sat; it_sat != _sats.end(); it_sat++)
			{
				list<pair<double, double>> wx;
				double mean, sigma, mean_sig;
				it_site = _sites.begin();
				for (it_site; it_site != _sites.end(); it_site++)
				{
					if (double_eq(_sat_site_amb[*it_site][*it_sat], 9.90)) continue;
					wx.push_back(make_pair(_sat_site_amb[*it_site][*it_sat], 1.0));
				}
				for (int i = 0; i < min(10.0, (double)(wx.size()) / 2.0); i++)
				{
					int wx_size = wx.size();
					getMeanWgt(true, wx, mean, sigma, mean_sig);
					if (wx_size == wx.size() || mean_sig == 0.0) break;
				}
				if (wx.size() > 0)
				{
					_cal_UPD->reset_upd(_updmode, epoch, *it_sat, mean, sigma, wx.size());
				}
				else
				{
					_cal_UPD->re_init_upd(_updmode, epoch, *it_sat);
				}
			}

			epoch.from_mjd(epoch.mjd(), epoch.sod() + intv, epoch.dsec());
		}
		
		// a zero-mean priori constraint
		for (auto itsat = _sats.begin(); itsat != _sats.end(); itsat++)
		{
			int nepo = ((end.mjd() - beg.mjd()) * 86400 + end.sod() - beg.sod()) / intv + 1;
			double upd_init = 0;
			epoch = beg;
			for (int iepo = 1; iepo < nepo; iepo++)
			{
				epoch.from_mjd(epoch.mjd(), epoch.sod() + intv, epoch.dsec());
				double value = _cal_UPD->get_upd_value(_updmode, epoch, *itsat);
				upd_init += (nepo - iepo)*value;
			}
			_cal_UPD->reset_upd(_updmode, beg, *itsat, -upd_init / nepo, 0.0, 999);
		}

		epoch = beg;
		while (epoch < end)
		{
			epoch.from_mjd(epoch.mjd(), epoch.sod() + intv, epoch.dsec());
			auto itsat = _sats.begin();
			double value1, value2;
			for (itsat; itsat != _sats.end(); itsat++)
			{
				value1 = _cal_UPD->get_upd_value(_updmode, epoch - intv, *itsat);
				value2 = _cal_UPD->get_upd_value(_updmode, epoch, *itsat);
				_cal_UPD->reset_upd_value(_updmode, epoch, *itsat, value1 + value2);
			}
		}

		for (set<string>::iterator it_site = _sites.begin(); it_site != _sites.end(); ++it_site)
		{
			string site = *it_site;
			_cal_UPD->delete_upd(_updmode, epoch, site);
		}

		return true;
	}

	//string t_gupdlsq::UPDType2String(UPDTYPE type)
	//{
	//	switch (type) {
	//	case UPDTYPE::EWL:  return "EWL";
	//	case UPDTYPE::WL:  return "WL";
	//	case UPDTYPE::NL:  return "NL";
	//	case UPDTYPE::EWL: 
	//	default:             return "UNDEF";
	//	}
	//	return "UNDEF";
	//}

	bool t_gupdlsq::_setWavelength()
	{
		t_gobsgnss* gnss = new t_gobsgnss();
		GSYS sys = t_gsys::str2gsys(_sys);
		GOBSBAND b1, b2, b3;
		b1 = _band_index[sys][FREQ_1];
		b2 = _band_index[sys][FREQ_2];
		string sys_short;
		switch (sys)
		{
		case GPS:
			sys_short = "G"; break;
		case GAL:
			sys_short = "E"; break;
		case BDS:
			sys_short = "C"; break;
		case GLO:
			sys_short = "R"; break;
		case QZS:
			sys_short = "J"; break;
		default:
			sys_short = "NONE"; break;
		}

		if (sys_short == "NONE")
		{
			write_log_info(_glog, 1, " ERROR[t_gupdlsq::_setWavelength] : Unrecognized system (GCREJ) : " + _sys);
			return false;
		}

		if (_sys == "GLO")
		{
			for (auto sat_iter = _glo_freq_num.begin(); sat_iter != _glo_freq_num.end(); sat_iter++)
			{
				string sat = sat_iter->first;
				gnss->sat(sat);
				gnss->channel(sat_iter->second);
				_sys_wavelen[sat]["L1"] = gnss->wavelength(b1);
				_sys_wavelen[sat]["L2"] = gnss->wavelength(b2);
				_sys_wavelen[sat]["WL"] = gnss->wavelength_WL(b1, b2);
				_sys_wavelen[sat]["NL"] = gnss->wavelength_NL(b1, b2);
			}
		}
		else
		{
			gnss->sat(sys_short);
			_sys_wavelen[sys_short]["L1"] = gnss->wavelength(b1);
			_sys_wavelen[sys_short]["L2"] = gnss->wavelength(b2);
			_sys_wavelen[sys_short]["WL"] = gnss->wavelength_WL(b1, b2);
			_sys_wavelen[sys_short]["NL"] = gnss->wavelength_NL(b1, b2);

			// if (_updtype != UPDTYPE::NL && _updtype != UPDTYPE::WL)
			if (_updmode != UPDTYPE::WL) // jqwu: for fixing L3 ambiguity in NL UPD estimation
			{
				b3 = _band_index[sys][FREQ_3];
				_sys_wavelen[sys_short]["L3"] = gnss->wavelength(b3);
				_sys_wavelen[sys_short]["EWL"] = gnss->wavelength_WL(b2, b3);
			}

		}

		delete gnss; gnss = nullptr;
		return true;
	}

	bool t_gupdlsq::_combineObsIF(string site, const t_gtime& T_beg, const t_gtime& T_end, double intv,
		t_epo_sat_amb& one_site_amb, t_all_sats_OW& OW)
	{
		//local variables
		int num_epo;
		string site_lower;
		vector<t_gsatdata> obsdata;
		t_vec_ambflag one_sat_flag;
		shared_ptr<t_oneway_ambiguity> temp_single_ow;
		t_gambflag one_site_flag, one_site_flag13;

		double mw_obs;
		double correct;
		GOBSBAND b1, b2, b3;
		t_gobs gobs1, gobs2, gobs3, gobs4, gobs5;
		try
		{
			//get number of epochs
			num_epo = round((T_end - T_beg) / intv);

			// get one site ambflag
			site_lower = "";
			transform(site.begin(), site.end(), back_inserter(site_lower), ::tolower);
			one_site_flag = _allambflag->getOneAmbFlag(site_lower);
			if (_updmode == UPDTYPE::EWL || _updmode == UPDTYPE::EWL_epoch) one_site_flag13 = _allambflag13->getOneAmbFlag(site_lower);
			if (_updmode == UPDTYPE::EWL24) one_site_flag13 = _allambflag14->getOneAmbFlag(site_lower);
			if (_updmode == UPDTYPE::EWL25) one_site_flag13 = _allambflag15->getOneAmbFlag(site_lower);

			t_gtime it_time = T_beg;
			///** Remove the impact of  TAI-tsys difference */
			//t_gtime it_gpstime = it_time + it_time.tai_tsys(t_gtime::GPS);
			//loop over epoch
			for (int iepo = 1; iepo <= num_epo; iepo++)
			{
				obsdata = _allobs->obs(site, it_time);

				//_outliers_proc.excludeBadObs(obsdata);
				if (_sys == "BDS")
				{
					t_gtriple BDS_tmp(0.0,0.0,0.0);
					_bds_codebias_cor.apply_IGSO_MEO(site, BDS_tmp, _allnav, obsdata);
					if (obsdata.size() > 0 && BDS_tmp.zero())
					{
						BDS_tmp = _allobs->getsitecrd(site);
						_bds_codebias_cor.apply_IGSO_MEO(site, BDS_tmp, _allnav, obsdata);
					}
				}

				//loop over sats
				for (auto it_satdata = obsdata.begin(); it_satdata != obsdata.end(); it_satdata++)
				{
					string sat = it_satdata->sat();
					one_site_amb[iepo][sat] = 0.0;
					if (sat.substr(0, 1) == "R" && it_satdata->channel() >= DEF_CHANNEL) it_satdata->channel(_glo_freq_num[sat]);

					// apply DCB on observation ~~
					it_satdata->apply_bias(*_allbias);

					//put amb arc into OW
					int AMB_pos;
					if (_updmode == UPDTYPE::WL)
					{
						if (!one_site_flag.isValid(sat, it_time, AMB_pos)) continue;
						one_sat_flag = one_site_flag.getAmbFlagData()[sat];
						if (one_sat_flag[AMB_pos]->iflag != "0")
						{
							temp_single_ow = make_shared<t_oneway_ambiguity>();
							(*temp_single_ow).beg_epo = iepo;
							(*temp_single_ow).sat = sat;
							OW[sat].push_back(temp_single_ow);
							one_site_flag.reset_iflag(sat, "0", AMB_pos);
							_lastSat = sat;
						}
						if (OW[sat].size() != 0) OW[sat].back()->end_epo = iepo;
						//b1 = t_gsys::band_priority(it_satdata->gsys(), FREQ_1);
						//b2 = t_gsys::band_priority(it_satdata->gsys(), FREQ_2);
						b1 =  _band_index[it_satdata->gsys()][FREQ_1];
						b2 =  _band_index[it_satdata->gsys()][FREQ_2];

						gobs1 = t_gobs(it_satdata->select_range(b1));
						gobs2 = t_gobs(it_satdata->select_range(b2));
						gobs3 = t_gobs(it_satdata->select_phase(b1));
						gobs4 = t_gobs(it_satdata->select_phase(b2));

						mw_obs = it_satdata->MW_cycle(gobs3, gobs4, gobs1, gobs2);
						cout << "MW" << " " << iepo << " " << sat << " " << mw_obs << endl;
					}
					else
					{
						int AMB_pos13;
						t_vec_ambflag one_sat_flag13;

						if (!one_site_flag.isValid(sat, it_time, AMB_pos) ||
							!one_site_flag13.isValid(sat, it_time, AMB_pos13)) continue;
						one_sat_flag = one_site_flag.getAmbFlagData()[sat];
						one_sat_flag13 = one_site_flag13.getAmbFlagData()[sat];

						if (one_sat_flag[AMB_pos]->iflag != "0" || one_sat_flag13[AMB_pos13]->iflag != "0")
						{
							temp_single_ow = make_shared<t_oneway_ambiguity>();
							(*temp_single_ow).beg_epo = iepo;
							(*temp_single_ow).sat = sat;
							OW[sat].push_back(temp_single_ow);
							one_site_flag.reset_iflag(sat, "0", AMB_pos);
							one_site_flag13.reset_iflag(sat, "0", AMB_pos13);
							_lastSat = sat;
						}
						if (OW[sat].size() != 0) OW[sat].back()->end_epo = iepo;
						//b1 = t_gsys::band_priority(it_satdata->gsys(), FREQ_1);
						//b2 = t_gsys::band_priority(it_satdata->gsys(), FREQ_2);
						//b3 = t_gsys::band_priority(it_satdata->gsys(), FREQ_3);
						b1 = _band_index[it_satdata->gsys()][FREQ_1];
						b2 = _band_index[it_satdata->gsys()][FREQ_2];
						b3 = _band_index[it_satdata->gsys()][FREQ_3];

						gobs1 = t_gobs(it_satdata->select_range(b1));
						gobs2 = t_gobs(it_satdata->select_range(b2));
						gobs3 = t_gobs(it_satdata->select_phase(b1));
						gobs4 = t_gobs(it_satdata->select_phase(b2));
						gobs5 = t_gobs(it_satdata->select_phase(b3));

						mw_obs = it_satdata->EWL_cycle(gobs3, gobs4, gobs5, gobs1, gobs2);
					}

					if (double_eq(mw_obs, 0.0) || std::isnan(mw_obs)) continue;

					//dcb correct
					correct = 0.0;
					bool isSuccess;
					if (_updmode == UPDTYPE::WL)
					{
						//isSuccess = apply_DCB(_allbias, *it_satdata, correct, OBSCOMBIN::MW_COMBIN, &gobs1, &gobs2);
						//isSuccess = correct_DCB("CAS", _allbias, *it_satdata, correct, OBSCOMBIN::MW_COMBIN, &gobs1, &gobs2);
					}
					else if (_updmode == UPDTYPE::EWL || _updmode == UPDTYPE::EWL_epoch)
					{
						//isSuccess = apply_DCB(_allbias, *it_satdata, correct, OBSCOMBIN::EWL_COMBIN, &gobs1, &gobs2);
						//isSuccess = correct_DCB("CAS", _allbias, *it_satdata, correct, OBSCOMBIN::EWL_COMBIN, &gobs1, &gobs2);
						if (sat.substr(0, 1) == "G")
						{
							one_epoch_ifcb epoch_ifcb = _gifcb->get_epo_ifcb(it_satdata->epoch());
							if (epoch_ifcb.find(sat) == epoch_ifcb.end()) continue;

							if (_gifcb->get_ifcb_sigma(it_satdata->epoch(), sat) <= 0.2 &&
								_gifcb->get_ifcb_npoint(it_satdata->epoch(), sat) > 2)
							{
								double factor = -3.113;
								correct -= factor * _gifcb->get_ifcb_value(it_satdata->epoch(), sat);
							}
						}
					}
				//	if (isSuccess)
					//{
#ifdef DEBUG
					std::cout << " WL/EWL Value " << " " << iepo << " " << sat << " " << fixed << setw(10)
						<< setprecision(6) << correct << setw(15) << setprecision(8) << mw_obs << endl;
#endif
						mw_obs += correct;
				//	}

					//赋值
						//赋值
					one_site_amb[iepo][sat] = mw_obs;

				}//end sat

				it_time = it_time + intv;
			}//end epo

			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 1, " ERROR[t_gupdlsq::_combineObsIF] : unknown mistake !");
			return false;
		}
	}

	bool t_gupdlsq::_combineObsALL(string site, const t_gtime& T_beg, const t_gtime& T_end, double intv,
		t_epo_sat_amb& one_site_amb, t_all_sats_OW& OW)
	{
		//local variables
		int num_epo;
		double correct = 0;
		shared_ptr<t_oneway_ambiguity> temp_single_ow;

		try
		{
			//get number of epochs
			num_epo = round((T_end - T_beg) / intv);

			t_gtime it_time = T_beg;

			//loop over epoch
			for (int iepo = 1; iepo <= num_epo; iepo++)
			{
				//loop over sats
				auto itsat = _sats.begin();
				for (itsat; itsat != _sats.end(); itsat++)
				{
					string sat = *itsat;
					one_site_amb[iepo][sat] = 0.0;

					if (!_allambupd->isAmbUpdValid(it_time, site, sat)) continue;

					//put amb arc into OW
					if (OW.find(sat) == OW.end())
					{
						temp_single_ow = make_shared<t_oneway_ambiguity>();
						(*temp_single_ow).beg_epo = iepo;
						(*temp_single_ow).end_epo = iepo;
						(*temp_single_ow).sat = sat;
						OW[sat].push_back(temp_single_ow);
					}
					else
					{
						if (OW[sat].back()->end_epo + 1 < iepo)
						{
							temp_single_ow = make_shared<t_oneway_ambiguity>();
							(*temp_single_ow).beg_epo = iepo;
							(*temp_single_ow).end_epo = iepo;
							(*temp_single_ow).sat = sat;
							OW[sat].push_back(temp_single_ow);
						}
						else
						{
							OW[sat].back()->end_epo = iepo;
						}
					}

					double l1, l2, N1, N2;
					string sys_idx = sat;
					if (_sys != "GLO") sys_idx = sat.substr(0, 1);

					l1 = _sys_wavelen[sys_idx]["L2"];
					l2 = _sys_wavelen[sys_idx]["L3"];

					if (_updmode == UPDTYPE::EWL || _updmode == UPDTYPE::EWL_epoch)
					{
						N1 = _allambupd->getOneSiteOneEpoOneSatAmbUpd(it_time, site, sat)->amb_wl / l1;
						N2 = _allambupd->getOneSiteOneEpoOneSatAmbUpd(it_time, site, sat)->amb_l3 / l2;
					}
					else if (_updmode == UPDTYPE::EWL24)
					{
						N1 = _allambupd->getOneSiteOneEpoOneSatAmbUpd(it_time, site, sat)->amb_wl / l1;
						N2 = _allambupd->getOneSiteOneEpoOneSatAmbUpd(it_time, site, sat)->amb_l4 / l2;
					}
					else if (_updmode == UPDTYPE::EWL25)
					{
						N1 = _allambupd->getOneSiteOneEpoOneSatAmbUpd(it_time, site, sat)->amb_wl / l1;
						N2 = _allambupd->getOneSiteOneEpoOneSatAmbUpd(it_time, site, sat)->amb_l5 / l2;
					}
					else
					{
						l1 = _sys_wavelen[sys_idx]["L1"];
						l2 = _sys_wavelen[sys_idx]["L2"];
						N1 = _allambupd->getOneSiteOneEpoOneSatAmbUpd(it_time, site, sat)->amb_if / l1;
						N2 = _allambupd->getOneSiteOneEpoOneSatAmbUpd(it_time, site, sat)->amb_wl / l2;
					}
									

					if (N1 != 0 && N2 != 0)
					{
						one_site_amb[iepo][*itsat] = N1 - N2;
#ifdef DEBUG
						std::cout << " WL/EWL Value " << " " << iepo << " " << sat << " " << fixed << setw(15) << setprecision(8) << N1 - N2 << endl;
#endif
					}
				}//end sat

				it_time = it_time + intv;
			}//end epo

			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 1, " ERROR[t_gupdlsq::_combineObsALL] : unknown mistake !");
			return false;
		}
	}

	bool t_gupdlsq::_combineObsGFIF(string site, const t_gtime& T_beg, const t_gtime& T_end, double intv,
		t_epo_sat_amb& one_site_GFIF, t_all_sats_OW& OW)
	{
		int num_epo;
		string site_lower;
		vector<t_gsatdata> obsdata;
		t_vec_ambflag one_sat_flag, one_sat_flag13;
		shared_ptr<t_oneway_ambiguity> temp_single_ow;
		t_gambflag one_site_flag, one_site_flag13;

		double GFIF_obs;
		GOBSBAND b1, b2, b3;
		t_gobs gobs1, gobs2, gobs3;
		try
		{
			num_epo = round((T_end - T_beg) / intv);
			site_lower = "";
			transform(site.begin(), site.end(), back_inserter(site_lower), ::tolower);
			one_site_flag = _allambflag->getOneAmbFlag(site_lower);
			one_site_flag13 = _allambflag13->getOneAmbFlag(site_lower);

			t_gtime it_time = T_beg;
			for (int iepo = 1; iepo <= num_epo; iepo++)
			{
				obsdata = _allobs->obs(site, it_time);

				//_outliers_proc.excludeBadObs(obsdata);
				if (_sys == "BDS")
				{
					t_gtriple BDS_tmp(0.0, 0.0, 0.0);
					_bds_codebias_cor.apply_IGSO_MEO(site, BDS_tmp, _allnav, obsdata);
					if (obsdata.size() > 0 && BDS_tmp.zero())
					{
						BDS_tmp = _allobs->getsitecrd(site);
						_bds_codebias_cor.apply_IGSO_MEO(site, BDS_tmp, _allnav, obsdata);
					}
				}

				for (auto it_satdata = obsdata.begin(); it_satdata != obsdata.end(); it_satdata++)
				{
					it_satdata->apply_bias(*_allbias);
					
					string sat = it_satdata->sat();
					one_site_GFIF[iepo][sat] = 0.0;

					int AMB_pos;
					int AMB_pos13;

					if (!one_site_flag.isValid(sat, it_time, AMB_pos) ||
						!one_site_flag13.isValid(sat, it_time, AMB_pos13)) continue;
					one_sat_flag = one_site_flag.getAmbFlagData()[sat];
					one_sat_flag13 = one_site_flag13.getAmbFlagData()[sat];

					if (one_sat_flag[AMB_pos]->iflag != "0" || one_sat_flag13[AMB_pos13]->iflag != "0")
					{
						temp_single_ow = make_shared<t_oneway_ambiguity>();
						(*temp_single_ow).beg_epo = iepo;
						(*temp_single_ow).sat = sat;
						OW[sat].push_back(temp_single_ow);
						one_site_flag.reset_iflag(sat, "0", AMB_pos);
						one_site_flag13.reset_iflag(sat, "0", AMB_pos13);
						_lastSat = sat;
					}
					if (OW[sat].size() != 0) OW[sat].back()->end_epo = iepo;
					//b1 = t_gsys::band_priority(it_satdata->gsys(), FREQ_1);
					//b2 = t_gsys::band_priority(it_satdata->gsys(), FREQ_2);
					//b3 = t_gsys::band_priority(it_satdata->gsys(), FREQ_3);
					b1 = _band_index[it_satdata->gsys()][FREQ_1];
					b2 = _band_index[it_satdata->gsys()][FREQ_2];
					b3 = _band_index[it_satdata->gsys()][FREQ_3];

					gobs1 = t_gobs(it_satdata->select_phase(b1));
					gobs2 = t_gobs(it_satdata->select_phase(b2));
					gobs3 = t_gobs(it_satdata->select_phase(b3));

					GFIF_obs = it_satdata->GFIF_meter(gobs1, gobs2, gobs3);

					if (double_eq(GFIF_obs, 0.0)) continue;

					one_site_GFIF[iepo][sat] = GFIF_obs;
				}
				it_time = it_time + intv;
			}
			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 1, " ERROR[t_gupdlsq::_combineObsGFIF] : unknown mistake !");
			return false;
		}
	}

	bool t_gupdlsq::_calculateAmbWL(t_epo_sat_amb& one_site_amb, t_all_sats_OW& OW)
	{
		list<pair<double, double>> wx;
		double mean, sigma, mean_sig;
		int npos = 0;
#ifdef DEBUG
		std::cout << "      ista isat nobs mean_value     sigma  iepo1  iepo2" << endl;
#endif
		try
		{
			auto itsat = OW.begin();
			// loop sats
			for (itsat; itsat != OW.end(); itsat++)
			{
				auto itepo = itsat->second.begin();
				// each sat has several arcs, loop arcs
				for (itepo; itepo != itsat->second.end(); )
				{
					wx.clear();
					npos = distance(itsat->second.begin(), itepo);  // index of arc in one sat
					OW[itsat->first][npos]->rwl = 0.0;
					// each arc has several epoch, loop epo
					for (int iepo = (*itepo)->beg_epo; iepo <= (*itepo)->end_epo; iepo++)
					{
						// one_site_amb means wl/ewl
						//cerr << one_site_amb[iepo][itsat->first] << endl;
						if (double_eq(one_site_amb[iepo][itsat->first], 0.0))  continue;
						wx.push_back(make_pair(one_site_amb[iepo][itsat->first], 1.0));
					}
					// Calculate mean as well as sigma
					getMeanWgt(false, wx, mean, sigma, mean_sig);
					// Select OK wl/ewl ambiguity according to threshold
					if ((_updmode == UPDTYPE::WL && (wx.size() <= 15 || mean_sig > 0.15))
						|| (_updmode == UPDTYPE::EWL && (wx.size() <= 15 || mean_sig > 0.15))
						|| (_updmode == UPDTYPE::EWL24 && (wx.size() <= 15 || mean_sig > 0.15))
						|| (_updmode == UPDTYPE::EWL25 && (wx.size() <= 15 || mean_sig > 0.15))
						|| (_updmode == UPDTYPE::EWL_epoch && wx.size() <= 0)
						|| (_updmode == UPDTYPE::IFCB && wx.size() <= 0))
					{
						itepo = itsat->second.erase(itepo);
					}
					else
					{
						OW[itsat->first][npos]->rwl = mean;
						OW[itsat->first][npos]->srwl = mean_sig;
						if (_obs_type == OBSCOMBIN::RAW_ALL) OW[itsat->first][npos]->srwl += 0.001;
#ifdef DEBUG
						std::cout <<updmode2str(_updmode) << " " << itsat->first << setw(16) << fixed << setprecision(2) << (*itepo)->rwl
							<< setw(16) << fixed << setprecision(4) << (*itepo)->srwl << setw(16) << (*itepo)->beg_epo << setw(16) << (*itepo)->end_epo << endl;
#endif
						itepo++;
					}
				}
			}
			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 1, " ERROR[t_gupdlsq::_calculateAmbWL] : unknown mistake !");
			return false;
		}
	}

	bool t_gupdlsq::_calculateAmbNL(string site, const t_gtime& t, t_epoch_amb& AMB)
	{
		map<string, double> xwl;
		map<string, double> wgt;
		list<pair<double, double>> wl_sig;
		list<pair<double, double>> ewl_sig;
		map<string, bool> isfixed;
		double sigx = 0.0;
		double mean = 0.0;
		double sig = 0.0;
		int npnt = 0;
		int ndel = 0;
		try
		{
			//get wide-lane ambiguities which are corrected with satellite UPDs
			auto itsat = AMB[site].begin();
			for (itsat; itsat != AMB[site].end(); itsat++)
			{
				// skip sat that unvalid in Wl UPD glfeng
				if (_allWLUPD->get_upd_sigma(UPDTYPE::WL, WL_IDENTIFY, itsat->first) > 1000.0) continue;
				/*skip past ambiguities(300 ago).Be aware ambiguities are kept if no new one comes in for each satellite*/
				double dt = (itsat->second.epo.dmjd() - t.dmjd()) * 86400;
				if (abs(dt) > 30000.0) continue;
				if (_site_ambupd_end[site] < t) continue;   // add for end of ambupd file glfeng

				xwl[itsat->first] = itsat->second.bw0 + itsat->second.bw - _allWLUPD->get_upd_value(UPDTYPE::WL, WL_IDENTIFY, itsat->first)
					- _allWLUPD->get_upd_value(UPDTYPE::WL, WL_IDENTIFY, site);
				wgt[itsat->first] = 1.0;
				wl_sig.push_back(make_pair(xwl[itsat->first], wgt[itsat->first]));
				npnt++;
			}
			if (npnt > 0)
			{
				getMeanFract(wl_sig, mean, sig, sigx);
				if (sigx <= 0.1)
				{
					double wl_rec = _allWLUPD->get_upd_value(UPDTYPE::WL, WL_IDENTIFY, site);
					if (!double_eq(wl_rec, 0.0)) mean = getFraction(mean, -0.5);
					_allWLUPD->reset_upd(UPDTYPE::WL, WL_IDENTIFY, site, wl_rec + mean, sigx, wl_sig.size());

#ifdef DEBUG
					std::cout << setprecision(2) << fixed << " Rec.WL - UPD:(Mean&Sigma) :" << site
						<< setw(4) << npnt << setw(4) << ndel << setw(10) << wl_rec + mean
						<< setw(10) << sigx << setw(10) << mean << endl;
#endif
				}
			}

			//get extra-wide-lane ambiguities which are corrected with satellite UPDs
			if (_allambflag13 && _carrier_range) {
				npnt = 0;
				for (auto itsat = AMB[site].begin(); itsat != AMB[site].end(); ++itsat) {
					// skip sat that unvalid in EWl UPD glfeng
					if (_allWLUPD->get_upd_sigma(UPDTYPE::EWL, EWL_IDENTIFY, itsat->first) > 1000.0) continue;
					/*skip past ambiguities(300 ago).Be aware ambiguities are kept if no new one comes in for each satellite*/
					double dt = (itsat->second.epo.dmjd() - t.dmjd()) * 86400;
					if (abs(dt) > 30000.0) continue;
					if (_site_ambupd_end[site] < t) continue;   // add for end of ambupd file glfeng
					if (double_eq(itsat->second.bewi, 0)) continue; // no float ewl ambiguity
					xwl[itsat->first] = itsat->second.bewi - _allWLUPD->get_upd_value(UPDTYPE::EWL, EWL_IDENTIFY, itsat->first)
						- _allWLUPD->get_upd_value(UPDTYPE::EWL, EWL_IDENTIFY, site);
					wgt[itsat->first] = 1.0;
					ewl_sig.push_back(make_pair(xwl[itsat->first], wgt[itsat->first]));
					npnt++;
				}
				if (npnt > 0) {
					getMeanFract(ewl_sig, mean, sig, sigx);
					if (sigx <= 0.1) {
						double wl_rec = _allWLUPD->get_upd_value(UPDTYPE::EWL, EWL_IDENTIFY, site);
						if (!double_eq(wl_rec, 0.0)) mean = getFraction(mean, -0.5);
						_allWLUPD->reset_upd(UPDTYPE::EWL, EWL_IDENTIFY, site, wl_rec + mean, sigx, ewl_sig.size());
					}
				}
				for (auto itsat = AMB[site].begin(); itsat != AMB[site].end(); ++itsat) {
					double dt = (itsat->second.epo.dmjd() - t.dmjd()) * 86400;
					if (abs(dt) > 30000.0) continue;
					if (_site_ambupd_end[site] < t) continue;
					if (double_eq(itsat->second.bewi, 0)) continue; // no float ewl ambiguity
					string sat = itsat->first;
					double ebw = AMB[site][sat].bewi - _allWLUPD->get_upd_value(UPDTYPE::EWL, EWL_IDENTIFY, sat) -
						_allWLUPD->get_upd_value(UPDTYPE::EWL, EWL_IDENTIFY, site);
					if (fabs(ebw - round(ebw)) < 0.30 && _allWLUPD->get_upd_sigma(UPDTYPE::EWL, EWL_IDENTIFY, site) <= 0.150 &&
						_allWLUPD->get_upd_sigma(UPDTYPE::EWL, EWL_IDENTIFY, sat) <= 0.150) {
						AMB[site][sat].bewi = round(ebw);
						_AMB_RES[site][sat].ewl_fixed = true;
					}
					else {
						AMB[site][sat].bewi = 0.0;
						_AMB_RES[site][sat].ewl_fixed = false;
					}
				}
			}


			//fix wide-lane using the estimated rec. UPD
#ifdef DEBUG
			std::cout << " FIX_ID PRN  #OBS    Nw+UPD      Sigma_Nw      Nl  " << endl;
#endif

			for (auto it_sat = AMB[site].begin(); it_sat != AMB[site].end(); ++it_sat)
			{
				// add for end of ambupd file glfeng
				double dt = (it_sat->second.epo.dmjd() - t.dmjd()) * 86400;
				if (abs(dt) > 30000.0) continue;

				if (_site_ambupd_end[site] < t) continue;

				string sat = it_sat->first;

				double bw = AMB[site][sat].bw0 + AMB[site][sat].bw - _allWLUPD->get_upd_value(UPDTYPE::WL, WL_IDENTIFY, sat) -
					_allWLUPD->get_upd_value(UPDTYPE::WL, WL_IDENTIFY, site);
				double sigma = AMB[site][sat].sbw;

				// define tmp freq
				string str_idx;
				str_idx = sat.substr(0, 1) == "R" ? sat : sat.substr(0, 1);
				if (fabs(bw - round(bw)) < 0.30 && sigma < 0.130 &&
					_allWLUPD->get_upd_sigma(UPDTYPE::WL, WL_IDENTIFY, site) <= 0.150 &&
					_allWLUPD->get_upd_sigma(UPDTYPE::WL, WL_IDENTIFY, sat) <= 0.150)
				{

					_AMB_RES[site][sat].wl_fixed = true;
					_AMB_RES[site][sat].wl_res = bw - round(bw);

					//compute narrow-lane ambiguties
					AMB[site][sat].bn = -AMB[site][sat].bc / _sys_wavelen[str_idx]["NL"] +
						_sys_wavelen[str_idx]["WL"] / _sys_wavelen[str_idx]["L2"] * round(bw);
					AMB[site][sat].bwi = round(bw);

#ifdef DEBUG
					// cout .......
					std::cout << fixed << " FIX  " + site + "  " + sat << setprecision(2) << setw(15)
						<< it_sat->second.bw << setw(15) << it_sat->second.bw0 +
						AMB[site][sat].bw - _allWLUPD->get_upd_value(UPDTYPE::WL, WL_IDENTIFY, sat) <<
						setw(15) << bw << setprecision(4) << setw(10) << sigma <<
						setprecision(2) << setw(15) << it_sat->second.bn << setw(15) <<
						it_sat->second.bc / _sys_wavelen[str_idx]["NL"] << endl;
#endif

					_sat_site_amb[site][sat] = getFraction(it_sat->second.bn, -0.50);
				}

				else
				{
					_AMB_RES[site][sat].wl_fixed = false;
					_AMB_RES[site][sat].wl_res = bw - round(bw);
#ifdef DEBUG
					// cout .......
					std::cout << fixed << " REAL " + site + "  " + sat << setprecision(2) << setw(15)
						<< AMB[site][sat].bw << setw(15) << AMB[site][sat].bw0 +
						AMB[site][sat].bw - _allWLUPD->get_upd_value(UPDTYPE::WL, WL_IDENTIFY, sat) <<
						setw(15) << bw << setprecision(4) << setw(10) << sigma <<
						setprecision(2) << setw(15) << AMB[site][sat].bn << setw(15) <<
						AMB[site][sat].bc / _sys_wavelen[str_idx]["NL"] << endl;
#endif
				}
			}
			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 1, " ERROR[t_gupdlsq::_calculateAmbNL] : unknown mistake !");
			return false;
		}
	}

	bool t_gupdlsq::_connectAmb(string site, t_all_sats_OW& OW)
	{
		struct OneArc
		{
			bool is_ref;   // ref arc or not
			bool is_con;   // connect to ref arc or not
			double rnw;
			double srnw;
			int nepos;
		};
		try
		{

			vector<OneArc> ARC;
			// Get _lastSat's position in OW
			// For example, if you process all GPS sats, _lastSat is G32, its index is 32.
			// May be it is only useful in BDS UPD estimation because of GEO's poor quality.
			//int pos;
			//for (auto iter = _sats.begin(); iter != _sats.end(); iter++)
			//{
			//	if (*iter == _lastSat) {
			//		pos = distance(_sats.begin(), iter) + 1;
			//	}
			//}
#ifdef DEBUG
			std::cout << "" << endl;
			if (_updmode == UPDTYPE::EWL) std::cout << " Connecting Extra-Wide-lane Data " << endl;
			if (_updmode == UPDTYPE::EWL24) std::cout << " Connecting Extra-Wide-lane(24) Data " << endl;
			if (_updmode == UPDTYPE::EWL25) std::cout << " Connecting Extra-Wide-lane(25) Data " << endl;
			if (_updmode == UPDTYPE::WL) std::cout << " Connecting Wide-lane Data " << endl;
#endif
			auto itsat_ow = OW.begin();
			int  j, k;
			double dif, idif, sigma;
			// loop sats
			for (itsat_ow; itsat_ow != OW.end(); itsat_ow++)
			{
				ARC.clear();
				auto itepo = (itsat_ow->second).begin();
				// loop arcs and put them into vector structure
				for (itepo; itepo != itsat_ow->second.end(); itepo++)
				{
					OneArc tmp_arc;
					tmp_arc.is_ref = false;
					tmp_arc.is_con = false;
					tmp_arc.rnw = (*itepo)->rwl;
					tmp_arc.srnw = (*itepo)->srwl;
					tmp_arc.nepos = (*itepo)->end_epo - (*itepo)->beg_epo;
					ARC.push_back(tmp_arc);
				}
				// connect ambiguities to reference ambiguities (reference ambiguities may be more than one)
				if (ARC.size() == 0) continue;
				if (ARC.size() == 1)
				{
					ARC[0].is_ref = true;
				}
				else
				{
					for (j = 0; j < ARC.size() - 1; j++)
					{
						if (ARC[j].is_con) continue;
						for (k = j + 1; k < ARC.size(); k++)
						{
							if (ARC[k].is_con) continue;
							dif = ARC[k].rnw - ARC[j].rnw;
							idif = round(dif);
							dif = dif - idif;
							sigma = sqrt(pow(ARC[k].srnw, 2) + pow(ARC[j].srnw, 2));
							if (sigma < 0.10) sigma = 0.10;
							if (abs(dif) > 3 * sigma || sigma > 0.15) continue;
							ARC[j].is_ref = true;
							ARC[k].is_con = true;
							ARC[j].nepos = ARC[j].nepos + ARC[k].nepos;
							ARC[j].rnw = (ARC[j].rnw * pow(ARC[k].srnw, 2) + (ARC[k].rnw - idif) * pow(ARC[j].srnw, 2))
								/ (pow(ARC[k].srnw, 2) + pow(ARC[j].srnw, 2));
							ARC[j].srnw = ARC[k].srnw * ARC[j].srnw / sqrt(ARC[k].srnw / pow(ARC[j].srnw, 2) + ARC[j].srnw / pow(ARC[k].srnw, 2));
						}
					}
				}

				//selected the longst ref ambiguity 
				k = 0;
				for (j = 0; j < ARC.size(); j++)
				{
					if (ARC[j].is_con) continue;
#ifdef DEBUG
					std::cout << "  ..connected " << right << setprecision(2) << setw(20)
						<< ARC[j].rnw << setprecision(3) << setw(10) << ARC[j].srnw <<
						setw(6) << ARC[j].nepos << endl;
#endif
					if (k > ARC[j].nepos) continue;
					if (ARC[j].rnw == 0.0) continue;
					// if (ARC[j].nepos < 1000 && pos < 6) continue;
					// selected ambiguity of corresponding site-sat pair
					_sat_site_amb[site][itsat_ow->first] = ARC[j].rnw - round(ARC[j].rnw);
					//std::cout << "XX " << OW.size() << OW.rbegin()->first << right << setw(20) << fixed << setprecision(16) << _sat_site_amb[site][itsat_ow->first] << endl;
					k = ARC[j].nepos;
				}
			}
			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 1, " ERROR[t_gupdlsq::_connectAmb] : unknown mistake !");
			return false;
		}
	}

	bool t_gupdlsq::_setRefStation(const t_gtime& t)
	{
		int msat_max;
		int msat;

		_iref_site = "";
		_jref_site = "";

		try
		{
			//search station with max.observed satellites with valid satelite upds
			msat_max = 0;
			for (set<string>::iterator it_site = _sites.begin(); it_site != _sites.end(); ++it_site)
			{
				string site = *it_site;
				msat = 0;
				for (set<string>::iterator it_sat = _sats.begin(); it_sat != _sats.end(); ++it_sat)
				{
					string sat = *it_sat;
					if (double_eq(_sat_site_amb[site][sat], 9.90)) continue;
					if (!_cal_UPD->upd_usable(t, sat)) continue;
					msat++;
				}

				if (msat < msat_max) continue;
				msat_max = msat;
				_iref_site = site;
				if (msat == _sats.size()) break;
			}

			if (msat_max < 5) _iref_site = "";

			// search station with max. observed satellites
			// init all sat/site upd
			if (_iref_site == "")
			{
				msat_max = 0;
				for (set<string>::iterator it_site = _sites.begin(); it_site != _sites.end(); ++it_site)
				{
					string site = *it_site;
					_cal_UPD->re_init_upd(_updmode, t, site);
					msat = 0;
					for (set<string>::iterator it_sat = _sats.begin(); it_sat != _sats.end(); ++it_sat)
					{
						string sat = *it_sat;
						_cal_UPD->re_init_upd(_updmode, t, sat);
						if (double_eq(_sat_site_amb[site][sat], 9.90)) continue;
						msat++;
					}

					if (msat <= msat_max) continue;
					msat_max = msat;
					_jref_site = site;
					if (msat == _sats.size()) break;
				}
			}
			else
			{
				// init all site upds & not usable sat upds 
				for (set<string>::iterator it_site = _sites.begin(); it_site != _sites.end(); ++it_site)
				{
					string site = *it_site;
					_cal_UPD->re_init_upd(_updmode, t, site);
				}
				for (set<string>::iterator it_sat = _sats.begin(); it_sat != _sats.end(); ++it_sat)
				{
					string sat = *it_sat;
					if (_cal_UPD->upd_usable(t, sat)) continue;
					_cal_UPD->re_init_upd(_updmode, t, sat);
				}
			}

			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 1, " ERROR[t_gupdlsq::_setRefStation] : unknown mistake !");
			return false;
		}
	}

	bool t_gupdlsq::_initUPD(const t_gtime& epoch)
	{
		int iIter;
		double mean, sigma, sigx;
		t_updrec one_upd;
		list<pair<double, double>> dxx; // initial value + wgt --> mean value/sig

#ifdef DEBUG
		std::cout << " iref, jref " << setw(10) << _iref_site << "," << setw(10) << _jref_site << endl;
#endif
		// take the fractional part of the reference station as satellite upd (assume
		// receiver upd is zero.It is possible that a satellite is not tracked at the
		// station.That means no upd for these satellites.This must be solved after
		// receiver upds are determined
		try
		{
			// 设立卫星upd初值(iref缺失)
			if (_iref_site.empty())
			{
				for (auto it_sat = _sats.begin(); it_sat != _sats.end(); it_sat++)
				{
					if (double_eq(_sat_site_amb[_jref_site][*it_sat], 9.90) || _cal_UPD->upd_usable(epoch, *it_sat))
						continue;

					one_upd.value = _sat_site_amb[_jref_site][*it_sat];
					one_upd.sigma = 1e-4;
					one_upd.npoint = 0;

					_cal_UPD->add_sat_upd(_updmode, epoch, *it_sat, one_upd);
				}
			}

#ifdef DEBUG
			std::cout << "" << endl;
			for (auto it_sat = _sats.begin(); it_sat != _sats.end(); it_sat++)
			{
				std::cout << " satupd " << right << scientific << setprecision(12) << setw(20)
					<< _cal_UPD->get_upd_value(_updmode, epoch, *it_sat) << setw(20) <<
					_cal_UPD->get_upd_sigma(_updmode, epoch, *it_sat) << "   " << left << boolalpha
					<< _cal_UPD->upd_usable(epoch, *it_sat) << endl;
			}
#endif
			//iteratively
			for (iIter = 1; iIter <= 3; iIter++)
			{
#ifdef DEBUG
				std::cout << "" << endl;
				std::cout << " Writing dxx of station " << endl;
#endif

				//计算除参考站以外的测站的upd
				// estimate the upds for receivers by removing satellite upds from the ud-amb.
				// and then taking the mean.
				for (auto it_site = _sites.begin(); it_site != _sites.end(); it_site++)
				{

					if (_cal_UPD->upd_usable(epoch, *it_site)) continue;

					//dxx清零
					dxx.clear();

					for (auto it_sat = _sats.begin(); it_sat != _sats.end(); it_sat++)
					{
						if (double_eq(_sat_site_amb[*it_site][*it_sat], 9.90) || !_cal_UPD->upd_usable(epoch, *it_sat))
						{
							continue;
						}

						dxx.push_back(make_pair(_sat_site_amb[*it_site][*it_sat] - _cal_UPD->get_upd(_updmode)[epoch][*it_sat]->value, 1.0));
					}

					if (dxx.size() > 1)
					{
						//计算该测站的upd(非参考站)
						getMeanFract(dxx, mean, sigma, sigx);

						if (sigma > 999.0 || double_eq(sigma, 999.0)) continue;
						one_upd.value = mean;
						one_upd.sigma = sigx;
						one_upd.npoint = dxx.size();
						_cal_UPD->add_sat_upd(_updmode, epoch, *it_site, one_upd);

#ifdef DEBUG
						std::cout << right << *it_site + "     " << fixed << setprecision(3) << setw(8) <<
							one_upd.value << setw(8) << one_upd.sigma << setw(5) << one_upd.npoint << endl;
#endif
						//set UPD for satellites observed at this station but not yet has an intial 
						//计算该测站对应的初值仍为0的卫星
						for (auto it_sat = _sats.begin(); it_sat != _sats.end(); it_sat++)
						{
							if (double_eq(_sat_site_amb[*it_site][*it_sat], 9.90) || _cal_UPD->upd_usable(epoch, *it_sat) ||
								!_cal_UPD->upd_usable(epoch, *it_site))
							{
								continue;
							}

							one_upd.value = _sat_site_amb[*it_site][*it_sat] - _cal_UPD->get_upd(_updmode)[epoch][*it_site]->value;
							one_upd.sigma = 1e-4;
							one_upd.npoint = _cal_UPD->get_upd(_updmode)[epoch][*it_sat]->npoint;
							_cal_UPD->add_sat_upd(_updmode, epoch, *it_sat, one_upd);
#ifdef DEBUG
							std::cout << " " + *it_sat + "     " << fixed << setprecision(3) << setw(8) <<
								one_upd.value << setw(8) << one_upd.sigma << setw(5) << one_upd.npoint << endl;
#endif
						}
					}

				}//end sites upd


#ifdef DEBUG
				std::cout << "" << endl;
				std::cout << " Writing dxx of satellite " << endl;
#endif
				// check missing satellites at the reference stations
				// --------计算卫星的upd
				for (auto it_sat = _sats.begin(); it_sat != _sats.end(); it_sat++)
				{

					//解算接收机upd时临时算出的卫星upd的point都设的0
					if (_cal_UPD->upd_usable(epoch, *it_sat) && _cal_UPD->get_upd(_updmode)[epoch][*it_sat]->npoint != 0)
					{
						continue;
					}

					//dxx清零
					dxx.clear();

					for (auto it_site = _sites.begin(); it_site != _sites.end(); it_site++)
					{
						if (double_eq(_sat_site_amb[*it_site][*it_sat], 9.90) || !_cal_UPD->upd_usable(epoch, *it_site))
						{
							continue;
						}

						dxx.push_back(make_pair(_sat_site_amb[*it_site][*it_sat] - _cal_UPD->get_upd(_updmode)[epoch][*it_site]->value, 1.0));

#ifdef DEBUG
						std::cout << *it_sat + "   " + *it_site + "   " << fixed << setprecision(3) << setw(8) <<
							_sat_site_amb[*it_site][*it_sat] - _cal_UPD->get_upd(_updmode)[epoch][*it_site]->value << setw(8) <<
							_sat_site_amb[*it_site][*it_sat] << setw(8) << _cal_UPD->get_upd(_updmode)[epoch][*it_site]->value << endl;
#endif
					}

					if (dxx.size() > 1)
					{
						//计算该卫星的upd
						getMeanFract(dxx, mean, sigma, sigx);

						one_upd.value = mean;
						one_upd.sigma = sigx;
						one_upd.npoint = dxx.size();
						_cal_UPD->add_sat_upd(_updmode, epoch, *it_sat, one_upd);
#ifdef DEBUG
						//std::cout
						std::cout << " " << *it_sat << "  " << fixed << setprecision(3) << setw(10) <<
							_cal_UPD->get_upd_value(_updmode, epoch, *it_sat) << setw(10) << _cal_UPD->get_upd_sigma(_updmode, epoch, *it_sat)
							<< setw(4) << _cal_UPD->get_upd(_updmode)[epoch][*it_sat]->npoint << right << setw(6) <<
							_cal_UPD->upd_usable(epoch, *it_sat) << endl;
#endif
					}
				}//end sats upd
			}// end Iteration

#ifdef DEBUG
			std::cout << "" << endl;
			for (set<string>::iterator it_sat = _sats.begin(); it_sat != _sats.end(); ++it_sat)
			{
				string sat = *it_sat;
				std::cout << " " << sat << "  " << right << setprecision(3) << setw(10) <<
					_cal_UPD->get_upd_value(_updmode, epoch, sat) << setw(10) << _cal_UPD->get_upd_sigma(_updmode, epoch, sat)
					<< setw(4) << _cal_UPD->get_upd(_updmode)[epoch][sat]->npoint << "  " << left << boolalpha << setw(5) <<
					_cal_UPD->upd_usable(epoch, sat) << endl;
			}

			for (set<string>::iterator it_site = _sites.begin(); it_site != _sites.end(); ++it_site)
			{
				string site = *it_site;
				std::cout << " " << site << "  " << right << setprecision(3) << setw(9) <<
					_cal_UPD->get_upd_value(_updmode, epoch, site) << setw(10) << _cal_UPD->get_upd_sigma(_updmode, epoch, site)
					<< setw(4) << _cal_UPD->get_upd(_updmode)[epoch][site]->npoint << "  " << left << boolalpha << setw(5) <<
					_cal_UPD->upd_usable(epoch, site) << endl;

			}

			std::cout << "" << endl;
			std::cout << " Ambiguities After Remove Initial receiver&satelite UPDs" << endl;
			for (set<string>::iterator it_site = _sites.begin(); it_site != _sites.end(); ++it_site)
			{
				string site = *it_site;
				std::cout << " " << site;

				for (set<string>::iterator it_sat = _sats.begin(); it_sat != _sats.end(); ++it_sat)
				{
					string sat = *it_sat;

					if (double_eq(_sat_site_amb[site][sat], 9.90))
					{
						std::cout << "  .. ";
						continue;
					}
					else
					{
						double cc_double = _sat_site_amb[site][sat] - _cal_UPD->get_upd_value(_updmode, epoch, site)
							- _cal_UPD->get_upd_value(_updmode, epoch, sat);
						std::cout << " " << fixed << setprecision(1) << setw(4)
							<< cc_double;
					}
				}
				std::cout << endl;
			}
#endif
			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 1, " ERROR[t_gupdlsq::_initUPD] : unknown mistake !");
			return false;
		}
	}

	bool t_gupdlsq::_updEstimate(const t_gtime& t)
	{
		// init par num and deleted par num
		int npar = _sites.size() + _sats.size();
		int npar_delete;
		// store idx of par in NEQ
		map<string, int> par_idx;
		// store used obs num in par estimate, default is 0
		map<string, int> par_nobs;
		// store site-sat pair amb whether is fixed, default is false
		map<pair<string, string>, bool> par_fixed;
		// ambiguities num
		int ntot_save_ambiguity;   //number fixed last time
		int ntot_ambiguity;        //number fixable
		int nall_ambiguity;        //number of all ambiguities
		int usable_ambiguity;     //number with both sat. and rec. upds
		//for sigma calculate
		double ltpl;
		//index
		int i, j, iter;
		// distinct NL/WL/EWL
		double threshold = 0.25;

		try
		{
			////select a refernece satellite with max. nsite
			//int msite_max = 0;
			//int msite;		
			//string iref_sat = "";
			//for (set<string>::iterator it_sat = _sats.begin(); it_sat != _sats.end(); ++it_sat)
			//{
			//	string sat = *it_sat;
			//	msite = 0;
			//	for (set<string>::iterator it_site = _sites.begin(); it_site != _sites.end(); ++it_site)
			//	{
			//		string site = *it_site;
			//		if (!double_eq(_sat_site_amb[site][sat], 9.90)) msite++;
			//		if (msite < msite_max) continue;
			//		msite_max = msite;
			//		iref_sat = sat;
			//		if (msite == _sites.size()) break;
			//	}
			//}

			// NL/EWL/WL may only eliminate a ref value
			/*if (_updtype == UPDTYPE::NL)
			{
				updref = _cal_UPD->get_upd_value(t, iref_sat);
				for (set<string>::iterator it_sat = _sats.begin(); it_sat != _sats.end(); ++it_sat)
				{
					_cal_UPD->reset_upd_value(t, *it_sat, _cal_UPD->get_upd_value(t, *it_sat) - updref);
				}
				for (set<string>::iterator it_site = _sites.begin(); it_site != _sites.end(); ++it_site)
				{
					_cal_UPD->reset_upd_value(t, *it_site, _cal_UPD->get_upd_value(t, *it_site) + updref);
				}
			}*/

			// Combine equation
			ntot_save_ambiguity = 0;
			for (iter = 1; iter <= 10; iter++)
			{
				// init par_nobs & par_fixed
				for (set<string>::iterator it_sat = _sats.begin(); it_sat != _sats.end(); ++it_sat)
				{
					string sat = *it_sat;
					par_nobs[sat] = 0;
					for (set<string>::iterator it_site = _sites.begin(); it_site != _sites.end(); ++it_site)
					{
						string site = *it_site;
						if (it_sat == _sats.begin()) par_nobs[site] = 0;
						par_fixed[make_pair(site, sat)] = false;
					}
				}

				// prepare NEQ/W
				Matrix NEQ(npar, npar); ColumnVector W(npar);
				ColumnVector xpar(npar); ColumnVector sxpar(npar);
				NEQ = 0.0; W = 0.0; xpar = 0.0; sxpar = 0.0;
				ltpl = 0.0; ntot_ambiguity = 0; nall_ambiguity = 0; usable_ambiguity = 0;
				i = 0;
				// loop sites
				for (set<string>::iterator it_site = _sites.begin(); it_site != _sites.end(); ++it_site)
				{
					i++;
					string site = *it_site;
					par_idx[site] = i + _sats.size();
					// loop sats
					j = 0;
					for (set<string>::iterator it_sat = _sats.begin(); it_sat != _sats.end(); ++it_sat)
					{
						j++;
						string sat = *it_sat;
						// first sats, then sites
						if (i == 1) par_idx[sat] = j;
						// cycle site-sat pair without data
						if (double_eq(_sat_site_amb[site][sat], 9.90)) continue;
						nall_ambiguity++;
						// cycle site-sat pair with unvalid upds
						if (!(_cal_UPD->upd_usable(t, site) && _cal_UPD->upd_usable(t, sat))) continue;
						usable_ambiguity++;
						// cycle site-sat pair whose fractional part not in threshold
						double xamb;
						int iamb;
						xamb = _sat_site_amb[site][sat] - _cal_UPD->get_upd_value(_updmode, t, sat) - _cal_UPD->get_upd_value(_updmode, t, site);
						iamb = round(xamb);
						if (fabs(xamb - iamb) > threshold) continue;
						// record site-sat pair used 
						par_fixed[make_pair(site, sat)] = true;
						// NEQ prepare
						NEQ(j, j) += 1.0;
						NEQ(j, i + _sats.size()) += 1.0;
						NEQ(i + _sats.size(), i + _sats.size()) += 1.0;
						// W prepare
						W(j) = W(j) - xamb + iamb;
						W(i + _sats.size()) = W(i + _sats.size()) - xamb + iamb;
						// record num of obs of obj (sat/site)
						par_nobs[sat]++;
						par_nobs[site]++;
						ntot_ambiguity++;
						// later become vtpv
						ltpl += (xamb - iamb) * (xamb - iamb);
					}
				}

				//select a refernece satellite with max. nsite
				//avoid select a ref sat without valid value ^_^
				int msite_max = 0;
				string iref_sat = "";
				for (set<string>::iterator it_sat = _sats.begin(); it_sat != _sats.end(); ++it_sat)
				{
					string sat = *it_sat;
					if (par_nobs[sat] > msite_max)
					{
						iref_sat = sat;
						msite_max = par_nobs[sat];
					}
				}
				if (iref_sat == "")
				{
					cerr << " Warning : Not enough sats in iter " << iter << " during the epoch : " << t.str_ymdhms() << endl;
					break;
				}

#ifdef DEBUG
				std::cout << "" << endl;
				std::cout << " niter , nfix_last, nfix " << right << setw(6) << iter <<
					setw(6) << ntot_save_ambiguity << setw(6) << ntot_ambiguity << endl;
				if (ntot_save_ambiguity >= ntot_ambiguity)
				{
					std::cout << " # all amb, # with init # fixable " << right << setw(6) << nall_ambiguity <<
						setw(6) << usable_ambiguity << setw(6) << ntot_ambiguity << endl;
				}
#endif

				// reference receiver upd, add strong constraint to avoid rank loss
				// consider obj without obs
				NEQ(par_idx[iref_sat], par_idx[iref_sat]) += 1e5;
				npar_delete = 0;
				for (auto it = _parlist.begin(); it != _parlist.end(); it++)
				{
					if (par_nobs[*it] == 0)
					{
						NEQ(par_idx[*it], par_idx[*it]) = 1e3;
						npar_delete++;
					}
					for (j = 1; j <= par_idx[*it] - 1; j++)
					{
						NEQ(par_idx[*it], j) = NEQ(j, par_idx[*it]);
					}
				}

				//if (double_eq(NEQ.determinant(), 0.0))
				//{
				//	cerr << " Warning : NEQ Matrix singularity" << " during the epoch : " << t.str_ymdhms() << endl;
				//	return false;
				//}
				// inv
				Matrix invNEQ;
				invNEQ << NEQ.i();

				// get solve
				for (i = 1; i <= npar; i++)
				{
					xpar(i) = 0.0;
					for (j = 1; j <= npar; j++)
					{
						xpar(i) += invNEQ(i, j) * W(j);
					}
					ltpl -= xpar(i) * W(i);  // vtpv
				}

				//cal sigma
				double sigma0 = 0.0;

				if ((ntot_ambiguity - (npar - npar_delete)*1.0) > 0.0)
				{
					sigma0 = sqrt(ltpl / (ntot_ambiguity - (npar - npar_delete)*1.0));
				}

				// add solve value to init value
				for (auto it = _parlist.begin(); it != _parlist.end(); it++)
				{
					if (par_nobs[*it] != 0)
					{
						sxpar(par_idx[*it]) = sigma0 * sqrt(invNEQ(par_idx[*it], par_idx[*it]));
					}
					if (*it != iref_sat && par_nobs[*it] == 0) continue;

#ifdef DEBUG
					std::cout << right << setw(5) << " " + *it << setprecision(3) << setw(10) <<
						_cal_UPD->get_upd_value(_updmode, t, *it) << setw(10) << xpar(par_idx[*it]) <<
						setw(10) << sxpar(par_idx[*it]) << setw(4) << par_nobs[*it] << endl;
#endif

					double tmp = _cal_UPD->get_upd_value(_updmode, t, *it) - xpar(par_idx[*it]);
					_cal_UPD->reset_upd(_updmode, t, *it, tmp, sxpar(par_idx[*it]), par_nobs[*it]);
				}

				//!! now fixed to 3 iteration
				if (ntot_save_ambiguity >= ntot_ambiguity) break;
				ntot_save_ambiguity = ntot_ambiguity;
			}


			// AMB_RES
			for (set<string>::iterator it_site = _sites.begin(); it_site != _sites.end(); ++it_site)
			{
				string site = *it_site;

				for (set<string>::iterator it_sat = _sats.begin(); it_sat != _sats.end(); ++it_sat)
				{
					string sat = *it_sat;
					if (double_eq(_sat_site_amb[site][sat], 9.90)) continue;
					if (!par_fixed[make_pair(site, sat)])
					{
						_AMB_RES[site][sat].nl_fixed = false;
						_AMB_RES[site][sat].nl_res = _sat_site_amb[site][sat] - _cal_UPD->get_upd_value(_updmode, t, sat)
							- _cal_UPD->get_upd_value(_updmode, t, site) - round(_sat_site_amb[site][sat]
								- _cal_UPD->get_upd_value(_updmode, t, sat) - _cal_UPD->get_upd_value(_updmode, t, site));

					}
					else
					{
						_AMB_RES[site][sat].nl_fixed = true;
						_AMB_RES[site][sat].nl_res = _sat_site_amb[site][sat] - _cal_UPD->get_upd_value(_updmode, t, sat)
							- _cal_UPD->get_upd_value(_updmode, t, site) - round(_sat_site_amb[site][sat]
								- _cal_UPD->get_upd_value(_updmode, t, sat) - _cal_UPD->get_upd_value(_updmode, t, site));
					}
				}
			}
#if DEBUG
			// Output to screen
			std::cout << "" << endl;
			std::cout << " Ambiguities After Remove estimated receiver&satelite UPDs (0.1 cycle)" << endl;
			for (set<string>::iterator it_site = _sites.begin(); it_site != _sites.end(); ++it_site)
			{
				string site = *it_site;
				std::cout << " " << site << "   ";
				for (set<string>::iterator it_sat = _sats.begin(); it_sat != _sats.end(); ++it_sat)
				{
					string sat = *it_sat;
					if (double_eq(_sat_site_amb[site][sat], 9.90))
					{
						std::cout << " .. ";
						continue;
					}
					else
					{
						double cc_double = (_sat_site_amb[site][sat] - _cal_UPD->get_upd_value(_updmode, t, site)
							- _cal_UPD->get_upd_value(_updmode, t, sat)) * 10.0;
						int cc_tmp = int(round(cc_double)) % 10;
						if (!par_fixed[make_pair(site, sat)]) {
							std::cout << right << setw(3) << cc_tmp << "*";
						}
						else
						{
							std::cout << right << setw(3) << cc_tmp << " ";
						}
					}
				}
				std::cout << endl;
			}
#endif
			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 1, " ERROR[t_gupdlsq::updEstimate] : unknown mistake !");
			return false;
		}
	}

	bool t_gupdlsq::_writeAmbResidual(const t_gtime& t)
	{
		ostringstream os;
		string wl_fix, nl_fix;

		//write Time
		os << " EPOCH-TIME" << setw(8) << t.mjd() << fixed << setprecision(1)
			<< setw(10) << t.sod() << endl;

		//write Data
		auto it_site = _AMB_RES.begin();
		for (it_site; it_site != _AMB_RES.end(); it_site++)
		{
			auto it_sat = it_site->second.begin();
			for (it_sat; it_sat != it_site->second.end(); it_sat++)
			{

				if (it_sat->second.wl_fixed)
				{
					wl_fix = "yes";
				}
				else
				{
					wl_fix = "no";
				}
				if (it_sat->second.nl_fixed)
				{
					nl_fix = "yes";
				}
				else
				{
					nl_fix = "no";
				}

				os << "  " << fixed << it_site->first << "  " << it_sat->first << right << setw(10) << setprecision(3) <<
					it_sat->second.wl_res << " " << left << setw(3) << wl_fix << " " << right << setw(10) << setprecision(3)
					<< it_sat->second.nl_res << " " << left << setw(3) << nl_fix  << endl;
			}
		}


		// Print amb res results
		if (_fambres) {
			_fambres->write(os.str().c_str(), os.str().size());
			_fambres->flush();
		}
		else {
			write_log_info(_glog, 1, " Warning[t_gupdlsq::_writeAmbResidual] : have no output file!");
		}
		return true;
	}


	void t_gupdlsq::_readOfile(set<string>::iterator it_site,t_gallobs& obs) {
		string rec = *it_site;
		//start reading o files one by one
		auto recs = dynamic_cast<t_gsetgen*>(_set)->recs();
		t_gdata* gdata = nullptr;

		t_gtime beg_reading = t_gtime::current_time(t_gtime::GPS);
		multimap<IFMT, string> inp = dynamic_cast<t_gsetinp*>(_set)->inputs_all();
		multimap<IFMT, string>::const_iterator itINP = inp.begin();
		for (size_t i = 0; i < inp.size() && itINP != inp.end(); ++i, ++itINP)
		{
			// Get the file format/path, which will be used in decoder
			IFMT   ifmt(itINP->first);
			string path(itINP->second);
			string id("ID" + int2str(i));

			t_gio* gio = nullptr;
			t_gcoder* gcoder = nullptr;

			// For different file format, we prepare different data container and decoder for them.
			if (ifmt == RINEXO_INP) { gdata = &obs; gcoder = new t_rinexo(_set, "", 5000); }
			else
			{
				continue;
			}

			// Check the file path
			if (path.substr(0, 7) != "file://")
			{
				_glog->comment(0, "main", "path is not file (skipped)!");
				continue;
			}

			// check the obs file valid
			auto pos = path.find_last_of("o");
			string name = path.substr(pos - 11, 4);
			transform(name.begin(), name.end(), name.begin(), ::toupper);
			if (name != rec) continue;
			if (recs.count(name) <= 0)
			{
				_glog->comment(0, "main", "no need to read this file (skipped)!");
				continue;
			}

			// runepoch for the time costed each epoch (i guess)
			t_gtime runepoch(t_gtime::GPS);
			t_gtime lstepoch(t_gtime::GPS);
			// READ DATA FROM FILE
			if (gcoder) {

				// prepare gio for the file
				gio = new t_gfile;
				gio->glog(_glog);
				gio->path(path);

				// cerr << LV + " Now read the file : " + path << endl;
				//glog.comment(t_glog::LOG_LV::LOG_INFO, "Main", "Main", " Now read the file : " + path);

				// Put the file into gcoder
				gcoder->clear();
				gcoder->path(path);
				gcoder->glog(_glog);
				// Put the data container into gcoder
				gcoder->add_data(id, gdata);
				//gcoder->add_data("OBJ", &obj);
				// Put the gcoder into the gio
				// Note, gcoder contain the gdata and gio contain the gcoder
				gio->coder(gcoder);

				runepoch = t_gtime::current_time(t_gtime::GPS);
				// Read the data from file here
				gio->run_read();
				lstepoch = t_gtime::current_time(t_gtime::GPS);
				// Write the information of reading process to log file
				_glog->comment(t_glog::LOG_LV::LOG_INFO, "main", "main", "READ: " + path + " time: "
					+ dbl2str(lstepoch.diff(runepoch)) + " sec");
				delete gio; gio = NULL;
				delete gcoder; gcoder = NULL;
			}

		}
		
	}


	bool t_gupdlsq::_writeCarrierRange(const t_gtime &t_beg, const double &intv)
	{
		// set output ambflag path
		string strambflagpath = dynamic_cast<t_gsetout*>(_set)->outputs("ambflag_dir");
		if (strambflagpath.empty())   strambflagpath = "file://log_fix";
		if (ACCESS(strambflagpath.substr(7).c_str(), 0) != 0)
		{
			MKDIR(strambflagpath.substr(7).c_str());
		}
		// loop all sites
		for (const auto &site : _sites) {
			string lower_site(site);
			transform(site.begin(), site.end(), lower_site.begin(), ::tolower);
			// check if there is ambflag file
			if (_allambflag->getAllAmbFlag().find(lower_site) == _allambflag->getAllAmbFlag().end())
				continue;
			auto ambflag_site = _allambflag->getOneAmbFlag(lower_site);
			auto amb_site = _allfixedamb->getOneSiteAmbUpd(site);
			for (const auto &ambflag : ambflag_site.getAmbFlagData()) {
				if (_sats.find(ambflag.first) == _sats.end()) continue;
				// loop over each ambflag record of one satellite
				for (const auto &iter : ambflag.second) {
					// skip DEL and BAD records
					if (iter->iflag != "1" && iter->iflag != "0" 
						|| iter->identify != "AMB" && iter->identify != "IAM") continue;
					t_gtime beg = ambflag_site.epoch2time(iter->beg_epo);
					t_gtime end = ambflag_site.epoch2time(iter->end_epo);
					if (beg > end) continue;
					unsigned nepo = (end - beg) / intv + 1;
					// get all fixed ambiguities of this arc
					beg = beg - intv;
					end = end + intv;
					t_satepo_ambupd satepo_amb = amb_site.getSatepoAmbUpd(beg, end, ambflag.first);
					// at least 50% ambiguity in a arc should be fixed
					if (satepo_amb.size() < round(nepo*0.5)) continue;
					// make sure that all ambiguity fixed to the same integer value
					double n1_first = satepo_amb.begin()->second->amb_if;
					double n2_first = satepo_amb.begin()->second->amb_wl;
					double n1_last = satepo_amb.rbegin()->second->amb_if;
					double n2_last = satepo_amb.rbegin()->second->amb_wl;
					if (fabs(n1_first - n1_last) < 0.05 && fabs(n2_first - n2_last) < 0.05) {
						iter->identify = "IAM";
						iter->iflag = "1";
						iter->C1 = n1_last;
						iter->C2 = n2_last;
					}
				}
			}
			string logsuffix = ".log";
			stringstream str;
			str << strambflagpath << PATH_SEPARATOR << lower_site << setw(3) << setfill('0') << t_beg.doy() << "0."
				<< setw(2) << setfill('0') << t_beg.yr() << "o" + logsuffix;
			t_gio* gout = new t_gfile;
			//set I/O
			gout->path(str.str()); gout->glog(_glog);
			//set coder
			t_gcoder* gcoder = new t_ambflag(nullptr, "", 4096);
			gcoder->clear(); gcoder->glog(_glog); gcoder->path(str.str());
			gcoder->add_data("IDX", _allambflag);
			gout->coder(gcoder);
			// write
			gout->run_write();
			delete gout; gout = nullptr;
			delete gcoder; gcoder = nullptr;
		}

		// for ambflag13
		if (_allambflag13) {
			for (const auto& site : _sites) {
				string lower_site(site);
				transform(site.begin(), site.end(), lower_site.begin(), ::tolower);
				// check if there is ambflag file
				if (_allambflag13->getAllAmbFlag().find(lower_site) == _allambflag13->getAllAmbFlag().end())
					continue;
				auto ambflag_site = _allambflag13->getOneAmbFlag(lower_site);
				auto amb_site = _allfixedamb->getOneSiteAmbUpd(site);
				for (const auto& ambflag : ambflag_site.getAmbFlagData()) {
					if (_sats.find(ambflag.first) == _sats.end()) continue;
					// loop over each ambflag record of one satellite
					for (const auto& iter : ambflag.second) {
						// skip DEL and BAD records
						if (iter->iflag != "1" && iter->iflag != "0"
							|| iter->identify != "AMB" && iter->identify != "IAM") continue;
						t_gtime beg = ambflag_site.epoch2time(iter->beg_epo);
						t_gtime end = ambflag_site.epoch2time(iter->end_epo);
						if (beg > end) continue;
						unsigned nepo = (end - beg) / intv + 1;
						// get all fixed ambiguities of this arc
						beg = beg - intv;
						end = end + intv;
						t_satepo_ambupd satepo_amb = amb_site.getSatepoAmbUpd(beg, end, ambflag.first);
						// at least 50% ambiguity in a arc should be fixed
						if (satepo_amb.size() < round(nepo * 0.5)) continue;
						// make sure that all ambiguity fixed to the same integer value
						double n1_first = satepo_amb.begin()->second->amb_if;
						double n3_first = satepo_amb.begin()->second->amb_l3;
						double n1_last = satepo_amb.rbegin()->second->amb_if;
						double n3_last = satepo_amb.rbegin()->second->amb_l3;
						if (fabs(n1_first - n1_last) < 0.05 && fabs(n3_first - n3_last) < 0.05) {
							iter->identify = "IAM";
							iter->iflag = "1";
							iter->C1 = n1_last;
							iter->C2 = n3_last;
						}
					}
				}
				string logsuffix = ".log13";
				stringstream str;
				str << strambflagpath << PATH_SEPARATOR << lower_site << setw(3) << setfill('0') << t_beg.doy() << "0."
					<< setw(2) << setfill('0') << t_beg.yr() << "o" + logsuffix;
				t_gio* gout = new t_gfile;
				//set I/O
				gout->path(str.str()); gout->glog(_glog);
				//set coder
				t_gcoder* gcoder = new t_ambflag(nullptr, "", 4096);
				gcoder->clear(); gcoder->glog(_glog); gcoder->path(str.str());
				gcoder->add_data("IDX", _allambflag13);
				gout->coder(gcoder);
				// write
				gout->run_write();
				delete gout; gout = nullptr;
				delete gcoder; gcoder = nullptr;
			}
		}
		return true;
	}

}
