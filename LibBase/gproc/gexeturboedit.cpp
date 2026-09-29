/**
*
* @verbatim
	History
	 -1.0 Guolong Feng  2020-03-20 created
  @endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file		    gexeturboedit.cpp
* @brief		main control class for turboedit
* @author       Guolong Feng, Wuhan University
* @date		    2020-03-20
*
*/

#include "gproc/gexeturboedit.h"

namespace great
{

	t_gexeturboedit::t_gexeturboedit(t_gsetbase* set, t_glog* glog, t_gallobs* gobs, t_gallnav* gnav):
		_set(set),
		_log(glog),
		_allobs(gobs),
		_allnav(gnav),
		_allobj(nullptr),
		_check_sf(false),
		_check_pc(true),
		_amb_output(true),
		_simulation(false),
		_check_short(true),
		_allambflag(nullptr),
		_site_kinematic(false),
		_check_statistics(true),
		_outliers_proc(set, glog)
	{
		// common intialization
		_initCommonProcess();

		// init post-processing / real-time-processing
		if (_liteMode || _realtime) _initRealTimeProcess();
		else                        _initPostProcess();        

		// output control infomation
		_outputControlInfo();
	
	}

	t_gexeturboedit::t_gexeturboedit(t_gsetbase* set, t_glog* glog, t_gallproc* gdata) :
		_set(set),
		_log(glog),
		_check_sf(false),
		_check_pc(true),
		_amb_output(true),
		_simulation(false),
		_check_short(true),
		_allambflag(nullptr),
		_site_kinematic(false),
		_check_statistics(true),
		_outliers_proc(set, glog)
	{
		// get data 
		_allobs = dynamic_cast<t_gallobs*>((*gdata)[t_gdata::ALLOBS]);
		_allnav = dynamic_cast<t_gallnav*>((*gdata)[t_gdata::GRP_EPHEM]);
		_allobj = dynamic_cast<t_gallobj*>((*gdata)[t_gdata::ALLOBJ]);

		// common intialization
		_initCommonProcess();

		// init post-processing / real-time-processing
		if (_liteMode || _realtime) _initRealTimeProcess();
		else                        _initPostProcess();

		// output control infomation
		_outputControlInfo();

	}

	t_gexeturboedit::~t_gexeturboedit()
	{
		if (_allambflag)
		{
			delete _allambflag;
			_allambflag = nullptr;
		}

		for (auto it : _rec_list)
		{
			if (_sppflt[it]) { delete _sppflt[it]; _sppflt[it] = nullptr; }
		}

		if (_debug_turbo) 
		{ 
			if (_debug_turbo->is_open()) _debug_turbo->close();
			delete _debug_turbo; 
			_debug_turbo = nullptr;
		}

	}

	bool t_gexeturboedit::ProcessBatch(string site, const t_gtime& begT, const t_gtime& endT, double sampling)
	{
		_crt_rec = site;
		_interval = sampling;

		// Check Time Duration
		if (begT > endT)//lvhb begT>endT
		{
			_log->comment(1, _crt_rec, "Error in time settings : beg time >= end time !");
			_log->comment(1, _crt_rec, "The beg time is : " + begT.str_ymdhms());
			_log->comment(1, _crt_rec, "The end time is : " + endT.str_ymdhms());
			return false;
		}

		// Real_time Turboedit
		if (_realtime || _liteMode) 
		{
			if(!_realtimeProcess(begT, endT)) return false;
		}
		else
		{
			// =================================================
			if (!_amb_output) return true;   //用于不覆盖原有log文件
			// =================================================

			if (_rec_list.find(_crt_rec) == _rec_list.end()) return false;

			if (_allambflag)
			{
				set<string> sites = _allambflag->getSiteList();
				if (sites.find(_crt_rec) != sites.end()) return false;
			}
			
			// Check simulation and kinematic
			if (_simulation && _site_kinematic)
			{
				_log->comment(1, _crt_rec, "Not support kinematic solution for simulation data! ");
				return false;
			}

			_epo_num = floor((_end_time - _beg_time) / sampling) + 1;

			// loop band combination
			for (unsigned int i = 1; i <= _post_all_band.size(); i++)
			{
				_post_crt_band = _post_all_band[i];
				_cycleslip[_crt_rec].push_back(shared_ptr<t_gcycleslip>(new t_gturboedit()));
				
				if (!_postProcess()) return false;
			}
		}
		return true;
	}

	void t_gexeturboedit::setGLOFreqID(t_gallsatinfo* gsatpars)
	{
		if (gsatpars)
		{
			t_gsatparam sat_param;
			for (auto sat : _sat_list)
			{
				if (sat.substr(0, 1) != "R") continue;
				if (gsatpars->getSatinfoUsed(sat, _beg_time, _end_time, sat_param))
				{
					_glonass_freq_id[sat] = sat_param.fid();
				}
			}
		}
	}

	void t_gexeturboedit::set_site_crd(string site, t_gtriple crd)
	{
		map<string, t_gtriple>::iterator it;
		it = _rec_crds.find(site);
		if (it == _rec_crds.end())
		{
			_rec_crds[site] = crd;
		}
	}

	void t_gexeturboedit::setSingleEpoData(vector<t_gsatdata>* input_data)
	{
		if (!input_data) return;

		if (_satdata_npp) _satdata_npp = nullptr;
		_satdata_npp = input_data;

		int nsize = _satdata_npp->size();

		_inputEpoData.clear();
		_map_nppdata_idx.clear();
		for (unsigned int i = 0; i < nsize; i++)
		{
			shared_ptr<t_gobsgnss> gobs(new t_gsatdata((*_satdata_npp)[i]));
			_inputEpoData.push_back(gobs);
			_map_nppdata_idx[gobs->site()][gobs->sat()] = i;
		}
		return;
	}

	bool t_gexeturboedit::checkSlipWithDoppler(string site, const t_gtime& begT, const t_gtime& endT, double sampling)
	{
		_crt_rec = site;
		_interval = sampling;

		// Check Time Duration
		if (begT > endT)//lvhb begT>endT
		{
			_log->comment(1, _crt_rec, "Error in time settings : beg time >= end time !");
			_log->comment(1, _crt_rec, "The beg time is : " + begT.str_ymdhms());
			_log->comment(1, _crt_rec, "The end time is : " + endT.str_ymdhms());
			return false;
		}

		_crt_time = begT;
		vector<shared_ptr<t_gobsgnss>>       epoData;
		while (_crt_time <= endT)
		{
			epoData = _allobs->obs_pt(_crt_rec, _crt_time);

			if (epoData.size() == 0)
			{
				_crt_time = _crt_time + _interval;
				continue;
			}

			for (auto &oneObs : epoData)
			{
				string sat = oneObs->sat();
				GSYS gsys = oneObs->gsys();

				// set glonass Frequency ID
				if (sat.substr(0, 1) == "R" && oneObs->channel() >= DEF_CHANNEL)
				{
					oneObs->channel(_glonass_freq_id[sat]);
				}

				// Not check in the first epoch
				if (_beg_time == _crt_time)
				{
					_epoDataPre_npp[_crt_rec][sat] = oneObs;
					continue;
				}

				if (_epoDataPre_npp[_crt_rec].find(sat) != _epoDataPre_npp[_crt_rec].end())
				{
					// Range / Doppler outliers Detection
					if (_frequency >= 2) _outliers_proc.flagRangeOutliers(_epoDataPre_npp[_crt_rec][sat], oneObs, _interval);

					if (_frequency == 1 || _observ == RAW_MIX) _checkSlipWithDoppler(_epoDataPre_npp[_crt_rec][sat], oneObs);
				}

				_epoDataPre_npp[_crt_rec][sat] = oneObs;
			}
			// Next
			_crt_time = _crt_time + _interval;

		}
		return true;
	}


	void  t_gexeturboedit::_outputControlInfo()
	{
		stringstream str2;
		if (_log) {
			_log->comment(0, "NOTE", "***********************************************");
			_log->comment(0, "NOTE", "+++  Control Parameters from XML file  +++");

			// Realtime
			_log->comment(0, "NOTE", "RealtimeMode :  " + bl2str(_realtime));
			_log->comment(0, "NOTE", "LiteMode     :  " + bl2str(_liteMode));
		}
		if (!_realtime && !_liteMode)
		{
			// SYS & Bands
			str2.str("");
			for (auto sys : _sys_list)
			{
				GSYS gsys = t_gsys::str2gsys(sys);
				vector<GOBSBAND> v_band = dynamic_cast<t_gsetgnss*>(_set)->band(gsys);

				str2 << sys + "->";
				if (_simulation && gsys == GSYS::LEO) continue;
				for (unsigned int i = 1; i <= _post_all_band.size(); i++)
				{
					if (_post_all_band[i].find(gsys) == _post_all_band[i].end()) continue;
					str2 << "BAND_" + gobsband2str(_post_all_band[i][gsys].first) +
						"+BAND_" + gobsband2str(_post_all_band[i][gsys].second) + "  ";
				}
				str2 << "  ";
			}
			if (_log) {
				_log->comment(0, "NOTE", "System&Bands :  " + str2.str());

				// Kin/Simu/Eph/CutoffElevation
				_log->comment(0, "NOTE", "Kinematic    :  " + bl2str(_site_kinematic));
				_log->comment(0, "NOTE", "Simulation   :  " + bl2str(_simulation));
				_log->comment(0, "NOTE", "Use Ephemeris:  " + bl2str(_ephemeris));
			}
			if (_ephemeris)
			{
				str2.str("  ");
				vector<string> inp = dynamic_cast<t_gsetinp*>(_set)->inputs("rinexn");
				for (auto iter : inp)  str2 << iter << "  ";
				if (_log)_log->comment(0, "NOTE", "Ephemeris    :  " + str2.str());
				str2.str(""); str2 << fixed << left << setw(7) << setprecision(2) << _min_elev;
				if (_log)_log->comment(0, "NOTE", "ElevationMask:  " + str2.str() + " degrees ");
			}
			// PC
			if (_log)_log->comment(0, "NOTE", "Check PC     :  " + bl2str(_check_pc));
			if (_check_pc)
			{
				str2.str(""); str2 << fixed << left << setw(7) << setprecision(2) << _pc_limit;
				if (_log)_log->comment(0, "NOTE", "PCLimitation :  " + str2.str() + " meters");
			}
			// LX-CX
			if (_log)_log->comment(0, "NOTE", "Check LX-CX  :  " + bl2str(_check_sf));
			if (_check_sf)
			{
				str2.str(""); str2 << fixed << left << setw(7) << setprecision(2) << _sf_limit;
				if (_log)_log->comment(0, "NOTE", "LXLimitation :  " + str2.str() + " meters");
			}
		}
		// MW
		if (_log)_log->comment(0, "NOTE", "TurboEdit(MW):  " + bl2str(_check_mw));
		if (_check_mw)
		{
			str2.str(""); str2 << fixed << left << setw(7) << setprecision(2) << _mw_limit;
			if (_log)_log->comment(0, "NOTE", "MWLimitation :  " + str2.str() + " cycles");
		}
		// GF
		if (_log)_log->comment(0, "NOTE", "TurboEdit(GF):  " + bl2str(_check_gf));
		if (_check_gf)
		{
			if (!_realtime && !_liteMode)
			{
				str2.str(""); str2 << fixed << left << setw(7) << setprecision(2) << _gf_limit;
				if (_log)_log->comment(0, "NOTE", "GFLimitation :  " + str2.str() + " cycles");
				str2.str(""); str2 << fixed << left << setw(7) << setprecision(2) << _gf_rms_limit;
				if (_log)_log->comment(0, "NOTE", "GFRmslimit   :  " + str2.str() + " cycles");
			}
			else
			{
				str2.str(""); str2 << fixed << left << setw(7) << setprecision(2) << _gf_limit;
				if (_log)_log->comment(0, "NOTE", "GFLimitation :  " + str2.str() + " meters");
			}
		}

		// check gap
		if (_log)_log->comment(0, "NOTE", "Check Gap    :  " + bl2str(_check_gap));
		if (_check_gap)
		{
			str2.str(""); str2 << fixed << left << setw(7) << _gap_limit;
			if (_log)_log->comment(0, "NOTE", "GapLimit     :  " + str2.str() + " epochs ");
		}

		// check short
		if (!_realtime && !_liteMode)
		{
			if (_log)_log->comment(0, "NOTE", "Check Short  :  " + bl2str(_check_short));
			if (_check_short)
			{
				str2.str(""); str2 << fixed << left << setw(7) << _short_limit;
				if (_log)_log->comment(0, "NOTE", "ShortLimit   :  " + str2.str() + " epochs ");
			}
			// check statistics
			if (_check_statistics)
			{
				str2.str(""); str2 << fixed << left << setw(7) << setprecision(2) << _max_mean_namb;
				if (_log)_log->comment(0, "NOTE", "MaxMeanAmbs  :  " + str2.str());
				str2.str(""); str2 << fixed << left << setw(7) << setprecision(2) << _min_percent;
				if (_log)_log->comment(0, "NOTE", "MinPercent>4 :  " + str2.str() + " %");
				str2.str(""); str2 << fixed << left << setw(7) << setprecision(2) << _min_mean_nprn;
				if (_log)_log->comment(0, "NOTE", "MinMeanPrns  :  " + str2.str());
			}
		}

	}

	void t_gexeturboedit::_initCommonProcess()
	{
		// Global Settings
		_beg_time = dynamic_cast<t_gsetgen*>(_set)->beg();
		_end_time = dynamic_cast<t_gsetgen*>(_set)->end();
		_sys_list = dynamic_cast<t_gsetgen*>(_set)->sys();        // processing systems
		_rec_list = dynamic_cast<t_gsetgen*>(_set)->recs();       // processing receivers
		for (auto it : _sys_list)                                 // processing satellites
		{
			GSYS gsys = t_gsys::str2gsys(it);
			std::set<string> sat_temp = dynamic_cast<t_gsetgnss*>(_set)->sat(gsys);
			_sat_list.insert(sat_temp.begin(), sat_temp.end());  // insert satellites
		}
		// exclude satellites mainly For simulation
		// real data can control processed sats by module <gps>/<bds>/<gal>/..., also useful for real data
		std::set<string> sat_rm = dynamic_cast<t_gsetgen*>(_set)->sat_rm();
		for (auto iter : sat_rm) { _sat_list.erase(iter); }

		// TurboEdit Settings
		_realtime = dynamic_cast<t_gsetproc*>(_set)->realtime();
		_min_elev = dynamic_cast<t_gsetproc*>(_set)->minimum_elev();
		_liteMode = dynamic_cast<t_gsetturboedit*>(_set)->liteMode();
		_check_mw = dynamic_cast<t_gsetturboedit*>(_set)->checkMW(_mw_limit);                // cycle for post/real_time
		_check_gf = dynamic_cast<t_gsetturboedit*>(_set)->checkGF(_gf_limit, _gf_rms_limit); // meter for real_time/cycle for post
		_check_gap = dynamic_cast<t_gsetturboedit*>(_set)->checkGap(_gap_limit);
		_ephemeris = dynamic_cast<t_gsetturboedit*>(_set)->isEphemeris();    // 适用于年份较老的动态数据，提供GLONASS_freq_ID
		if (!_allnav) _ephemeris = false;

		if (!_check_gap)  _gap_limit = 999999;

		// SPPFLT initialization
		for (auto it : _rec_list)
		{
			if (_sppflt[it]) { delete _sppflt[it]; _sppflt[it] = nullptr; }
			_sppflt[it] = new t_gsppflt(it, _set);
			_sppflt[it]->setDAT(_allobs, _allnav);
			_sppflt[it]->setOBJ(_allobj);
		}

		_observ = dynamic_cast<t_gsetproc*>(_set)->obs_combin();
	}

	void t_gexeturboedit::_initRealTimeProcess()
	{
		// 平滑窗口
		_smooth_window = dynamic_cast<t_gsetturboedit*>(_set)->smoothWindows();
		// GF设置 参考 张小红《构建阈值模型改善TurboEdit实时周跳探测》
		// 对于每一组数据，GF组合历元差与高度角成反比，主要原因是卫星高度角低时，多路径效应和观测噪声的影响较大，因此，需要在上述初始阈
		// 值中引入一个与卫星高度角有关的加权因子，目前在前面乘以一个比较大的因子来替代
		// MW 采用自适应

		double sample = dynamic_cast<t_gsetgen*>(_set)->sampling();
		if (sample <= 1)
		{
			if (double_eq(_gf_limit, 0.0)) _gf_limit = 0.05;                     ///< meters
			if (double_eq(_mw_limit, 0.0)) _mw_limit = 3.0;                      ///< cycles
		}
		else if (sample <= 20)
		{
			if (double_eq(_gf_limit, 0.0)) _gf_limit = 0.1 / 20.0 * sample + 0.05;   ///< meters
			if (double_eq(_mw_limit, 0.0)) _mw_limit = 2.5 / 20.0 * sample + 2.5;    ///< cycles
		}
		else if (sample <= 60)
		{
			if (double_eq(_gf_limit, 0.0)) _gf_limit = 0.15;                    ///< meters
			if (double_eq(_mw_limit, 0.0)) _mw_limit = 4.0;                     ///< cycles
		}
		else if (sample <= 100)
		{
			if (double_eq(_gf_limit, 0.0)) _gf_limit = 0.25;                    ///< meters
			if (double_eq(_mw_limit, 0.0)) _mw_limit = 5.0;                     ///< cycles
		}
		else
		{
			if (double_eq(_gf_limit, 0.0)) _gf_limit = 0.3;                   ///< meters
			if (double_eq(_mw_limit, 0.0)) _mw_limit = 6.0;                   ///< cycles
		}

		_frequency = dynamic_cast<t_gsetproc*>(_set)->frequency();

		if (_ephemeris)
		{
			_glonass_freq_id = _allnav->glo_freq_num();
		}

		if (1) // lvhb shut down the outputting file in 202011
		{
			if (_debug_turbo)
			{
				if (_debug_turbo->is_open()) _debug_turbo->close();
				delete _debug_turbo;
				_debug_turbo = nullptr;
			}
		}
		else {
			_debug_turbo = new t_giof;
			_debug_turbo->tsys(t_gtime::GPS);
			_debug_turbo->mask("debug_turbo.log");
			_debug_turbo->append(dynamic_cast<t_gsetout*>(_set)->append());
		}

	}

	bool t_gexeturboedit::_realtimeProcess(const t_gtime& begT, const t_gtime& endT)
	{
		_crt_time = begT;
		vector<shared_ptr<t_gobsgnss>>       epoData;
		while (_crt_time <= endT)
		{
			// 如果设为false，需要在每个历元周跳探测前传入测站坐标
			if (_ephemeris)
			{
				if (_siteSPP(_crt_rec, _crt_time))
					_rec_crds[_crt_rec] = _sppflt[_crt_rec]->getCrd(_crt_time);
			}

			if(_inputEpoData.size() == 0) epoData = _allobs->obs_pt(_crt_rec, _crt_time);
			else                          epoData = _inputEpoData;

			if (epoData.size() == 0)
			{
				_crt_time = _crt_time + _interval;
				continue;
			}

			for (auto &oneObs : epoData)
			{
				string sat = oneObs->sat();
				GSYS gsys = oneObs->gsys();
				string site = oneObs->site();

				if (site != _crt_rec) continue;  // For Npp, maybe no use

				// set glonass Frequency ID
				if (sat.substr(0, 1) == "R" && oneObs->channel() >= DEF_CHANNEL)
				{
					oneObs->channel(_glonass_freq_id[sat]);
				}

				if (_frequency > 1)
				{
					set<GOBSBAND> bands = oneObs->band_avail(true);
					vector<GOBSBAND> sorted_bands = sort_band(gsys, bands);

					for (unsigned int i = 1; i < sorted_bands.size(); ++i)
					{
						if (sorted_bands[0] == sorted_bands[i]) continue;

						//_realTimeMWGF(oneObs, sorted_bands[0], sorted_bands[i]);
						_realTimeMWGF_new(oneObs, sorted_bands[0], sorted_bands[i]);
						//_realTimeMWGF_test(oneObs, sorted_bands[0], sorted_bands[i]);
					}
				}

				// Not check in the first epoch
				if (_beg_time == _crt_time)
				{
					_epoDataPre[_crt_rec][sat] = oneObs;
					continue;
				}

				if (_epoDataPre[_crt_rec].find(sat) != _epoDataPre[_crt_rec].end() && _epoDataPre_npp.size() == 0)
				{
					// Range / Doppler outliers Detection
					if (_frequency >= 2) _outliers_proc.flagRangeOutliers(_epoDataPre[_crt_rec][sat], oneObs, _interval);

					if (_frequency == 1 || _observ == RAW_MIX) _checkSlipWithDoppler(_epoDataPre[_crt_rec][sat], oneObs);
				}

				_epoDataPre[_crt_rec][sat] = oneObs;

			}
			// Next
			_crt_time = _crt_time + _interval;

		}
		return true;
	}

	void t_gexeturboedit::_realTimeMWGF(shared_ptr<t_gobsgnss> oneObs, GOBSBAND b1, GOBSBAND b2)
	{
		string sat = oneObs->sat();
		//// get the set<GOBS>
		//set<GOBS> obs_vec1 = oneObs->obs_phase(b1);
		//set<GOBS> obs_vec2 = oneObs->obs_phase(b2);
        // add select Phase
        // 避免 L1C+L2L的周跳对L1C+L2W造成影响（即 L1C L2L判定为周跳，但L2W无周跳， 实际上选用的为L1C L2W）
		GOBS obs_L1 = oneObs->select_phase(b1, true);
		GOBS obs_L2 = oneObs->select_phase(b2, true);
		set<GOBS> obs_vec1, obs_vec2;
		if (obs_L1 != GOBS::X) obs_vec1.insert(obs_L1);
		if (obs_L2 != GOBS::X) obs_vec2.insert(obs_L2);

		t_gobs gobs_P1 = t_gobs(oneObs->select_range(b1, true));
		t_gobs gobs_P2 = t_gobs(oneObs->select_range(b2, true));


		double distance, sat_clk;
		double sat_elev = _min_elev;
		double GF_factor = 1.0;
		if (this->_getSatInfo(sat, distance, sat_elev, sat_clk))
		{
			if (sat_elev < _min_elev) sat_elev = _min_elev;
			if (sat_elev < 30) GF_factor = sqrt(sin(D2R * 30) / sin(D2R * sat_elev));
		}
		else
		{
			GF_factor = 5.0;
		}

		ostringstream os;
		os.str("");

		for (auto obs1 : obs_vec1)
		{
			t_gobs gobs_L1 = t_gobs(obs1);

			for (auto obs2 : obs_vec2)
			{
				t_gobs gobs_L2 = t_gobs(obs2);

				double crt_MW_val = oneObs->MW_cycle(gobs_L1, gobs_L2, gobs_P1, gobs_P2);     // MW combination [cycles]
				double crt_GF_val = oneObs->L4(gobs_L1, gobs_L2);                             // Geometry-free combination [meters]

				os << fixed << setw(10) << " MW[c] " << oneObs->sat() << "  " << setw(6) << _crt_time.sod()
					<< "  " << gobs2str(obs1) << "  " << gobs2str(obs2)
					<< "  " << gobs2str(gobs_P1.gobs()) << "  " << gobs2str(gobs_P2.gobs()) 
					<< setw(15) << setprecision(4) << crt_MW_val << endl;

				os << fixed << setw(10) << " GF[m] " << oneObs->sat() << "  " << setw(6) << _crt_time.sod()
					<< "  " << gobs2str(obs1) << "  " << gobs2str(obs2)
					<< setw(15) << setprecision(4) << crt_GF_val << endl;

				// 如果某卫星某参考GOBS/第二个GOBS不存在（开始），重新再计算
				if (_rt_data[_crt_rec][sat][obs1].epo_num.find(obs2) == _rt_data[_crt_rec][sat][obs1].epo_num.end()
			        || _rt_data[_crt_rec][sat][obs1].mwslip[obs2] || _rt_data[_crt_rec][sat][obs1].gfslip[obs2])
				{
					_rt_data[_crt_rec][sat][obs1].epo_num[obs2] = 1;
					_rt_data[_crt_rec][sat][obs1].mwslip[obs2] = false;
					_rt_data[_crt_rec][sat][obs1].gfslip[obs2] = false;
					_rt_data[_crt_rec][sat][obs1].last_time[obs2] = _crt_time;
					_rt_data[_crt_rec][sat][obs1].GF[obs2] = crt_GF_val;
					_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2] = make_pair(crt_MW_val, 0.0);

					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].clear();
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(0.0);
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(crt_MW_val);

					if (_satdata_npp)
					{
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs1, 1);
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs2, 1);
					}
					else
					{
						oneObs->addlli(obs1, 1);
						oneObs->addlli(obs2, 1);
					}
					continue;
				}

				// 数据中断与前后历元的比较，存在BUG
				// 如  epoch->(1,2,3,4,20,21,22,...),MW->(-3.22,-3.24,-3.26,-3.25,-33.52,-33.42,-33.47,...)
				double dt = _crt_time.diff(_rt_data[_crt_rec][sat][obs1].last_time[obs2]);
				double dif_MW = crt_MW_val - _rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].back();
				// 数据中断且MW前后时刻变化大视为模糊度重置
				if ((dt > _interval && abs(dif_MW) > 4) || _check_gap && dt > _gap_limit * _interval)
				{
					_rt_data[_crt_rec][sat][obs1].epo_num[obs2] = 1;
					_rt_data[_crt_rec][sat][obs1].mwslip[obs2] = false;
					_rt_data[_crt_rec][sat][obs1].gfslip[obs2] = false;
					_rt_data[_crt_rec][sat][obs1].last_time[obs2] = _crt_time;
					_rt_data[_crt_rec][sat][obs1].GF[obs2] = crt_GF_val;
					_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2] = make_pair(crt_MW_val, 0.0);

					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].clear();
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(0.0);
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(crt_MW_val);

					ostringstream msg;
					msg << fixed << setw(10) << " Data Gap or Max_diff_MW : " << _crt_time.str_hms() << "   " << oneObs->sat() << "   " << gobs2str(obs1)
						<< "   " << gobs2str(obs2) << setw(16) << setprecision(3) << dif_MW << setw(16) << setprecision(3) << dt
						<< setw(16) << setprecision(3) << _gap_limit * _interval;
					_log->comment(1, _crt_rec, msg.str());

					if (_satdata_npp)
					{
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs1, 1);
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs2, 1);
					}
					else
					{
						oneObs->addlli(obs1, 1);
						oneObs->addlli(obs2, 1);
					}
					continue;
				}
		
				_rt_data[_crt_rec][sat][obs1].epo_num[obs2] += 1;
				_rt_data[_crt_rec][sat][obs1].last_time[obs2] = _crt_time;

				int epo_num = _rt_data[_crt_rec][sat][obs1].epo_num[obs2];

				if (_check_gf)
				{
					double pre_GF_val = _rt_data[_crt_rec][sat][obs1].GF[obs2];

					_rt_data[_crt_rec][sat][obs1].GF[obs2] = crt_GF_val;

					if (abs(crt_GF_val - pre_GF_val) > _gf_limit * GF_factor)
					{
						_rt_data[_crt_rec][sat][obs1].gfslip[obs2] = true;
						ostringstream msg;
						msg << fixed << setw(10) << " GF cycle slip : " << _crt_time.str_hms() << "   " << oneObs->sat() << "   " << gobs2str(obs1)
							<< "   " << gobs2str(obs2) << setw(16) << setprecision(3) << pre_GF_val << setw(16) << setprecision(3) << crt_GF_val
							<< setw(16) << setprecision(3) << (crt_GF_val - pre_GF_val) << setw(16) << setprecision(3) << _gf_limit * GF_factor;
						_log->comment(1, _crt_rec, msg.str());
					}

				}

				if (_check_mw)
				{
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(crt_MW_val);

					double sigma;
					double pre_MW_val = _rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].first;
					double pre_MW_std = _rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].second;
					sigma = sqrt(pre_MW_std);

					if (epo_num <= _smooth_window)
					{
						double smooth_MW_val = pre_MW_val * (epo_num - 1) / epo_num + crt_MW_val / epo_num;
						_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].first = smooth_MW_val;
						double smooth_MW_std = pre_MW_std * (epo_num - 1) / epo_num + pow(crt_MW_val - pre_MW_val, 2) / epo_num;
						_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].second = smooth_MW_std;
					}
					else
					{
						double smooth_MW_val = 0.0;
						for (int i = epo_num - _smooth_window; i <= epo_num - 1; i++)
						{
							smooth_MW_val += _rt_data[_crt_rec][sat][obs1].origin_MWs[obs2][i] / _smooth_window;
						}
						double smooth_MW_std = 0.0;
						for (int i = epo_num - _smooth_window + 1; i <= epo_num; i++)
						{
							smooth_MW_std += pow(_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2][i] - smooth_MW_val, 2) / _smooth_window;
						}
						_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].first = smooth_MW_val;
						_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].second = smooth_MW_std;
					}

					// 自适应阈值
					double k;
					if (sigma >= 0.0 && sigma <= 0.4) k = 8.0 - 5.0 * sigma;
					else if (sigma > 0.4 && sigma <= 0.6) k = 6.0;
					else if (sigma > 0.6 && sigma <= 1.0) k = 6.0 - 2.5 * (sigma - 0.6);
					else                                  k = 5.0;

					// outliers	
					if (oneObs->getobs(gobs_P1.gobs()) < 1.9e7 || oneObs->getobs(gobs_P2.gobs()) < 1.9e7)
					{
						k = 1.e9;
					}
					else if (oneObs->getoutliers(gobs_P1.gobs()) > 0 || oneObs->getoutliers(gobs_P2.gobs()) > 0)
					{
						k *= 2;
					}
					
					_mw_limit = k * sigma;

					if (abs(crt_MW_val - pre_MW_val) > _mw_limit && !double_eq(_mw_limit, 0.0))
					{
						_rt_data[_crt_rec][sat][obs1].mwslip[obs2] = true;
						if (_log)
						{
							ostringstream msg;
							msg << fixed << setw(10) << " MW cycle slip : " << _crt_time.str_hms() << "   " << oneObs->sat() << "   " << gobs2str(obs1)
								<< "   " << gobs2str(obs2) << setw(16) << setprecision(3) << pre_MW_val << setw(16) << setprecision(3) << crt_MW_val
								<< setw(16) << setprecision(3) << (crt_MW_val - pre_MW_val) << setw(16) << setprecision(3) << _mw_limit;
							_log->comment(1, _crt_rec, msg.str());
						}
					}
				}

				// 任一不合格就要设为周跳
				if (_rt_data[_crt_rec][sat][obs1].gfslip[obs2] || _rt_data[_crt_rec][sat][obs1].mwslip[obs2])
				{
					if (_satdata_npp)
					{
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs1, 1);
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs2, 1);
					}
					else
					{
						oneObs->addlli(obs1, 1);
						oneObs->addlli(obs2, 1);
					}
				}

			}
		}

		// Print flt results
		if (_debug_turbo) {
			_debug_turbo->write(os.str().c_str(), os.str().size());
			_debug_turbo->flush();
		}

	}

	void t_gexeturboedit::_realTimeMWGF_test(shared_ptr<t_gobsgnss> oneObs, GOBSBAND b1, GOBSBAND b2)
	{
		string sat = oneObs->sat();
		// get the set<GOBS>
		set<GOBS> obs_vec1 = oneObs->obs_phase(b1);
		set<GOBS> obs_vec2 = oneObs->obs_phase(b2);
		t_gobs gobs_P1 = t_gobs(oneObs->select_range(b1, true));
		t_gobs gobs_P2 = t_gobs(oneObs->select_range(b2, true));

		//xiongyun test
		double distance, elev, sat_clk;
		double fact = 1;
		if (_getSatInfo(sat, distance, elev, sat_clk))
		{
			//if (_interval >= 29.5)
			//{
			//	if (_interval <= 2.0)		fact = 1.0;
			//	else if (_interval <= 4.0)	fact = 1.25;
			//	else if (_interval <= 6.0)	fact = 1.5;
			//	else					fact = 2.0;
			//}

			if (elev < _min_elev)
			{
				elev = _min_elev;
			}
		}

		for (auto obs1 : obs_vec1)
		{
			t_gobs gobs_L1 = t_gobs(obs1);

			for (auto obs2 : obs_vec2)
			{
				t_gobs gobs_L2 = t_gobs(obs2);

				double crt_MW_val = oneObs->MW_cycle(gobs_L1, gobs_L2, gobs_P1, gobs_P2);     // MW combination [cycles]
				double crt_GF_val = oneObs->L4(gobs_L1, gobs_L2);                             // Geometry-free combination [meters]

				// 如果某卫星某参考GOBS以及第二个GOBS不存在（开始）+ 数据中断，重新再计算
				if (_rt_data[_crt_rec][sat][obs1].epo_num.find(obs2) == _rt_data[_crt_rec][sat][obs1].epo_num.end()
					|| _check_gap && _crt_time.diff(_rt_data[_crt_rec][sat][obs1].last_time[obs2]) > _gap_limit* _interval
					|| _rt_data[_crt_rec][sat][obs1].mwslip[obs2] || _rt_data[_crt_rec][sat][obs1].gfslip[obs2])
				{
					_rt_data[_crt_rec][sat][obs1].epo_num[obs2] = 1;
					_rt_data[_crt_rec][sat][obs1].mwslip[obs2] = false;
					_rt_data[_crt_rec][sat][obs1].gfslip[obs2] = false;
					_rt_data[_crt_rec][sat][obs1].last_time[obs2] = _crt_time;
					_rt_data[_crt_rec][sat][obs1].GF[obs2] = crt_GF_val;
					_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2] = make_pair(crt_MW_val, 0.0);

					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].clear();
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(crt_MW_val);

					if (_satdata_npp)
					{
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs1, 1);
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs2, 1);
					}
					else
					{
						oneObs->addlli(obs1, 1);
						oneObs->addlli(obs2, 1);
					}
					continue;
				}

				double dt = _crt_time.diff(_rt_data[_crt_rec][sat][obs1].last_time[obs2]);
				_rt_data[_crt_rec][sat][obs1].epo_num[obs2] += 1;
				_rt_data[_crt_rec][sat][obs1].last_time[obs2] = _crt_time;

				int epo_num = _rt_data[_crt_rec][sat][obs1].epo_num[obs2];

				if (_check_gf)
				{
					double pre_GF_val = _rt_data[_crt_rec][sat][obs1].GF[obs2];
					//double pre_GF_std = _rt_data[_crt_site][sat][obs1].smooth_GF[obs2].second;

					//double smooth_GF_val = pre_GF_val * (epo_num - 1) / epo_num + crt_GF_val / epo_num;
					//_rt_data[_crt_site][sat][obs1].smooth_GF[obs2].first = smooth_GF_val;
					//double smooth_GF_std = pre_GF_std * (epo_num - 1) / epo_num + pow(crt_GF_val - pre_GF_val, 2) / epo_num;
					//_rt_data[_crt_site][sat][obs1].smooth_GF[obs2].second = smooth_GF_std;

					_rt_data[_crt_rec][sat][obs1].GF[obs2] = crt_GF_val;

					//xiongyun test
					if (elev < 15)
					{
						_gf_limit = -1 * elev * 1.0 / 15.0 + 2;
					}
					else
					{
						_gf_limit = 1;
					}
					int gapepo = round(dt / _interval);
					if (gapepo >= 3 && dt >= 90.0)
					{
						_gf_limit = min(_gf_limit * gapepo, 0.35);
					}

					//if (abs(crt_GF_val - pre_GF_val) > _gf_limit * 5)
					if (abs(crt_GF_val - pre_GF_val) > _gf_limit)
					{
						_rt_data[_crt_rec][sat][obs1].gfslip[obs2] = true;
						ostringstream msg;
						msg << fixed << setw(10) << " GF cycle slip : " << _crt_time.str_hms() << "   " << oneObs->sat() << "   " << gobs2str(obs1)
							<< "   " << gobs2str(obs2) << setw(16) << setprecision(3) << pre_GF_val << setw(16) << setprecision(3) << crt_GF_val
							<< setw(16) << setprecision(3) << (crt_GF_val - pre_GF_val) << setw(16) << setprecision(3) << _gf_limit * 5;
						_log->comment(1, _crt_rec, msg.str());
					}

				}

				if (_check_mw)
				{
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(crt_MW_val);

					double sigma;
					double pre_MW_val = _rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].first;
					double pre_MW_std = _rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].second;
					sigma = sqrt(pre_MW_std);

					if (epo_num < _smooth_window)
					{
						double smooth_MW_val = pre_MW_val * (epo_num - 1) / epo_num + crt_MW_val / epo_num;
						_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].first = smooth_MW_val;
						double smooth_MW_std = pre_MW_std * (epo_num - 1) / epo_num + pow(crt_MW_val - pre_MW_val, 2) / epo_num;
						_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].second = smooth_MW_std;
					}
					else
					{
						double smooth_MW_val = 0.0;
						for (int i = epo_num - _smooth_window; i < epo_num - 1; i++)
						{
							smooth_MW_val += _rt_data[_crt_rec][sat][obs1].origin_MWs[obs2][i] / _smooth_window;
						}
						double smooth_MW_std = 0.0;
						for (int i = epo_num - _smooth_window; i < epo_num; i++)
						{
							smooth_MW_std += pow(_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2][i] - smooth_MW_val, 2) / _smooth_window;
						}
						_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].first = smooth_MW_val;
						_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].second = smooth_MW_std;

					}

					// 自适应阈值
					double k;
					if (sigma >= 0.0 && sigma <= 0.4) k = 8.0 - 5.0 * sigma;
					else if (sigma > 0.4 && sigma <= 0.6) k = 6.0;
					else if (sigma > 0.6 && sigma <= 1.0) k = 6.0 - 2.5 * (sigma - 0.6);
					else                                  k = 5.0;

					// outliers	
					if (oneObs->getoutliers(gobs_P1.gobs()) > 0 || oneObs->getoutliers(gobs_P2.gobs()) > 0)
					{
						k *= 2;
					}
					_mw_limit = k * sigma;

					//xiongyun test
					if (elev >= 20) _mw_limit = 1;
					else		   _mw_limit = -1 * 0.1 * elev + 3;

					//xiongyun test
					//if (abs(crt_MW_val - pre_MW_val) > _mw_limit && !double_eq(_mw_limit, 0.0))
					if (abs(crt_MW_val - pre_MW_val) > _mw_limit* fact && !double_eq(_mw_limit, 0.0))
					{
						_rt_data[_crt_rec][sat][obs1].mwslip[obs2] = true;
						if (_log)
						{
							ostringstream msg;
							msg << fixed << setw(10) << " MW cycle slip : " << _crt_time.str_hms() << "   " << oneObs->sat() << "   " << gobs2str(obs1)
								<< "   " << gobs2str(obs2) << setw(16) << setprecision(3) << pre_MW_val << setw(16) << setprecision(3) << crt_MW_val
								<< setw(16) << setprecision(3) << (crt_MW_val - pre_MW_val) << setw(16) << setprecision(3) << _mw_limit;
							_log->comment(1, _crt_rec, msg.str());
						}
					}
				}

				// 任一不合格就要设为周跳
				if (_rt_data[_crt_rec][sat][obs1].gfslip[obs2] || _rt_data[_crt_rec][sat][obs1].mwslip[obs2])
				{
					if (_satdata_npp)
					{
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs1, 1);
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs2, 1);
					}
					else
					{
						oneObs->addlli(obs1, 1);
						oneObs->addlli(obs2, 1);
					}
				}

			}
		}
	}

	
	void t_gexeturboedit::_realTimeMWGF_new(shared_ptr<t_gobsgnss> oneObs, GOBSBAND b1, GOBSBAND b2)
	{
		string sat = oneObs->sat();
		//// get the set<GOBS>
		//set<GOBS> obs_vec1 = oneObs->obs_phase(b1);
		//set<GOBS> obs_vec2 = oneObs->obs_phase(b2);
		// add select Phase
        // 避免 L1C+L2L的周跳对L1C+L2W造成影响（即 L1C L2L判定为周跳，但L2W无周跳， 实际上选用的为L1C L2W）
		GOBS obs_L1 = oneObs->select_phase(b1, true);
		GOBS obs_L2 = oneObs->select_phase(b2, true);
		set<GOBS> obs_vec1, obs_vec2;
		if(obs_L1 != GOBS::X) obs_vec1.insert(obs_L1);
		if(obs_L2 != GOBS::X) obs_vec2.insert(obs_L2);

		t_gobs gobs_P1 = t_gobs(oneObs->select_range(b1, true));
		t_gobs gobs_P2 = t_gobs(oneObs->select_range(b2, true));

		double distance, sat_clk;
		double sat_elev = _min_elev;
		double scale_factor = 1;
		double thres_MW = 0, thres_GF = 0;
		if (this->_getSatInfo(sat, distance, sat_elev, sat_clk))
		{
			if (sat_elev < _min_elev) sat_elev = _min_elev;
		}

		ostringstream os;
		os.str("");

		for (auto obs1 : obs_vec1)
		{
			t_gobs gobs_L1 = t_gobs(obs1);

			for (auto obs2 : obs_vec2)
			{
				t_gobs gobs_L2 = t_gobs(obs2);

				double crt_MW_val = oneObs->MW_cycle(gobs_L1, gobs_L2, gobs_P1, gobs_P2);     // MW combination [cycles]
				double crt_GF_val = oneObs->L4(gobs_L1, gobs_L2);                             // Geometry-free combination [meters]

				os << fixed << setw(10) << " MW[c] " << oneObs->sat() << "  " << setw(6) << _crt_time.sod()
					<< "  " << gobs2str(obs1) << "  " << gobs2str(obs2)
					<< "  " << gobs2str(gobs_P1.gobs()) << "  " << gobs2str(gobs_P2.gobs())
					<< setw(15) << setprecision(4) << crt_MW_val << endl;

				os << fixed << setw(10) << " GF[m] " << oneObs->sat() << "  " << setw(6) << _crt_time.sod()
					<< "  " << gobs2str(obs1) << "  " << gobs2str(obs2)
					<< setw(15) << setprecision(4) << crt_GF_val << endl;

				// 如果某卫星某参考GOBS以及第二个GOBS不存在（开始）+ 数据中断，重新再计算
				if (_rt_data[_crt_rec][sat][obs1].epo_num.find(obs2) == _rt_data[_crt_rec][sat][obs1].epo_num.end()
					|| _rt_data[_crt_rec][sat][obs1].mwslip[obs2] || _rt_data[_crt_rec][sat][obs1].gfslip[obs2])
				{
					_rt_data[_crt_rec][sat][obs1].epo_num[obs2] = 1;
					_rt_data[_crt_rec][sat][obs1].mwslip[obs2] = false;
					_rt_data[_crt_rec][sat][obs1].gfslip[obs2] = false;
					_rt_data[_crt_rec][sat][obs1].last_time[obs2] = _crt_time;
					_rt_data[_crt_rec][sat][obs1].GF[obs2] = crt_GF_val;
					_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2] = make_pair(crt_MW_val, 0.0);

					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].clear();
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(0.0);
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(crt_MW_val);

					if (_satdata_npp)
					{
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs1, 1);
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs2, 1);
					}
					else
					{
						oneObs->addlli(obs1, 1);
						oneObs->addlli(obs2, 1);
					}
					continue;
				}

				// 数据中断与前后历元的比较，存在BUG
				// 如  epoch->(1,2,3,4,20,21,22,...),MW->(-3.22,-3.24,-3.26,-3.25,-33.52,-33.42,-33.47,...)
				double dt = _crt_time.diff(_rt_data[_crt_rec][sat][obs1].last_time[obs2]);
				double dif_MW = crt_MW_val - _rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].back();
				// 数据中断且MW前后时刻变化大视为模糊度重置
				if ((dt > _interval && abs(dif_MW) > 8) || _check_gap && dt > _gap_limit * _interval)
				{
					_rt_data[_crt_rec][sat][obs1].epo_num[obs2] = 1;
					_rt_data[_crt_rec][sat][obs1].mwslip[obs2] = false;
					_rt_data[_crt_rec][sat][obs1].gfslip[obs2] = false;
					_rt_data[_crt_rec][sat][obs1].last_time[obs2] = _crt_time;
					_rt_data[_crt_rec][sat][obs1].GF[obs2] = crt_GF_val;
					_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2] = make_pair(crt_MW_val, 0.0);

					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].clear();
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(0.0);
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(crt_MW_val);


					ostringstream msg;
					msg << fixed << setw(10) << " Data Gap or Max_diff_MW : " << _crt_time.str_hms() << "   " << oneObs->sat() << "   " << gobs2str(obs1)
						<< "   " << gobs2str(obs2) << setw(16) << setprecision(3) << dif_MW << setw(16) << setprecision(3) << dt
						<< setw(16) << setprecision(3) << _gap_limit * _interval;
					_log->comment(1, _crt_rec, msg.str());

					if (_satdata_npp)
					{
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs1, 1);
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs2, 1);
					}
					else
					{
						oneObs->addlli(obs1, 1);
						oneObs->addlli(obs2, 1);
					}
					continue;
				}

				int gapNum = dt / _interval - 1;

				_rt_data[_crt_rec][sat][obs1].epo_num[obs2] += 1;
				_rt_data[_crt_rec][sat][obs1].last_time[obs2] = _crt_time;

				int epo_num = _rt_data[_crt_rec][sat][obs1].epo_num[obs2];

				if (_check_gf)
				{
					double pre_GF_val = _rt_data[_crt_rec][sat][obs1].GF[obs2];

					_rt_data[_crt_rec][sat][obs1].GF[obs2] = crt_GF_val;

					if (sat_elev < 15)
					{
						thres_GF = -_gf_limit * 1.0 / 15.0 * sat_elev + 2 * _gf_limit;
					}
					else
					{
						thres_GF = _gf_limit;
					}

					if (dt >= 90 && gapNum >= 3)
					{
						thres_GF = min(thres_GF * gapNum, 0.35);
					}

					if (abs(crt_GF_val - pre_GF_val) > thres_GF)
					{
						_rt_data[_crt_rec][sat][obs1].gfslip[obs2] = true;
						ostringstream msg;
						msg << fixed << setw(10) << " GF cycle slip : " << _crt_time.str_hms() << "   " << oneObs->sat() << "   " << gobs2str(obs1)
							<< "   " << gobs2str(obs2) << setw(16) << setprecision(3) << pre_GF_val << setw(16) << setprecision(3) << crt_GF_val
							<< setw(16) << setprecision(3) << (crt_GF_val - pre_GF_val) << setw(16) << setprecision(3) << thres_GF;
						_log->comment(1, _crt_rec, msg.str());
					}

				}

				if (_check_mw)
				{
					_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2].push_back(crt_MW_val);

					double pre_MW_val = _rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].first;
					double pre_MW_std = _rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].second;
					
					if (epo_num <= _smooth_window)
					{
						double smooth_MW_val = pre_MW_val * (epo_num - 1) / epo_num + crt_MW_val / epo_num;
						_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].first = smooth_MW_val;
						double smooth_MW_std = pre_MW_std * (epo_num - 1) / epo_num + pow(crt_MW_val - pre_MW_val, 2) / epo_num;
						_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].second = smooth_MW_std;
					}
					else
					{
						double smooth_MW_val = 0.0;
						for (int i = epo_num - _smooth_window; i <= epo_num - 1; i++)
						{
							smooth_MW_val += _rt_data[_crt_rec][sat][obs1].origin_MWs[obs2][i] / _smooth_window;
						}
						double smooth_MW_std = 0.0;
						for (int i = epo_num - _smooth_window + 1; i <= epo_num; i++)
						{
							smooth_MW_std += pow(_rt_data[_crt_rec][sat][obs1].origin_MWs[obs2][i] - smooth_MW_val, 2) / _smooth_window;
						}
						_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].first = smooth_MW_val;
						_rt_data[_crt_rec][sat][obs1].smooth_MW[obs2].second = smooth_MW_std;

					}

					if (sat_elev < 20) thres_MW = -_mw_limit * 0.1 * sat_elev + 3 * _mw_limit;
					else               thres_MW = _mw_limit;

					if (gapNum <= 2)		scale_factor = 1.0;
					else if (gapNum <= 4)	scale_factor = 1.25;
					else if (gapNum <= 8)	scale_factor = 1.5;
					else					scale_factor = 2.0;

					if (abs(crt_MW_val - pre_MW_val) > thres_MW * scale_factor && !double_eq(thres_MW, 0.0))
					{
						_rt_data[_crt_rec][sat][obs1].mwslip[obs2] = true;
						if (_log)
						{
							ostringstream msg;
							msg << fixed << setw(10) << " MW cycle slip : " << _crt_time.str_hms() << "   " << oneObs->sat() << "   " << gobs2str(obs1)
								<< "   " << gobs2str(obs2) << setw(16) << setprecision(3) << pre_MW_val << setw(16) << setprecision(3) << crt_MW_val
								<< setw(16) << setprecision(3) << (crt_MW_val - pre_MW_val) << setw(16) << setprecision(3) << thres_MW;
							_log->comment(1, _crt_rec, msg.str());
						}
					}
				}

				// 任一不合格就要设为周跳
				if (_rt_data[_crt_rec][sat][obs1].gfslip[obs2] || _rt_data[_crt_rec][sat][obs1].mwslip[obs2])
				{
					if(_satdata_npp)
					{
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs1, 1);
						(*_satdata_npp)[_map_nppdata_idx[_crt_rec][sat]].addlli(obs2, 1);
					}
					else
					{
						oneObs->addlli(obs1, 1);
						oneObs->addlli(obs2, 1);
					}
				}

			}
		}

		// Print flt results
		if (_debug_turbo) {
			_debug_turbo->write(os.str().c_str(), os.str().size());
			_debug_turbo->flush();
		}
	}


	void t_gexeturboedit::_checkSlipWithDoppler(shared_ptr<t_gobsgnss> ObsPre, shared_ptr<t_gobsgnss> Obs)
	{
		// 受到速度的影响 Todo 修改
		// 参考彭文杰博士论文，后续进行具体值的修改
		vector<GOBS> obs_vec = Obs->obs();

		double diff_time = _crt_time.diff(ObsPre->epoch());

		double factor = 1.0, factor_P = 1.0;

		double nepoch = diff_time / _interval;
		if (nepoch <= 3.0)       factor = 1.0;
		else if (nepoch <= 6.0)  factor = 1.5;
		else if (nepoch <= 12.0) factor = 2.0;
		else                     factor = 2.5;

		ostringstream os;
		os.str("");
		for (auto obs_phase : obs_vec)
		{
			// skip not phase obs
			if (!t_gobs(obs_phase).is_phase()) continue;

			// epoch Gap Judge
			if (_check_gap && diff_time > _gap_limit * _interval)
			{
				Obs->addlli(obs_phase, 1);

				ostringstream msg;
				msg << fixed << setw(10) << " cycle slip : " << _crt_time.str_hms() << "   " << Obs->sat() << "   " << gobs2str(obs_phase)
					<< "   epoch gap !!! ";
				_log->comment(1, _crt_rec, msg.str());
				continue;
			}

			double crt_phase = Obs->getobs(obs_phase);      // cycle
			double pre_phase = ObsPre->getobs(obs_phase);   // cycle

			if (double_eq(pre_phase, NULL_GOBS) || double_eq(crt_phase, NULL_GOBS)) continue;
	
			string strGOBS = gobs2str(obs_phase);
			strGOBS.replace(0, 1, "D");
			GOBS obs_doppler = str2gobs(strGOBS);
			strGOBS.replace(0, 1, "C");
			GOBS obs_range = str2gobs(strGOBS);

			//// 通常多普勒观测值存在粗差时，相位观测值也会出现周跳
			//if (Obs->getoutliers(obs_doppler) > 0)
			//{
			//	Obs->addlli(obs_phase, 1);
			//	ostringstream msg;
			//	msg << fixed << setw(10) << " cycle slip : " << _crt_time.str_hms() << "   " << Obs->sat() << "   " << gobs2str(obs_phase)
			//		<< "   doppler observation outliers !!! ";
			//	_log->comment(1, _crt_rec, msg.str());
			//	continue;
			//}

			double crt_doppler = Obs->getobs(obs_doppler);      // cycle / s
			double pre_doppler = ObsPre->getobs(obs_doppler);   // cycle / s
			if (double_eq(crt_doppler, NULL_GOBS) || double_eq(pre_doppler, NULL_GOBS)) continue;

			double delta_LD = (crt_phase - pre_phase) + 0.5 * (crt_doppler + pre_doppler) * diff_time;

			os << fixed << setw(10) << " LD[c] " << Obs->sat() << "  " << setw(6) << ObsPre->epoch().sod() << setw(6) << _crt_time.sod() <<
				"  " << gobs2str(obs_phase) << "  " << gobs2str(obs_doppler) << setw(15) << setprecision(4) << delta_LD << endl;

			double crt_range = Obs->getobs(obs_range);      // m
			double pre_range = ObsPre->getobs(obs_range);   // m
			if (double_eq(crt_range, NULL_GOBS) || double_eq(pre_range, NULL_GOBS)) continue;

			crt_phase *= Obs->wavelength(t_gobs(obs_doppler).band());  // m 
			pre_phase *= Obs->wavelength(t_gobs(obs_doppler).band());  // m
			double delta_PL = (crt_range - crt_phase) - (pre_range - pre_phase);

			os << fixed << setw(10) << " PL[m] " << Obs->sat() << "  " << setw(6) << ObsPre->epoch().sod() << setw(6) << _crt_time.sod() <<
				"  " << gobs2str(obs_phase) << "  " << gobs2str(obs_range) << setw(15) << setprecision(4) << delta_PL << endl;

			if (Obs->getoutliers(obs_range) > 0 || ObsPre->getoutliers(obs_range) > 0) factor_P *= 2;

			// 为了避免极少数情况下多普勒观测值存在粗差并将相位观测值误判为存在周跳,阈值设为5周
			if (abs(delta_LD) > 5 * factor && abs(delta_PL) > 3 * factor * factor_P)
			{
				Obs->addlli(obs_phase, 1);
				ostringstream msg;
				msg << fixed << setw(10) << " cycle slip : " << _crt_time.str_hms() << "   " << Obs->sat() << "   " << gobs2str(obs_phase) << "  " << gobs2str(obs_doppler) << "  " << gobs2str(obs_range)
					<< setw(16) << setprecision(3) << delta_LD << " > 5 cycles[delta_LD],  " << setw(16) << setprecision(3) << delta_PL << "  > 3 meters[delta_PL] ";
				_log->comment(1, _crt_rec, msg.str());
			}
		}

		// Print flt results
		if (_debug_turbo) {
			_debug_turbo->write(os.str().c_str(), os.str().size());
			_debug_turbo->flush();
		}

	}

	void t_gexeturboedit::_initPostProcess()
	{
		// TurboEdit Settings for post-processing
		_site_kinematic = dynamic_cast<t_gsetproc*>(_set)->pos_kin();
		_simulation = dynamic_cast<t_gsetproc*>(_set)->simulation();
		_amb_output = dynamic_cast<t_gsetturboedit*>(_set)->isAmbOutput();
		//_simulation = dynamic_cast<t_gsetturboedit*>(_set)->isSimulation();
		_check_sf = dynamic_cast<t_gsetturboedit*>(_set)->checkSingleFreq(_sf_limit);
		_check_pc = dynamic_cast<t_gsetturboedit*>(_set)->checkPC(_pc_limit);
		_check_short = dynamic_cast<t_gsetturboedit*>(_set)->checkShort(_short_limit);
		_check_statistics = dynamic_cast<t_gsetturboedit*>(_set)->checkStatistics(_min_percent, _min_mean_nprn, _max_mean_namb);

		// Single frequency and combination are not compatible
		if (_check_mw || _check_gf) _check_sf = false;
		if (!_check_short) _short_limit = 0;

		// set GNSS bands
		for (auto it : _sys_list)
		{
			GSYS gsys = t_gsys::str2gsys(it);
			if (_simulation && gsys == GSYS::LEO) continue;
			std::set<GOBSBAND> s_band;                           // insert bands
			vector<GOBSBAND> v_band = dynamic_cast<t_gsetgnss*>(_set)->band(gsys);
			for (auto it_band : v_band) s_band.insert(it_band);
			v_band = sort_band(gsys, s_band);

			_freq_index[gsys] = dynamic_cast<t_gsetgnss*>(_set)->freq_index(gsys);

			if (_check_sf)
			{
				_post_all_band[1][gsys] = make_pair(v_band[0], v_band[0]);
			}
			else
			{
				for (unsigned int i = 1; i < v_band.size(); ++i)
				{
					_post_all_band[i][gsys] = make_pair(v_band[0], v_band[i]);
				}
			}
		}

		// Get site coordinates from XML/rinex O files
		set<string> rec_rm;
		for (const auto& it : _rec_list)
		{
			// XML 
			_rec_crds[it] = dynamic_cast<t_gsetrec*>(_set)->get_crd_xyz(it);
			// SNX
			if (_rec_crds[it].zero() && _allobj)
			{
				shared_ptr<t_gobj> _grec = _allobj->obj(it);
				t_gtriple crd(0.0, 0.0, 0.0);
				t_gtriple std(0.0, 0.0, 0.0);
				if (_grec)
				{
					_grec->get_recent_crd(_beg_time, 0.1, crd, std);
					_rec_crds[it] = crd;
				}
			}
			// Rinex O
			if (_rec_crds[it].zero() && _allobs)
			{
				_rec_crds[it] = _allobs->getsitecrd(it);
			}
			// delete it
			if (_rec_crds[it].zero())
			{
				rec_rm.insert(it);
			}
		}

		for (auto it : rec_rm) 
		{
			_rec_list.erase(it);
			if (_sppflt[it]) { delete _sppflt[it]; _sppflt[it] = nullptr; } // SPPFLT delete
		}

		// Renew satlist according brdm files
		if (_ephemeris)
		{
			// clean invalid/duplicate NAV messages
			_allnav->clean_invalid();
			_allnav->clean_duplicit();
			_glonass_freq_id = _allnav->glo_freq_num();
			set<string> sat_xml = _sat_list;
			_sat_list.clear();
			set<string> sat_temp = _allnav->satellites();
			set_intersection(sat_xml.begin(), sat_xml.end(), sat_temp.begin(), sat_temp.end(), inserter(_sat_list, _sat_list.begin()));
		}

		// Set scaling factor for large cycle checks
		_scaling_factor["LX"] = 1.0;  // single frequency
		_scaling_factor["MW"] = 5.0;
		if (_simulation)
		{
			_ephemeris = false;       // if simulation, noy use ephemeris
			_scaling_factor["GF"] = 20.0;
		}
		else
		{
			_scaling_factor["GF"] = 5.0;
		}

		// init ambflag
		if (_allambflag) { delete _allambflag; _allambflag = nullptr; }
		_allambflag = new t_gallambflag(t_gdata::AMBFLAG);

	}

	bool t_gexeturboedit::_postProcess()
	{   
		 _log->comment(1, _crt_rec, " --- Start TURBOEDIT --- ");
		 
		 // Set log suffix
		 string log_suffix = ".log";
		 string cerr_out = " Process Band : ";
		 for (auto it : _sys_list)
		 {
			 GSYS gsys = t_gsys::str2gsys(it);
			 if (_post_crt_band.find(gsys) == _post_crt_band.end()) continue;
			 if (_simulation && gsys == GSYS::LEO) continue;
			 string str_freq1 = gfreqseq2str(_freq_index[gsys][_post_crt_band[gsys].first]);
			 string str_freq2 = gfreqseq2str(_freq_index[gsys][_post_crt_band[gsys].second]);
			 string str_freq = "L" + str_freq1 + str_freq2;
			 if (str_freq1 > "2" || str_freq2 > "2") log_suffix = ".log" + str_freq1 + str_freq2;
			 
			 cerr_out += " " + it + "->BAND_" + gobsband2str(_post_crt_band[gsys].first) + " ";

			 if (!_check_sf) cerr_out += "+ BAND_" + gobsband2str(_post_crt_band[gsys].second) + " ";			 
		 }
		 cerr << cerr_out << endl;


		 // Init turboEdit Information for each site
		 cerr << " --- Init Turboedit Information --- " << endl;
		 _initPostTurboInfo();

		 //  Basic check
		 cerr << " --- Basic Check --- " << endl;
		 _postBasicCheck();

		 // Remove short piece and flag gap
		 cerr << " --- Check Gap and Remove Short --- " << endl;
		 _log->comment(1, _crt_rec, " --- First remove short, short_limit = " + dbl2str(_short_limit * _interval / 2.0) + "s --- ");
		 for (auto sat : _sat_list)
		 {
			 _log->comment(1, _crt_rec, "## remove short for satellite, " + sat);
			 while (true)
			 {
				 if (!_postCheckGapShort(_short_limit * _interval / 2.0, sat)) break;
			 }
		 }

     	 // check_lwlg check WN and LG for all satellites
		 if (!_simulation && (_check_gf || _check_mw))
		 {
			 cerr << " --- Turboedit Check --- " << endl;
			 _log->comment(1, _crt_rec, " --- check MW or GF for all satellites --- ");
			 for (auto sat : _sat_list)
			 {
				 _log->comment(1, _crt_rec, "## turboedit for satellite " + sat);
				 _postCheckMWGF(sat);
			 }
			 // remove short piece and flag gap after turboedit check
			 cerr << " --- Check Gap and Remove Short --- " << endl;
			 _log->comment(1, _crt_rec, " --- Second remove short, short_limit = " + dbl2str(_short_limit * _interval / 2.0) + "s --- ");
			 for (auto sat : _sat_list)
			 {
				 _log->comment(1, _crt_rec, "## remove short for satellite, " + sat);
				 while (true)
				 {
					 if (!_postCheckGapShort(_short_limit * _interval / 2.0, sat)) break;
				 }
			 }		 
		 }

     	 // make the decision whether the staion is usable and write the report
         // in case there are two many ambiguities, try to remove some of the short pieces by incresing `ength_short`		
     	 string decision = "";
     	 double new_length_short = _short_limit * _interval;
     	 while (new_length_short < 7200)
     	 {
			 cerr << " --- Print Process Report --- " << endl;
     	     _postPrintReport(decision);
     	 
     	     if (decision != "BAD too many AMBs") break;
     	 
     	     cerr << " --- Check Gap and Remove Short --- " << endl;
     	     new_length_short *= 1.2;
     	     _log->comment(1, _crt_rec, "Too many ambiguities and try to remove some short piece. short_length="
     	 	  + dbl2str(new_length_short) + "s");
     	 
     	     for (auto sat : _sat_list)
     	     {
				 _log->comment(1, _crt_rec, "## remove short for satellite, " + sat);
     	 	      while (true)
     	 	      {
				  if (!_postCheckGapShort(new_length_short, sat)) break;
     	 	      }
     	     }
     	 }
     	 	
     	 // set log-files
     	 _allambflag->clearSiteAmbFlag(_crt_rec);
     	 cerr << " --- Set Single Ambflag --- " << endl;
     	 if (!_postSetAmbflag()) return false;
     
     	 // write log-files
		 if (_amb_output)
		 {
			 cerr << " --- Encoder Product: " + _crt_rec + " --- " << endl;
     	 	 if (!_postEncoderProduct(log_suffix)) return false;
		 }
		 return true;
	}

	void t_gexeturboedit::_initPostTurboInfo()
	{
		// re-initialize _post_data/_post_time
		int iepo;
		int flag = _setFlag(0, "NODATA");
		for (auto sat : _sat_list)
		{		
			for (iepo = 1; iepo <= _epo_num; iepo++)
			{
				_post_data[sat].obs[iepo][1] = 0.0;    // geometry-free combination  cycle
				_post_data[sat].obs[iepo][2] = 0.0;    // Melbourne-Wuebenna combination cycle
				_post_data[sat].obs[iepo][3] = 0.0;    // ionosphere-free combination for phase cycle
				_post_data[sat].obs[iepo][4] = 0.0;    // ionosphere-free combination for range m
				_post_data[sat].obs[iepo][5] = 0.0;    // residuals of geometry-free combination
				_post_data[sat].obs[iepo][6] = 999.0;  // sigma of geometry-free combination
				_post_data[sat].elev[iepo] = 0.0;
				_post_data[sat].iflag[iepo] = flag;
			}
			// init last valid epoch for every satellite
			_post_data[sat].last_epo = 0;
		}

		for (iepo = 1; iepo <= _epo_num; iepo++)
		{

			_post_time[iepo] = _beg_time + _interval * (iepo - 1);
		}

		_log->comment(1, _crt_rec, "Turboedit Information Initialization Success For Post-Processing!");
	}

	bool t_gexeturboedit::_siteSPP(const string& site, const t_gtime& now)
	{
		// SPPFlt for station
		if (_sppflt[site]->processBatch(now, now) < 0)
		{
			_log->comment(1, site, "Kinematic SPPFlt Fail at " + now.str_ymdhms());
			return false;
		}
		return true;
	}

	void t_gexeturboedit::_postBasicCheck()
	{
		_crt_time = _beg_time;
		_crt_epo = 1;
		int flag_tmp;
		while (_crt_time <= _end_time)
		{	
			// cerr << _crt_time.str_ymdhms() << endl;

			// For Kinematic Stations
			if ((_site_kinematic || _rec_crds[_crt_rec].zero()) && _ephemeris)
			{
				if (_siteSPP(_crt_rec, _crt_time)) 
					_rec_crds[_crt_rec] = _sppflt[_crt_rec]->getCrd(_crt_time);
			}

			/* get observation */
			vector<t_gsatdata> crt_obs = _allobs->obs(_crt_rec, _crt_time);

			// loop observations
			for (auto &it_obs : crt_obs)
			{
				string sat = it_obs.sat();
				GSYS sys = it_obs.gsys();

				// exclude sats not in _sat_list
				if (_sat_list.find(sat) == _sat_list.end()) continue;

				// with data, set its flag to good at first
				_post_data[sat].iflag[_crt_epo] = DATAFLAG::GOOD;

				// set glonass Frequency ID
				if (sat.substr(0, 1) == "R" && it_obs.channel() >= DEF_CHANNEL)
				{
					it_obs.channel(_glonass_freq_id[sat]);
				}

				//Todo
				//check wether after shadow for GPS of BlockType II / IIA

				// check observation exist or not, and obs combination
				if (!_postPrepareObs(it_obs)) continue;

				if (!_simulation && _ephemeris)
				{
					double distance, elev, sat_clk;

					// compute satellite position, distance, elevation
					if (!_getSatInfo(sat, distance, elev, sat_clk))
					{
						flag_tmp = _post_data[sat].iflag[_crt_epo];
						_post_data[sat].iflag[_crt_epo] = _setFlag(flag_tmp, "NOBRD");
						continue;
					}

					_post_data[sat].obs[_crt_epo][4] -= (distance - sat_clk);
					_post_data[sat].elev[_crt_epo] = elev;

					if (elev < _min_elev)
					{
						flag_tmp = _post_data[sat].iflag[_crt_epo];
						_post_data[sat].iflag[_crt_epo] = _setFlag(flag_tmp, "LOWELE");
						_log->comment(2, _crt_rec, sat + " has low elevation = " + dbl2str(elev) + " at Epoch " + int2str(_crt_epo));
						continue;
					}
				}
				
				//check widelane and ionosphere for detecting large cycle slips
				_postCheckLargeSlips(sat);

				// check loss-of-lock indicator
				if (it_obs.getlli(_post_phase_obs[sys][0]) >= 1 || it_obs.getlli(_post_phase_obs[sys][1]) >= 1)
				{
					flag_tmp = _post_data[sat].iflag[_crt_epo];
					_post_data[sat].iflag[_crt_epo] = _setFlag(flag_tmp, "LLI");
					_log->comment(2, _crt_rec, "Exist loss-of-lock For " + sat + " in Epoch " + int2str(_crt_epo));
				}
		
			}  // end loop observations

			// check range
			if (!_simulation && _ephemeris && _check_pc) _postCheckRange();

			// Update time & epoch
			_crt_time = _crt_time + _interval;
			_crt_epo = floor((_crt_time - _beg_time) / _interval) + 1;
		}
	}

	bool t_gexeturboedit::_postPrepareObs(t_gobsgnss& obsdata)
	{
		GSYS gsys = obsdata.gsys();
		string sat = obsdata.sat();

		// For GLONASS when 3 frequency 
		if (_post_crt_band.find(gsys) == _post_crt_band.end()) return false;

		GOBS obsL1 = obsdata.select_phase(_post_crt_band[gsys].first);
		GOBS obsP1 = obsdata.select_range(_post_crt_band[gsys].first);
		double P1 = obsdata.getobs(obsP1);
		double L1 = obsdata.getobs(obsL1);
		_post_phase_obs[gsys].clear();
		_post_phase_obs[gsys].push_back(obsL1);
		_post_phase_obs[gsys].push_back(obsL1);

		// check single frequency
		if (_check_sf)
		{
			if (double_eq(P1, 0.0) || double_eq(L1, 0.0))
			{
				int flag_tmp = _post_data[sat].iflag[_crt_epo];
				_post_data[sat].iflag[_crt_epo] = _setFlag(flag_tmp, "NO4");
				_log->comment(2, _crt_rec, sat + " lack of full data at Epoch " + int2str(_crt_epo));
				return false;
			}
			_post_data[sat].obs[_crt_epo][3] = L1 * obsdata.frequency(_post_crt_band[gsys].first);  // L1 [m]
			_post_data[sat].obs[_crt_epo][4] = P1;                                                  // P1 [m]
			_post_data[sat].obs[_crt_epo][1] = _post_data[sat].obs[_crt_epo][3] - P1;               // L1-P1 [m]
		}
		// check MW+GF
		if (_check_mw || _check_gf)
		{
			GOBS obsL2 = obsdata.select_phase(_post_crt_band[gsys].second);
			GOBS obsP2 = obsdata.select_range(_post_crt_band[gsys].second);
			double P2 = obsdata.getobs(obsP2);
			double L2 = obsdata.getobs(obsL2);
			_post_phase_obs[gsys].back() = obsL2;

			if (double_eq(P1, 0.0) || double_eq(L1, 0.0) || double_eq(P2, 0.0) || double_eq(L2, 0.0))
			{
				int flag_tmp = _post_data[sat].iflag[_crt_epo];
				_post_data[sat].iflag[_crt_epo] = _setFlag(flag_tmp, "NO4");
				_log->comment(2, _crt_rec, sat + " lack of full data at Epoch " + int2str(_crt_epo));
				return false;
			}

			t_gobs gobsL1 = t_gobs(obsL1);
			t_gobs gobsL2 = t_gobs(obsL2);
			t_gobs gobsP1 = t_gobs(obsP1);
			t_gobs gobsP2 = t_gobs(obsP2);

			_post_data[sat].obs[_crt_epo][1] = obsdata.L4_cycle(gobsL1, gobsL2);                   // Geometry-free combination [cycle]
			_post_data[sat].obs[_crt_epo][2] = obsdata.MW_cycle(gobsL1, gobsL2, gobsP1, gobsP2);   // MW combination [cycle]
			_post_data[sat].obs[_crt_epo][3] = obsdata.L3_cycle(gobsL1, gobsL2);                   // LC combination [cycle]
			_post_data[sat].obs[_crt_epo][4] = obsdata.P3(gobsP1, gobsP2);                         // PC combination [m]

		}

		return true;

	}

	bool t_gexeturboedit::_getSatInfo(string sat, double& distance, double& elev, double& sat_clk)
	{
		bool valid = false;
		double sat_crd[3];
		t_gtriple sat_crd_triple, site_sat_vector;

		double delay = 0.0;
		t_gtime sat_time;
		shared_ptr<t_geph> geph;
		while (true)
		{
			sat_time = _crt_time - delay;
			geph = _allnav->find(sat, sat_time, false);  // add false to not mind healthy status, test by glfeng
			// no broadcast ephemeris
			if (!geph) break;
			// invalid broadcast ephemeris
			if (fabs(sat_time.diff(geph->epoch()) / 3600.0) > 4.0) break;

			// can't get sat_crd from broadcast ephemeris
			if (geph->pos(sat_time, sat_crd, NULL, NULL, false) < 0 ) break; // add false to not mind healthy status, test by glfeng

			sat_crd_triple = t_gtriple(sat_crd);

			if (_rec_crds[_crt_rec].zero()) break;

			site_sat_vector = sat_crd_triple - _rec_crds[_crt_rec];
			distance = site_sat_vector.norm();

			if (abs(distance / CLIGHT - delay) < 1E-9)
			{
				valid = true;
				break;
			}

			delay = distance / CLIGHT;
		}

		// Compute elevation & sat clks
		if (valid)
		{

			if (geph->clk(sat_time, &sat_clk, NULL, NULL, false) >= 0) // add false to not mind healthy status, test by glfeng
			{
				sat_clk *= CLIGHT;
				// CRS
				elev = (_rec_crds[_crt_rec][0] * site_sat_vector[0] + _rec_crds[_crt_rec][1] * site_sat_vector[1]
					+ _rec_crds[_crt_rec][2] * site_sat_vector[2]) / _rec_crds[_crt_rec].norm() / distance;
				elev = 90.0 - acos(elev) * 180.0 / G_PI;
				
				return true;
			}
		}

		distance = 0.0; sat_clk = 0.0, elev = 0.0;
		_log->comment(2, _crt_rec, sat + " lack of suitable broadcast ephemeris at Time " + _crt_time.str_ymdhms());
		return false;
	}

	void t_gexeturboedit::_postCheckLargeSlips(string sat)
	{
		// check widelane and ionosphere for detecting large cycle slips
		if (_post_data[sat].last_epo >= 1)
		{
			int flag_tmp;

			double dummy;
			if (_check_mw)   // widelane check [cycle]
			{				
				dummy = fabs(_post_data[sat].obs[_crt_epo][2] - _post_data[sat].obs[_post_data[sat].last_epo][2]);
				if (dummy > _mw_limit * _scaling_factor["MW"])
				{
					flag_tmp = _post_data[sat].iflag[_crt_epo];
					_post_data[sat].iflag[_crt_epo] = _setFlag(flag_tmp, "MWJUMP");
					_log->comment(1, _crt_rec, "Bad Widelane at Epoch " + int2str(_crt_epo) + " for "
						+ sat + ", it residuals : " + dbl2str(dummy));
				}
			}
			if (_check_gf)
			{
				// ionosphere check [cycle]
				dummy = fabs(_post_data[sat].obs[_crt_epo][1] - _post_data[sat].obs[_post_data[sat].last_epo][1])
					/ _crt_time.diff(_post_time[_post_data[sat].last_epo]);
				if (!_simulation && _post_data[sat].elev[_crt_epo] < 45.0)
				{
					dummy *= sin(_post_data[sat].elev[_crt_epo] / 180.0 * G_PI);
				}
				if (dummy > _gf_limit * _scaling_factor["GF"])
				{
					flag_tmp = _post_data[sat].iflag[_crt_epo];
					_post_data[sat].iflag[_crt_epo] = _setFlag(flag_tmp, "GFJUMP");
					_log->comment(1, _crt_rec, "Bad Ionosphere at Epoch " + int2str(_crt_epo) + " for "
						+ sat + ", its residuals and ele: " + dbl2str(dummy) + "  " +
						dbl2str(_post_data[sat].elev[_crt_epo]));
				}
			}
			if (_check_sf)
			{
				dummy = fabs(_post_data[sat].obs[_crt_epo][1] - _post_data[sat].obs[_post_data[sat].last_epo][1]);  // [m] 
				if (dummy > _sf_limit * _scaling_factor["LX"])
				{
					flag_tmp = _post_data[sat].iflag[_crt_epo];
					_post_data[sat].iflag[_crt_epo] = _setFlag(flag_tmp, "MWJUMP");
					_log->comment(1, _crt_rec, "Bad LX-CX at Epoch " + int2str(_crt_epo) + " "
						+ sat + " residuals : " + dbl2str(dummy));
				}
			}
		}
		_post_data[sat].last_epo = _crt_epo;
	}

	void t_gexeturboedit::_postCheckRange()
	{		
		vector<double> pc_res;
		vector<string> useful_sats;
		for (auto it_sat : _sat_list)
		{
			if (!_isTrue(_post_data[it_sat].iflag[_crt_epo], "OK")) continue;
			pc_res.push_back(_post_data[it_sat].obs[_crt_epo][4]);
			useful_sats.push_back(it_sat);
		}

		int flag_tmp;
		if (pc_res.size() == 1 || (pc_res.size() == 2 && fabs(pc_res[0] - pc_res[1]) > 40.0))
		{
			for (auto sat : useful_sats)
			{
				flag_tmp = _post_data[sat].iflag[_crt_epo];
				_post_data[sat].iflag[_crt_epo] = _setFlag(flag_tmp, "PCBAD");
			}
			return;
		}

		double mean = 0.0;
		double sig = 0.0;
		vector<double> res_tmp;

		for (auto it : pc_res) { mean += it; }
		mean /= pc_res.size();

		for (auto it : pc_res)
		{
			sig += pow(it - mean, 2);
			res_tmp.push_back(fabs(it - mean));
		}
		sig = sqrt(sig / (pc_res.size() - 1));

		double sig1 = 2000;
		double mean1 = 0.0;
		bool isFound = true;
		set<string> bad_sats;
		while (isFound && sig1 > 1000 && pc_res.size() >= 3)
		{
			// find out the max.residial
			auto maxRes = max_element(res_tmp.begin(), res_tmp.end());
			int maxPos = distance(res_tmp.begin(), maxRes);
			double MaxValue = pc_res[maxPos];

			// get mean and sigma after removal of the largest residual.
			mean1 = 0.0;
			sig1 = 0.0;
			res_tmp.clear();

			for (auto it : pc_res)
			{
				if (double_eq(it, MaxValue)) continue;
				mean1 += it / (pc_res.size() - 1);
			}

			for (auto it : pc_res)
			{
				if (double_eq(it, MaxValue)) continue;
				sig1 += (it - mean1) * (it - mean1);
				res_tmp.push_back(fabs(it - mean1));
			}
			sig1 = sqrt(sig1 / (pc_res.size() - 2));

			isFound = sig / sig1 > 2.5 || fabs(mean - mean1) > 200.0
				|| fabs(MaxValue - mean1) > (3 * (sig1 > 100.0 ? sig1 : 100.0));
			if (!isFound) continue;
			mean = mean1;
			sig = sig1;
			bad_sats.insert(useful_sats[maxPos]);
			pc_res.erase(pc_res.begin() + maxPos);
			useful_sats.erase(useful_sats.begin() + maxPos);
		}

		int nbad = 0;
		for (auto sat : _sat_list)
		{
			if (bad_sats.find(sat) == bad_sats.end()) continue;
			if (fabs(_post_data[sat].obs[_crt_epo][4]) <= _pc_limit) continue;
			flag_tmp = _post_data[sat].iflag[_crt_epo];
			_post_data[sat].iflag[_crt_epo] = _setFlag(flag_tmp, "PCBAD");
			_log->comment(1, _crt_rec, "Flagged Bad PC " + sat + "  at Epoch "
				+ int2str(_crt_epo) + "  pc_residuals and threshold: " +
				dbl2str(_post_data[sat].obs[_crt_epo][4]) + dbl2str(_pc_limit));
			nbad++;
		}

		if (nbad >= 3)
		{
			_log->comment(1, _crt_rec, "Too many bad ranges, be sure the coordinates used are good enough!");
		}

		if (pc_res.size() != 0) _post_time[_crt_epo] = _post_time[_crt_epo] - mean / CLIGHT;

		return;
	}

	bool t_gexeturboedit::_postCheckGapShort(double length_short, string  sat)
	{
		// remove short piece
		bool is_removed = false;

		vector<int> search_flag_data;
		search_flag_data.push_back(-999);   // no use
		int last_valid_epo = 0;
		int iepo = 0;
		// mark large gap and first epoch
		for (iepo = 1; iepo <= _epo_num; iepo++)
		{
			search_flag_data.push_back(_post_data[sat].iflag[iepo]);
			if (!_isTrue(_post_data[sat].iflag[iepo], "OK")) continue;  // not usable data

			if (last_valid_epo == 0)
			{
				_post_data[sat].iflag[iepo] = _setFlag(_post_data[sat].iflag[iepo], "GAP");
			}
			else
			{
				if (_post_time[iepo].diff(_post_time[last_valid_epo]) > _gap_limit * _interval)
				{
					_post_data[sat].iflag[iepo] = _setFlag(_post_data[sat].iflag[iepo], "GAP");
					_log->comment(1, _crt_rec, "Flagged gap at " + int2str(iepo)
						+ ", last good epoch  " + int2str(last_valid_epo));
				}
			}
			last_valid_epo = iepo;
			search_flag_data[iepo] = _post_data[sat].iflag[iepo];  // update iflag
		}

		// no usable data
		if (last_valid_epo <= 0)
		{
			_log->comment(1, _crt_rec, "No data to remove: " + sat);
			return is_removed;
		}

		int j = 1;
		int i, k;
		while (j <= _epo_num) {
			// find the first ambiguity
			i = _findFlag(j, _epo_num, search_flag_data, "AMB", "+");
			if (i == -999)  break;

			// find the next ambiguity or end of the data set
			j = _findFlag(i + 1, _epo_num, search_flag_data, "AMB", "+");
			if (j == -999) j = _epo_num + 1;

			// find backwards the last usable epoch
			k = _findFlag(i, j - 1, search_flag_data, "OK", "-");

			// count the number of usable
			int ngood = 1;  // for epoch i
			for (iepo = i + 1; iepo <= k; iepo++)
			{
				if (_isTrue(search_flag_data[iepo], "OK")) ngood++;
			}

			// Removed short
			if (_post_time[k].diff(_post_time[i + 1]) < length_short
				|| ngood * _interval <= length_short / 2)
			{
				_log->comment(1, _crt_rec, "Removed short, from epoch "
					+ int2str(i) + " to " + int2str(k));
				for (iepo = i; iepo <= k; iepo++)
				{
					if (_isTrue(search_flag_data[iepo], "OK"))
					{
						_post_data[sat].iflag[iepo] = _setFlag(_post_data[sat].iflag[iepo], "SHRT");
						search_flag_data[iepo] = _post_data[sat].iflag[iepo];
					}
				}
				is_removed = true;
			}

			j = k + 1;
		}
		return is_removed;
	}

	void t_gexeturboedit::_postCheckMWGF(string sat)
	{
		int nused = 0;
		for (int iepo = 1; iepo <= _epo_num; iepo++)
		{
			if (!_isTrue(_post_data[sat].iflag[iepo], "OK")) continue;
			nused++;
		}

		if (nused > 0)
		{
			if(_check_mw)  _postEditWidelane(sat);
			if(_check_gf)  _postCheckIonosphere(sat);
		}

		_log->comment(1, _crt_rec, " --- After MW+GF check, short_limit = " + dbl2str(_short_limit * _interval) + "s --- ");
		while (true)
		{
			if (!_postCheckGapShort(_short_limit * _interval, sat)) break;
		}
		return;
	}

	void t_gexeturboedit::_postEditWidelane(string sat)
	{
		vector<int> search_flg;
		search_flg.push_back(999);   // no use
		int iepo = 0;
		// simply check for jump
		for (iepo = 1; iepo <= _epo_num; iepo++)
		{
			search_flg.push_back(3);
			if (!_isTrue(_post_data[sat].iflag[iepo], "OK")) continue; // OK = GOOD+AMB
			search_flg[iepo] = 1;
			if (!_isTrue(_post_data[sat].iflag[iepo], "MWJUMP") && !_isTrue(_post_data[sat].iflag[iepo], "GAP"))       // 排除上述情况，剩下的基本就是正常的MW值
			{
				search_flg[iepo] = 0;
			}
		}

		//  sequential mean of NW and flag possible are apply to detect detect cycle slips 
		//          ( by testing  dabs(nwi-meani) > limit cycles)
		// if a data point is flagged, mean calculation is re-initialized by finding the
		// next two points with a dif.less than limit.
		int kobs = 0, last_epo = 0;
		double mean_mw, mw_diff, mw_adjacent_epo;
		for (iepo = 1; iepo <= _epo_num; iepo++)
		{
			if (search_flg[iepo] > 1) continue;    // data with bad point
			if (search_flg[iepo] == 1) kobs = 0;   // possible exsit LWJUMP/GAP, piece start

			if (last_epo >= 1)
			{
				mw_adjacent_epo = _post_data[sat].obs[last_epo][2] - _post_data[sat].obs[iepo][2];
			}

			if (kobs == 0)                              // first point. piece start
			{
				mean_mw = _post_data[sat].obs[iepo][2];
				search_flg[iepo] = 1;
				kobs = 1;
			}
			else if (kobs == 1)                        // second point
			{
				mw_diff = _post_data[sat].obs[iepo][2] - mean_mw;
				if (fabs(mw_diff) <= _mw_limit)        // two points are ok
				{
					kobs++;
					mean_mw = mean_mw + mw_diff / kobs;
				}
				else
				{
					_log->comment(1, _crt_rec, "Widelane bad point, " + int2str(iepo)
						+ "  " + dbl2str(mw_adjacent_epo));
					search_flg[last_epo] = 2;                 // bad point = 2
					search_flg[iepo] = 1;                     // piece start = 1
					mean_mw = _post_data[sat].obs[iepo][2];
				}
			}
			else                                       // third or later point
			{
				mw_diff = _post_data[sat].obs[iepo][2] - mean_mw;
				if (fabs(mw_diff) <= _mw_limit ||
					fabs(mw_diff) <= 1.4 * _mw_limit && fabs(mw_adjacent_epo) <= 0.5)
				{
					kobs++;
					mean_mw = mean_mw + mw_diff / kobs;
					if (search_flg[last_epo] == 1)
					{
						_log->comment(1, _crt_rec, "Widelane bad point, " + int2str(last_epo)
							+ "  " + dbl2str(mw_adjacent_epo));
						search_flg[last_epo] = 2;
					}
				}
				else                                  // possible cycle slip
				{
					if (search_flg[last_epo] == 1)           // set crt epo as second point
					{
						_log->comment(1, _crt_rec, "Widelane possible jump, " + int2str(last_epo));
						kobs = 1;
						mean_mw = _post_data[sat].obs[last_epo][2];
						mw_diff = _post_data[sat].obs[iepo][2] - mean_mw;
						if (fabs(mw_diff) <= _mw_limit)
						{
							kobs++;
							mean_mw = mean_mw + mw_diff / kobs;
						}
						else                         // set crt epo as bad point
						{
							_log->comment(1, _crt_rec, "Widelane bad point, " + int2str(iepo)
								+ "  " + dbl2str(mw_adjacent_epo));
							search_flg[last_epo] = 2;
							search_flg[iepo] = 1;
							mean_mw = _post_data[sat].obs[iepo][2];
						}
					}
					else
					{
						search_flg[iepo] = 1;
					}					
				}
			}
			last_epo = iepo;
		}

		// set flag
		for (iepo = 1; iepo <= _epo_num; iepo++)
		{
			if (search_flg[iepo] == 1)
			{
				if (!_isTrue(_post_data[sat].iflag[iepo], "MWJUMP")
					&& !_isTrue(_post_data[sat].iflag[iepo], "GAP"))
				{
					_post_data[sat].iflag[iepo] = _setFlag(_post_data[sat].iflag[iepo], "MWJUMP");
				}
			}
			else if (search_flg[iepo] == 2)
			{
				_post_data[sat].iflag[iepo] = _setFlag(_post_data[sat].iflag[iepo], "MWBAD");
			}
		}

	}

	void t_gexeturboedit::_postCheckIonosphere(string sat)
	{
		vector<int> search_flg;
		search_flg.push_back(999);  // no use, for counts
		int iepo = 0;
		// simply check for jump
		for (iepo = 1; iepo <= _epo_num; iepo++)
		{
			search_flg.push_back(2);
			if (!_isTrue(_post_data[sat].iflag[iepo], "OK")) continue;
			search_flg[iepo] = 1;
			if (!_isTrue(_post_data[sat].iflag[iepo], "GFJUMP")) search_flg[iepo] = 0;   // Not LGJUMP
		}

		// loop epochs
		int ibeg, iend, kobs, i;        // first/last usable epoch, useful epo_num, loop useful epochs
		vector<double> inp_lg, inp_ti, inp_res;
		vector<int> inp_flg;
		for (iepo = 1; iepo <= _epo_num; iepo++)
		{
			if (search_flg[iepo] == 2) continue;     // exclude bad point

			// search backward 5 epochs flagged as ok
			if (_searchEpoch(search_flg, iepo, 0, "-", ibeg) == 0)
			{
				_post_data[sat].iflag[iepo] = _setFlag(_post_data[sat].iflag[iepo], "GFJUMP");
				continue;
			}

			// search forwards, be aware that we have at least one usable point iepo.
			_searchEpoch(search_flg, iepo, _epo_num + 1, "+", iend);

			// total epochs
			kobs = iend - ibeg + 1;
			// insert relavant elements
			inp_lg.clear(); inp_ti.clear(); inp_flg.clear(); inp_res.clear();
			for (i = ibeg; i <= iend; i++)
			{
				inp_lg.push_back(_post_data[sat].obs[i][1]);
				inp_res.push_back(_post_data[sat].obs[i][5]);
				inp_flg.push_back(search_flg[i]);
				inp_ti.push_back(_post_time[i].diff(_beg_time));
			}

			// check for jump within the selected piece
			double rms = 0.0;
			bool MaybeNotjump = _checkForJump("LGC", inp_ti, inp_lg, inp_flg, inp_res, rms);
			// update search_flg
			for (i = ibeg; i <= iend; i++)
			{
				search_flg[i] = inp_flg[i - ibeg];
			}
			// true, valid check , may not exsit LG jump
			if (MaybeNotjump)
			{
				// update turbo data
				_post_data[sat].obs[iepo][5] = inp_res[iepo - ibeg];
				_post_data[sat].obs[iepo][6] = rms;
				// residual control only
				if (double_eq(_gf_rms_limit, 0.0) && fabs(_post_data[sat].obs[iepo][5]) <= _gf_limit) {
					continue;
				}
				// with rms control
				if (!double_eq(_gf_rms_limit, 0.0) && rms <= _gf_rms_limit &&
					fabs(_post_data[sat].obs[iepo][5]) <= (3 * rms > _gf_limit ? 3 * rms : _gf_limit)) {
					continue;
				}
			}

			_post_data[sat].iflag[iepo] = _setFlag(_post_data[sat].iflag[iepo], "GFJUMP");
		}
	}

	bool t_gexeturboedit::_checkForJump(string strflag, const vector<double>& ti,
		const vector<double>& values, vector<int>& flag, vector<double>& res, double& rms)
	{
		// with coefficients of polynomial, degree = 2
		vector<double> coeff;
		coeff.push_back(0.0); coeff.push_back(0.0);
		coeff.push_back(0.0); coeff.push_back(0.0);


		vector<int> ipt;           // with information about used epoch, default is -1
		int kobs = ti.size();      // number of obs
		int epo0 = round(fmod(ti[0], 86400.0) / _interval) + 1;  // sod in a single day
		double pre_rms = 1.e10;    // previous rms
		int degree = 0, nflag = 0, iter = 0, iepo = 0;

		bool iteration = true;
		while (iteration)
		{
			if (!_polynomialFit(degree, ti, values, flag, coeff, res, ipt, rms))
			{
				_log->comment(1, _crt_rec, "Ploydf error!");
				return false;
			}

			// detecting flagged number
			nflag = 0;
			for (iepo = 0; iepo < kobs; iepo++)
			{
				if (ipt[iepo] == -1) continue;

				if (fabs(res[iepo]) > (3 * rms > 0.30 ? 3 * rms : 0.30))// residuals with 3 times RMS
				{
					if (flag[iepo] == 0)
					{
						_log->comment(1, _crt_rec, strflag + " flagged(b) : at Epoch " + int2str(epo0 + iepo)
							+ "   " + dbl2str(ti[iepo]) + "  [max(0.30,3rms)] -- 3rms=" + dbl2str(3 * rms)
							+ " -- residuals=" + dbl2str(res[iepo]));
						nflag++;
						flag[iepo] = 1;
					}
				}
			}
			// Number of iterations
			iter++;
			iteration = false;

			if (nflag > 0)
			{
				iteration = true;
			}
			else if (rms > 0.3)   // rms threshold (unit error)
			{
				iteration = true;
				if (rms / pre_rms > 0.95)  // have not obvious change, add degree
				{
					degree++;
					iter = 1;
				}
				if (degree > 2) iteration = false;
			}

			pre_rms = rms;
			if (iteration && iter <= 2) iteration = true;
		}

		return true;
	}

	bool t_gexeturboedit::_polynomialFit(int degree, const vector<double>& X, const vector<double>& Y,
		const vector<int>& L, vector<double>& C, vector<double>& V, vector<int>& ipt, double& rms)
	{
		ipt.clear();        // ipt  re-initialization for next time
		// normalisation of independent argument
		double Xmax = -1.e20;
		double Xmin = 1.e20;
		int N = X.size();
		int m, i, j, k;
		// searching the Max.  and  min.  Independent argument values
		for (i = 0; i < N; i++)
		{
			ipt.push_back(-1);    // flag used epoch
			// L: 0: best - NO LG/LC/.. JUMP data; 1: Followed by best - OK data[Jump + GOOD]; 2: Not OK
			// Preset Independent variable search range
			if (L[i] <= 1)        
			{
				Xmax = X[i] > Xmax ? X[i] : Xmax;
				Xmin = X[i] < Xmin ? X[i] : Xmin;
				V[i] = 0.0;
			}
		}
		// No way to search
		if (Xmax == Xmin) return false;

		double dx = (Xmax + Xmin) / 2.0;
		double ft = 2.0 / (Xmax - Xmin);

		// initial importent parameter
		Matrix Q; ColumnVector F; ColumnVector A;
		if (degree > 0)
		{
			Q.ReSize(degree, degree); Q = 0.0;
			F.ReSize(degree);         F = 0.0;
			A.ReSize(degree);         A = 0.0;
		}

		int neffective_obs = 0;   // number of obs

		// loop over all functional values	
		double Xnorm, Xnorj;
		for (j = N - 1; j >= 1; j--)   // From the last to the second
		{
			if (L[j] == 0)  // get first usable functional value, must best
			{
				// searching next usable functional value to form difference
				m = j - 1;
				while (m >= 0 && L[m] > 1)  m--;
				if (m != -1)     // get two obs
				{
					Xnorj = (X[j] - dx)*ft;
					Xnorm = (X[m] - dx)*ft;
					neffective_obs++;                           // Number of observations
					if (degree <= 0) continue;
					// computaion of a-matrix (row number j)
					for (k = 1; k <= degree; k++)
					{
						A(k) = pow(Xnorm, k) - pow(Xnorj, k);   // equal to coefficient matrix
					}
					// contribution to q - matrix
					Q += A * A.t();                             // equal to BTPB
					// contribution to right hand side
					for (i = 1; i <= degree; i++)
					{
						F(i) += A(i) * (Y[m] - Y[j]);           // equal to BTPL
					}
				}
			}
		}

		// not enough obs < number of parameters
		if (neffective_obs < degree) return false;

		// solution vector (polynomial coefficients) equal to inv(BTPB)BTPL
		// Check Matrix singularity
		if (degree > 0)
		{
			Matrix solution;
			if (double_eq(Q.maximum_absolute_value(), 0.0)) return false;
			
			try
			{
				solution << Q.i() * F;
			}
			catch (exception e)
			{
				return false;
			}

			for (i = 1; i <= degree; i++)
			{
				C[i - 1] = solution(i, 1);
			}
		}

		//  residuals, rms
		rms = 0.0;
		double Vmax = 0.0, ys = 0.0;
		for (j = N - 1; j >= 1; j--)
		{
			if (L[j] <= 1)  // same ways ..... like before
			{
				// searching next usable functional value to form difference
				m = j - 1;
				while (m >= 0 && L[m] > 1)  m--;
				if (m != -1)
				{
					Xnorj = (X[j] - dx)*ft;
					Xnorm = (X[m] - dx)*ft;

					ys = 0.0;
					if (degree > 0)
					{
						for (k = 1; k <= degree; k++)
						{
							ys += C[k - 1] * (pow(Xnorm, k) - pow(Xnorj, k));
						}
					}

					V[j] = ys - (Y[m] - Y[j]);  // get residuals

					if (L[j] == 0)    // calculate rms with best data
					{
						if (fabs(V[j]) > Vmax) Vmax = fabs(V[j]);
						rms += V[j] * V[j];
					}

					ipt[j] = m;
				}
			}
		}

		double rms1;
		if (neffective_obs - degree != 0 && rms > 0.0)   // judge Degrees of freedom to compute Unit weight error
		{
			Vmax = rms - Vmax * Vmax;
			rms = sqrt(rms / (neffective_obs - degree));   // Unit weight error
			if (neffective_obs - degree - 1 != 0 && Vmax > 0.0)
			{
				rms1 = sqrt(Vmax / (neffective_obs - degree - 1));
				if (rms / rms1 > 2.0) rms = rms1;
			}

			return true;
		}
		else
		{
			_log->comment(2, _crt_rec, "BAD RMS, " + dbl2str(rms) +
				"  number of effective obs and degree  " + int2str(neffective_obs)
				+ "  " + int2str(degree));
			return false;
		}
	}

	void t_gexeturboedit::_postPrintReport(string& decision)
	{
		map<int, map <string, int> > sat_state;  // five state/sat/number
		map<int, int> state_sum;                 // five sate/number
		map<string, bool> lfirst_sat;            // first_sats
		int useful_sat = 0;                      // used sat number
		int iepo;
		for (auto sat : _sat_list)
		{
			lfirst_sat[sat] = true;
			for (iepo = 1; iepo <= _epo_num; iepo++)
			{
				if (_isTrue(_post_data[sat].iflag[iepo], "NODATA")) continue;

				sat_state[1][sat] += 1;                            // 1: Total observatons     
				if (_isTrue(_post_data[sat].iflag[iepo], "NO4"))
				{
					sat_state[4][sat] += 1;                        // 4: Incomplete observations   
				}
				else if (_isTrue(_post_data[sat].iflag[iepo], "LOWELE"))
				{
					sat_state[5][sat] += 1;                        // 5: Low elevation observations             
				}
				else if (!_isTrue(_post_data[sat].iflag[iepo], "OK"))
				{
					sat_state[2][sat] += 1;                        // 2: No-brd+Bad-brd+Short+PC/LW/LG-Bad+Shadow            
				}
				else if (!_isTrue(_post_data[sat].iflag[iepo], "GOOD"))
				{
					if (lfirst_sat[sat])  lfirst_sat[sat] = false;
					else                  sat_state[3][sat] += 1;  // 3: OK-GOOD[Gap+LWJump+LGJump+LLI], means new ambiguities
				}
			}

			// how many observations for TOTAL,NO4,LOWELE,OK,GOOD
			if (sat_state[1][sat] == 0) continue;
			useful_sat++;
			state_sum[1] += sat_state[1][sat];
			state_sum[2] += sat_state[2][sat];
			state_sum[3] += sat_state[3][sat];
			state_sum[4] += sat_state[4][sat];
			state_sum[5] += sat_state[5][sat];
		}

		map<int, int> nobs_epo;    // nbos each epoch
		map<int, int> nepo_nsv;    // e.g 1 obs have 3 epochs, 2 obs have 10 epochs, ...


		int kepo4 = 0;
		int kepo4_nobs = 0;
		int nepo_input = 1;
		bool linput;
		for (int i = 1; i <= _epo_num; i++)
		{
			linput = false;
			nobs_epo[i] = 0;
			for (auto sat : _sat_list)
			{
				if (!_isTrue(_post_data[sat].iflag[i], "NODATA"))
				{
					linput = true;
				}

				if (_isTrue(_post_data[sat].iflag[i], "OK"))
				{
					nobs_epo[i] += 1;
				}
			}
			if (linput) nepo_input++;
			nepo_nsv[nobs_epo[i]] += 1;
			if (nobs_epo[i] >= 4)
			{
				kepo4++;
				kepo4_nobs += nobs_epo[i];
			}
		}

		double mean_namb_per_prn = 0.0;
		if (useful_sat != 0) mean_namb_per_prn = state_sum[3] / (useful_sat*1.0);

		double mean_nprn_per_epo = 0.0;
		if (kepo4 != 0) mean_nprn_per_epo = kepo4_nobs / (kepo4 * 1.0);

		double percent_kepo4 = 0.0;
		// i.e. useful_sat != 0
		if (state_sum[1] != 0) percent_kepo4 = kepo4 / (nepo_input * 1.0);


		decision = "GOOD";
		if (_check_statistics)
		{
			if (percent_kepo4 < _min_percent*1.e-2)
			{
				decision = "BAD too many removed";
			}
			else if (mean_namb_per_prn > _max_mean_namb)
			{
				decision = "BAD too many AMBs";
			}
			else if (mean_nprn_per_epo < _min_mean_nprn)
			{
				decision = "BAD too small mean nprn";
			}
		}


		_log->comment(1, _crt_rec, "+++  Summary  +++");
		_log->comment(1, _crt_rec, " Site and Evaluation : " + _crt_rec + "   " + decision);
		_log->comment(1, _crt_rec, " #EPOCH > 4_SAT      : " + dbl2str(percent_kepo4*100) + "%");
		_log->comment(1, _crt_rec, " MEAN_#AMB_PER_PRN   : " + dbl2str(mean_namb_per_prn));
		_log->comment(1, _crt_rec, " MEAN_#PRN_PER_EPOCH : " + dbl2str(mean_nprn_per_epo));
		_log->comment(1, _crt_rec, " NUMBER_OF_OBS       :    #TOTAL:   " + int2str(state_sum[1]) + 
			"   #REMOVED:   "+ int2str(state_sum[2]) + "   #LOWELEV:   " + int2str(state_sum[5]) +
			"   #AMB:   " + int2str(state_sum[3]));
		stringstream str_prn, str_total, str_removed, str_no4, str_lowelev, str_amb;
		str_prn     << " PRN      TOTAL";
		str_total   << " #TOTAL   " << setw(5) << state_sum[1];
		str_removed << " #REMOVED " << setw(5) << state_sum[2];
		str_no4     << " #NO_4TYP " << setw(5) << state_sum[4];
		str_lowelev << " #LOWELEV " << setw(5) << state_sum[5];
		str_amb     << " #AMB     " << setw(5) << state_sum[3];
		for (auto sat : _sat_list)	
		{
			if (sat_state[1][sat] == 0) continue;
			str_prn     << "  " + sat;
			str_total   << setw(5) << sat_state[1][sat];
			str_removed << setw(5) << sat_state[2][sat];
			str_no4     << setw(5) << sat_state[4][sat];
			str_lowelev << setw(5) << sat_state[5][sat];
			str_amb     << setw(5) << sat_state[3][sat];
		}
		_log->comment(1, _crt_rec, str_prn.str());
		_log->comment(1, _crt_rec, str_total.str());
		_log->comment(1, _crt_rec, str_removed.str());
		_log->comment(1, _crt_rec, str_no4.str());
		_log->comment(1, _crt_rec, str_lowelev.str());
		_log->comment(1, _crt_rec, str_amb.str());
		
		stringstream str_svobs, str_epoch;
		str_svobs << " #SV obs. ";
		str_epoch << " #epoch   ";
		for (int j = 0; j <= 36; j++)
		{
			str_svobs << setw(5) << j;
			str_epoch << setw(5) << nepo_nsv[j];
		}
		_log->comment(1, _crt_rec, str_svobs.str());
		_log->comment(1, _crt_rec, str_epoch.str());
		_log->comment(1, _crt_rec, "---  Summary  ---");
	}

	bool t_gexeturboedit::_postSetAmbflag()
	{	
		map<string, map<int, int>> local_flag;
		map<string, vector<int>> search_flag_data;
		map<int,int> namb_epo;
		int iepo, jepo, kepo;
		for (auto sat : _sat_list)
		{
			search_flag_data[sat].push_back(999);
			for (iepo = 1; iepo <= _epo_num; iepo++)
			{
				search_flag_data[sat].push_back(_post_data[sat].iflag[iepo]);				
				local_flag[sat][iepo] = 0;
				namb_epo[iepo] = 0;
			}
		}

		// Ambiguities infomation flagged at local_flag
		// local_flag[sat][iepo]==iepo BAD data
		// local_flag[sat][iepo] < 0   AMB data
		// local_flag[sat][iepo] >iepo DEL data
		ambflag_hd ambflag_head;
		ambflag_head.new_rm_obs = 0;
		ambflag_head.new_amb = 0;
		ambflag_head.avaiable_obs = 0;
		for (auto sat : _sat_list)
		{
			iepo = 0;
			while (iepo < _epo_num)
			{
				iepo++;
				// check obs or not
				if (_isTrue(_post_data[sat].iflag[iepo], "NODATA")) continue;
				// select bad point 
				if (!_isTrue(_post_data[sat].iflag[iepo], "OK"))
				{
					jepo = iepo;
					// Traverse from the current epoch to the end, if the data is OK, jump out
					while (jepo <= _epo_num && !_isTrue(_post_data[sat].iflag[jepo], "OK"))
					{
						// Delete an observation if it exists
						if (!_isTrue(_post_data[sat].iflag[jepo], "NODATA"))
						{
							ambflag_head.new_rm_obs += 1;
							local_flag[sat][iepo] = jepo;   // iepo->jepo  bad duration
						}
						jepo++;
					}
					iepo = jepo - 1;  // Minus 1 because jepo adds 1
				}
				else                  // OK data
				{
					if (!_isTrue(_post_data[sat].iflag[iepo], "GOOD"))  // OK-GOOD：insert AMB because of GAP/JUMP
					{
						//  find the next ambiguity or end of the data set
						jepo = _findFlag(iepo + 1, _epo_num, search_flag_data[sat], "AMB", "+");
						if (jepo == -999)  jepo = _epo_num + 1;

						// find backwards the last usable epoch
						jepo = _findFlag(iepo + 1, jepo - 1, search_flag_data[sat], "GOOD", "-");
						if (jepo == -999) jepo = iepo;

						// accumulate namb per epoch
						for (kepo = iepo; kepo <= jepo; kepo++)
						{
							namb_epo[kepo] += 1;
						}
						// new amb num
						ambflag_head.new_amb += 1;

						local_flag[sat][iepo] = -jepo;   // New amb 是负数
					}
					ambflag_head.avaiable_obs += 1;
				}
			}
		}

		// update ambflag header
		string sat = *(_sat_list.begin());
		ambflag_head.intv = _interval;
		ambflag_head.duration = _end_time.diff(_beg_time);
		ambflag_head.beg_mjd = _beg_time.mjd();
		ambflag_head.beg_sod = _beg_time.sod() + _beg_time.dsec();
		int max_nobs = -1;
		for (iepo = 1; iepo <= _epo_num; iepo++)
		{
			max_nobs = max_nobs > namb_epo[iepo] ? max_nobs : namb_epo[iepo];
		}
		ambflag_head.max_amb_1epo = max_nobs;

		// add ambflag_head
		_allambflag->addAmbFlagHead(_crt_rec, ambflag_head);
		_cycleslip[_crt_rec].back()->set_active_amb(_crt_rec, ambflag_head.max_amb_1epo);
		// add ambflag data
		if (!_postSetAmbflagData(local_flag)) return false;

		else return true;
	}

	bool t_gexeturboedit::_postSetAmbflagData(map<string, map<int, int>>& sat_epo_beg_end)
	{
		ambflag_data AmbFlag_data;
		t_gtime beg, end;
		for (int iepo = 1; iepo <= _epo_num; iepo++)
		{
			for (std::string sat : _sat_list)
			{	
				AmbFlag_data.reason = "RN_     ";
				if (sat_epo_beg_end[sat][iepo] == iepo)  // End epoch is stored
				{
					AmbFlag_data.identify = "BAD";
					AmbFlag_data.iflag = "5";
					AmbFlag_data.reason = _strFlagDescription(sat, "BAD", iepo);
				}
				else if (sat_epo_beg_end[sat][iepo] > iepo)
				{
					AmbFlag_data.identify = "DEL";
					AmbFlag_data.iflag = "6";
					AmbFlag_data.reason = _strFlagDescription(sat, "DEL", iepo);
				}
				else if (sat_epo_beg_end[sat][iepo] < 0)
				{
					AmbFlag_data.identify = "AMB";
					AmbFlag_data.iflag = "1";
					AmbFlag_data.reason = _strFlagDescription(sat, "AMB", iepo);
					sat_epo_beg_end[sat][iepo] = -sat_epo_beg_end[sat][iepo];
				}

				if (AmbFlag_data.reason == "RN_Unkown_Amb" || AmbFlag_data.reason == "RN_Unkown_Bad")
				{
					_log->comment(1, _crt_rec, AmbFlag_data.reason + " at Epoch  " + int2str(iepo) +
							+ "  " + sat + ", iflag= " + int2str(_post_data[sat].iflag[iepo]));
					return false;
				}

				if (AmbFlag_data.reason != "RN_     ")
				{
					AmbFlag_data.beg_epo = iepo;
					AmbFlag_data.end_epo = sat_epo_beg_end[sat][iepo];
					AmbFlag_data.C1 = 0.0;
					AmbFlag_data.C2 = 0.0;
					_allambflag->addAmbFlagData(_crt_rec, sat, AmbFlag_data);
					beg = _beg_time + (AmbFlag_data.beg_epo - 1) * _interval;
					end = _beg_time + (AmbFlag_data.end_epo - 1) * _interval;
					_cycleslip[_crt_rec].back()->add_ambflag(_crt_rec, sat, AmbFlag_data.identify, beg, end);
				}
			}
		}
		return true;
	}

	string t_gexeturboedit::_strFlagDescription(string sat, string str3, int epo)
	{
		string str_back = "     ";
		if (str3 == "BAD")
		{
			if      (_isTrue(_post_data[sat].iflag[epo], "LOWELE"))     str_back = "lowelevation";
			else if (_isTrue(_post_data[sat].iflag[epo], "SHRT"))       str_back = "shortpiece";
			else if (_isTrue(_post_data[sat].iflag[epo], "MWBAD"))      str_back = "badwidelane";
			else if (_isTrue(_post_data[sat].iflag[epo], "GFBAD"))      str_back = "badionosphere";
			else if (_isTrue(_post_data[sat].iflag[epo], "NO4"))        str_back = "lessthan4obs";
			else if (_isTrue(_post_data[sat].iflag[epo], "PCBAD"))      str_back = "badrange";		
			else if (_isTrue(_post_data[sat].iflag[epo], "ASHADOW"))    str_back = "aftershadow";
			else if (_isTrue(_post_data[sat].iflag[epo], "NOBRD"))      str_back = "nobrd";
			else                                                        str_back = "Unkown_Bad";
		}
		else if (str3 == "AMB")
		{
			if      (_isTrue(_post_data[sat].iflag[epo], "GAP"))    str_back = "biggap";
			else if (_isTrue(_post_data[sat].iflag[epo], "MWJUMP")) str_back = "widelanejump";
			else if (_isTrue(_post_data[sat].iflag[epo], "GFJUMP")) str_back = "iono.jump";
			else if (_isTrue(_post_data[sat].iflag[epo], "LLI"))    str_back = "flaginrinex";
			else                                                    str_back = "Unkown_Amb";
		}
		else   // DEL
		{
			str_back = "mixed";
		}
		str_back = "RN_" + str_back;
		return str_back;
	}

	bool t_gexeturboedit::_postEncoderProduct(string logsuffix)
	{
		if (_allambflag->getAllAmbFlag().find(_crt_rec) == _allambflag->getAllAmbFlag().end())
		{
			_log->comment(1, _crt_rec, "No ambflag info : " + _crt_rec);
			return false;
		}
		
		string strambflagpath = dynamic_cast<t_gsetout*>(_set)->outputs("ambflag_dir");
		if (strambflagpath.empty())   strambflagpath = "file://log_tb";
		if (ACCESS(strambflagpath.substr(7).c_str(), 0) != 0)
		{
			MKDIR(strambflagpath.substr(7).c_str());
		}

		string site = _crt_rec;
		transform(site.begin(), site.end(), site.begin(), ::tolower);
		stringstream str;
		str << strambflagpath << PATH_SEPARATOR << site << setw(3) << setfill('0') << _beg_time.doy() << "0."
			<< setw(2) << setfill('0') << _beg_time.yr() << "o" + logsuffix;

		t_gio* gout = new t_gfile;
		//set I/O
		gout->path(str.str()); gout->glog(_log);
		//set coder
		t_gcoder* gcoder = new t_ambflag(nullptr, "", 4096);
		gcoder->clear(); gcoder->glog(_log); gcoder->path(str.str());
		gcoder->add_data("IDX", _allambflag);
		gout->coder(gcoder);
		// write
		gout->run_write();
		delete gout; gout = nullptr;
		delete gcoder; gcoder = nullptr;
		return true;
	}

	int t_gexeturboedit::_str2DataFlag(string data_flag)
	{
		string tmp = data_flag;
		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
		if      (tmp == "GOOD")    return DATAFLAG::GOOD;
		else if (tmp == "GFJUMP")  return DATAFLAG::GFJUMP;
		else if (tmp == "MWJUMP")  return DATAFLAG::MWJUMP;
		else if (tmp == "GAP")     return DATAFLAG::GAP;
		else if (tmp == "LLI")     return DATAFLAG::LLI;
		else if (tmp == "PCBAD")   return DATAFLAG::PCBAD;
		else if (tmp == "NOBRD")   return DATAFLAG::NOBRD;
		else if (tmp == "GFBAD")   return DATAFLAG::GFBAD;
		else if (tmp == "MWBAD")   return DATAFLAG::MWBAD;
		else if (tmp == "ASHADOW") return DATAFLAG::ASHADOW;
		else if (tmp == "SHRT")    return DATAFLAG::SHRT;
		else if (tmp == "LOWELE")  return DATAFLAG::LOWELE;
		else if (tmp == "NO4")     return DATAFLAG::NO4;
		else if (tmp == "NODATA")  return DATAFLAG::NODATA;
		else                       return -1;
	}

	int t_gexeturboedit::_setFlag(int iflag, string type)
	{
		string trim_type = trim(type);
		transform(trim_type.begin(), trim_type.end(), trim_type.begin(), ::toupper);

		int set_flag = iflag;

		if (trim_type == "OK")
		{
			for (int i = 16; i <= 31; i++)
			{
				set_flag = set_flag & ~(1 << (i + 1 - 1));  // 把右数第i + 1位变成0
			}
		}
		else if (trim_type == "GOOD")
		{
			for (int i = 0; i <= 31; i++)
			{
				set_flag = set_flag & ~(1 << (i + 1 - 1));  // 把右数第i + 1位变成0
			}
		}
		else
		{
			set_flag = set_flag | (1 << (_str2DataFlag(trim_type) + 1 - 1)); // 把右数第str2DataFlag(trim_type) + 1位变成1
		}
		return set_flag;
	}

	int t_gexeturboedit::_findFlag(int istart, int istop, const vector<int>& flg, string flag_type, string direction)
	{
		int  i1, i2, istep;

		if (direction == "+") {
			i1 = istart;
			i2 = istop;
			istep = 1;
		}
		else {
			i1 = istop;
			i2 = istart;
			istep = -1;
		}

		int iepo = -999;   // if can't find, return the -1
		for (int i = i1; i != i2 + istep; i += istep) {
			if (!_isTrue(flg[i], flag_type)) {
				continue;
			}
			iepo = i;
			break;
		}

		return iepo;
	}

	bool t_gexeturboedit::_isTrue(int iflag, string type)
	{
		bool istrue = false;
		string trim_type = trim(type);
		transform(trim_type.begin(), trim_type.end(), trim_type.begin(), ::toupper);

		if (trim_type == "OK")
		{
			// note: z'FFFF0000' = -65536
			istrue = (iflag & -65536) == 0 ? true : false;
		}
		else if (trim_type == "GOOD")
		{
			istrue = iflag == 0 ? true : false;
		}
		else if (trim_type == "AMB")
		{
			if ((iflag & -65536) == 0 && iflag != 0)
			{
				istrue = true;
			}
		}
		else
		{
			int tmp = iflag | (1 << (_str2DataFlag(trim_type) + 1 - 1));
			istrue = tmp == iflag ? true : false;
		}

		return istrue;
	}

	int t_gexeturboedit::_searchEpoch(const vector<int>& flg, int crt_epo, int limit_epo, string mode, int& use_epo)
	{
		int jepo, prefix;
		int max_epo_num;
		if (mode == "-")  // backward "-"
		{
			jepo = crt_epo - 1;
			max_epo_num = 30;
			prefix = -1;
		}
		if (mode == "+")  // forward "-"
		{
			jepo = crt_epo;
			max_epo_num = 29;
			prefix = 1;
		}

		int npoint = 0;
		use_epo = crt_epo;
		while (true)
		{
			if (npoint >= 5) break;
			if (abs(crt_epo - jepo) > max_epo_num) break;
			if (jepo == limit_epo) break;
			if (flg[jepo] < 2) {
				npoint++;
				use_epo = jepo;
			}
			jepo = jepo + prefix;
		}

		return npoint;
	}

}


