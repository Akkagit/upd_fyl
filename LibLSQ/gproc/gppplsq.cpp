// history : 
//	jdhuang 2019-05-12 add the iono pars init function in _InitLsq

#include "gproc/gppplsq.h"
#include "gproc/gspplsq.h"
#include "gproc/gupdateparIF.h"
#include "gproc/gupdateparALL.h"
#include "gset/gsetpar.h"
#include "gutils/ginfolog.h"
#include "gutils/gturboedit.h"
namespace great
{
	t_gppplsq::t_gppplsq(string site, t_gsetbase* set, t_gallproc* data, t_glog* log ) :
		t_glsqproc(set, data, log),
		t_gxml("kml"),
		_kml(false),
		_fltfile(0),
		_enufile(0),
		_isFirstFix(true),
		_recclk_file("recclk"),   //for iterater
		_isClient(false),
		_quality_control(_gset, nullptr)
	{
		// set site
		_crt_rec = site;
		_rec_list.clear();
		_rec_list.insert(site);

		// init ppp settings
		_init_ppp_setings();
		_gallbias = dynamic_cast<t_gallbias*>((*data)[t_gdata::ALLBIAS]);
		// out file pointer
		_setOut();
	}

	t_gppplsq::~t_gppplsq()
	{
		// delete _result
		if (_fltfile)
		{
			if(_fltfile->is_open()) _fltfile->close();
			delete _fltfile;
			_fltfile = nullptr;
		}

		// delete _enufile
		if (_enufile)
		{
			if (_enufile->is_open()) _enufile->close();
			delete _enufile;
			_enufile = nullptr;
		}

		if (_ambfix)
		{
			delete _ambfix;
			_ambfix = nullptr;
		}

		if (_spplsq)
		{
			delete _spplsq;
			_spplsq = nullptr;
		}

		// recover for every file
		if (_gall_rcv)
		{
			delete _gall_rcv;
			_gall_rcv = nullptr;
		}

		// delete tropfile
		if (_tropfile)
		{
			if (_tropfile->is_open()) _tropfile->close();
			delete _tropfile;
			_tropfile = nullptr;
		}

		if (_kml)
		{
			this->write(_kml_name);
		}

#ifdef DEBUG_PPPLSQIF
		if (!_result) _writeLogInfo(_log, 0, "DEBUG : the result file pointer is delete.");
		if (!_enufile) _writeLogInfo(_log, 0, "DEBUG : the enu file pointer is delete.");
		if (!_ambfix) _writeLogInfo(_log, 0, "DEBUG : the amb file pointer is delete.");
#endif // DEBUG_PPPLSQIF
	}

	bool t_gppplsq::GenerateProduct()
	{
		bool prod_valid = _initLsqProdData(nullptr);
		if (!prod_valid) return false;
		return true;
	}

	bool t_gppplsq::_init_ppp_setings()
	{
		///< set reference coordinate
		if (_crd_est == CONSTRPAR::SIMU_KIN || _crd_est == CONSTRPAR::EST)
		{
			xyz_standard = dynamic_cast<t_gsetrec*>(_gset)->get_crd_xyz(_crt_rec);
			neu_sum_nepo = make_pair(0, t_gtriple(0.0, 0.0, 0.0));
		}

		///< Judge whether use simulation data (PPPLSQ/PPPLEO)
		_simulation = dynamic_cast<t_gsetproc*>(_gset)->simulation();

		///< PPPLEO settings
		if (_simulation)
		{
			// ZTD STO
			_recclk_threshold = 1e-6;
			_bias_model = shared_ptr<t_gbiasmodel>(new t_simuleobias(_gall_proc, _glog, _gset));
			_base_model = shared_ptr<t_gbasemodel>(new t_gcombIF(_gset, _glog, _bias_model, _gall_proc));
		}

		if (_crd_est == KIN)
		{
			_recclk_threshold = 1e-6;
		}

		///< ambiguity fixing settings
		if (_ambfix)
		{
			delete _ambfix;
			_ambfix = nullptr;
		}

		if (_fix_mode != t_gsetamb::FIX_MODE::NO)
		{
			_ambfix = new t_gambiguity(_crt_rec, _gset);
			_ambfix->setActiveAmb(_sat_list.size());
		}

		///< recover info
		if (_gall_rcv)
		{
			delete _gall_rcv;
			_gall_rcv = nullptr;
		}
		_gall_rcv = new t_gallrecover();
		_gall_rcv->glog(_glog);

		///< real time
		_wl_Upd_time = WL_IDENTIFY;
		_ewl_Upd_time = EWL_IDENTIFY;
		_ewl24_Upd_time = EWL24_IDENTIFY;
		_ewl25_Upd_time = EWL25_IDENTIFY;

		return true;
	}

	bool t_gppplsq::_generateAmbupdIF(t_gallambupd* allambupd)
	{
		if (_lsq->mode() == gnut::EPO) return false;
		map<string, map<t_gtime, map<int, double> > > mw_used;
		list<pair<double, double>> wx;
		t_gtime it;
		double mean, sigma, mean_sig;
		// get wl ambiguities and sigma <remove outlier>
		auto iter_sat = _sat_list.begin();
		for (; iter_sat != _sat_list.end(); iter_sat++)
		{
			if (_AmbValidTime.find(*iter_sat) == _AmbValidTime.end()) continue;
			auto iter_time = _AmbValidTime[*iter_sat].begin();
			for (; iter_time != _AmbValidTime[*iter_sat].end(); iter_time++)
			{
				it = iter_time->first;
				while (it <= iter_time->second) {
					if (!double_eq(_MW[it][*iter_sat][4], 0.0)) {
						wx.push_back(make_pair(_MW[it][*iter_sat][4], 1.0));
					}
					it = it + _obs_intv;
				}
				mean = 0.0, sigma = 0.0, mean_sig = 0.0;
				getMeanWgt(true, wx, mean, sigma, mean_sig);
				if (double_eq(mean_sig, 0.0)) mean_sig = 999.0;

				it = iter_time->first;
				while (it <= iter_time->second) {
					mw_used[*iter_sat][it][1] = mean;
					mw_used[*iter_sat][it][2] = mean_sig;
					it = it + _obs_intv;
				}
				wx.clear();
			}
		}

		// get IF ambiguities [cycle]
		t_gambupd gambupd;
		t_ambupd_data ambupd_data;
		// get all amb par
		vector<t_grecover_par> par_amb_IF = _gall_rcv->get_recover_par(t_gpar::t_type::AMB_IF);
		int parnum = _lsq->_x_solve.parNumber();
		for (int i = 0; i < parnum; i++)
		{
			if (_lsq->_x_solve[i].parType != t_gpar::t_type::AMB_IF) continue;
			t_grecover_par rcv_par(_lsq->_x_solve[i], _lsq->dx()(i + 1));			
			par_amb_IF.push_back(rcv_par);
		}

		for (auto iter = par_amb_IF.begin(); iter != par_amb_IF.end(); iter++)
		{
			if (iter->par.site != _crt_rec) continue;
			double value = iter->par.value() + iter->correct_value;
			string sat = iter->par.prn;
			if (!double_eq(value, 0.0) && mw_used.find(sat) != mw_used.end())
			{
				it = iter->par.beg;
				while (it <= iter->par.end)
				{
					if (mw_used[sat].find(it) != mw_used[sat].end())
					{
						ambupd_data.amb_if = value;
						ambupd_data.amb_wl = mw_used[sat][it][1];
						ambupd_data.amb_wl_sigma = mw_used[sat][it][2];
						gambupd.addAmbUpd(it, sat, ambupd_data);
					}
					it = it + _obs_intv;
				}
			}
		}
		allambupd->addAmbUpd(_crt_rec, gambupd);
		return true;
	}

	bool t_gppplsq::_generateAmbupdALL(t_gallambupd* allambupd)
	{
		if (_lsq->mode() == gnut::EPO) return false;

		t_gtime it;
		// get UDUC ambiguities [meter]
		t_gambupd gambupd;
		t_ambupd_data ambupd_data;
		// get all amb par
		vector<t_grecover_par> par_amb_ALL, par_amb_tmp;
		par_amb_tmp = _gall_rcv->get_recover_par(t_gpar::AMB_L1);
		par_amb_ALL.insert(par_amb_ALL.end(), par_amb_tmp.begin(), par_amb_tmp.end());
		par_amb_tmp = _gall_rcv->get_recover_par(t_gpar::AMB_L2);
		par_amb_ALL.insert(par_amb_ALL.end(), par_amb_tmp.begin(), par_amb_tmp.end());
		par_amb_tmp = _gall_rcv->get_recover_par(t_gpar::AMB_L3);
		par_amb_ALL.insert(par_amb_ALL.end(), par_amb_tmp.begin(), par_amb_tmp.end());
		par_amb_tmp = _gall_rcv->get_recover_par(t_gpar::AMB_L4);
		par_amb_ALL.insert(par_amb_ALL.end(), par_amb_tmp.begin(), par_amb_tmp.end());
		par_amb_tmp = _gall_rcv->get_recover_par(t_gpar::AMB_L5);
		par_amb_ALL.insert(par_amb_ALL.end(), par_amb_tmp.begin(), par_amb_tmp.end());

		int parnum = _lsq->_x_solve.parNumber();
		for (int i = 0; i < parnum; i++)
		{
			if (_lsq->_x_solve[i].parType == t_gpar::AMB_L1
				|| _lsq->_x_solve[i].parType == t_gpar::AMB_L2 || _lsq->_x_solve[i].parType == t_gpar::AMB_L3
				|| _lsq->_x_solve[i].parType == t_gpar::AMB_L4 || _lsq->_x_solve[i].parType == t_gpar::AMB_L5)
			{
				t_grecover_par rcv_par(_lsq->_x_solve[i], _lsq->dx()(i + 1));
				par_amb_ALL.push_back(rcv_par);
			}
		}

		for (auto iter = par_amb_ALL.begin(); iter != par_amb_ALL.end(); iter++)
		{
			if (iter->par.site != _crt_rec) continue;
			double value = iter->par.value() + iter->correct_value;
			string sat = iter->par.prn;
			if (!double_eq(value, 0.0))
			{
				it = iter->par.beg;
				while (it <= iter->par.end)
				{
					//amb_if: ambiguity ionosphere-free value meter  /N1 in UCUD mode (m)
					//amb_wl: ambiguity Wide-Lane value       cycle  /N2 in UCUD mode (m)
					if (iter->par.parType == t_gpar::AMB_L1) ambupd_data.amb_if = value;
					if (iter->par.parType == t_gpar::AMB_L2) ambupd_data.amb_wl = value;
					if (iter->par.parType == t_gpar::AMB_L3) ambupd_data.amb_l3 = value;
					if (iter->par.parType == t_gpar::AMB_L4) ambupd_data.amb_l4 = value;
					if (iter->par.parType == t_gpar::AMB_L5) ambupd_data.amb_l5 = value;
					ambupd_data.amb_wl_sigma = 0.0;
					gambupd.addAmbUpd(it, sat, ambupd_data);
					it = it + _obs_intv;
				}

			}
			if (iter->par.parType == t_gpar::AMB_L1) ambupd_data.amb_if = 0;
			if (iter->par.parType == t_gpar::AMB_L2) ambupd_data.amb_wl = 0;
			if (iter->par.parType == t_gpar::AMB_L3) ambupd_data.amb_l3 = 0;
			if (iter->par.parType == t_gpar::AMB_L4) ambupd_data.amb_l4 = 0;
			if (iter->par.parType == t_gpar::AMB_L5) ambupd_data.amb_l5 = 0;
		}
		allambupd->addAmbUpd(_crt_rec, gambupd);
		return true;
	}

	bool t_gppplsq::InitProc(t_gallproc* data,const t_gtime& beg, const t_gtime& end)
	{
		if (!data ||
			beg > end)
		{
			write_log_info(_glog, 0, "ERROR : beg > end or no proc data");
			return false;
		}

		// get the time
		_beg_time = (_beg_time < beg) ? beg : _beg_time;
		_end_time = (_end_time > end) ? end : _end_time;
		_crt_time = beg;


		bool init_data = _initLsqProcData(data);
		if (!init_data)
		{
			_glog->comment(0, "ERROR : can not init data");
			return false;
		}

		bool init_pars = _initLsqProcPars(_lsq);
		if (!init_pars)
		{
			_glog->comment(0, "ERROR : can not init pars");
			return false;
		}

		// quality control
		_quality_control.setNav(_gall_nav);

		return true;
	}

	int t_gppplsq::ProcessOneEpoch(const t_gtime& now, vector<t_gsatdata>* obs)
	{
		_crt_time = now;
		vector<t_gsatdata> crt_obs;
		

		if (obs) crt_obs = *obs;
		else crt_obs = _gall_obs->obs(_rec_list, now);//maybe bug
		
		_crt_obs_list = crt_obs;
		
		if (crt_obs.size() == 0) {
			if (_glog && _glog->verb() >= 1) {
				_glog->comment(1, "gppplsq", _crt_rec + now.str_ymdhms(" no observation found at epoch: "));
			}
			return -1;
		}

		/*set glonass frequency number and correct dcb*/
		_single_freq_sats.clear();
		for (auto& iter : crt_obs)
		{
			if (iter.sys() == "R" && iter.channel() >= DEF_CHANNEL) iter.channel(_glofrq_num[iter.sat()]);
			iter.apply_bias(*_gallbias);
			if (_obs_mode == RAW_MIX) {
				set<GOBSBAND> band_avail = iter.band_avail(true);
				if (band_avail.size() < 2) _single_freq_sats.insert(iter.sat());
			}
		}

		/*set kinematic coordinate*/
		if (_crd_est == KIN)
		{
			if (_spplsq)
			{
				_spplsq->processBatch(now, now);
				//_spplsq->ProcessOneEpoch(now, &crt_obs);
				_kin_xyz[_crt_rec][now] = _spplsq->getCrd(now);
				cout << _crt_time.sod() << "  " << fixed << setw(16) << setprecision(4) << _kin_xyz[_crt_rec][now][0]
					 << setw(16) << setprecision(4) << _kin_xyz[_crt_rec][now][1]
					 << setw(16) << setprecision(4) << _kin_xyz[_crt_rec][now][2] << endl;
			}
			_quality_control.processOneEpoch(now, _crt_rec, _kin_xyz[_crt_rec][now], crt_obs);
		}
		else
		{
			_quality_control.processOneEpoch(now, _crt_rec, xyz_standard, crt_obs);
		}

		/*process one epoch*/
		bool epoch_valid = _processOneEpoch(now, crt_obs);
		if (!epoch_valid)
		{
			write_log_info(_glog, 0, "Epoch process failed" + now.str_ymdhms());
			return -1;
		}

		if (_lsq->mode() == gnut::EPO && now < _end_time)
		{
			int solve_valid = t_gppplsq::_SolveEpoch();
			if (solve_valid == -1)
			{
				write_log_info(_glog, 0, "PPP Falied!" + now.str_ymdhms());
				return -1;
			}
			else if (solve_valid == 0)
			{
				write_log_info(_glog, 0, "PPP-Float" + now.str_ymdhms());
				return 0;
			}
			else if (solve_valid == 1)
			{
				write_log_info(_glog, 0, "PPP-AR" + now.str_ymdhms());
				return 1;
			}
		}
		return 1;
	}

	bool t_gppplsq::ExtractResfile()
	{
		try
		{
			_lsq->solve_NEQ();
		}
		catch (exception e)
		{
			write_log_info(_glog, 1, e.what(), "Solve Equation Fail!");
			return false;
		}
		catch (NPDException e)
		{
			writeLogInfo(_glog, 0, "t_gppplsq::processBatch ", "Can't Solve NEQ! Too Few Observations!");
			return false;
		}

		/* recover equation */
		try
		{
			_lsq->get_result_parameter(*_gall_rcv);
			_lsq->recover_parameters(*_gall_rcv);
			_extract_resfile(*_gall_rcv, _crt_rec + "_resfile");
			_extractTropDelay(*_lsq);
			_extract_clkfile(*_gall_rcv, t_gallprec::clk_type::AR);
		}
		catch (exception e)
		{
			write_log_info(_glog, 1, e.what(), "Recover Parameter Fail!");
			return false;
		}
		return true;
	}

	void t_gppplsq::resetLsqUpdate(t_gsetbase* set, t_gallproc* data, t_glog* log)
	{
		_init_lsq_settings(set, data, log);
	}

	t_gallpar& t_gppplsq::getAllPar()
	{
		/*return _lsq->_x_solve;*/
		if (_isClient) _ambfix->getFinalParams();
		else return _ambfix->updateUCAmb();
	}

	void t_gppplsq::setWLParams(t_gallpar& params)
	{
		_wl_params = params; 
	}

	bool t_gppplsq::ProcessBatch(t_gallproc * data, const t_gtime & beg, const t_gtime & end)
	{
		if (!InitProc(data, beg, end)) return false;

		while (_crt_time <= _end_time)
		{
			cerr << _crt_time.str_ymdhms() << ":" << "\r";
			write_log_info(_glog, 1, "t_gppplsq::processBatch", _crt_time.str_ymdhms("Processing epoch: "));

			ProcessOneEpoch(_crt_time);
			_crt_time = _crt_time + _obs_intv;
		}

		if (_lsq->mode() == gnut::EPO && _crd_est != KIN && _enufile)
		{
			neu_sum_nepo.second /= neu_sum_nepo.first;
			ostringstream enu_info;
			enu_info << " RMS = " << fixed << setw(12) << setprecision(4) << sqrt(neu_sum_nepo.second.crd(1))
				<< setw(12) << setprecision(4) << sqrt(neu_sum_nepo.second.crd(0))
				<< setw(12) << setprecision(4) << sqrt(neu_sum_nepo.second.crd(2)) << endl;
			_enufile->write(enu_info.str().c_str(), enu_info.str().size());
			_enufile->flush();
		}

		if (!ExtractResfile())return false;

		return true;
	}

	bool t_gppplsq::_initLsqProdData(t_gallprod * data)
	{
		t_gallambupd* allambupd = new t_gallambupd();

		if (_obs_mode == OBSCOMBIN::IONO_FREE && !_generateAmbupdIF(allambupd) ||
			_obs_mode == OBSCOMBIN::RAW_ALL && !_generateAmbupdALL(allambupd))
		{
			delete allambupd; allambupd = NULL;
			return false;
		}

		_extract_ambupd(allambupd);
		// Encoder ambupd files
		//string path = "file://" + _crt_rec + "_ambupd_" + _beg_time.str_yyyydoy();

		//t_gfile* gout = new t_gfile;
		//gout->glog(_glog);
		//gout->path(path);

		//t_gcoder* gcoder = new t_ambupd(_gset, "", 4096);
		//gcoder->clear();
		//gcoder->glog(_glog);
		//gcoder->path(path);
		//gcoder->add_data("IDX", allambupd);

		//gout->coder(gcoder);

		//_glog->comment(1, "t_gppplsqALL::_initLsqProdData", " Writing AMBUPD file : " + path);
		//gout->run_write();
		//delete gout; gout = NULL;
		//delete gcoder; gcoder = NULL;
		//delete allambupd; allambupd = NULL;

		return true;
	}
	
	bool t_gppplsq::_initLsqProcPars(t_glsq* lsq)
	{
		if (_rec_list.empty() || _sat_list.empty())
		{
			write_log_info(_glog, 4, "ERROR : no rec or no sat in PPP_LSQ.");
			return false;
		}

		// Crd pars
		bool crd_valid = _init_rec_crd_pars(lsq);
		bool clk_valid = _init_rec_clk_pars(lsq);
		bool trop_valid =  _init_rec_trop_pars(lsq);

		bool sys_valid = _init_sys_bia_pars(lsq);
		if (_frequency >= 3)
		{
			if (_obs_mode == OBSCOMBIN::IONO_FREE)
			{
				if (!_init_rec_clk13_pars(lsq) || 
					!_init_IF_amb13_pars(lsq))return false;
			}
			else if(_obs_mode == OBSCOMBIN::RAW_ALL)
			{
				//if (!_init_ifb_bia_pars(lsq))return false; _init_sys_rec_ifb_pars		
				bool ifb_valid = false;
				if ( _ifb_model == IFB_MODEL::EST_IFB)          ifb_valid = _init_sys_ifb_pars(lsq);
				else if ( _ifb_model == IFB_MODEL::EST_REC_IFB) ifb_valid = _init_sys_rec_ifb_pars(lsq);
				else if ( _ifb_model == IFB_MODEL::EST_SAT_IFB) ifb_valid = _init_sys_sat_ifb_pars(lsq);
			
				if (!ifb_valid) return false;
			}
		}

		if (crd_valid && 
			clk_valid && 
			trop_valid && 
			sys_valid)
		{
			write_log_info(_glog, 4, "INIT : init the pars estimated in PPP_LSQ.");
			return true;
		}
		else
		{
			write_log_info(_glog, 4, "ERROR : can not init the pars estimated in PPP_LSQ.");
			return false;
		}
	}

	bool t_gppplsq::_initLsqProcData(t_gallproc * data)
	{
		gtrace("t_gppplsq::_initLsqProcData");

		// prepare obs data
		_gall_bias = dynamic_cast<t_gallbias*>((*data)[t_gdata::ALLBIAS]);
		_gall_obs = dynamic_cast<t_gallobs*>((*data)[t_gdata::ALLOBS]);
		_gall_nav = dynamic_cast<t_gallnav*>((*data)[t_gdata::GRP_EPHEM]);
		_gall_obj = dynamic_cast<t_gallobj*>((*data)[t_gdata::ALLOBJ]);
		if (_fix_mode == t_gsetamb::FIX_MODE::SEARCH || _fix_mode == t_gsetamb::FIX_MODE::ROUND)
		{
			_gupd = dynamic_cast<t_gupd*>((*data)[t_gdata::UPD]);
			//if (_frequency == 3) 
			//{ 
			//	_gupd_epoch = dynamic_cast<t_gupd_epoch*>((*data)[t_gdata::UPD_EPOCH]); 
			//}
		}

		// Whether data is valid lvhb 20201021
		if (/*!_gall_bias ||*/ !_gall_obs || !_gall_obj || (!_gall_nav && _crd_est == KIN) ||
			(!_isClient && _fix_mode != t_gsetamb::FIX_MODE::NO && _upd_mode != t_gsetamb::UPD_MODE::OSB && !_gupd))
		{
			write_log_info(_glog, 0, "ERROR ：can not init lsq process data");
			return false;
		}

		// init SPP for kinematic station
		if (_crd_est == KIN && _kin_xyz.empty())
		{
			_spplsq = new t_gspplsq(_crt_rec, _gset, _gall_obj);
			_spplsq->setDAT(_gall_obs, _gall_nav);
			//_spplsq->processBatch(_beg_time, _end_time);
		}

		// set glofrq number 
		_glofrq_num = _gall_obs->glo_freq_num();
		if (_glofrq_num.empty() && _gall_nav)
		{
			_glofrq_num = _gall_nav->glo_freq_num();
		}

		//// special for simulation PPP
		//if (_simulation)
		//{
		//	raw_CODE_bias(_beg_time, _allbias, _glofrq_num, _simu_raw_dcb);
		//}

		write_log_info(_glog, 4, "INIT : process data in lsq is prepared.");
		return true;
	}

	bool t_gppplsq::_init_sys_rec_ifb_pars(t_glsq* lsq)
	{
		if (_ifb_model != IFB_MODEL::EST_REC_IFB) return false;

		// if raw all mode, init IFB par for each satellite
		// one sat one ifb
		map<GSYS, set<t_gpar::t_type>> par_list;

		// 最多支持处理五频数据
		switch (_frequency)
		{
		case 5:
			par_list[GAL].insert(t_gpar::t_type::IFB_GAL_3);
			par_list[BDS].insert(t_gpar::t_type::IFB_BDS_3);
		case 4:
			par_list[GAL].insert(t_gpar::t_type::IFB_GAL_2);
			par_list[BDS].insert(t_gpar::t_type::IFB_BDS_2);
		case 3:
			par_list[GPS].insert(t_gpar::t_type::IFB_GPS);
			par_list[GAL].insert(t_gpar::t_type::IFB_GAL);
			par_list[BDS].insert(t_gpar::t_type::IFB_BDS);
			par_list[QZS].insert(t_gpar::t_type::IFB_QZS);
		case 2:
			break;
		default:
			return false;
		}

		// one rec one ifb
		for (const auto& rec : _rec_list)
		{
			t_gpar par_IFB;

			for (const auto& sys : _sys_list)
			{
				auto gsys = t_gsys::str2gsys(sys);
				switch (gsys)
				{
				case GPS:
				{
					for (const auto& par_type : par_list[GPS])
					{
						par_IFB = t_gpar(rec, par_type, lsq->_x_solve.parNumber() + 1, "");
						_set_isb_par_setting(lsq, par_IFB, _sysbias_model, 0.0, 3000.0);
					}
					break;
				}
				case GAL:
				{
					for (const auto& par_type : par_list[GAL])
					{
						par_IFB = t_gpar(rec, par_type, lsq->_x_solve.parNumber() + 1, "");

						_set_isb_par_setting(lsq, par_IFB, _sysbias_model, 0.0, 3000.0);
					}
					break;
				}
				case BDS:
					for (const auto& par_type : par_list[BDS])
					{
						par_IFB = t_gpar(rec, par_type, lsq->_x_solve.parNumber() + 1, "");
						_set_isb_par_setting(lsq, par_IFB, _sysbias_model, 0.0, 3000.0);
					}
					break;
				case QZS:
					for (const auto& par_type : par_list[QZS])
					{
						par_IFB = t_gpar(rec, par_type, lsq->_x_solve.parNumber() + 1, "");
						_set_isb_par_setting(lsq, par_IFB, _sysbias_model, 0.0, 3000.0);
					}
					break;
				default:
					throw logic_error("can not support sys : " + sys);
				}
			}
		}

		return true;
	}
	

	bool t_gppplsq::_processOneRecOneSat(const t_gtime & crt_epoch, const std::string & rec, const std::string & sat, t_gsatdata & crt_obs)
	{
		bool proc_valid = t_glsqproc::_processOneRecOneSat(crt_epoch, rec, sat,crt_obs);

		if (!proc_valid)
		{
			write_log_info(_glog, 0, "ERROR", "can not process the sat of " +
				sat + " and the site of " + rec + " in the time " + crt_epoch.str_ymdhms());
			return false;
		}

		if (_obs_mode == OBSCOMBIN::IONO_FREE)
		{
			//if (_simulation)
			//{
			//	/*_simu_raw_range.clear();
			//	for (auto obs : _simu_raw_dcb[sat])
			//	{
			//		double raw_obs = crd_obs.getobs(obs.first);
			//		crd_obs.resetobs(obs.first, raw_obs + obs.second);
			//		_simu_raw_range[obs.first] = raw_obs;
			//	}*/

			//	_combineMW(crd_obs);
			//	
			//	//for (auto obs : _simu_raw_range)
			//	//{
			//	//	crd_obs.resetobs(obs.first, obs.second);
			//	//}
			//}
			//else
			//{

			_combineMW(crt_obs);

			if (_frequency >= 3)
			{
				_combineEWL(crt_obs);
				_combineLW(crt_obs);
				_combineLE(crt_obs);
				_combineLWL(crt_obs);
			}
			if (_frequency >= 4)
			{
				_combineEWL24(crt_obs);
				_combineLE24(crt_obs);
				_combineLWL124(crt_obs);
			}
			if (_frequency == 5)
			{
				_combineEWL25(crt_obs);
				_combineLE25(crt_obs);
				_combineLWL125(crt_obs);
			}
		}
		else if (_obs_mode == OBSCOMBIN::RAW_ALL && !_isClient)
		{
			if (_frequency >= 3)
			{
				_combineLW(crt_obs);
				_combineLE(crt_obs);
				_combineLWL(crt_obs);
			}
			if (_frequency >= 4)
			{
				_combineLE24(crt_obs);
				_combineLWL124(crt_obs);
			}
			if (_frequency == 5)
			{
				_combineLE25(crt_obs);
				_combineLWL125(crt_obs);
			}
		}
				
			//}
		

		// Changed by wangbo: add for PPP-RTK
		_crt_ele[crt_obs.sat()] = crt_obs.ele_deg();
		
		return true;
	}

	bool t_gppplsq::_processOneRec(const t_gtime& crt_epoch, const std::string& crt_rec, std::vector<t_gsatdata>& crt_obs)
	{
		bool proc_site = t_glsqproc::_processOneRec(crt_epoch, crt_rec, crt_obs);
		if (!proc_site)
		{
			cout << crt_epoch.str_ymdhms(crt_rec + " has no equations ") << endl;
			write_log_info(_glog, 1, "NOTE", crt_rec + ": havo no equations");
			return false;
		}

		if (_frequency >= 3 && _obs_mode == IONO_FREE && _fix_mode == t_gsetamb::FIX_MODE::SEARCH)
		{
			_addFixConstraintAFIF(crt_epoch, crt_obs, "EWL");
		}
		if (_frequency >= 4 && _obs_mode == IONO_FREE && _fix_mode == t_gsetamb::FIX_MODE::SEARCH)
		{
			_addFixConstraintAFIF(crt_epoch, crt_obs, "EWL24");
		}
		if (_frequency == 5 && _obs_mode == IONO_FREE && _fix_mode == t_gsetamb::FIX_MODE::SEARCH)
		{
			_addFixConstraintAFIF(crt_epoch, crt_obs, "EWL25");
		}

		return true;
	}

	void t_gppplsq::_IFConstrainedByWL(t_glsq* lsq)
	{		
		double wl_wavelength = 0;
		t_gobsgnss* gnss = new t_gobsgnss();

		// reset MW by WL
		_MW[_crt_time].clear();
		set<string> amb_prns = _wl_params.amb_prns();
		for (auto iter = amb_prns.begin(); iter != amb_prns.end(); iter++)
		{
			t_gobsgnss gnss(*iter);
			GSYS gsys = gnss.gsys();
			wl_wavelength = gnss.wavelength_WL(_band_index[gsys][FREQ_1], _band_index[gsys][FREQ_2]);

			int index = _wl_params.getParam(_crt_rec, t_gpar::AMB_WL, *iter);
			if (index < 0)continue;
			_MW[_crt_time][*iter][1] = 10;
			_MW[_crt_time][*iter][2] = _wl_params.getAllPar()[index].value()/wl_wavelength;
			_MW[_crt_time][*iter][3] = pow(_wl_params.getAllPar()[index].apriori(), 2) * _MW[_crt_time][*iter][1];
		}

		// add WL constraint to IF
		vector<t_gpar::t_type> constraint_par;
		constraint_par.push_back(t_gpar::CRD_X);
		constraint_par.push_back(t_gpar::CRD_Y);
		constraint_par.push_back(t_gpar::CRD_Z);
		//constraint_par.push_back(t_gpar::CLK);
		double p0 = 1E4;
		for (auto iter = constraint_par.begin(); iter != constraint_par.end(); iter++)
		{
			int index_if = lsq->_x_solve.getParam(_crt_rec, *iter, "");
			int index_wl = _wl_params.getParam(_crt_rec, *iter, "");

			double dl = _wl_params[index_wl].value() - lsq->_x_solve[index_if].value();

			vector<pair<int, double>> B;
			B.push_back(make_pair(index_if + 1, 1));
			t_glsqEquationMatrix virtual_equ;
			t_gobscombtype type;
			virtual_equ.add_equ(B, p0, dl, _crt_rec, "", type, false);
			lsq->add_equation(virtual_equ);
		}
	}

	void t_gppplsq::_setOut()
	{
		// set output file
		string flt_tmp;
		flt_tmp = dynamic_cast<t_gsetout*>(_gset)->outputs("flt");
		if (!flt_tmp.empty()) {
			substitute(flt_tmp, "$(rec)", _crt_rec, false);
			if (_fltfile) {
				if (_fltfile->is_open()) _fltfile->close();
				delete _fltfile;
				_fltfile = nullptr;
			}
			_fltfile = new t_giof;
			_fltfile->tsys(t_gtime::GPS);
			_fltfile->mask(flt_tmp);
			_fltfile->append(dynamic_cast<t_gsetout*>(_gset)->append());
			_printOutHeader();
		}

		// set enu file
		string enu_tmp;
		enu_tmp = dynamic_cast<t_gsetout*>(_gset)->outputs("enu");
		if (!enu_tmp.empty()) {
			substitute(enu_tmp, "$(rec)", _crt_rec, false);
			if (_enufile) {
				if (_enufile->is_open()) _enufile->close();
				delete _enufile;
				_enufile = nullptr;
			}
			_enufile = new t_giof;
			_enufile->mask(enu_tmp);
			_enufile->append(dynamic_cast<t_gsetout*>(_gset)->append());
		}

		// rec clk file
		string clk_tmp = dynamic_cast<t_gsetout*>(_gset)->outputs("recclk");
		if (clk_tmp.empty())
		{
			clk_tmp = "rec_" + int2str(_beg_time.year()) + int2str(_beg_time.doy(), 3);
		}
		else
		{
			substitute(clk_tmp, "$(date)", int2str(_beg_time.year()) + int2str(_beg_time.doy(), 3), false);
		}
		_recclk_file = clk_tmp;

		// set trop file
		string trop_tmp = dynamic_cast<t_gsetout*>(_gset)->outputs("trop");
		if (!trop_tmp.empty()) {
			substitute(trop_tmp, "$(rec)", _crt_rec, false);
			if (_tropfile) {
				if (_tropfile->is_open()) _tropfile->close();
				delete _tropfile;
				_tropfile = nullptr;
			}
			_tropfile = new t_giof;
			_tropfile->tsys(t_gtime::GPS);
			_tropfile->mask(trop_tmp);
			_tropfile->append(dynamic_cast<t_gsetout*>(_gset)->append());
		}

		string kml_tmp = dynamic_cast<t_gsetout*>(_gset)->outputs("kml");
		if (!kml_tmp.empty()) {
			substitute(kml_tmp, "$(rec)", _crt_rec, false);
			_kml_name = kml_tmp;
			_kml = true;

			xml_node node, document, placemark, multigeometry;
			xml_node root = _doc;
			node = this->_default_node(root, _root.c_str());

			document = this->_default_node(node, "Document");

			string mode = "_LSQ";
			string desc = _crt_rec + mode;
			this->_default_node(document, "name", desc.c_str(), false);

			// for style in document
			xml_node doc_style1 = this->_default_node(document, "Style");
			string id1 = "P1";
			this->_default_attr(doc_style1, "id", id1);
			xml_node IconStyle = this->_default_node(doc_style1, "IconStyle");
			this->_default_node(IconStyle, "color", "ff008800");
			this->_default_node(IconStyle, "scale", "0.3");
			xml_node Icon = this->_default_node(IconStyle, "Icon");
			this->_default_node(Icon, "href", "http://maps.google.com/mapfiles/kml/pal2/icon18.png");

			xml_node doc_style2 = document.append_copy(doc_style1);
			string id2 = "P2";
			this->_default_attr(doc_style2, "id", id2, true);
			IconStyle = this->_default_node(doc_style2, "IconStyle");
			this->_default_node(IconStyle, "color", "ff00aaff", true);

			xml_node doc_style3 = document.append_copy(doc_style1);
			string id3 = "P3";
			this->_default_attr(doc_style3, "id", id3, true);
			IconStyle = this->_default_node(doc_style3, "IconStyle");
			this->_default_node(IconStyle, "color", "ff0000ff", true);

			// for style in placemark
			placemark = this->_default_node(document, "Placemark");
			this->_default_node(placemark, "name", "Trajection", false);
			xml_node style = this->_default_node(placemark, "Style");
			xml_node LineStyle = this->_default_node(style, "LineStyle");
			this->_default_node(LineStyle, "color", "ff00ffff", false);
			this->_default_node(LineStyle, "width", "2", false);
			// multigeometry = this->_default_node(placemark, "MultiGeometry");
			_kml_line = this->_default_node(placemark, "LineString");
			this->_default_node(_kml_line, "coordinates");
		}

	}

	int t_gppplsq::_SolveEpoch()
	{
		t_glsq lsqtmp(*_lsq);
		//_lsq->print_result_parameter(true);
		if (_isClient && _obs_mode == IONO_FREE && _frequency == 2)_IFConstrainedByWL(&lsqtmp);
		try
		{
			lsqtmp.solve_NEQ();
		}
		catch (NPDException e)
		{
			writeLogInfo(_glog, 0, "t_gppplsq::_SolveEpoch[ERROR] ", "Can't Solve NEQ! Too Few Observations!");
			return -1;
		}


	//	cout << "NEQ Matrix" << endl;
	//	cout << fixed << setw(30) << setprecision(3) << _lsq->NEQ() << endl;
		//_lsq->print_result_parameter(true);
	
		cout << _crt_time.str_ymdhms() << endl;
		for (int i = 0; i < lsqtmp.NEQ().Nrows(); i++) {
			cout << fixed << setw(20) << "Float EPO  " << setw(20) << lsqtmp._x_solve[i].str_type() + "  " <<
				setw(20) << setprecision(5) << lsqtmp._x_solve[i].value() <<
				setw(15) << setprecision(5) << lsqtmp.dx()(i + 1) <<
				setw(20) << setprecision(5) << lsqtmp._x_solve[i].value() + lsqtmp.dx()(i + 1) <<
				setw(12) << setprecision(5) << lsqtmp.stdx(i + 1) << endl;
		}

		int fix_valid = 0;   // 0 Float; 1 Fixed
		// PPP AR after 1 epoch
		if (_fix_mode != t_gsetamb::FIX_MODE::NO)
		{
			if (_realtime == true) {
				_ewl_Upd_time = _crt_time;
				_ewl24_Upd_time = _crt_time;
				_ewl25_Upd_time = _crt_time;
				_wl_Upd_time = _crt_time;
			}
			if (_isFirstFix)
			{
				_ambfix->setUPD(_gupd);
				//_ambfix->setUPD_epoch(_gupd_epoch);
				_ambfix->setLOG(_glog);
				_ambfix->setWaveLength(_glofrq_num);
				_ambfix->setObsType(_obs_mode);
				if (_slip12) _ambfix->setActiveAmb(_slip12->get_active_amb(_crt_rec));
				_isFirstFix = false;
			}

			// getting the reference satellite
			if (dynamic_cast<t_gsetamb*>(_gset)->isSetRefSat())
			{
				bool ref_valid = _getSatRef();
				if (ref_valid)_ambfix->setSatRef(_sat_ref);
				else return 0;
			}

			//t_gamb_cmn amb_cmn(_crt_time, &lsqtmp);
			if (_obs_mode == OBSCOMBIN::IONO_FREE)
			{
				_ambfix->setMW(_MW[_crt_time]);
				fix_valid = _ambfix->processBatch(_crt_time, &lsqtmp, _slip12.get(), "NL");
			}
			else if (_obs_mode == OBSCOMBIN::RAW_ALL || _obs_mode == OBSCOMBIN::RAW_MIX)
			{
				_ambfix->setELE(_crt_ele);
				if (!(_sys_list.find("BDS") != _sys_list.end() && _sys_list.size() == 1))
				{
					if (_frequency >= 3) { _ambfix->processBatch(_crt_time, &lsqtmp, _slip13.get(), "EWL"); }
					if (_frequency >= 4) { _ambfix->processBatch(_crt_time, &lsqtmp, _slip14.get(), "EWL24"); }
					if (_frequency == 5) { _ambfix->processBatch(_crt_time, &lsqtmp, _slip15.get(), "EWL25"); }
				}
				_ambfix->processBatch(_crt_time, &lsqtmp, _slip12.get(), "WL");
				if (_frequency >= 3 && !_isClient) { _addFixConstraintAFIF_RAW(_crt_time, &lsqtmp, "EWL"); }
				if (_frequency >= 4 && !_isClient) { _addFixConstraintAFIF_RAW(_crt_time, &lsqtmp, "EWL24"); }
				if (_frequency == 5 && !_isClient) { _addFixConstraintAFIF_RAW(_crt_time, &lsqtmp, "EWL25"); }
			    fix_valid = _ambfix->processBatch(_crt_time, &lsqtmp, _slip12.get(), "NL");
				
			}
			else if (_isClient && _obs_mode == OBSCOMBIN::WL_COMBIN)
			{
				fix_valid = _ambfix->processBatch(_crt_time, &lsqtmp, _slip12.get(), "NL");
			}

			if (fix_valid == 1)
			{
				//cout << _crt_time.str_ymdhms() << endl;
				cout << " After Ambiguity Fix: " << _crt_time.str_ymdhms() << endl;
				for (int i = 0; i < lsqtmp.NEQ().Nrows(); i++) {
					cout << fixed << setw(20) << "Fix EPO  " << setw(20) << lsqtmp._x_solve[i].str_type() + "  " <<
						setw(20) << setprecision(5) << lsqtmp._x_solve[i].value() <<
						setw(15) << setprecision(5) << lsqtmp.dx()(i + 1) <<
						setw(20) << setprecision(5) << lsqtmp._x_solve[i].value() + lsqtmp.dx()(i + 1) <<
						setw(12) << setprecision(5) << lsqtmp.stdx(i + 1) << endl;
					_lsq->_x_solve[i].amb_ini = false;
				}
			}
			else
			{
				fix_valid = 0;
				//if ambiguity can't be fixed in PPP-RTK, the ambiguity in the next epoch will be initialize
				if (_isClient) 
				{
					for (int i = 0; i < lsqtmp.NEQ().Nrows(); i++) _lsq->_x_solve[i].amb_ini = true;
				}
			}
			
			if (_obs_mode == WL_COMBIN)
			{
				return fix_valid;
			}
		}

		//LX changed: if fix_valid is false, print the float solutions
		if (fix_valid==1)
		{
			_printOut(lsqtmp, "Fixed");
		}
		else
		{
			_printOut(lsqtmp, "Float");

		}
		
		
		return fix_valid;
	}

	void t_gppplsq::_printOutHeader()
	{
		ostringstream header_os;

		// output formaat
		header_os << "%EPO " << fixed << setprecision(4)
			<< setw(24) << "DATE[Y-M-D h:m:s]"
			<< setw(17) << "Sec of Week[s]"
			<< setw(16) << "Sec of Day[s]"
			// write data
			<< setw(15) << "CRD_CRDX[m]"
			<< setw(15) << "CRD_CRDY[m]"
			<< setw(15) << "CRD_CRDZ[m]"
			<< setw(15) << "SIGMA_X[m]"
			<< setw(15) << "SIGMA_Y[m]"
			<< setw(15) << "SIGMA_Z[m]"
			<< setw(15) << "REC_CLK[m]"
			<< setw(15) << "SIGMA_CLK[m]"
			<< setw(15) << "NSAT"
			<< setw(15) << "SIGMA0[m]"
			<< setw(15) << "GDOP"
			<< setw(15) << "PDOP"
			<< setw(15) << "AmbStatus"
			<< endl;

		if (_fltfile) {
			_fltfile->write(header_os.str().c_str(), header_os.str().size());
			_fltfile->flush();
		}

	}

	void t_gppplsq::_printOut(const t_glsq& lsqestimator, string mode)
	{
		// print ztd/grd
		_extractTropDelay(lsqestimator);

		ostringstream os;

		string _site = *_rec_list.begin();
		t_gallpar par = lsqestimator._x_solve;
		ColumnVector dx = lsqestimator.dx();
		ColumnVector stdx = lsqestimator.stdx();
		string str_dsec = dbl2str(_crt_time.dsec());
		int nsat = (par.amb_prns()).size();

		// get XYZ and XYZ_sig & output enu coordinates
		t_gtriple xyz(0.0, 0.0, 0.0);
		t_gtriple xyz_sigma(0.0, 0.0, 0.0);
		if (_crd_est != CONSTRPAR::FIX)
		{
			int idx_x, idx_y, idx_z;
			idx_x = par.getParam(_site, t_gpar::CRD_X, "");
			idx_y = par.getParam(_site, t_gpar::CRD_Y, "");
			idx_z = par.getParam(_site, t_gpar::CRD_Z, "");
			xyz[0] = par[idx_x].value() + dx(idx_x + 1);
			xyz[1] = par[idx_y].value() + dx(idx_y + 1);
			xyz[2] = par[idx_z].value() + dx(idx_z + 1);
			xyz_sigma[0] = stdx(idx_x + 1);
			xyz_sigma[1] = stdx(idx_y + 1);
			xyz_sigma[2] = stdx(idx_z + 1);

			// write enu
			if (_enufile)
			{
				t_gtriple dxyz, ell, neu;

				ostringstream enu_info;
				double crt_t = _crt_time.sod() + _crt_time.dsec();
				

				if (_crd_est == CONSTRPAR::KIN)
				{
					// just for now
					enu_info <<  _crt_rec << " " << fixed << setw(8) << _crt_time.mjd() << 
						fixed << setw(10) << setprecision(1) << crt_t; 
					//enu_info << fixed << setw(10) << setprecision(1) << crt_t;
					enu_info << setw(16) << setprecision(4) << xyz[0]
						<< setw(16) << setprecision(4) << xyz[1]
						<< setw(16) << setprecision(4) << xyz[2] << "  " << nsat << endl;
				}
				else
				{
					enu_info << fixed << setw(10) << setprecision(1) << crt_t;
					xyz2ell(xyz_standard, ell, false);
					dxyz = xyz - xyz_standard;
					xyz2neu(ell, dxyz, neu);

					// RMS
					neu_sum_nepo.first++;   // nepoch
					neu_sum_nepo.second += t_gtriple(pow(neu.crd(0), 2), pow(neu.crd(1), 2), pow(neu.crd(2), 2)); // sum

					enu_info << fixed << setw(16) << setprecision(4) << neu[1]
						<< setw(16) << setprecision(4) << neu[0]
						<< setw(16) << setprecision(4) << neu[2] << endl;
				}
				_enufile->write(enu_info.str().c_str(), enu_info.str().size());
				_enufile->flush();
			}

			// write kml
			if (_kml) {
				t_gtriple ell;
				xyz2ell(xyz, ell, false);
				ostringstream out;
				out << fixed << setprecision(11) << " " << setw(0) << ell[1] * R2D << ',' << ell[0] * R2D;
				string val = out.str();

				xml_node root = _doc;
				_kml_line = this->_default_node(root, _root.c_str());
				xml_node document = _kml_line.child("Document");
				xml_node last_child = document.last_child();
				xml_node placemark = document.insert_child_after("Placemark", last_child);
				if (mode == "Fixed")this->_default_node(placemark, "styleUrl", "#P1"); // value: #P1 (according to precision)
				else this->_default_node(placemark, "styleUrl", "#P3");
				//this->_default_node(placemark, "name", int2str(_crt_time.sow()).c_str());
				xml_node point = this->_default_node(placemark, "Point");
				this->_default_node(point, "coordinates", val.c_str()); // for point

				xml_node Placemark = document.child("Placemark");
				xml_node LineString = Placemark.child("LineString");
				this->_default_node(LineString, "coordinates", val.c_str(), false);  // for line
			}
		}
		
		// get CLK param
		double clk, clk_sigma;
		clk = 0.0; clk_sigma = 0.0;
		int idx_clk = par.getParam(_site, t_gpar::CLK, "");
		if (idx_clk >= 0)
		{
			clk = par[idx_clk].value() + dx(idx_clk + 1);
			clk_sigma = stdx(idx_clk + 1);
		}

		// get unit sigma0
		double sigma0 = 0.0;
		sigma0 = lsqestimator.sigma0();

		// get GDOP PDOP
		double GDOP, PDOP;
		double qxx, qyy, qzz, qtt;
		qxx = pow(xyz_sigma[0] / sigma0, 2);
		qyy = pow(xyz_sigma[1] / sigma0, 2);
		qzz = pow(xyz_sigma[2] / sigma0, 2);
		qtt = pow(clk_sigma / sigma0, 2);
		GDOP = sqrt(qxx + qyy + qzz + qtt);
		PDOP = sqrt(qxx + qyy + qzz);

		// output formaat
		os << "=EPO " << fixed << setprecision(4)
			<< _crt_time.str_ymdhms("", false, true) << str_dsec.substr(2) 
			<< setw(17) << _crt_time.sow() + _crt_time.dsec()
			<< setw(16) << _crt_time.sod() + _crt_time.dsec() << setprecision(4);

		// write data
		os << setw(15) << xyz[0]           //[m]
			<< setw(15) << xyz[1]           //[m]
			<< setw(15) << xyz[2]           //[m]
			<< setw(15) << xyz_sigma[0]     //[m]
			<< setw(15) << xyz_sigma[1]     //[m]
			<< setw(15) << xyz_sigma[2]     //[m]
			<< setw(15) << clk              //[m]
			<< setw(15) << clk_sigma        //[m]
			<< setw(15) << nsat
			<< setw(15) << sigma0           //[m]
			<< setw(15) << GDOP
			<< setw(15) << PDOP
			<< setw(15) << mode
			<< endl;

		// Print flt results
		if (_fltfile) 
		{
			_fltfile->write(os.str().c_str(), os.str().size());
			_fltfile->flush();
		}
		else 
		{
			_glog->comment(1, "t_gspplsq::_printout", "have no output file!");
		}
	}

	void t_gppplsq::_extractTropDelay(const t_glsq& lsqestimator)
	{
		if (!_tropfile) return;

		ostringstream os;
		map<t_gtime, map<t_gpar::t_type, double> > TropDelay;
		map<t_gtime, int> nsat;
		if (_lsq->mode() == LSQMODE::LSQ)
		{
			t_gtime it_time;
			// get recover trop delay
			vector<t_grecover_par> par_Trop, par_tmp;
			par_Trop = _gall_rcv->get_recover_par(t_gpar::TRP);
			par_tmp = _gall_rcv->get_recover_par(t_gpar::GRD_N);
			par_Trop.insert(par_Trop.end(), par_tmp.begin(), par_tmp.end());
			par_tmp = _gall_rcv->get_recover_par(t_gpar::GRD_E);
			par_Trop.insert(par_Trop.end(), par_tmp.begin(), par_tmp.end());
			// store trop according time
			for (auto iter = par_Trop.begin(); iter != par_Trop.end(); iter++)
			{
				if (iter->par.site != _crt_rec) continue;
				double value = iter->par.value() + iter->correct_value;
				if (iter->par.parType == t_gpar::TRP) value += iter->par.zhd;
				if (double_eq(value, 0.0)) continue;
				it_time = iter->par.beg;
				for (; it_time <= iter->par.end; it_time = it_time + _obs_intv)
				{
					if (it_time > _end_time) continue;
					TropDelay[it_time][iter->par.parType] = value;
				}
			}
			// get number of sats
			nsat = _gall_rcv->get_sat_number_per_epo();
		}
		// current epo tropospheric delays
		t_gallpar par = lsqestimator._x_solve;
		ColumnVector dx = lsqestimator.dx();
		_crt_time = _crt_time < _end_time ? _crt_time : _end_time;
		TropDelay[_crt_time][t_gpar::TRP] = 0.0;
		TropDelay[_crt_time][t_gpar::GRD_N] = 0.0;
		TropDelay[_crt_time][t_gpar::GRD_E] = 0.0;
		nsat[_crt_time] = par.amb_prns().size();
		int idx_ztd, idx_zgn, idx_zge;
		idx_ztd = par.getParam(_crt_rec, t_gpar::TRP, "");
		idx_zgn = par.getParam(_crt_rec, t_gpar::GRD_N, "");
		idx_zge = par.getParam(_crt_rec, t_gpar::GRD_E, "");
		if (idx_ztd >= 0) TropDelay[_crt_time][t_gpar::TRP] = par[idx_ztd].value() + dx(idx_ztd + 1) + par[idx_ztd].zhd;
		if (idx_zgn >= 0) TropDelay[_crt_time][t_gpar::GRD_N] = par[idx_zgn].value() + dx(idx_zgn + 1);
		if (idx_zge >= 0) TropDelay[_crt_time][t_gpar::GRD_E] = par[idx_zge].value() + dx(idx_zge + 1);
		
		// output  unit mm
		for (auto it : TropDelay)
		{
			os.str("");
			os << fixed << it.first.str_yyyy() << " " << setw(5) << it.first.doy()
				<< setw(8) << it.first.sod() << setprecision(3)
				<< setw(12) << it.second[t_gpar::TRP] * 1000.0
				<< setw(12) << it.second[t_gpar::GRD_N] * 1000.0
				<< setw(12) << it.second[t_gpar::GRD_E] * 1000.0
				<< setw(12) << nsat[it.first] << endl;
			_tropfile->write(os.str().c_str(), os.str().size());
			_tropfile->flush();
		}
		return;
	}


	void t_gppplsq::_combineMW(t_gsatdata& satdata)
	{
		string sat = satdata.sat();

		if (_crt_time == _beg_time)
		{
			_MW[_beg_time - _obs_intv][sat][1] = 0.0;
			_MW[_beg_time - _obs_intv][sat][2] = 0.0;
			_MW[_beg_time - _obs_intv][sat][3] = 0.0;
			_MW[_beg_time - _obs_intv][sat][4] = 0.0;
			_MW[_beg_time - _obs_intv][sat][5] = 0.0; // ele
		}

		GOBSBAND b1 = _band_index[satdata.gsys()][FREQ_1];
		GOBSBAND b2 = _band_index[satdata.gsys()][FREQ_2];
		t_gobs gobsP1 = t_gobs(satdata.select_range(b1));
		t_gobs gobsP2 = t_gobs(satdata.select_range(b2));
		t_gobs gobsL1 = t_gobs(satdata.select_phase(b1));
		t_gobs gobsL2 = t_gobs(satdata.select_phase(b2));
		double mw_obs = satdata.MW_cycle(gobsL1, gobsL2, gobsP1, gobsP2);

		if (!double_eq(mw_obs, 0.0))
		{
			// calculate smooth MW
			t_gtime pre_time = _crt_time - _obs_intv;
			if (_MW[pre_time].find(sat) != _MW[pre_time].end() && !double_eq(_MW[pre_time][sat][4], 0.0))
			{
				_MW[_crt_time][sat][1] = _MW[pre_time][sat][1] + 1;
				_MW[_crt_time][sat][4] = mw_obs;
				if (!double_eq(_MW[_crt_time][sat][1], 1.0) && _MW[_crt_time][sat][1] > 1.0)
				{
					_MW[_crt_time][sat][3] = _MW[pre_time][sat][3] * (_MW[_crt_time][sat][1] - 2) /
						(_MW[_crt_time][sat][1] - 1) + pow(mw_obs - _MW[pre_time][sat][2], 2)
						/ _MW[_crt_time][sat][1];
				}
				else
				{
					_MW[_crt_time][sat][3] = _MW[pre_time][sat][3];
				}
				_MW[_crt_time][sat][2] = _MW[pre_time][sat][2] + (mw_obs - _MW[pre_time][sat][2])
					/ _MW[_crt_time][sat][1];
				_MW[_crt_time][sat][5] = satdata.ele_deg();

				// for output ambupd 
				_AmbValidTime[sat].rbegin()->second = _crt_time;
			}
			// exclude mw begining ,it is different from last all
			else
			{
				_MW[_crt_time][sat][1] = 1;
				_MW[_crt_time][sat][4] = mw_obs;
				_MW[_crt_time][sat][3] = 0.0;
				_MW[_crt_time][sat][2] = mw_obs / _MW[_crt_time][sat][1];
				_MW[_crt_time][sat][5] = satdata.ele_deg();

				// for output ambupd
				if(_AmbValidTime.find(sat) == _AmbValidTime.end() 
					|| _AmbValidTime[sat].rbegin()->second != _crt_time)
				{
					_AmbValidTime[sat].push_back(make_pair(_crt_time, _crt_time));
				}
			}
			if (_frequency >= 3 && _fix_mode!= t_gsetamb::FIX_MODE::NO && _upd_mode== t_gsetamb::UPD_MODE::UPD)
			{
				one_epoch_upd epoch_upd = _gupd->get_epo_upd(UPDTYPE::WL, WL_IDENTIFY);
				if (epoch_upd.find(sat) != epoch_upd.end() && epoch_upd[sat]->sigma <= 0.2 &&
					epoch_upd[sat]->npoint > 2)
				{
					_rMW[_crt_time][sat][1] = _MW[_crt_time][sat][2] - epoch_upd[sat]->value;
				}
			}
			else if (_frequency >= 3 && _fix_mode != t_gsetamb::FIX_MODE::NO && _upd_mode != t_gsetamb::UPD_MODE::UPD)
			{
				_rMW[_crt_time][sat][1] = _MW[_crt_time][sat][2];
			}
		}
		// exclude zero wl
		else
		{
			_MW[_crt_time][sat][1] = 0.0;
			_MW[_crt_time][sat][2] = 0.0;
			_MW[_crt_time][sat][3] = 0.0;
			_MW[_crt_time][sat][4] = 0.0;
			_MW[_crt_time][sat][5] = 0.0;
		}
	//	satdata.apply_bias_tmp2(*_gallbias);
	}

	void t_gppplsq::_combineEWL(t_gsatdata& satdata)
	{
		string sat = satdata.sat();
		if (_crt_time == _beg_time)
		{
			_EWL[_beg_time - _obs_intv][sat][1] = 0.0;
			_EWL[_beg_time - _obs_intv][sat][2] = 0.0;
			_EWL[_beg_time - _obs_intv][sat][3] = 0.0;
			_EWL[_beg_time - _obs_intv][sat][4] = 0.0;
			_EWL[_beg_time - _obs_intv][sat][5] = 0.0; // ele
		}

		//GOBSBAND b1 = t_gsys::band_priority(satdata.gsys(), FREQ_1);
		//GOBSBAND b2 = t_gsys::band_priority(satdata.gsys(), FREQ_2);
		//GOBSBAND b3 = t_gsys::band_priority(satdata.gsys(), FREQ_3);
		GOBSBAND b1 = _band_index[satdata.gsys()][FREQ_1];
		GOBSBAND b2 = _band_index[satdata.gsys()][FREQ_2];
		GOBSBAND b3 = _band_index[satdata.gsys()][FREQ_3];
		t_gobs gobsP1 = t_gobs(satdata.select_range(b1));
		t_gobs gobsP2 = t_gobs(satdata.select_range(b2));
		t_gobs gobsL1 = t_gobs(satdata.select_phase(b1));
		t_gobs gobsL2 = t_gobs(satdata.select_phase(b2));
		t_gobs gobsL3 = t_gobs(satdata.select_phase(b3));
		double ewl_obs = satdata.EWL_cycle(gobsL1, gobsL2, gobsL3, gobsP1, gobsP2);
		////dcb correct
		//OBSCOMBIN type = OBSCOMBIN::EWL_COMBIN;

		if (!double_eq(ewl_obs, 0.0))
		{
			//apply_DCB(_allbias, satdata, ewl_obs, type, &gobsP1, &gobsP2);
			//correct_DCB("CAS", _allbias, satdata, ewl_obs, type, &gobsP1, &gobsP2);
			// calculate smooth MW
			if (_EWL[_crt_time - _obs_intv].find(sat) != _EWL[_crt_time - _obs_intv].end())
			{
				_EWL[_crt_time][sat][1] = _EWL[_crt_time - _obs_intv][sat][1] + 1;
				_EWL[_crt_time][sat][4] = ewl_obs;
				if (!double_eq(_EWL[_crt_time][sat][1], 1.0) && _EWL[_crt_time][sat][1] > 1.0)
				{
					_EWL[_crt_time][sat][3] = _EWL[_crt_time - _obs_intv][sat][3] * (_EWL[_crt_time][sat][1] - 2) /
						(_EWL[_crt_time][sat][1] - 1) + pow(ewl_obs - _EWL[_crt_time - _obs_intv][sat][2], 2)
						/ _EWL[_crt_time][sat][1];
				}
				else
				{
					_EWL[_crt_time][sat][3] = _EWL[_crt_time - _obs_intv][sat][3];
				}
				_EWL[_crt_time][sat][2] = _EWL[_crt_time - _obs_intv][sat][2] + (ewl_obs - _EWL[_crt_time - _obs_intv][sat][2])
					/ _EWL[_crt_time][sat][1];
				_EWL[_crt_time][sat][5] = satdata.ele_deg();
			}
			// exclude ewl begining ,it is different from last all
			else
			{
				_EWL[_crt_time][sat][1] = 1;
				_EWL[_crt_time][sat][4] = ewl_obs;
				_EWL[_crt_time][sat][3] = 0.0;
				_EWL[_crt_time][sat][2] = ewl_obs / _EWL[_crt_time][sat][1];
				_EWL[_crt_time][sat][5] = satdata.ele_deg();
			}
			if (_fix_mode != t_gsetamb::FIX_MODE::NO && _upd_mode == t_gsetamb::UPD_MODE::UPD)
			{
				one_epoch_upd epoch_upd = _gupd->get_epo_upd(UPDTYPE::EWL_epoch, _crt_time);
				if (epoch_upd.find(sat) != epoch_upd.end())
				{
					if (epoch_upd[sat]->sigma <= 0.2 && epoch_upd[sat]->npoint >= 2)
					{
						_rEWL[_crt_time][sat][1] = _EWL[_crt_time][sat][2] - epoch_upd[sat]->value;
					}
				}
				else
				{
					epoch_upd = _gupd->get_epo_upd(UPDTYPE::EWL, EWL_IDENTIFY);
					if (epoch_upd.find(sat) != epoch_upd.end() && epoch_upd[sat]->sigma <= 0.2 && epoch_upd[sat]->npoint >= 2)
					{
						_rEWL[_crt_time][sat][1] = _EWL[_crt_time][sat][2] - epoch_upd[sat]->value;
					}

				}
			}
			else if (_fix_mode != t_gsetamb::FIX_MODE::NO && _upd_mode != t_gsetamb::UPD_MODE::UPD)
			{
				_rEWL[_crt_time][sat][1] = _EWL[_crt_time][sat][2];
			}
		}
		// exclude zero wl
		else
		{
			_EWL[_crt_time][sat][1] = 0.0;
			_EWL[_crt_time][sat][2] = 0.0;
			_EWL[_crt_time][sat][3] = 0.0;
			_EWL[_crt_time][sat][4] = 0.0;
			_EWL[_crt_time][sat][5] = 0.0;
		}
	}

	void t_gppplsq::_combineEWL24(t_gsatdata& satdata)
	{
		string sat = satdata.sat();
		if (_crt_time == _beg_time)
		{
			_EWL24[_beg_time - _obs_intv][sat][1] = 0.0;
			_EWL24[_beg_time - _obs_intv][sat][2] = 0.0;
			_EWL24[_beg_time - _obs_intv][sat][3] = 0.0;
			_EWL24[_beg_time - _obs_intv][sat][4] = 0.0;
			_EWL24[_beg_time - _obs_intv][sat][5] = 0.0; // ele
		}

		//GOBSBAND b1 = t_gsys::band_priority(satdata.gsys(), FREQ_1);
		//GOBSBAND b2 = t_gsys::band_priority(satdata.gsys(), FREQ_2);
		//GOBSBAND b3 = t_gsys::band_priority(satdata.gsys(), FREQ_3);
		GOBSBAND b1 = _band_index[satdata.gsys()][FREQ_1];
		GOBSBAND b2 = _band_index[satdata.gsys()][FREQ_2];
		GOBSBAND b3 = _band_index[satdata.gsys()][FREQ_4];
		t_gobs gobsP1 = t_gobs(satdata.select_range(b1));
		t_gobs gobsP2 = t_gobs(satdata.select_range(b2));
		t_gobs gobsL1 = t_gobs(satdata.select_phase(b1));
		t_gobs gobsL2 = t_gobs(satdata.select_phase(b2));
		t_gobs gobsL3 = t_gobs(satdata.select_phase(b3));
		double ewl24_obs = satdata.EWL_cycle(gobsL1, gobsL2, gobsL3, gobsP1, gobsP2);
		////dcb correct
		//OBSCOMBIN type = OBSCOMBIN::EWL_COMBIN;

		if (!double_eq(ewl24_obs, 0.0))
		{
			//apply_DCB(_allbias, satdata, ewl_obs, type, &gobsP1, &gobsP2);
			//correct_DCB("CAS", _allbias, satdata, ewl_obs, type, &gobsP1, &gobsP2);
			// calculate smooth MW
			if (_EWL24[_crt_time - _obs_intv].find(sat) != _EWL24[_crt_time - _obs_intv].end())
			{
				_EWL24[_crt_time][sat][1] = _EWL24[_crt_time - _obs_intv][sat][1] + 1;
				_EWL24[_crt_time][sat][4] = ewl24_obs;
				if (!double_eq(_EWL24[_crt_time][sat][1], 1.0) && _EWL24[_crt_time][sat][1] > 1.0)
				{
					_EWL24[_crt_time][sat][3] = _EWL24[_crt_time - _obs_intv][sat][3] * (_EWL24[_crt_time][sat][1] - 2) /
						(_EWL24[_crt_time][sat][1] - 1) + pow(ewl24_obs - _EWL24[_crt_time - _obs_intv][sat][2], 2)
						/ _EWL24[_crt_time][sat][1];
				}
				else
				{
					_EWL24[_crt_time][sat][3] = _EWL24[_crt_time - _obs_intv][sat][3];
				}
				_EWL24[_crt_time][sat][2] = _EWL24[_crt_time - _obs_intv][sat][2] + (ewl24_obs - _EWL24[_crt_time - _obs_intv][sat][2])
					/ _EWL24[_crt_time][sat][1];
				_EWL24[_crt_time][sat][5] = satdata.ele_deg();
			}
			// exclude ewl begining ,it is different from last all
			else
			{
				_EWL24[_crt_time][sat][1] = 1;
				_EWL24[_crt_time][sat][4] = ewl24_obs;
				_EWL24[_crt_time][sat][3] = 0.0;
				_EWL24[_crt_time][sat][2] = ewl24_obs / _EWL24[_crt_time][sat][1];
				_EWL24[_crt_time][sat][5] = satdata.ele_deg();
			}
			if (_fix_mode != t_gsetamb::FIX_MODE::NO && _upd_mode == t_gsetamb::UPD_MODE::UPD)
			{
				one_epoch_upd epoch_upd = _gupd->get_epo_upd(UPDTYPE::EWL24, EWL24_IDENTIFY);
				if (epoch_upd.find(sat) != epoch_upd.end() && epoch_upd[sat]->sigma <= 0.2 && epoch_upd[sat]->npoint >= 2)
				{
					_rEWL24[_crt_time][sat][1] = _EWL24[_crt_time][sat][2] - epoch_upd[sat]->value;
				}
			}
			else if (_fix_mode != t_gsetamb::FIX_MODE::NO && _upd_mode != t_gsetamb::UPD_MODE::UPD)
			{
				_rEWL24[_crt_time][sat][1] = _EWL24[_crt_time][sat][2];
			}
		}
		// exclude zero wl
		else
		{
			_EWL24[_crt_time][sat][1] = 0.0;
			_EWL24[_crt_time][sat][2] = 0.0;
			_EWL24[_crt_time][sat][3] = 0.0;
			_EWL24[_crt_time][sat][4] = 0.0;
			_EWL24[_crt_time][sat][5] = 0.0;
		}
	}

	void t_gppplsq::_combineEWL25(t_gsatdata& satdata)
	{
		string sat = satdata.sat();
		if (_crt_time == _beg_time)
		{
			_EWL25[_beg_time - _obs_intv][sat][1] = 0.0;
			_EWL25[_beg_time - _obs_intv][sat][2] = 0.0;
			_EWL25[_beg_time - _obs_intv][sat][3] = 0.0;
			_EWL25[_beg_time - _obs_intv][sat][4] = 0.0;
			_EWL25[_beg_time - _obs_intv][sat][5] = 0.0; // ele
		}

		//GOBSBAND b1 = t_gsys::band_priority(satdata.gsys(), FREQ_1);
		//GOBSBAND b2 = t_gsys::band_priority(satdata.gsys(), FREQ_2);
		//GOBSBAND b3 = t_gsys::band_priority(satdata.gsys(), FREQ_3);
		GOBSBAND b1 = _band_index[satdata.gsys()][FREQ_1];
		GOBSBAND b2 = _band_index[satdata.gsys()][FREQ_2];
		GOBSBAND b3 = _band_index[satdata.gsys()][FREQ_5];
		t_gobs gobsP1 = t_gobs(satdata.select_range(b1));
		t_gobs gobsP2 = t_gobs(satdata.select_range(b2));
		t_gobs gobsL1 = t_gobs(satdata.select_phase(b1));
		t_gobs gobsL2 = t_gobs(satdata.select_phase(b2));
		t_gobs gobsL3 = t_gobs(satdata.select_phase(b3));
		double ewl25_obs = satdata.EWL_cycle(gobsL1, gobsL2, gobsL3, gobsP1, gobsP2);
		////dcb correct
		//OBSCOMBIN type = OBSCOMBIN::EWL_COMBIN;

		if (!double_eq(ewl25_obs, 0.0))
		{
			//apply_DCB(_allbias, satdata, ewl_obs, type, &gobsP1, &gobsP2);
			//correct_DCB("CAS", _allbias, satdata, ewl_obs, type, &gobsP1, &gobsP2);
			// calculate smooth MW
			if (_EWL25[_crt_time - _obs_intv].find(sat) != _EWL25[_crt_time - _obs_intv].end())
			{
				_EWL25[_crt_time][sat][1] = _EWL25[_crt_time - _obs_intv][sat][1] + 1;
				_EWL25[_crt_time][sat][4] = ewl25_obs;
				if (!double_eq(_EWL25[_crt_time][sat][1], 1.0) && _EWL25[_crt_time][sat][1] > 1.0)
				{
					_EWL25[_crt_time][sat][3] = _EWL25[_crt_time - _obs_intv][sat][3] * (_EWL25[_crt_time][sat][1] - 2) /
						(_EWL25[_crt_time][sat][1] - 1) + pow(ewl25_obs - _EWL25[_crt_time - _obs_intv][sat][2], 2)
						/ _EWL25[_crt_time][sat][1];
				}
				else
				{
					_EWL25[_crt_time][sat][3] = _EWL25[_crt_time - _obs_intv][sat][3];
				}
				_EWL25[_crt_time][sat][2] = _EWL25[_crt_time - _obs_intv][sat][2] + (ewl25_obs - _EWL25[_crt_time - _obs_intv][sat][2])
					/ _EWL25[_crt_time][sat][1];
				_EWL25[_crt_time][sat][5] = satdata.ele_deg();
			}
			// exclude ewl begining ,it is different from last all
			else
			{
				_EWL25[_crt_time][sat][1] = 1;
				_EWL25[_crt_time][sat][4] = ewl25_obs;
				_EWL25[_crt_time][sat][3] = 0.0;
				_EWL25[_crt_time][sat][2] = ewl25_obs / _EWL25[_crt_time][sat][1];
				_EWL25[_crt_time][sat][5] = satdata.ele_deg();
			}
			if (_fix_mode != t_gsetamb::FIX_MODE::NO && _upd_mode == t_gsetamb::UPD_MODE::UPD)
			{
				one_epoch_upd epoch_upd = _gupd->get_epo_upd(UPDTYPE::EWL25, EWL25_IDENTIFY);
				if (epoch_upd.find(sat) != epoch_upd.end() && epoch_upd[sat]->sigma <= 0.2 && epoch_upd[sat]->npoint >= 2)
				{
					_rEWL25[_crt_time][sat][1] = _EWL25[_crt_time][sat][2] - epoch_upd[sat]->value;
				}
			}
			else if (_fix_mode != t_gsetamb::FIX_MODE::NO && _upd_mode != t_gsetamb::UPD_MODE::UPD)
			{
				_rEWL25[_crt_time][sat][1] = _EWL25[_crt_time][sat][2];
			}
		}
		// exclude zero wl
		else
		{
			_EWL25[_crt_time][sat][1] = 0.0;
			_EWL25[_crt_time][sat][2] = 0.0;
			_EWL25[_crt_time][sat][3] = 0.0;
			_EWL25[_crt_time][sat][4] = 0.0;
			_EWL25[_crt_time][sat][5] = 0.0;
		}
	}

	void t_gppplsq::_combineLW(t_gsatdata& satdata)
	{
		string sat = satdata.sat();
		if (_crt_time == _beg_time)
		{
			_LW[_beg_time - _obs_intv][sat][1] = 0.0;
		}

		//GOBSBAND b1 = t_gsys::band_priority(satdata.gsys(), FREQ_1);
		//GOBSBAND b2 = t_gsys::band_priority(satdata.gsys(), FREQ_2);
		GOBSBAND b1 = _band_index[satdata.gsys()][FREQ_1];
		GOBSBAND b2 = _band_index[satdata.gsys()][FREQ_2];
		t_gobs gobsL1 = t_gobs(satdata.select_phase(b1));
		t_gobs gobsL2 = t_gobs(satdata.select_phase(b2));

		double LW_obs = satdata.LW_meter(gobsL1, gobsL2);

		if (!double_eq(LW_obs, 0.0))
		{
			_LW[_crt_time][sat][1] = LW_obs;
		}
		else
		{
			_LW[_crt_time][sat][1] = 0.0;
		}
	}

	void t_gppplsq::_combineLE(t_gsatdata& satdata)
	{
		string sat = satdata.sat();
		if (_crt_time == _beg_time)
		{
			_LE[_beg_time - _obs_intv][sat][1] = 0.0;
		}

		//GOBSBAND b2 = t_gsys::band_priority(satdata.gsys(), FREQ_2);
		//GOBSBAND b3 = t_gsys::band_priority(satdata.gsys(), FREQ_3);
		GOBSBAND b2 = _band_index[satdata.gsys()][FREQ_2];
		GOBSBAND b3 = _band_index[satdata.gsys()][FREQ_3];
		t_gobs gobsL2 = t_gobs(satdata.select_phase(b2));
		t_gobs gobsL3 = t_gobs(satdata.select_phase(b3));

		double LE_obs = satdata.LE_meter(gobsL2, gobsL3);

		if (!double_eq(LE_obs, 0.0))
		{
			_LE[_crt_time][sat][1] = LE_obs;
		}
		else
		{
			_LE[_crt_time][sat][1] = 0.0;
		}
	}

	void t_gppplsq::_combineLE24(t_gsatdata& satdata)
	{
		string sat = satdata.sat();
		if (_crt_time == _beg_time)
		{
			_LE24[_beg_time - _obs_intv][sat][1] = 0.0;
		}

		//GOBSBAND b2 = t_gsys::band_priority(satdata.gsys(), FREQ_2);
		//GOBSBAND b3 = t_gsys::band_priority(satdata.gsys(), FREQ_3);
		GOBSBAND b2 = _band_index[satdata.gsys()][FREQ_2];
		GOBSBAND b4 = _band_index[satdata.gsys()][FREQ_4];
		t_gobs gobsL2 = t_gobs(satdata.select_phase(b2));
		t_gobs gobsL4 = t_gobs(satdata.select_phase(b4));

		double LE24_obs = satdata.LE_meter(gobsL2, gobsL4);

		if (!double_eq(LE24_obs, 0.0))
		{
			_LE24[_crt_time][sat][1] = LE24_obs;
		}
		else
		{
			_LE24[_crt_time][sat][1] = 0.0;
		}
	}

	void t_gppplsq::_combineLE25(t_gsatdata& satdata)
	{
		string sat = satdata.sat();
		if (_crt_time == _beg_time)
		{
			_LE25[_beg_time - _obs_intv][sat][1] = 0.0;
		}

		//GOBSBAND b2 = t_gsys::band_priority(satdata.gsys(), FREQ_2);
		//GOBSBAND b3 = t_gsys::band_priority(satdata.gsys(), FREQ_3);
		GOBSBAND b2 = _band_index[satdata.gsys()][FREQ_2];
		GOBSBAND b5 = _band_index[satdata.gsys()][FREQ_5];
		t_gobs gobsL2 = t_gobs(satdata.select_phase(b2));
		t_gobs gobsL5 = t_gobs(satdata.select_phase(b5));

		double LE25_obs = satdata.LE_meter(gobsL2, gobsL5);

		if (!double_eq(LE25_obs, 0.0))
		{
			_LE25[_crt_time][sat][1] = LE25_obs;
		}
		else
		{
			_LE25[_crt_time][sat][1] = 0.0;
		}
	}

	void t_gppplsq::_combineLWL(t_gsatdata& satdata)
	{
		string sat = satdata.sat();
		if (_crt_time == _beg_time)
		{
			_LWL[_beg_time - _obs_intv][sat][1] = 0.0;
		}

		//GOBSBAND b1 = t_gsys::band_priority(satdata.gsys(), FREQ_1);
		//GOBSBAND b3 = t_gsys::band_priority(satdata.gsys(), FREQ_3);
		GOBSBAND b1 = _band_index[satdata.gsys()][FREQ_1];
		GOBSBAND b3 = _band_index[satdata.gsys()][FREQ_3];
		t_gobs gobsL1 = t_gobs(satdata.select_phase(b1));
		t_gobs gobsL3 = t_gobs(satdata.select_phase(b3));

		double factor13 = satdata.LWL_factor13(gobsL1, gobsL3);
		double LWL_obs;
		if (!double_eq(factor13, 0.0) && !double_eq(_LW[_crt_time][sat][1], 0.0) && !double_eq(_LE[_crt_time][sat][1], 0.0))
		{
			LWL_obs = 1 / (1 - factor13)*_LW[_crt_time][sat][1] - 1 / (1 / factor13 - 1)*_LE[_crt_time][sat][1];
		}
		else
		{
			LWL_obs = 0.0;
		}

		if (!double_eq(LWL_obs, 0.0))
		{
			_LWL[_crt_time][sat][1] = LWL_obs;
		}
		else
		{
			_LWL[_crt_time][sat][1] = 0.0;
		}
	}

	void t_gppplsq::_combineLWL124(t_gsatdata& satdata)
	{
		string sat = satdata.sat();
		if (_crt_time == _beg_time)
		{
			_LWL124[_beg_time - _obs_intv][sat][1] = 0.0;
		}

		//GOBSBAND b1 = t_gsys::band_priority(satdata.gsys(), FREQ_1);
		//GOBSBAND b3 = t_gsys::band_priority(satdata.gsys(), FREQ_3);
		GOBSBAND b1 = _band_index[satdata.gsys()][FREQ_1];
		GOBSBAND b4 = _band_index[satdata.gsys()][FREQ_4];
		t_gobs gobsL1 = t_gobs(satdata.select_phase(b1));
		t_gobs gobsL4 = t_gobs(satdata.select_phase(b4));

		double factor14 = satdata.LWL_factor13(gobsL1, gobsL4);
		double LWL124_obs;
		if (!double_eq(factor14, 0.0) && !double_eq(_LW[_crt_time][sat][1], 0.0) && !double_eq(_LE24[_crt_time][sat][1], 0.0))
		{
			LWL124_obs = 1 / (1 - factor14)*_LW[_crt_time][sat][1] - 1 / (1 / factor14 - 1)*_LE24[_crt_time][sat][1];
		}
		else
		{
			LWL124_obs = 0.0;
		}

		if (!double_eq(LWL124_obs, 0.0))
		{
			_LWL124[_crt_time][sat][1] = LWL124_obs;
		}
		else
		{
			_LWL124[_crt_time][sat][1] = 0.0;
		}
	}

	void t_gppplsq::_combineLWL125(t_gsatdata& satdata)
	{
		string sat = satdata.sat();
		if (_crt_time == _beg_time)
		{
			_LWL125[_beg_time - _obs_intv][sat][1] = 0.0;
		}

		//GOBSBAND b1 = t_gsys::band_priority(satdata.gsys(), FREQ_1);
		//GOBSBAND b3 = t_gsys::band_priority(satdata.gsys(), FREQ_3);
		GOBSBAND b1 = _band_index[satdata.gsys()][FREQ_1];
		GOBSBAND b5 = _band_index[satdata.gsys()][FREQ_5];
		t_gobs gobsL1 = t_gobs(satdata.select_phase(b1));
		t_gobs gobsL5 = t_gobs(satdata.select_phase(b5));

		double factor15 = satdata.LWL_factor13(gobsL1, gobsL5);
		double LWL125_obs;
		if (!double_eq(factor15, 0.0) && !double_eq(_LW[_crt_time][sat][1], 0.0) && !double_eq(_LE25[_crt_time][sat][1], 0.0))
		{
			LWL125_obs = 1 / (1 - factor15)*_LW[_crt_time][sat][1] - 1 / (1 / factor15 - 1)*_LE25[_crt_time][sat][1];
		}
		else
		{
			LWL125_obs = 0.0;
		}

		if (!double_eq(LWL125_obs, 0.0))
		{
			_LWL125[_crt_time][sat][1] = LWL125_obs;
		}
		else
		{
			_LWL125[_crt_time][sat][1] = 0.0;
		}
	}

	map<GSYS,double> t_gppplsq::_recWidelaneUPD(string site, std::vector<t_gsatdata>& crt_obs)
	{
		try
		{
			map<GSYS, double> result;
			set<GSYS> gsys = { GSYS::GPS,GSYS::GAL,GSYS::BDS,GSYS::QZS };
			for (auto sys = gsys.begin(); sys != gsys.end(); sys++)
			{
				double upd_rec = 0.0;
				map<string, double> xwl;
				map<string, double> wgt;
				list<pair<double, double>> wl_sig;
				double sigx = 0.0, mean = 0.0, sig = 0.0;
				int npnt = 0;

				std::vector<t_gsatdata>::iterator it_beg = crt_obs.begin();
				std::vector<t_gsatdata>::iterator it_end = crt_obs.end();
				string crt_sat("NONE");
				for (auto it = it_beg; it != it_end; it++)
				{
					if (_sat_list.empty() ||
						_sat_list.find(it->sat()) == _sat_list.end())
					{
						continue;
					}
					crt_sat = it->sat();
					if (it->gsys() != *sys) continue;

					if (_rMW[_crt_time][crt_sat][1] == 0.0 || _LWL[_crt_time][crt_sat][1] == 0.0) continue;

					xwl[crt_sat] = _rMW[_crt_time][crt_sat][1];
					wgt[crt_sat] = 1.0;
					wl_sig.push_back(make_pair(xwl[crt_sat], wgt[crt_sat]));
					npnt++;
				}
				if (npnt > 2)
				{
					getMeanFract(wl_sig, mean, sig, sigx);
					mean = getFraction(mean, -0.5);
					upd_rec = sigx <= 0.10 ? mean : 0.0;
				}
				result.insert({ *sys, upd_rec });
			}
			
			return result;
		}
		catch (...)
		{
			cout << "ERROR : t_gambiguity::_rec_wl_upd throw exception" << endl;
			throw(-1); //glfeng
		}
	}

	map<GSYS, double> t_gppplsq::_recExtraWidelaneUPD(string site, std::vector<t_gsatdata>& crt_obs, string mode)
	{
		try
		{
			map<GSYS, double> result;
			set<GSYS> gsys = { GSYS::GPS,GSYS::GAL,GSYS::BDS,GSYS::QZS };
			for (auto sys = gsys.begin(); sys != gsys.end(); sys++)
			{
				double upd_rec = 0.0;
				map<string, double> xewl;
				map<string, double> wgt;
				list<pair<double, double>> ewl_sig;
				double sigx = 0.0, mean = 0.0, sig = 0.0;
				int npnt = 0;
				std::vector<t_gsatdata>::iterator it_beg = crt_obs.begin();
				std::vector<t_gsatdata>::iterator it_end = crt_obs.end();
				string crt_sat("NONE");
				for (auto it = it_beg; it != it_end; it++)
				{
					if (_sat_list.empty() ||
						_sat_list.find(it->sat()) == _sat_list.end())
					{
						continue;
					}
					crt_sat = it->sat();
					if (it->gsys() != *sys) continue;

					if (mode == "EWL" && (_rEWL[_crt_time][crt_sat][1] == 0.0 || _LWL[_crt_time][crt_sat][1] == 0.0)) continue;
					if (mode == "EWL24" && (_rEWL24[_crt_time][crt_sat][1] == 0.0 || _LWL124[_crt_time][crt_sat][1] == 0.0)) continue;
					if (mode == "EWL25" && (_rEWL25[_crt_time][crt_sat][1] == 0.0 || _LWL125[_crt_time][crt_sat][1] == 0.0)) continue;

					if (mode == "EWL") xewl[crt_sat] = _rEWL[_crt_time][crt_sat][1];
					else if(mode=="EWL24") xewl[crt_sat] = _rEWL24[_crt_time][crt_sat][1];
					else if(mode == "EWL25") xewl[crt_sat] = _rEWL25[_crt_time][crt_sat][1];

					wgt[crt_sat] = 1.0;
					ewl_sig.push_back(make_pair(xewl[crt_sat], wgt[crt_sat]));
					npnt++;
				}
				if (npnt > 2)
				{
					getMeanFract(ewl_sig, mean, sig, sigx);
					mean = getFraction(mean, -0.5);
					upd_rec = sigx <= 0.10 ? mean : 0.0;
				}
				result.insert({ *sys, upd_rec });
			}

			return result;
		}
		catch (...)
		{
			cout << "ERROR : t_gambiguity::_rec_ewl_upd throw exception" << endl;
			throw(-1); //glfeng
		}
	}

	map<GSYS, double> t_gppplsq::_recWidelaneUPD_RAW(string site, std::vector<t_gsatdata>& crt_obs)
	{
		try
		{
			map<GSYS, double> result;
			set<GSYS> gsys = { GSYS::GPS,GSYS::GAL,GSYS::BDS,GSYS::QZS };
			for (auto sys = gsys.begin(); sys != gsys.end(); sys++)
			{
				double upd_rec = 0.0;
				map<string, double> xwl;
				map<string, double> wgt;
				list<pair<double, double>> wl_sig;
				double sigx = 0.0, mean = 0.0, sig = 0.0;
				int npnt = 0;

				std::vector<t_gsatdata>::iterator it_beg = crt_obs.begin();
				std::vector<t_gsatdata>::iterator it_end = crt_obs.end();
				string crt_sat("NONE");
				for (auto it = it_beg; it != it_end; it++)
				{
					if (_sat_list.empty() ||
						_sat_list.find(it->sat()) == _sat_list.end())
					{
						continue;
					}
					crt_sat = it->sat();
					if (it->gsys() != *sys) continue;

					GSYS sys = it->gsys();
					GOBSBAND b1 = _band_index[sys][FREQ_1];
					GOBSBAND b2 = _band_index[sys][FREQ_2];
					double wlength1 = it->wavelength(b1);
					double wlength2 = it->wavelength(b2);
					if (wlength1 == 0.0 || wlength2 == 0.0) continue;

					double amb_L1 = 0.0, amb_L2 = 0.0;
					t_gallpar pars = _ambfix->getFinalParams();
					int idx_L1 = pars.getParam(_crt_rec, t_gpar::AMB_L1, crt_sat);
					int idx_L2 = pars.getParam(_crt_rec, t_gpar::AMB_L2, crt_sat);
					if (idx_L1 != -1) amb_L1 = pars[idx_L1].value() / wlength1;
					if (idx_L2 != -1) amb_L2 = pars[idx_L2].value() / wlength2;

					double sat_wl_upd = 0.0;
					one_epoch_upd epoch_upd_wl = _gupd->get_epo_upd(UPDTYPE::WL, WL_IDENTIFY);
					if (epoch_upd_wl.find(crt_sat) == epoch_upd_wl.end()) continue;

					if (epoch_upd_wl[crt_sat]->sigma > 0.2 || epoch_upd_wl[crt_sat]->npoint <= 2)
					{
						continue;
					}
					sat_wl_upd = epoch_upd_wl[crt_sat]->value;

					if (amb_L1 == 0.0 || amb_L2 == 0.0 || _LWL[_crt_time][crt_sat][1] == 0.0) continue;

					xwl[crt_sat] = amb_L1 - amb_L2 - sat_wl_upd;
					wgt[crt_sat] = 1.0;
					wl_sig.push_back(make_pair(xwl[crt_sat], wgt[crt_sat]));
					npnt++;
				}
				if (npnt > 2)
				{
					getMeanFract(wl_sig, mean, sig, sigx);
					mean = getFraction(mean, -0.5);
					upd_rec = sigx <= 0.10 ? mean : 0.0;
				}
				result.insert({ *sys, upd_rec });
			}

			return result;
		}
		catch (...)
		{
			cout << "ERROR : t_gambiguity::_rec_wl_upd throw exception" << endl;
			throw(-1);
		}
	}

	map<GSYS, double> t_gppplsq::_recExtraWidelaneUPD_RAW(string site, std::vector<t_gsatdata>& crt_obs, string mode)
	{
		try
		{
			map<GSYS, double> result;
			set<GSYS> gsys = { GSYS::GPS,GSYS::GAL,GSYS::BDS,GSYS::QZS };
			for (auto sys = gsys.begin(); sys != gsys.end(); sys++)
			{
				double upd_rec = 0.0;
				map<string, double> xewl;
				map<string, double> wgt;
				list<pair<double, double>> ewl_sig;
				double sigx = 0.0, mean = 0.0, sig = 0.0;
				int npnt = 0;
				std::vector<t_gsatdata>::iterator it_beg = crt_obs.begin();
				std::vector<t_gsatdata>::iterator it_end = crt_obs.end();
				string crt_sat("NONE");
				for (auto it = it_beg; it != it_end; it++)
				{
					if (_sat_list.empty() ||
						_sat_list.find(it->sat()) == _sat_list.end())
					{
						continue;
					}
					crt_sat = it->sat();
					if (it->gsys() != *sys) continue;

					GSYS sys = it->gsys();
					GOBSBAND b2 = _band_index[sys][FREQ_2];
					GOBSBAND b3;
					if (mode == "EWL") b3 = _band_index[sys][FREQ_3];
					else if (mode == "EWL24") b3 = _band_index[sys][FREQ_4];
					else if (mode == "EWL25") b3 = _band_index[sys][FREQ_5];
					double wlength2 = it->wavelength(b2);
					double wlength3 = it->wavelength(b3);
					if (wlength2 == 0.0 || wlength3 == 0.0) continue;

					double amb_L2 = 0.0, amb_L3 = 0.0;
					t_gallpar pars = _ambfix->getFinalParams();
					int idx_L2 = pars.getParam(_crt_rec, t_gpar::AMB_L2, crt_sat);
					int idx_L3;
					if (mode == "EWL") idx_L3 = pars.getParam(_crt_rec, t_gpar::AMB_L3, crt_sat);
					else if (mode == "EWL24") idx_L3 = pars.getParam(_crt_rec, t_gpar::AMB_L4, crt_sat);
					else if (mode == "EWL25") idx_L3 = pars.getParam(_crt_rec, t_gpar::AMB_L5, crt_sat);
					if (idx_L2 != -1) amb_L2 = pars[idx_L2].value() / wlength2;
					if (idx_L3 != -1) amb_L3 = pars[idx_L3].value() / wlength3;

					double sat_ewl_upd = 0.0;
					one_epoch_upd epoch_upd_ewl;
					if (mode == "EWL") epoch_upd_ewl = _gupd->get_epo_upd(UPDTYPE::EWL, EWL_IDENTIFY);
					else if (mode == "EWL24") epoch_upd_ewl = _gupd->get_epo_upd(UPDTYPE::EWL24, EWL24_IDENTIFY);
					else if (mode == "EWL25") epoch_upd_ewl = _gupd->get_epo_upd(UPDTYPE::EWL25, EWL25_IDENTIFY);
					if (epoch_upd_ewl.find(crt_sat) == epoch_upd_ewl.end()) continue;

					if (epoch_upd_ewl[crt_sat]->sigma > 0.2 || epoch_upd_ewl[crt_sat]->npoint <= 2)
					{
						continue;
					}

					sat_ewl_upd = epoch_upd_ewl[crt_sat]->value;

					double LWL;
					if (mode == "EWL") LWL = _LWL[_crt_time][crt_sat][1];
					else if (mode == "EWL24") LWL = _LWL124[_crt_time][crt_sat][1];
					else if (mode == "EWL25") LWL = _LWL125[_crt_time][crt_sat][1];

					if (amb_L2 == 0.0 || amb_L3 == 0.0 || LWL == 0.0) continue;

					xewl[crt_sat] = amb_L2 - amb_L3 - sat_ewl_upd;
					wgt[crt_sat] = 1.0;
					ewl_sig.push_back(make_pair(xewl[crt_sat], wgt[crt_sat]));
					npnt++;
				}
				if (npnt > 2)
				{
					getMeanFract(ewl_sig, mean, sig, sigx);
					mean = getFraction(mean, -0.5);
					upd_rec = sigx <= 0.10 ? mean : 0.0;
				}
				result.insert({ *sys, upd_rec });
			}

			return result;
		}
		catch (...)
		{
			cout << "ERROR : t_gambiguity::_rec_ewl_upd throw exception" << endl;
			throw(-1); //glfeng
		}
	}

	int t_gppplsq::_addFixConstraintAFIF(const t_gtime& crt_epoch, std::vector<t_gsatdata>& crt_obs, string mode)
	{
		vector<double> codeOmc = _crt_equ.get_codeomc();
		map<GSYS, double> rec_wl_upd = _recWidelaneUPD(_crt_rec, crt_obs);
		map<GSYS, double> rec_ewl_upd = _recExtraWidelaneUPD(_crt_rec, crt_obs, mode);
		std::vector<t_gsatdata>::iterator it_beg = crt_obs.begin();
		std::vector<t_gsatdata>::iterator it_end = crt_obs.end();
		string crt_sat("NONE");
		int num = 0;
		auto sat_list = _crt_equ.get_sitelist(_crt_rec);
		for (auto it = it_beg; it != it_end; it++)
		{			
			crt_sat = it->sat();
			num = num + 1;
			
			if (_sat_list.empty() ||
				_sat_list.find(it->sat()) == _sat_list.end() || sat_list.find(it->sat()) == sat_list.end())
			{
				continue;
			}

			if (rec_wl_upd.find(it->gsys())->second == 0.0 || rec_ewl_upd.find(it->gsys())->second == 0.0)
				continue;

			double sat_wl_upd = 0.0, sat_ewl_upd = 0.0;

			if (_upd_mode == t_gsetamb::UPD_MODE::UPD)
			{
				one_epoch_upd epoch_upd_wl = _gupd->get_epo_upd(UPDTYPE::WL, WL_IDENTIFY);
				one_epoch_upd epoch_upd_ewl;
				if (mode == "EWL") epoch_upd_ewl = _gupd->get_epo_upd(UPDTYPE::EWL, EWL_IDENTIFY);
				else if (mode == "EWL24") epoch_upd_ewl = _gupd->get_epo_upd(UPDTYPE::EWL24, EWL24_IDENTIFY);
				else if (mode == "EWL25") epoch_upd_ewl = _gupd->get_epo_upd(UPDTYPE::EWL25, EWL25_IDENTIFY);

				if (epoch_upd_wl.find(crt_sat) == epoch_upd_wl.end()) continue;
				if (epoch_upd_ewl.find(crt_sat) == epoch_upd_ewl.end()) continue;

				if (epoch_upd_wl[crt_sat]->sigma > 0.2 || epoch_upd_wl[crt_sat]->npoint <= 2 ||
					epoch_upd_ewl[crt_sat]->sigma > 0.2 || epoch_upd_ewl[crt_sat]->npoint <= 2)
				{
					continue;
				}
				sat_wl_upd = epoch_upd_wl[crt_sat]->value;
				sat_ewl_upd = epoch_upd_ewl[crt_sat]->value;

				//epoch_upd_ewl = _gupd->get_epo_upd(UPDTYPE::EWL_epoch, crt_epoch);
				//if (epoch_upd_ewl.find(crt_sat) != epoch_upd_ewl.end())
				//{
				//	if (epoch_upd_ewl[crt_sat]->sigma <= 0.2 && epoch_upd_ewl[crt_sat]->npoint >= 2)
				//	{
				//		sat_ewl_upd = epoch_upd_ewl[crt_sat]->value;
				//	}
				//}
			}
			GSYS sys = it->gsys();
			double ratio_a = 1000;
			if (sys == GSYS::GAL) ratio_a = 5000;
			GOBSBAND b1 = _band_index[sys][FREQ_1];
			GOBSBAND b2 = _band_index[sys][FREQ_2];
			GOBSBAND b3;
			if (mode == "EWL") b3 = _band_index[sys][FREQ_3];
			else if (mode == "EWL24") b3 = _band_index[sys][FREQ_4];
			else if (mode == "EWL25") b3 = _band_index[sys][FREQ_5];
			double wlength1 = it->wavelength(b1);
			double wlength2 = it->wavelength(b2);
			double wlength3 = it->wavelength(b3);
			if (wlength1 == 0.0 || wlength2 == 0.0 || wlength3 == 0.0) continue;

			double LWL, rMW, rEWL;
			rMW = _rMW[_crt_time][crt_sat][1];
			if (mode == "EWL")
			{
				LWL = _LWL[_crt_time][crt_sat][1];
				rEWL = _rEWL[_crt_time][crt_sat][1];
			}
			else if (mode == "EWL24")
			{
				LWL = _LWL124[_crt_time][crt_sat][1];
				rEWL = _rEWL24[_crt_time][crt_sat][1];
			}
			else if (mode == "EWL25")
			{
				LWL = _LWL125[_crt_time][crt_sat][1];
				rEWL = _rEWL25[_crt_time][crt_sat][1];
			}

			if (LWL != 0.0 && rMW != 0.0 && rEWL != 0.0)
			{
				double res_wl = rMW - rec_wl_upd.find(it->gsys())->second;
				double res_ewl = rEWL - rec_ewl_upd.find(it->gsys())->second;
				if (abs(res_wl - round(res_wl)) < 0.25 && abs(res_ewl - round(res_ewl)) < 0.25)
				{
					_IMW[_crt_time][crt_sat][1] = round(res_wl);
					_IEWL[_crt_time][crt_sat][1] = round(res_ewl);
					_PW[_crt_time][crt_sat][1] = LWL +
						1 / (wlength3 / wlength1 - 1)*(wlength2*wlength3) / (wlength3 - wlength2)
						*(_IEWL[_crt_time][crt_sat][1] + rec_ewl_upd.find(it->gsys())->second + sat_ewl_upd)
						- 1 / (1 / wlength1 - 1 / wlength2)*(1 / (1 - wlength1 / wlength3))
						*(_IMW[_crt_time][crt_sat][1] + rec_wl_upd.find(it->gsys())->second + sat_wl_upd);
					GOBSBAND b1 = _band_index[sys][FREQ_1];
					GOBSBAND b2 = _band_index[sys][FREQ_2];
					t_gobs obsP1 = crt_obs[(int64_t)num - 1].select_range(b1);
					t_gobs obsP2 = crt_obs[(int64_t)num - 1].select_range(b2);
					double obs = crt_obs[(int64_t)num - 1].P3(obsP1, obsP2);
					//seek for crt_sat idx
					int idx;
					for (idx = 0; idx < _crt_equ.B.size(); idx++)
					{
						if (_crt_equ.get_satname(idx) == crt_sat)
							break;
					}
					_AFIF[_crt_time][crt_sat][1] = _PW[_crt_time][crt_sat][1] - obs + codeOmc[idx / 2];
					//add AFIF omc to equation
					vector<pair<int, double> > coef_AFIF = _crt_equ.B[idx];
					t_gobscombtype type(obsP1, _obs_mode);
					_crt_equ.add_equ(coef_AFIF, _crt_equ.P[idx], _AFIF[_crt_time][crt_sat][1], _crt_rec, crt_sat, type, false);
				}
			}
		}
		return 0;
	}

	int t_gppplsq::_addFixConstraintAFIF_RAW(const t_gtime& crt_epoch, t_glsq*  glsq, string mode)
	{
		//vector<double> codeOmc = _crt_equ.get_codeomc();
		vector<double> phaseOmc = _crt_equ.get_phaseomc();

		//double rec_wl_upd = _recWidelaneUPD(_crt_rec, _crt_obs_list);
		//double rec_ewl_upd = _recExtraWidelaneUPD(_crt_rec, _crt_obs_list);
		map<GSYS, double> rec_wl_upd = _recWidelaneUPD_RAW(_crt_rec, _crt_obs_list);
		map<GSYS, double> rec_ewl_upd = _recExtraWidelaneUPD_RAW(_crt_rec, _crt_obs_list, mode);
		std::vector<t_gsatdata>::iterator it_beg = _crt_obs_list.begin();
		std::vector<t_gsatdata>::iterator it_end = _crt_obs_list.end();
		string crt_sat("NONE");
		int num = 0;
		auto sat_list = _crt_equ.get_sitelist(_crt_rec);
		for (auto it = it_beg; it != it_end; it++)
		{
			crt_sat = it->sat();
			num = num + 1;
			
			if (_sat_list.empty() ||
				_sat_list.find(it->sat()) == _sat_list.end() || sat_list.find(it->sat()) == sat_list.end())
			{
				continue;
			}
			
			if (rec_wl_upd.find(it->gsys())->second == 0.0 || rec_ewl_upd.find(it->gsys())->second == 0.0)
				continue;

			double sat_wl_upd = 0.0, sat_ewl_upd = 0.0;

			if (_upd_mode == t_gsetamb::UPD_MODE::UPD&&!_isClient)
			{
				one_epoch_upd epoch_upd_wl = _gupd->get_epo_upd(UPDTYPE::WL, WL_IDENTIFY);
				one_epoch_upd epoch_upd_ewl;
				if (mode == "EWL") epoch_upd_ewl = _gupd->get_epo_upd(UPDTYPE::EWL, EWL_IDENTIFY);
				else if (mode == "EWL24") epoch_upd_ewl = _gupd->get_epo_upd(UPDTYPE::EWL24, EWL24_IDENTIFY);
				else if (mode == "EWL25") epoch_upd_ewl = _gupd->get_epo_upd(UPDTYPE::EWL25, EWL25_IDENTIFY);
				if (epoch_upd_wl.find(crt_sat) == epoch_upd_wl.end()) continue;
				if (epoch_upd_ewl.find(crt_sat) == epoch_upd_ewl.end()) continue;

				if (epoch_upd_wl[crt_sat]->sigma > 0.2 || epoch_upd_wl[crt_sat]->npoint <= 2 ||
					epoch_upd_ewl[crt_sat]->sigma > 0.2 || epoch_upd_ewl[crt_sat]->npoint <= 2)
				{
					continue;
				}
				sat_wl_upd = epoch_upd_wl[crt_sat]->value;
				sat_ewl_upd = epoch_upd_ewl[crt_sat]->value;

				epoch_upd_ewl = _gupd->get_epo_upd(UPDTYPE::EWL_epoch, crt_epoch);
				if (epoch_upd_ewl.find(crt_sat) != epoch_upd_ewl.end())
				{
					if (epoch_upd_ewl[crt_sat]->sigma <= 0.2 && epoch_upd_ewl[crt_sat]->npoint >= 2)
					{
						sat_ewl_upd = epoch_upd_ewl[crt_sat]->value;
					}
				}


			}
			GSYS sys = it->gsys();
			double ratio_a = 1000;
			if (sys == GSYS::GAL) ratio_a = 5000;
			GOBSBAND b1 = _band_index[sys][FREQ_1];
			GOBSBAND b2 = _band_index[sys][FREQ_2];
			GOBSBAND b3;
			if (mode == "EWL") b3 = _band_index[sys][FREQ_3];
			else if (mode == "EWL24") b3 = _band_index[sys][FREQ_4];
			else if (mode == "EWL25") b3 = _band_index[sys][FREQ_5];
			double wlength1 = it->wavelength(b1);
			double wlength2 = it->wavelength(b2);
			double wlength3 = it->wavelength(b3);
			if (wlength1 == 0.0 || wlength2 == 0.0 || wlength3 == 0.0) continue;

			double amb_L1 = 0.0, amb_L2 = 0.0, amb_L3 = 0.0;
			t_gallpar pars = _ambfix->getFinalParams();
			int idx_L1 = pars.getParam(_crt_rec, t_gpar::AMB_L1, crt_sat);
			int idx_L2 = pars.getParam(_crt_rec, t_gpar::AMB_L2, crt_sat);
			int idx_L3;
			if (mode == "EWL") idx_L3 = pars.getParam(_crt_rec, t_gpar::AMB_L3, crt_sat);
			else if (mode == "EWL24") idx_L3 = pars.getParam(_crt_rec, t_gpar::AMB_L4, crt_sat);
			else if (mode == "EWL25") idx_L3 = pars.getParam(_crt_rec, t_gpar::AMB_L5, crt_sat);
			if (idx_L1 != -1) amb_L1 = pars[idx_L1].value() / wlength1;
			if (idx_L2 != -1) amb_L2 = pars[idx_L2].value() / wlength2;
			if (idx_L3 != -1) amb_L3 = pars[idx_L3].value() / wlength3;

			double LWL;
			if (mode == "EWL") LWL = _LWL[_crt_time][crt_sat][1];
			else if (mode == "EWL24") LWL = _LWL124[_crt_time][crt_sat][1];
			else if (mode == "EWL25") LWL = _LWL125[_crt_time][crt_sat][1];

			if (amb_L1 != 0.0 && amb_L2 != 0.0 && amb_L3 != 0.0 && LWL != 0.0 )
			{
				double res_wl = amb_L1 - amb_L2 - rec_wl_upd.find(it->gsys())->second - sat_wl_upd;
				double res_ewl = amb_L2 - amb_L3 - rec_ewl_upd.find(it->gsys())->second - sat_ewl_upd;
				if (abs(res_wl - round(res_wl)) < 0.25 && abs(res_ewl - round(res_ewl)) < 0.25)
				{
					_IMW[_crt_time][crt_sat][1] = round(res_wl);
					_IEWL[_crt_time][crt_sat][1] = round(res_ewl);
					_PW[_crt_time][crt_sat][1] = LWL +
						1 / (wlength3 / wlength1 - 1)*(wlength2*wlength3) / (wlength3 - wlength2)
						*(_IEWL[_crt_time][crt_sat][1] + rec_ewl_upd.find(it->gsys())->second + sat_ewl_upd)
						- 1 / (1 / wlength1 - 1 / wlength2)*(1 / (1 - wlength1 / wlength3))
						*(_IMW[_crt_time][crt_sat][1] + rec_wl_upd.find(it->gsys())->second + sat_wl_upd);
					GOBSBAND b1 = _band_index[sys][FREQ_1];
					GOBS obsP1 = _crt_obs_list[(int64_t)num - 1].select_range(b1);
					//double C1 = it->getobs(obsP1);
					GOBS obsL1 = _crt_obs_list[(int64_t)num - 1].select_phase(b1);
					double L1 = it->getobs(obsL1);
				
					//seek for crt_sat idx
					int idx;
					for (idx = 0; idx < _crt_equ.B.size(); idx++)
					{
						if (_crt_equ.get_satname(idx) == crt_sat)
							break;
					}
					//_AFIF[_crt_time][crt_sat][1] = _PW[_crt_time][crt_sat][1] - C1 + codeOmc[idx / 2];
					_AFIF[_crt_time][crt_sat][1] = _PW[_crt_time][crt_sat][1] - L1 * wlength1 + phaseOmc[idx / 2] + pars[idx_L1].value();
				
				//	cerr << "AFIF value" << _AFIF[_crt_time][crt_sat][1] << endl;
					//add AFIF omc to equation
					vector<pair<int, double> > coef_AFIF = _crt_equ.B[idx];
					t_glsqEquationMatrix virtual_equ;
					t_gobscombtype type(obsP1, _obs_mode);
					virtual_equ.add_equ(coef_AFIF, _crt_equ.P[idx], _AFIF[_crt_time][crt_sat][1], _crt_rec, crt_sat, type, false);
					glsq->add_equation(virtual_equ);
				}
			}
		}
		try
		{
			for (int i = 0; i < glsq->NEQ().Nrows(); i++)
			{
				glsq->_x_solve[i].apriori(0.0);
			}
			glsq->solve_NEQ();
		}
		catch (exception e)
		{
			write_log_info(_glog, 1, e.what(), "Solve Equation Fail!");
			return false;
		}
		auto tmp_dx = glsq->dx();
		auto tmp_stdx = glsq->stdx();
		_ambfix->updateFixParam(glsq->_x_solve, tmp_dx, &tmp_stdx);

		cout << _crt_time.str_ymdhms() << endl;
		for (int i = 0; i < glsq->_x_solve.parNumber(); i++) {
			cout << setw(20) << "EPO  " << setw(20) << glsq->_x_solve[i].str_type() + "  " <<
				setw(20) << setprecision(5) << glsq->_x_solve[i].value() <<
				setw(15) << setprecision(5) << glsq->dx()(i + 1) <<
				setw(20) << setprecision(5) << glsq->_x_solve[i].value() + glsq->dx()(i + 1) <<
				setw(12) << setprecision(5) << glsq->stdx(i + 1) << endl;
		}
		return 0;
	}

	bool t_gppplsq::_getSatRef()
	{
		_sat_ref.clear();

		one_epoch_upd wl_upd, nl_upd, ewl_upd;
		if (_gupd)
		{
			wl_upd = _gupd->get_epo_upd(UPDTYPE::WL, _wl_Upd_time);
			ewl_upd = _gupd->get_epo_upd(UPDTYPE::EWL, _ewl_Upd_time);
			if(_upd_mode == t_gsetamb::UPD_MODE::UPD) nl_upd = _gupd->get_epo_upd(UPDTYPE::NL, _crt_time);
		}

		vector<t_gpar> params_ALL = _lsq->_x_solve.getAllPar();

		for (auto itsys = _sys_list.begin(); itsys != _sys_list.end(); itsys++)
		{
			GSYS sys = t_gsys::str2gsys(*itsys);
			double max_mw = 0.0;
			string ref;

			for (auto it_par = params_ALL.begin(); it_par != params_ALL.end(); it_par++)
			{
				if (it_par->str_type().find("AMB") == string::npos)continue;

				if (t_gsys::sat2gsys(it_par->prn) != sys)continue;

				if (sys == GSYS::BDS && stoi(it_par->prn.substr(1, 3)) <= 5)continue; // exclude BDS C01~C05

				// reference satellite must have WL/NL UPD
				//LX 注释窄巷因为整数钟时没有窄巷upd
				// elev >= 15 angle && continuous observation epoch > 600s
				// modified 0->300 for Server !!!!!
				if (_upd_mode == t_gsetamb::UPD_MODE::UPD)
				{
					if (_gupd && _frequency >= 3 && (ewl_upd.find(it_par->prn) == ewl_upd.end() || ewl_upd[it_par->prn]->npoint <= 2)) continue;
					if (_gupd && !(wl_upd.find(it_par->prn) != wl_upd.end() && nl_upd.find(it_par->prn) != nl_upd.end())) continue;
					if (_gupd && !(wl_upd[it_par->prn]->npoint > 2)) continue;
				}
				else if (_upd_mode == t_gsetamb::UPD_MODE::IRC)
				{
					if (_gupd && !(wl_upd.find(it_par->prn) != wl_upd.end()))continue;
				}

				if(_crt_ele[it_par->prn] < 15 || !_isClient && (_crt_time - it_par->beg <= 0)) continue;

				if (_obs_mode == RAW_MIX && _single_freq_sats.find(it_par->prn) != _single_freq_sats.end()) continue;

				double temp = (_crt_time - it_par->beg + 1) * _crt_ele[it_par->prn];
				// max  elev * continuous observation epoch 
				if (temp > max_mw)
				{
					max_mw = temp;
					ref = it_par->prn;
				}
			}
		//	cerr << "ref_sat    " << (_gupd->get_epo_upd(EWL_IDENTIFY).find(ref) != _gupd->get_epo_upd(EWL_IDENTIFY).end()) << endl;
		//	cerr << "ref_sat    " << ref << endl;
			if (!ref.empty()){
				_sat_ref.insert(ref);
			}
		}
		if (_sat_ref.empty())return false;
		else return true;
	}
}