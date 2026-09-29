/**
*
* @verbatim
History
-1.0 zhangwei  2019-03-18  creat
-1.1 Qinyujie  2019-03-20  optim
-1.2 Qinyujie  2019-04-11  add Doxygen style code
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gorblsq_leo.cpp
* @brief		lsq processing for leo
*
* @author       zhangwei, Wuhan University
* @version		1.0.0
* @date		    2019-03-18
*
*/
#include "gproc/gpodlsqIF.h"
#include "gutils/gleo.h"
#include "gset/gsetpar.h"
#include "gset/gsetout.h"
#include "gset/gsetinp.h"
#include "gutils/ginfolog.h"
#include "gutils/gturboedit.h"
#include "gio/gfile.h"
#include "gcoders/icsnav.h"
#include "gcoders/icsleo.h"
#include "gutils/gstring.h"
#include "gproc/gupdateparSLR.h"
#include <algorithm>
#include "gcoders/rinexo.h"


namespace great
{
	t_gpodlsqIF::t_gpodlsqIF()
	{

	}

	/** @brief default constructor. */
	t_gpodlsqIF::t_gpodlsqIF(t_gsetbase* set, t_gallproc* data, t_glog* log) :
		t_glsqproc(set, data, log),
		_result(0)
	{

		string tmp;
		tmp = dynamic_cast<t_gsetout*>(set)->outputs("flt");
		if (!tmp.empty())
		{
			substitute(tmp, "$(rec)", _crt_rec, false);
			_result = new t_giof;
			_result->tsys(t_gtime::GPS);
			_result->mask(tmp);
			_result->append(dynamic_cast<t_gsetout*>(set)->append());
		}

		if (_slr)
		{
			if (_slr_model)
			{
				delete _slr_model;
				_slr_model = nullptr;
			}
			_slr_model = new t_gslrmodel(data, _glog, _gset);

			_set_correlation_out();
		}
	}

	/** @brief default destructor. */
	t_gpodlsqIF::~t_gpodlsqIF()
	{

		if (_result)
		{
			if (_result->is_open()) _result->close();
			delete _result;
			_result = nullptr;
		}
		if (_slr_model)
		{
			delete _slr_model;
			_slr_model = nullptr;
		}
	}

	bool t_gpodlsqIF::ProcessBatch(t_gallproc * data, const t_gtime & beg, const t_gtime & end)
	{
		if (data == nullptr)
		{
			_glog->logError("t_gpodlsqIF", "ProcessBatch", "t_gallproc is nullptr");
			return false;
		}

		if (beg > end)
		{
			_glog->logError("t_gpodlsqIF", "ProcessBatch", "beg > end");
			_glog->logError("t_gpodlsqIF", "ProcessBatch", beg.str_mjdsod("beg is ==> "));
			_glog->logError("t_gpodlsqIF", "ProcessBatch", end.str_mjdsod("end is ==> "));
			return false;
		}

		// get the time
		_beg_time = (_beg_time > beg) ? beg : _beg_time;
		_end_time = (_end_time < end) ? end : _end_time;
		_crt_time = beg;
		_glog->logInfo("t_gpodlsqIF", "ProcessBatch", _beg_time.str_mjdsod("_beg_time is ==> "));
		_glog->logInfo("t_gpodlsqIF", "ProcessBatch", _end_time.str_mjdsod("_end_time is ==> "));
		_glog->logInfo("t_gpodlsqIF", "ProcessBatch", _crt_time.str_mjdsod("_crt_time is ==> "));
		//Init data (include SITE/GNSS/LEO/SLR data)
		bool init_data = _initLsqProcData(data);
		if (!init_data)
		{
			_glog->logError("t_gpodlsqIF", "ProcessBatch", "can not init data");
			return false;
		}
		//Init pars (include SITE/GNSS/LEO/SLR pars)
		bool init_pars = _initLsqProcPars(_lsq);
		if (!init_pars)
		{
			_glog->logError("t_gpodlsqIF", "ProcessBatch", "can not init pars");
			return false;
		}

		READ_OFILE_MODE rom = dynamic_cast<t_gsetproc*>(_gset)->read_ofile_mode();
		int step = dynamic_cast<t_gsetproc*>(_gset)->read_ofile_step();
		int cmlstep = step;
		t_gallobs*  obs=new t_gallobs();  obs->glog(_glog); obs->gset(_gset);
		while (_crt_time <= _end_time)
		{
			chrono::high_resolution_clock::time_point beg_t = chrono::high_resolution_clock::now();
			_glog->logInfo("t_gpodlsqIF", "ProcessBatch", _crt_time.str_mjdsod("PROCESS ==> "));
			
			if ((rom == READ_OFILE_MODE::REALTIME) && (_crt_time > _beg_time + cmlstep))
			{
				_readOfiles(step, obs);
				cmlstep += step;
			}
			_glog->logInfo("t_gpodlsqIF", "ProcessBatch", _crt_time.str_mjdsod("PROCESS ==> "));
			cerr << _crt_time.str_mjdsod("PROCESS ==> ") << endl;

			/* get observation */
			vector<t_gsatdata> crt_obs = _gall_obs->obs(_rec_list, _crt_time);
			
			bool epoch_valid = _processOneEpoch(_crt_time, crt_obs);
			if (!epoch_valid)
			{
				write_log_info(_glog, 1, "ERROR" + _crt_time.str_ymdhms("Processing failed in epoch "));
				_crt_time = _crt_time + _obs_intv;
				continue;
			}
			if (_lsq->mode() == gnut::EPO && _crt_time < _end_time)
			{
				_get_epoch_solution();
			}
			chrono::high_resolution_clock::time_point end_t = chrono::high_resolution_clock::now();
			writeLogInfo(_glog, 0, "NOTE", "Process Epoch ==>" + _crt_time.str_ymdhms() + "successed.");
			cerr << "Totally Spends time is " << chrono::duration_cast<chrono::milliseconds>(end_t - beg_t).count() / 1000.0 << " sec " << endl;
			_crt_time = _crt_time + _obs_intv;
		}

		// SLR parts
		if (_slr)
		{
			vector<t_gsatdata> crt_obs_none;
			_change_update_mode(data,_lsq->_x_solve);

			set<string>::iterator rec_beg = _rec_list.begin();
			set<string>::iterator rec_end = _rec_list.end();
			map<t_gtime, t_gpreprocslr> allonesat;
			map<string, map<t_gtime, t_gpreprocslr>>* all = _slr_model->_allpreproc->getallSLROBS();
			for (set<string>::iterator it = rec_beg; it != rec_end; it++)
			{
				string crt_rec = *it;
				auto key = all->find(strs2strl(crt_rec));
				if (key == all->end())
				{
					write_log_info(_glog, 1, "ERROR", "  NO SLR obs for " + strs2strl(crt_rec));
				}
				else
				{
					allonesat = (*all)[strs2strl(crt_rec)];
				}

				for (auto iter = allonesat.begin(); iter != allonesat.end(); iter++)
				{
					//bool updata_valid = _lsq->update_parameter(iter->first, crt_obs_none, false);
					//if (!updata_valid)
					//{
						//write_log_info(_glog, 1, iter->first.str(), ": updating slr parameter failed ");
						//return false;
					//}

					_slr_model->_cmb_equ(iter->first, _lsq->_x_solve, crt_rec, _crt_equ, iter->second._sta_name);

					_lsq->add_equation(_crt_equ, iter->first);
				}
			}
		}

		try
		{
			_lsq->lsq_sysbias_constraint();
			_lsq->solve_NEQ();
		}
		catch (exception e)
		{
			write_log_info(_glog, 1, e.what(), "Solve Equation Fail!");
			return false;
		}

		/* recover equation */
		try
		{
			//_get_correlation_matrix();
			_lsq->get_result_parameter(*_gall_rcv);
			_lsq->recover_parameters(*_gall_rcv);
			/*bool prod_valid = _initLsqProdData(nullptr);
			if (!prod_valid) return false;*/
		}
		catch (exception e)
		{
			write_log_info(_glog, 1, e.what(), "Recover Parameter Fail!");
			return false;
		}

		return true;
	}

	void t_gpodlsqIF::_get_correlation_matrix()
	{
		SymmetricMatrix Q_temp;
		Q_temp = _lsq->Qx();
		vector<string> namelist;
		vector<t_gpar> allpars = _lsq->_x_solve.getAllPar();
		int i = 1;

		auto par_end = allpars.end();
		for (auto iter = allpars.begin(); iter != par_end; iter++)
		{
			if (iter->parType == t_gpar::AMB_IF || iter->parType == t_gpar::CLK)
			{
				Matrix_remRC(Q_temp, i, i);
				iter = allpars.erase(iter);
				i--;
				iter--;
			}
			else if (iter->parType == t_gpar::DUT1)
			{
				string type_temp = iter->str_type();
				namelist.push_back(iter->str_type());
			}
			else if (iter->parType != t_gpar::RB && (iter + 1) != par_end && iter->parType == (iter + 1)->parType)
			{
				Matrix_remRC(Q_temp, i, i);
				iter = allpars.erase(iter);
				iter--;
				i--;
			}
			else
			{
				string type_temp = iter->str_type();
				if (iter->parType == t_gpar::RB || iter->parType == t_gpar::CRD_X || iter->parType == t_gpar::CRD_Y || iter->parType == t_gpar::CRD_Z)
				{
					namelist.push_back(iter->site + "_" + iter->str_type());
				}
				else
				{
					namelist.push_back(iter->str_type());
				}

			}
			i++;
		}

		SymmetricMatrix  _Qnew(Q_temp.Nrows());
		for (int i = 1; i <= Q_temp.Nrows(); i++)
		{
			for (int j = 1; j <= Q_temp.Ncols(); j++)
			{
				_Qnew(i, j) = Q_temp(i, j) / sqrt(Q_temp(i, i)*Q_temp(j, j));
			}
		}

		//cout << _Qnew;
		if (_corr_map) {
			_extract_corr_map(_Qnew, namelist);
			//return true;
		}

	}

	void t_gpodlsqIF::_extract_corr_map(SymmetricMatrix& Qnew, vector<string>& name)
	{
		ostringstream os;
		for (int i = 1; i <= Qnew.Nrows(); i++)
		{
			os << setw(12) << left << name[i - 1];
			for (int j = 1; j <= Qnew.Ncols(); j++)
			{
				os << setw(15) << setiosflags(ios::fixed) << setprecision(6) << right << Qnew(i, j);
			}
			os << endl;

		}

		if (_corr_map) {
			_corr_map->write(os.str().c_str(), os.str().size());
			_corr_map->flush();
		}
		else {
			_glog->comment(2, "t_gpodlsqIF::_printout", "have no output file!");
		}

		return;

	}


	void t_gpodlsqIF::_change_update_mode(t_gallproc* eop_data,t_gallpar& pars)
	{
		shared_ptr<t_gupdatepar> updateSLR(new t_gupdateparSLR(eop_data));
		//updateSLR->set_interval(_obs_intv);
		_lsq->set_update_par(updateSLR);

		int ix = pars.getParam("", t_gpar::XPOLE, "");
		int iy = pars.getParam("", t_gpar::YPOLE, "");
		int iu = pars.getParam("", t_gpar::UT1, "");
		int idx = pars.getParam("", t_gpar::DXPOLE, "");
		int idy = pars.getParam("", t_gpar::DYPOLE, "");
		int idu = pars.getParam("", t_gpar::DUT1, "");

		bool judge = true;
		pars.setParRemove(ix, judge);
		pars.setParRemove(iy, judge);
		pars.setParRemove(iu, judge);
		pars.setParRemove(idx, judge);
		pars.setParRemove(idy, judge);
		pars.setParRemove(idu, judge);

	}

	bool t_gpodlsqIF::GenerateProduct()
	{
		bool prod_valid = _initLsqProdData(nullptr);
		if (!prod_valid) return false;
		return true;
	}

	bool t_gpodlsqIF::_initLsqProdData(t_gallprod * data)
	{
		// The res file
		_extract_resfile(*_gall_rcv);
		_extract_stafile(*_gall_rcv);

		// update rec clk file
		_extract_clkfile(*_gall_rcv, t_gallprec::AR);

		//LEO parts
		if (_mode.count("D") > 0)
		{
			_update_leo_ics();
			// Update ICS-LEO file
			t_gio*    gout       = new t_gfile; 
			t_gcoder* gics_coder = new t_ics_leo(_gset, "", 4096);

			string ics_path = dynamic_cast<t_gsetout*>(_gset)->outputs("icsleo");
			if (ics_path == "")
			{
				ics_path = (dynamic_cast<t_gsetinp*>(_gset)->inputs("icsleo"))[0];
			}
			gout->path(ics_path); gout->glog(_glog);
			gics_coder->clear(); gics_coder->glog(_glog); gics_coder->path(ics_path);
			gics_coder->add_data("ID1", _gall_leo);
			gics_coder->add_data("OBJ", _gall_obj);
			gout->coder(gics_coder);
			gout->run_write();
			if (gout)
			{
				delete gout;
				gout = nullptr;
			}

			if (gics_coder)
			{
				delete gics_coder;
				gics_coder = nullptr;
			}
		}

		//GNSS parts
		if (_mode.count("S") > 0)
		{
			//sat clk
			_extract_clkfile(*_gall_rcv, t_gallprec::AS);

			_update_nav_ics();
			// Update ICS-GNSS file
			t_gio*          gout = new t_gfile; 
			t_gcoder* gics_coder = new t_ics_nav(_gset, "", 4096);
			string ics_path = dynamic_cast<t_gsetout*>(_gset)->outputs("ics");
			if (ics_path == "")
			{
				ics_path = (dynamic_cast<t_gsetinp*>(_gset)->inputs("ics"))[0] + "_out";
			}
			gout->path(ics_path); gout->glog(_glog);
			gics_coder->clear(); gics_coder->glog(_glog); gics_coder->path(ics_path);
			gics_coder->add_data("ID1", _gall_gns);
			gics_coder->add_data("OBJ", _gall_obj);
			gout->coder(gics_coder);
			gout->run_write();

			if (gout)
			{
				delete gout;
				gout = nullptr;
			}

			if (gics_coder)
			{
				delete gics_coder;
				gics_coder = nullptr;
			}
		}

		//SLR parts
		if (_slr)
		{
			_extract_rbfile(*_gall_rcv);
			_extract_stafile(*_gall_rcv);
			_extract_eopfile(*_gall_rcv);
			_extract_polefile(*_gall_rcv);
		}

		// update xml
		//_extract_xmlcrd(*_allrcv);

		return true;
	}

	bool t_gpodlsqIF::_initLsqProcPars(t_glsq* lsq)
	{
		if (_rec_list.empty() || _sat_list.empty())
		{
			write_log_info(_glog, 1, "ERROR", "no rec or no sat in POD_LSQ.");
			return false;
		}
		//Common parts
		bool rec_valid = _init_rec_clk_pars(lsq);
		bool sys_valid = _init_sys_bia_pars(lsq);

		bool ifb_valid = false;
		     if (_frequency != 2 && _ifb_model == IFB_MODEL::EST_IFB)     ifb_valid = _init_sys_ifb_pars(lsq);
		else if (_frequency != 2 && _ifb_model == IFB_MODEL::EST_REC_IFB) ifb_valid = _init_sys_rec_ifb_pars(lsq);
		else if (_frequency != 2 && _ifb_model == IFB_MODEL::EST_SAT_IFB) ifb_valid = _init_sys_sat_ifb_pars(lsq);
		else if (_frequency == 3 && _obs_mode == IONO_FREE && _ifb_model == IFB_MODEL::NONE)        ifb_valid = _init_rec_clk13_pars(lsq);
		else ifb_valid = true;

		if (!rec_valid)
		{
			write_log_info(_glog, 1, "ERROR", "Init the rec clk pars failed");
			return false;
		}
		if (!ifb_valid)
		{
			write_log_info(_glog, 1, "ERROR", "Init the ifb pars failed");
			return false;
		}
		if (!sys_valid)
		{
			write_log_info(_glog, 1, "ERROR", "Init the sys bias pars failed");
			return false;
		}
		//SLR parts
		if (_slr)
		{
			bool init_slrpars = _initSLRProcPars(_lsq);
			if (!init_slrpars)
			{
				write_log_info(_glog, 1, "ERROR", "Init the SLR pars failed");
				return false;
			}
		}
		//LEO parts
		if (_mode.count("D") > 0)
		{
			bool icsleo_valid = _init_leo_ics_pars(lsq);
			if (!icsleo_valid)
			{
				write_log_info(_glog, 1, "ERROR", "Init the LEO ICS pars failed");
				return false;
			}
		}
		//GNSS parts
		if (_mode.count("S") > 0)
		{
			bool ics_valid = _init_gns_ics_pars(lsq);
			bool erp_valid = _init_crs_erp_pars(lsq);
			bool crd_valid = _init_rec_crd_pars(lsq, false);
			bool trop_valid = _init_rec_trop_pars(lsq);
			bool clk_valid = _init_sat_clk_pars(lsq);
			if (!ics_valid)
			{
				write_log_info(_glog, 1, "ERROR", "Init the GNSS ICS pars failed");
				return false;
			}
			if (!erp_valid)
			{
				write_log_info(_glog, 1, "ERROR", "Init the ERP pars failed");
				return false;
			}
			if (!crd_valid)
			{
				write_log_info(_glog, 1, "ERROR", "Init the CRD pars failed");
				return false;
			}
			if (!trop_valid)
			{
				write_log_info(_glog, 1, "ERROR", "Init the TROP pars failed");
				return false;
			}
			if (!clk_valid)
			{
				write_log_info(_glog, 1, "ERROR", "Init the sat clk pars failed");
				return false;
			}
		}
		return true;
	}

	bool t_gpodlsqIF::_initLsqProcData(t_gallproc * data)
	{
		gtrace("t_gpodlsqIF::_initLsqProcData");

		// prepare obs data
		_gall_bias = dynamic_cast<t_gallbias*>((*data)[t_gdata::ALLBIAS]);
		_gall_obs = dynamic_cast<t_gallobs*>((*data)[t_gdata::ALLOBS]);
		_gall_obj = dynamic_cast<t_gallobj*>((*data)[t_gdata::ALLOBJ]);
		_gall_rcv = dynamic_cast<t_gallrecover*>((*data)[t_gdata::ALLRECOVER]);
		_gdata_erp = dynamic_cast<t_gpoleut1*>((*data)[t_gdata::ALLPOLEUT1]);
		_gdata_ambcon = dynamic_cast<t_gambcon*>((*data)[t_gdata::AMBCON]);
		//lijie for glonass fid
		_gdata_satinfo= dynamic_cast<t_gallsatinfo*>((*data)[t_gdata::SATPARS]); 
		if (!_gall_bias)
		{
			_glog->logInfo("t_gpodlsqIF", "_initLsqProcData", "ALLBIAS is empty!");
			// jdhuang
			// 1 + 7 无电离层的gal我们不改dcb
			cerr << "t_gpodlsqIF  _initLsqProcData  ALLBIAS is empty!" << endl;
			return false;
		}
		if (!_gall_obs)
		{
			_glog->logInfo("t_gpodlsqIF", "_initLsqProcData", "ALLOBS is empty!");
			return false;
		}
		if (!_gall_obj)
		{
			_glog->logInfo("t_gpodlsqIF", "_initLsqProcData", "ALLOBJ is empty!");
			return false;
		}
		if (!_gall_rcv)
		{
			_glog->logInfo("t_gpodlsqIF", "_initLsqProcData", "ALLRECOVER is empty!");
			return false;
		}
		if (!_gdata_erp)
		{
			_glog->logInfo("t_gpodlsqIF_leo", "_initLsqProcData", "ALLPOLEUT1 is empty!");
			return false;
		}
		if (_mode.count("S") > 0) // GPS parts (GPS POD or Integrated POD)
		{
			_gall_gns = dynamic_cast<t_gallics_nav*>((*data)[t_gdata::ALLICS]);
			if (!_gall_gns)
			{
				_glog->logInfo("t_gpodlsqIF", "_initLsqProcData", "ALLICS is empty!");
				return false;
			}
		}
		if (_mode.count("D") > 0) // LEO parts (LEO POD or Integrated POD)
		{
			_gall_leo = dynamic_cast<t_gallics_leo*>((*data)[t_gdata::ALLICSLEO]);
			if (!_gall_leo)
			{
				_glog->logInfo("t_gpodlsqIF", "_initLsqProcData", "ALLICSLEO is empty!");
				return false;
			}
		}
		if (_slr)  // With SLR observations
		{
			_gdata_slrf = dynamic_cast<t_gslrf*>((*data)[t_gdata::SLRF]);
			_gall_preproc = dynamic_cast<t_gallpreproc*>((*data)[t_gdata::ALLPREPROC]);
			_gall_rb = dynamic_cast<t_gallrb*>((*data)[t_gdata::ALLRB]);
			if (_gall_rb == nullptr)
			{
				_gall_rb = new t_gallrb();
			}
			
			if (!_gdata_slrf)
			{
				_glog->logInfo("t_gpodlsqIF", "_initLsqProcData", "SLRF is empty!");
				return false;
			}
			if (!_gall_preproc)
			{
				_glog->logInfo("t_gpodlsqIF", "_initLsqProcData", "ALLPREPROC is empty!");
				return false;
			}
			/*if (!_allrb)
			{
				_glog->logInfo("t_gpodlsqIF", "_initLsqProcData", "ALLRB is empty!");
				return false;
			}*/
			bool init_slrdata = _InitSLRFData(_slr_model->_allslrf);
			if (!init_slrdata)
			{
				_glog->logInfo("t_gpodlsqIF", "_initLsqProcData", "Init slr data");
				return false;
			}
		}

		// GPS parts
		if (_mode.count("S") > 0)
		{
			set<string> sat_list;
			set<string>::const_iterator sat_beg = _sat_list.begin();
			set<string>::const_iterator sat_end = _sat_list.end();
			set<string>::const_iterator ics_beg = _gall_gns->sat_list.begin();
			set<string>::const_iterator ics_end = _gall_gns->sat_list.end();
			set_intersection(sat_beg, sat_end, ics_beg, ics_end,
				inserter(sat_list, sat_list.begin()));
			_sat_list = sat_list;
			if (_sat_list.empty())
			{
				_glog->logInfo("t_gpodlsqIF", "_initLsqProcData", "_sat_list is empty!");
				return false;
			}
		}
		// LEO parts
		if (_mode.count("D") > 0)
		{
			set<string> sites; // ground sites
			set_difference(_rec_list.begin(), _rec_list.end(), _leo_list.begin(), _leo_list.end(), inserter(sites, sites.begin()));
			
			set<string> longsite = _gall_leo->get_leo_sat_list();
			set<string> short_list;
			for (auto site = longsite.begin(); site != longsite.end(); site++)
			{
				short_list.insert(strl2strs(*site));
			}
			set<string> leo_list;
			set<string>::const_iterator leo_beg = _leo_list.begin();
			set<string>::const_iterator leo_end = _leo_list.end();
			set<string>::const_iterator ics_beg = short_list.begin();
			set<string>::const_iterator ics_end = short_list.end();
			set_intersection(leo_beg, leo_end, ics_beg, ics_end,
				inserter(leo_list, leo_list.begin()));
			_leo_list = leo_list;
			//update _rec_list
			set_intersection(_leo_list.begin(), _leo_list.end(), sites.begin(), sites.end(), inserter(_rec_list, _rec_list.begin()));
			if (_leo_list.size() == 0) {
				write_log_info(_glog, 0, "ERROR: ", "The LEO ics file has no data!");
				return false;
			}
		}

		if (_mode.count("F") > 0)
		{
			set<string> sites; // ground sites
			set_difference(_rec_list.begin(), _rec_list.end(), _leo_list.begin(), _leo_list.end(), inserter(sites, sites.begin()));

			/*set<string> longsite;
			for (auto leo = _leo_list.begin(); leo != _leo_list.end(); leo++)
			{
				longsite.insert(strs2strl(*leo));
			}
			set<string> short_list;
			for (auto site = longsite.begin(); site != longsite.end(); site++)
			{
				short_list.insert(strl2strs(*site));
			}
			set<string>::const_iterator leo_beg = _leo_list.begin();
			set<string>::const_iterator leo_end = _leo_list.end();*/

			//update _rec_list
			set_intersection(_leo_list.begin(), _leo_list.end(), sites.begin(), sites.end(), inserter(_rec_list, _rec_list.begin()));
			if (_leo_list.size() == 0) {
				write_log_info(_glog, 0, "ERROR: ", "The LEO ics file has no data!");
				return false;
			}
		}

		auto tb_slip = dynamic_cast<t_gturboedit*>(_slip12.get());
		if (tb_slip) {
			set<string> rec_list_log = tb_slip->get_sitelist_of_logfile();
			for (auto iter = _rec_list.begin(); iter != _rec_list.end();)
			{
				if (rec_list_log.count(*iter) == 0) {
					_glog->comment(t_glog::LOG_LV::LOG_ERROR, "t_gpodlsqIF", "_initLsqProcData", *iter + " logfile don't exist! Erase the site!");
					iter = _rec_list.erase(iter);
				}
				else {
					iter++;
				}
			}
		}

		write_log_info(_glog, 4, "INIT : process data in lsq is prepared.");
		return true;
	}

	void t_gpodlsqIF::_update_nav_ics()
	{
		// update
		for (auto sat_iter = _sat_list.begin(); sat_iter != _sat_list.end(); sat_iter++)
		{
			auto ics = _gall_gns->get_nav_ics(*sat_iter);
			auto& par_list = ics->get_par_list();

			// cycle par list
			for (auto par_iter = par_list.begin(); par_iter != par_list.end(); par_iter++)
			{
				t_gpar::t_type par_type = t_gics_nav::_icspar2gpar(par_iter->first);

				int par_idx = _lsq->_x_solve.getParam("", par_type, *sat_iter);

				ColumnVector dx = _lsq->dx();
				ColumnVector stdx = _lsq->stdx();

				// update value && sig
				par_iter->second.value += dx(par_idx + 1);
				par_iter->second.sigma = stdx(par_idx + 1);

			}

		}
	}
	
	void t_gpodlsqIF::_update_leo_ics()
	{
		// update
		for (auto sat_iter = _leo_list.begin(); sat_iter != _leo_list.end(); sat_iter++)
		{
			string sat = strs2strl(*sat_iter);
			auto ics = _gall_leo->get_leo_ics(sat);
			map<ICS_PAR, vector<t_icsleo_record>>* par_list = ics->leoics();

			for (auto par_iter = par_list->begin(); par_iter != par_list->end(); par_iter++)
				{
					t_gpar::t_type par_type = t_gics_leo::icspar2gpar(par_iter->first);
					vector<t_icsleo_record> &pars = par_iter->second;
					for (int k = 0; k < pars.size(); k++)
					{
						t_gtime begtime = pars[k].beg;
						t_gtime endtime = pars[k].end;
						int par_idx = _lsq->_x_solve.getParam(*sat_iter, par_type, "", begtime, endtime);
						if (par_idx >= 0)
						{
							ColumnVector dx = _lsq->dx();
							ColumnVector stdx = _lsq->stdx();

							// update value && sig
							pars[k].value += dx(par_idx + 1);
							pars[k].sig = stdx(par_idx + 1);

						}

					}

				}

			}
		}


	void t_gpodlsqIF::_get_epoch_solution()
	{
		t_glsq lsqtmp(*_lsq);
		lsqtmp.solve_NEQ();
		_print_out_content(lsqtmp);

	}
	//	

	void t_gpodlsqIF::_print_out_header()
	{
		ostringstream header_os;

		// output formaat
		header_os << "LSQ:= " << fixed << setprecision(4)
			<< setw(24) << "TIME"
			// write data
			<< setw(15) << "PXSAT[m]"
			<< setw(15) << "PYSAT[m]"
			<< setw(15) << "PZSAT[m]"
			<< setw(15) << "REC_CLK[m]"
			<< setw(15) << "SIGMA0[m]"
			<< setw(15) << "SIGMA_PXSAT[m]"
			<< setw(15) << "SIGMA_PYSAT[m]"
			<< setw(15) << "SIGMA_PZSAT[m]"
			<< setw(15) << "GDOP"
			<< setw(15) << "PDOP"
			<< endl;

		if (_result) {
			_result->write(header_os.str().c_str(), header_os.str().size());
			_result->flush();
		}

	}

	void t_gpodlsqIF::_set_correlation_out()
	{
		//t_gppplsqIF::_setOut();

		string corr_tmp;
		corr_tmp = dynamic_cast<t_gsetout*>(_gset)->outputs("correlation");
		if (corr_tmp.empty()) {
			write_log_info(_glog, 0, "Don't have correlation out file,Please Check out the xml file!!");
			corr_tmp = "correlation";
		}

		substitute(corr_tmp, "$(date)", int2str(_beg_time.year()) + int2str(_beg_time.doy()), false);
		if (_corr_map) {
			if (_corr_map->is_open()) {
				_corr_map->close();
			}
			delete _corr_map;
		}
		_corr_map = new t_giof;
		_corr_map->tsys(t_gtime::GPS);
		_corr_map->mask(corr_tmp);
		_corr_map->append(dynamic_cast<t_gsetout*>(_gset)->append());


	}

	void t_gpodlsqIF::_print_out_content(const t_glsq& lsqestimator)
	{
		ostringstream os;
		t_gallpar par = lsqestimator._x_solve;
		ColumnVector dx = lsqestimator.dx();
		//SymmetricMatrix Qx = lsqestimator.Qx();
		vector<t_gpar> parm = par.getAllPar();
		int i = 0;
		for (auto ipar = parm.begin(); ipar != parm.end(); ipar++)
		{
			i++;
			string site = "NONE";
			string sat = "NONE";
			string isite = ipar->site;
			string isat = ipar->prn;
			if (isite != "")site = isite;
			if (isat != "")sat = isat;
			t_gtime tbeg = ipar->beg;
			string str_dsec_beg = dbl2str(tbeg.dsec());
			t_gtime tend = ipar->end;
			string str_dsec_end = dbl2str(tbeg.dsec());
			os << setw(5) << left << i <<setw(10) << left << ipar->str_type()
				<< setw(10) << right << site << setw(10) << right << sat
				<< setprecision(8) << setw(15) << right << ipar->value()
				<< setprecision(8) << setw(15) << right << dx(i) << setprecision(8) << setw(15) << right << dx(i)
				<< setprecision(8) << setw(15) << right << lsqestimator.stdx(i)
				<< setprecision(4) << setw(25) << right << tbeg.str_ymdhms() << setprecision(4) << str_dsec_beg.substr(2)
				<< setprecision(4) << setw(25) << right << tend.str_ymdhms() << setprecision(4) << str_dsec_end.substr(2)
				<< endl;
			
		}
		if (_result)
		{
			_result->write(os.str().c_str(), os.str().size());
			_result->flush();
		}
		else
		{
			_glog->comment(1, "t_gpodlsq_leo::_printout", "have no output file!");
		}
		
	}

	void t_gpodlsqIF::_readOfiles(int step ,t_gallobs*& obs)
	{
		if (obs) delete obs;
		obs = new t_gallobs();  obs->glog(_glog); obs->gset(_gset);
		t_gdata* gdata = nullptr;

		t_gtime beg_reading = t_gtime::current_time(t_gtime::GPS);
		multimap<IFMT, string> inp = dynamic_cast<t_gsetinp*>(_gset)->inputs_all();
		multimap<IFMT, string>::const_iterator itINP = inp.begin();
		for (size_t i = 0; i < inp.size() && itINP != inp.end(); ++i, ++itINP)
		{
			// Get the file format/path, which will be used in decoder
			IFMT   ifmt(itINP->first);
			if (ifmt != RINEXO_INP) continue;
			string path(itINP->second);
			string id("ID" + int2str(i));

			t_gio* gio = nullptr;
			t_gcoder* gcoder = nullptr;

			// For different file format, we prepare different data container and decoder for them.
			gdata = obs; gcoder = new t_rinexo(_crt_time, _crt_time + step, _gset, "", 5000);

			// Check the file path
			if (path.substr(0, 7) != "file://")
			{
				_glog->comment(0, "main", "path is not file (skipped)!");
				continue;
			}

			// check the obs file valid

			// runepoch for the time costed each epoch (i guess)
			t_gtime runepoch(t_gtime::GPS);
			t_gtime lstepoch(t_gtime::GPS);
			// READ DATA FROM FILE
			if (gcoder) {

				// prepare gio for the file
				gio = new t_gfile;
				gio->glog(_glog);
				gio->path(path);


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
			_gall_obs = obs;
		}
	}

}
