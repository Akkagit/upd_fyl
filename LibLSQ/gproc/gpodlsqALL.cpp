// history : 
//	jdhuang 2019-05-12 add the iono pars init function in _InitLsq

#include "gproc/gpodlsqALL.h"
#include "gset/gsetpar.h"
#include "gset/gsetout.h"
#include "gset/gsetinp.h"
#include "gutils/ginfolog.h"
#include "gutils/gturboedit.h"
#include "gio/gfile.h"
#include "gcoders/icsnav.h"
#include "gdata/gion.h"
#include "gutils/gstring.h"
#include <algorithm>
#include "gcoders/rinexo.h"

namespace great
{
	t_gpodlsqALL::t_gpodlsqALL()
	{
	}

	t_gpodlsqALL::t_gpodlsqALL(t_gsetbase* set, t_gallproc* data, t_glog* log ) :
		t_glsqprocALL(set, data, log)
	{
		// write log
		log->logInfo("", "", "creat t_gpodlsqALL");


		_dd_model = dynamic_cast<t_gsetamb*>(_gset)->dd_mode();


	}

	t_gpodlsqALL::~t_gpodlsqALL()
	{
	}

	//void t_gpodlsqALL::setUseIon(bool useIon)
	//{
	//	_useIon_cmd = useIon;
	//}

	//void t_gpodlsqALL::setAmbfix(bool ambfix)
	//{
	//	_ambfix_cmd = ambfix;
	//}

	bool t_gpodlsqALL::ProcessBatch(t_gallproc * data, const t_gtime & beg, const t_gtime & end)
	{
		// check the glog
		this->_isCreatLog = this->_create_log(this->_glog);

		string class_id = "t_gpodlsqALL";
		string funct_id = "ProcessBatch";

		// check info
		if (data == nullptr)
		{
			this->_glog->logError(class_id, funct_id, "t_gallproc is nullptr");
			return false;
		}
		if (beg > end)
		{		
			this->_glog->logError(class_id, funct_id, "The beg > end, beg is : " + beg.str_mjdsod() + " end is : " + end.str_mjdsod());
			return false;
		}
		if (!_check())
		{
			this->_glog->logError(class_id, funct_id, "can not pass the check function");
			return false;
		}


		// get the time
		if (_beg_time < beg) _beg_time = beg;
		if (_end_time > end) _end_time = end;
		_crt_time = _beg_time;    
		this->_glog->logInfo(class_id, funct_id, "beg time is " + _beg_time.str_mjdsod());
		this->_glog->logInfo(class_id, funct_id, "end time is " + _end_time.str_mjdsod());
		this->_glog->logInfo(class_id, funct_id, "1st time is " + _crt_time.str_mjdsod());
	    
		// init the process data 
		this->_glog->logInfo(class_id, funct_id, "init the process date used in lsq process");
		if (!_initLsqProcData(data))
		{
			this->_glog->logError(class_id, funct_id, "can not init data");
			return false;
		}

		// init the process pars 
		this->_glog->logInfo(class_id, funct_id, "init the process pars ested in lsq process");
		if (!_initLsqProcPars(_lsq))
		{
			this->_glog->logError(class_id, funct_id, "can not init pars");
			return false;
		}

		// loop while time
		this->_glog->logInfo(class_id, funct_id, "start loop process batch !");
		READ_OFILE_MODE rom = dynamic_cast<t_gsetproc*>(_gset)->read_ofile_mode();
		int step = dynamic_cast<t_gsetproc*>(_gset)->read_ofile_step();
		int cmlstep = step;
		t_gallobs* obs = new t_gallobs();  obs->glog(_glog); obs->gset(_gset);
		while (_crt_time <= _end_time)
		{
			if ((rom == READ_OFILE_MODE::REALTIME) && (_crt_time > _beg_time + cmlstep))
			{
				_readOfiles(step, obs);
				cmlstep += step;
			}
			/* get observation */
			vector<t_gsatdata> crt_obs = _gall_obs->obs(_rec_list, _crt_time);
			if (!_processOneEpoch(_crt_time, crt_obs))
			{
				this->_glog->logInfo(class_id, funct_id, _crt_time.str_ymdhms("Processing failed in epoch "));
				_crt_time = _crt_time + _obs_intv;
				continue;
			}

			this->_glog->logInfo(class_id, funct_id, _crt_time.str_ymdhms("Processing successed in epoch "));
			_crt_time = _crt_time + _obs_intv;
		}

		try
		{
			this->_glog->logInfo(class_id, funct_id, "add lsq_sysbias_constraint");
			_lsq->lsq_sysbias_constraint();

			this->_glog->logInfo(class_id, funct_id, "solve_NEQ");
			_lsq->solve_NEQ();
		}
		catch (exception e)
		{
			this->_glog->logInfo(class_id, funct_id, string("Solve Equation Fail!") + e.what());
			return false;
		}

		/* recover equation */
		this->_glog->logInfo(class_id, funct_id, "recover parameters.");
		_lsq->get_result_parameter(*_gall_rcv);
		_lsq->recover_parameters(*_gall_rcv);
		this->_glog->logInfo(class_id, funct_id, "creat recover files.");
		_extract_resfile(*_gall_rcv);

		return true;
	}

	bool t_gpodlsqALL::GenerateProduct()
	{
		bool prod_valid = t_gpodlsqALL::_initLsqProdData(nullptr);
		if (!prod_valid) return false;
		return true;
	}


	bool t_gpodlsqALL::_initLsqProdData(t_gallprod * data)
	{
		// The clk file
		// update clk file
		_extract_clkfile(*_gall_rcv, t_gallprec::AR);
		_extract_clkfile(*_gall_rcv, t_gallprec::AS);

		// 放在外面，用一个extrion的app专门实现这个功能，因为电离层参数不需要更新的这么频繁。
		// t_glsqprocALL::_extract_ion_file(*_allrcv);

		_update_icsfile();

		// Update ICS file
		t_gfile gout; 
		t_ics_nav gics_coder(_gset, "", 4096);
		string ics_path = (dynamic_cast<t_gsetout*>(_gset)->outputs("ics"));
		if (ics_path.empty()) {
			ics_path = "grt$(date).ics";
		}
		if (ics_path.find("$(gwkd)") != string::npos)
		{
			substitute(ics_path, "$(gwkd)", _beg_time.str_gwkd(), false);
		}
		else if (ics_path.find("$(date)") != string::npos)
		{
			substitute(ics_path, "$(date)", _beg_time.str_yyyydoy(), false);
		}

		gout.path(ics_path); 
		gout.glog(_glog);
		gics_coder.clear(); 
		gics_coder.glog(_glog); 
		gics_coder.path(ics_path);
		gics_coder.add_data("ID1", _gall_gns);
		gics_coder.add_data("OBJ", _gall_obj);
		gout.coder(&gics_coder);
		gout.run_write();

		return true;
	}

	bool t_gpodlsqALL::_initLsqProcPars(t_glsq* lsq)
	{
		string class_id = "t_gpodlsqALL";
		string funct_id = "_initLsqProcPars";

		if (_rec_list.empty()){this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "no rec");return false;}
		if (_sat_list.empty()){this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "no sat");return false;}

		// Crd pars
		bool crd_valid  =  _init_rec_crd_pars(lsq, false);
		bool trop_valid = _init_rec_trop_pars(lsq);
		bool clk_valid  = _init_rec_clk_pars(lsq) && _init_sat_clk_pars(lsq);
		//bool ion_valid  = _init_gns_sion_pars(lsq);
		bool isb_valid  = _init_sys_bia_pars(lsq);
		bool ics_valid  = _init_gns_ics_pars(lsq);
		bool erp_valid  = _init_crs_erp_pars(lsq);

		bool ifb_valid = false;
		     if (_frequency  != 2 && _ifb_model == IFB_MODEL::EST_IFB)     ifb_valid = _init_sys_ifb_pars(lsq);
		else if (_frequency  != 2 && _ifb_model == IFB_MODEL::EST_REC_IFB) ifb_valid = _init_sys_rec_ifb_pars(lsq);
		else if (_frequency  != 2 && _ifb_model == IFB_MODEL::EST_SAT_IFB) ifb_valid = _init_sys_sat_ifb_pars(lsq);
		else ifb_valid = true;

		if (!crd_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "crd_valid  is false"); return false; }
		if (!clk_valid ){ this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "clk_valid  is false"); return false;}
		if (!trop_valid){ this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "trop_valid is false"); return false; }
		if (!erp_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "erp_valid  is false"); return false; }
		if (!ics_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "ics_valid  is false"); return false; }
		//if (!ion_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "ion_valid  is false"); return false; }
		if (!isb_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "isb_valid  is false"); return false;}
		if (!ifb_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "ifb_valid  is false"); return false;}


		for (int i = 0; i < lsq->_x_solve.parNumber(); i++)
		{
			cout << setw(15) << lsq->_x_solve[i].str_type()
				<< " prn   " << setw(6) << right << lsq->_x_solve[i].prn
				<< " site  " << setw(6) << right << lsq->_x_solve[i].site
				<< " value " << setw(16) << setprecision(8) << right << lsq->_x_solve[i].value()
				<< " sigma " << setw(16) << setprecision(8) << right << lsq->_x_solve[i].apriori() << endl;
				 
		}
		this->_glog->comment(t_glog::LOG_LV::LOG_INFO, class_id, funct_id, "finish init process pars.");
		return true;
	}

	bool t_gpodlsqALL::_initLsqProcData(t_gallproc * data)
	{
		string class_id = "t_gpodlsqALL";
		string funct_id = "_initLsqProcData";

		// prepare obs data
		_gall_bias = dynamic_cast<t_gallbias*>((*data)[t_gdata::ALLBIAS]);
		_gall_obs = dynamic_cast<t_gallobs*>((*data)[t_gdata::ALLOBS]);
		_gall_obj = dynamic_cast<t_gallobj*>((*data)[t_gdata::ALLOBJ]);
		_gall_rcv = dynamic_cast<t_gallrecover*>((*data)[t_gdata::ALLRECOVER]);
		_gall_gns = dynamic_cast<t_gallics_nav*>((*data)[t_gdata::ALLICS]);
		_gdata_erp = dynamic_cast<t_gpoleut1*>((*data)[t_gdata::ALLPOLEUT1]);
		_gdata_ion = dynamic_cast<t_gion*>((*data)[t_gdata::ION]);
		_gdata_ambcon= dynamic_cast<t_gambcon*>((*data)[t_gdata::AMBCON]);

		if (!_gall_bias){ this->_glog->logError(class_id, funct_id, "can not init lsq process data : _allbias!"); return false; }
		if (!_gall_obs) { this->_glog->logError(class_id, funct_id, "can not init lsq process data : _allobs !"); return false; }
		if (!_gall_obj) { this->_glog->logError(class_id, funct_id, "can not init lsq process data : _allobj !"); return false; }
		if (!_gall_rcv) { this->_glog->logError(class_id, funct_id, "can not init lsq process data : _allrcv !"); return false; }
		if (!_gall_gns) { this->_glog->logError(class_id, funct_id, "can not init lsq process data : _allgns !"); return false; }
		if (!_gdata_erp) { this->_glog->logError(class_id, funct_id, "can not init lsq process data : _allerp !"); return false; }

		// check the ics sat list and xml sat list
		set<string> sat_list;
		set_intersection(_sat_list.begin(), _sat_list.end(), _gall_gns->sat_list.begin(), _gall_gns->sat_list.end(), inserter(sat_list, sat_list.begin()));
		_sat_list = sat_list;
		if (_sat_list.empty())
		{
			this->_glog->logInfo(class_id, funct_id, "sat list is empty");
			return false;
		}

		this->_glog->comment(t_glog::LOG_LV::LOG_INFO, class_id, funct_id, "Finish init process data!");
		return true;
	}


	void t_gpodlsqALL::_update_icsfile()
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
				// jdhuang : 有时候ics和光压模型不匹配，会内存泄漏
				try
				{
					par_iter->second.value += dx(par_idx + 1);
					par_iter->second.sigma = stdx(par_idx + 1);
				}
				catch (...)
				{
					throw logic_error("check your xml and ics file. ");
				}
			}

		}
	}

	bool t_gpodlsqALL::_check()
	{
		return true;
	}

	void t_gpodlsqALL::_readOfiles(int step, t_gallobs*& obs)
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