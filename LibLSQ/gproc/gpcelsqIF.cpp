// history : 
//	jdhuang 2019-05-12 add the iono pars init function in _InitLsq

#include "gproc/gpcelsqIF.h"
#include "gset/gsetpar.h"
#include "gset/gsetout.h"
#include "gset/gsetgen.h"
#include "gutils/ginfolog.h"
#include "gutils/gturboedit.h"
#include "gproc/gupdateparIRC.h"
#include "gproc/gupdateparIF.h"

#include "gproc/gpreproc.h"
#include "gcoders/bnccorr.h"
#include "gcoders/bncobs.h"  //added by liyuhao
#include <gcoders/recover.h>
#include <gio/gfile.h>
#include "gutils/gstring.h"
#include <thread>
#include <chrono>

namespace great
{

	t_gpcelsqIF::t_gpcelsqIF(t_gsetbase* set, t_gallproc* data, t_glog* log ) :
		t_glsqproc(set, data, log)
	{
		// prepare model
		if (_bias_pcemodel)
		{
			delete _bias_pcemodel;
			_bias_pcemodel = nullptr;
		}
		//if (_carrier_range) {
		//	shared_ptr<t_gupdatepar> updateIRC(new t_gupdateparIRC(_slip12));
		//	_lsq->set_update_par(updateIRC);
		//}
		_bias_pcemodel = new t_gpceprecisemodel(data, _glog, _gset);
		//if (_bias_pcemodel)
		//{
		//	delete _bias_pcemodel;
		//	_bias_pcemodel = nullptr;
		//}
		//_bias_pcemodel = new t_gpceprecisemodel(data, _glog, _gset);
		//if (_frequency == 3) {
		//	shared_ptr<t_gupdatepar> update_if13(new t_gupdateparIF13(_slip12, _slip13));
		//	update_if13->set_interval(_obs_intv);
		//	_lsq->set_update_par(update_if13);
		//}


		if (_lsq->mode() == gnut::EPO)
		{
			string satclk_file = "clk_" + _beg_time.str_yyyydoy() + "_epo";
			_satclkfile = new t_giof;
			_satclkfile->tsys(t_gtime::GPS);
			_satclkfile->mask(satclk_file);
			_satclkfile->append(false);
			//if (_frequency == 3)
			//{
			//	satclk_file = "clk13_" + _beg_time.str_yyyydoy() + "_epo";
			//	_satclk13file = new t_giof;
			//	_satclk13file->tsys(t_gtime::GPS);
			//	_satclk13file->mask(satclk_file);
			//	_satclk13file->append(false);
			//}
		}

		//if (_realtime || _lsq->mode() == gnut::EPO)
		//if (_realtime)
		//{
		//	string satclk_file = "clk_" + _beg_time.str_yyyydoy() + "_rtnet";
		//	_rtnetclkfile = new t_giof;
		//	_rtnetclkfile->tsys(t_gtime::GPS);
		//	_rtnetclkfile->mask(satclk_file);
		//	_rtnetclkfile->append(false);
		//}
		_clk_log = new t_glog();
		_clk_log->mask("clk.log");
		_clk_log->cache_size(99);
		_clk_log->tsys(t_gtime::GPS);
		_clk_log->time_stamp(true);

		t_gsetinp* temp_input = dynamic_cast<t_gsetinp*>(set);
		if (temp_input)
		{
			multimap<IFMT, string> inp = temp_input->inputs_all();
			multimap<IFMT, string>::const_iterator itINP = inp.begin();
			for (size_t i = 0; i < inp.size() && itINP != inp.end(); ++i, ++itINP)
			{
				if (itINP->first == BNCSEND_INP)
				{
					string path(itINP->second);
					if (path.substr(0, 6) == "tcp://" && _realtime) {
						_gioout = new t_gtcp(path, true);
						_gioout->glog(log);
						_upload = true;
					}
				}
				if (itINP->first == SP3_INP)
				{
					std::size_t found = (itINP->second).find_last_of("/\\");
					_sp3dir = (itINP->second).substr(0, found);
					_sp3name = (itINP->second).substr(found + 1);
				}
				if (itINP->first != RECOVER_INP) continue;

				this->set_recover_data_by_resfile(itINP->second);
			}
		}
	}

	t_gpcelsqIF::~t_gpcelsqIF()
	{
		// prepare model
		if (_bias_pcemodel)
		{
			delete _bias_pcemodel;
			_bias_pcemodel = nullptr;
		}
		//if (_bias_pcemodel)
		//{
		//	delete _bias_pcemodel;
		//	_bias_pcemodel = nullptr;
		//}
		if (_clk_log)
		{
			delete _clk_log;
			_clk_log = nullptr;
		}
		if (_lsq->mode() == gnut::EPO)
		{
			if (_satclkfile)
			{
				delete _satclkfile;
				_satclkfile = nullptr;
			}
			//if (_frequency == 3 && _satclk13file)
			//{
			//	delete _satclk13file;
			//	_satclk13file = nullptr;
			//}
		}
		//if (_realtime || _lsq->mode() == gnut::EPO)
		//{
		//	if (_rtnetclkfile)
		//	{
		//		delete _rtnetclkfile;
		//		_rtnetclkfile = nullptr;
		//	}
		//}
	}

	bool t_gpcelsqIF::ProcessBatch(t_gallproc * data, const t_gtime & beg, const t_gtime & end)
	{
		int beg_mjd, crt_mjd;
		chrono::high_resolution_clock::time_point beg_epo_time;
		chrono::high_resolution_clock::time_point end_epo_time;
		if (!data ||
			(beg > end))
		{
			write_log_info(_glog, 0, "ERROR", "beg > end or no proc data");
			return false;
		}

		// get the time
		_beg_time = (_beg_time < beg) ? beg : _beg_time;
		_end_time = (_end_time > end) ? end : _end_time;
		_crt_time = beg;
		beg_mjd = _crt_time.mjd();

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

#ifdef USE_OPENMP
		dynamic_cast<t_gprecisebias*>(_bias_model.get())->set_multi_thread(_rec_list);
		omp_set_num_threads(_num_threads);
#endif

		// init gturbo
		if (dynamic_cast<t_gsetproc*>(_gset)->slip_model() == SLIPMODEL::TURBO_EDIT)
		{
			for (auto it = _rec_list.begin(); it != _rec_list.end();) {
				// get site crd at first time
				t_gtriple site_crd;
				if (_get_site_crd(*it, site_crd))
				{
					_gturbo[*it] = make_shared<t_gexeturboedit>(_gset, _glog, data);
					_gturbo[*it]->set_site_crd(*it, site_crd);
					it++;
				}
				else
				{
					it = _rec_list.erase(it);
					continue;
				}

			}
		}

		string sat_ref = dynamic_cast<t_gsetproc*>(_gset)->ref_clk();
		while (_crt_time <= _end_time)
		{
			beg_epo_time = chrono::high_resolution_clock::now();
			//vector<t_gsatdata> crt_obs1 = _allobs->obs(_rec_list, _crt_time);
			//crt_time = _crt_time;
	
			//update sp3 file in realtime
			if (_ultrasp3)
			{
				int time_gap = 1;
				_update_ultrasp3(time_gap);

				//string sp3type, str_yeardoyhour, str_gwkdhour, sp3path;
				//t_gtime sp3time(t_gtime::GPS);
				////判断长名或短名文件
				//if (_sp3name.length() == 15){
				//	sp3type = "short";
				//}
				//else if (_sp3name.size() == 38){
				//	sp3type = "long";
				//}
				////get current sp3 time
				//if (sp3type == "short"){
				//	int gwk = str2int(_sp3name.substr(3, 4));
				//	int wkd = str2int(_sp3name.substr(7, 1));
				//	int sod = str2int(_sp3name.substr(9, 2)) * 3600;
				//	sp3time.from_gwd(gwk, wkd, sod);
				//}
				//else if (sp3type == "long"){
				//	int year = str2int(_sp3name.substr(11, 4));
				//	int doy = str2int(_sp3name.substr(15, 3));
				//	int hour = str2int(_sp3name.substr(18, 2));
				//	int mon = 0;
				//	int day = 0;
				//	sp3time.yeardoy2monthday(year, doy, &mon, &day);
				//	sp3time.from_ymdhms(year, mon, day, hour, 0, 0);
				//}

				//clock_t beg_read_time = clock();
				//t_gio* gio = 0;
				//t_gcoder* gcoder = 0;
				//gcoder = new t_sp3(_gset, "", 8172);

				//int hour_gap;
				//sp3time = sp3time + 3600;
				//while (_crt_time >= sp3time)
				//{
				//	hour_gap = _crt_time.hour() - sp3time.hour();
				//	if (hour_gap > 7)
				//	{
				//		sp3time = sp3time + 7;
				//	}
				//	//update sp3 filename
				//	string name = _sp3name;
				//	if (sp3type == "long"){
				//		str_yeardoyhour = sp3time.str_yyyy() + sp3time.str_doy() + sp3time.str_hms().substr(1, 2);
				//		sp3path = _sp3dir + PATH_SEPARATOR + name.replace(11, 9, str_yeardoyhour);
				//	}
				//	else if (sp3type == "short"){
				//		str_gwkdhour = sp3time.str_gwkd() + "_" + sp3time.str_hms().substr(1, 2);
				//		sp3path = _sp3dir + PATH_SEPARATOR + name.replace(3, 8, str_gwkdhour);
				//	}
				//	substitute(sp3path, GFILE_PREFIX, "");
				//	if (access(sp3path.c_str(), 0) != -1)
				//	{
				//		gio = new t_gfile;
				//		gio->glog(_glog);
				//		gio->path(sp3path);

				//		// Put the file into gcoder
				//		gcoder->clear();
				//		gcoder->path(sp3path);
				//		gcoder->glog(_glog);
				//		// Put the data container into gcoder
				//		string id = "id_sp3";
				//		gcoder->add_data(id, _allnav);
				//		// Put the gcoder into the gio
				//		// Note, gcoder contain the gdata and gio contain the gcoder
				//		gio->coder(gcoder);
				//		// Read the data from file here
				//		gio->run_read();
				//		// get read time 
				//		clock_t end_read_time = clock();
				//		cout << "Time for read" + sp3path + " is " << (end_read_time - beg_read_time) / CLOCKS_PER_SEC << " sec " << endl;

				//		if (sp3type == "long"){
				//			_sp3name.replace(11, 9, str_yeardoyhour);
				//		}
				//		else if (sp3type == "short"){
				//			_sp3name.replace(3, 8, str_gwkdhour);
				//		}
				//		break;
				//	}
				//	sp3time = sp3time + 3600;
				//}
				//if (gio) delete gio; gio = NULL;
				//if (gcoder) delete gcoder; gcoder = NULL;
			}

			// for not realtime but want to use realtime turboedit
			if (!_realtime && _lite_turboedit)
			{
				vector<string> rec_list(_rec_list.begin(),_rec_list.end());
				bool turboedit = (dynamic_cast<t_gsetproc*>(_gset)->slip_model() == SLIPMODEL::TURBO_EDIT);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
				for (int i=0; i < rec_list.size();i++)
				{
					string site = rec_list[i];
					// change _crd_est to EST temporary for SPP
					//t_gsppflt* sppflt = new t_gsppflt(site, _gset);
					//sppflt->setDAT(_allobs, _allnav);
					//sppflt->setOBJ(_allobj);

					//// SPPFlt for kinematic station
					//if (sppflt->processBatch(_crt_time, _crt_time) < 0)
					//{
					//	_glog->comment(1, site, "Kinematic SPPFlt Fail ! Erase!");
					//	_allobs->erase(site, _crt_time);
					//}
					//_outlier_preres[site] = sppflt->get_outlier_sat();
					//if (sppflt) { delete sppflt; sppflt = nullptr; }
					if (turboedit)
					{
						_gturbo.at(site)->ProcessBatch(site, _crt_time, _crt_time, _obs_intv);
					}
					else
					{
						t_gpreproc prep(_gall_obs, _gset);
						prep.glog(_glog);
						t_gtime epo_prev = _crt_time - 300.0;
						prep.setNav(_gall_nav);
						prep.preprocess(site, epo_prev, _crt_time, _obs_intv, false);

						//_allobs->erase(site, _crt_time);
					}
				}
			}

			if (dynamic_cast<t_gsetproc*>(_gset)->realtime())
			{
				while (!((t_bncobs*)_gcoder[0])->available(_crt_time) || (!_ultrasp3 ? !((t_bnccorr*)_gcoder[1])->available(_crt_time) : false))
					t_gtime::gmsleep(1);

				vector<string> rec_list(_rec_list.begin(), _rec_list.end());
				bool turboedit = (dynamic_cast<t_gsetproc*>(_gset)->slip_model() == SLIPMODEL::TURBO_EDIT);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
				for (int i=0;i<rec_list.size();i++)
				{
					string site = rec_list[i];
					// change _crd_est to EST temporary for SPP
					//t_gsppflt* sppflt = new t_gsppflt(site, _gset);
					//sppflt->setDAT(_allobs, _allnav);
					//sppflt->setOBJ(_allobj);

					//// SPPFlt for kinematic station
					//if (sppflt->processBatch(_crt_time, _crt_time) < 0)
					//{
					//	_glog->comment(1, site, "Kinematic SPPFlt Fail ! Erase!");
					//	_allobs->erase(site, _crt_time);
					//}
					//_outlier_preres[site] = sppflt->get_outlier_sat();
					//if (sppflt) { delete sppflt; sppflt = nullptr; }

					if (turboedit)
					{
						//gturbo = new t_gexeturboedit(&gset, &glog, gobs, gorb);
						_gturbo[site]->ProcessBatch(site, _crt_time, _crt_time, _obs_intv);

					}
					else
					{
						t_gpreproc prep(_gall_obs, _gset);
						prep.glog(_glog);
						t_gtime epo_prev = _crt_time - 300.0;
						prep.setNav(_gall_nav);
						prep.preprocess(site, epo_prev, _crt_time, _obs_intv, false);

						//_allobs->erase(site, _crt_time);
					}
				}
			}

			cerr << _crt_time.str_ymdhms() << ":" << endl;
			write_log_info(_glog, 0, "NOTE", _crt_time.str_ymdhms("Processing epoch "));

			/* get observation */
			vector<t_gsatdata> crt_obs = _gall_obs->obs(_rec_list, _crt_time);

			bool epoch_valid;
			if (_realtime || _lsq->mode() == gnut::EPO) {
				_check_out = false;
				epoch_valid = _processOneEpoch(_crt_time, crt_obs);
				if (!epoch_valid)
				{
					write_log_info(_glog, 1, "ERROR" + _crt_time.str_ymdhms("Processing failed in epoch "));
					_crt_time = _crt_time + _obs_intv;
					continue;
				}

				do {
					t_glsq lsqtmp(*_lsq);
					_check_out = outlierDetect(lsqtmp, _map_all_equ);
					if (_outlier_data.size() == 0) continue;
					//epoch_valid = _processOneEpoch(_crt_time, crt_obs);
					outilerDelete();
				} while (!_check_out);

				if (!epoch_valid)
				{
					write_log_info(_glog, 1, "ERROR" + _crt_time.str_ymdhms("Processing failed in epoch "));
					_crt_time = _crt_time + _obs_intv;
					continue;
				}
			}
			else {
				epoch_valid = _processOneEpoch(_crt_time, crt_obs);
				if (!epoch_valid)
				{
					write_log_info(_glog, 1, "ERROR" + _crt_time.str_ymdhms("Processing failed in epoch "));
					_crt_time = _crt_time + _obs_intv;
					continue;
				}
				if (sat_ref == "") _lsq->lsq_clk_constraint();
			}
				
			//bool epoch_valid = _processOneEpoch(_crt_time, crt_obs);
			//if (!epoch_valid)
			//{
			//	write_log_info(_glog, 1, "ERROR" + _crt_time.str_ymdhms("Processing failed in epoch "));
			//	_crt_time = _crt_time + _obs_intv;
			//	continue;
			//}

			if (_realtime || _lsq->mode() == gnut::EPO)
			{
				//临时存储_lsq用于计算dx和Qx
				//t_glsq lsqtmp(*_lsq);
				set<string>::iterator rec_beg = _rec_list.begin();
				set<string>::iterator rec_end = _rec_list.end();
				//for (set<string>::iterator it = rec_beg; it != rec_end; it++)
				//{
				//	// if total obs of the site is zero, skip the site
				//	auto it_site = find(_invalid_site.begin(), _invalid_site.end(), *it);
				//	if (it_site != _invalid_site.end())
				//	{
				//		continue;
				//	}
				//	lsqtmp.add_equation(_IRC, _map_all_equ[*it], _crt_time);
				//}
				//if (!_IRC)
				//{
				//	lsqtmp.lsq_sysbias_constraint();
				//}
				//lsqtmp.solve_NEQ();
				////check post-fit residuals
				//outlierDetect(lsqtmp, _map_all_equ);

				//update observation equations
				for (set<string>::iterator it = rec_beg; it != rec_end; it++)
				{
					// if total obs of the site is zero, skip the site
					auto it_site = find(_invalid_site.begin(), _invalid_site.end(), *it);
					if (it_site != _invalid_site.end())
					{
						continue;
					}
					_lsq->add_equation(_IRC, _map_all_equ[*it], _crt_time);
				}
				_lsq->get_equ_obs_total_num();

				if (_lsq->mode() == gnut::EPO && _crt_time < _end_time)
				{
					crt_mjd = _crt_time.mjd();
					t_glsq lsqepo(*_lsq);
					if (!_IRC)
					{
						lsqepo.lsq_sysbias_constraint();
					}
					if (sat_ref == "") lsqepo.lsq_clk_constraint();
					lsqepo.solve_x();

					if (crt_mjd > beg_mjd)
					{
						//delete last day
						delete _satclkfile;
						_satclkfile = nullptr;
						//creat current day
						string satclk_file = "clk_" + _crt_time.str_yyyydoy() + "_epo";
						_satclkfile = new t_giof;
						_satclkfile->tsys(t_gtime::GPS);
						_satclkfile->mask(satclk_file);
						_satclkfile->append(false);

						beg_mjd = crt_mjd;
					}
					if (_satclkfile) {
						if (_frequency == 3 && _satclk13file)
						{
							_print_clk_solved(&lsqepo, _crt_time, true);
						}
						else
						{
							_print_clk_solved(&lsqepo, _crt_time);
						}
					}

					//output rtnet clk
					if (_realtime && _upload)
					{
						_print_rtnet_clk(&lsqepo, _crt_time);
					}
				}
			}

			write_log_info(_glog, 1, "NOTE", _crt_time.str_ymdhms("End Processing epoch "));
			end_epo_time = chrono::high_resolution_clock::now();
			cout << "Time for " + _crt_time.str_ymdhms("Processing epoch ") + " is " << chrono::duration_cast<chrono::milliseconds>(end_epo_time-beg_epo_time).count() / 1000.0 << " sec " << endl;
			_crt_time = _crt_time + _obs_intv;
		}

		try
		{
			if (sat_ref == "") _lsq->lsq_clk_constraint();
			if (!_IRC)
			{
				_lsq->lsq_sysbias_constraint();
			}
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
			_lsq->get_result_parameter(*_gall_rcv);
			_lsq->recover_parameters(*_gall_rcv);
			_extract_resfile(*_gall_rcv);
		}
		catch (exception e)
		{
			write_log_info(_glog, 1, e.what(), "Recover Parameter Fail!");
			return false;
		}

		return true;
	}

	bool t_gpcelsqIF::GenerateProduct()
	{
		// update clk file
		if (_frequency == 2)
		{
			_extract_clkfile(*_gall_rcv, t_gallprec::AS);
			_extract_clkfile(*_gall_rcv, t_gallprec::AR);
		}
		else if (_frequency == 3)
		{
			_extract_clkfile(*_gall_rcv, t_gallprec::AS, true);
			_extract_clkfile(*_gall_rcv, t_gallprec::AR);
		}
		return true;
	}

	bool t_gpcelsqIF::_initLsqProdData(t_gallprod * data)
	{
		return false;
	}

	bool t_gpcelsqIF::_initLsqProcPars(t_glsq* lsq)
	{
		string class_id = "t_gpcelsqIF";
		string funct_id = "_initLsqProcPars";

		if (_rec_list.empty() || _sat_list.empty())
		{
			this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "no rec or sat");
			return false;
		}

		// Crd pars
		bool rec_crd_valid = _init_rec_crd_pars(lsq);	
		bool sat_clk_valid = _init_sat_clk_pars(lsq);
		bool sat_ifcb_valid = false;
		if (_frequency == 3)
		{
			sat_ifcb_valid = _init_sat_ifcb_pars(lsq);
		}
		bool rec_clk_valid;
		if (_IRC)
		{
			rec_clk_valid = _init_rec_clk_pars_irc(lsq);
		}
		else
		{
			rec_clk_valid = _init_rec_clk_pars(lsq);
		}
		bool rec_clk13_valid = false;
		if (_frequency == 3)
		{
			rec_clk13_valid = _init_rec_clk13_pars(lsq);
		}
		bool trop_valid = _init_rec_trop_pars(lsq);

		bool isb_valid = false;
		if (!_IRC)
		{
			isb_valid = _init_sys_bia_pars_pce(lsq);
		}

		if (!rec_crd_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "rec_crd_valid  is false"); return false; }
		if (!sat_clk_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "sat_clk_valid  is false"); return false; }
		if (!rec_clk_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "rec_clk_valid  is false"); return false; }
		if (!trop_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "trop_valid  is false"); return false; }
		if (_frequency == 3 && !sat_ifcb_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "sat_ifcb_valid  is false"); return false; }
		if (_frequency == 3 && !rec_clk13_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "rec_clk13_valid  is false"); return false; }
		if (!_IRC && !isb_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "isb_valid  is false"); return false; }

		write_log_info(_glog, 1, "INIT", "init the pars estimated in PCE_LSQ_IF.");
		return true;
		//if ((rec_clk_valid && sat_clk_valid && trop_valid) &&
		//	//&& (amb_valid || _IRC) && 
		//	(_frequency == 2 || _frequency == 3 && sat_ifcb_valid && clk13_valid))
		//{
		//	write_log_info(_glog, 1, "INIT","init the pars estimated in PPP_LSQ.");
		//	return true;
		//}
		//else
		//{
		//	write_log_info(_glog, 1, "ERROR","can not init the pars estimated in PPP_LSQ.");
		//	return false;
		//}
	}

	bool t_gpcelsqIF::_update_ultrasp3(int time_gap)
	{
		string sp3type, str_yeardoyhour, str_gwkdhour, sp3path;
		t_gtime sp3time(t_gtime::GPS);
		//判断长名或短名文件
		if (_sp3name.length() == 15) {
			sp3type = "short";
		}
		else if (_sp3name.size() == 38) {
			sp3type = "long";
		}
		//get current sp3 time
		if (sp3type == "short") {
			int gwk = str2int(_sp3name.substr(3, 4));
			int wkd = str2int(_sp3name.substr(7, 1));
			int sod = str2int(_sp3name.substr(9, 2)) * 3600;
			sp3time.from_gwd(gwk, wkd, sod);
		}
		else if (sp3type == "long") {
			int year = str2int(_sp3name.substr(11, 4));
			int doy = str2int(_sp3name.substr(15, 3));
			int hour = str2int(_sp3name.substr(18, 2));
			int mon = 0;
			int day = 0;
			sp3time.yeardoy2monthday(year, doy, &mon, &day);
			sp3time.from_ymdhms(year, mon, day, hour, 0, 0);
		}

		clock_t beg_read_time, end_read_time;
		t_gio* gio = 0;
		t_gcoder* gcoder = 0;
		gcoder = new t_sp3(_gset, "", 8172);

		int hour_gap;
		sp3time = sp3time + 3600 * time_gap;
		while (_crt_time >= sp3time)
		{
			hour_gap = _crt_time.hour() - sp3time.hour();
			if (hour_gap > 7)
			{
				sp3time = sp3time + 7;
			}
			//update sp3 filename
			string name = _sp3name;
			if (sp3type == "long") {
				str_yeardoyhour = sp3time.str_yyyy() + sp3time.str_doy() + sp3time.str_hms().substr(1, 2);
				sp3path = _sp3dir + PATH_SEPARATOR + name.replace(11, 9, str_yeardoyhour);
			}
			else if (sp3type == "short") {
				str_gwkdhour = sp3time.str_gwkd() + "_" + sp3time.str_hms().substr(1, 2);
				sp3path = _sp3dir + PATH_SEPARATOR + name.replace(3, 8, str_gwkdhour);
			}
			substitute(sp3path, GFILE_PREFIX, "");
			if (access(sp3path.c_str(), 0) != -1)
			{
				gio = new t_gfile;
				gio->glog(_glog);
				gio->path(sp3path);
				beg_read_time = clock();

				// Put the file into gcoder
				gcoder->clear();
				gcoder->path(sp3path);
				gcoder->glog(_glog);
				// Put the data container into gcoder
				string id = "id_sp3";
				gcoder->add_data(id, _gall_nav);
				// Put the gcoder into the gio
				// Note, gcoder contain the gdata and gio contain the gcoder
				gio->coder(gcoder);
				// Read the data from file here
				gio->run_read();
				// get read time 
				end_read_time = clock();
				cout << "Time for read" + sp3path + " is " << (end_read_time - beg_read_time) / CLOCKS_PER_SEC << " sec " << endl;

				if (sp3type == "long") {
					_sp3name.replace(11, 9, str_yeardoyhour);
				}
				else if (sp3type == "short") {
					_sp3name.replace(3, 8, str_gwkdhour);
				}
				break;
			}
			sp3time = sp3time + 3600 * time_gap;
		}
		if (gio) delete gio; gio = NULL;
		if (gcoder) delete gcoder; gcoder = NULL;

		return true;
	}

	bool t_gpcelsqIF::_initLsqProcData(t_gallproc * data)
	{
		gtrace("t_gppplsqIF::_initLsqProcData");

		// prepare obs data
		_gall_bias = dynamic_cast<t_gallbias*>((*data)[t_gdata::ALLBIAS]);
		_gall_obs = dynamic_cast<t_gallobs*>((*data)[t_gdata::ALLOBS]);
		_gall_obj = dynamic_cast<t_gallobj*>((*data)[t_gdata::ALLOBJ]);
		_gall_rcv = dynamic_cast<t_gallrecover*>((*data)[t_gdata::ALLRECOVER]);
		_gall_nav = dynamic_cast<t_gallnav*>((*data)[t_gdata::GRP_EPHEM]);
		_gall_nav->overwrite(true);


		if (!_gall_bias ||
			!_gall_obs ||
			!_gall_obj ||
			!_gall_rcv)
		{
			write_log_info(_glog, 0, "ERROR ：can not init lsq process data");
			return false;
		}

		write_log_info(_glog, 4, "INIT : process data in lsq is prepared.");
		return true;
	}

	bool t_gpcelsqIF::_processOneRecOneSat(const t_gtime & crt_epoch,
		const std::string & rec, const std::string & sat, t_gsatdata & crt_obs)
	{
		t_gtime epoch = crt_epoch;
		//bool proc_valid = _bias_model->_cmb_equ(crt_epoch, crt_obs, _lsq->_x_solve, _crt_equ);
		bool proc_valid = _base_model->cmb_equ(epoch, _lsq->_x_solve, crt_obs, _crt_equ);
		string class_id = "t_gpcelsqIF";
		string funct_id = "_processOneRecOneSat";

		if (!proc_valid)
		{
			//write_log_info(_glog, 2, "ERROR", crt_epoch.str_ymdhms("can not process the sat of " +
			//	sat + " and the site of " + rec + " in the time "));
			_glog->comment(t_glog::LOG_LV::LOG_WARN, class_id, funct_id, 
				crt_epoch.str_mjdsod("can not combine equ for code in "  + sat + rec + "at "));
			return false;
		}

		return true;
	}

	bool t_gpcelsqIF::_processOneEpoch(const t_gtime& crt_epoch, std::vector<t_gsatdata>& crt_obs)
	{
		string class_id = "t_gpcelsqIF";
		string funct_id = "_processOneEpoch";
		_invalid_site.clear();
		_map_all_equ.clear();

		// updata pars
		if (!_slip12 || crt_obs.empty()) return false;
		//// for debug
		ColumnVector W = _lsq->W();					///< BTPL Matrix
		SymmetricMatrix NEQ = _lsq->NEQ();

		bool updata_valid;
		if (!_check_out)
		{
			updata_valid = _lsq->update_parameter(_IRC, crt_epoch, crt_obs, false);
			if (!updata_valid ||
				crt_obs.empty())
			{
				write_log_info(_glog, 1, crt_epoch.str(), ": no observation found ");
				return false;
			}
		}

		// combine obs equations
		//t_glsqEquationMatrix equ;

		int iObs = 0;
		vector<t_gsatdata>::iterator it_obs;
		vector<string>::iterator it_outlier1;
		vector<pair<string, string> >::iterator it_outlier2;
		set<string>::iterator rec_beg = _rec_list.begin();
		set<string>::iterator rec_end = _rec_list.end();
		vector<t_gsatdata>::iterator obs_beg = crt_obs.begin();
		vector<t_gsatdata>::iterator obs_end = crt_obs.end();

		// prepare Obs
		map < string, vector<t_gsatdata> >map_site_obs;
		for (set<string>::iterator it = rec_beg; it != rec_end; it++)
		{
			vector<t_gsatdata> crt_rec_obs;
			string crt_rec = *it;

			// get all the obs for the rec
			bool select_valid = _select_rec_obs(crt_rec, crt_epoch, crt_obs, crt_rec_obs);
			if (!select_valid)
			{
				_glog->logInfo("t_gpcelsqIF", "_processOneEpoch", crt_epoch.str_mjdsod("no useful data : " + crt_rec));
				_invalid_site.push_back(crt_rec);
				continue;
			}
			it_obs = crt_rec_obs.begin();
			while (it_obs != crt_rec_obs.end()) {
				string sat = it_obs->sat();
				//spp res outlier
				it_outlier1 = find(_outlier_preres[crt_rec].begin(), _outlier_preres[crt_rec].end(), sat);
				if (it_outlier1 != _outlier_preres[crt_rec].end()) {
					it_obs = crt_rec_obs.erase(it_obs);
					continue;
				}
				//post-fit res outlier
				it_outlier2 = find(_outlier_data.begin(), _outlier_data.end(), make_pair(crt_rec, sat));
				if (it_outlier2 != _outlier_data.end()) {
					it_obs = crt_rec_obs.erase(it_obs);
					continue;
				}
				it_obs++;
			}
			//number of sat less than 4 continue
			if (crt_rec_obs.size() < 4) {
				_invalid_site.push_back(crt_rec);
				continue;
			}
			map_site_obs[crt_rec] = crt_rec_obs;
		}

		
		vector<string> vec_sites;
		for (auto item : map_site_obs) {
			vec_sites.push_back(item.first);
		}
		t_gmutex add_mtx;
		// multi-thread  Process REC [Use OpenMP]
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
		for (int site_i = 0; site_i < vec_sites.size(); site_i++)
		{
			t_glsqEquationMatrix equ_temp;
			string rec_temp = vec_sites[site_i];
			// Process this site and get the equations
			//_glog->logInfo("t_gpcelsqIF", "_processOneEpoch", "Begin process " + rec_temp);
			bool proc_site = _processOneRec_thread_safe(crt_epoch, rec_temp,map_site_obs[rec_temp],equ_temp);
			if (!proc_site)
			{
				cout << crt_epoch.str_ymdhms(rec_temp + " has no equations ", false, false) << endl;
				_glog->logInfo("t_gpcelsqIF", "_processOneEpoch", crt_epoch.str_mjdsod("no useful equations : " + rec_temp));
				continue;
			}

			// add the new equations

			add_mtx.lock();
			//TODO:EPO模式下先做残差检验再添加法方程
			if (_lsq->mode() == gnut::EPO)
			{
				_map_all_equ[rec_temp] = equ_temp;
			}
			else
			{
				_lsq->add_equation(_IRC, equ_temp, crt_epoch);
			}
			add_mtx.unlock();
			//_lsq->add_equation(_IRC, _crt_equ, crt_epoch);
			


#ifdef DEBUG_BPL_LOG
			//Matrix tmp = B;
			//for (int i = 1; i <= tmp.Nrows(); i++)
			//{
			//	char row[4]; sprintf(row, "%4d", i); string rows(row);
			//	string tmp_str;
			//	for (int j = 1; j <= tmp.Ncols(); j++)
			//	{
			//		tmp_str.clear();
			//		char col[4];        sprintf(col, "%4d", j); string cols(col);
			//		char tmp_array[20]; sprintf(tmp_array, "%18.8f", tmp(i, j));

			//		tmp_str = tmp_str + rows + cols + tmp_array;
			//		write_log_info(_B_log, 0, "DEBUG", tmp_str);
			//	}
			//}

			//DiagonalMatrix tmp2 = P;
			//for (int i = 1; i <= tmp2.Nrows(); i++)
			//{
			//	char row[4]; sprintf(row, "%4d", i); string rows(row);
			//	string tmp_str;

			//	char col[4];        sprintf(col, "%4d", 1); string cols(col);
			//	char tmp_array[20]; sprintf(tmp_array, "%18.8f", tmp2(i));
			//	tmp_str = tmp_str + rows + cols + tmp_array;
			//	write_log_info(_P_log, 0, "DEBUG", tmp_str);
			//}

			/*ColumnVector tmp3 = l;
			for (int i = 1; i <= tmp3.Nrows(); i++)
			{
				char row[4]; sprintf(row, "%4d", i); string rows(row);
				string tmp_str;
				char col[4];        sprintf(col, "%4d", 1); string cols(col);
				char tmp_array[20]; sprintf(tmp_array, "%18.8f", tmp3(i));
				tmp_str = tmp_str + rows + cols + tmp_array;
				write_log_info(_L_log, 0, "DEBUG", tmp_str);
			}*/
#endif // DEBUG_BPL_LOG

		}

		return true;
	}

	bool t_gpcelsqIF::_processOneRec(const t_gtime& crt_epoch, const std::string& crt_rec, std::vector<t_gsatdata>& crt_obs)
	{
		string class_id = "t_gpcelsqIF";
		string funct_id = "_processOneRec";

		if (crt_obs.empty() ||
			crt_rec.empty())
		{
			if (_glog) _glog->comment(2, "t_ppplsq::_processSite", "This epoch have no useful data: " + crt_rec);
			return false;
		}
		_crt_rec = crt_rec;

		vector<double> codeOmc;

		// loop for Combine equation
		bool while_valid = false;
		do
		{
			// process one rec one sat and get equ
			_crt_equ.clear_allequ();

			std::vector<t_gsatdata>::iterator it_beg = crt_obs.begin();
			std::vector<t_gsatdata>::iterator it_end = crt_obs.end();
			string crt_sat("NONE");
			for (auto it = it_beg; it != it_end; it++)
			{
				if (_IRC && _carrier_range)
				{
					bool exist = false;
					map<GOBS, int> lli = it->lli();
					for (auto it_lli = lli.begin(); it_lli != lli.end(); it_lli++)
					{
						if (it_lli->second == 9)
						{
							exist = true;
						}
					}
					if (!exist)
					{
						continue;
					}
				}
				if (_sat_list.empty() ||
					_sat_list.find(it->sat()) == _sat_list.end())
				{
					continue;
				}
				crt_sat = it->sat();
				bool proc_valid = _processOneRecOneSat(crt_epoch, crt_rec, crt_sat, *it);
				if (!proc_valid)
				{
					//cout << crt_epoch.str_ymdhms(crt_epoch.str_ymdhms("can not process the rec : " + crt_rec + " and the sat : " + crt_sat + " "), false, false) << endl;
					_glog->comment(t_glog::LOG_LV::LOG_WARN, class_id, funct_id, crt_epoch.str_mjdsod("can not process the rec : " + crt_rec + " and the sat : " + crt_sat + " "));
					continue;
				}

			}

			codeOmc = _crt_equ.get_codeomc();	

			//reset recclk omc, check range
			double recclk;
			double sigma;
			if (_IRC)
			{
				int k = Check_Range(codeOmc, recclk, sigma);

				bool allvalid = true;
				vector<int> idx_irc;
				set<string>::const_iterator itsys = _sys_list.begin();
				t_gpar::t_type clktype;
				for (itsys; itsys != _sys_list.end(); itsys++)
				{
					GSYS crt_sys = t_gsys::str2gsys(*itsys);
					switch (crt_sys)
					{
						case GSYS::GPS:
							clktype = t_gpar::CLK_G;
							break;
						case GSYS::BDS:
							clktype = t_gpar::CLK_C;
							break;
						case GSYS::GAL:
							clktype = t_gpar::CLK_E;
							break;
						case GSYS::GLO:
							clktype = t_gpar::CLK_R;
							break;
						case GSYS::QZS:
							clktype = t_gpar::CLK_J;
							break;
						default:
							break;
					}
					idx_irc.push_back(_lsq->_x_solve.getParam(crt_rec, clktype, ""));
				}

				// reset recclk
				if (idx_irc[0] >= 0 && !while_valid && (fabs(recclk / CLIGHT) > 1e-6 ||
					(k >= 1 && fabs(_lsq->_x_solve[idx_irc[0]].value() / CLIGHT) < 1e-15)))
				//if (idx_irc[0] >= 0 && !while_valid && (fabs(recclk / CLIGHT) > 1e-6 ||
				//	(k >= 1 && fabs(_crt_time == _beg_time ? _lsq->_x_solve[idx_irc[0]].value() / CLIGHT 
				//		: _bias_pcemodel->get_rec_clk(crt_rec)) < 1e-15)))
				{
					double recclk0;
					for (int isys = 0; isys < idx_irc.size(); isys++)
					{
						recclk0 = recclk;
						if (_crt_time == _beg_time)
						{
							recclk0 += _lsq->_x_solve[idx_irc[isys]].value();
						}
						else
						{
							recclk0 += _bias_model->get_rec_clk(crt_rec) * CLIGHT;
						}
						//recclk0 += _lsq->_x_solve[idx_irc[isys]].value();
						_lsq->_x_solve[idx_irc[isys]].value(recclk0);
					}
					_bias_model->update_obj_clk(crt_rec, crt_epoch, recclk0 / CLIGHT);

					codeOmc.clear();
					while_valid = true;
				}
				else
				{
					while_valid = false;
				}
			}
			else
			{
				int k = Check_Range(codeOmc, recclk, sigma);
				int idx = _lsq->_x_solve.getParam(crt_rec, t_gpar::CLK, "");
				if (idx >= 0 && !while_valid && (fabs(recclk / CLIGHT) > 1e-6 ||
					(k >= 1 && fabs(_lsq->_x_solve[idx].value() / CLIGHT) < 1e-15)))
				{
					// reset recclk
					recclk += _lsq->_x_solve[idx].value();
					_lsq->_x_solve[idx].value(recclk);
					_bias_model->update_obj_clk(crt_rec, crt_epoch, recclk / CLIGHT);

					codeOmc.clear();
					while_valid = true;
				}
				else
				{
					while_valid = false;
				}
			}
		} while (while_valid);

		if (!_init_sat_clk())
		{
			_glog->logInfo("t_gpcelsqIF", "_processOneRec", crt_epoch.str_mjdsod("_init_sat_clk failed" + _crt_sat));
			return false;
		}

		return true;
	}

	bool t_gpcelsqIF::_processOneRec_thread_safe(const t_gtime & crt_epoch, const std::string & crt_rec,const std::vector<t_gsatdata>& crt_obs, t_glsqEquationMatrix & equ_result)
	{
		string class_id = "t_gpcelsqIF";
		string funct_id = "_processOneRec";

		if (crt_obs.empty() ||
			crt_rec.empty())
		{
			if (_glog) _glog->comment(2, "t_ppplsq::_processSite", "This epoch have no useful data: " + crt_rec);
			return false;
		}

		vector<double> codeOmc;

		// loop for Combine equation
		bool while_valid = false;
		do
		{
			// process one rec one sat and get equ
			equ_result.clear_allequ();

			auto it_beg = crt_obs.begin();
			auto it_end = crt_obs.end();
			string crt_sat("NONE");
			for (auto it:crt_obs)
			{
				if (_IRC && _carrier_range)
				{
					bool exist = false;
					map<GOBS, int> lli = it.lli();
					for (auto it_lli = lli.begin(); it_lli != lli.end(); it_lli++)
					{
						if (it_lli->second == 9)
						{
							exist = true;
						}
					}
					if (!exist)
					{
						continue;
					}
				}
				if (_sat_list.empty() ||
					_sat_list.find(it.sat()) == _sat_list.end())
				{
					continue;
				}
				crt_sat = it.sat();
				t_gtime epo(crt_epoch);
				if (!_base_model->cmb_equ(epo, _lsq->_x_solve, it, equ_result))
				{
					//cout << crt_epoch.str_ymdhms(crt_epoch.str_ymdhms("can not process the rec : " + crt_rec + " and the sat : " + crt_sat + " "), false, false) << endl;
					_glog->comment(t_glog::LOG_LV::LOG_WARN, class_id, funct_id, crt_epoch.str_mjdsod("can not process the rec : " + crt_rec + " and the sat : " + crt_sat + " "));
					continue;
				}

			}

			codeOmc = equ_result.get_codeomc();

			//reset recclk omc, check range
			double recclk;
			double sigma;
			if (_IRC)
			{
				int k = Check_Range(codeOmc, recclk, sigma);

				bool allvalid = true;
				vector<int> idx_irc;
				set<string>::const_iterator itsys = _sys_list.begin();
				t_gpar::t_type clktype;
				for (itsys; itsys != _sys_list.end(); itsys++)
				{
					GSYS crt_sys = t_gsys::str2gsys(*itsys);
					switch (crt_sys)
					{
					case GSYS::GPS:
						clktype = t_gpar::CLK_G;
						break;
					case GSYS::BDS:
						clktype = t_gpar::CLK_C;
						break;
					case GSYS::GAL:
						clktype = t_gpar::CLK_E;
						break;
					case GSYS::GLO:
						clktype = t_gpar::CLK_R;
						break;
					case GSYS::QZS:
						clktype = t_gpar::CLK_J;
						break;
					default:
						break;
					}
					idx_irc.push_back(_lsq->_x_solve.getParam(crt_rec, clktype, ""));
				}

				// reset recclk
				if (idx_irc[0] >= 0 && !while_valid && (fabs(recclk / CLIGHT) > 1e-6 ||
					(k >= 1 && fabs(_lsq->_x_solve[idx_irc[0]].value() / CLIGHT) < 1e-15)))
					//if (idx_irc[0] >= 0 && !while_valid && (fabs(recclk / CLIGHT) > 1e-6 ||
					//	(k >= 1 && fabs(_crt_time == _beg_time ? _lsq->_x_solve[idx_irc[0]].value() / CLIGHT 
					//		: _bias_pcemodel->get_rec_clk(crt_rec)) < 1e-15)))
				{
					double recclk0;
					for (int isys = 0; isys < idx_irc.size(); isys++)
					{
						recclk0 = recclk;
						if (crt_epoch == _beg_time)
						{
							recclk0 += _lsq->_x_solve[idx_irc[isys]].value();
						}
						else
						{
							recclk0 += _bias_model->get_rec_clk(crt_rec) * CLIGHT;
						}
						//recclk0 += _lsq->_x_solve[idx_irc[isys]].value();
						_lsq->_x_solve[idx_irc[isys]].value(recclk0);
					}
					_bias_model->update_obj_clk(crt_rec, crt_epoch, recclk0 / CLIGHT);

					codeOmc.clear();
					while_valid = true;
				}
				else
				{
					while_valid = false;
				}
			}
			else
			{
				int k = Check_Range(codeOmc, recclk, sigma);
				int idx = _lsq->_x_solve.getParam(crt_rec, t_gpar::CLK, "");
				if (idx >= 0 && !while_valid && (fabs(recclk / CLIGHT) > 1e-6 ||
					(k >= 1 && fabs(_lsq->_x_solve[idx].value() / CLIGHT) < 1e-15)))
				{
					// reset recclk
					recclk += _lsq->_x_solve[idx].value();
					_lsq->_x_solve[idx].value(recclk);
					_bias_model->update_obj_clk(crt_rec, crt_epoch, recclk / CLIGHT);

					codeOmc.clear();
					while_valid = true;
				}
				else
				{
					while_valid = false;
				}
			}
		} while (while_valid);

		//TODO:How to init_sat_clk when multi-thread????????
		//if (!_init_sat_clk())
		//{
		//	_glog->logInfo("t_gpcelsqIF", "_processOneRec", crt_epoch.str_mjdsod("_init_sat_clk failed" + _crt_sat));
		//	return false;
		//}

		return true;
	}


	bool t_gpcelsqIF::outlierDetect(t_glsq& lsq, map<string, t_glsqEquationMatrix>& all_equ)
	{
		Matrix B_value;
		DiagonalMatrix P_value;
		ColumnVector l_value;
		map<string, map<string, double>> res;
		//vector<pair<string, string> > erase_data;
		vector<pair<string, string> >::iterator it_erase;
		//idx,site,sat,obs_type
		map<int, pair<pair<string, string>, string>> obs_sitesat;
		map<int, pair<pair<string, string>, string>>::iterator it_sitesat;
		vector<string>::iterator it_site;
		pair<string, string> site_sat;

		set<string>::iterator rec_beg = _rec_list.begin();
		set<string>::iterator rec_end = _rec_list.end();

		int i, num, ipar;
		int num_equ = 0;
		string site, sat, type;
		GSYS gs;

		//solve equation
		for (set<string>::iterator it = rec_beg; it != rec_end; it++)
		{
			// if total obs of the site is zero, skip the site
			auto it_site = find(_invalid_site.begin(), _invalid_site.end(), *it);
			if (it_site != _invalid_site.end())
			{
				continue;
			}
			lsq.add_equation(_IRC, _map_all_equ[*it], _crt_time);
		}
		if (!_IRC)
		{
			lsq.lsq_sysbias_constraint();
		}
		lsq.solve_x();

		//ColumnVector v_orig, v_norm;
		ColumnVector dx = lsq.dx();

		//clear old outlier data
		_outlier_data.clear();

		//get totoal obs and par number of current epoch
		int nobs = lsq.get_equ_obs_total_num();
		int npar = lsq._x_solve.parNumber();

		B_value.ReSize(nobs, npar); B_value = 0.0;
		l_value.ReSize(nobs);       l_value = 0.0;
		P_value.ReSize(nobs);       P_value = 0.0;
		//check post-fit residuals
		for (set<string>::iterator it = rec_beg; it != rec_end; it++)
		{
			site = *it;
			res.clear();
			// if total obs of the site is zero, skip the site
			it_site = find(_invalid_site.begin(), _invalid_site.end(), site);
			if (it_site != _invalid_site.end())
			{
				continue;
			}

			for (num = 0; num < all_equ[site].num_equ(); num++)
			{
				sat = all_equ[site].get_satname(num);
				type = all_equ[site].get_obscombtype2str(num);
				gs = t_gsys::str2gsys(sat.substr(0, 1));
				site_sat = make_pair(site, sat);
				for (ipar = 0; ipar < all_equ[site].B[num].size(); ipar++)
				{
					//B matrix
					//B_value(num + 1 + num_equ, all_equ[site].B[num][ipar].first) = all_equ[site].B[num][ipar].second;
					//B * dx
					res[sat][type] += all_equ[site].B[num][ipar].second * dx(all_equ[site].B[num][ipar].first);
					//cout << all_equ[*it].B[num][ipar].second << " "  << dx(all_equ[*it].B[num][ipar].first) << " " << res[sat][type] << endl;
				}
				//post-fit residuals
				//cout << all_equ[*it].l[num] << endl;
				res[sat][type] -= all_equ[site].l[num];
				//P matrix
				//P_value(num + 1 + num_equ) = all_equ[site].P[num];
				//l matrix
				//l_value(num + 1 + num_equ) = all_equ[site].l[num];
				//save position of erase site_sat
				//obs_sitesat[num + 1 + num_equ].first = site_sat;
				//obs_sitesat[num + 1 + num_equ].second = type;

				//check outlier post-fit residual and save bad site_sat
				if (type == "LC")
				{
					if (fabs(res[sat][type]) > _maxres_L[gs])
					{
						it_erase = find(_outlier_data.begin(), _outlier_data.end(), site_sat);
						if (it_erase == _outlier_data.end())
						{
							_outlier_data.push_back(site_sat);
							_glog->logInfo("t_gpcelsqIF", "outlierDetect:", " outlier (v_orig: " + type + " " + site + " " + sat + ":" + dbl2str(fabs(res[sat][type])) + ")");
							//cout << site_sat.first << " " << site_sat.second << endl;
						}
					}
				}
				else if (type == "PC")
				{
					if (fabs(res[sat][type]) > _maxres_C[gs])
					{
						it_erase = find(_outlier_data.begin(), _outlier_data.end(), site_sat);
						if (it_erase == _outlier_data.end())
						{
							_outlier_data.push_back(site_sat);
							_glog->logInfo("t_gpcelsqIF", "outlierDetect:", " outlier (v_orig: " + type + " " + site + " " + sat + ":" + dbl2str(fabs(res[sat][type])) + ")");
							//cout << site_sat.first << " " << site_sat.second << endl;
						}
					}
				}
			}
			num_equ += all_equ[site].num_equ();
		}

		// check normalized post-fit residuals
		//v_orig = l_value - B_value * dx;
		//Matrix Qv = P_value.i() - B_value * Qx * B_value.t();
		//v_norm.ReSize(v_orig.Nrows());
		//for (int i = 1; i <= v_norm.Nrows(); i++) {
		//	v_norm(i) = sqrt(1 / Qv(i, i)) * v_orig(i);
		//	if (fabs(v_norm(i)) > _maxres_norm)
		//	{
		//		it_erase = find(erase_data.begin(), erase_data.end(), obs_sitesat[i].first);
		//		if (it_erase == erase_data.end())
		//		{
		//			erase_data.push_back(obs_sitesat[i].first);
		//			_glog->logInfo("t_gpcelsqIF", "outlierDetect:", " outlier (v_norm: " + obs_sitesat[i].second + " " + obs_sitesat[i].first.first + " " + obs_sitesat[i].first.second + ")");
		//		}
		//	}
		//	//if (obs_sitesat[i].second == "LC")
		//	//{
		//	//	if (fabs(v_orig(i)) > _maxres_L[gs])
		//	//	{
		//	//		it_erase = find(erase_data.begin(), erase_data.end(), obs_sitesat[i].first);
		//	//		if (it_erase == erase_data.end())
		//	//		{
		//	//			erase_data.push_back(obs_sitesat[i].first);
		//	//			_glog->logInfo("t_gpcelsqIF", "outlierDetect:", " outlier ( v_norm: " + obs_sitesat[i].second + " " + obs_sitesat[i].first.first + " " + obs_sitesat[i].first.second);
		//	//		}
		//	//	}
		//	//}
		//	//else if (obs_sitesat[i].second == "PC")
		//	//{
		//	//	if (fabs(v_orig(i)) > _maxres_C[gs])
		//	//	{
		//	//		it_erase = find(erase_data.begin(), erase_data.end(), obs_sitesat[i].first);
		//	//		if (it_erase == erase_data.end())
		//	//		{
		//	//			erase_data.push_back(obs_sitesat[i].first);
		//	//			_glog->logInfo("t_gpcelsqIF", "outlierDetect:", " outlier ( v_norm: " + obs_sitesat[i].second + " " + obs_sitesat[i].first.first + " " + obs_sitesat[i].first.second);
		//	//		}
		//	//	}
		//	//}
		//}

		//num_equ = 0;
		//erase bad site_sat in obs equ
		//for (set<string>::iterator it = rec_beg; it != rec_end; it++)
		//{
		//	site = *it;
		//	it_site = find(_invalid_site.begin(), _invalid_site.end(), site);
		//	if (it_site != _invalid_site.end())
		//	{
		//		continue;
		//	}
		//	num_equ = all_equ[site].num_equ();
		//	for (num = 0; num < num_equ; num++)
		//	{
		//		sat = all_equ[site].get_satname(num);
		//		type = all_equ[site].get_obscombtype2str(num);
		//		pair<string, string> site_sat(site, sat);
		//		it_erase = find(erase_data.begin(), erase_data.end(), site_sat);
		//		if (it_erase != erase_data.end())
		//		{
		//			all_equ[site].remove(num + 1);
		//			num -= 1;
		//			num_equ -= 1;
		//		}
		//	}

		//	//if number of sat less than 4, erase the site
		//	if (all_equ[site].num_equ() < 8)
		//	{
		//		_invalid_site.push_back(site);
		//		_glog->logInfo("t_gpcelsqIF", "outlierDetect:", site + " number of sat obs less than 4");
		//	}
		//}

		return true;
	}

	void t_gpcelsqIF::outilerDelete()
	{
		for (auto iter = _map_all_equ.begin(); iter != _map_all_equ.end();iter++)
		{
			string rec = iter->first;
			for (int i = 0; i < iter->second.num_equ();) {
				string sat = iter->second.get_satname(i);

				//spp res outlier
				auto it_outlier1 = find(_outlier_preres[rec].begin(), _outlier_preres[rec].end(), sat);
				if (it_outlier1 != _outlier_preres[rec].end()) {
					iter->second.remove(i + 1);
					continue;
				}

				//post-fit res outlier
				auto it_outlier2 = find(_outlier_data.begin(), _outlier_data.end(), make_pair(rec, sat));
				if (it_outlier2 != _outlier_data.end()) {
					iter->second.remove(i + 1);
					continue;
				}
				i++;
			}
		}
	}

	void t_gpcelsqIF::_print_clk_solved(t_glsq* lsq, t_gtime& crt_epoch, bool isTriple)
	{
		//set<string>::const_iterator rec_end = _rec_list.end();
		set<string>::const_iterator sat_end = _sat_list.end();
		//set<string>::const_iterator it_rec = _rec_list.begin();
		set<string>::const_iterator it_sat = _sat_list.begin();

		t_gallpar par = lsq->_x_solve;
		ColumnVector dx = lsq->dx();

		ostringstream os;
		t_gtime creat_t(t_gtime::UTC);
		creat_t.add_secs(28800);

		string sat_ref = dynamic_cast<t_gsetproc*>(_gset)->ref_clk();
		// write sp3 clk file header
		if (!_satclkfile->is_open()) {
			os << "     2.00           C                                       RINEX VERSION / TYPE" << endl;

			os << left << setw(20) << "G-Nut" << setw(20) << "GREAT" << setw(20) << creat_t.str("%Y%m%d %H%M%S") << "PGM / RUN BY / DATE" << endl;

			os << right << setw(9) << fixed << setprecision(2) << _obs_intv << setw(51) << " " << "INTERVAL" << endl;

			// write # / TYPES OF DATA
			os << "    " << right << setw(2) << "1" << "    " << "AS" << setw(48) << "" << "# / TYPES OF DATA" << endl;

			// write # OF CLK REF
			os << "    " << right << setw(2) << "1" << setw(54) << "" << "# OF CLK REF" << endl;

			// write ANALYSIS CLK REF
			os << setw(4) << sat_ref << setw(56) << " " << "ANALYSIS CLK REF" << endl;

			// write END OF HEADER
			os << setw(60) << " " << "END OF HEADER" << endl;
		}

		// Print sp3 clk file header
		if (_satclkfile) {
			_satclkfile->write(os.str().c_str(), os.str().size());
			_satclkfile->flush();
			//if (isTriple)
			//{
			//	_satclk13file->write(os.str().c_str(), os.str().size());
			//	_satclk13file->flush();
			//}
			os.str("");
		}
		else {
			_glog->comment(2, "t_gpcelsqIF::_write_clk_header", "have no output file!");
		}

		// write sp3 clk file data
		double satclk;
		for (it_sat; it_sat != sat_end; it_sat++) {
			int idx = par.getParam("", t_gpar::CLK_SAT, *it_sat, crt_epoch, crt_epoch);
			if (idx < 0) continue;
			int idx_ifcb;

			double x = par[idx].value();
			double dclk = dx(idx + 1);
			double x_ifcb;
			double difcb;
			if (isTriple)
			{
				idx_ifcb = par.getParam("", t_gpar::IFCB_F3, *it_sat, crt_epoch, crt_epoch);
				if (idx_ifcb > -1)
				{
					x_ifcb = par[idx_ifcb].value();
					difcb = dx(idx_ifcb + 1);
				}
			}
			double clk;
			double ifcb = 0.0;
			if (dclk != 0) {
				//continue;
				clk = x + dclk;
				satclk = clk / CLIGHT;
				if (idx_ifcb > -1)
				{
					ifcb = (x_ifcb + difcb) / CLIGHT;
				}

				if (!isTriple || ifcb == 0.0)
				{
					os << "AS" << " " << setw(3) << *it_sat << "  " <<
						setw(26) << crt_epoch.str("%Y %C %L %K %O %P") << " " << setw(2)
						<< "1" << "   " << setw(19) << fixed << scientific << setprecision(12) << satclk << endl;
				}
				else
				{
					os << "AS" << " " << setw(3) << *it_sat << "  " <<
						setw(26) << crt_epoch.str("%Y %C %L %K %O %P") << " " << setw(2)
						<< "4" << "   " << setw(19) << fixed << scientific << setprecision(12) << satclk
						<< " " << setw(19) << fixed << scientific << setprecision(12) << 0.0
						<< " " << setw(19) << fixed << scientific << setprecision(12) << 0.0
						<< " " << setw(19) << fixed << scientific << setprecision(12) << ifcb << endl;
				}
			}

			//if (isTriple)
			//{
			//	t_gpar::t_type clktype;
			//	GSYS crt_sys = t_gsys::str2gsys((*it_sat).substr(0, 1));
			//	switch (crt_sys)
			//	{
			//	case GSYS::GPS:
			//		clktype = t_gpar::G15CLK_SAT;
			//		break;
			//	case GSYS::BDS:
			//		clktype = t_gpar::C15CLK_SAT;
			//		break;
			//	case GSYS::GAL:
			//		clktype = t_gpar::E15CLK_SAT;
			//		break;
			//	case GSYS::QZS:
			//		clktype = t_gpar::J15CLK_SAT;
			//		break;
			//	default:
			//		break;
			//	}
			//	idx = par.getParam("", clktype, *it_sat, crt_epoch, crt_epoch);

			//	x = par[idx].value();
			//	dclk = dx(idx + 1);
			//	if (dclk != 0) {
			//		clk = x + dclk;
			//		satclk = clk / CLIGHT;

			//		os13 << "AS" << " " << setw(3) << *it_sat << "  " <<
			//			setw(26) << crt_epoch.str("%Y %C %L %K %O %P") << " " << setw(2)
			//			<< "1" << "   " << setw(19) << fixed << scientific << setprecision(12) << satclk << endl;
			//	}
			//}

			//double clk = x + dclk;
			//satclk = clk / CLIGHT;

			//os << "AS" << " " << setw(3) << *it_sat << "  " <<
			//	setw(26) << crt_epoch.str("%Y %C %L %K %O %P") << " " << setw(2)
			//	<< "1" << "   " << setw(19) << fixed << scientific << setprecision(12) << satclk << endl;

		}

		// Print sp3 clk file data
		if (_satclkfile) {
			_satclkfile->write(os.str().c_str(), os.str().size());
			_satclkfile->flush();
			os.str("");
			//if (isTriple)
			//{
			//	_satclk13file->write(os13.str().c_str(), os13.str().size());
			//	_satclk13file->flush();
			//	os13.str("");
			//}
		}
		else {
			_glog->comment(2, "t_gpcelsqIF::_printout", "have no output file!");
		}
	}

	void t_gpcelsqIF::_print_rtnet_clk(t_glsq* lsq, t_gtime& crt_epoch)
	{
		set<string>::const_iterator sat_end = _sat_list.end();
		set<string>::const_iterator it_sat = _sat_list.begin();

		t_gallpar par = lsq->_x_solve;
		ColumnVector dx = lsq->dx();
		SymmetricMatrix Qx = lsq->Qx();

		ostringstream os;

		// write rtnet clk time
		os << setw (2) << "* " << setw(25) << crt_epoch.str("%Y %C %L %K %O %P") << endl;
		dynamic_cast<t_gtcp*>(_gioout)->run_send(os.str());
		//test
		cout << os.str() << endl;
		os.str("");

		// write rtnet clk data
		for (it_sat; it_sat != sat_end; it_sat++) {
			int idx = par.getParam("", t_gpar::CLK_SAT, *it_sat, crt_epoch, crt_epoch);
			if (idx < 0) continue;
			int idx_ifcb;

			double x = par[idx].value();
			double dclk = dx(idx + 1);

			double clk;
			double ifcb = 0.0;
			if (dclk != 0) {
				clk = x + dclk;
				os << setw(3) << *it_sat << " " << setw(6) << "CLK 1 " << setw(15) << right << fixed << setprecision(6) << clk << endl;
				dynamic_cast<t_gtcp*>(_gioout)->run_send(os.str());
				//test
				cout << os.str() << endl;
				os.str("");
			}
		}

		// write rtnet clk end string
		os << setw(3) << "EOE" << endl;
		dynamic_cast<t_gtcp*>(_gioout)->run_send(os.str());
		//test
		cout << os.str() << endl;
		os.str("");

	}

	void t_gpcelsqIF::set_recover_data_by_resfile(string resfile_name)
	{
		int i = resfile_name.find("_resfile");
		string site = resfile_name.substr(i - 4, 4);
		site = str2upper(site);
		if (_siteres[site]) {
			delete _siteres[site];
			_siteres[site] = nullptr;
		}
		_siteres[site] = new t_gallrecover();

		_glog->comment(t_glog::LOG_LV::LOG_INFO, "t_gpcelsqIF", "set_recover_data_by_resfile", "Reading resfile: " + resfile_name);

		t_gio* gin = new t_gfile;
		gin->path(resfile_name); gin->glog(_glog);
		t_gcoder* read_coder = new t_resfile(nullptr, "", 4096);
		read_coder->clear(); read_coder->glog(_glog); read_coder->path(resfile_name);
		read_coder->add_data("ID1", _siteres[site]);
		gin->coder(read_coder);
		gin->run_read();
		delete gin;
		delete read_coder;

		_glog->comment(t_glog::LOG_LV::LOG_INFO, "t_gpcelsqIF", "set_recover_data_by_resfile", "Reading resfile: Finish!");
	}
}