// history : 
//	jdhuang 2019-05-12 add the iono pars init function in _InitLsq

#include "gproc/gpcelsqALL.h"
#include "gset/gsetpar.h"
#include "gset/gsetout.h"
#include "gset/gsetgen.h"
#include "gutils/ginfolog.h"
#include "gutils/gturboedit.h"
namespace great
{
	t_gpcelsqALL::t_gpcelsqALL()
	{
	}

	t_gpcelsqALL::t_gpcelsqALL(t_gsetbase* set, t_gallproc* data, t_glog* log ) :
		t_glsqproc(set, data, log)
	{

		_bias_pcemodel = new t_gpceprecisemodel(data, _glog, _gset);

		_clk_log = new t_glog();
		_clk_log->mask("clk.log");
		_clk_log->cache_size(99);
		_clk_log->tsys(t_gtime::GPS);
		_clk_log->time_stamp(true);
	}

	t_gpcelsqALL::~t_gpcelsqALL()
	{
		if (_bias_pcemodel)
		{
			delete _bias_pcemodel;
			_bias_pcemodel = nullptr;
		}
		if (_clk_log)
		{
			delete _clk_log;
			_clk_log = nullptr;
		}
	}

	bool t_gpcelsqALL::ProcessBatch(t_gallproc * data, const t_gtime & beg, const t_gtime & end)
	{
		string func_id = "t_gpcelsqALL";
		if (!data ||
			beg > end)
		{
			write_log_info(_glog, 0, "ERROR", "beg > end or no proc data");
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

		string sat_ref = dynamic_cast<t_gsetproc*>(_gset)->ref_clk();
		while (_crt_time <= _end_time)
		{
			cerr << _crt_time.str_ymdhms() << ":" << endl;
			write_log_info(_glog, 1, "NOTE", _crt_time.str_ymdhms("Processing epoch "));

			/* get observation */
			vector<t_gsatdata> crt_obs = _gall_obs->obs(_rec_list, _crt_time);
			
			bool epoch_valid = _processOneEpoch(_crt_time, crt_obs);
			if (!epoch_valid)
			{
				write_log_info(_glog, 1, "ERROR" + _crt_time.str_ymdhms("Processing failed in epoch "));
				_crt_time = _crt_time + _obs_intv;
				continue;
			}
			if (sat_ref == "") _lsq->lsq_clk_constraint();
			try
			{
				if (_lsq->mode() == gnut::EPO && _crt_time < _end_time)
				{
					t_glsq* lsqtmp = new t_glsq(*_lsq);
					if (!_IRC)
					{
						lsqtmp->lsq_sysbias_constraint();
					}
					lsqtmp->solve_NEQ();
					if (_satclkfile) {
						if (_frequency == 3 && _satclk13file)
						{
							_print_clk_solved(lsqtmp, _crt_time, true);
						}
						else
						{
							_print_clk_solved(lsqtmp, _crt_time);
						}
					}
				}
			}
			catch (...)
			{
				_crt_time = _crt_time + _obs_intv;
				continue;
			}
			_crt_time = _crt_time + _obs_intv;
		}

		try
		{
			if (sat_ref == "") _lsq->lsq_clk_constraint();
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


	bool t_gpcelsqALL::_initLsqProdData(t_gallprod * data)
	{
		return false;
	}

	bool t_gpcelsqALL::_initLsqProcPars(t_glsq* lsq)
	{
		string class_id = "t_gpcelsqALL";
		string funct_id = "_initLsqProcPars";

		if (_rec_list.empty() || _sat_list.empty())
		{
			this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "no rec or sat");
			return false;
		}

		// Crd pars
		bool rec_crd_valid = _init_rec_crd_pars(lsq);
		bool sat_clk_vaild = _init_sat_clk_pars(lsq);
		bool sat_ifcb_valid = false;
		if (_frequency == 3)
		{
			sat_ifcb_valid = _init_sat_ifcb_pars(lsq);
		}
		bool rec_clk_valid = _init_rec_clk_pars(lsq);
		bool rec_clk13_valid = false;
		if (_frequency == 3)
		{
			rec_clk13_valid = _init_rec_clk13_pars(lsq);
		}
		bool trop_valid = _init_rec_trop_pars(lsq);
		bool ion_valid = _init_gns_sion_pars(lsq);
		bool isb_valid = _init_sys_bia_pars(lsq);

		if (!rec_crd_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "rec_crd_valid  is false"); return false; }
		if (!sat_clk_vaild) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "sat_clk_valid  is false"); return false; }
		if (_frequency == 3 && !sat_ifcb_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "sat_ifcb_valid  is false"); return false; }
		if (!rec_clk_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "rec_clk_valid  is false"); return false; }
		if (!trop_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "trop_valid  is false"); return false; }
		if (!ion_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "ion_valid  is false"); return false; }
		if (!isb_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "isb_valid  is false"); return false; }

		write_log_info(_glog, 1, "INIT", "init the pars estimated in PCE_LSQ_ALL.");
		return true;

		//if (clk_valid && 
		//	trop_valid)
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

	bool t_gpcelsqALL::_initLsqProcData(t_gallproc * data)
	{
		gtrace("t_gppplsqIF::_initLsqProcData");

		// prepare obs data
		_gall_bias = dynamic_cast<t_gallbias*>((*data)[t_gdata::ALLBIAS]);
		_gall_obs = dynamic_cast<t_gallobs*>((*data)[t_gdata::ALLOBS]);
		_gall_obj = dynamic_cast<t_gallobj*>((*data)[t_gdata::ALLOBJ]);
		_gall_rcv = dynamic_cast<t_gallrecover*>((*data)[t_gdata::ALLRECOVER]);


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

	//bool t_gpcelsqALL::_processOneRecOneSat(const t_gtime & crt_epoch,
	//	const std::string & rec, const std::string & sat, 
	//	t_gsatdata & crd_obs, int & obs_num)
	//{
	//	bool proc_valid = t_gpcelsqALL::_processOneRecOneSat(crt_epoch, rec, sat,crd_obs, obs_num);

	//	if (!proc_valid)
	//	{
	//		write_log_info(_log, 2, "ERROR", crt_epoch.str_ymdhms("can not process the sat of " +
	//			sat + " and the site of " + rec + " in the time ",false));
	//		return false;
	//	}

	//	return true;
	//}

	bool t_gpcelsqALL::_processOneEpoch(const t_gtime& crt_epoch, std::vector<t_gsatdata>& crt_obs)
	{
		string class_id = "t_gpcelsqALL";
		string funct_id = "_processOneEpoch";

		// updata pars
		if (!_slip12 || crt_obs.empty()) return false;
		//// for debug
		ColumnVector W = _lsq->W();					///< BTPL Matrix
		SymmetricMatrix NEQ = _lsq->NEQ();

		bool updata_valid;
		updata_valid = _lsq->update_parameter(_IRC, crt_epoch, crt_obs, false);
		if (!updata_valid ||
			crt_obs.empty())
		{
			write_log_info(_glog, 1, crt_epoch.str(), ": no observation found ");
			return false;
		}

		int iObs = 0;
		set<string>::iterator rec_beg = _rec_list.begin();
		set<string>::iterator rec_end = _rec_list.end();
		vector<t_gsatdata>::iterator obs_beg = crt_obs.begin();
		vector<t_gsatdata>::iterator obs_end = crt_obs.end();
		for (set<string>::iterator it = rec_beg; it != rec_end; it++)
		{
			vector<t_gsatdata> crt_rec_obs;
			string crt_rec = *it;

			// claer equ before
			_crt_equ.clear_allequ();

			// reset init value of coord parameters for kin mode
			if (_crd_est == KIN && _kin_xyz[_crt_rec].find(crt_epoch) != _kin_xyz[_crt_rec].end())
			{
				int idx = _lsq->_x_solve.getParam(crt_rec, t_gpar::CRD_X, "");
				_lsq->_x_solve.setParValue(idx, _kin_xyz[crt_rec][crt_epoch].crd(0));
				idx = _lsq->_x_solve.getParam(crt_rec, t_gpar::CRD_Y, "");
				_lsq->_x_solve.setParValue(idx, _kin_xyz[crt_rec][crt_epoch].crd(1));
				idx = _lsq->_x_solve.getParam(crt_rec, t_gpar::CRD_Z, "");
				_lsq->_x_solve.setParValue(idx, _kin_xyz[crt_rec][crt_epoch].crd(2));
			}

			// get all the obs for the rec
			bool select_valid = _select_rec_obs(crt_rec, crt_epoch, crt_obs, crt_rec_obs);
			if (!select_valid)
			{
				_glog->comment(t_glog::LOG_LV::LOG_INFO, class_id, funct_id, crt_epoch.str_mjdsod(crt_rec + " has no useful data"));
				continue;
			}

			// Process this site and get the equations
			bool proc_site = _processOneRec(crt_epoch, crt_rec, crt_rec_obs);
			if (!proc_site)
			{
				cout << crt_epoch.str_ymdhms(crt_rec + " has no equations ", false, false) << endl;
				_glog->comment(t_glog::LOG_LV::LOG_INFO, class_id, funct_id, crt_epoch.str_mjdsod(crt_rec + " has no equations"));
				continue;
			}

			// add the new equations
			try
			{
				_lsq->add_equation(_IRC, _crt_equ, crt_epoch);
			}
			catch (...)
			{
				_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, crt_epoch.str_mjdsod("add_equation throw error"));
				return false;
			}

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

	bool t_gpcelsqALL::_processOneRec(const t_gtime& crt_epoch, const std::string& crt_rec, std::vector<t_gsatdata>& crt_obs)
	{
		string class_id = "t_gpcelsqALL";
		string funct_id = "_processOneRec";

		if (crt_obs.empty() ||
			crt_rec.empty())
		{
			_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, crt_epoch.str_mjdsod("This epoch have no useful data: " + crt_rec));
			return false;
		}

		// update crt_rec in priviate
		_crt_rec = crt_rec;

		// loop for Combine equation
		vector<double> codeOmc;
		bool while_valid = false;
		do
		{
			_crt_equ.clear_allequ();
			for (auto& it : crt_obs)
			{
				string crt_sat = it.sat();
				if (_sat_list.empty() ||
					_sat_list.find(crt_sat) == _sat_list.end())
				{
					continue;
				}

				if (!t_gpcelsqALL::_processOneRecOneSat(crt_epoch, crt_rec, crt_sat, it))
				{
					_glog->comment(t_glog::LOG_LV::LOG_WARN, class_id, funct_id, crt_epoch.str_mjdsod("can not process the rec : " + crt_rec + " and the sat : " + crt_sat + " "));
					continue;
				}
			}

			//reset recclk omc
			codeOmc = _crt_equ.get_codeomc();
			int idx = _lsq->_x_solve.getParam(crt_rec, t_gpar::CLK, "");
			double recclk = 0.0;
			double sigma = 0.0;
			int k = Check_Range(codeOmc, recclk, sigma);

			// reset recclk
			if (idx >= 0 && !while_valid && (fabs(recclk / CLIGHT) > 1e-6 ||
				(k >= 1 && fabs(_lsq->_x_solve[idx].value() / CLIGHT) < 1e-15)))
			{
				recclk += _lsq->_x_solve[idx].value();
				_lsq->_x_solve[idx].value(recclk);
				_bias_pcemodel->update_obj_clk(crt_rec, crt_epoch, recclk / CLIGHT);

				codeOmc.clear();
				while_valid = true;
			}
			else
			{
				while_valid = false;
			}
		} while (while_valid);

		// give the satclk initial value
		if (!_init_sat_clk()) return false;

		return true;
	}

	bool t_gpcelsqALL::_processOneRecOneSat(const t_gtime& crt_epoch,
		const std::string& rec,
		const std::string& sat,
		t_gsatdata& crt_obs)
	{
		string class_id = "t_gpcelsqALL";
		string funct_id = "_processOneRecOneSat";

#ifdef DEBUG_LSQ_ALL
		_glog->comment(t_glog::LOG_LV::LOG_INFO, class_id, funct_id, crt_epoch.str_mjdsod("process the " + rec + " + " + sat));
#endif

		bool code_valid = _bias_pcemodel->_cmb_equ(crt_epoch, crt_obs, _lsq->_x_solve, _crt_equ);
		if (!code_valid)
		{
			_glog->comment(t_glog::LOG_LV::LOG_WARN, class_id, funct_id, crt_epoch.str_mjdsod("can not combine equ for code in "));
			return false;
		}
		return true;
	}

	bool t_gpcelsqALL::GenerateProduct()
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

	void t_gpcelsqALL::_print_clk_solved(t_glsq* lsq, t_gtime& crt_epoch, bool isTriple)
	{
		set<string>::const_iterator sat_end = _sat_list.end();;
		set<string>::const_iterator it_sat = _sat_list.begin();

		t_gallpar par = lsq->_x_solve;
		ColumnVector dx = lsq->dx();
		SymmetricMatrix Qx = lsq->Qx();

		ostringstream os;
		ostringstream os13;
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
			_glog->comment(2, "t_gpcelsqALL::_write_clk_header", "have no output file!");
		}

		// write sp3 clk file data
		double satclk;
		for (it_sat; it_sat != sat_end; it_sat++) {
			int idx = par.getParam("", t_gpar::CLK_SAT, *it_sat, crt_epoch, crt_epoch);
			int idx_ifcb = par.getParam("", t_gpar::IFCB_F3, *it_sat, crt_epoch, crt_epoch);

			double x = par[idx].value();
			double dclk = dx(idx + 1);
			double x_ifcb;
			double difcb;
			if (idx_ifcb > -1)
			{
				x_ifcb = par[idx_ifcb].value();
				difcb = dx(idx_ifcb + 1);
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
		}

		// Print sp3 clk file data
		if (_satclkfile) {
			_satclkfile->write(os.str().c_str(), os.str().size());
			_satclkfile->flush();
			os.str("");
		}
		else {
			_glog->comment(2, "t_gpcelsqALL::_printout", "have no output file!");
		}
	}


}