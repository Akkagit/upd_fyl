#ifdef USE_LIDAR
#include "gmultign_lidar.h"
#include "gproc/gpreproc.h"
#include "gutils/gtimesync.h"
#include <iostream>
using namespace gsins;
using namespace std;


///******************************************* t_gmultign ***********************************************/
t_gmultign_lidar::t_gmultign_lidar(string site, string site_base, t_gsetbase* gset, t_gallproc *allproc) :
	t_gintegration(site, site_base, gset, allproc),
	t_gsinskf(gset),
	t_gpvtflt(site, site_base, gset, allproc),
	t_glidar(gset),
	t_gspp(site, gset)
{
	mIsFirstLidar = true;
	_lidardata = nullptr;
	_framecount = 0;
	_lidarOdo.scanPeriod_ = dynamic_cast<t_gsetlidar*>(gset)->scan_period();
}

void gsins::t_gmultign_lidar::setgallparuse(bool status)
{
	is_allTC = status;
}

void gsins::t_gmultign_lidar::Add_LIDAR(t_glidardata * lidar)
{
	_lidardata = lidar;
}


bool gsins::t_gmultign_lidar::_lidar_time_valid(t_gtime lidart, t_gtime inst)
{
	gtrace("t_gintegration::_time_valid");

	bool res_valid = false;
	//out << "diff:" << abs(inst.diff(lidart)) << endl;
	if (abs(inst.diff(lidart)) < _shm.delay) res_valid = true;
	else res_valid = false;

	return res_valid;
}


int gsins::t_gmultign_lidar::processBatch(const t_gtime & beg, const t_gtime & end)
{
	gtrace("t_gintegration::processBatch");
#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif 
	//_gmutex.lock();

	if (_grec == nullptr) {
		ostringstream os;
		os << "ERROR: No object found (" << _site << "). Processing terminated!!! " << beg.str_ymdhms() << " -> " << end.str_ymdhms() << endl;
		if (_log) _log->comment(0, "t_gintegration", os.str());
		//	_gmutex.unlock(); 
		return -1;
	}

	int sign = 1; if (!_beg_end) sign = -1;
	t_gpvtflt::InitProc(beg, end);

	if (!this->_init()) {
		if (_log) _log->comment(0, "t_gintegration", "init failed");
		return -1;
	}

	// need execute frequently in real-time mode
	_gobs->setepoches(_site);
	bool realtime = dynamic_cast<t_gsetproc*>(_set)->realtime();
	// time align
	if (!realtime) {
		if (_ins_beg < _gnss_beg) {
			_ins_crt = _imudata->erase_bef(_gnss_beg);
		}
		else
		{
			_gmutex.unlock();
			t_gpvtflt::processBatch(_gnss_crt, _ins_crt, false);
			_gmutex.lock();
			// next gnss time 
			if (_gnss_crt < _ins_crt) {
				int nEpo = round(_ins_crt.diff(_gnss_crt) / (sign * _sampling) + 0.5);
				if (_sampling > 1) _gnss_crt.add_secs(int(sign * _sampling * nEpo));  // =<1Hz data
				else               _gnss_crt.add_dsec(sign * _sampling * nEpo);       //  >1Hz data
			}
			//while (_gnss_crt < _ins_crt)
			//{
			//	if (_sampling > 1) _gnss_crt.add_secs(int(sign * _sampling));  // =<1Hz data
			//	else               _gnss_crt.add_dsec(sign * _sampling);       //  >1Hz data
			//}
		}
	}

	std::cerr << "\n" << _site << ": Start GNSS/INS Integration Processing: " << _ins_beg.str_ymdhms() << " " << _ins_end.str_ymdhms() << endl;

	t_gposdata::data_pos posdata;
	//to do remove this one
	bool time_loop = true;
	// this bool is to eliminate the variance if mode is VIO
	if (_ign_type == IGN_TYPE::MULTIGN_LIDAR_TC) this->setgallparuse(true);
	
	bool is_GNSS = ((_ign_type != IGN_TYPE::LIDAR_INS) && (_ign_type != IGN_TYPE::PURE_INS));
	bool is_LIDAR = ((_ign_type == IGN_TYPE::LIDAR_INS) || (_ign_type == IGN_TYPE::MULTIGN_LIDAR_LC)
		|| (_ign_type == IGN_TYPE::MULTIGN_LIDAR_TC) || (_ign_type == IGN_TYPE::MULTIGN_LIDAR_STC));
	//is_GNSS = false;
	//cout << "is_GNSS:" << is_GNSS << endl;
	//cout << "is_LIDAR:" << is_LIDAR << endl;
	while (time_loop)
	{
		//lsy
		if (_ins_crt > _ins_end || _gnss_crt > _gnss_end)
		{
			cout << " Processing Finish!" << endl;
			if (_log) _log->comment(0, "t_gintegration", "Processing Finish!");
			time_loop = false; break;
		}
		if (!_imudata->load(_wm, _vm, _shm.t, _shm.ts, _shm.nSamples)) { cerr << "\n Processing Finished!" << endl; break; }
		_ins_crt = t_gtime(_gnss_crt.gwk(), _shm.t);
		if (!_aligned)
		{
			if (_time_valid(_gnss_crt, _ins_crt))Flag = _getPOS(posdata);
			_aligned = cascaded_align(posdata.pos, posdata.vn);
		}
		else
		{
			if ((_ign_type == IGN_TYPE::TCI || _ign_type==IGN_TYPE::MULTIGN_LIDAR_TC) && !_initial_merge) _merge_init();
			sins.Update(_wm, _vm, _shm);
			if (_ign_type == PURE_INS) {
				t_gsinskf::write();
				continue;
			}
			time_update(sins.nts, 0.0);

			
			if (_ign_type == IGN_TYPE::LIDAR_INS || _ign_type == IGN_TYPE::MULTIGN_LIDAR_LC || 
				_ign_type==IGN_TYPE::MULTIGN_LIDAR_TC || _ign_type==IGN_TYPE::MULTIGN_LIDAR_STC)
			{
				_lidardata->load(_shm.t, _lidar_t, _lidar_path);
				//cout << "imu_t:" << setprecision(8) << _shm.t << "   lidar_t:" << setprecision(8) << _lidar_t << endl;
				_lidar_crt = t_gtime(_gnss_crt.gwk(), _lidar_t);
				//lidarproc->load_imuobs(_shm.t, _wm, _vm);
			}
			// todo need more judge
			bool flag_gnss = (is_GNSS && (_time_valid(_gnss_crt, _ins_crt)));
			bool flag_lidar = (is_LIDAR && (_lidar_time_valid(_lidar_crt, _ins_crt)));
			//cout << "ins_time:" << setprecision(8) << (_ins_crt.sow() + _ins_crt.dsec()) << endl;
			//cout << "lidar_time:" << setprecision(8) << (_lidar_crt.sow() + _lidar_crt.dsec())<<" is_LIDAR:"<<is_LIDAR << endl;
			//std::cout << "flag_lidar:" << flag_lidar << endl;
			if (flag_lidar || flag_gnss)
			{
				std::cout << endl;
				std::cout << "imu_time:" << setprecision(10) << _shm.t << endl;
				std::cout << "lidar_now:" << setprecision(10) << _lidar_t << "  " << flag_lidar << endl;
				std::cout << "gns_now:" << setprecision(10) << (_gnss_crt.sow() + _gnss_crt.dsec()) << "  " << flag_gnss << endl;
			}
			
			if (flag_gnss)
			{
				if (_ign_type == IGN_TYPE::LCI)
				{
					Flag = _getPOS(posdata);
					if (Flag != NO_MEAS) _meas_update();
				}
				else if (_ign_type == IGN_TYPE::TCI)
				{

					_timeUpdate(_gnss_crt);
					Flag = MEAS_TYPE(this->_processGNSSEpoch(_gnss_crt));
					if (Flag == NO_MEAS) { _reset(); }

				}
				else if (_ign_type == IGN_TYPE::STCI)
				{
					if (_valid_ins_constraint())
					{
						Eigen::Vector3d XYZ_INS = sins.pos_ecef + sins.Ceb*lever;
						t_gtriple XYZ(XYZ_INS(0), XYZ_INS(1), XYZ_INS(2));
						double ins_error = SQRT(Pk(6, 6) + Pk(7, 7) + Pk(8, 8));
						//t_gtriple XYZ_RMS(SQRT(Pk(6, 6)), SQRT(Pk(7, 7)), SQRT(Pk(8, 8)));
						//t_gtriple XYZ_RMS(ins_error, ins_error, ins_error);
						t_gtriple XYZ_RMS(0.1, 0.1, 0.1);
						_external_pos(XYZ, XYZ_RMS);
					}

					Flag = _getPOS(posdata);
					if (Flag != NO_MEAS) {
						_meas_update();
						//if (outlier_detect()) { _reset(); Flag = NO_MEAS; }
					}
				}
				else if (_ign_type == IGN_TYPE::MULTIGN_LIDAR_LC)
				{
					Flag = _getPOS(posdata);
					if (Flag != NO_MEAS) _meas_update();
				}
				else if (_ign_type == IGN_TYPE::MULTIGN_LIDAR_TC)
				{
					cout << "gnss update" << endl;
					_timeUpdate(_gnss_crt);

					//whd test
					//Eigen::Vector3d before_p = imu_state.position;
					Flag = MEAS_TYPE(this->_processGNSSEpoch(_gnss_crt));
					//Eigen::Vector3d after_p = imu_state.position;
					/*if ((before_p - after_p).norm() > 0.2)
					{
						pos_jump = true;
						cout << (before_p - after_p).norm() << endl;
						getchar();
					}
					else
						pos_jump = false;*/

					if (Flag == NO_MEAS) { _reset(); }
				}
				else if (_ign_type == IGN_TYPE::MULTIGN_LIDAR_STC)
				{
					if (_valid_ins_constraint())
					{
						Eigen::Vector3d XYZ_INS = sins.pos_ecef + sins.Ceb*lever;
						t_gtriple XYZ(XYZ_INS(0), XYZ_INS(1), XYZ_INS(2));
						double ins_error = SQRT(Pk(6, 6) + Pk(7, 7) + Pk(8, 8));
						//t_gtriple XYZ_RMS(SQRT(Pk(6, 6)), SQRT(Pk(7, 7)), SQRT(Pk(8, 8)));
						//t_gtriple XYZ_RMS(ins_error, ins_error, ins_error);
						t_gtriple XYZ_RMS(0.1, 0.1, 0.1);
						_external_pos(XYZ, XYZ_RMS);
					}

					Flag = _getPOS(posdata);
					if (Flag != NO_MEAS) {
						_mult_meas_update();
						//if (outlier_detect()) { _reset(); Flag = NO_MEAS; }
					}
				}
				else
				{
					cerr << "flag_gnss=true is a wrong mode!" << endl;
				}

				if (Flag != NO_MEAS) {
					feedback(glv.INF);
					sins.pure_ins_time = 0.0;
					cerr << "\r" << _ins_crt.str_ymdhms();
					if (_log) _log->comment(0, "t_gintegration",
						_ins_crt.str_ymdhms(" epoch ") + "Processing GNSS/INS Integration Success!");
				}

			}
			if (flag_lidar)
			{
				//cout << flag_lidar << "," << _ign_type << endl;
				// todo  merge to one function
				if (_ign_type == IGN_TYPE::LIDAR_INS)
				{
					//cout << flag_lidar << "," << _ign_type << endl;
					//cout << _lidar_path.t << "," << _lidar_path.lidar_path << endl;
				    lidarproc->load_lidarobs(_lidar_t, _lidar_path);
					Flag = MEAS_TYPE(this->_processLidarEpoch(_lidar_crt));
				}
				else if (_ign_type == IGN_TYPE::MULTIGN_LIDAR_LC)
				{
					lidarproc->load_lidarobs(_lidar_t, _lidar_path);
					Flag = MEAS_TYPE(this->_processLidarEpoch(_lidar_crt));
				}
				else if (_ign_type == IGN_TYPE::MULTIGN_LIDAR_TC)
				{
					lidarproc->load_lidarobs(_lidar_t, _lidar_path);
					Flag = MEAS_TYPE(this->_processLidarEpoch(_lidar_crt));
				}
				else if (_ign_type == IGN_TYPE::MULTIGN_LIDAR_STC)
				{
					lidarproc->load_lidarobs(_lidar_t, _lidar_path);
					Flag = MEAS_TYPE(this->_processLidarEpoch(_lidar_crt));
				}
				else
				{
					cerr << "flag_lidar=true is a wrong mode!" << endl;
				}

				if (_ins_crt <= _ins_end&&use_map&&_lidarmap.systemInited_)
				{
					imu_state.time = _ins_crt.sow() + _ins_crt.dsec();
					publisher.UpdateNewState(imu_state, lidar_states, _lidarmap.first_t_l_e, _lidarmap.submap_corner);
				}
				else if (_ins_crt <= _ins_end && use_scan&&_lidarOdo.systemInited_)
				{
					imu_state.time = _ins_crt.sow() + _ins_crt.dsec();
					publisher.UpdateNewState(imu_state, lidar_states,Eigen::Vector3d::Zero(), _lidarOdo.test_lastcorner);
				}
			}

			

			if(flag_gnss||flag_lidar)
			{
				t_gsinskf::write();
				Eigen::Vector3d Car_pos = Geod2Cart(sins.pos, false);
				Eigen::Vector3d Ve = t_gbase::Cen(sins.pos)*(sins.vn);
				std::cout << fixed << setprecision(4) <<sins.t;
				std::cout << fixed << setprecision(3) <<
					" " << Car_pos(0) << " "
					" " << Car_pos(1) << " "
					" " << Car_pos(2) << " ";
				std::cout << fixed << setprecision(3) <<
					" " << sins.vn(0) << " "
					" " << sins.vn(1) << " "
					" " << sins.vn(2) << " " << endl;
			}


			//use for glfw 
			imu_state.orientation = Eigen::Quaterniond(t_gbase::q2mat(sins.qnb));
			imu_state.position = Geod2Cart(sins.pos, false);

			/*if (lidar_states.size() > 0)
			{
				cout << "lidar state:" <<lidar_states.size()<< endl;
				getchar();
			}*/
			if (_ins_crt <= _ins_end&& _ign_type==IGN_TYPE::TCI)
			{
				//only show the odometry result got from imu time update
				imu_state.time = _ins_crt.sow()+_ins_crt.dsec();
				publisher.UpdateNewState(imu_state, lidar_states);
			}
			
		}
		while (_ins_crt.diff(_gnss_crt) > _sampling)
		{
			if (_sampling > 1) _gnss_crt.add_secs(int(_sampling));  // =<1Hz data
			else               _gnss_crt.add_dsec(_sampling);       //  >1Hz data
		}

	}


	_gmutex.unlock();
	return 0;
}

int gsins::t_gmultign_lidar::_processGNSSEpoch(const t_gtime & runEpoch)
{
	gtrace("t_gintegration::_processEpoch");
#ifdef BMUTEX   0
	boost::mutex::scoped_lock lock(_mutex);
#endif 
	//	_gmutex.lock();

#ifdef DEBUG
	cout << "t_gintegration: Epoch - full data size: " << _data.size() << " "
		<< runEpoch.str_hms() << endl;
#endif

	vector<t_gsatdata>::iterator ItException = _data.end();
	ColumnVector Threshold;
	if (_grec == nullptr) {
		if (_log) _log->comment(1, "t_gintegration", "No receiver settings available!!!");
		return NO_MEAS;
	}

	_epoch = runEpoch;
	t_gtriple crdapr = _grec->crd_arp(_epoch);

	if (double_eq(crdapr[0], 0.0) && double_eq(crdapr[1], 0.0) &&
		double_eq(crdapr[2], 0.0)) {
		_valid_crd_xml = false;
	}
	else {
		_valid_crd_xml = true;
	}

	if (!_valid_crd_xml) _sig_init_crd = 100.0;

	Matrix A;
	SymmetricMatrix P;
	ColumnVector l, dx;
	SymmetricMatrix Qsav, QsavBP;
	vector<t_gsatdata>::iterator it;
	t_gallpar XsavBP;
	ColumnVector v_norm, v_orig;
	string outlier = "";
	double vtpv;
	int nobs_total, npar_number;

	_cntrep = 0;    // number of iterations caused by outliers

	do
	{

		// erase sat satellite because of outlier
		it = _data.begin();
		while (it != _data.end()) {
			if (it->sat() == outlier) {
				it = _data.erase(it); continue;
			}
			it++;
		}

		if (_valid_ins_constraint()) {
			Eigen::Vector3d XYZ_INS = sins.pos_ecef + sins.Ceb*lever;
			t_gtriple XYZ(XYZ_INS(0), XYZ_INS(1), XYZ_INS(2));
			_external_pos(XYZ, t_gtriple());
		}



		if (_prepareData() < 0) {
			_predict();//todo:add time gap
			return NO_MEAS;
		}
		_pos_constrain = false;
		QsavBP = _Qx; XsavBP = _param;
		_predict();

		_initialized = true;

		if (_data.size() < _minsat) {
			if (_log) _log->comment(2, "t_gintegration", "Not enough visible satellites!");
			t_gsppflt::_restore(QsavBP, XsavBP);
			return NO_MEAS;
		}

		// define a number of measurements 
		unsigned int nObs = _data.size();
		unsigned int mult = 1;
		if (_observ == RAW_DOUBLE) { mult = 2; nObs *= 2; }
		if (_observ == RAW_ALL) { mult = 2; nObs *= 5; } // reservation for 5 freq - not used raws will be removed
		if (_phase) { mult *= 2; nObs *= 2; }  // code + phase

		unsigned int nPar = _param.parNumber();

		A.ReSize(nObs, nPar); A = 0.0;
		l.ReSize(nObs);       l = 0.0;
		P.ReSize(nObs);       P = 0.0;
		dx.ReSize(nPar);     dx = 0.0;
		unsigned int iobs = 1;
		_frqNum.clear(); _obs_index.clear();

		Matrix A_new;
		SymmetricMatrix P_new;
		ColumnVector l_new, dx_new;
		dx_new.ReSize(nPar);    dx_new = 0.0;
		unsigned int iobs_new = 1;
		_frqNum.clear(); _obs_index.clear();
		if (_isBase)
		{
			dynamic_cast<t_gcombDD*>(&(*_base_model))->set_base_data(&_data_base);
			dynamic_cast<t_gcombDD*>(&(*_base_model))->set_rec_info(_gallobj->obj(_site_base)->crd(_epoch), _vBanc[3], _vBanc_base[3]);
		}
		// use combmodel
		t_glsqEquationMatrix equ;
		it = _data.begin();
		for (it = _data.begin(); it != _data.end(); ) {
			int irc = _base_model->cmb_equ(_epoch, _param, *it, equ);
			if (irc < 0) { it = _data.erase(it); continue; }
			if (_observ == IONO_FREE) _combineMW(*it);
			++it;
		}
		iobs_new = equ.num_equ();
		// tranform equ into A, P, l
		_equationToMatrix(equ, A_new, P_new, l_new, nPar);

		// generate obs_index
		_obs_index.clear();
		_generateObsIndex(equ);

		if (iobs_new < _minsat * mult) {
			if (_log) _log->comment(2, "gpvtflt", "Not enough processable observations!");
			t_gsppflt::_restore(QsavBP, XsavBP);
			return NO_MEAS;
		}

		A = A_new;
		P = P_new;
		l = l_new;
		dx = dx_new;

		if (_isBase || _sd_sat) {
			if (_combineDD(A, P, l) < 0) {
				return NO_MEAS;
			}
		}
		Qsav = _Qx;

		if (_merge(A) < 0) {
			dx = 0.0; _Qx = Qsav;
			if (_log) _log->comment(1, "gintegration", "merge failed!");
			return NO_MEAS;
		}

		try
		{
			_filter->update(A, P, l, dx, _Qx);
		}
		catch (...)
		{
			if (_log) _log->comment(1, "gintegration", "meas_update failed!");
			dx = 0.0; _Qx = Qsav;
			Flag = NO_MEAS;
			_reset();
			return Flag;
		}

		// increasing variance after update in case of introducing new ambiguity
		for (size_t iPar = 0; iPar < _param.parNumber(); iPar++) {
			if (_param[iPar].parType == t_gpar::AMB_IF) {
				string sat = _param[iPar].prn;
				if (_newAMB.find(sat) != _newAMB.end() && _cntrep == 1) {
					if (_newAMB[sat] == 1) _Qx(iPar + 1, iPar + 1) += 10;
					if (_newAMB[sat] == 2 && _Qx(iPar + 1, iPar + 1) > 0.01) _Qx(iPar + 1, iPar + 1) += 1;
					_newAMB[sat]++;
				}
			}
		}

		//if (_isBase)
		//{
		//	if (_observ == RAW_ALL || _observ == RAW_SINGLE)Flag = DD_UC_OBS_MEAS;
		//	if (_observ == IONO_FREE)Flag = DD_IF_OBS_MEAS;
		//}
		//else
		//{
		//	if (_observ == RAW_ALL || _observ == RAW_SINGLE)Flag = UD_UC_OBS_MEAS;
		//	if (_observ == IONO_FREE)Flag = UD_IF_OBS_MEAS;
		//}

		Flag = GNSS_MEAS;

		// post-fit residuals
		v_orig = l - A * dx;
		ColumnVector vtPv = v_orig.t() * P * v_orig;
		vtpv = vtPv(1);
		nobs_total = A.Nrows();
		npar_number = A.ncols();
		int freedom = A.Nrows() - A.Ncols();
		if (freedom < 0)freedom = 1;
		_sig_unit = vtPv(1) / freedom;

		// normalized post-fit residuals
		Matrix Qv = A * _Qx * A.t() + P.i();
		v_norm.ReSize(v_orig.Nrows());
		for (int i = 1; i <= v_norm.Nrows(); i++) {
			v_norm(i) = sqrt(1 / Qv(i, i)) * v_orig(i);
		}


	} while (_outlierDetect(v_norm, Qsav, outlier) != 0);

	_filter->add_data(_param, dx, _Qx, _sig_unit, Qsav);
	_filter->add_data(A, P, l);
	_filter->add_data(vtpv, nobs_total, npar_number);

	if (_gnss_feedback(dx) < 0) {
		if (_log) _log->comment(1, "gintegration", "meas_update failed!");
		dx = 0.0; _Qx = Qsav;
		Flag = NO_MEAS;
		return Flag;
	}
	//cam parameter change
	//cout << "GNSS_Meas_UP:" << endl;
#ifdef USE_LIDAR
	if (_ign_type == IGN_TYPE::MULTIGN_LIDAR_TC)
	{
		if (_lidar_feedback(dx) < 0)
		{
			cerr << "GNSS MeasUpdate: lidar feed back wrong" << endl;
		}
	}
#endif
	Phik = Eigen::MatrixXd::Identity(15, 15);

	_amb_resolution();

	//	int ins_size = dynamic_cast<t_gsetign*>(_set)->nq();
	//	if (_amb_state)
	//	{
	//		Xk = Columns2VectorXd(_filter->dx()).block(0, 0, ins_size, 1);
	//	}
	//	else
	//	{
	//		Xk = Columns2VectorXd(dx).block(0, 0, ins_size, 1);
	//	}
		// todo: the fixed-dx must be transmited to t_gsinskf::write
	if (_isBase || dynamic_cast<t_gsetamb*>(_set)->fix_mode() != t_gsetamb::FIX_MODE::NO)
	{
		t_gsins sins_out = sins;
		int ins_size = dynamic_cast<t_gsetign*>(_set)->nq();
		// output fixed
		if (_amb_state)
		{
			Xk = Columns2VectorXd(_filter->dx()).block(0, 0, ins_size, 1);
			feedback(glv.INF);
			t_gsinskf::write();
			_prt_ins_kml();
		}
		// transmit float
		sins = sins_out;
		Xk = Columns2VectorXd(dx).block(0, 0, ins_size, 1);


		//Xk = Columns2VectorXd(_filter->dx()).block(0, 0, ins_size, 1);
		//Pk = BaseMatrix2Eigen(_filter->Qx()).block(0, 0, ins_size, ins_size);
		//_param = _param_fixed;
		//_Qx = _filter->Qx();
	}
	//	_gmutex.unlock();
	return Flag;
}
	

int gsins::t_gmultign_lidar::_processLidarEpoch(const t_gtime & runEpoch)
{
	double tmp = runEpoch.sow() + runEpoch.dsec();
	_lidarframe = lidarproc->ProcessBatch(tmp);

	//cout << "test:" << _lidarframe.LessSharp.size() << endl;
	//bool use_scan,use_map;
	//int map_gap;

	if (!use_scan && !use_map)
	{
		std::cout << "Warning: if lidar measurement is needed, you should chose at least one mode from LIDAR_SCAN or LIDAR_MAP!" << std::endl;
		getchar();
		return LIDAR_MEAS;
	}
	if (_lidarframe.fullCloud.points.size() <= 0)
	{
		std::cout << "Error: There are no lidar points read in!" << endl;
		getchar();
		return LIDAR_MEAS;
	}

	imu_state_id = t_glidar::next_id++; //put this ahead
	if (mIsFirstLidar)
	{
		double run_epoch = runEpoch.sow() + runEpoch.dsec();
		double dt = sins.t - run_epoch;
		
		//compensate due to time 
		Eigen::Vector3d BLH = sins.pos - sins.eth.v2dp(sins.vn, dt);//BLH
		Eigen::Vector3d XYZ = gsins::Geod2Cart(BLH, false);
		Eigen::Vector3d imu_vel = sins.vn - sins.an*dt;//ENU
		const Eigen::Matrix3d& R_i_l = R_lidar_imu.transpose();
		const Eigen::Vector3d& t_l_i = t_lidar_imu;
		const Eigen::Matrix3d& R_i_n = t_gbase::q2mat(sins.qnb);
		Eigen::Matrix3d R_n_i = R_i_n.transpose();
		Eigen::Matrix3d R_n_e = t_gbase::Cen(BLH);
		Eigen::Matrix3d R_e_i = R_n_i * R_n_e.transpose();
		Eigen::Vector3d t_l_e = XYZ + R_n_e * R_i_n*t_l_i;
		Eigen::Matrix3d R_e_l = R_i_l * R_e_i;

		_lidarframe.R_l_e = R_e_l.transpose();
		_lidarframe.t_l_e = t_l_e;
		_lidarframe.empty = false;

		_lastframe = _lidarframe;
		mIsFirstLidar = false;
		return LIDAR_MEAS;
	}

	// TODO checkstaticmotion
	if (true)
	{
		_framecount++;
		stateAugmentation(runEpoch);
		FrameBuffer.push_back(imu_state_id);

		auto iter_lidar = --lidar_states.end();
		_lidarframe.R_l_e = iter_lidar->second.orientation;
		_lidarframe.t_l_e = iter_lidar->second.position;
		_lidarframe.empty = false;


		if (use_corrdistort)
		{
			_lidarOdo.removeDistortion(_lastframe, _lidarframe);
		}
		if (use_scan)
		{
			if (!_lidarOdo.systemInited_)
			{
				_lidarOdo.systemInited_ = true;
			}
			else
			{
				Eigen::MatrixXd H;
				Eigen::VectorXd r;
				if (FrameBuffer.size() > 2)
				{
					_lidarOdo.process(_lastframe, _lidarframe);
					lidarMeasurementJacobian(_lidarOdo.lidarOdoObs, H, r, false, 1);
					cout << "odo update:" << endl;
					lidar_meas_update(H, r);
					//getchar();
				}
			}
		}

		if (use_map)
		{
			if (!_lidarmap.systemInited_)
			{
				_lidarmap.addPointcloudToMap(_lastframe);
				_lidarmap.last_R_l_e = _lidarframe.R_l_e;
				_lidarmap.last_t_l_e = _lidarframe.t_l_e;
			}
			else if (FrameBuffer.size() > 2&&_framecount%map_gap == 0)
			{
				R_w_e = _lidarmap.first_R_l_e;
				t_w_e = _lidarmap.first_t_l_e;
				Eigen::MatrixXd H;
				Eigen::VectorXd r;
				_lidarmap.process(_lidarframe, proj_lastmap);
				if(proj_lastmap)
					lidarMeasurementJacobian_test(_lidarmap.lidarMapObs, H, r, false, 1);
				else
					lidarMeasurementJacobian(_lidarmap.lidarMapObs, H, r, false, 1);

				std::cout << "map update:" << endl;
				lidar_meas_update(H, r);

				_lidarframe.R_l_e = (--lidar_states.end())->second.orientation;
				_lidarframe.t_l_e = (--lidar_states.end())->second.position;
				_lastframe.R_l_e = (++lidar_states.begin())->second.orientation;
				_lastframe.t_l_e = (++lidar_states.begin())->second.position;

				//_lidarmap.first_R_l_e = lidar_states.begin()->second.orientation;
				//_lidarmap.first_t_l_e = lidar_states.begin()->second.position;
				Eigen::Matrix3d new_firstR(lidar_states.begin()->second.orientation);
				Eigen::Vector3d new_firstt= lidar_states.begin()->second.position;
				_lidarmap.renewWframe(new_firstR, new_firstt);

				_lidarmap.last_R_l_e = _lidarframe.R_l_e;
				_lidarmap.last_t_l_e = _lidarframe.t_l_e;

				_lidarmap.addPointcloudToMap(_lastframe);

			/*	Eigen::Vector3d corr_diff =_lidarframe.t_l_e - _lastframe.t_l_e;
				cout << "corr_diff:" << corr_diff.transpose() << endl;
				cout << "corr_diff:" << corr_diff.norm() << endl;
				if (corr_diff.norm() < 1&&pos_jump==false)
				{
					_lidarmap.addPointcloudToMap(_lastframe);
				}*/
			}
			else
			{
				/*Eigen::Vector3d corr_diff = _lidarframe.t_l_e - _lastframe.t_l_e;
				cout << "corr_diff:" << corr_diff.transpose() << endl;
				cout << "corr_diff:" << corr_diff.norm() << endl;
				if (corr_diff.norm() < 1&& pos_jump==false)
				{
					_lidarmap.addPointcloudToMap(_lastframe);
				}*/

				_lidarmap.addPointcloudToMap(_lastframe);
				_lidarmap.last_R_l_e = _lidarframe.R_l_e;
				_lidarmap.last_t_l_e = _lidarframe.t_l_e;
			}
		}
		vector<StateIDType> rm_lidar_state_ids(0);
		if (FrameBuffer.size()>2)
		{
			//remove the last lidarframe
			rm_lidar_state_ids.push_back(*(++FrameBuffer.begin()));
		}
		
		//sort(rm_lidar_state_ids.begin(), rm_lidar_state_ids.end());

		//renew the pose stored in _lidarframe
		//auto iter_lidar = --lidar_states.end();
		_lidarframe.R_l_e = (--lidar_states.end())->second.orientation;
		_lidarframe.t_l_e = (--lidar_states.end())->second.position;

		removeLidarState(rm_lidar_state_ids);

		vector<StateIDType> vTemp;
		bool isRemove;
		for (int i = 0; i < FrameBuffer.size(); i++)
		{
			isRemove = false;

			for (int j = 0; j < rm_lidar_state_ids.size(); j++)
			{
				if (FrameBuffer[i] == rm_lidar_state_ids[j])
					isRemove = true;
			}

			if (!isRemove)
				vTemp.push_back(FrameBuffer[i]);
		}
		FrameBuffer = vTemp;

		_lastframe = _lidarframe;
	}
	return LIDAR_MEAS;
}

void gsins::t_gmultign_lidar::removeLidarState(vector<StateIDType>& rm_lidar_state_ids)
{
	if (rm_lidar_state_ids.size() <= 0)
		return;
	if (is_allTC)
	{
		for (const auto& lidar_id : rm_lidar_state_ids)
		{
			StateIDType _lidar_id = lidar_id;
			string string_id = state_id2str(_lidar_id);
			int idx_att_x = _param.getParam(_site, t_gpar::LIDAR_ATT_X, string_id);
			Matrix_remRC(_Qx, _param[idx_att_x].index, _param[idx_att_x].index);
			_param.delParam(idx_att_x);
			_param.reIndex();

			int idx_att_y = _param.getParam(_site, t_gpar::LIDAR_ATT_Y, string_id);
			Matrix_remRC(_Qx, _param[idx_att_y].index, _param[idx_att_y].index);
			_param.delParam(idx_att_y);
			_param.reIndex();

			int idx_att_z = _param.getParam(_site, t_gpar::LIDAR_ATT_Z, string_id);
			Matrix_remRC(_Qx, _param[idx_att_z].index, _param[idx_att_z].index);
			_param.delParam(idx_att_z);
			_param.reIndex();

			int idx_crd_x = _param.getParam(_site, t_gpar::LIDAR_CRD_X, string_id);
			Matrix_remRC(_Qx, _param[idx_crd_x].index, _param[idx_crd_x].index);
			_param.delParam(idx_crd_x);
			_param.reIndex();

			int idx_crd_y = _param.getParam(_site, t_gpar::LIDAR_CRD_Y, string_id);
			Matrix_remRC(_Qx, _param[idx_crd_y].index, _param[idx_crd_y].index);
			_param.delParam(idx_crd_y);
			_param.reIndex();

			int idx_crd_z = _param.getParam(_site, t_gpar::LIDAR_CRD_Z, string_id);
			Matrix_remRC(_Qx, _param[idx_crd_z].index, _param[idx_crd_z].index);
			_param.delParam(idx_crd_z);
			_param.reIndex();

			lidar_states.erase(lidar_id);

		}
	}
	else
	{
		for (const auto& lidar_id : rm_lidar_state_ids)
		{
			int lidar_sequence = std::distance(lidar_states.begin(),
				lidar_states.find(lidar_id));
			int lidar_state_start = 15 + 6 * lidar_sequence;
			int lidar_state_end = lidar_state_start + 6;
			Eigen::MatrixXd state_cov = get_statecov();
			if (lidar_state_end < state_cov.rows())
			{
				state_cov.block(lidar_state_start, 0,
					state_cov.rows() - lidar_state_end,
					state_cov.cols()) =
					state_cov.block(lidar_state_end, 0,
						state_cov.rows() - lidar_state_end,
						state_cov.cols());
			
			
				state_cov.block(0, lidar_state_start,
					state_cov.rows(),
					state_cov.cols() - lidar_state_end) =
					state_cov.block(0, lidar_state_end,
						state_cov.rows(),
						state_cov.cols() - lidar_state_end);
			
				state_cov.conservativeResize(
					state_cov.rows() - 6, state_cov.cols() - 6);
			}
			else
			{
				state_cov.conservativeResize(
					state_cov.rows() - 6, state_cov.cols() - 6);
			}
			lidar_states.erase(lidar_id);
			P_LL.conservativeResize(6 * lidar_states.size(), 6 * lidar_states.size());
			P_IL.conservativeResize(15, 6 * lidar_states.size());
			Eigen::MatrixXd state_cov_fixed = (state_cov +
				state_cov.transpose()) / 2.0;
			assign_statecov(state_cov_fixed);
		}
	}
}

int gsins::t_gmultign_lidar::_lidar_feedback(const ColumnVector & dx)
{
	if (lidar_states.size() > 0)
	{
		auto lidar_state_iter = lidar_states.begin();
		for (int i = 0; i < lidar_states.size();
			++i, ++lidar_state_iter)
		{
			StateIDType cur_id = lidar_state_iter->first;
			string _lidarid = state_id2str(cur_id);
		
			int idx = _param.getParam(_site, t_gpar::LIDAR_ATT_X, _lidarid);
			if (idx < 0) cerr << "lidar_feedback() idx wrong!" << endl;

			Eigen::Vector3d att_dx = Eigen::Vector3d(dx(idx + 1), dx(idx + 2), dx(idx + 3));
			Eigen::Vector3d crd_dx = Eigen::Vector3d(dx(idx + 4), dx(idx + 5), dx(idx + 6));

			t_gquat dq_lidar = t_gbase::rv2q(-att_dx);
			Eigen::Quaterniond _dq_lidar;
			_dq_lidar.w() = dq_lidar.q0; _dq_lidar.x() = dq_lidar.q1; _dq_lidar.y() = dq_lidar.q2; _dq_lidar.z() = dq_lidar.q3;
			_dq_lidar.normalize();
			lidar_state_iter->second.orientation = lidar_state_iter->second.orientation * _dq_lidar;
			lidar_state_iter->second.position -= crd_dx;
		}
	}
	return 1;
}

int gsins::t_gmultign_lidar::_lidar_feedback(const Eigen::VectorXd & dx)
{
	//cout << "lidar_feedaback====lidar_statesize:" << lidar_states.size() << endl;
	//cout<< "lidar_dx.segment:" << (dx.segment<6>(15)).transpose() << endl;
	int ins_size = 15;
	auto lidar_state_iter = lidar_states.begin();
	for (int i = 0; i < lidar_states.size();
		++i, ++lidar_state_iter)
	{
		const Eigen::VectorXd& delta_x_lidar = dx.segment<6>(ins_size + i * 6);
		t_gquat dq_lidar = t_gbase::rv2q(-delta_x_lidar.head<3>());
		Eigen::Quaterniond _dq_lidar;
		_dq_lidar.w() = dq_lidar.q0; _dq_lidar.x() = dq_lidar.q1; _dq_lidar.y() = dq_lidar.q2; _dq_lidar.z() = dq_lidar.q3;
		_dq_lidar.normalize();
		lidar_state_iter->second.orientation = lidar_state_iter->second.orientation * _dq_lidar;
		lidar_state_iter->second.position -= delta_x_lidar.tail<3>();
	}
	return 1;
}

string gsins::t_gmultign_lidar::state_id2str(const StateIDType& id)
{
	//get string camera state id
	ostringstream os;
	os << id;
	string lidar_id;
	istringstream is(os.str());
	is >> lidar_id;
	return lidar_id;
	
}



void gsins::t_gmultign_lidar::lidar_meas_update(const Eigen::MatrixXd& H, const Eigen::VectorXd& r)
{
	if (H.rows() == 0 || r.rows() == 0)
		return;
	if (is_allTC)
	{
		assert(H.cols() == 6 * lidar_states.size());
		int par_size = _param.parNumber();
		int ins_size = 15;
		int obs_size = H.rows();
		int N = lidar_states.size();
		Eigen::MatrixXd All_Variance = BaseMatrix2Eigen(_Qx);
		All_Variance.block(0, 0, ins_size, ins_size) = Pk;
		Eigen::MatrixXd All_H = Eigen::MatrixXd::Zero(H.rows(), par_size);

		assert(N == int(H.cols() / 6));
		assert(obs_size == r.rows());


		auto lidar_state_iter = lidar_states.begin();
		for (int i = 0; i < lidar_states.size();
			++i, ++lidar_state_iter)
		{
			string lidar_id = state_id2str(lidar_state_iter->first);
			int idx = _param.getParam(_site, t_gpar::LIDAR_ATT_X, lidar_id);
			All_H.block(0, idx, obs_size, 6) = H.block(0, 6 * i, obs_size, 6);

		}

		int ig_size = _param.parNumber() - 6 * lidar_states.size();
		Eigen::MatrixXd R = Eigen::MatrixXd::Identity(obs_size, obs_size)*observation_noise;

		ColumnVector dx;

		Eigen::MatrixXd H_thin;
		Eigen::VectorXd r_thin;
		if (All_H.rows() > All_H.cols())
		{

			Eigen::HouseholderQR<Eigen::MatrixXd> qr_helper(All_H);
			Eigen::MatrixXd Q = qr_helper.householderQ();
			Eigen::MatrixXd Q1;
			Q1 = Q.leftCols(par_size);
			H_thin = Q1.transpose() * All_H;
			r_thin = Q1.transpose() * r;
		}
		else
		{

			H_thin = All_H;
			r_thin = r;
		}

		const Eigen::MatrixXd& P = All_Variance;
		Eigen::MatrixXd S = H_thin * P*H_thin.transpose() +
			observation_noise * Eigen::MatrixXd::Identity(
				H_thin.rows(), H_thin.rows());

		Eigen::MatrixXd K_transpose = S.ldlt().solve(H_thin*P);
		Eigen::MatrixXd K = K_transpose.transpose();

		Eigen::VectorXd delta_x = K * r_thin;
		Eigen::MatrixXd I_KH = Eigen::MatrixXd::Identity(K.rows(), H_thin.cols()) - K * H_thin;
		All_Variance = I_KH * All_Variance*I_KH.transpose() + K * K.transpose()*observation_noise;
		dx = VectorXd2Columns(delta_x);
		_Qx = Eigen2BaseMatrix(All_Variance);


		int ign_size = _param.parNumber();
		Xk = Columns2VectorXd(dx).block(0, 0, ins_size, 1);
		Pk = BaseMatrix2Eigen(_Qx).block(0, 0, ins_size, ins_size);

		if (_gnss_feedback(dx) < 0) {
			cerr << "gns¡ªmeas-update-fail" << endl;
		}

		cout << "imu_xk:" << setprecision(8) << Xk.block(0, 0, 9, 1).transpose() << endl;

		sins.qeb = sins.qeb - Xk.block(0, 0, 3, 1);
		sins.ve = sins.ve - Xk.block(3, 0, 3, 1);
		sins.pos_ecef = sins.pos_ecef - Xk.block(6, 0, 3, 1);
		sins.eb = sins.eb + Xk.block(9, 0, 3, 1);
		sins.db = sins.db + Xk.block(12, 0, 3, 1);

		sins.qnb = t_gbase::m2qua(sins.eth.Cne) * sins.qeb;
		sins.vn = sins.eth.Cne*sins.ve;
		sins.pos = Cart2Geod(sins.pos_ecef, false);
		sins.att = t_gbase::q2att(sins.qnb);
		sins.Ceb = t_gbase::q2mat(sins.qeb);
		sins.Cnb = t_gbase::q2mat(sins.qnb);

		if (_lidar_feedback(dx) < 0) cerr << "lidar_update_error:" << endl;
	}
	else
	{
		int ins_size = 15;
		Eigen::MatrixXd All_H = Eigen::MatrixXd::Zero(H.rows(), H.cols() + ins_size);
		All_H.block(0, 15, H.rows(), H.cols()) = H;
		
		Eigen::MatrixXd H_thin;
		Eigen::VectorXd r_thin;
	
		if (All_H.rows() > All_H.cols())
		{
	
			Eigen::HouseholderQR<Eigen::MatrixXd> qr_helper(All_H);
			Eigen::MatrixXd Q = qr_helper.householderQ();
			Eigen::MatrixXd Q1 = Q.leftCols(15 + lidar_states.size() * 6);
	
			H_thin = Q1.transpose() * All_H;
			r_thin = Q1.transpose() * r;
		}
		else
		{
	
			H_thin = All_H;
			r_thin = r;
		}
		Eigen::MatrixXd state_cov = get_statecov();
		const Eigen::MatrixXd& P = state_cov;
		Eigen::MatrixXd S = H_thin * P*H_thin.transpose() +
			observation_noise * Eigen::MatrixXd::Identity(
				H_thin.rows(), H_thin.rows());
	
		Eigen::MatrixXd K_transpose = S.ldlt().solve(H_thin*P);
		Eigen::MatrixXd K = K_transpose.transpose();
	
		Eigen::VectorXd delta_x = K * r_thin;
	
		const Eigen::VectorXd& delta_x_imu = delta_x.head<15>();
	
		//cout << "res:" << r.transpose() << endl;

		sins.qeb = sins.qeb - delta_x_imu.head<3>();
		sins.ve = sins.ve - delta_x_imu.segment<3>(3);
		sins.pos_ecef = sins.pos_ecef - delta_x_imu.segment<3>(6);
		sins.eb = sins.eb + delta_x_imu.segment<3>(9);
		sins.db = sins.db + delta_x_imu.segment<3>(12);
	
		sins.qnb = t_gbase::m2qua(sins.eth.Cne) * sins.qeb;
		sins.vn = sins.eth.Cne*sins.ve;
		sins.pos = Cart2Geod(sins.pos_ecef, false);
		sins.att = t_gbase::q2att(sins.qnb);
		sins.Ceb = t_gbase::q2mat(sins.qeb);
		sins.Cnb = t_gbase::q2mat(sins.qnb);
		cout << "imu_xk:" <<setprecision(8)<< delta_x.block(0, 0, 15, 1).transpose() << endl;
	
		if(!_lidar_feedback(delta_x)) cerr<< "lidar_update_error:" << endl;
	
	
		Eigen::MatrixXd I_KH = Eigen::MatrixXd::Identity(K.rows(), H_thin.cols()) - K * H_thin;
	
		state_cov = I_KH * state_cov*I_KH.transpose() +
			K * K.transpose()*observation_noise;
		Eigen::MatrixXd state_cov_fixed = (state_cov + state_cov.transpose()) / 2.0;
	
		assign_statecov(state_cov_fixed);
	}
}


void gsins::t_gmultign_lidar::stateAugmentation(const t_gtime& runEpoch)
{
	double run_epoch = runEpoch.sow() + runEpoch.dsec();
	double dt = sins.t - run_epoch;

	//compensate due to time 
	Eigen::Vector3d BLH = sins.pos - sins.eth.v2dp(sins.vn, dt);//BLH
	Eigen::Vector3d XYZ = gsins::Geod2Cart(BLH, false);

	Eigen::Vector3d imu_vel = sins.vn - sins.an*dt;//ENU
	const Eigen::Matrix3d& R_i_l = R_lidar_imu.transpose();
	const Eigen::Vector3d& t_l_i = t_lidar_imu;
	const Eigen::Matrix3d& R_i_n = t_gbase::q2mat(sins.qnb);
	Eigen::Matrix3d R_n_i = R_i_n.transpose();
	Eigen::Matrix3d R_n_e = t_gbase::Cen(BLH);
	Eigen::Matrix3d R_e_i = R_n_i * R_n_e.transpose();
	Eigen::Vector3d t_l_e = XYZ + R_n_e * R_i_n*t_l_i;
	Eigen::Matrix3d R_e_l = R_i_l * R_e_i;
	//maintain the numerical stability of derivatives

	lidar_states[imu_state_id] = LIDARState(imu_state_id);
	LIDARState& lidar_state = lidar_states[imu_state_id];

	lidar_state.time = run_epoch;
	lidar_state.orientation = R_e_l.transpose();
	lidar_state.position = t_l_e;
	
	if (is_allTC)
	{
		//get string lidar state id
		string lidar_id = state_id2str(lidar_state.id);
		t_gpar att_x_par(_site, t_gpar::LIDAR_ATT_X, _param.parNumber() + 1, lidar_id);
		_param.addParam(att_x_par);
		t_gpar att_y_par(_site, t_gpar::LIDAR_ATT_Y, _param.parNumber() + 1, lidar_id);
		_param.addParam(att_y_par);
		t_gpar att_z_par(_site, t_gpar::LIDAR_ATT_Z, _param.parNumber() + 1, lidar_id);
		_param.addParam(att_z_par);

		t_gpar crd_x_par(_site, t_gpar::LIDAR_CRD_X, _param.parNumber() + 1, lidar_id);
		_param.addParam(crd_x_par);
		t_gpar crd_y_par(_site, t_gpar::LIDAR_CRD_Y, _param.parNumber() + 1, lidar_id);
		_param.addParam(crd_y_par);
		t_gpar crd_z_par(_site, t_gpar::LIDAR_CRD_Z, _param.parNumber() + 1, lidar_id);
		_param.addParam(crd_z_par);

		int ins_size = 15;
		Eigen::MatrixXd All_Variance = BaseMatrix2Eigen(_Qx);
		//time update
		size_t old_rows = All_Variance.rows();
		size_t old_cols = All_Variance.cols();
		All_Variance.block(0, 0, ins_size, ins_size) = Pk;
		All_Variance.block(ins_size, 0, old_rows - ins_size, ins_size) = All_Variance.block(ins_size, 0, old_rows - ins_size, ins_size)*Phik.transpose();
		All_Variance.block(0, ins_size, ins_size, old_rows - ins_size) = Phik * All_Variance.block(0, ins_size, ins_size, old_rows - ins_size);
		Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, old_cols);
		//lidar att to imu att
		J.block<3, 3>(0, 0) = -R_e_l;
		//lidar pos to imu att
		J.block<3, 3>(3, 0) = t_gbase::askew(R_n_e*R_i_n*t_l_i);
		//lidar pos to imu pos
		J.block<3, 3>(3, 6) = Eigen::Matrix3d::Identity();
		All_Variance.conservativeResize(old_rows + 6, old_cols + 6);
		const Eigen::MatrixXd P11 = All_Variance.block(0, 0, old_rows, old_cols);
		//matrix augment
		All_Variance.block(old_rows, 0, 6, old_cols) = J * P11;
		All_Variance.block(0, old_cols, old_rows, 6) = All_Variance.block(old_rows, 0, 6, old_cols).transpose();
		All_Variance.block<6, 6>(old_rows, old_cols) = J * P11 * J.transpose();
		//Eigen::EigenSolver<Eigen::MatrixXd> es(All_Variance);
		//cout << "The eigenvalues of A are:" << endl<<setprecision(12) << es.eigenvalues() << endl;
		//keep positive definite of matrix
		All_Variance.block<6, 6>(old_rows, old_cols) += Eigen::MatrixXd::Identity(6, 6)*1e-10;
		_Qx = Eigen2BaseMatrix(All_Variance);
		assert(_Qx.rows() == _param.parNumber());
	}
	else
	{
		Eigen::Matrix<double, 6, 15> J = Eigen::Matrix<double, 6, 15>::Zero();
		//lidar att to imu att
		J.block<3, 3>(0, 0) = -R_e_l;
		//lidar pos to imu att
		J.block<3, 3>(3, 0) = t_gbase::askew(R_n_e*R_i_n*t_l_i);
		//lidar pos to imu pos
		J.block<3, 3>(3, 6) = Eigen::Matrix3d::Identity();
		
		//covariance timeupdate
		if (lidar_states.size() > 0)
		{
			P_IL = Phik * P_IL;
		}

		Eigen::MatrixXd state_cov;
		state_cov = get_statecov();
		size_t old_rows = state_cov.rows();
		size_t old_cols = state_cov.cols();
		state_cov.conservativeResize(old_rows + 6, old_cols + 6);
	
		const Eigen::Matrix<double, 15, 15>& P11 =
			state_cov.block<15, 15>(0, 0);
		const Eigen::MatrixXd& P12 =
			state_cov.block(0, 15, 15, old_cols - 15);
	
		state_cov.block(old_rows, 0, 6, old_cols) << J * P11, J*P12;
		state_cov.block(0, old_cols, old_rows, 6) =
			state_cov.block(old_rows, 0, 6, old_cols).transpose();
		state_cov.block<6, 6>(old_rows, old_cols) =
			J * P11 * J.transpose();
	
		//Eigen::EigenSolver<Eigen::MatrixXd> es(state_cov);
		//cout << "The eigenvalues of A are:" << endl <<setprecision(12)<< es.eigenvalues() << endl;
		//state_cov.block<6, 6>(old_rows, old_cols) += state_cov.block<6, 6>(old_rows, old_cols) + Eigen::MatrixXd::Identity(6,6)*1e-10;
		
		P_LL.conservativeResize(P_LL.rows() + 6, P_LL.cols() + 6);
		P_IL.conservativeResize(P_IL.rows(), P_IL.cols() + 6);
		Eigen::MatrixXd state_cov_fixed = (state_cov +
			state_cov.transpose()) / 2.0;
	
		assign_statecov(state_cov_fixed);
		
	}

	int ins_size = 15;
	Phik = Eigen::MatrixXd::Identity(ins_size, ins_size);

}



//void gsins::t_gmultign_lidar::pruneLidarStateBuffer()
//{
//	lidar_meas_update(H_x, r);
//
//	if (is_allTC)
//	{
//		for (const auto& lidar_id : rm_lidar_state_ids)
//		{
//			StateIDType _lidar_id = lidar_id;
//			string string_id = state_id2str(_lidar_id);
//			int idx_att_x = _param.getParam(_site, t_gpar::LIDAR_ATT_X, string_id);
//			Matrix_remRC(_Qx, _param[idx_att_x].index, _param[idx_att_x].index);
//			_param.delParam(idx_att_x);
//			_param.reIndex();
//
//			int idx_att_y = _param.getParam(_site, t_gpar::LIDAR_ATT_Y, string_id);
//			Matrix_remRC(_Qx, _param[idx_att_y].index, _param[idx_att_y].index);
//			_param.delParam(idx_att_y);
//			_param.reIndex();
//
//			int idx_att_z = _param.getParam(_site, t_gpar::LIDAR_ATT_Z, string_id);
//			Matrix_remRC(_Qx, _param[idx_att_z].index, _param[idx_att_z].index);
//			_param.delParam(idx_att_z);
//			_param.reIndex();
//
//			int idx_crd_x = _param.getParam(_site, t_gpar::LIDAR_CRD_X, string_id);
//			Matrix_remRC(_Qx, _param[idx_crd_x].index, _param[idx_crd_x].index);
//			_param.delParam(idx_crd_x);
//			_param.reIndex();
//
//			int idx_crd_y = _param.getParam(_site, t_gpar::LIDAR_CRD_Y, string_id);
//			Matrix_remRC(_Qx, _param[idx_crd_y].index, _param[idx_crd_y].index);
//			_param.delParam(idx_crd_y);
//			_param.reIndex();
//
//			int idx_crd_z = _param.getParam(_site, t_gpar::LIDAR_CRD_Z, string_id);
//			Matrix_remRC(_Qx, _param[idx_crd_z].index, _param[idx_crd_z].index);
//			_param.delParam(idx_crd_z);
//			_param.reIndex();
//
//			lidar_states.erase(lidar_id);
//
//		}
//	}
//	else
//	{
//		for (const auto& lidar_id : rm_lidar_state_ids)
//		{
//			int lidar_sequence = std::distance(lidar_states.begin(),
//				lidar_states.find(lidar_id));
//			int lidar_state_start = 15 + 6 * lidar_sequence;
//			int lidar_state_end = lidar_state_start + 6;
//			Eigen::MatrixXd state_cov = get_statecov();
//			if (lidar_state_end < state_cov.rows())
//			{
//				state_cov.block(lidar_state_start, 0,
//					state_cov.rows() - lidar_state_end,
//					state_cov.cols()) =
//					state_cov.block(lidar_state_end, 0,
//						state_cov.rows() - lidar_state_end,
//						state_cov.cols());
//			
//			
//				state_cov.block(0, lidar_state_start,
//					state_cov.rows(),
//					state_cov.cols() - lidar_state_end) =
//					state_cov.block(0, lidar_state_end,
//						state_cov.rows(),
//						state_cov.cols() - lidar_state_end);
//			
//				state_cov.conservativeResize(
//					state_cov.rows() - 6, state_cov.cols() - 6);
//			}
//			else
//			{
//				state_cov.conservativeResize(
//					state_cov.rows() - 6, state_cov.cols() - 6);
//			}
//			lidar_states.erase(lidar_id);
//			P_LL.conservativeResize(6 * lidar_states.size(), 6 * lidar_states.size());
//			P_IL.conservativeResize(15, 6 * lidar_states.size());
//			Eigen::MatrixXd state_cov_fixed = (state_cov +
//				state_cov.transpose()) / 2.0;
//			assign_statecov(state_cov_fixed);
//		}
//	}
//	return;
//}

void gsins::t_gmultign_lidar::assign_statecov(const Eigen::MatrixXd &P)
{
	size_t rows = P.rows();
	size_t cols = P.cols();
	assert(rows == cols);
	int N = lidar_states.size();
	Pk = P.block(0, 0, 15, 15);
	P_LL = P.block(15, 15, 6 * N, 6 * N);
	P_IL = P.block(0, 15, 15, 6 * N);
}

Eigen::MatrixXd gsins::t_gmultign_lidar::get_statecov()
{
	size_t imu_rows = Pk.rows();
	size_t lidar_rows = P_LL.rows();
	Eigen::MatrixXd state_cov = Eigen::MatrixXd::Zero(imu_rows + lidar_rows, imu_rows + lidar_rows);
	state_cov.block(0, 0, imu_rows, imu_rows) = Pk;
	state_cov.block(imu_rows, imu_rows, lidar_rows, lidar_rows) = P_LL;
	state_cov.block(0, imu_rows, imu_rows, lidar_rows) = P_IL;
	state_cov.block(imu_rows, 0, lidar_rows, imu_rows) = P_IL.transpose();
	return state_cov;
}

int gsins::t_gmultign_lidar::_mult_meas_update(double fading)
{
	gtrace("t_gintegration::_meas_update");

	//Eigen::VectorXd Pxz = Eigen::VectorXd::Zero(nq),
	//	Kk = Eigen::VectorXd::Zero(nq), Hi = Eigen::VectorXd::Zero(nq);
	

	if (!_set_meas()) {
		if (_log) _log->comment(1, "t_gintegration", _ins_crt.str_ymdhms(" _set_meas failed: "));
		_reset();
		Flag = NO_MEAS;
		return -1;
	}

	int new_nq = nq + 6 * lidar_states.size();
	Eigen::VectorXd Pxz = Eigen::VectorXd::Zero(new_nq),
		Kk = Eigen::VectorXd::Zero(new_nq), Hi = Eigen::VectorXd::Zero(new_nq);
	Eigen::MatrixXd new_Pk = get_statecov();
	assert(new_Pk.rows() == new_nq);
	Eigen::MatrixXd new_Hk = Eigen::MatrixXd::Zero(Hk.rows(), new_nq);
	new_Hk.block(0, 0, Hk.rows(), Hk.cols()) = Hk;
	Eigen::VectorXd new_Xk = Eigen::VectorXd::Zero(new_nq);

	// before estimating 
	if (outlier_detect()) { _reset(); Flag = NO_MEAS; return -1; }

	for (int i = 0; i < nr; i++)
	{
		Hi = new_Hk.block(i, 0, 1, new_nq).transpose();
		Pxz = new_Pk * Hi;
		double Pz0 = (Hi.transpose()*Pxz), r = Zk(i) - (Hi.transpose()*new_Xk);
		double Pzz = Pz0 + Rk(i, i);
		Kk = Pxz * (1.0 / Pzz);
		new_Xk += Kk * r;
		new_Pk = new_Pk - Kk * Pxz.transpose();
	}
	t_gbase::symmetry(new_Pk);
	Pk = (new_Pk).block(0, 0, Pk.rows(), Pk.cols());
	Xk = new_Xk.block(0, 0, 15, 1);

	if (fading > 1.0) Pk = Pk * fading;
	t_gbase::symmetry(Pk);
	Phik = Eigen::MatrixXd::Identity(15, 15);
	set_meas_flag(0);

	// after estimating
	if (outlier_detect()) { _reset(); Flag = NO_MEAS; return -1; }

	if (Flag != NO_MEAS)
	{
		_lidar_feedback(new_Xk);
		assign_statecov(new_Pk);
	}

	return 1;
}


#endif // USE_VIS