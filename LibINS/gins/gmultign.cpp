#ifdef USE_VIS
#include "gmultign.h"
#include "gproc/gpreproc.h"
#include "gutils/gtimesync.h"
#include <iostream>
using namespace gsins;
using namespace std;
using namespace gvis;

///******************************************* t_gmultign ***********************************************/
t_gmultign::t_gmultign(string site, string site_base, t_gsetbase* gset, t_gallproc *allproc) :
	t_gintegration(site, site_base, gset, allproc),
	t_gsinskf(gset),
	t_gpvtflt(site, site_base, gset, allproc),
	t_gvis(gset),
	t_gspp(site, gset)
{

	_imgdata = nullptr;
}

void gsins::t_gmultign::setgallparuse(bool status)
{
	is_allTC = status;
}

void gsins::t_gmultign::Add_IMG(t_gimgdata * img)
{
	_imgdata = img;
}


bool gsins::t_gmultign::_img_time_valid(t_gtime imgt, t_gtime inst)
{
	gtrace("t_gintegration::_time_valid");

	bool res_valid = false;
	//cout << "diff:" << abs(inst.diff(imgt)) << endl;
	if (abs(inst.diff(imgt)) < _shm.delay) res_valid = true;
	else res_valid = false;

	return res_valid;
}


int gsins::t_gmultign::processBatch(const t_gtime & beg, const t_gtime & end, bool prtOut)
{
	gtrace("t_gintegration::processBatch");
#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif 
	_gmutex.lock();

	if (_grec == nullptr) {
		ostringstream os;
		os << "ERROR: No object found (" << _site << "). Processing terminated!!! " << beg.str_ymdhms() << " -> " << end.str_ymdhms() << endl;
		if (_log) _log->comment(0, "t_gintegration", os.str());
		 _gmutex.unlock(); 
		return -1;
	}

	// get the time
	_beg_time = _ins_beg; _end_time = _ins_end;
	_beg_time = (_beg_time < beg) ? beg : _beg_time;
	_end_time = (_end_time > end) ? end : _end_time;

	if (_isBase) {
		if (_setCrd() < 0) {
			if (_log) _log->comment(0, "t_gintegration", "Base station have no coordinate");
			_gmutex.unlock(); 
			return -1;
		}
	}

	_dop.set_data(_gnav, _gobs, _site);
	_running = true;
	int sign = 1;

	if (!this->_init()) {
		if (_log) _log->comment(0, "t_gintegration", "init failed");
		 _gmutex.unlock(); 
		return -1;
	}

	t_gtime tmp = t_gtime::current_time(t_gtime::GPS);
	if (!_beg_end) {
		tmp = _ins_end; _ins_end = _ins_beg; _ins_beg = tmp;
		tmp = _gnss_end; _gnss_end = _gnss_beg; _gnss_beg = tmp;
		_ins_crt = _ins_beg;
		sign = -1;
		if (_log) _log->comment(2, "t_gintegration", "Filtering in end -> begin direction!");
	}
	else {
		sign = 1;
		if (_log) _log->comment(2, "t_gintegration", "Filtering in begin -> end direction!");
	}

	// need execute frequently in real-time mode
	_gobs->setepoches(_site);
	bool is_GNSS = ((_ign_type != IGN_TYPE::VIO) && (_ign_type != IGN_TYPE::PURE_INS));
	bool is_VISION = ((_ign_type == IGN_TYPE::VIO) || (_ign_type == IGN_TYPE::MULTIGN_LC)
		|| (_ign_type == IGN_TYPE::MULTIGN_TC) || (_ign_type == IGN_TYPE::MULTIGN_STC));
	// time align
	if (_ins_beg < _gnss_beg && is_GNSS) {
		_ins_crt = _imudata->erase_bef(_gnss_beg);
	}
	else if (((_ign_type != IGN_TYPE::VIO) && (_ign_type != IGN_TYPE::PURE_INS)))
	{
		_gmutex.unlock();
		t_gpvtflt::processBatch(_gnss_crt, _ins_crt, false);
		_gmutex.lock();

		// next gnss time 
		while (_gnss_crt < _ins_crt)
		{
			if (_sampling > 1) _gnss_crt.add_secs(int(sign * _sampling));  // =<1Hz data
			else               _gnss_crt.add_dsec(sign * _sampling);       //  >1Hz data
		}
	}

	std::cerr << "\n" << _site << ": Start GNSS/INS Integration Processing: " << _ins_beg.str_ymdhms() << " " << _ins_end.str_ymdhms() << endl;

	_extrinsic_init();

	t_gposdata::data_pos posdata; posdata = 0.0;
	//to do remove this one
	bool time_loop = true;
	// this bool is to eliminate the variance if mode is VIO
	if (_ign_type == IGN_TYPE::MULTIGN_TC) this->setgallparuse(true);
	
	//bool is_VISION = false;
	while (time_loop)
	{
		if (_ins_crt > _ins_end || _gnss_crt > _gnss_end)
		{
			std::cout << " Processing Finish!" << endl;
			if (_log) _log->comment(0, "t_gintegration", "Processing Finish!");
			time_loop = false; break;
		}
		_imudata->load(_wm, _vm, _shm.t, _shm.ts, _shm.nSamples);
		_ins_crt = t_gtime(_gnss_crt.gwk(), _shm.t);
		if (!_aligned)
		{
			if (_time_valid(_gnss_crt, _ins_crt))Flag = _getPOS(posdata);
			_aligned = cascaded_align(posdata.pos, posdata.vn);
		}
		else
		{
			if ((_ign_type == IGN_TYPE::TCI || _ign_type == IGN_TYPE::MULTIGN_TC) && !_initial_merge)
				_merge_init();

			sins.Update(_wm, _vm, _shm);
			if (_ign_type == PURE_INS) {
				t_gsinskf::write();
				continue;
			}
			time_update(sins.nts, 0.0);

			// all mode that need image data
			if (_ign_type == IGN_TYPE::VIO || _ign_type == IGN_TYPE::MULTIGN_LC || 
				_ign_type==IGN_TYPE::MULTIGN_TC || _ign_type==IGN_TYPE::MULTIGN_STC)
			{
				_imgdata->load(_shm.t, _img_t, _img_path);
				_img_crt = t_gtime(_gnss_crt.gwk(), _img_t);
				imgproc->load_imuobs(_shm.t, _wm, _vm);
			}
			// todo need more judge
			bool flag_gnss = (is_GNSS && (_time_valid(_gnss_crt, _ins_crt)));
			bool flag_img = (is_VISION && (_img_time_valid(_img_crt, _ins_crt)));
			
			if (flag_img || flag_gnss)
			{
				cout << endl;
				cout << "imu_time:" << setprecision(10) << _shm.t << endl;
				cout << "img_now:" << setprecision(10) << _img_t << "  " << flag_img << endl;
				cout << "gns_now:" << setprecision(10) << (_gnss_crt.sow() + _gnss_crt.dsec()) << "  " << flag_gnss << endl;
			}
		
			if (flag_img)
			{
				// todo  merge to one function
				if (_ign_type == IGN_TYPE::VIO)
				{
				    imgproc->load_imgobs(_img_t, _img_path);
					Flag = MEAS_TYPE(this->_processCamEpoch(_img_crt));
				}
				else if (_ign_type == IGN_TYPE::MULTIGN_LC)
				{
					imgproc->load_imgobs(_img_t, _img_path);
					Flag = MEAS_TYPE(this->_processCamEpoch(_img_crt));
				}
				else if (_ign_type == IGN_TYPE::MULTIGN_TC)
				{
					imgproc->load_imgobs(_img_t, _img_path);
					Flag = MEAS_TYPE(this->_processCamEpoch(_img_crt));
				}
				else if (_ign_type == IGN_TYPE::MULTIGN_STC)
				{
					imgproc->load_imgobs(_img_t, _img_path);
					Flag = MEAS_TYPE(this->_processCamEpoch(_img_crt));
				}
				else
				{
					cerr << "flag_img=true is a wrong mode!" << endl;
				}
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
					// if (!_isBase) {
						_timeUpdate(_gnss_crt);
						Flag = MEAS_TYPE(this->_processEpoch(_gnss_crt));
						if (Flag == NO_MEAS) { _reset(); }
					// }
					// else
					// { // bugs
					// 	Flag = _getFixedDDObs();
					// 	if (Flag != NO_MEAS) {
					// 		_meas_update();
					// 	}
					// }
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
				else if (_ign_type == IGN_TYPE::MULTIGN_LC)
				{
					Flag = _getPOS(posdata);
					if (Flag != NO_MEAS) _mult_meas_update();
				}
				else if (_ign_type == IGN_TYPE::MULTIGN_TC)
				{
					_timeUpdate(_gnss_crt);
					Flag = MEAS_TYPE(this->_processGNSSEpoch(_gnss_crt));
					if (Flag == NO_MEAS) { _reset(); }
				}
				else if (_ign_type == IGN_TYPE::MULTIGN_STC)
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

			//if (flag_gnss)
			if (flag_gnss || flag_img)
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
					" " << Ve(0) << " "
					" " << Ve(1) << " "
					" " << Ve(2) << " " << endl;
			}


			//use for glfw 
			imu_state.orientation = Eigen::Quaterniond(t_gbase::q2mat(sins.qnb));
			imu_state.position = Geod2Cart(sins.pos, false);

			if(_ins_crt <= _ins_end)
			publisher.UpdateNewState(imu_state, cam_states, initCamPos);
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

int gsins::t_gmultign::processBatch(const t_gtime& beg, const t_gtime& end)
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


	// this bool is to eliminate the variance if mode is VIO
	if (_ign_type == IGN_TYPE::MULTIGN_TC) this->setgallparuse(true);
	is_GNSS = ((_ign_type != IGN_TYPE::VIO) && (_ign_type != IGN_TYPE::PURE_INS));
	is_VISION = ((_ign_type == IGN_TYPE::VIO) || (_ign_type == IGN_TYPE::MULTIGN_LC)
		|| (_ign_type == IGN_TYPE::MULTIGN_TC) || (_ign_type == IGN_TYPE::MULTIGN_STC));


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
	cout << "EPOCH " << _ins_crt.str_ymdhms() << endl;
	if (!_extrinsic_init()) {
		if (_log) _log->comment(0, "t_gintegration", "init extrinsic parameters failed");
		return -1;
	}
	t_gposdata::data_pos posdata;
	bool time_loop = true;
	// ins time loop
	while (time_loop)
	{
		if (_ins_crt > _ins_end || _gnss_crt > _gnss_end)
		{
			cerr << " Processing Finished!" << endl;
			if (_log) _log->comment(0, "t_gintegration", "Processing Finished!");
			time_loop = false; break;
		}
		if (realtime) {
			t_gtime crt = t_gtime::current_time(t_gtime::GPS); this->_avaliable(_ins_crt);
		}
		//if (!_imudata->load(_wm, _vm, _shm.t, _shm.ts, _shm.nSamples, true)) { cerr << "\n Processing Finished!" << endl; break; }
		if (!_imudata->load(_wm, _vm, _shm.t, _shm.ts, _shm.nSamples)) { cerr << "\n Processing Finished!" << endl; break; }
		_ins_crt = t_gtime(_gnss_crt.gwk(), _shm.t);
		if (!_aligned)
		{
			if (_time_valid(_gnss_crt, _ins_crt)) Flag = _getPOS(posdata);
			_aligned = cascaded_align(posdata.pos, posdata.vn);
		}
		else
		{
			if ((_ign_type == IGN_TYPE::TCI || _ign_type == IGN_TYPE::MULTIGN_TC) && !_initial_merge)
				_merge_init();

			sins.Update(_wm, _vm, _shm);
			if (_ign_type == PURE_INS) {
				cerr << "\r" << _ins_crt.str_ymdhms();
				t_gsinskf::write_tmp(); continue;
			}

			time_update(sins.nts, 0.0);

			//fill imu data in imgproc
			//all mode that need image data
			if (_ign_type == IGN_TYPE::VIO || _ign_type == IGN_TYPE::MULTIGN_LC ||
				_ign_type == IGN_TYPE::MULTIGN_TC || _ign_type == IGN_TYPE::MULTIGN_STC)
			{
				_imgdata->load(_shm.t, _img_t, _img_path);
				_img_crt = t_gtime(_gnss_crt.gwk(), _img_t);
				imgproc->load_imuobs(_shm.t, _wm, _vm);
			}

			int imeas = _getMeas(), irc = -1;
			set<MEAS_TYPE>::const_iterator it = _Meas_Type.begin();
			//cout sensor time
			while (it != _Meas_Type.end())
			{
				cout << "imu_time:" << setprecision(10) << _shm.t << endl;
				switch (*it)
				{
				case GNSS_MEAS:
				{
					cout << "gns_now:" << setprecision(10) << (_gnss_crt.sow() + _gnss_crt.dsec()) << endl;
				}break;
				case VIS_MEAS:
				{
					cout << "img_now:" << setprecision(10) << _img_t << endl;
				}break;
				default:break;
				}it++;
				cout << endl;
			}
			it = _Meas_Type.begin();
			while (it != _Meas_Type.end())
			{
				switch (*it)
				{
				case GNSS_MEAS:
				{
					if (irc = _GNSS_Update() > 0)t_gsinskf::_restore(Xk, Pk);
					else t_gsinskf::_recover(Xk, Pk);

					//new_add
					if (Flag != NO_MEAS) {
						feedback(glv.INF);
						sins.pure_ins_time = 0.0;
					}
				}break;
				case VIS_MEAS:
				{
					if (irc = _Cam_Update() > 0) t_gsinskf::_restore(Xk, Pk);
					else t_gsinskf::_recover(Xk, Pk);
					//no need feedback
				}break;
				case ZUPT_MEAS:
				{
					if (irc = _ZUPT_Update() > 0) t_gsinskf::_restore(Xk, Pk);
					else t_gsinskf::_recover(Xk, Pk);

					//new_add
					if (Flag != NO_MEAS) {
						feedback(glv.INF);
						sins.pure_ins_time = 0.0;
					}
				}break;
				case NHC_MEAS:
				{
					if (irc = _NHC_Update() > 0) t_gsinskf::_restore(Xk, Pk);
					else t_gsinskf::_recover(Xk, Pk);

					//new_add
					if (Flag != NO_MEAS) {
						feedback(glv.INF);
						sins.pure_ins_time = 0.0;
					}
				} break;
				default:break;
				}it++;
			}

			if (Flag != NO_MEAS) {
				//feedback(glv.INF);
				//sins.pure_ins_time = 0.0;
				_prt_port();
				cerr << "\r" << _ins_crt.str_ymdhms();
				if (_log) _log->comment(0, "t_gintegration",
					_ins_crt.str_ymdhms(" epoch ") + "Processing GNSS/INS Integration Success!");
			}

			if (!(_isBase && _ign_type == IGN_TYPE::TCI)) {
				t_gsinskf::write_tmp();
			}
			//use for glfw 
			imu_state.orientation = Eigen::Quaterniond(t_gbase::q2mat(sins.qnb));
			imu_state.position = Geod2Cart(sins.pos, false);

			if (_ins_crt <= _ins_end)
				publisher.UpdateNewState(imu_state, cam_states, initCamPos);
		}
		while (_ins_crt.diff(_gnss_crt) > _sampling)
		{
			if (_sampling > 1) _gnss_crt.add_secs(int(_sampling));  // =<1Hz data
			else               _gnss_crt.add_dsec(_sampling);       //  >1Hz data
		}

	}

	//	_gmutex.unlock();
	return 0;
}

bool gsins::t_gmultign::_extrinsic_init()
{
	if (!imgproc->estimate_extrinsic || (_ign_type != IGN_TYPE::MULTIGN_TC && _ign_type != IGN_TYPE::VIO)) return true;
	if (_ign_type == IGN_TYPE::MULTIGN_TC)
	{
		try
		{
			t_gallpar param_extended;
			SymmetricMatrix Qx_extended;

			for (int i = 0; i < _param.parNumber(); i++)
			{
				param_extended.addParam(_param[i]);
			}
			
			// attitude	    
			Eigen::Vector3d att = imgproc->T_cam0_imu.linear().eulerAngles(2, 1, 0);
			for (int ipar = t_gpar::EXTRINSIC_ATT_X, i = 0; ipar <= t_gpar::EXTRINSIC_ATT_Z; ipar++, i++)
			{
				t_gpar attr(_site, t_gpar::t_type(ipar), ipar, "");
				param_extended.addParam(attr);
			}
			// position 
			Eigen::Vector3d trans = imgproc->t_cam0_imu;
			for (int ipar = t_gpar::EXTRINSIC_CRD_X, i = 0; ipar <= t_gpar::EXTRINSIC_CRD_Z; ipar++, i++)
			{
				t_gpar pos(_site, t_gpar::t_type(ipar), ipar, "");
				param_extended.addParam(pos);
			}
			int ibefore = param_extended.getParam(_site, t_gpar::EXTRINSIC_ATT_X, "");
			int iend = param_extended.parNumber();

			Qx_extended.ReSize(param_extended.parNumber()); Qx_extended = 0.0;


			for (int i = ibefore, j = 0; j < 3; i++, j++)
			{
				Qx_extended(i + 1, i + 1) = imgproc->initial_extrinsic_rotation_cov(j);
			}
			for (int i = ibefore + 3, j = 0; j < 3; i++, j++)
			{
				Qx_extended(i + 1, i + 1) = imgproc->initial_extrinsic_translation_cov(j);
			}

			_param = param_extended;
			_Qx = Qx_extended;

			_param.reIndex();
		}
		catch (...)
		{
			return false;
		}
	}
	else if (_ign_type == IGN_TYPE::VIO)
	{
		try
		{
			CAMState cam(0);
			cam.orientation = Eigen::Quaterniond(imgproc->T_cam0_imu.linear());
			cam.position = imgproc->t_cam0_imu;
			t_gvis::extrinsic = cam;

			const int r = P_CC.rows(), c = P_CC.cols();
			assert(r != 0 || c != 0);
			Eigen::MatrixXd P_CC_new = Eigen::MatrixXd::Zero(6, 6);
			for (int i = 0; i < 3; i++)
			{
				P_CC_new(i, i) = imgproc->initial_extrinsic_rotation_cov(i);
			}
			for (int i = 3, j = 0; j < 3; i++, j++)
			{
				P_CC_new(i, i) = imgproc->initial_extrinsic_translation_cov(j);
			}
			P_CC = P_CC_new;
			Eigen::Matrix<double, 15, 6> P_IC_new = Eigen::Matrix<double, 15, 6>::Zero();
			P_IC = P_IC_new;
		}
		catch (...)
		{
			return false;
		}
	}
	return true;
}

int gsins::t_gmultign::_processGNSSEpoch(const t_gtime & runEpoch)
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

	//	double tmp_gnss = runEpoch.sow() + runEpoch.dsec();
	//	//if (tmp_gnss >= 185956.3 && tmp_gnss <=717 185976.3)
	//	if(tmp_gnss>=442961.6 && tmp_gnss<=443061.6)	
	//{
	//	while (_data.size() > 5)
	//		_data.pop_back();
	//
	//	cout << _data.size() << endl;
	//}
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
		//这里可能有问题
		//if (_valid_ins_constraint())
		//{
		//	Eigen::Vector3d XYZ_INS = sins.pos_ecef + sins.Ceb*lever;
		//	t_gtriple XYZ(XYZ_INS(0), XYZ_INS(1), XYZ_INS(2));
		//	_external_pos(XYZ, t_gtriple());
		//}
		Eigen::Vector3d XYZ_INS = sins.pos_ecef + sins.Ceb*lever;
		t_gtriple XYZ(XYZ_INS(0), XYZ_INS(1), XYZ_INS(2));
		_external_pos(XYZ, t_gtriple());

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
			dynamic_cast<t_gcombDD*>(&(*_base_model))->set_rec_info(_gallobj->obj(_site_base)->crd(_epoch), _vBanc(4), _vBanc_base(4));
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
			if (_log) _log->comment(2, "gintegration", "merge failed!");
			return NO_MEAS;
		}

		try
		{
			////从update函数中拿出来
			//Matrix K;            // Kalman gain
			//Matrix NN;
			//SymmetricMatrix Pli = P.i();

			//NN = Pli + A * _Qx * A.t();
			//K = _Qx * A.t() * NN.i();
			////cout << dx << endl;
			////  cout  << "Qx" << endl << Qx << endl << "KAQx" << endl << K*A*Qx << endl << "K" << K << endl;
			//IdentityMatrix I = IdentityMatrix(_Qx.nrows());
			//Matrix KA = K * A;
			//Matrix I_KA = I - K * A;

			//// ZhouYuxuan modified. "K * ( l - A * dx )" => "K * l"
			//dx = K * l;               // update state vector
			//_Qx << I_KA * _Qx * I_KA.t() + K * Pli *K.t();
			_filter->update(A, P, l, dx, _Qx);

		}
		catch (...)
		{
			if (_log) _log->comment(2, "gintegration", "meas_update failed!");
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
	cout << "GNSS_Meas_UP:" << endl;
#ifdef USE_VIS
	if (_ign_type == IGN_TYPE::MULTIGN_TC)
	{
		if (_cam_feedback(dx) < 0)
		{
			cerr << "GNSS MeasUpdate: camera feed back wrong" << endl;
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


int gsins::t_gmultign::_getMeas()
{
	gtrace("t_gintegration::_getMeas");
#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif 
	//	_gmutex.lock();

	_Meas_Type.clear(); Flag = NO_MEAS;
	t_gsinskf::_restore(Xk, Pk);

	if (is_GNSS && _time_valid(_gnss_crt, _ins_crt))
		_Meas_Type.insert(MEAS_TYPE::GNSS_MEAS);
	if (is_VISION && _img_time_valid(_img_crt, _ins_crt))
		_Meas_Type.insert(MEAS_TYPE::VIS_MEAS);

	// 2Hz
	if (fmod(sins.pure_ins_time, 0.5) < 1e-6 || fmod(sins.pure_ins_time, 0.5) > 0.5 - 1e-6)
	{
		MEAS_TYPE meas_type = meas_state();
		//if (meas_type != NO_MEAS)_Meas_Type.insert(meas_type);
	}

	if (_mm.size() != 0)
	{
		//_Meas_Type.insert(MEAS_TYPE::ATT_MEAS);
	}

	//	_gmutex.unlock();
	return _Meas_Type.size();
}

int gsins::t_gmultign::_GNSS_Update()
{
	t_gposdata::data_pos posdata;

	if (_ign_type == IGN_TYPE::LCI)
	{
		Flag = _getPOS(posdata);
		if (Flag != NO_MEAS) {
			_meas_update();
		}
	}
	else if (_ign_type == IGN_TYPE::MULTIGN_LC)
	{
		Flag = _getPOS(posdata);
		if (Flag != NO_MEAS) {
			_mult_meas_update();
		}
	}
	else if (_ign_type == IGN_TYPE::STCI)
	{
		if (_valid_ins_constraint())
		{
			Eigen::Vector3d XYZ_INS = sins.pos_ecef + sins.Ceb*lever;
			t_gtriple XYZ(XYZ_INS(0), XYZ_INS(1), XYZ_INS(2));
			double ins_error = SQRT(Pk(6, 6) + Pk(7, 7) + Pk(8, 8));
			t_gtriple XYZ_RMS(0.5, 0.5, 0.5);
			_external_pos(XYZ, XYZ_RMS);
		}
		Flag = _getPOS(posdata);
		if (Flag != NO_MEAS) {
			_meas_update();
		}
	}
	else if (_ign_type == IGN_TYPE::MULTIGN_STC)
	{
		if (_valid_ins_constraint())
		{
			Eigen::Vector3d XYZ_INS = sins.pos_ecef + sins.Ceb*lever;
			t_gtriple XYZ(XYZ_INS(0), XYZ_INS(1), XYZ_INS(2));
			double ins_error = SQRT(Pk(6, 6) + Pk(7, 7) + Pk(8, 8));
			t_gtriple XYZ_RMS(0.5, 0.5, 0.5);
			_external_pos(XYZ, XYZ_RMS);
		}
		Flag = _getPOS(posdata);
		if (Flag != NO_MEAS) {
			_mult_meas_update();
		}
	}
	else if (_ign_type == IGN_TYPE::TCI)
	{
		_timeUpdate(_gnss_crt);
		Flag = MEAS_TYPE(this->_processEpoch(_gnss_crt));
		if (Flag == NO_MEAS) { _reset(); }
	}
	else if (_ign_type == IGN_TYPE::MULTIGN_TC)
	{
		_timeUpdate(_gnss_crt);
		Flag = MEAS_TYPE(this->_processGNSSEpoch(_gnss_crt));
		if (Flag == NO_MEAS) { _reset(); }
	}
	

	if (Flag == NO_MEAS)return -1;

	return 1;
}

int gsins::t_gmultign::_Cam_Update()
{
	if (_ign_type == IGN_TYPE::MULTIGN_LC)
	{
		imgproc->load_imgobs(_img_t, _img_path);
		Flag = MEAS_TYPE(this->_processCamEpoch(_img_crt));
	}
	else if (_ign_type == IGN_TYPE::MULTIGN_STC)
	{
		imgproc->load_imgobs(_img_t, _img_path);
		Flag = MEAS_TYPE(this->_processCamEpoch(_img_crt));
	}
	else if (_ign_type == IGN_TYPE::MULTIGN_TC)
	{
		imgproc->load_imgobs(_img_t, _img_path);
		Flag = MEAS_TYPE(this->_processCamEpoch(_img_crt));
	}
	else if (_ign_type == IGN_TYPE::VIO)
	{
		imgproc->load_imgobs(_img_t, _img_path);
		Flag = MEAS_TYPE(this->_processCamEpoch(_img_crt));
	}
	
	if (Flag == NO_MEAS)return -1;

	return 1;
}

int gsins::t_gmultign::_processCamEpoch(const t_gtime & runEpoch)
{
	double tmp = runEpoch.sow() + runEpoch.dsec();
	_pointcloud = imgproc->ProcessBatch(tmp);

	// TODO checkstaticmotion
	// need more test,this should be used when GNSS is avaliable
	if (!checkStaticMotion())
	{
		imu_state_id = t_gvis::next_id++; //put this ahead

		//not used function
		keyframeCheck(_pointcloud);
		cout << "keyframe:" << KeyFrameBuffer.size() << endl;
		cout << "nonkeyframe:" << nonKeyFrameBuffer.size() << endl;


		stateAugmentation(runEpoch);

		addFeatureObservations(_pointcloud);

		removeLostFeatures();

		pruneCamStateBuffer();

	}
	return VIS_MEAS;

}

int gsins::t_gmultign::_cam_feedback(const ColumnVector & dx)
{
	if (cam_states.size() > 0)
	{
	
		auto cam_state_iter = cam_states.begin();
		for (int i = 0; i < cam_states.size();
			++i, ++cam_state_iter)
		{
			StateIDType cur_id = cam_state_iter->first;
			string _camid = state_id2str(cur_id);
		
			int idx = _param.getParam(_site, t_gpar::CAM_ATT_X, _camid);
			if (idx < 0) cerr << "cam_feedback() idx wrong!" << endl;

			Eigen::Vector3d att_dx = Eigen::Vector3d(dx(idx + 1), dx(idx + 2), dx(idx + 3));
			Eigen::Vector3d crd_dx = Eigen::Vector3d(dx(idx + 4), dx(idx + 5), dx(idx + 6));
			//cout << "cam_datt:" <<setprecision(12)<< att_dx.transpose() << endl;
			//cout << "cam_dx:" <<setprecision(12) <<crd_dx.transpose() << endl;
			t_gquat dq_cam = t_gbase::rv2q(-att_dx);
			Eigen::Quaterniond _dq_cam;
			_dq_cam.w() = dq_cam.q0; _dq_cam.x() = dq_cam.q1; _dq_cam.y() = dq_cam.q2; _dq_cam.z() = dq_cam.q3;
			_dq_cam.normalize();
			cam_state_iter->second.orientation = cam_state_iter->second.orientation * _dq_cam;
			cam_state_iter->second.position -= crd_dx;
		}
	}
	return 1;
}

int gsins::t_gmultign::_cam_feedback(const Eigen::VectorXd & dx)
{
	cout << "cam_feedaback====cam_statesize:" << cam_states.size() << endl;
	cout<< "cam_feedaback====dx.rows:" << dx.rows() << endl;
	//cerr << dx.transpose() << endl;
	int ins_size = 15;
	if (imgproc->estimate_extrinsic) ins_size = 21;

	auto cam_state_iter = cam_states.begin();
	//cerr << (*cam_state_iter).first << endl;

	for (int i = 0; i < cam_states.size();
		++i, ++cam_state_iter)
	{
		const Eigen::VectorXd& delta_x_cam = dx.segment<6>(ins_size + i * 6);
		t_gquat dq_cam = t_gbase::rv2q(-delta_x_cam.head<3>());
		Eigen::Quaterniond _dq_cam;
		_dq_cam.w() = dq_cam.q0; _dq_cam.x() = dq_cam.q1; _dq_cam.y() = dq_cam.q2; _dq_cam.z() = dq_cam.q3;
		_dq_cam.normalize();
		cam_state_iter->second.orientation = cam_state_iter->second.orientation * _dq_cam;
		cam_state_iter->second.position -= delta_x_cam.tail<3>();
		//cout << "STATE:" << (*cam_state_iter).first << delta_x_cam.transpose() << endl;
	}
	
	return 1;
}

string gsins::t_gmultign::state_id2str(const StateIDType& id)
{
	//get string camera state id
	ostringstream os;
	os << id;
	string cam_id;
	istringstream is(os.str());
	is >> cam_id;
	return cam_id;
	
}



void gsins::t_gmultign::cam_meas_update(const Eigen::MatrixXd& H, const Eigen::VectorXd& r)
{
	if (H.rows() == 0 || r.rows() == 0)
		return;
	if (is_allTC)
	{
		if (imgproc->estimate_extrinsic)
			assert(H.cols() == 6 * cam_states.size() + 6);
		else
			assert(H.cols() == 6 * cam_states.size());
		int par_size = _param.parNumber();
		int ins_size = 15;
		int obs_size = H.rows();
		int N = cam_states.size();
		
		Eigen::MatrixXd All_Variance = BaseMatrix2Eigen(_Qx);
		All_Variance.block(0, 0, ins_size, ins_size) = Pk;
		Eigen::MatrixXd All_H = Eigen::MatrixXd::Zero(H.rows(), par_size);

		assert(N == int(H.cols() / 6));
		assert(obs_size == r.rows());


		auto cam_state_iter = cam_states.begin();
		for (int i = 0; i < cam_states.size();
			++i, ++cam_state_iter)
		{
			string cam_id = state_id2str(cam_state_iter->first);
			int idx = _param.getParam(_site, t_gpar::CAM_ATT_X, cam_id);
			All_H.block(0, idx, obs_size, 6) = H.block(0, 6 * i, obs_size, 6);

		}

		int ig_size = _param.parNumber() - 6 * cam_states.size();
		Eigen::MatrixXd R = Eigen::MatrixXd::Identity(obs_size, obs_size)*observation_noise;

		ColumnVector dx;

		Eigen::MatrixXd H_thin;
		Eigen::VectorXd r_thin;
		string cam_id = state_id2str((*cam_states.begin()).first);
		int idx = _param.getParam(_site, t_gpar::CAM_ATT_X, cam_id);
		//cerr << All_H.block(0, 0, obs_size, idx) << endl;
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
		All_Variance = I_KH * All_Variance * I_KH.transpose() + K * K.transpose() * observation_noise;
		dx = VectorXd2Columns(delta_x);
		_Qx = Eigen2BaseMatrix(All_Variance);


		int ign_size = _param.parNumber();
		Xk = Columns2VectorXd(dx).block(0, 0, ins_size, 1);
		Pk = BaseMatrix2Eigen(_Qx).block(0, 0, ins_size, ins_size);

		if (_gnss_feedback(dx) < 0) {
			cerr << "gns—meas-update-fail" << endl;
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

		if (imgproc->estimate_extrinsic)
		{
			int index = _param.getParam(_site, t_gpar::EXTRINSIC_ATT_X, "");
			if (index < 0) cerr << "Extrinsic Update Error!" << endl;
			Eigen::Vector3d att_dx = Eigen::Vector3d(dx(index + 1), dx(index + 2), dx(index + 3));
			Eigen::Vector3d crd_dx = Eigen::Vector3d(dx(index + 4), dx(index + 5), dx(index + 6));
			cout << "index" << index << "\next_datt:" << setprecision(12) << att_dx.transpose() << endl;
			cout << "ext_dx:" << setprecision(12) << crd_dx.transpose() << endl;
			t_gquat dq_cam = t_gbase::rv2q(-att_dx);
			Eigen::Quaterniond _dq_ext;
			_dq_ext.w() = dq_cam.q0; _dq_ext.x() = dq_cam.q1; _dq_ext.y() = dq_cam.q2; _dq_ext.z() = dq_cam.q3;
			_dq_ext.normalize();
			this->R_cam0_imu = R_cam0_imu * _dq_ext;
			this->t_cam0_imu -= crd_dx;
			this->T_cam0_imu.linear() = R_cam0_imu;
			this->T_cam0_imu.translation() = t_cam0_imu;
			imgproc->T_cam0_imu = this->T_cam0_imu;
			imgproc->t_cam0_imu = t_cam0_imu;
		}

		if (_cam_feedback(dx) < 0) cerr << "cam_update_error:" << endl;
	}
	else
	{
		int ins_size = 15;
		if (imgproc->estimate_extrinsic) ins_size = 21;

		Eigen::MatrixXd All_H = Eigen::MatrixXd::Zero(H.rows(), H.cols() + 15);
		All_H.block(0, 15, H.rows(), H.cols()) = H;
		
		Eigen::MatrixXd H_thin;
		Eigen::VectorXd r_thin;
	
		if (All_H.rows() > All_H.cols())
		{
	
			Eigen::HouseholderQR<Eigen::MatrixXd> qr_helper(All_H);
			Eigen::MatrixXd Q = qr_helper.householderQ();
			Eigen::MatrixXd Q1 = Q.leftCols(ins_size + cam_states.size() * 6);
	
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
	
		if (delta_x_imu.segment<3>(6).norm() > 0.5 ||
			delta_x_imu.segment<3>(12).norm() > 1.0)
		{
			printf("delta velocity: %f\n", delta_x_imu.segment<3>(6).norm());
			printf("delta position: %f\n", delta_x_imu.segment<3>(12).norm());
		}

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
		if (imgproc->estimate_extrinsic)
		{
			const Eigen::VectorXd& delta_x_cam = delta_x.segment<6>(15);
			std::cout << delta_x.size() << " dx:" << delta_x_cam.transpose() << "size" << delta_x_cam.size() << endl;
			t_gquat dq_cam = t_gbase::rv2q(-delta_x_cam.head<3>());
			Eigen::Quaterniond _dq_cam;
			_dq_cam.w() = dq_cam.q0; _dq_cam.x() = dq_cam.q1; _dq_cam.y() = dq_cam.q2; _dq_cam.z() = dq_cam.q3;
			_dq_cam.normalize();
			extrinsic.orientation = extrinsic.orientation * _dq_cam;
			extrinsic.position -= delta_x_cam.tail<3>();

			imgproc->T_cam0_imu.linear() = extrinsic.orientation.toRotationMatrix();
			imgproc->T_cam0_imu.translation() = extrinsic.position;
			imgproc->t_cam0_imu = extrinsic.position;
			this->R_cam0_imu = extrinsic.orientation.toRotationMatrix();;
			this->t_cam0_imu = extrinsic.position;
			this->T_cam0_imu = imgproc->T_cam0_imu;
		}

	
		if(!_cam_feedback(delta_x)) cerr<< "cam_update_error:" << endl;
	
	
		Eigen::MatrixXd I_KH = Eigen::MatrixXd::Identity(K.rows(), H_thin.cols()) - K * H_thin;
	
		state_cov = I_KH * state_cov*I_KH.transpose() +
			K * K.transpose()*observation_noise;
		Eigen::MatrixXd state_cov_fixed = (state_cov + state_cov.transpose()) / 2.0;
	
		assign_statecov(state_cov_fixed);
	}
}


bool gsins::t_gmultign::checkStaticMotion()
{
	double fx = cam0_intrinsics(0); double fy = cam0_intrinsics(1);
	double cx = cam0_intrinsics(2); double cy = cam0_intrinsics(3);
	// only uesd for mono
	if (num_of_cam != 1) return false;
	double mean_motion = 0.0;
	int common_feature_size = 0;
	for (auto pts : _pointcloud.features)
	{
		int feature_id = pts.id;
		auto iter = map_server.find(feature_id);

		if (iter != map_server.end())
		{
			double x1, y1, x2, y2;
			x1 = fx * pts.cam0_point.x + cx;
			y1 = fy * pts.cam0_point.y + cy;
			x2 = fx * iter->second.observations[imu_state_id](0) + cx;
			y2 = fy * iter->second.observations[imu_state_id](1) + cy;
			//cout << x1 <<" "<< x2 << endl;
			double dx = x2 - x1;
			double dy = y2 - y1;
			mean_motion += sqrt(dx * dx + dy * dy);
			common_feature_size++;
		}
	}
	mean_motion = mean_motion / common_feature_size;

	if (mean_motion < 1.2)
	{
		cout << " Static Motion" << endl;
		return true;
		cin.get();
	}
	return false;
}

void gsins::t_gmultign::stateAugmentation(const t_gtime& runEpoch)
{
	double run_epoch = runEpoch.sow() + runEpoch.dsec();
	double dt = sins.t - run_epoch;

	//compensate due to time 
	Eigen::Vector3d BLH = sins.pos - sins.eth.v2dp(sins.vn, dt);//BLH
	Eigen::Vector3d XYZ = gsins::Geod2Cart(BLH, false);

	Eigen::Vector3d imu_vel = sins.vn - sins.an*dt;//ENU
	const Eigen::Matrix3d& R_i_c = R_cam0_imu.transpose();
	const Eigen::Vector3d& t_c_i = t_cam0_imu;
	const Eigen::Matrix3d& R_i_n = t_gbase::q2mat(sins.qnb);
	Eigen::Matrix3d R_n_i = R_i_n.transpose();
	Eigen::Matrix3d R_n_e = t_gbase::Cen(BLH);
	Eigen::Matrix3d R_e_i = R_n_i * R_n_e.transpose();
	Eigen::Vector3d t_c_e = XYZ + R_n_e * R_i_n*t_c_i;
	Eigen::Matrix3d R_e_c = R_i_c * R_e_i;
	//maintain the numerical stability of derivatives
	if (mIsFirstImg)
	{
		initCamPos = t_c_e;
		t_c_e = Eigen::Vector3d::Zero();
		mIsFirstImg = false;
	}
	else
	{
		t_c_e = t_c_e - initCamPos;
	}

	cam_states[imu_state_id] = CAMState(imu_state_id);
	
	CAMState& cam_state = cam_states[imu_state_id];

	if (isKeyFrame)
		cam_state.isKeyFrame = true;

	cam_state.time = run_epoch;
	cam_state.orientation = R_e_c.transpose();
	cam_state.position = t_c_e;
	
	if (is_allTC)
	{
		//get string camera state id
		string cam_id = state_id2str(cam_state.id);
		t_gpar att_x_par(_site, t_gpar::CAM_ATT_X, _param.parNumber() + 1, cam_id);
		_param.addParam(att_x_par);
		t_gpar att_y_par(_site, t_gpar::CAM_ATT_Y, _param.parNumber() + 1, cam_id);
		_param.addParam(att_y_par);
		t_gpar att_z_par(_site, t_gpar::CAM_ATT_Z, _param.parNumber() + 1, cam_id);
		_param.addParam(att_z_par);

		t_gpar crd_x_par(_site, t_gpar::CAM_CRD_X, _param.parNumber() + 1, cam_id);
		_param.addParam(crd_x_par);
		t_gpar crd_y_par(_site, t_gpar::CAM_CRD_Y, _param.parNumber() + 1, cam_id);
		_param.addParam(crd_y_par);
		t_gpar crd_z_par(_site, t_gpar::CAM_CRD_Z, _param.parNumber() + 1, cam_id);
		_param.addParam(crd_z_par);

		int ins_size = 15;
		Eigen::MatrixXd All_Variance = BaseMatrix2Eigen(_Qx);
		//time update
		size_t old_rows = All_Variance.rows();
		size_t old_cols = All_Variance.cols();
		All_Variance.block(0, 0, ins_size, ins_size) = Pk;
		All_Variance.block(ins_size, 0, old_rows - ins_size, ins_size) = All_Variance.block(ins_size, 0, old_rows - ins_size, ins_size) * Phik.transpose();
		All_Variance.block(0, ins_size, ins_size, old_rows - ins_size) = Phik * All_Variance.block(0, ins_size, ins_size, old_rows - ins_size);
		Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, old_cols);
		//cam att to imu att
		J.block<3, 3>(0, 0) = -R_e_c;
		//cam pos to imu att
		J.block<3, 3>(3, 0) = t_gbase::askew(R_n_e * R_i_n * t_c_i);
		//cam pos to imu pos
		J.block<3, 3>(3, 6) = Eigen::Matrix3d::Identity();
		if (imgproc->estimate_extrinsic)
		{
			int index = _param.getParam(_site, t_gpar::EXTRINSIC_ATT_X, "");
			J.block<3, 3>(0, index + 1) = Eigen::Matrix3d::Identity();
			//J.block<3, 3>(15, 3) = Eigen::Matrix3d::Zero();
			//J.block<3, 3>(18, 0) = Eigen::Matrix3d::Zero();
			J.block<3, 3>(3, index + 4) = R_e_i.transpose();
		}

		All_Variance.conservativeResize(old_rows + 6, old_cols + 6);
		const Eigen::MatrixXd P11 = All_Variance.block(0, 0, old_rows, old_cols);
		//matrix augment
		All_Variance.block(old_rows, 0, 6, old_cols) = J * P11;
		All_Variance.block(0, old_cols, old_rows, 6) = All_Variance.block(old_rows, 0, 6, old_cols).transpose();
		All_Variance.block<6, 6>(old_rows, old_cols) = J * P11 * J.transpose();
		//Eigen::EigenSolver<Eigen::MatrixXd> es(All_Variance);
		//cout << "The eigenvalues of A are:" << endl<<setprecision(12) << es.eigenvalues() << endl;
		//keep positive definite of matrix
		All_Variance.block<6, 6>(old_rows, old_cols) += Eigen::MatrixXd::Identity(6, 6) * 1e-10;
		_Qx = Eigen2BaseMatrix(All_Variance);
		assert(_Qx.nrows() == _param.parNumber());
	}
	else
	{
		Eigen::MatrixXd J;
		int ins_size;
		if (imgproc->estimate_extrinsic)
		{
			J = Eigen::Matrix<double, 6, 21>::Zero(); ins_size = 21;
		}
		else
		{
			J = Eigen::Matrix<double, 6, 15>::Zero(); ins_size = 15;
		}
		//cam att to imu att
		J.block<3, 3>(0, 0) = -R_e_c;
		//cam pos to imu att
		J.block<3, 3>(3, 0) = t_gbase::askew(R_n_e * R_i_n * t_c_i);
		//cam pos to imu pos
		J.block<3, 3>(3, 6) = Eigen::Matrix3d::Identity();
		if (imgproc->estimate_extrinsic)
		{
			J.block<3, 3>(0, 15) = Eigen::Matrix3d::Identity();
			//J.block<3, 3>(15, 3) = Eigen::Matrix3d::Zero();
			//J.block<3, 3>(18, 0) = Eigen::Matrix3d::Zero();
			J.block<3, 3>(3, 18) = R_e_i.transpose();
		}
		//cout << get_statecov() << endl;

		//covariance timeupdate
		if (cam_states.size() > 0)
		{
			P_IC = Phik * P_IC;
		}

		Eigen::MatrixXd state_cov;
		state_cov = get_statecov();
		size_t old_rows = state_cov.rows();
		size_t old_cols = state_cov.cols();
		state_cov.conservativeResize(old_rows + 6, old_cols + 6);
		//cout << state_cov << endl;

		const Eigen::MatrixXd& P11 =
			state_cov.block(0, 0, ins_size, ins_size);
		const Eigen::MatrixXd& P12 =
			state_cov.block(0, ins_size, ins_size, old_cols - ins_size);

		state_cov.block(old_rows, 0, 6, old_cols) << J * P11, J* P12;
		state_cov.block(0, old_cols, old_rows, 6) =
			state_cov.block(old_rows, 0, 6, old_cols).transpose();
		state_cov.block<6, 6>(old_rows, old_cols) =
			J * P11 * J.transpose();

		//Eigen::EigenSolver<Eigen::MatrixXd> es(state_cov);
		//cout << "The eigenvalues of A are:" << endl <<setprecision(12)<< es.eigenvalues() << endl;
		//state_cov.block<6, 6>(old_rows, old_cols) += state_cov.block<6, 6>(old_rows, old_cols) + Eigen::MatrixXd::Identity(6,6)*1e-10;

		P_CC.conservativeResize(P_CC.rows() + 6, P_CC.cols() + 6);
		P_IC.conservativeResize(P_IC.rows(), P_IC.cols() + 6);
		Eigen::MatrixXd state_cov_fixed = (state_cov +
			state_cov.transpose()) / 2.0;

		assign_statecov(state_cov_fixed);

	}

	int ins_size = 15;
	Phik = Eigen::MatrixXd::Identity(ins_size, ins_size);

}


bool gsins::t_gmultign::gatingTest(const Eigen::MatrixXd& H, const Eigen::VectorXd& r, const int& dof)
{ 
	if (is_allTC)
	{
		int par_size = _param.parNumber();
		int ins_size = 15;
		int obs_size = H.rows();
		int N = cam_states.size();
		Eigen::MatrixXd All_Variance = BaseMatrix2Eigen(_Qx);
		All_Variance.block(0, 0, 15, 15) = Pk;
		Eigen::MatrixXd All_H = Eigen::MatrixXd::Zero(H.rows(), par_size);

		assert(obs_size == r.rows());

		auto cam_state_iter = cam_states.begin();
		for (int i = 0; i < cam_states.size();
			++i, ++cam_state_iter)
		{
			string cam_id = state_id2str(cam_state_iter->first);
			int idx = _param.getParam(_site, t_gpar::CAM_ATT_X, cam_id);
			All_H.block(0, idx, obs_size, 6) = H.block(0, 6 * i, obs_size, 6);
		}
		Eigen::MatrixXd state_cov = All_Variance;
		Eigen::MatrixXd P1 = All_H * state_cov*All_H.transpose();
		Eigen::MatrixXd P2 = observation_noise * Eigen::MatrixXd::Identity(All_H.rows(), All_H.rows());
		// calculate weighted average of the residuals
		double gamma = r.transpose() * (P1 + P2).ldlt().solve(r);
		//cout << r.rows() << " " << gamma << " " << chi_squared_test_table[dof] << " " << endl;
		if (gamma < chi_squared_test_table[dof]) {
			return true;
		}
		else {
			return false;
		}
	}
	else
	{
		int ins_size = 15;
		if (imgproc->estimate_extrinsic) ins_size = 21;

		Eigen::MatrixXd state_cov = get_statecov();
		//cerr << "H: " << H.rows() << " " << H.cols() << " Cam: " << cam_states.size() * 6 << endl;
		Eigen::MatrixXd All_H = Eigen::MatrixXd::Zero(H.rows(), H.cols() + 15);
		All_H.block(0, 15, H.rows(), H.cols()) = H;
		Eigen::MatrixXd P1 = All_H * state_cov*All_H.transpose();
		Eigen::MatrixXd P2 = observation_noise * Eigen::MatrixXd::Identity(All_H.rows(), All_H.rows());
		double gamma = r.transpose() * (P1 + P2).ldlt().solve(r);
		//cout << dof << " " << gamma << " " << chi_squared_test_table[dof] << " " << endl;
		if (gamma < chi_squared_test_table[dof]){
			return true;}
		else{
			return false;}
	}
	
}

void gsins::t_gmultign::removeLostFeatures()
{
	int pointNum = 0;
	int pointInitialfail = 0;
	int pointgatingTestfail = 0;
	int pointUpdate = 0;

	int jacobian_row_size = 0;
	vector<FeatureIDType> invalid_feature_ids(0);
	vector<StateIDType> processed_feature_ids(0);

	for (auto iter = map_server.begin();
		iter != map_server.end(); ++iter)
	{

		auto& feature = iter->second;

		if (!feature.is_KeyFrame)
			continue;
		// lost or not determined in addFeatureObservations()
		if (!feature.isLost)
			continue;

		// change  3->5
		if (feature.observations.size() < 5)
		{
			invalid_feature_ids.push_back(feature.id);
			continue;
		}

		pointNum++;

		if ((!feature.is_initialized) || ((feature.is_initialized) && (feature.is_initialized_NonKey)))
		{
			// TODO not only used in mono lk
			// uniform velocity change of the feature points to remove outliers
			//if (!stereo && !feature.checkMotion(cam_states))
			if (!feature.checkMotion(cam_states))
			{
				invalid_feature_ids.push_back(feature.id);
				cout << "checkmotionfail:" << endl;
				pointInitialfail++;
				continue;


			}
			if (!feature.triangulatePoint(cam_states))
			{
				invalid_feature_ids.push_back(feature.id);
				cout << "triangulatefail:" << endl;

				pointInitialfail++;
				continue;
			}
			else
			{
				#ifdef USE_CERES
				if (!feature.LMOptimizatePositionAndCheck(cam_states))
				{
					invalid_feature_ids.push_back(feature.id);
					//cout << "initposfail:" << endl;
					pointInitialfail++;
					continue;
				}
				else
				{
					feature.is_initialized_NonKey = false;
				}
				#else
				if (!feature.initializePosition(cam_states))
				{
					invalid_feature_ids.push_back(feature.id);
					//cout << "initposfail:" << endl;
					pointInitialfail++;
					continue;
				}
				else
					feature.is_initialized_NonKey = false;
				#endif // USE_CERES
			}
		}


		int observationSize = 0;

		auto it = feature.observations.begin();
		while (it != feature.observations.end())
		{

			if (cam_states[it->first].isKeyFrame)
				observationSize++;

			it++;
		}
		// this situation only appears in keyframe strategy
		if (!stereo && observationSize <= 2)
		{
			invalid_feature_ids.push_back(feature.id);
			
			pointInitialfail++;
			continue;
		}
		if (stereo) jacobian_row_size += 4 * observationSize - 3;
		else jacobian_row_size += 2 * observationSize - 3;
		processed_feature_ids.push_back(feature.id);
	}

	for (const auto& feature_id : invalid_feature_ids)
		map_server.erase(feature_id);

	cout << "processed_feature_ids.size():" << processed_feature_ids.size() << endl;

	if (processed_feature_ids.size() == 0)
		return;

	Eigen::MatrixXd H_x;
	Eigen::VectorXd r;
	if (imgproc->estimate_extrinsic)
	{
		H_x = Eigen::MatrixXd::Zero(jacobian_row_size, 6 * cam_states.size() + 6);
		r = Eigen::VectorXd::Zero(jacobian_row_size);
	}
	else
	{
		H_x = Eigen::MatrixXd::Zero(jacobian_row_size, 6 * cam_states.size());
		r = Eigen::VectorXd::Zero(jacobian_row_size);
	}


	int stack_cntr = 0;


	for (const auto& feature_id : processed_feature_ids)
	{
		auto& feature = map_server[feature_id];

		vector<StateIDType> cam_state_ids(0);

		for (const auto& measurement : feature.observations)
		{
			if (cam_states[measurement.first].isKeyFrame)
				cam_state_ids.push_back(measurement.first);
		}


		Eigen::MatrixXd H_xj;
		Eigen::VectorXd r_j;

		

		if (!featureJacobian(feature.id, cam_state_ids, H_xj, r_j))
		{
			continue;
		}

		if (gatingTest(H_xj, r_j, cam_state_ids.size() - 1))
		{
			H_x.block(stack_cntr, 0, H_xj.rows(), H_xj.cols()) = H_xj;
			r.segment(stack_cntr, r_j.rows()) = r_j;
			stack_cntr += H_xj.rows();
		}
		else
		{
		
			cout << "not_gatingtest:" << endl;
			pointgatingTestfail++;
		}
		if (stack_cntr > 1500)
			break;
	}

	
	pointUpdate = pointNum - pointInitialfail - pointgatingTestfail;
	cout << "remove_jacobian_size:" << stack_cntr << endl;
	H_x.conservativeResize(stack_cntr, H_x.cols());
	r.conservativeResize((Eigen::Index)stack_cntr);

	cam_meas_update(H_x, r);


	for (const auto& feature_id : processed_feature_ids)
	{
		map_server.erase(feature_id);
	}

	cout << "!!!remove_finshed!!!" << endl;
	return;
}


void gsins::t_gmultign::pruneCamStateBuffer()
{
	cout << "prune_start:" << endl;
	if (KeyFrameBuffer.size() < max_camstate_size)
		return;

	int pointNum = 0;
	int pointInitialfail = 0;
	int pointgatingTestfail = 0;
	int pointUpdate = 0;


	vector<StateIDType> rm_cam_state_ids(0);

	findRedundantCamStates(rm_cam_state_ids);


	int jacobian_row_size = 0;
	for (auto& item : map_server)
	{
		auto& feature = item.second;


		vector<StateIDType> involved_cam_state_ids(0);
		for (const auto& cam_id : rm_cam_state_ids)
		{
			if (feature.observations.find(cam_id) !=
				feature.observations.end())
				involved_cam_state_ids.push_back(cam_id);
		}

		if (involved_cam_state_ids.size() == 0)
			continue;
		if (involved_cam_state_ids.size() == 1)
		{
			feature.observations.erase(involved_cam_state_ids[0]);
			continue;
		}
		if ((involved_cam_state_ids.size() == 2) && (!feature.is_KeyFrame))
		{
			feature.observations.erase(involved_cam_state_ids[0]);
			feature.observations.erase(involved_cam_state_ids[1]);
			continue;
		}


		pointNum++;

		if ((!feature.is_initialized) || ((feature.is_initialized) && (feature.is_initialized_NonKey)))
		{
			if (!stereo && !feature.checkMotion(cam_states))
			{

				for (const auto& cam_id : involved_cam_state_ids)
					feature.observations.erase(cam_id);
				cout << "checkmotionfaile:" << endl;
				pointInitialfail++;
				continue;

			}
			if (!feature.triangulatePoint(cam_states))
			{

				for (const auto& cam_id : involved_cam_state_ids)
					feature.observations.erase(cam_id);
				cout << "triangulatefail:" << endl;
				pointInitialfail++;
				continue;
			}
			else
			{
#ifdef USE_CERES
				if (!feature.LMOptimizatePositionAndCheck(cam_states))
				{
					for (const auto& cam_id : involved_cam_state_ids)
						feature.observations.erase(cam_id);
					//cout << "initposfail:" << endl;
					pointInitialfail++;

					continue;
				}
				else
					feature.is_initialized_NonKey = false;
#else
				if (!feature.initializePosition(cam_states))
				{

					for (const auto& cam_id : involved_cam_state_ids)
						feature.observations.erase(cam_id);
					//cout << "initposfail:" << endl;
					pointInitialfail++;

					continue;
				}
				else
					feature.is_initialized_NonKey = false;
#endif // USE_CERES			
			}
		}

		if (stereo) jacobian_row_size += 4 * involved_cam_state_ids.size() - 3;
		else	jacobian_row_size += 2 * involved_cam_state_ids.size() - 3;
	}

	/*Eigen::MatrixXd H_x = Eigen::MatrixXd::Zero(jacobian_row_size,
		15 + 6 * cam_states.size());*/
	Eigen::MatrixXd H_x;
	if (imgproc->estimate_extrinsic)
		H_x = Eigen::MatrixXd::Zero(jacobian_row_size,
			6 * cam_states.size() + 6);
	else
		H_x = Eigen::MatrixXd::Zero(jacobian_row_size,
			6 * cam_states.size());
	
	Eigen::VectorXd r = Eigen::VectorXd::Zero(jacobian_row_size);

	int stack_cntr = 0;
	for (auto& item : map_server)
	{
		auto& feature = item.second;


		vector<StateIDType> involved_cam_state_ids(0);
		for (const auto& cam_id : rm_cam_state_ids)
		{
			if (feature.observations.find(cam_id) !=
				feature.observations.end())
				involved_cam_state_ids.push_back(cam_id);
		}

		if (involved_cam_state_ids.size() == 0)
			continue;


		Eigen::MatrixXd H_xj;
		Eigen::VectorXd r_j;

		if (!featureJacobian(feature.id, involved_cam_state_ids, H_xj, r_j))
		{
			for (const auto& cam_id : involved_cam_state_ids)
				feature.observations.erase(cam_id);
			continue;
		}



		if (gatingTest(H_xj, r_j, involved_cam_state_ids.size()))
		{
			H_x.block(stack_cntr, 0, H_xj.rows(), H_xj.cols()) = H_xj;
			r.segment(stack_cntr, r_j.rows()) = r_j;
			stack_cntr += H_xj.rows();
		}
		else
		{
			cout << "not_gatingtest:" << endl;
			pointgatingTestfail++;
		}
		for (const auto& cam_id : involved_cam_state_ids)
			feature.observations.erase(cam_id);
	}
	//cout << "H_X.rows()"<< stack_cntr << endl;

	H_x.conservativeResize(stack_cntr, H_x.cols());
	r.conservativeResize(stack_cntr);
	cout << "prune_jacobian_size:" << stack_cntr << endl;

	pointUpdate = pointNum - pointInitialfail - pointgatingTestfail;

	cam_meas_update(H_x, r);

	if (is_allTC)
	{
		for (const auto& cam_id : rm_cam_state_ids)
		{
			StateIDType _cam_id = cam_id;
			string string_id = state_id2str(_cam_id);
			int idx_att_x = _param.getParam(_site, t_gpar::CAM_ATT_X, string_id);
			Matrix_remRC(_Qx, _param[idx_att_x].index, _param[idx_att_x].index);
			_param.delParam(idx_att_x);
			_param.reIndex();

			int idx_att_y = _param.getParam(_site, t_gpar::CAM_ATT_Y, string_id);
			Matrix_remRC(_Qx, _param[idx_att_y].index, _param[idx_att_y].index);
			_param.delParam(idx_att_y);
			_param.reIndex();

			int idx_att_z = _param.getParam(_site, t_gpar::CAM_ATT_Z, string_id);
			Matrix_remRC(_Qx, _param[idx_att_z].index, _param[idx_att_z].index);
			_param.delParam(idx_att_z);
			_param.reIndex();

			int idx_crd_x = _param.getParam(_site, t_gpar::CAM_CRD_X, string_id);
			Matrix_remRC(_Qx, _param[idx_crd_x].index, _param[idx_crd_x].index);
			_param.delParam(idx_crd_x);
			_param.reIndex();

			int idx_crd_y = _param.getParam(_site, t_gpar::CAM_CRD_Y, string_id);
			Matrix_remRC(_Qx, _param[idx_crd_y].index, _param[idx_crd_y].index);
			_param.delParam(idx_crd_y);
			_param.reIndex();

			int idx_crd_z = _param.getParam(_site, t_gpar::CAM_CRD_Z, string_id);
			Matrix_remRC(_Qx, _param[idx_crd_z].index, _param[idx_crd_z].index);
			_param.delParam(idx_crd_z);
			_param.reIndex();

			cam_states.erase(cam_id);

		}
	}
	else
	{
		for (const auto& cam_id : rm_cam_state_ids)
		{
			int cam_sequence = std::distance(cam_states.begin(),
				cam_states.find(cam_id));
			if (imgproc->estimate_extrinsic) cam_sequence++;

			int cam_state_start = 15 + 6 * cam_sequence;
			int cam_state_end = cam_state_start + 6;
			Eigen::MatrixXd state_cov = get_statecov();
			if (cam_state_end < state_cov.rows())
			{
				state_cov.block(cam_state_start, 0,
					state_cov.rows() - cam_state_end,
					state_cov.cols()) =
					state_cov.block(cam_state_end, 0,
						state_cov.rows() - cam_state_end,
						state_cov.cols());
			
			
				state_cov.block(0, cam_state_start,
					state_cov.rows(),
					state_cov.cols() - cam_state_end) =
					state_cov.block(0, cam_state_end,
						state_cov.rows(),
						state_cov.cols() - cam_state_end);
			
				state_cov.conservativeResize(
					state_cov.rows() - 6, state_cov.cols() - 6);
			}
			else
			{
				state_cov.conservativeResize(
					state_cov.rows() - 6, state_cov.cols() - 6);
			}
			cam_states.erase(cam_id);

			P_CC.conservativeResize(P_CC.rows() - 6, P_CC.cols() - 6);
			P_IC.conservativeResize(15, P_IC.cols() - 6);

			Eigen::MatrixXd state_cov_fixed = (state_cov +
				state_cov.transpose()) / 2.0;
			assign_statecov(state_cov_fixed);
		}
	}
	return;
}

void gsins::t_gmultign::assign_statecov(const Eigen::MatrixXd &P)
{
	size_t rows = P.rows();
	size_t cols = P.cols();
	assert(rows == cols);
	int N = cam_states.size();
	Pk = P.block(0, 0, 15, 15);

	P_CC = P.block(15, 15, rows - 15, cols - 15);
	P_IC = P.block(0, 15, 15, cols - 15);

}

Eigen::MatrixXd gsins::t_gmultign::get_statecov()
{
	size_t imu_rows = Pk.rows();
	size_t cam_rows = P_CC.rows();
	Eigen::MatrixXd state_cov = Eigen::MatrixXd::Zero(imu_rows + cam_rows, imu_rows + cam_rows);
	state_cov.block(0, 0, imu_rows, imu_rows) = Pk;
	state_cov.block(imu_rows, imu_rows, cam_rows, cam_rows) = P_CC;
	state_cov.block(0, imu_rows, imu_rows, cam_rows) = P_IC;
	state_cov.block(imu_rows, 0, cam_rows, imu_rows) = P_IC.transpose();
	return state_cov;
}

int gsins::t_gmultign::_mult_meas_update(double fading)
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

	int new_nq = nq + 6 * cam_states.size();
	Eigen::VectorXd Pxz = Eigen::VectorXd::Zero(new_nq),
		Kk = Eigen::VectorXd::Zero(new_nq), Hi = Eigen::VectorXd::Zero(new_nq);
	Eigen::MatrixXd new_Pk = get_statecov();
	assert(new_Pk.rows() == new_nq);
	Eigen::MatrixXd new_Hk = Eigen::MatrixXd::Zero(Hk.rows(), new_nq);
	new_Hk.block(0, 0, Hk.rows(), Hk.cols()) = Hk;
	Eigen::VectorXd new_Xk = Eigen::VectorXd::Zero(new_nq);

	// before estimating 
	//if (outlier_detect()) { _reset(); Flag = NO_MEAS; return -1; }

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
	//if (outlier_detect()) { _reset(); Flag = NO_MEAS; return -1; }

	if (Flag != NO_MEAS)
	{
		_cam_feedback(new_Xk);
		assign_statecov(new_Pk);
	}

	return 1;
}












//int gsins::t_gmultign::test_processBatch(const t_gtime & beg, const t_gtime & end, bool prtOut)
//{
//
//	gtrace("t_gintegration::processBatch");
//#ifdef BMUTEX   
//	boost::mutex::scoped_lock lock(_mutex);
//#endif 
//	_gmutex.lock();
//
//	if (_grec == nullptr) {
//		ostringstream os;
//		os << "ERROR: No object found (" << _site << "). Processing terminated!!! " << beg.str_ymdhms() << " -> " << end.str_ymdhms() << endl;
//		if (_log) _log->comment(0, "t_gintegration", os.str());
//		return -1;
//	}
//
//	// get the time
//	_beg_time = (_beg_time < beg) ? beg : _beg_time;
//	_end_time = (_end_time > end) ? end : _end_time;
//
//	if (_isBase) {
//		if (_setCrd() < 0) {
//			if (_log) _log->comment(0, "t_gintegration", "Base station have no coordinate");
//			return -1;
//		}
//	}
//
//	_dop.set_data(_gnav, _gobs, _site);
//	_running = true;
//	_param.setSite(_site);
//	int sign = 1;
//	double subint = 0.1;
//	if (_scale > 0) subint = 1.0 / _scale;
//	if (_sampling > 1) subint = pow(10, floor(log10(_sampling)));
//
//	if (!this->_init()) {
//		if (_log) _log->comment(0, "t_gintegration", "init failed");
//		return -1;
//	}
//	t_gtime tmp = t_gtime::current_time(t_gtime::GPS);
//	if (!_beg_end) {
//		tmp = _ins_end; _ins_end = _ins_beg; _ins_beg = tmp;
//		tmp = _gnss_end; _gnss_end = _gnss_beg; _gnss_beg = tmp;
//		_ins_crt = _ins_beg;
//		sign = -1;
//		if (_log) _log->comment(2, "t_gintegration", "Filtering in end -> begin direction!");
//	}
//	else {
//		sign = 1;
//		if (_log) _log->comment(2, "t_gintegration", "Filtering in begin -> end direction!");
//	}
//
//	// need execute frequently in real-time mode
//	_gobs->setepoches(_site);
//
//	// time align
//	if (_ins_beg < _gnss_beg) {
//		_ins_crt = _imudata->erase_bef(_gnss_beg);
//	}
//	else
//	{
//		//t_out("Qx:", _Qx);
//		t_gpvtflt_newmodel::processBatch(_gnss_crt, _ins_crt, false);
//		// next gnss time 
//		while (_gnss_crt < _ins_crt)
//		{
//			if (_sampling > 1) _gnss_crt.add_secs(int(sign * _sampling));  // =<1Hz data
//			else               _gnss_crt.add_dsec(sign * _sampling);       //  >1Hz data
//		}
//	}
//	//t_out("Qx:", _Qx);
//
//	std::cout << _site << ": Start Processing: " << " " << beg.str_ymdhms() << " " << end.str_ymdhms() << endl;
//
//	t_gposdata::data_pos posdata;
//	bool time_loop = true;
//	//lsy add
//	bool first_loop = true;
//	// ins time loop
//	while (time_loop)
//	{
//		if (_ins_crt > _ins_end || _gnss_crt > _gnss_end)
//		{
//			cout << " Processing Finish!" << endl;
//			if (_log) _log->comment(0, "t_gintegration", "Processing Finish!");
//			time_loop = false; break;
//		}
//		_imudata->load(_wm, _vm, _shm.t, _shm.ts, _shm.nSamples, true);
//
//		_ins_crt = t_gtime(_gnss_crt.gwk(), _shm.t);
//		//_img_crt = t_gtime(_gnss_crt.gwk(), _img_t);
//
//		if (!_aligned)
//		{
//			if (_time_valid(_gnss_crt, _ins_crt))
//				Flag = _getPOS(posdata);
//			// todo: modified 
//			_aligned = fsiam(posdata.pos, false);
//		}
//		else
//		{
//
//			_imgdata->load(_shm.t, _img_t, _img_path);
//			_img_crt = t_gtime(_gnss_crt.gwk(), _img_t);
//
//			sins.Update(_wm, _vm, _shm);
//
//			if (_ign_type == PURE_INS) {
//				t_gsinskf::write();
//				continue;
//			}
//
//			time_update(sins.nts, 0.0);
//			Phi= Eigen::MatrixXd::Identity(nq, nq) + (Ft*sins.nts);
//			test_cam_time_update();
//
//
//			if (_ign_type == IGN_TYPE::VIO || _ign_type == IGN_TYPE::MULTIGN)
//			{
//				imgproc->load_imuobs(_shm.t, _wm, _vm);
//			}
//
//		     //  bool flag_gnss = false;
//				bool flag_gnss = _time_valid(_gnss_crt, _ins_crt);
//			//	bool flag_img = false;
//			    bool flag_img = _img_time_valid(_img_crt, _ins_crt);
//
//
//			if (flag_img || flag_gnss)
//			{
//				cout << endl;
//				cout << "imu_time:" << setprecision(10) << _shm.t << endl;
//				cout << "img_now:" << setprecision(10) << _img_t << "  " << flag_img << endl;
//				cout << "gns_now:" << setprecision(10) << (_gnss_crt.sow() + _gnss_crt.dsec()) << "  " << flag_gnss << endl;
//
//
//			}
//			if (flag_gnss && flag_img)
//			{
//				if (_ign_type == IGN_TYPE::MULTIGN)
//				{
//			
//					imgproc->load_imgobs(_img_t, _img_path);
//					Flag = MEAS_TYPE(this->test_processCamEpoch(_img_crt));
//
//					//double tmp_gnss = _gnss_crt.sow() + _gnss_crt.dsec();
//					//if (tmp_gnss <= 442961.6 || tmp_gnss >= 443061.6)
//					//{
//						_timeUpdate(_gnss_crt);
//						Flag = MEAS_TYPE(this->test_processGNSSEpoch(_gnss_crt));
//						//Flag = MEAS_TYPE(this->_processEpoch(_gnss_crt));
//						if (Flag == NO_MEAS) { _reset(); }
//						else { feedback(glv.INF); }
//					//}
//				}
//				else
//				{
//					cout << "IGN_TYPE_ERROR" << endl;
//					cin.get();
//				}
//			}
//			else if (flag_gnss && !flag_img)
//			{
//				if (_ign_type == IGN_TYPE::TCI || _ign_type == IGN_TYPE::MULTIGN)
//				{
//					//double tmp_gnss = _gnss_crt.sow() + _gnss_crt.dsec();
//					//if (tmp_gnss <= 442961.6 || tmp_gnss >= 443061.6)
//					//{
//						_timeUpdate(_gnss_crt);
//						Flag = MEAS_TYPE(this->test_processGNSSEpoch(_gnss_crt));
//						//Flag = MEAS_TYPE(this->_processEpoch(_gnss_crt));
//						if (Flag == NO_MEAS) { _reset(); }
//						else { feedback(glv.INF); }
//					//}
//				}
//				else
//				{
//					cout << "IGN_TYPE_ERROR" << endl;
//					cin.get();
//				}
//			}
//			else if (!flag_gnss && flag_img)
//			{
//				if (_ign_type == IGN_TYPE::VIO || _ign_type == IGN_TYPE::MULTIGN)
//				{
//					imgproc->load_imgobs(_img_t, _img_path);
//					Flag = MEAS_TYPE(this->test_processCamEpoch(_img_crt));
//
//				}
//				else
//				{
//					cout << "IGN_TYPE_ERROR" << endl;
//					cin.get();
//				}
//
//
//			}
//			if (flag_gnss || flag_img)
//			{
//				t_gsinskf::write();
//				if (_log) _log->comment(0, "t_gintegration",
//					_ins_crt.str_ymdhms(" epoch ") + "Processing GNSS/INS Integration Success!");
//
//				Eigen::Vector3d Car_pos = Geod2Cart(sins.pos, false);
//				Eigen::Vector3d Ve = t_gbase::Cen(sins.pos)*(sins.vn);
//				std::cout << fixed << setprecision(4) << sins.t;
//				std::cout << fixed << setprecision(3) <<
//					" " << Car_pos(0) << " "
//					" " << Car_pos(1) << " "
//					" " << Car_pos(2) << " ";
//				std::cout << fixed << setprecision(3) <<
//					" " << Ve(0) << " "
//					" " << Ve(1) << " "
//					" " << Ve(2) << " " << endl;
//			}
//
//
//			//use for glfw 
//
//			imu_state.orientation = Eigen::Quaterniond(t_gbase::q2mat(sins.qnb));
//			imu_state.position = Geod2Cart(sins.pos, false);
//
//			if (_ins_crt <= _ins_end)
//				publisher.UpdateNewState(imu_state, cam_states, initCamPos);
//		}
//		while (_ins_crt.diff(_gnss_crt) > _sampling)
//		{
//			if (_sampling > 1) _gnss_crt.add_secs(int(_sampling));  // =<1Hz data
//			else               _gnss_crt.add_dsec(_sampling);       //  >1Hz data
//		}
//
//	}
//
//
//	_gmutex.unlock();
//	return 0;
//}



//int gsins::t_gmultign::test_processGNSSEpoch(const t_gtime & runEpoch)
//{
//	gtrace("t_gintegration::_processEpoch");
//#ifdef BMUTEX   
//	boost::mutex::scoped_lock lock(_mutex);
//#endif 
//	_gmutex.lock();
//
//#ifdef DEBUG
//	cout << "t_gintegration: Epoch - full data size: " << _data.size() << " "
//		<< runEpoch.str_hms() << endl;
//#endif
//
//	vector<t_gsatdata>::iterator ItException = _data.end();
//
//	double tmp_gnss = runEpoch.sow() + runEpoch.dsec();
//	//if (tmp_gnss >= 185956.3 && tmp_gnss <=717 185976.3)
//	//if (tmp_gnss >= 442961.6 && tmp_gnss <= 443061.6)
//	//{
//	//	while (_data.size() > 5)
//	//		_data.pop_back();
//
//	//	cout << _data.size() << endl;
//	//}
//
//
//	ColumnVector Threshold;
//	if (_grec == nullptr) {
//		if (_log) _log->comment(1, "t_gintegration", "No receiver settings available!!!");
//		_gmutex.unlock(); return NO_MEAS;
//	}
//
//	_epoch = runEpoch;
//	t_gtriple crdapr = _grec->crd_arp(_epoch);
//
//	//if (double_eq(crdapr[0], 0.0) && double_eq(crdapr[1], 0.0) &&
//	//	double_eq(crdapr[2], 0.0)) {
//	//	_valid_crd_xml = false;
//	//}
//	//else {
//	//	_valid_crd_xml = true;
//	//}
//
//	if (!_valid_crd_xml) _sig_init_crd = 100.0;
//
//	Matrix A;
//	SymmetricMatrix P;
//	ColumnVector l, dx, dx_fix;
//	SymmetricMatrix Qsav, QsavBP;
//	vector<t_gsatdata>::iterator it;
//	t_gallpar XsavBP;
//	ColumnVector v_norm, v_orig;
//	string outlier = "";
//	double vtpv;
//	int nobs_total, npar_number;
//
//	//做仿真 在这里进行
//
//	_cntrep = 0;    // number of iterations caused by outliers
//
//	P_IC_sav = P_IC;
//	P_CC_sav = P_CC;
//	do
//	{
//
//		// erase sat satellite because of outlier
//		it = _data.begin();
//		while (it != _data.end()) {
//			if (it->sat() == outlier) {
//				it = _data.erase(it); continue;
//			}
//			it++;
//		}
//
//		Eigen::Vector3d XYZ_INS = sins.pos_ecef + sins.Ceb*lever;
//		t_gtriple XYZ(XYZ_INS(0), XYZ_INS(1), XYZ_INS(2));
//		_external_pos(XYZ, t_gtriple());
//
//		if (_prepareData() < 0) {
//			_predict();//todo:add time gap
//			_gmutex.unlock(); return NO_MEAS;
//		}
//
//		QsavBP = _Qx; XsavBP = _param;
//		_predict();
//
//		_initialized = true;
//
//		if (_data.size() < _minsat) {
//			if (_log) _log->comment(2, "t_gintegration", "Not enough visible satellites!");
//			_restore(QsavBP, XsavBP);
//			_gmutex.unlock(); return NO_MEAS;
//		}
//
//		// define a number of measurements 
//		unsigned int nObs = _data.size();
//		unsigned int mult = 1;
//		if (_observ == RAW_DOUBLE) { mult = 2; nObs *= 2; }
//		if (_observ == RAW_ALL) { mult = 2; nObs *= 5; } // reservation for 5 freq - not used raws will be removed
//		if (_phase) { mult *= 2; nObs *= 2; }  // code + phase
//
//		unsigned int nPar = _param.parNumber();
//
//		A.ReSize(nObs, nPar); A = 0.0;
//		l.ReSize(nObs);       l = 0.0;
//		P.ReSize(nObs);       P = 0.0;
//		dx.ReSize(nPar);     dx = 0.0;
//		unsigned int iobs = 1;
//		_frqNum.clear(); _obs_index.clear();
//
//		Matrix A_new;
//		SymmetricMatrix P_new;
//		ColumnVector l_new, dx_new;
//		unsigned int iobs_new = 1;
//		_frqNum.clear(); _obs_index.clear();
//		if (_isBase)
//		{
//			dynamic_cast<t_gcombDD*>(&(*_base_model))->set_base_data(&_data_base);
//			dynamic_cast<t_gcombDD*>(&(*_base_model))->set_rec_info(_vBanc[3], _vBanc_base[3]);
//		}
//		// use combmodel
//		t_glsqEquationMatrix equ;
//		it = _data.begin();
//		for (it = _data.begin(); it != _data.end(); ) {
//			int irc = _base_model->cmb_equ(_epoch, _param, *it, equ);
//			if (irc < 0) { it = _data.erase(it); continue; }
//			if (_observ == IONO_FREE) _combineMW(*it);
//			++it;
//		}
//		iobs_new = equ.num_equ();
//		// tranform equ into A, P, l
//		_equationToMatrix(equ, A_new, P_new, l_new, nPar);
//
//		// generate obs_index
//		_obs_index.clear();
//		_generateObsIndex(equ);
//
//		if (iobs_new <= _minsat * mult) {
//			if (_log) _log->comment(2, "gpvtflt", "Not enough processable observations!");
//			_restore(QsavBP, XsavBP);
//			_gmutex.unlock(); return NO_MEAS;
//		}
//
//		A = A_new;
//		P = P_new;
//		l = l_new;
//		dx = dx_new;
//
//		if (_isBase || _sd_sat) {
//			if (_combineDD(A, P, l) < 0) {
//				_gmutex.unlock(); return NO_MEAS;
//			}
//		}
//		Qsav = _Qx;
//
//		//if (_merge(A, P, l) < 0) {
//		//	dx = 0.0; _Qx = Qsav;
//		//	if (_log) _log->comment(2, "gintegration", "merge failed!");
//		//	_gmutex.unlock(); return NO_MEAS;
//		//}
//
//
//		//if (_multmerge(A) < 0) {
//		//	dx = 0.0; _Qx = Qsav;
//		//	if (_log) _log->comment(2, "gintegration", "merge failed!");
//		//	_gmutex.unlock(); return NO_MEAS;
//		//}
//		if (_merge(A) < 0) {
//			dx = 0.0; _Qx = Qsav;
//			if (_log) _log->comment(2, "gintegration", "merge failed!");
//			_gmutex.unlock(); return NO_MEAS;
//		}
//
//		try
//		{
//			// 把视觉参数想办法塞进去
//			int ins_size = 15;
//			int ig_size = _param.parNumber();
//			int N = 6 * cam_states.size();
//			int obs_size = A.Nrows();
//
//			//H阵扩维 
//			Eigen::MatrixXd H = NewMat2Eigen(A);
//			H.conservativeResize(obs_size, ig_size + N);
//			H.block(0, ig_size, obs_size, N) = Eigen::MatrixXd::Zero(obs_size, N);
//
//			//P阵扩维
//			Eigen::MatrixXd Variance = BaseMatrix2Eigen(_Qx);
//			Variance.conservativeResize(ig_size + N, ig_size + N);
//			Variance.block(ig_size, ig_size, N, N) = P_CC_sav;
//			Variance.block(0, ig_size, ig_size, N) = Eigen::MatrixXd::Zero(ig_size, N);
//			Variance.block(ig_size, 0, N, ig_size) = Eigen::MatrixXd::Zero(N, ig_size);
//			Variance.block(0, ig_size, ins_size, N) = P_IC_sav;
//			Variance.block(ig_size, 0, N, ins_size) = P_IC_sav.transpose();
//
//			//if (cam_states.size() > 0) cout << "P_IC_sav:" << P_IC_sav << endl;
//			//l
//			Eigen::VectorXd res = Columns2VectorXd(l);
//
//			Eigen::MatrixXd R = BaseMatrix2Eigen(P).inverse();
//			Eigen::MatrixXd NN = R + H * Variance*H.transpose();
//			Eigen::MatrixXd K = Variance * H.transpose()*NN.inverse();
//			Eigen::MatrixXd KH = K * H;
//			Eigen::MatrixXd I_KH = Eigen::MatrixXd::Identity(K.rows(),H.cols()) - KH;
//			Eigen::VectorXd DX = K * res;
//			//cout << "DX:" << DX.transpose() << endl;
//			Variance = I_KH * Variance*I_KH.transpose() + K * R*K.transpose();
//
//			dx = VectorXd2Columns(DX.block(0, 0, ig_size, 1));
//			_Qx = Eigen2BaseMatrix(Variance.block(0, 0, ig_size, ig_size));
//
//
//			Dx_sav = DX;
//			Dx_row = ig_size;
//			//用Eigen来解算
//
//
//
//			////从update函数中拿出来
//			//Matrix K;            // Kalman gain
//			//Matrix NN;
//			//SymmetricMatrix Pli = P.i();
//
//			//NN = Pli + A * _Qx * A.t();
//			//K = _Qx * A.t() * NN.i();
//			////cout << dx << endl;
//			////  cout  << "Qx" << endl << Qx << endl << "KAQx" << endl << K*A*Qx << endl << "K" << K << endl;
//			//IdentityMatrix I = IdentityMatrix(_Qx.nrows());
//			//Matrix KA = K * A;
//			//Matrix I_KA = I - K * A;
//
//			//// ZhouYuxuan modified. "K * ( l - A * dx )" => "K * l"
//			//dx = K * l;               // update state vector
//			//_Qx << I_KA * _Qx * I_KA.t() + K * Pli *K.t();
//			//_filter->update(A, P, l, dx, _Qx);
//
//		}
//		catch (...)
//		{
//			if (_log) _log->comment(2, "gintegration", "meas_update failed!");
//			dx = 0.0; _Qx = Qsav;
//			Flag = NO_MEAS;
//			_reset();
//			return Flag;
//		}
//
//		// increasing variance after update in case of introducing new ambiguity
//		for (size_t iPar = 0; iPar < _param.parNumber(); iPar++) {
//			if (_param[iPar].parType == t_gpar::AMB_IF) {
//				string sat = _param[iPar].prn;
//				if (_newAMB.find(sat) != _newAMB.end() && _cntrep == 1) {
//					if (_newAMB[sat] == 1) _Qx(iPar + 1, iPar + 1) += 10;
//					if (_newAMB[sat] == 2 && _Qx(iPar + 1, iPar + 1) > 0.01) _Qx(iPar + 1, iPar + 1) += 1;
//					_newAMB[sat]++;
//				}
//			}
//		}
//
//		if (_isBase)
//		{
//			if (_observ == RAW_ALL)Flag = DD_UC_OBS_MEAS;
//			if (_observ == IONO_FREE)Flag = DD_IF_OBS_MEAS;
//		}
//		else if (_sd_sat)
//		{
//			if (_observ == IONO_FREE)Flag = SD_IF_OBS_MEAS;
//			if (_observ == RAW_ALL)Flag = SD_UC_OBS_MEAS;
//		}
//		else
//		{
//			if (_observ == RAW_ALL)Flag = UD_UC_OBS_MEAS;
//			if (_observ == IONO_FREE)Flag = UD_IF_OBS_MEAS;
//		}
//
//		//
//		//if (_meas_update() < 0)
//		//{
//		//	if (_log) _log->comment(2, "gintegration", "meas_update failed!");
//		//	Flag = NO_MEAS;
//		//	_reset();
//		//	return Flag;
//		//}
//
//		// post-fit residuals
//		//cout << "pos_fit_residual:" << endl;
//		//t_out("dl", l);
//		//t_out("dx", dx);
//		v_orig = l - A * dx;
//
//		ColumnVector vtPv = v_orig.t() * P * v_orig;
//		vtpv = vtPv(1);
//		nobs_total = A.Nrows();
//		npar_number = A.ncols();
//		int freedom = A.Nrows() - A.Ncols();
//		_sig_unit = vtPv(1) / freedom;
//
//		// normalized post-fit residuals
//		Matrix Qv = A * _Qx * A.t() + P.i();
//
//		v_norm.ReSize(v_orig.Nrows());
//		for (int i = 1; i <= v_norm.Nrows(); i++) {
//			v_norm(i) = sqrt(1 / Qv(i, i)) * v_orig(i);
//		}
//
//
//	} while (_outlierDetect(v_norm, Qsav, outlier) != 0);
//
//	//_gnss_feedback();
//	if (_gnss_feedback(dx) < 0) {
//		if (_log) _log->comment(2, "gintegration", "meas_update failed!");
//		dx = 0.0; _Qx = Qsav;
//		Flag = NO_MEAS;
//		return Flag;
//	}
//
//	//cam parameter change
//	cout << "GNSS_Meas_UP:" << endl;
//	//if (_cam_feedback(dx) < 0)
//	//{
//	//	cerr << "GNSS MeasUpdate: camera feed back wrong" << endl;
//	//}
//	if (test_cam_feedback() < 0)
//	{
//		cerr << "GNSS MeasUpdate: camera feed back wrong" << endl;
//	}
//
//	Phik = Eigen::MatrixXd::Identity(15, 15);
//
//	_filter->add_data(_param, dx, _Qx, _sig_unit, Qsav);
//	_filter->add_data(A, P, l);
//	_filter->add_data(vtpv, nobs_total, npar_number);
//	//_amb_resolution(dx_fix);
//
//	// todo: the fixed-dx must be transmited to t_gsinskf::write
//
//
//	//ColumnVector dx_tmp = dx;
//	//SymmetricMatrix   Qx_tmp = _Qx;
//	//t_gallpar _param_tmp = _param;
//	//// ambfix added by zhshen
//	//if (_fix_mode != t_gsetamb::FIX_MODE::NO)
//	//{
//	//	if (_isFirstFix)
//	//	{
//	//		_ambfix->setUPD(_gupd);
//	//		_ambfix->setLOG(_log);
//	//		_ambfix->setFreq(_gnav);
//	//		_isFirstFix = false;
//	//	}
//	//
//	//	_filter->add_data(_param, dx, _Qx, _sig_unit, Qsav);
//	//	_filter->add_data(A, P, l);
//	//	_filter->add_data(vtpv, nobs_total, npar_number);
//	//	_ambfix->setObsType(_observ);
//	//	_ambfix->setActiveAmb(npar_number);
//	//	if (_observ == IONO_FREE)	_ambfix->setMW(_MW[_epoch]);
//	//
//	//	if (_observ == RAW_ALL && !_isBase&&_frequency >= 3)
//	//	{
//	//		_ambfix->processBatch(_epoch, _filter, nullptr, "EWL");
//	//	}
//	//	if (_observ == RAW_ALL && !_isBase)
//	//	{
//	//		_ambfix->processBatch(_epoch, _filter, nullptr, "WL");
//	//	}
//	//	if (_ambfix->processBatch(_epoch, _filter, nullptr, "NL") < 0)
//	//	{
//	//		_amb_state = false;
//	//		if (_log) _log->comment(2, "gpvtflt", _site + _epoch.str_ymdhms(" epoch ")
//	//			+ ": fix ambiguity failed !");
//	//	}
//	//	else _amb_state = _ambfix->amb_fixed();
//	//
//	//	dx_tmp = _filter->dx();
//	//	Qx_tmp = _filter->Qx();
//	//
//	//}
//
//	_gmutex.unlock();
//	return Flag;
//}
//
//
//
//int gsins::t_gmultign::test_processCamEpoch(const t_gtime & runEpoch)
//{
//	cout << "!!!TEST_PROCESS_CAM_EPOCH!!!" << endl;
//
//	double tmp = runEpoch.sow() + runEpoch.dsec();
//	_pointcloud = imgproc->ProcessBatch(tmp);
//
//	//cout << "before_PK:" << endl;
//	//for (int i = 0; i < 15; i++) cout << Pk(i, i) << endl;
//	//cout << "before_PCC:" << endl;
//	//for (int i = 0; i < P_CC.rows(); i++)cout << P_CC(i, i) << endl;
//
//	bool little0 = false;
//	for (int i = 0; i < 15; i++)
//	{
//		if (Pk(i, i) < 0)
//		{
//			cout << "before_cam_update  Pk<0" << endl;
//			
//			little0 = true;
//		}
//	}
//	if (little0 == true)
//	{
//		cout << setprecision(8) << Pk << endl;
//		cin.get();
//	}
//
//	//for (int i = 0; i < P_CC.rows(); i++)
//	//{
//	//	if (P_CC(i, i) < 0)
//	//	{
//	//		cout << "before_cam_update  Pcc<0" << endl;
//	//		cin.get();
//	//	}
//	//}
//
//
//	// TODO checkstaticmotion
//	// need more test,this should be used when GNSS is avaliable
//	if (!checkStaticMotion())
//	{
//		imu_state_id = t_gvis::next_id++; //put this ahead
//
//		keyframeCheck(_pointcloud);
//		cout << "keyframe:" << KeyFrameBuffer.size() << endl;
//		cout << "nonkeyframe:" << nonKeyFrameBuffer.size() << endl;
//
//
//		test_stateAugmentation(runEpoch);
//
//		addFeatureObservations(_pointcloud);
//
//		test_removeLostFeatures();
//
//		test_pruneCamStateBuffer();
//
//	}
//	//cout << "ater_PK:" << endl;
//	//for (int i = 0; i < 15; i++) cout << Pk(i, i) << endl;
//	//cout << "after_PCC:" << endl;
//	//for (int i = 0; i < P_CC.rows(); i++)cout << P_CC(i, i) << endl;
//
//
//	//for (int i = 0; i < 15; i++)
//	//{
//	//	if (Pk(i, i) < 0)
//	//	{
//	//		cout << "after_cam_update  Pk<0" << endl;
//	//		cin.get();
//	//	}
//	//}
//
//	//for (int i = 0; i < P_CC.rows(); i++)
//	//{
//	//	if (P_CC(i, i) < 0)
//	//	{
//	//		cout << "after_cam_update  Pcc<0" << endl;
//	//		cin.get();
//	//	}
//	//}
//
//
//
//	return VIS_MEAS;
//
//}
//
//
//
//int gsins::t_gmultign::test_cam_feedback()
//{
//	//另外这里需要对cam的参数进行修改
//	auto cam_state_iter = cam_states.begin();
//	for (int i = 0; i < cam_states.size();
//		++i, ++cam_state_iter)
//	{
//		const Eigen::VectorXd& delta_x_cam = Dx_sav.segment<6>(Dx_row + i * 6);
//		t_gquat dq_cam = t_gbase::rv2q(-delta_x_cam.head<3>());
//		Eigen::Quaterniond _dq_cam;
//		_dq_cam.w() = dq_cam.q0; _dq_cam.x() = dq_cam.q1; _dq_cam.y() = dq_cam.q2; _dq_cam.z() = dq_cam.q3;
//		_dq_cam.normalize();
//		cam_state_iter->second.orientation = cam_state_iter->second.orientation * _dq_cam;
//		cam_state_iter->second.position -= delta_x_cam.tail<3>();
//
//		cout << "cam_dq:" << setprecision(12) << delta_x_cam.head<3>().transpose() << endl;
//		cout << "cam_dpos_ecef:" << setprecision(12) << delta_x_cam.tail<3>() << endl;
//	}
//	return 1;
//}
//
//void gsins::t_gmultign::test_stateAugmentation(const t_gtime& runEpoch)
//{
//	double run_epoch = runEpoch.sow() + runEpoch.dsec();
//	double dt = sins.t - run_epoch;
//	cout << "dt:" << setprecision(4) << dt << endl;
//	//把惯导向前推0.0025s
//	Eigen::Vector3d BLH = sins.pos - sins.eth.v2dp(sins.vn, dt);//BLH
//	Eigen::Vector3d XYZ = gsins::Geod2Cart(BLH, false);
//
//	Eigen::Vector3d imu_vel = sins.vn - sins.an*dt;//ENU
//	const Eigen::Matrix3d& R_i_c = R_cam0_imu.transpose();
//	const Eigen::Vector3d& t_c_i = t_cam0_imu;
//	const Eigen::Matrix3d& R_i_n = t_gbase::q2mat(sins.qnb);
//	Eigen::Matrix3d R_n_i = R_i_n.transpose();
//	Eigen::Matrix3d R_n_e = t_gbase::Cen(BLH);
//	Eigen::Matrix3d R_e_i = R_n_i * R_n_e.transpose();
//	Eigen::Vector3d t_c_e = XYZ + R_n_e * R_i_n*t_c_i;
//	Eigen::Matrix3d R_e_c = R_i_c * R_e_i;
//	//第一帧影像位置赋0
//	if (mIsFirstImg)
//	{
//		initCamPos = t_c_e;
//		t_c_e = Eigen::Vector3d::Zero();
//		mIsFirstImg = false;
//	}
//	else
//	{
//		t_c_e = t_c_e - initCamPos;
//	}
//
//	cam_states[imu_state_id] = CAMState(imu_state_id);
//	CAMState& cam_state = cam_states[imu_state_id];
//
//	if (isKeyFrame)
//		cam_state.isKeyFrame = true;
//
//	cam_state.time = run_epoch;
//	cam_state.orientation = R_e_c.transpose();
//
//	cam_state.position = t_c_e;
//	Eigen::Matrix<double, 6, 15> J = Eigen::Matrix<double, 6, 15>::Zero();
//
//	////相机姿态对惯导姿态
//	//J.block<3, 3>(0, 0) = -R_i_c * (R_i_n.transpose());
//	////相机位置对惯导姿态
//	//J.block<3, 3>(3, 0) = R_n_e * t_gbase::askew(R_i_n*t_c_i);
//	////相机位置对惯导位置
//	//J.block<3, 3>(3, 6) = t_gbase::dGeod2Cart(sins.eth, BLH);
//
//	//相机姿态对惯导姿态 //R_i_c->R_e_c;
//	J.block<3, 3>(0, 0) = -R_e_c;
//	//相机位置对惯导姿态
//	J.block<3, 3>(3, 0) = t_gbase::askew(R_n_e*R_i_n*t_c_i);
//	//相机位置对惯导位置
//	J.block<3, 3>(3, 6) = Eigen::Matrix3d::Identity();
//
//	Eigen::MatrixXd state_cov;
//	state_cov = test_get_statecov();
//	size_t old_rows = state_cov.rows();
//	size_t old_cols = state_cov.cols();
//	state_cov.conservativeResize(old_rows + 6, old_cols + 6);
//
//	const Eigen::Matrix<double, 15, 15>& P11 =
//		state_cov.block<15, 15>(0, 0);
//	const Eigen::MatrixXd& P12 =
//		state_cov.block(0, 15, 15, old_cols - 15);
//
//	state_cov.block(old_rows, 0, 6, old_cols) << J * P11, J*P12;
//	state_cov.block(0, old_cols, old_rows, 6) =
//		state_cov.block(old_rows, 0, 6, old_cols).transpose();
//	state_cov.block<6, 6>(old_rows, old_cols) =
//		J * P11 * J.transpose();
//
//	//Eigen::EigenSolver<Eigen::MatrixXd> es(state_cov);
//	//cout << "The eigenvalues of A are:" << endl <<setprecision(12)<< es.eigenvalues() << endl;
//	//state_cov.block<6, 6>(old_rows, old_cols) += state_cov.block<6, 6>(old_rows, old_cols) + Eigen::MatrixXd::Identity(6,6)*1e-10;
//	
//	P_CC.conservativeResize(P_CC.rows() + 6, P_CC.cols() + 6);
//	P_IC.conservativeResize(P_IC.rows(), P_IC.cols() + 6);
//	Eigen::MatrixXd state_cov_fixed = (state_cov +
//		state_cov.transpose()) / 2.0;
//
//	test_assign_statecov(state_cov_fixed);
//	//Pk = state_cov.block(0, 0, 15, 15);
//
//
//}
//
//

//

//void gsins::t_gmultign::test_removeLostFeatures()
//{
//	cout << "remove_start:" << endl;
//	int pointNum = 0;
//	int pointInitialfail = 0;
//	int pointgatingTestfail = 0;
//	int pointUpdate = 0;
//
//	int jacobian_row_size = 0;
//	vector<FeatureIDType> invalid_feature_ids(0);
//	vector<StateIDType> processed_feature_ids(0);
//
//	for (auto iter = map_server.begin();
//		iter != map_server.end(); ++iter)
//	{
//
//		auto& feature = iter->second;
//
//		if (!feature.is_KeyFrame)
//			continue;
//		// lost or not determined in addFeatureObservations()
//		if (!feature.isLost)
//			continue;
//
//		//lsy change  3->5
//		if (feature.observations.size() < 5)
//		{
//			invalid_feature_ids.push_back(feature.id);
//			continue;
//		}
//
//		pointNum++;
//
//		if ((!feature.is_initialized) || ((feature.is_initialized) && (feature.is_initialized_NonKey)))
//		{
//			// TODO not only used in mono lk
//			// uniform velocity change of the feature points to remove outliers
//			//if (!stereo && !feature.checkMotion(cam_states))
//			if (!feature.checkMotion(cam_states))
//			{
//				invalid_feature_ids.push_back(feature.id);
//				cout << "checkmotionfail:" << endl;
//				pointInitialfail++;
//				continue;
//
//
//			}
//			if (!feature.triangulatePoint(cam_states))
//			{
//				invalid_feature_ids.push_back(feature.id);
//				cout << "triangulatefail:" << endl;
//
//				pointInitialfail++;
//				continue;
//			}
//			else
//			{
//				//if (!feature.LMOptimizatePositionAndCheck(cam_states))
//				//{
//				//	invalid_feature_ids.push_back(feature.id);
//				//	//cout << "initposfail:" << endl;
//				//	pointInitialfail++;
//				//	continue;
//				//}
//				if (!feature.initializePosition(cam_states))
//				{
//					invalid_feature_ids.push_back(feature.id);
//					//cout << "initposfail:" << endl;
//					pointInitialfail++;
//					continue;
//				}
//				else
//					feature.is_initialized_NonKey = false;
//			}
//		}
//
//
//		int observationSize = 0;
//
//		auto it = feature.observations.begin();
//		while (it != feature.observations.end())
//		{
//
//			if (cam_states[it->first].isKeyFrame)
//				observationSize++;
//
//			it++;
//		}
//		// this situation only appears in keyframe strategy
//		if (!stereo && observationSize <= 2)
//		{
//			invalid_feature_ids.push_back(feature.id);
//			// 深度恢复失败的点
//			pointInitialfail++;
//			continue;
//		}
//		if (stereo) jacobian_row_size += 4 * observationSize - 3;
//		else jacobian_row_size += 2 * observationSize - 3;
//		processed_feature_ids.push_back(feature.id);
//	}
//
//	for (const auto& feature_id : invalid_feature_ids)
//		map_server.erase(feature_id);
//
//	cout << "processed_feature_ids.size():" << processed_feature_ids.size() << endl;
//
//	if (processed_feature_ids.size() == 0)
//		return;
//
//
//	Eigen::MatrixXd H_x = Eigen::MatrixXd::Zero(jacobian_row_size,
//		15 + 6 * cam_states.size());
//	//Eigen::MatrixXd H_x = Eigen::MatrixXd::Zero(jacobian_row_size,
//	//	6 * cam_states.size());
//	Eigen::VectorXd r = Eigen::VectorXd::Zero(jacobian_row_size);
//	int stack_cntr = 0;
//
//
//	for (const auto& feature_id : processed_feature_ids)
//	{
//		auto& feature = map_server[feature_id];
//
//		vector<StateIDType> cam_state_ids(0);
//
//		for (const auto& measurement : feature.observations)
//		{
//			if (cam_states[measurement.first].isKeyFrame)
//				cam_state_ids.push_back(measurement.first);
//		}
//
//
//		Eigen::MatrixXd H_xj;
//		Eigen::VectorXd r_j;
//
//		test_featureJacobian(feature.id, cam_state_ids, H_xj, r_j);
//
//		if (test_gatingTest(H_xj, r_j, cam_state_ids.size() - 1))
//		{
//			H_x.block(stack_cntr, 0, H_xj.rows(), H_xj.cols()) = H_xj;
//			r.segment(stack_cntr, r_j.rows()) = r_j;
//			stack_cntr += H_xj.rows();
//		}
//		else
//		{
//			// 通不过卡方检测的点数
//			cout << "not_gatingtest:" << endl;
//			pointgatingTestfail++;
//		}
//		if (stack_cntr > 1500)
//			break;
//	}
//
//	//实际更新点数
//	pointUpdate = pointNum - pointInitialfail - pointgatingTestfail;
//	cout << "remove_jacobian_size:" << stack_cntr << endl;
//	H_x.conservativeResize(stack_cntr, H_x.cols());
//	r.conservativeResize(stack_cntr);
//
//	test_cam_meas_update(H_x, r);
//
//
//	for (const auto& feature_id : processed_feature_ids)
//	{
//		map_server.erase(feature_id);
//	}
//
//	cout << "!!!remove_finshed!!!" << endl;
//	return;
//}
//
//bool gsins::t_gmultign::test_gatingTest(const Eigen::MatrixXd& H, const Eigen::VectorXd& r, const int& dof)
//{
//	Eigen::MatrixXd state_cov = test_get_statecov();
//	//Eigen::MatrixXd state_cov = BaseMatrix2Eigen(_Qx);
//	//Eigen::MatrixXd state_cov = get_statecov();
//	Eigen::MatrixXd P1 = H * state_cov*H.transpose();
//	//double tmp_observation_noise = (5.0 / 1000) * (5.0 / 1000) ;
//	//MatrixXd P1 = H * state_server.state_cov * H.transpose();
//	Eigen::MatrixXd P2 = observation_noise * Eigen::MatrixXd::Identity(H.rows(), H.rows());
//
//	// 残差加权平均数
//	double gamma = r.transpose() * (P1 + P2).ldlt().solve(r);
//
//	//cout << r << endl;
//	//cout << dof << " " << gamma << " " << chi_squared_test_table[dof] << " " << endl;
//	if (gamma < chi_squared_test_table[dof])
//	{
//
//		return true;
//	}
//	else
//	{
//
//		return false;
//	}
//}
//
//
//void gsins::t_gmultign::test_cam_meas_update(const Eigen::MatrixXd& H, const Eigen::VectorXd& r)
//{
//	if (H.rows() == 0 || r.rows() == 0)
//		return;
//
//	Eigen::MatrixXd H_thin;
//	Eigen::VectorXd r_thin;
//
//
//	//lsy
//	//cout << "r:" << r.transpose() << endl;
//	if (H.rows() > H.cols())
//	{
//
//		Eigen::HouseholderQR<Eigen::MatrixXd> qr_helper(H);
//		Eigen::MatrixXd Q = qr_helper.householderQ();
//		Eigen::MatrixXd Q1 = Q.leftCols(15 + cam_states.size() * 6);
//
//		H_thin = Q1.transpose() * H;
//		r_thin = Q1.transpose() * r;
//	}
//	else
//	{
//
//		H_thin = H;
//		r_thin = r;
//	}
//
//	Eigen::MatrixXd state_cov = test_get_statecov();
//	const Eigen::MatrixXd& P = state_cov;
//	Eigen::MatrixXd S = H_thin * P*H_thin.transpose() +
//		observation_noise * Eigen::MatrixXd::Identity(
//			H_thin.rows(), H_thin.rows());
//
//	Eigen::MatrixXd K_transpose = S.ldlt().solve(H_thin*P);
//	Eigen::MatrixXd K = K_transpose.transpose();
//
//	Eigen::VectorXd delta_x = K * r_thin;
//	//cout << "delta_x_rows" << delta_x.rows() << endl;
//
//	const Eigen::VectorXd& delta_x_imu = delta_x.head<15>();
//
//	if (delta_x_imu.segment<3>(6).norm() > 0.5 ||
//		delta_x_imu.segment<3>(12).norm() > 1.0)
//	{
//		printf("delta velocity: %f\n", delta_x_imu.segment<3>(6).norm());
//		printf("delta position: %f\n", delta_x_imu.segment<3>(12).norm());
//	}
//	//sins.qnb = sins.qnb - delta_x_imu.head<3>();
//	//sins.vn = sins.vn - delta_x_imu.segment<3>(3);
//	//sins.pos = sins.pos - delta_x_imu.segment<3>(6);
//	//sins.eb = sins.eb + delta_x_imu.segment<3>(9);
//	//sins.db = sins.db + delta_x_imu.segment<3>(12);
//
//	sins.qeb = sins.qeb - delta_x_imu.head<3>();
//	sins.ve = sins.ve - delta_x_imu.segment<3>(3);
//	sins.pos_ecef = sins.pos_ecef - delta_x_imu.segment<3>(6);
//	sins.eb = sins.eb + delta_x_imu.segment<3>(9);
//	sins.db = sins.db + delta_x_imu.segment<3>(12);
//
//	sins.qnb = t_gbase::m2qua(sins.eth.Cne) * sins.qeb;
//	sins.vn = sins.eth.Cne*sins.ve;
//	sins.pos = Cart2Geod(sins.pos_ecef, false);
//	sins.att = t_gbase::q2att(sins.qnb);
//	sins.Ceb = t_gbase::q2mat(sins.qeb);
//	sins.Cnb = t_gbase::q2mat(sins.qnb);
//	cout << "CAM_MEAS_UP:" << endl;
//	cout << "imu_dq:" << setprecision(12) << delta_x_imu.head<3>().transpose() << endl;
//	cout << "imu_dpos_ecef:" << setprecision(12) << delta_x_imu.segment<3>(6).transpose() << endl;
//
//	auto cam_state_iter = cam_states.begin();
//	for (int i = 0; i < cam_states.size();
//		++i, ++cam_state_iter)
//	{
//
//		const Eigen::VectorXd& delta_x_cam = delta_x.segment<6>(15 + i * 6);
//		t_gquat dq_cam = t_gbase::rv2q(-delta_x_cam.head<3>());
//		Eigen::Quaterniond _dq_cam;
//		_dq_cam.w() = dq_cam.q0; _dq_cam.x() = dq_cam.q1; _dq_cam.y() = dq_cam.q2; _dq_cam.z() = dq_cam.q3;
//		_dq_cam.normalize();
//
//		cam_state_iter->second.orientation = cam_state_iter->second.orientation * _dq_cam;
//		cam_state_iter->second.position -= delta_x_cam.tail<3>();
//
//
//	}
//
//
//	Eigen::MatrixXd I_KH = Eigen::MatrixXd::Identity(K.rows(), H_thin.cols()) - K * H_thin;
//
//	state_cov = I_KH * state_cov*I_KH.transpose() +
//		K * K.transpose()*observation_noise;
//	Eigen::MatrixXd state_cov_fixed = (state_cov + state_cov.transpose()) / 2.0;
//
//	test_assign_statecov(state_cov_fixed);
//}
//
//void gsins::t_gmultign::test_pruneCamStateBuffer()
//{
//	if (KeyFrameBuffer.size() < max_camstate_size)
//		return;
//
//	int pointNum = 0;
//	int pointInitialfail = 0;
//	int pointgatingTestfail = 0;
//	int pointUpdate = 0;
//
//
//	vector<StateIDType> rm_cam_state_ids(0);
//
//	findRedundantCamStates(rm_cam_state_ids);
//
//
//	int jacobian_row_size = 0;
//	for (auto& item : map_server)
//	{
//		auto& feature = item.second;
//
//
//		vector<StateIDType> involved_cam_state_ids(0);
//		for (const auto& cam_id : rm_cam_state_ids)
//		{
//			if (feature.observations.find(cam_id) !=
//				feature.observations.end())
//				involved_cam_state_ids.push_back(cam_id);
//		}
//
//		if (involved_cam_state_ids.size() == 0)
//			continue;
//		if (involved_cam_state_ids.size() == 1)
//		{
//			feature.observations.erase(involved_cam_state_ids[0]);
//			continue;
//		}
//		if ((involved_cam_state_ids.size() == 2) && (!feature.is_KeyFrame))
//		{
//			feature.observations.erase(involved_cam_state_ids[0]);
//			feature.observations.erase(involved_cam_state_ids[1]);
//			continue;
//		}
//
//
//		pointNum++;
//
//		if ((!feature.is_initialized) || ((feature.is_initialized) && (feature.is_initialized_NonKey)))
//		{
//			if (!stereo && !feature.checkMotion(cam_states))
//			{
//
//				for (const auto& cam_id : involved_cam_state_ids)
//					feature.observations.erase(cam_id);
//				cout << "checkmotionfaile:" << endl;
//				pointInitialfail++;
//				continue;
//
//			}
//			if (!feature.triangulatePoint(cam_states))
//			{
//
//				for (const auto& cam_id : involved_cam_state_ids)
//					feature.observations.erase(cam_id);
//				cout << "triangulatefail:" << endl;
//				pointInitialfail++;
//				continue;
//			}
//			else
//			{
//				//if (!feature.LMOptimizatePositionAndCheck(cam_states))
//				//{
//				//	for (const auto& cam_id : involved_cam_state_ids)
//				//		feature.observations.erase(cam_id);
//				//	//cout << "initposfail:" << endl;
//				//	pointInitialfail++;
//
//				//	continue;
//				//}
//				if (!feature.initializePosition(cam_states))
//				{
//
//					for (const auto& cam_id : involved_cam_state_ids)
//						feature.observations.erase(cam_id);
//					//cout << "initposfail:" << endl;
//					pointInitialfail++;
//
//					continue;
//				}
//				else
//					feature.is_initialized_NonKey = false;
//			}
//		}
//
//		if (stereo) jacobian_row_size += 4 * involved_cam_state_ids.size() - 3;
//		else	jacobian_row_size += 2 * involved_cam_state_ids.size() - 3;
//	}
//
//	Eigen::MatrixXd H_x = Eigen::MatrixXd::Zero(jacobian_row_size,
//		15 + 6 * cam_states.size());
//	Eigen::VectorXd r = Eigen::VectorXd::Zero(jacobian_row_size);
//
//	int stack_cntr = 0;
//	for (auto& item : map_server)
//	{
//		auto& feature = item.second;
//
//
//		vector<StateIDType> involved_cam_state_ids(0);
//		for (const auto& cam_id : rm_cam_state_ids)
//		{
//			if (feature.observations.find(cam_id) !=
//				feature.observations.end())
//				involved_cam_state_ids.push_back(cam_id);
//		}
//
//		if (involved_cam_state_ids.size() == 0)
//			continue;
//
//
//		Eigen::MatrixXd H_xj;
//		Eigen::VectorXd r_j;
//		test_featureJacobian(feature.id, involved_cam_state_ids, H_xj, r_j);
//
//
//		if (test_gatingTest(H_xj, r_j, involved_cam_state_ids.size()))
//		{
//			H_x.block(stack_cntr, 0, H_xj.rows(), H_xj.cols()) = H_xj;
//			r.segment(stack_cntr, r_j.rows()) = r_j;
//			stack_cntr += H_xj.rows();
//		}
//		else
//		{
//			cout << "not_gatingtest:" << endl;
//			pointgatingTestfail++;
//		}
//		for (const auto& cam_id : involved_cam_state_ids)
//			feature.observations.erase(cam_id);
//	}
//	H_x.conservativeResize(stack_cntr, H_x.cols());
//	r.conservativeResize(stack_cntr);
//	cout << "prune_jacobian_size:" << stack_cntr << endl;
//	//实际更新点数
//	pointUpdate = pointNum - pointInitialfail - pointgatingTestfail;
//
//
//	test_cam_meas_update(H_x, r);
//
//
//
//	for (const auto& cam_id : rm_cam_state_ids)
//	{
//
//		int cam_sequence = std::distance(cam_states.begin(),
//			cam_states.find(cam_id));
//		int cam_state_start = 15 + 6 * cam_sequence;
//		int cam_state_end = cam_state_start + 6;
//
//		//协方差消除 做好
//		Eigen::MatrixXd state_cov = test_get_statecov();
//		if (cam_state_end < state_cov.rows())
//		{
//			state_cov.block(cam_state_start, 0,
//				state_cov.rows() - cam_state_end,
//				state_cov.cols()) =
//				state_cov.block(cam_state_end, 0,
//					state_cov.rows() - cam_state_end,
//					state_cov.cols());
//
//
//			state_cov.block(0, cam_state_start,
//				state_cov.rows(),
//				state_cov.cols() - cam_state_end) =
//				state_cov.block(0, cam_state_end,
//					state_cov.rows(),
//					state_cov.cols() - cam_state_end);
//
//			state_cov.conservativeResize(
//				state_cov.rows() - 6, state_cov.cols() - 6);
//		}
//		else
//		{
//			state_cov.conservativeResize(
//				state_cov.rows() - 6, state_cov.cols() - 6);
//		}
//
//
//		cam_states.erase(cam_id);
//
//
//
//		//	//消除后 重新赋值
//		P_CC.conservativeResize(6 * cam_states.size(), 6 * cam_states.size());
//
//		P_IC.conservativeResize(15, 6 * cam_states.size());
//		Eigen::MatrixXd state_cov_fixed = (state_cov +
//			state_cov.transpose()) / 2.0;
//
//		test_assign_statecov(state_cov_fixed);
//		//assign_statecov(state_cov——);
//
//
//	}
//
//
//	return;
//}

#endif // USE_VIS