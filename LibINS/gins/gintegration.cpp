#include "gintegration.h"
#include "gproc/gpreproc.h"
#include "gutils/gtimesync.h"
#include <iostream>
using namespace gsins;
using namespace std;


/******************************************* t_gintegration ***********************************************/
gsins::t_gintegration::t_gintegration(string site, string site_base, t_gsetbase* gset, t_gallproc *allproc) :
	t_gsinskf(gset,&(/*site_base+*/site)),
	t_gspp(site, gset),
	t_gpvtflt(site, site_base, gset, allproc)
{
	kftk = 0.0;
	nq = dynamic_cast<t_gsetign*>(gset)->nq();
	nr = dynamic_cast<t_gsetign*>(gset)->nr();
	Ft = Pk = Eigen::MatrixXd::Zero(nq, nq);
	Hk = Eigen::MatrixXd::Zero(nr, nq);
	Qt = Pmin = Xk = Eigen::VectorXd::Zero(nq);  Pmax = Eigen::VectorXd::Ones(nq)* glv.INF;
	Rt = Zk = Eigen::VectorXd::Zero(nr);  rts = Eigen::VectorXd::Ones(nr);
	Rk = Eigen::MatrixXd::Zero(nr, nr);
	Rmax = Eigen::VectorXd::Ones(nr)* glv.INF; Rmin = Rb = Eigen::VectorXd::Zero(nr); Rbeta = Eigen::VectorXd::Ones(nr);
	FBTau = FBMax = Eigen::VectorXd::Ones(nq)* glv.INF; FBXk = FBTotal = Eigen::VectorXd::Zero(nq);
	Phik = Eigen::MatrixXd::Identity(nq, nq);
	measflag = 0;
	FilterFlag = dynamic_cast<t_gsetign*>(gset)->fltmode();
	_gnss_beg = dynamic_cast<t_gsetgen*>(gset)->beg();
	_gnss_end = dynamic_cast<t_gsetgen*>(gset)->end();
	double start = dynamic_cast<t_gsetins*>(gset)->start();
	double end = dynamic_cast<t_gsetins*>(gset)->end();
	_ins_beg = t_gtime(_gnss_beg.gwk(), start);
	_ins_end = t_gtime(_gnss_end.gwk(), end);
	_ign_type = dynamic_cast<t_gsetign*>(gset)->ign_type();
	_imudata = nullptr;
	_initial_merge = false;

	string port = dynamic_cast<t_gsetout*>(_setkf)->output_port("ins");
	if (port.size() > 0)
	{
		t_gtcp* ins = new t_gtcp(port);
		_maptcp.insert(make_pair(INS_OUT, ins));
	}
}

void gsins::t_gintegration::Add_IMU(t_gimudata * imu)
{
	_imudata = imu;
}

/************************************************************
** Function: Integration Continuous time measurement vector coefficient matrix
***********************************************************/
void gsins::t_gintegration::set_Hk()
{
	int col = _param.parNumber();
	Eigen::Vector3d vartheta = sins.Ceb*t_gbase::askew(lever)*sins.wib + t_gbase::askew(sins.eth.weie)*sins.Ceb*lever;
	switch (Flag)
	{
	case POS_VEL_MEAS:
		Hk.block(0, 3, 3, 3) = Eigen::Matrix3d::Identity();      // vel-vel
		Hk.block(0, 0, 3, 3) = -t_gbase::askew(vartheta);		 // vel-att
		Hk.block(0, 9, 3, 3) = -sins.Ceb*t_gbase::askew(lever);  // vel-eb
		Hk.block(3, 6, 3, 3) = Eigen::Matrix3d::Identity();      // pos-pos
		Hk.block(3, 0, 3, 3) = t_gbase::askew(sins.Ceb*lever);   // pos-att
		break;
	case POS_MEAS:
		Hk.block(0, 6, 3, 3) = Eigen::Matrix3d::Identity();
		Hk.block(0, 0, 3, 3) = t_gbase::askew(sins.Ceb*lever);	 // pos-att
		break;
	case VEL_MEAS:
		Hk.block(0, 3, 3, 3) = Eigen::Matrix3d::Identity();      // vel-vel
		Hk.block(0, 0, 3, 3) = -t_gbase::askew(vartheta);		 // vel-att
		Hk.block(0, 9, 3, 3) = -sins.Ceb*t_gbase::askew(lever);  // vel-eb
		break;
	case FIXED_DD_UC_OBS_MEAS:  // independent-amb for rtk 
		Hk.block(0, 6, nr, 3) = -_gnss_design_mat.block(0, 0, nr, 3);   // obs-pos
		//Hk.block(0, 0, nr, 3) = -_gnss_design_mat.block(0, 0, nr, 3)*t_gbase::askew(sins.Ceb*lever);	 // obs-att
		break;
	case ZUPT_MEAS:
		Hk.block(0, 3, 3, 3) = Eigen::Matrix3d::Identity();
		break;
	case NHC_MEAS:
		Hk.block(0, 3, 1, 3) = sins.Cbe.block(0, 0, 1, 3);
		Hk.block(1, 3, 1, 3) = sins.Cbe.block(2, 0, 1, 3);
		break;
	case NO_MEAS:
		Hk = Eigen::MatrixXd::Zero(nr, nq);
		break;
	case ATT_MEAS:
		Hk.block(0, 0, 3, 3) = Eigen::Matrix3d::Identity();
	case ZUPT_POS_MEAS:
		Hk.block(0, 6, 3, 3) = Eigen::Matrix3d::Identity();
	default:
		break;
	}
}

int gsins::t_gintegration::processBatch(const t_gtime & beg, const t_gtime & end, bool prtOut)
{
	gtrace("t_gintegration::processBatch");
#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif 
	_gmutex.lock();

	if (_grec == nullptr) {
		ostringstream os;
		os << "ERROR: No object found (" << _site << "). Processing terminated!!! " << beg.str_ymdhms() << " -> " << end.str_ymdhms() << endl;
		if (_log) _log->comment(0, "gintegration", os.str());
		_gmutex.unlock(); return -1;
	}

	int fltflag = 1; if (!_beg_end) fltflag = -1;
	t_gpvtflt::InitProc(beg, end);

	if (!this->_init()) {
		if (_log) _log->comment(0, "gintegration", "init failed");
		_gmutex.unlock(); return -1;
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
				int nEpo = round(_ins_crt.diff(_gnss_crt) / (fltflag * _sampling) + 0.5);
				if (_sampling > 1) _gnss_crt.add_secs(int(fltflag * _sampling * nEpo));  // =<1Hz data
				else               _gnss_crt.add_dsec(fltflag * _sampling * nEpo);       //  >1Hz data
			}
			while (_gnss_crt < _ins_crt)
			{
				if (_sampling > 1) _gnss_crt.add_secs(int(fltflag * _sampling));  // =<1Hz data
				else               _gnss_crt.add_dsec(fltflag * _sampling);       //  >1Hz data
			}
		}
	}

	std::cerr << "\n" << _site << ": Start GNSS/INS Integration Processing: " << _ins_beg.str_ymdhms() << " " << _ins_end.str_ymdhms() << endl;
	cout << "EPOCH " << _ins_crt.str_ymdhms() << endl;
	t_gposdata::data_pos posdata; posdata = 0.0;
	bool time_loop = true;
	while (time_loop)
	{
		if (_ins_crt >= _ins_end || _gnss_crt >= _gnss_end)
		{
			cerr << " Processing Finished!" << endl;
			if (_log) _log->comment(0, "gintegration", "Processing Finished!");
			time_loop = false; break;
		}

		if (realtime) {
			t_gtime crt = t_gtime::current_time(t_gtime::GPS); this->_avaliable(_ins_crt);
		}
		if (!_imudata->load(_wm, _vm, _shm.t, _shm.ts, _shm.nSamples)) { cerr << "\n Processing Finished!" << endl; break; }

		bool imu_gap = false;
		double inflation = 1.0;

		if (fabs(_shm.ts*dynamic_cast<t_gsetins*>(_set)->freq()) >= 2.0)
		{
			imu_gap = true;
			inflation = pow(_shm.ts*dynamic_cast<t_gsetins*>(_set)->freq(), 6);
			std::cerr << "Find IMU Gap! :  " << _ins_crt.str_hms() << " " << setprecision(5) << setiosflags(ios::fixed) << _shm.ts << std::endl;
		}

		_ins_crt = t_gtime(_gnss_crt.gwk(), _shm.t);
		if (!_aligned)
		{
			if (_time_valid(_gnss_crt, _ins_crt))Flag = _getPOS(posdata);
			_aligned = cascaded_align(posdata.pos, posdata.vn);
		}
		else
		{
			if (_ign_type == IGN_TYPE::TCI && !_initial_merge) _merge_init();
			sins.Update(_wm, _vm, _shm);
			if (_ign_type == PURE_INS) {
				cerr << "\r" << _ins_crt.str_ymdhms();
				t_gsinskf::write_tmp(); continue;
			}

			time_update(sins.nts, 0, inflation);

			if (_time_valid(_gnss_crt, _ins_crt))
			{
				if (_ign_type == IGN_TYPE::LCI)
				{
					Flag = _getPOS(posdata);
					if (Flag != NO_MEAS) _meas_update();
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
				else if (_ign_type == IGN_TYPE::TCI)
				{
					_timeUpdate(_gnss_crt);
					Flag = MEAS_TYPE(this->_processEpoch(_gnss_crt));
					if (Flag == NO_MEAS) { _reset(); }
				}

				if (Flag != NO_MEAS) {
					feedback(glv.INF);
					sins.pure_ins_time = 0.0;
					cerr << "\r" << _ins_crt.str_ymdhms();

					if (_log) _log->comment(0, "gintegration",
						_ins_crt.str_ymdhms(" epoch ") + "Processing GNSS/INS Integration Success!");
				}
			}

			if (!(_isBase && _ign_type == IGN_TYPE::TCI)) {
				t_gsinskf::write_tmp();
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

int gsins::t_gintegration::processBatch(const t_gtime & beg, const t_gtime & end)
{
	gtrace("t_gintegration::processBatch");
#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif 
	//_gmutex.lock();

	if (_grec == nullptr) {
		ostringstream os;
		os << "ERROR: No object found (" << _site << "). Processing terminated!!! " << beg.str_ymdhms() << " -> " << end.str_ymdhms() << endl;
		if (_log) _log->comment(0, "gintegration", os.str());
		//	_gmutex.unlock(); 
		return -1;
	}

	int fltflag = 1; if (!_beg_end) fltflag = -1;
	t_gpvtflt::InitProc(beg, end);

	if (!this->_init()) {
		if (_log) _log->comment(0, "gintegration", "init failed");
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
				int nEpo = round(_ins_crt.diff(_gnss_crt) / (fltflag * _sampling) + 0.5);
				if (_sampling > 1) _gnss_crt.add_secs(int(fltflag * _sampling * nEpo));  // =<1Hz data
				else               _gnss_crt.add_dsec(fltflag * _sampling * nEpo);       //  >1Hz data
			}
			//while (_gnss_crt < _ins_crt)
			//{
			//	if (_sampling > 1) _gnss_crt.add_secs(int(fltflag * _sampling));  // =<1Hz data
			//	else               _gnss_crt.add_dsec(fltflag * _sampling);       //  >1Hz data
			//}
		}
	}

	std::cerr << "\n" << _site << ": Start GNSS/INS Integration Processing: " << _ins_beg.str_ymdhms() << " " << _ins_end.str_ymdhms() << endl;
	cout << "EPOCH " << _ins_crt.str_ymdhms() << endl;
	t_gposdata::data_pos posdata;
	bool time_loop = true;
	// ins time loop
	while (time_loop)
	{
		if (_ins_crt > _ins_end || _gnss_crt > _gnss_end)
		{
			cerr << " Processing Finished!" << endl;
			if (_log) _log->comment(0, "gintegration", "Processing Finished!");
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
			if (_ign_type == IGN_TYPE::TCI && !_initial_merge) _merge_init();
			sins.Update(_wm, _vm, _shm);
			if (_ign_type == PURE_INS) {
				cerr << "\r" << _ins_crt.str_ymdhms();
				t_gsinskf::write_tmp(); continue;
			}

			time_update(sins.nts, 0.0);
			int imeas = _getMeas(), irc = -1;
			set<MEAS_TYPE>::const_iterator it = _Meas_Type.begin();
			while (it != _Meas_Type.end())
			{
				switch (*it)
				{
				case GNSS_MEAS:
				{
					if (irc = _GNSS_Update() > 0)t_gsinskf::_restore(Xk, Pk);
					else t_gsinskf::_recover(Xk, Pk);
				}break;
				case ZUPT_MEAS:
				{
					if (irc = _ZUPT_Update() > 0) t_gsinskf::_restore(Xk, Pk);
					else t_gsinskf::_recover(Xk, Pk);
				}break;
				case NHC_MEAS:
				{
					if (irc = _NHC_Update() > 0) t_gsinskf::_restore(Xk, Pk);
					else t_gsinskf::_recover(Xk, Pk);
				} break;
				default:break;
				}it++;
			}

			if (Flag != NO_MEAS) {
				feedback(glv.INF); _prt_port();
				sins.pure_ins_time = 0.0;
				cerr << "\r" << _ins_crt.str_ymdhms() << setw(5) << " Q = " << (_amb_state ? 1 : 2);
				if (_log) _log->comment(0, "gintegration",
					_ins_crt.str_ymdhms(" epoch ") + "Processing GNSS/INS Integration Success!");
			}

			if (!(_isBase && _ign_type == IGN_TYPE::TCI)) {
				t_gsinskf::write_tmp();
			}
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

int gsins::t_gintegration::processBatch2(const t_gtime & beg, const t_gtime & end, bool prtOut,bool beg_end)
{
	gtrace("t_gintegration::processBatch2");
#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif 

	if (_grec == nullptr) {
		ostringstream os;
		os << "ERROR: No object found (" << _site << "). Processing terminated!!! " << beg.str_ymdhms() << " -> " << end.str_ymdhms() << endl;
		if (_log) _log->comment(0, "gintegration", os.str());
		//_gmutex.unlock(); 
		return -1;
	}

	fltflag = 1;
	string flag = "(forward direction) ";
	string symbol="¨€¨€¨€¨€¨€¨€¨€¨€¨€¨€";
	_beg_end = beg_end;

	t_gpvtflt::InitProc(beg, end);

	if (!this->_init()) {
		if (_log) _log->comment(0, "gintegration", "init failed" + flag);
		//_gmutex.unlock();
		return -1;
	}

	if (!_beg_end)
	{
		flag = "(backward direction)";
		fltflag = -1;
		sins.eth.wie = -sins.eth.wie;
	}
	// need execute frequently in real-time mode
	_gobs->setepoches(_site);
	bool realtime = dynamic_cast<t_gsetproc*>(_set)->realtime();
	// time align
	if (!realtime) {
		if ((_beg_end && _ins_beg < _gnss_beg)|| (!_beg_end && _ins_beg > _gnss_beg)) {
			_ins_crt = _imudata->erase_bef(_gnss_beg, _beg_end);
		}
		else
		{
			//_gmutex.unlock();
			if (_beg_end) t_gpvtflt::processBatch(_gnss_crt, _ins_crt, false);
			else t_gpvtflt::processBatch(_ins_crt, _gnss_crt, false);
			//_gmutex.lock();
			// next gnss time 
			if ((_beg_end && _gnss_crt < _ins_crt)||(!_beg_end &&_gnss_crt> _ins_crt)){
				int nEpo = round(_ins_crt.diff(_gnss_crt) / (fltflag * _sampling) + 0.5);
				if (_sampling > 1) _gnss_crt.add_secs(int(fltflag * _sampling * nEpo));  //  < 1Hz data
				else               _gnss_crt.add_dsec(fltflag * _sampling * nEpo);       //  >=1Hz data
			}
		}
	}

	std::cerr << "\n" << _site << ": Start GNSS/INS Integration Processing:"<<flag << _ins_beg.str_ymdhms() << " " << _ins_end.str_ymdhms() << endl;
	cout << "EPOCH " << _ins_crt.str_ymdhms() << endl;
	t_gposdata::data_pos posdata;
	bool time_loop = true;
	// ins time loop

	while (time_loop)
	{
		if ((_beg_end && (_ins_crt > _ins_end || _gnss_crt > _gnss_end)) || (!_beg_end && (_ins_crt < _ins_end || _gnss_crt < _gnss_end)))
		{
			//cerr << " Processing Finished!" << flag << endl;
			if (_log) _log->comment(0, "gintegration", "Processing Finished!" + flag);
			time_loop = false; break;
		}
		if (realtime) {
			t_gtime crt = t_gtime::current_time(t_gtime::GPS); this->_avaliable(_ins_crt);
		}
		
		if (!_imudata->load(_wm, _vm, _shm.t, _shm.ts, _shm.nSamples,_beg_end)) 
		{ 
			//cerr << "\n Processing Finished!" << flag << endl; 
			break; 
		}
		_ins_crt = t_gtime(_gnss_crt.gwk(), _shm.t);
		if (!_aligned)
		{
			if (_time_valid(_gnss_crt, _ins_crt)) Flag = _getPOS(posdata);
			_aligned = cascaded_align(posdata.pos, posdata.vn);
		}
		else
		{
			if (_ign_type == IGN_TYPE::TCI && !_initial_merge) _merge_init();
			sins.Update(_wm, _vm, _shm);
			if (_ign_type == PURE_INS) {
				cerr << "\r" << "forward direction: " << left << setw(20) << symbol.substr(0, 2*int(*_fwd_per * 10))
					<< " backward direction: " << right << setw(20) << symbol.substr(0, 2 * int(*_bwd_per * 10));
				if(prtOut) t_gsinskf::write_tmp(_beg_end);
				continue;
			}

			time_update(sins.nts, 0.0); 

			//sins.int_sec();
			kftk = sins.t;
			int imeas = _getMeas(), irc = -1;
			set<MEAS_TYPE>::const_iterator it = _Meas_Type.begin();
			while (it != _Meas_Type.end())
			{
				switch (*it)
				{
				case GNSS_MEAS:
				{
					if (irc = _GNSS_Update() > 0)t_gsinskf::_restore(Xk, Pk);
					else t_gsinskf::_recover(Xk, Pk);
				}break;
				case ZUPT_MEAS:
				{
					if (irc = _ZUPT_Update() > 0) t_gsinskf::_restore(Xk, Pk);
					else t_gsinskf::_recover(Xk, Pk);
				}break;
				case NHC_MEAS:
				{
					if (irc = _NHC_Update() > 0) t_gsinskf::_restore(Xk, Pk);
					else t_gsinskf::_recover(Xk, Pk);
				} break;
				default:break;
				}it++;
			}

			if (_beg_end)
			{
				if (!prtOut)
				{
					Xk_forward.insert(make_pair(kftk, Xk));
					Pk_forward.insert(make_pair(kftk, Pk));
					sins_forward.insert(make_pair(kftk, sins));
				}
				*_fwd_crt = _ins_crt;
				*_fwd_per = (*_fwd_crt).diff(_ins_beg) / _ins_end.diff(_ins_beg);
			}
			else
			{
				if (!prtOut) 
				{
					Xk_backward.insert(make_pair(kftk, Xk));
					Pk_backward.insert(make_pair(kftk, Pk));
					sins_backward.insert(make_pair(kftk, sins));
				}			
				*_bwd_crt = _ins_crt;
				*_bwd_per = _ins_beg.diff(*_bwd_crt) / _ins_beg.diff(_ins_end);
			}

			if (Flag != NO_MEAS) {
				feedback(glv.INF); _prt_port();
				sins.pure_ins_time = 0.0;
				/*cerr << "\r" << "forward direction: "<< (*_fwd_crt).str_ymdhms() 
					<< " backward direction: " << (*_bwd_crt).str_ymdhms();*/
				cerr << "\r" << "forward direction: " << left << setw(20) << symbol.substr(0, 2 * int(*_fwd_per * 10))
					<< " backward direction: " << right << setw(20) << symbol.substr(0, 2 * int(*_bwd_per * 10));
				if (_log) _log->comment(0, "gintegration",
					_ins_crt.str_ymdhms(" epoch ") + "Processing GNSS/INS Integration Success!" + flag);
			}

			if (!(_isBase && _ign_type == IGN_TYPE::TCI)&& prtOut) {
				t_gsinskf::write_tmp(_beg_end);
			}
		}
		while (abs(_ins_crt.diff(_gnss_crt)) > _sampling)
		{
			if (_sampling > 1) _gnss_crt.add_secs(fltflag*int(_sampling));  // =<1Hz data
			else               _gnss_crt.add_dsec(fltflag*_sampling);       //  >1Hz data
		}

	}

	//_gmutex.unlock();
	return 0;
}

void gsins::t_gintegration::feedback(double kfts)
{
	gtrace("t_gintegration::feedback");
	try
	{
		t_gkalman::feedback(kfts);
	}
	catch (out_of_range)
	{
		cout << "t_gintegration::feedback" + to_string(kftk) + ": " + "FeedBack Failed!" << endl;
	}
	//cout << FBXk << endl;
	sins.qeb = sins.qeb - FBXk.block(0, 0, 3, 1);
	sins.ve = sins.ve - FBXk.block(3, 0, 3, 1);
	sins.pos_ecef = sins.pos_ecef - FBXk.block(6, 0, 3, 1);
	sins.eb = sins.eb + FBXk.block(9, 0, 3, 1);
	sins.db = sins.db + FBXk.block(12, 0, 3, 1);

	sins.qnb = t_gbase::m2qua(sins.eth.Cne) * sins.qeb;
	sins.vn = sins.eth.Cne*sins.ve;
	sins.pos = Cart2Geod(sins.pos_ecef, false);
	sins.att = t_gbase::q2att(sins.qnb);
	sins.Ceb = t_gbase::q2mat(sins.qeb);
	sins.Cnb = t_gbase::q2mat(sins.qnb);
}

int gsins::t_gintegration::_gnss_feedback(const ColumnVector & dx)
{
	gtrace("t_gintegration::_gnss_feedback");

	try {
		// parameter dx
		int ign_size = _param.parNumber();
		int ins_size = dynamic_cast<t_gsetign*>(_set)->nq();

		Xk = Columns2VectorXd(dx).block(0, 0, ins_size, 1);
		Pk = BaseMatrix2Eigen(_Qx).block(0, 0, ins_size, ins_size);


		// trophere
		for (unsigned int iPar = 0; iPar < _param.parNumber(); iPar++) {
			if (_param[iPar].parType == t_gpar::TRP)
			{
				string site = _param[iPar].site;
				t_gtriple Ell, XYZ;
				if (_param.getCrdParam(site, XYZ) > 0) {
				}
				else { XYZ = _gallobj->obj(site)->crd_arp(_epoch); }
				xyz2ell(XYZ, Ell, false);
				if (site == _site && _gModel->tropoModel() != 0)
					_param[iPar].apriori(_gModel->tropoModel()->getZHD(Ell, _epoch));
				else if (site == _site_base && _gModel_base->tropoModel() != 0)
					_param[iPar].apriori(_gModel_base->tropoModel()->getZHD(Ell, _epoch));
			}
		}

		// parameter transmit
		for (unsigned int iPar = 0; iPar < _param.parNumber(); iPar++) {
			if (_param[iPar].parType == t_gpar::CRD_X || _param[iPar].parType == t_gpar::CRD_Y || _param[iPar].parType == t_gpar::CRD_Z)
			{
				//cout << "EPO_INS        " << _param[iPar].str_type() << "         " << fixed << setw(12) << setprecision(3) << _param[iPar].value() << fixed << setw(12) << setprecision(3) << dx(_param[iPar].index) << fixed << setw(12) << setprecision(3) << _param[iPar].value() - dx(_param[iPar].index) << endl;
				_param[iPar].value(_param[iPar].value() - dx(_param[iPar].index));

			}
			else
			{

				//cout << "EPO_INS2        " << _param[iPar].str_type() << "         " << fixed << setw(12) << setprecision(3) << _param[iPar].value() << fixed << setw(12) << setprecision(3) << dx(_param[iPar].index) << fixed << setw(12) << setprecision(3) << _param[iPar].value() + dx(_param[iPar].index) << endl;
				_param[iPar].value(_param[iPar].value() + dx(_param[iPar].index));
			}
		}
		//ostringstream os;
		//_prtOut(_epoch, _param, _Qx, _data, os, line, true);
		//// Print flt results
		//if (_flt) {
		//	_flt->write(os.str().c_str(), os.str().size());
		//	_flt->flush();
		//}
		return 1;
	}
	catch (...)
	{
		if (_log) _log->comment(0, "gintegration", "GNSS feedback failed!");
		_reset();
		Flag = NO_MEAS;
		return -1;
	}
}

int gsins::t_gintegration::_reset()
{
	gtrace("t_gintegration::_reset");

	try
	{
		int ins_size = dynamic_cast<t_gsetign*>(_set)->nq();
		// reset Xk
		Xk.resize(ins_size); Xk = Eigen::VectorXd::Zero(ins_size);
		// reset Pk
		Pk = Eigen::MatrixXd::Zero(ins_size, ins_size); Pk = Pk_Sav;  // because of outlier loop
	}
	catch (...)
	{
		if (_log) _log->comment(0, "gintegration", "GNSS and INS filter Reset failed!");
		return -1;
	}

	if (_log) _log->comment(0, "gintegration", "GNSS and INS filter Reset.");
	return 1;
}

int gsins::t_gintegration::_prt_ins_kml()
{
	gtrace("t_gintegration::_prt_ins_kml");

	Eigen::Vector3d Geo_pos = Eigen::Vector3d(sins.pos(0) / glv.deg, sins.pos(1) / glv.deg, sins.pos(2));
	bool ins = dynamic_cast<t_gsetinp*>(_set)->input_size("imu") > 0 ? true : false;
	// write kml
	if (_kml && ins) {
		ostringstream out;
		out << fixed << setprecision(11) << " " << setw(0) << Geo_pos[1] << ',' << Geo_pos[0];
		string val = out.str();

		xml_node root = _doc;
		xml_node node = this->_default_node(root, _root.c_str());
		xml_node document = node.child("Document");
		xml_node last_child = document.last_child();
		xml_node placemark = document.insert_child_after("Placemark", last_child);
		this->_default_node(placemark, "styleUrl", "#P2");
		this->_default_node(placemark, "time", dbl2str(sins.t).c_str());
		xml_node point = this->_default_node(placemark, "Point");
		this->_default_node(point, "coordinates", val.c_str()); // for point

		xml_node Placemark = document.child("Placemark");
		xml_node LineString = Placemark.child("LineString");
		this->_default_node(LineString, "coordinates", val.c_str(), false);  // for line
	}

	return 1;
}

void gsins::t_gintegration::_prt_port()
{
	gtrace("t_gintegration::_prt_port");
	set<string> ambs = _param.amb_prns();
	int nsat = ambs.size();
	// get amb status
	string amb = "Float";
	if (_amb_state)amb = "Fixed";
	Eigen::Vector3d XYZ_INS = sins.pos_ecef + sins.Ceb*lever;
	ostringstream os;
	os << fixed << setprecision(4) << " "
		<< " " << _ins_crt.sow() + _ins_crt.dsec()
		<< fixed << setprecision(4)
		<< " " << setw(15) << XYZ_INS[0]          // [m]
		<< " " << setw(15) << XYZ_INS[1]          // [m]
		<< " " << setw(15) << XYZ_INS[2]          // [m]
		<< " " << setw(10) << sins.vn[0]          // [m]
		<< " " << setw(10) << sins.vn[1]          // [m]
		<< " " << setw(10) << sins.vn[2]          // [m]
		<< fixed << setprecision(4)
		<< " " << setw(10) << sins.att(0) / glv.deg
		<< " " << setw(10) << sins.att(1) / glv.deg
		<< " " << setw(10) << -sins.att(2) / glv.deg
		<< fixed << setprecision(4) 
		<< " " << setw(15) << sins.eb(0) / glv.dph
		<< " " << setw(15) << sins.eb(1) / glv.dph
		<< " " << setw(15) << sins.eb(2) / glv.dph
	    << fixed << setprecision(4)
		<< " " << setw(15) << sins.db(0) / glv.mg
		<< " " << setw(15) << sins.db(1) / glv.mg
		<< " " << setw(15) << sins.db(2) / glv.mg
		<< fixed << setprecision(0)
		<< " " << setw(5) << nsat            // nsat
		<< fixed << setprecision(1)
		<< fixed << setprecision(2)
		<< " " << setw(5) << _dop.pdop()            // pdop
		<< fixed << setprecision(2)
		<< " " << setw(8) << amb
		<< endl;
	if (_maptcp.find(FLT_OUT) != _maptcp.end())_maptcp[FLT_OUT]->run_send(os.str());
}

int gsins::t_gintegration::_init()
{
	gtrace("t_gintegration::_init");

	if (!_ins_init() || !_gnss_init())
	{
		if (_log) _log->comment(1, "gintegration", "initialize integration filter failed!!! ");
		return -1;
	}

	_ins_crt = _gnss_crt = t_gtime::current_time(t_gtime::GPS);

	//if (dynamic_cast<t_gsetinp*>(_set)->input_size("rinexo") == 0)
	//{
	//	if (_log) _log->comment(1, "gintegration", "GNSS observation is not existing!!! ");
	//	return -1;
	//}

	//_gnss_beg = dynamic_cast<t_gsetgen*>(_setkf)->beg();
	//_gnss_end = dynamic_cast<t_gsetgen*>(_setkf)->end();
	_gnss_beg = _beg_time;
	_gnss_end = _end_time;

	if (_imudata)
	{
		double start = _imudata->beg_obs(_beg_end);
		double end = _imudata->end_obs(_beg_end);
		double start_set = dynamic_cast<t_gsetins*>(_set)->start();
		double end_set = dynamic_cast<t_gsetins*>(_set)->end();
		start = (start < start_set) ? start_set : start;
		end = (end < end_set) ? end : end_set;

		if (_beg_end)
		{
			_gnss_crt = _gobs->beg_obs(_site);
			_ins_beg = t_gtime(_gnss_beg.gwk(), start);
			_ins_end = t_gtime(_gnss_end.gwk(), end);
		}
		else {
			_gnss_crt = _gobs->end_obs(_site);
			_ins_beg = t_gtime(_gnss_beg.gwk(), end);
			_ins_end = t_gtime(_gnss_end.gwk(), start);
		}
		_ins_crt = _ins_beg;
	}
	else
	{
		if (_log) _log->comment(1, "gintegration", "IMU observation is not existing!!! ");
		return -1;
	}
	if (dynamic_cast<t_gsetproc*>(_set)->realtime())_ins_end = LAST_TIME;

	return 1;
}

bool gsins::t_gintegration::_ins_init()
{
	gtrace("t_gintegration::_ins_init");

	kftk = sins.t; MeasVel = MeasPos = Eigen::Vector3d::Zero(), Flag = NO_MEAS, fltflag = 1;
	this->lever = dynamic_cast<t_gsetign*>(_setkf)->lever();
	Eigen::VectorXd v(nq), vmax(nq), vmin(nq);

	Eigen::Vector3d angle_tmp, v_tmp, pos_tmp, eb_tmp, db_tmp;
	angle_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_misalignment_std();
	v_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_vel_std();
	pos_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_pos_std();
	eb_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_gyro_std();
	db_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_acce_std();
	v << angle_tmp*glv.deg, v_tmp, pos_tmp, eb_tmp*glv.dph, db_tmp*glv.mg;
	v = v.array().abs2();
	Pk = v.array().matrix().asDiagonal();

	angle_tmp = dynamic_cast<t_gsetign*>(_setkf)->min_misalignment_std();
	v_tmp = dynamic_cast<t_gsetign*>(_setkf)->min_vel_std();
	pos_tmp = dynamic_cast<t_gsetign*>(_setkf)->min_pos_std();
	eb_tmp = dynamic_cast<t_gsetign*>(_setkf)->min_gyro_std();
	db_tmp = dynamic_cast<t_gsetign*>(_setkf)->min_acce_std();
	vmin << angle_tmp*glv.min, v_tmp, pos_tmp, eb_tmp*glv.dph, db_tmp*glv.mg;
	Pmin = vmin.array().abs2();
	vmax = vmin * 1e6;
	Pmax = vmax.array().abs2();

	angle_tmp = dynamic_cast<t_gsetign*>(_setkf)->misalignment_psd();
	v_tmp = dynamic_cast<t_gsetign*>(_setkf)->vel_psd();
	pos_tmp = dynamic_cast<t_gsetign*>(_setkf)->pos_psd();
	eb_tmp = dynamic_cast<t_gsetign*>(_setkf)->gyro_psd();
	db_tmp = dynamic_cast<t_gsetign*>(_setkf)->acce_psd();
	Qt << angle_tmp*glv.dpsh, v_tmp*glv.mgpsHz, pos_tmp, eb_tmp*glv.dphpsh, db_tmp*glv.mgpsh;
	Qt = Qt.array().abs2();

	v_tmp = dynamic_cast<t_gsetign*>(_setkf)->meas_vel_noise();
	pos_tmp = dynamic_cast<t_gsetign*>(_setkf)->meas_pos_noise();
	Rt << v_tmp, pos_tmp;
	Rt = Rt.array().abs2();

	Rmax = Rt * 1000; Rmin = Rt*0.01;
	Rb << 0.9, 0.9, 0.9, 0.9, 0.9, 0.9;
	FBTau << 1.0, 1.0, 10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0;
	FBMax << glv.INF, glv.INF, glv.INF, glv.INF, glv.INF, glv.INF, glv.INF, glv.INF, glv.INF,
		3000.0*glv.dph, 3000.0*glv.dph, 3000.0*glv.dph, 50.0*glv.mg, 50.0*glv.mg, 50.0*glv.mg;

	return true;
}

bool gsins::t_gintegration::_gnss_init()
{
	if ((_ign_type == IGN_TYPE::TCI || _ign_type == IGN_TYPE::MULTIGN_TC || _ign_type==IGN_TYPE::MULTIGN_LIDAR_TC))
	{
		try
		{
			t_gallpar param_extended;
			SymmetricMatrix Qx_extended;

			// attitude
			for (int ipar = t_gpar::ATT_X; ipar <= t_gpar::ATT_Z; ipar++)
				param_extended.addParam(t_gpar(_site, t_gpar::t_type(ipar), ipar, ""));
			// velocity
			for (int ipar = t_gpar::VEL_X; ipar <= t_gpar::VEL_Z; ipar++)
				param_extended.addParam(t_gpar(_site, t_gpar::t_type(ipar), ipar, ""));
			// position 
			for (int ipar = t_gpar::CRD_X; ipar <= t_gpar::CRD_Z; ipar++)
				param_extended.addParam(t_gpar(_site, t_gpar::t_type(ipar), ipar, ""));
			// bias
			for (int ipar = t_gpar::eb_X; ipar <= t_gpar::db_Z; ipar++)
				param_extended.addParam(t_gpar(_site, t_gpar::t_type(ipar), ipar, ""));

			int icrdz_gnss = _param.getParam(_site, t_gpar::CRD_Z, "");

			for (int i = icrdz_gnss + 1; i < _param.parNumber(); i++)
			{
				param_extended.addParam(_param[i]);
			}

			int idbz = param_extended.getParam(_site, t_gpar::db_Z, "");

			Qx_extended.ReSize(param_extended.parNumber()); Qx_extended = 0.0;

			if (idbz + 2 < param_extended.parNumber())
			{
				Qx_extended.SymSubMatrix(idbz + 2, param_extended.parNumber()) =
					_Qx.SymSubMatrix(icrdz_gnss + 2, _param.parNumber());
			}

			for (int i = 0; i < idbz + 1; i++)
			{
				Qx_extended(i + 1, i + 1) = 1.0;
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
	return true;
}

bool gsins::t_gintegration::_valid_ins_constraint()
{
	int sig_init_crd = dynamic_cast<t_gsetproc*>(_set)->sig_init_crd();
	double ins_error = SQRT(Pk(6, 6) + Pk(7, 7) + Pk(8, 8));
	// the pure ins time threshold is maybe shorter than this for low-grade imu
	if (ins_error > sig_init_crd || sins.pure_ins_time > 30)return false;
	return true;
}

int gsins::t_gintegration::_prepareData()
{
	int irc = t_gpvtflt::_prepareData();
	return irc;
}

int gsins::t_gintegration::_prepareData(const string & ssite, vector<t_gsatdata>& sdata)
{
	int irc = t_gpvtflt::_prepareData(ssite, sdata);
	return irc;
}

int gsins::t_gintegration::_processEpoch(const t_gtime & runEpoch)
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
		if (_log) _log->comment(1, "gintegration", "No receiver settings available!!!");
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
			if (_log) _log->comment(2, "gintegration", "Not enough visible satellites!");
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

int gsins::t_gintegration::_merge(Matrix & A)
{
	gtrace("t_gintegration::_merge");

	try
	{
		// reset _param_ign first of all(because of outlier loop)
		int size = _param.parNumber();
		int nobs = A.Nrows();
		// A 
		int icrdx = _param[_param.getParam(_site, t_gpar::CRD_X, "")].index;
		A.SubMatrix(1, nobs, icrdx, icrdx + 2) = -A.SubMatrix(1, nobs, icrdx, icrdx + 2);
		int iattx = _param[_param.getParam(_site, t_gpar::ATT_X, "")].index;
		A.SubMatrix(1, nobs, iattx, iattx + 2) = A.SubMatrix(1, nobs, icrdx, icrdx + 2)*Eigen2newMat(t_gbase::askew(sins.Ceb*lever));
		tmeas = _gnss_crt.sow() + _gnss_crt.dsec();
		if (fabs(sins.t - tmeas) > 0.5) throw("merge failed");
		// Qx
		int ins_size = dynamic_cast<t_gsetign*>(_set)->nq(); Pk = Pk_Sav;  // because of outlier loop
		Eigen::MatrixXd GNSS_Variance = NewMat2Eigen(_Qx);
		GNSS_Variance.block(0, 0, ins_size, ins_size) = Pk; Pk = GNSS_Variance;
		Pk.block(ins_size, 0, size - ins_size, ins_size) = Pk.block(ins_size, 0, size - ins_size, ins_size)*Phik.transpose();
		Pk.block(0, ins_size, ins_size, size - ins_size) = Phik * Pk.block(0, ins_size, ins_size, size - ins_size);
		_Qx = Eigen2BaseMatrix(Pk);
		Pk.conservativeResize(ins_size, ins_size);

		// maybe some questions in the way
		//_Qx.SubMatrix(1, ins_size, 1, ins_size) = Eigen2BaseMatrix(Pk);
		//_Qx.SubMatrix(1, ins_size, ins_size + 1, size) = Eigen2newMat(Phik)*_Qx.SubMatrix(1, ins_size, ins_size + 1, size);
		//_Qx.SubMatrix(ins_size + 1, size, 1, ins_size) = _Qx.SubMatrix(1, ins_size, ins_size + 1, size).t();

	}
	catch (...)
	{
		if (_log) _log->comment(1, "gintegration", "GNSS and INS infomation fusion failed!");
		_reset();
		Flag = NO_MEAS;
		return -1;
	}
	return 1;
}

int gsins::t_gintegration::_merge_init()
{
	t_gtriple gnss_pos;
	_param.getCrdParam(_site, gnss_pos);
	sins.eth.Update(Cart2Geod(Eigen::Vector3d(gnss_pos[0], gnss_pos[1], gnss_pos[2]), false), Eigen::Vector3d::Zero());
	sins.Cnb = t_gbase::q2mat(sins.qnb);
	sins.pos_ecef = Eigen::Vector3d(gnss_pos[0], gnss_pos[1], gnss_pos[2]) - sins.eth.Cen*sins.Cnb * lever;
	sins.pos = Cart2Geod(sins.pos_ecef, false);
	int icrdx = _param[_param.getParam(_site, t_gpar::CRD_X, "")].index;
	Pk.block(icrdx - 1, icrdx - 1, 3, 3) = NewMat2Eigen(_Qx).block(icrdx - 1, icrdx - 1, 3, 3);
	for (int i = icrdx + 3; i <= _param.parNumber(); i++)
		for (int j = icrdx; j <= icrdx + 2; j++)
			_Qx(i, j) = -_Qx(i, j);
	_initial_merge = true;
	return 1;
}

int gsins::t_gintegration::_GNSS_Update()
{
	t_gposdata::data_pos posdata;

	if (_ign_type == IGN_TYPE::LCI)
	{
		Flag = _getPOS(posdata);
		if (Flag != NO_MEAS) {
			_meas_update();
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
	else if (_ign_type == IGN_TYPE::TCI)
	{
		_timeUpdate(_gnss_crt);
		Flag = MEAS_TYPE(this->_processEpoch(_gnss_crt));
		if (Flag == NO_MEAS) { _reset(); }
	}

	if (Flag == NO_MEAS)return -1;

	return 1;
}

int gsins::t_gintegration::_ZUPT_Update()
{
	set_meas_ZUPT(_ins_crt.sow() + _ins_crt.dsec());
	_meas_update();
	if (outlier_detect()) { _reset(); Flag = NO_MEAS; }

	if (Flag == NO_MEAS)
		return -1;
	return 1;
}

int gsins::t_gintegration::_ZUPT_POS_Update()
{
	set_meas_ZUPT_POS(_ins_crt.sow() + _ins_crt.dsec());
	_meas_update();
	if (outlier_detect()) { _reset(); Flag = NO_MEAS; }

	if (Flag == NO_MEAS)
		return -1;
	return 1;
}


int gsins::t_gintegration::_NHC_Update()
{
	set_meas_NHC(_ins_crt.sow() + _ins_crt.dsec());
	_meas_update();
	if (outlier_detect()) { _reset(); Flag = NO_MEAS; }

	if (Flag == NO_MEAS)return -1;
	return 1;
}

int gsins::t_gintegration::_ATT_Update()
{
	set_meas_ATT(_ins_crt.sow() + _ins_crt.dsec());
	_meas_update();
	if (outlier_detect()) { _reset(); Flag = NO_MEAS; }

	if (Flag == NO_MEAS)
		return -1;
	return 1;
}

MEAS_TYPE gsins::t_gintegration::_getPOS(t_gposdata::data_pos& pos)
{
	gtrace("t_gintegration::_getPOS");

	MEAS_TYPE res_type;
	double crt = _ins_crt.sow() + _ins_crt.dsec();
	t_gtime runEpoch = _gobs->load(_site, crt);
	int irc = ProcessOneEpoch(runEpoch);
	if (irc < 0) return MEAS_TYPE::NO_MEAS;
	_get_result(runEpoch, pos);

	MeasVel = pos.vn*fltflag; MeasPos = pos.pos; tmeas = pos.t;
	_Cov_MeasVn = pos.Rvn; _Cov_MeasPos = pos.Rpos;

	if (!_aligned)
	{
		sins.pos = Cart2Geod(MeasPos, false);
		sins.vn = t_gbase::Cen(sins.pos).transpose()*  MeasVel;
	}
	res_type = MEAS_TYPE::POS_MEAS;
	if (!double_eq(MeasVel.norm(), 0.0))
		res_type = POS_VEL_MEAS;

	if (pos.PDOP > _shm.max_pdop)res_type = NO_MEAS;
	if (pos.nSat < _shm.min_sat)res_type = NO_MEAS;
	//if (_isBase && !pos.amb_state)res_type = NO_MEAS;

	//if (_sampling > 1) _gnss_crt.add_secs(int(_sampling));  // =<1Hz data
	//else               _gnss_crt.add_dsec(_sampling);       //  >1Hz data

	return res_type;

}

MEAS_TYPE gsins::t_gintegration::_get_DD_Fixed_Obs()
{
	if (_valid_ins_constraint() && 1)
	{
		Eigen::Vector3d XYZ_INS = sins.pos_ecef + sins.Ceb*lever;
		t_gtriple XYZ(XYZ_INS(0), XYZ_INS(1), XYZ_INS(2));
		double ins_error = SQRT(Pk(6, 6) + Pk(7, 7) + Pk(8, 8));
		t_gtriple XYZ_RMS(0.2, 0.2, 0.2);
		_external_pos(XYZ, XYZ_RMS);
	}

	MEAS_TYPE res_type;
	double crt = _ins_crt.sow() + _ins_crt.dsec();
	t_gtime runEpoch = _gobs->load(_site, crt);
	int irc = ProcessOneEpoch(runEpoch);
	if (irc == 0) return MEAS_TYPE::NO_MEAS; // Float

	_gnss_design_mat = NewMat2Eigen(_post_A);
	_gnss_meas_vector = Columns2VectorXd(_post_l);
	_gnss_meas_cov = BaseMatrix2Eigen(_post_P).inverse();

	tmeas = runEpoch.sow() + runEpoch.dsec();

	return FIXED_DD_UC_OBS_MEAS;
}

int gsins::t_gintegration::_getMeas()
{
	gtrace("t_gintegration::_getMeas");
#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif 
	//	_gmutex.lock();

	_Meas_Type.clear(); Flag = NO_MEAS;
	t_gsinskf::_restore(Xk, Pk);

	if (_time_valid(_gnss_crt, _ins_crt))
		_Meas_Type.insert(MEAS_TYPE::GNSS_MEAS);

	// 2Hz
	if (fmod(sins.pure_ins_time, 0.5) < 1e-6 || fmod(sins.pure_ins_time, 0.5) > 0.5 - 1e-6)
	{
		MEAS_TYPE meas_type = meas_state();
		if (meas_type != NO_MEAS)_Meas_Type.insert(meas_type);
	}

	if (_mm.size() != 0)
	{
		//_Meas_Type.insert(MEAS_TYPE::ATT_MEAS);
	}

	//	_gmutex.unlock();
	return _Meas_Type.size();
}


int gsins::t_gintegration::_getMeas_test()
{
	gtrace("t_gintegration::_getMeas");
#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif 
	//	_gmutex.lock();

	_Meas_Type.clear(); Flag = NO_MEAS;
	t_gsinskf::_restore(Xk, Pk);

	//if (_time_valid(_gnss_crt, _ins_crt))
		//_Meas_Type.insert(MEAS_TYPE::GNSS_MEAS);

	// 2Hz
	double nt = 0.5;
	if (fmod(sins.pure_ins_time, nt) < 1e-6 || fmod(sins.pure_ins_time, nt) > nt - 1e-6)
	{
		MEAS_TYPE meas_type = meas_state();
		/*&& meas_type!= NHC_MEAS*/
		if (meas_type != NO_MEAS)_Meas_Type.insert(meas_type);
	}

	if (_mm.size() != 0)
	{
		//_Meas_Type.insert(MEAS_TYPE::ATT_MEAS);
	}

	//	_gmutex.unlock();
	return _Meas_Type.size();
}

bool gsins::t_gintegration::_time_valid(t_gtime gt, t_gtime inst)
{
	gtrace("t_gintegration::_time_valid");

	bool res_valid = false;
	double crt = inst.sow() + inst.dsec();
	_data.erase(_data.begin(), _data.end());
	t_gtime obsEpo = _gobs->load(_site, crt);

	if ((abs(inst.diff(obsEpo)) < _shm.delay && inst >= obsEpo)
		|| abs(inst.diff(obsEpo)) < 1e-6)
	{
		if (_isClient) {//lvhb added for RTK/INS
			if (!_gInterpol->interpolAug(obsEpo,_data,_data_base)){
				return false;
			}
		}
		//detect cycle slip by lvhb
		if (dynamic_cast<t_gsetproc*>(_set)->slip_model() == SLIPMODEL::TURBO_EDIT) {
			if (_isClient) _gturboedit->setSingleEpoData(&_data);
			_gturboedit->ProcessBatch(_site, obsEpo, obsEpo, _sampling);
			if (_isBase && !_isClient) {
				_gturboedit->ProcessBatch(_site_base, obsEpo, obsEpo, _sampling);
			}
		}
		else {
			if (_isClient) _gpre->setSingleEpoData(&_data);
			_gpre->preprocess(_site, obsEpo, obsEpo, _sampling, false);
			if (_isBase && !_isClient) {
				_gpre->preprocess(_site_base, obsEpo, obsEpo, _sampling, false);
			}
		}
		if(!_isClient) _data = _gobs->obs(_site, obsEpo);
		if (_data.size() > 0) {
			res_valid = true;
			// apply dcb
			if (_gallbias) {
				for (auto& itdata : _data) {
					itdata.apply_bias(*_gallbias);
				}
			}
		}
		else {
			if (_log && _log->verb() >= 1) {
				_log->comment(1, "gintegration", _site + gt.str_ymdhms(" no observation found at epoch: "));
			}
			res_valid = false;
		}
		if (_isBase)
		{
			if (!_isClient) {
				_data_base.erase(_data_base.begin(), _data_base.end());
				_data_base = _gobs->obs(_site_base, obsEpo);
			}
			if (_data_base.size() > 0) {
				res_valid = true;
				// apply dcb
				if (_gallbias) {
					for (auto& itdata : _data_base) {
						itdata.apply_bias(*_gallbias);
					}
				}
			}
			else {
				if (_log && _log->verb() >= 1) {
					_log->comment(1, "gintegration", _site_base + gt.str_ymdhms(" no base observation found at epoch: "));
				}
				res_valid = false;
			}
		}
	}

	// judge ultimately(for rtk)
	if (_data.size() == 0)res_valid = false;
	if (res_valid)_gnss_crt = obsEpo;
	return res_valid;
}

bool gsins::t_gintegration::_avaliable(t_gtime crt)
{
	gtrace("t_gintegration::_avaliable");

	bool obs, corr, imu;
	corr = false;

	_gobs->erase(_site, _ins_crt);
	if (_isBase) _gobs->erase(_site_base, _ins_crt);
	_gobs->setepoches(_site);

	do {

		obs = (crt + 1.0 < _gobs->end_obs(_site));
		imu = _imudata->available(crt,_beg_end);
		if (_isBase)corr = true;
		else corr = dynamic_cast<t_gallprec*>(_gnav)->corr_avali(crt + 1.0);

		//obs = ((t_bncobs*)_gcoder[0])->available(crt + 1.0);
		//imu = ((t_imufile*)_gcoder[3])->available(crt);
		//if (_isBase)corr = true;
		//else 	corr = ((t_bnccorr*)_gcoder[1])->available(crt + 1.0);

		//t_gtime::gmsleep(_shm.ts * 100);

	} while (!obs || !imu || !corr);

	return true;
}

void gsins::t_gintegration::time_update(double kfts, int fback, double inflation)
{
	gtrace("t_gintegration::time_update");

	set_Ft1();

	Eigen::MatrixXd Fk = Eigen::MatrixXd::Identity(nq, nq) + (Ft*kfts);

	Phik = Fk * Phik;
	Xk = Fk * Xk;
	Eigen::MatrixXd Qk = (Qt*kfts*inflation).array().matrix().asDiagonal();
	Pk = Fk * Pk*(Fk.transpose());  Pk += Qk;
	Pk_Sav = Pk;
	if (fback)  feedback(kfts);

}

int gsins::t_gintegration::_meas_update(double fading)
{
	gtrace("t_gintegration::_meas_update");

	Eigen::VectorXd Pxz = Eigen::VectorXd::Zero(nq),
		Kk = Eigen::VectorXd::Zero(nq), Hi = Eigen::VectorXd::Zero(nq);


	if (!_set_meas()) {
		if (_log) _log->comment(1, "gintegration", _ins_crt.str_ymdhms(" _set_meas failed: "));
		_reset();
		Flag = NO_MEAS;
		return -1;
	}

#ifdef DEBUG
	t_out("Zk:", Zk.transpose());
	t_out("Xk:", Xk.transpose());
	t_out("Hk:", Hk);
	t_out("Rk:", Rk);
#endif
	// before estimating 
	//if (outlier_detect()) { _reset(); Flag = NO_MEAS; return -1; }

	for (int i = 0; i < nr; i++)
	{
		Hi = Hk.block(i, 0, 1, nq).transpose();
		Pxz = Pk*Hi;
		double Pz0 = (Hi.transpose()*Pxz), r = Zk(i) - (Hi.transpose()*Xk);
		double Pzz = Pz0 + Rk(i, i);
		Kk = Pxz*(1.0 / Pzz);
		Xk += Kk*r;
		Pk = Pk - Kk*Pxz.transpose();
	}
	if (fading > 1.0) Pk = Pk* fading;
	t_gbase::symmetry(Pk);
	Phik = Eigen::MatrixXd::Identity(15, 15);
	set_meas_flag(0);

#ifdef DEBUG
	t_out("Xk:", Xk.transpose());
#endif

	// after estimating
	if (outlier_detect()) { _reset(); Flag = NO_MEAS; return -1; }

	return 1;
}

bool gsins::t_gintegration::_set_meas()
{
	double dt = sins.t - tmeas;
	if (fabs(dt) > _shm.delay) return false;
	_resize();
	Eigen::Vector3d vartheta = sins.Ceb*t_gbase::askew(lever)*sins.wib + t_gbase::askew(sins.eth.weie)*sins.Ceb*lever;
	switch (Flag)
	{
	case POS_VEL_MEAS:
		Zk.block(0, 0, 3, 1) = sins.ve - MeasVel - sins.ae*dt - vartheta;
		Rk.block(0, 0, 3, 3) = _Cov_MeasVn.array().matrix().asDiagonal();
		set_meas_flag(007);
		Zk.block(3, 0, 3, 1) = sins.pos_ecef - MeasPos - sins.ve*dt + sins.Ceb*lever;
		Rk.block(3, 3, 3, 3) = _Cov_MeasPos.array().matrix().asDiagonal();
		set_meas_flag(070);
		MeasVel = MeasPos = Eigen::Vector3d::Zero();
		break;
	case POS_MEAS:
		Zk.block(0, 0, 3, 1) = sins.pos_ecef - MeasPos - sins.ve*dt + sins.Ceb*lever;
		Rk.block(0, 0, 3, 3) = _Cov_MeasPos.array().matrix().asDiagonal();
		set_meas_flag(070);
		MeasVel = MeasPos = Eigen::Vector3d::Zero();
		break;
	case VEL_MEAS:
		Zk.block(0, 0, 3, 1) = sins.ve - MeasVel  - sins.ae*dt + sins.Ceb*t_gbase::askew(sins.web)*lever;
		Rk.block(0, 0, 3, 3) = _Cov_MeasVn.array().matrix().asDiagonal();
		set_meas_flag(007);
		MeasVel = MeasPos = Eigen::Vector3d::Zero();
		break;
	case FIXED_DD_UC_OBS_MEAS:
		Zk = _gnss_meas_vector;
		Rk = _gnss_meas_cov;
		set_meas_flag(001);
		break;
	case ZUPT_MEAS:  // MeasVel=Ve
		Zk.block(0, 0, 3, 1) = MeasVel ;
		cout << Zk << endl;
		Rk.block(0, 0, 3, 3) = _Cov_MeasVn.array().matrix().asDiagonal();
		set_meas_flag(007);
		MeasVel = Eigen::Vector3d::Zero();
		break;
	case NHC_MEAS:  // MeasVel(0)=vb(0) and  MeasVel(2)=vb(2) 
		Zk(0) = MeasVel(0); Zk(1) = MeasVel(2);
		Rk(0, 0) = _Cov_MeasVn(0); Rk(1, 1) = _Cov_MeasVn(2);
		MeasVel = Eigen::Vector3d::Zero();
		set_meas_flag(005);
		break;
	case ODO_MEAS: // MeasVf=odo_vel
		Zk(0) = MeasVf;
		Rk(0, 0) = SQR(0.1);
		MeasVf = 0;
		set_meas_flag(001);
	case Hgt_MEAS:
		Zk(0) = MeasHgt;
		Rk(0, 0) = SQR(0.5);
		MeasHgt = 0;
		set_meas_flag(001);
	case YAW_MEAS:
		Zk(0) = MeasYaw;
		Rk(0, 0) = SQR(0.1*glv.deg);
		MeasYaw = 0;
		set_meas_flag(001);
		break;
	case ATT_MEAS:
		Zk.block(0, 0, 3, 1) = MeasAtt - sins.att;   ///? attµÄ×ø±êÏµ
		Rk.block(0, 0, 3, 3) = _Cov_MeasAtt.array().matrix().asDiagonal();
		set_meas_flag(007);
		MeasAtt = Eigen::Vector3d::Zero();
		_mm.clear();
		break;
	case ZUPT_POS_MEAS:
		Zk.block(0, 0, 3, 1) = sins.pos_ecef - MeasPos;
		cout << Zk << endl;
		Rk.block(0, 0, 3, 3) = _Cov_MeasPos.array().matrix().asDiagonal();
		set_meas_flag(070);
		MeasPos = Eigen::Vector3d::Zero();
		break;
	case NO_MEAS:
		set_meas_flag(0);
		break;
	default:
		break;
	}
	if (this->measflag)set_Hk();
	return true;
}

bool gsins::t_gintegration::_resize()
{
	switch (Flag)
	{
	case POS_VEL_MEAS:
		nr = 6; break;
	case POS_MEAS:
		nr = 3; break;
	case VEL_MEAS:
		nr = 3; break;
	case ZUPT_MEAS:
		nr = 3; break;
	case NHC_MEAS:
		nr = 2; break;
	case ODO_MEAS:
		nr = 1; break;
	case Hgt_MEAS:
		nr = 1; break;
	case YAW_MEAS:
		nr = 1; break;
	case FIXED_DD_UC_OBS_MEAS:
		nr = _gnss_meas_vector.rows(); break;
	case ZUPT_POS_MEAS:
		nr = 3; break;
	case NO_MEAS:
		break;
	default:
		break;
	}
	Hk = Eigen::MatrixXd::Zero(nr, nq); Zk = Eigen::VectorXd::Zero(nr);
	Rk = Eigen::MatrixXd::Zero(nr, nr);
	return true;
}

int gsins::t_gintegration::outlier_detect()
{
	gtrace("t_gintegration::outlier_detect");

	//return 0;
	double max_set = dynamic_cast<t_gsetproc*>(_set)->max_res_norm();
	// compute normalized measument res
	_v_post = Zk - Hk*Xk;
#ifdef DEBUG
	t_out("v_post:", _v_post.transpose());
#endif
	Eigen::MatrixXd Qz = Hk*Pk*Hk.transpose() + Rk;
	Eigen::VectorXd v_norm = Eigen::VectorXd::Zero(nr);
#ifdef DEBUG
	t_out("Qz:", Qz.diagonal().transpose());
#endif
	for (int i = 0; i < nr; i++)
	{
		v_norm(i) = sqrt(1.0 / Qz(i, i))*abs(_v_post(i));
	}
#ifdef DEBUG
	t_out("v_norm:", v_norm.transpose());
#endif
	// detect the outlier
	double MaxNormRes = 0.0;
	int MaxIndex = 1000;
	for (int i = 0; i < nr; i++)
	{
		if (v_norm[i] > max_set)
		{
			Flag = NO_MEAS;
			Pk = Pk_Sav;
			if (_log)
				_log->comment(1, "gintegration", "obs(" + to_string(i) + ")" + "GNSS and INS loosely coupled outlier!");
			return 1;
		}
	}
	return 0;

}

/// <summary>
/// ZUPT_Check for foot mount imu 
/// </summary>
/// <returns></returns>
bool gsins::t_gintegration::_ZUPT_Check()
{
	int slide_window = 20;
	auto all_tz = cal_Tz();

	for (double tz : all_tz)
	{
		if (_Tz.size() == 0)
		{
			_Tz.push_back(tz);
		}
		else if (_Tz.size() < slide_window)
		{
			double tmp = (_Tz[_Tz.size() - 1] * _Tz.size() + tz) / (_Tz.size() + 1);
			_Tz.push_back(tmp);
		}
		else
		{
			double tmp = _Tz[_Tz.size() - 1] + (tz - _Tz[_Tz.size() - slide_window]) / slide_window;
			_Tz.push_back(tmp);
			double min_tz = *min_element((_Tz.end() - slide_window), _Tz.end());
			double max_tz = *max_element((_Tz.end() - slide_window), _Tz.end());
			_thresh = min_tz + 0.25 * (max_tz - min_tz);

			if (_thresh > tz)
			{
				_Meas_Type.insert(ZUPT_MEAS);
			}
		}
	}
	return 1;

}

vector<double> gsins::t_gintegration::cal_Tz()
{
	vector<double> res;
	int slide_window = 20;
	double a_std = dynamic_cast<t_gsetign*>(_setkf)->initial_acce_std()[0];
	double g_std = dynamic_cast<t_gsetign*>(_setkf)->initial_gyro_std()[0];
	for (int i = 0; i < _wm.size(); i++)
	{
		Eigen::Vector3d f = _vm[i] / sins.nts;
		Eigen::Vector3d w = _wm[i] / sins.nts;
		if (_f.size() == 0)
		{
			_mf = f;  ///< initialize
		}
		else if (_f.size() < slide_window)
		{
			_mf = _mf * _f.size() + f;
			_mf = _mf / (_f.size() + 1);
		}
		else
		{
			_mf = _mf + (f - _f[_f.size() - slide_window] / slide_window);
		}
		_f.push_back(f);
		f = f - (_mf / _mf.norm() * sins.eth.gn.norm());
		res.push_back(f.norm() / (a_std * a_std) - w.norm() / (g_std * g_std));
	}
	return res;

}

//pdr test
int gsins::t_gintegration::processBatch_test(const t_gtime& beg, const t_gtime& end)
{
	t_gposdata::data_pos posdata;
	bool time_loop = true;
	_ins_crt = _ins_beg;
	_ins_init();
	// ins time loop
	while (time_loop)
	{
		if (_ins_crt > _ins_end)
		{
			cerr << " Processing Finished!" << endl;
			if (_log) _log->comment(0, "gintegration", "Processing Finished!");
			time_loop = false; break;
		}
		//if (!_imudata->load(_wm, _vm,_mm, _shm.t, _shm.ts, _shm.nSamples, true)) { cerr << "\n Processing Finished!" << endl; break; }
		if (!_imudata->load(_wm, _vm, _shm.t, _shm.ts, _shm.nSamples)) { cerr << "\n Processing Finished!" << endl; break; }
		_ins_crt = t_gtime(_gnss_beg.gwk(), _shm.t);
		_aligned = 1;
		if (0/*!_aligned*/)
		{
			//if (_time_valid(_gnss_crt, _ins_crt)) Flag = _getPOS(posdata);
			//_aligned = cascaded_align(posdata.pos, posdata.vn);
		}
		else
		{
			sins.Update(_wm, _vm, _shm);
			if (_ign_type == PURE_INS) {
				t_gsinskf::write();

				//wh test
				//continue;
			}

			time_update(sins.nts, 0.0);
			int imeas = _getMeas_test();
			int irc = -1;
			//_ZUPT_Check();
			set<MEAS_TYPE>::const_iterator it = _Meas_Type.begin();
			while (it != _Meas_Type.end())
			{
				switch (*it)
				{
				case GNSS_MEAS:
				{
					if (irc = _GNSS_Update() > 0)t_gsinskf::_restore(Xk, Pk);
					else t_gsinskf::_recover(Xk, Pk);
				}break;
				case ZUPT_MEAS:
				{
					if (irc = _ZUPT_Update() > 0) t_gsinskf::_restore(Xk, Pk);
					else t_gsinskf::_recover(Xk, Pk);
				}break;
				case NHC_MEAS:
				{
					if (irc = _NHC_Update() > 0) t_gsinskf::_restore(Xk, Pk);
					else t_gsinskf::_recover(Xk, Pk);
				}break;

				case ATT_MEAS:
				{
					/*
					if (irc = _ATT_Update() > 0) t_gsinskf::_restore(Xk, Pk);
					else t_gsinskf::_recover(Xk, Pk);
					*/
				}break;

				case ZUPT_POS_MEAS:
				{
					if (irc = _ZUPT_POS_Update() > 0) t_gsinskf::_restore(Xk, Pk);
					else t_gsinskf::_recover(Xk, Pk);
				}break;

				default:break;
				}it++;
			}

			if (Flag != NO_MEAS) {
				feedback(glv.INF);
				sins.pure_ins_time = 0.0;
				cerr << "\r" << _ins_crt.str_ymdhms();
				if (_log) _log->comment(0, "gintegration",
					_ins_crt.str_ymdhms(" epoch ") + "Processing GNSS/INS Integration Success!");
			}
		}
	}

	return 0;
}

void gsins::t_gintegration::Set_pro(t_gtime *fwd_crt, t_gtime *bwd_crt, double* fwd_per, double *bwd_per)
{
	_fwd_crt = fwd_crt;
	_bwd_crt = bwd_crt;
	_fwd_per = fwd_per;
	_bwd_per = bwd_per;
}
