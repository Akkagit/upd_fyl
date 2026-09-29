#include "gins.h"
#include <iostream>
#include <iomanip>
#include "gset/gsetout.h"
#include "gset/gsetins.h"
#include <Eigen/src/Core/util/DisableStupidWarnings.h>
using namespace gsins;
using namespace Eigen;
using namespace std;
using namespace gnut;



/******************************************* t_gsins ***********************************************/
t_gsins::t_gsins(const t_gquat& qnb0, const Vector3d& vn0, const Vector3d& pos0, double t0)
{
	qnb = qnb0; vn = vn0; pos = pos0; t = t0;
	eth.Update(pos0, vn0);
	Cnb = t_gbase::q2mat(qnb); att = t_gbase::q2att(qnb); Cbn = Cnb.transpose(); vb = Cbn*vn;
	Kg = Ka = Matrix3d::Identity();
	eb = db = _tauG = _tauA = Vector3d::Zero();
	wib = web = wnb = fb = fn = an = Vector3d::Zero();
	err_mat();

}

gsins::t_gsins::t_gsins(gnut::t_gsetbase * set)
{
	Vector3d Cart= dynamic_cast<t_gsetins*>(set)->pos();
	pos = Cart2Geod(Cart, false);
	vn = t_gbase::Cen(pos).transpose()*(dynamic_cast<t_gsetins*>(set)->vel());
	att = dynamic_cast<t_gsetins*>(set)->att();
	eth.Update(pos, vn);
	qnb = t_gbase::a2qua(att);
	Cnb = t_gbase::q2mat(qnb); Cbn = Cnb.transpose(); vb = Cbn*vn;
	Kg = Ka = Matrix3d::Identity();
	eb = db = _tauG = _tauA = Vector3d::Zero();
	wib = web = wnb = fb = fn = an = Vector3d::Zero();
	err_mat();
	_order = dynamic_cast<t_gsetins*>(set)->out_order();
	output_freq = dynamic_cast<t_gsetins*>(set)->out_freq();
	output_intsec = dynamic_cast<t_gsetins*>(set)->out_intsec();
	pure_ins_time = 1.0; 
}

void gsins::t_gsins::set_posvel(const Eigen::Vector3d & pos0, const Eigen::Vector3d & vn0)
{
	pos = pos0; vn = vn0;
}

void t_gsins::set_tau(const Vector3d& tauG, const Vector3d& tauA)
{
	_tauG(0) = tauG(0) > glv.INF / 2 ? 0.0 : 1.0 / tauG(0);
	_tauG(1) = tauG(1) > glv.INF / 2 ? 0.0 : 1.0 / tauG(1);
	_tauG(2) = tauG(2) > glv.INF / 2 ? 0.0 : 1.0 / tauG(2);
	_tauA(0) = tauA(0) > glv.INF / 2 ? 0.0 : 1.0 / tauA(0);
	_tauA(1) = tauA(1) > glv.INF / 2 ? 0.0 : 1.0 / tauA(1);
	_tauA(2) = tauA(2) > glv.INF / 2 ? 0.0 : 1.0 / tauA(2);
}
/************************************************************
** Function: SINS Initial Coarse Alignment on static base
** Input: Angular velocity and Velocity increment from m-1 to m and initial latitude
** See also P151 - P153
***********************************************************/
Vector3d t_gsins::align_coarse(const Vector3d& wmm,const Vector3d& vmm)
{
	double latitude = pos(0);
	double T11, T12, T13, T21, T22, T23, T31, T32, T33;
	double cb = cos(latitude), tb = tan(latitude), nn;
	Vector3d wbib = wmm / wmm.norm(), fb = vmm / vmm.norm();
	T31 = fb(0), T32 = fb(1), T33 = fb(2);
	T21 = wbib(0) / cb - T31*tb, T22 = wbib(1) / cb - T32*tb, T23 = wbib(2) / cb - T33*tb;       
	nn = sqrt(T21*T21 + T22*T22 + T23*T23);  T21 /= nn, T22 /= nn, T23 /= nn;
	T11 = T22*T33 - T23*T32, T12 = T23*T31 - T21*T33, T13 = T21*T32 - T22*T31;
	nn = sqrt(T11*T11 + T12*T12 + T13*T13);  T11 /= nn, T12 /= nn, T13 /= nn;
	Matrix3d Cnb;
	Cnb << T11, T12, T13, T21, T22, T23, T31, T32, T33;
	return t_gbase::m2att(Cnb);
}

Eigen::Vector3d gsins::t_gsins::acc_align(const Eigen::Vector3d & wmm)
{
	return Eigen::Vector3d();
}

/************************************************************
** Function: SINS Update ( Mechanical arrangement )
** Input: Angular velocity and Velocity increment from m-1 to m and Processing t_scheme
** See also P69 - P76
***********************************************************/
void t_gsins::Update(const vector<Vector3d>& wm, const vector<Vector3d>& vm, const t_scheme& scm)
{
	int n = wm.size();
	nts = scm.ts;
	t = scm.t;
	nts = abs(nts);
	double nts_2 = nts / 2.0;
	imu.Update(wm, vm, scm);
	imu.phim = Kg*imu.phim - eb*nts; imu.dvbm = Ka*imu.dvbm - db*nts;
	Vector3d vn1_2 = vn + an*nts_2, pos1_2 = pos + eth.v2dp(vn1_2, nts_2);
	eth.Update(pos1_2, vn1_2);
	wib = imu.phim / nts; fb = imu.dvbm / nts;
	web = wib - Cbn*eth.wnie;
	wnb = wib - (qnb*t_gbase::rv2q(imu.phim / 2))*eth.wnin;
	fn = qnb*fb;
	an = t_gbase::rv2q(-eth.wnin*nts_2)*fn + eth.gcc;
	Vector3d vn1 = vn + an*nts;
	pos = pos + eth.v2dp(vn + vn1, nts_2);	vn = vn1;
	qnb = t_gbase::rv2q(-eth.wnin*nts)*qnb*t_gbase::rv2q(imu.phim);
	Cnb = t_gbase::q2mat(qnb); att = t_gbase::m2att(Cnb); Cbn = Cnb.transpose(); vb = Cbn*vn;
	eth.Update(pos, vn); pos_ecef = Geod2Cart(pos, false);
	Ceb = eth.Cen*Cnb; Cbe = Ceb.transpose(); qeb = t_gbase::m2qua(Ceb); ve = eth.Cen*vn; ae = eth.Cen*an;
	pure_ins_time += nts;
}

/************************************************************
** Function: SINS Error Analysis
** Attitude, Velocity and Postion Error
** See also P78 - P82
***********************************************************/
void t_gsins::err_mat()
{
	double tb = eth.tb, secb = 1.0 / eth.cb, secb2 = secb*secb, scb = eth.sb*eth.cb,
		wN = eth.wnie(1), wU = eth.wnie(2), vE = vn(0), vN = vn(1);
	double f_RMh = eth.f_RMh, f_RNh = eth.f_RNh, f_cbRNh = eth.f_cbRNh,
		f_RMh2 = f_RMh*f_RMh, f_RNh2 = f_RNh*f_RNh;
	Matrix3d F1, F2, F3;
	F1 << 0, 0, 0, -wU, 0, 0, wN, 0, 0;
	F2 << 0, 0, vN*f_RMh2, 0, 0, -vE*f_RNh2, vE*secb2*f_RNh, 0, -vE*tb*f_RNh2;
	F3 << 0, 0, 0, 0, 0, 0, -glv.g0*(5.27094e-3 * 2 * scb + 2.32718e-5 * 4 * eth.sb2*scb), 0, 3.086e-6;
	Faa = t_gbase::askew(-eth.wnin);
	Fav << 0, -f_RMh, 0, f_RNh, 0, 0, tb*f_RNh, 0, 0;
	Fap = F1 + F2;
	Fva = t_gbase::askew(fn);
	Fvv = t_gbase::askew(vn)*Fav - t_gbase::askew(eth.wnie + eth.wnin);
	Fvp = t_gbase::askew(vn)*(2 * F1 + F2) + F3;
	Fpv << 0, f_RMh, 0, f_cbRNh, 0, 0, 0, 0, 1;
	Fpp << 0, 0, -vN*f_RMh2, vE*tb*f_cbRNh, 0, -vE*secb*f_RNh2, 0, 0, 0;
	//Fng << 0, eth.RMh, 0, eth.cbRNh, 0, 0, 0, 0, 1;
}

void gsins::t_gsins::prt_header(ostringstream & os)
{
	// the first line
	os << "% ";
	os << setw(15) << "Seconds of Week";
	switch (_order)
	{
		/// the att_out is positive if defination is from north to west[-Pi,Pi]
	case XYZ_XYZ_PRY_NW:
		os << setw(18) << "X-ECEF" <<
			  setw(18) << "Y-ECEF" <<
			  setw(18) << "Z-ECEF";
		os << setw(10) << "VX" <<
			  setw(10) << "VY" <<
			  setw(10) << "VZ";
		os << setw(10) << "Pitch" <<
			  setw(10) << "Roll" <<
			  setw(10) << "Yaw";
		break;

		/// the att_out is positive if defination is from north to west[-Pi,Pi]
	case BLH_ENU_PRY_NW:
		os << setw(18) << "Latitude" <<
			  setw(18) << "Longitude" <<
			  setw(18) << "Altitude";
		os << setw(10) << "Ve" <<
			  setw(10) << "Vn" <<
			  setw(10) << "Vu";
		os << setw(10) << "Pitch" <<
			  setw(10) << "Roll" <<
			  setw(10) << "Yaw";
		break;

		/// the att_out is positive if defination is from north to west[-Pi,Pi]
		/// HOLO format is Longitude,latitude,altitude,Roll,Pitch,Yaw,Ve,Vn,Vu
	case HOLO_ODO:
		os << setw(18) << "Longitude" <<
			  setw(18) << "Latitude" <<
			  setw(18) << "Altitude";
		os << setw(10) << "Pitch" <<
			  setw(10) << "Roll" <<
		 	  setw(10) << "Yaw";
		os << setw(10) << "Ve" <<
			  setw(10) << "Vn" <<
			  setw(10) << "Vu";
		break;

	default:
		break;
	}

	os << setw(15) << "GyroBiasX" <<
		setw(15) << "GyroBiasY" <<
		setw(15) << "GyroBiasZ" <<
		setw(15) << "AcceBiasX" <<
		setw(15) << "AcceBiasY" <<
		setw(15) << "AcceBiasZ";
	os << "         " << endl;


	// the second line
	os << "% ";
	os << setw(15) << "(s)";
	switch (_order)
	{
	case XYZ_XYZ_PRY_NW:
		os << setw(18) << "(m)" <<
			  setw(18) << "(m)" <<
			  setw(18) << "(m)";
		os << setw(10) << "(m/s)" <<
			  setw(10) << "(m/s)" <<
			  setw(10) << "(m/s)";
		os << setw(10) << "(deg)" <<
			  setw(10) << "(deg)" <<
			  setw(10) << "(deg)";
		break;

	case BLH_ENU_PRY_NW:
		os << setw(18) << "(deg)" <<
			  setw(18) << "(deg)" <<
			  setw(18) << "(m)";
		os << setw(10) << "(m/s)" <<
			  setw(10) << "(m/s)" <<
			  setw(10) << "(m/s)";
		os << setw(10) << "(deg)" <<
			  setw(10) << "(deg)" <<
			  setw(10) << "(deg)";
		break;

	case HOLO_ODO:
		os << setw(18) << "(deg)" <<
			  setw(18) << "(deg)" <<
			  setw(18) << "(m)";
		os << setw(10) << "(deg)" <<
			  setw(10) << "(deg)" <<
			  setw(10) << "(deg)";
		os << setw(10) << "(m/s)" <<
			  setw(10) << "(m/s)" <<
			  setw(10) << "(m/s)";
		break;
	default:
		break;
	}
	os << setw(15) << "(deg/h)";
	os << setw(15) << "(deg/h)";
	os << setw(15) << "(deg/h)";
	os << setw(15) << "(mg)";
	os << setw(15) << "(mg)";
	os << setw(15) << "(mg)";
	os << endl;

}

void gsins::t_gsins::prt_sins(ostringstream & os)
{
	Vector3d att1 = att / glv.deg;
	Vector3d eb1 = eb / glv.deg*glv.hur;
	Vector3d db1 = db / glv.mg;
	Vector3d Car_pos = Geod2Cart(pos, false);
	Vector3d Geo_pos = Vector3d(pos(0) / glv.deg, pos(1) / glv.deg, pos(2));
	Vector3d Ve = t_gbase::Cen(pos)*(vn);

	Vector3d pos_out, vel_out, att_out;
	switch (_order)
	{
	/// the att_out is positive if defination is from north to west[-Pi,Pi]
	case XYZ_XYZ_PRY_NW:
		pos_out = Car_pos;
		//pos_out = Car_pos + Ceb * Vector3d(-0.05, -0.325, 0.165);
		vel_out = ve;
		att_out = att1;
		att_out(2) = -att_out(2);
		break;

	/// the att_out is positive if defination is from north to west[-Pi,Pi]
	case BLH_ENU_PRY_NW:
		pos_out = Geo_pos;
		vel_out = vn;
		att_out = att1;
		att_out(2) = -att_out(2);
		break;

	/// the att_out is positive if defination is from north to west[-Pi,Pi]
	/// HOLO format is latitude,longitude,altitude,Roll,Pitch,Yaw,Ve,Vn,Vu
	case HOLO_ODO:
		pos_out = Vector3d(Geo_pos(0), Geo_pos(1), Geo_pos(2));
		// exchange vel_out and att_out for HOLO format
		vel_out = att1;
		att_out = vn;
		vel_out(2) = -vel_out(2);
		break;

	default:
		break;
	}

	os << fixed << setprecision(6) << setw(18)<< t;
	os << fixed << setprecision(3) <<
		setw(18) << pos_out(0) <<
		setw(18) << pos_out(1) <<
		setw(18) << pos_out(2);
	os << fixed << setprecision(3) <<
		setw(10) << vel_out(0) <<
		setw(10) << vel_out(1) <<
		setw(10) << vel_out(2);
	os << fixed << setprecision(4) <<
		setw(10) << att_out(0) <<
		setw(10) << att_out(1) <<
		setw(10) << att_out(2);
	os << fixed << setprecision(4) <<
		setw(15) << eb1(0) <<
		setw(15) << eb1(1) <<
		setw(15) << eb1(2);
	os << fixed << setprecision(4) <<
		setw(15) << db1(0) <<
		setw(15) << db1(1) <<
		setw(15) << db1(2);
	os << endl;
}

void gsins::t_gsins::int_sec(bool _beg_end)
{
	double out_nts = 1.0 / output_freq;
	double dt = fmod(t, out_nts);
	auto p = vec_qnb.begin();
	t_gquat q1; //q1 = vec_qnb[0];
	q1 = p->second;  p++;
	t_gquat q2; //q2 = vec_qnb[1];
	q2 = p->second;
	if (_beg_end)
	{
		q1 = q1 * ((nts - dt) / nts);
		q2 = q2 * (dt / nts);
		t = t - dt;
	}
	else
	{
		dt = out_nts - dt;
		q1 = q1 * (dt / nts);
		q2 = q2 * ((nts - dt) / nts);
		t = t + dt;
	}
	Vector3d vn1 = vn - an * dt;
	pos = pos - eth.v2dp(vn + vn1, dt / 2);	vn = vn1;
	//qnb = t_gbase::rv2q(-eth.wnin*nts)*qnb*t_gbase::rv2q(imu.phim);

	qnb = q1 + q2;
	qnb.normlize(qnb);
	Cnb = t_gbase::q2mat(qnb); att = t_gbase::m2att(Cnb); Cbn = Cnb.transpose(); vb = Cbn * vn;
	eth.Update(pos, vn); pos_ecef = Geod2Cart(pos, false);
	Ceb = eth.Cen*Cnb; Cbe = Ceb.transpose(); qeb = t_gbase::m2qua(Ceb); ve = eth.Cen*vn; ae = eth.Cen*an;
}

void gsins::t_gsins::_debug_ins_info()
{
	cout << "ins info : \n"
		<< "     time: " << setw(16) << t << "\n"
		<< "     att: " << setw(10) << att.transpose() / glv.deg << "\n"
		<< "     vel: " << setw(10) << ve.transpose() << "\n"
		<< "     pos: " << setw(10) << pos_ecef.transpose() << "\n"
		<< "     an: " << setw(10) << an.transpose() << "\n"
		<< "     ae: " << setw(10) << ae.transpose() << "\n"
		<< "     web: " << setw(10) << web.transpose() / glv.dps << "\n"
		<< "     wnb: " << setw(10) << wnb.transpose() / glv.dps << "\n"
		<< endl;
}


/******************************************* t_gsinskf ***********************************************/
gsins::t_gsinskf::t_gsinskf() :t_gkalman()
{
}

gsins::t_gsinskf::t_gsinskf(gnut::t_gsetbase * gset,string* filename) : t_gkalman(gset)
{
	this->sins = t_gsins(gset);
	this->gset(gset);
	if (filename) _ofilename = *filename;
	this->set_out();
	wmm = vmm = Eigen::Vector3d::Zero();
	_shm = t_scheme(gset);
	if (_shm.align)_aligned = false;
	else _aligned = true;
	_align_count = 0;
	_first_align = true;
	Eigen::VectorXd r = Eigen::VectorXd::Zero(5);
	//r << 0.1, 0.01, 1 * glv.deg;
	// an vn wnb fbib wbib
	r << 0.1, 0.1, 1 * glv.deg, 0.1, 0.1 * glv.dps;
	_avar.init(5, dynamic_cast<t_gsetins*>(_setkf)->ts(), r, Eigen::VectorXd::Ones(5));
}

gsins::t_gsinskf::~t_gsinskf()
{
	if (_fins)delete _fins; _fins = nullptr;
	if (_fkf)delete _fkf; _fkf = nullptr;
	if (_fpk)delete _fpk; _fpk = nullptr;
}


bool gsins::t_gsinskf::align_coarse(const vector<Eigen::Vector3d>& wm, const vector<Eigen::Vector3d>& vm, const t_scheme & scm)
{
	gtrace("t_gsinskf::align_coarse");
	sins.t = scm.t;
	for (int j = 0; j < scm.nSamples; j++) {
		wmm = wmm + wm[j];
		vmm = vmm + vm[j];
		_align_count++;
	}
	if (abs(_align_count - scm.align_time * scm.freq) < 1e-5)
	{
		wmm = wmm * fltflag / _align_count; vmm = vmm / _align_count;
		sins.qnb = t_gbase::a2qua(sins.align_coarse(wmm, vmm));
		return true;
	}
	return false;
}

bool gsins::t_gsinskf::align_pva(const Eigen::Vector3d& pos, const t_scheme & scm)
{
	gtrace("t_gsinskf::align_pva");

	if (_first_align == true)
	{
		//if amb fix
		// first ecef pos
		_first_pos = Geod2Cart(pos, false);
		_first_align = false;
		return false;
	}

	Vector3d endpos = Geod2Cart(pos, false);

	Vector3d baseline = XYZ2ENU(endpos, _first_pos);
	double dist = SQRT(SQR(baseline(0)) + SQR(baseline(1)));

	//if amb fix
	if (dist > 3)
	{
		double yaw = acos(fabs(baseline(1)) / dist);
		if (baseline(1) > 0 && baseline(0) > 0)yaw = yaw;
		if (baseline(1) > 0 && baseline(0) < 0)yaw = -yaw;
		if (baseline(1) < 0 && baseline(0) > 0)yaw = glv.PI - yaw;
		if (baseline(1) < 0 && baseline(0) < 0)yaw = -(glv.PI - yaw);
		sins.qnb = t_gbase::a2qua(Vector3d(0, 0, yaw));
		return true;
	}

	// if amb float
	if (dist > 10)
	{
	}


	return false;
}

bool gsins::t_gsinskf::align_vva(const Eigen::Vector3d& vel, const t_scheme & scm)
{
	gtrace("t_gsinskf::align_vva");
	if (vel.norm() > 2.0)
	{
		double yaw = atan2(vel(0), vel(1));
		if (vel(1) > 0 && vel(0) > 0)yaw = yaw;
		if (vel(1) > 0 && vel(0) < 0)yaw = -yaw;
		if (vel(1) < 0 && vel(0) > 0)yaw = glv.PI - yaw;
		if (vel(1) < 0 && vel(0) < 0)yaw = -(glv.PI - yaw);
		sins.qnb = t_gbase::a2qua(Vector3d(0, 0, yaw));
		return true;
	}
	return false;
}

bool gsins::t_gsinskf::cascaded_align(Eigen::Vector3d& pos, Eigen::Vector3d& vel)
{
	// todo: modified 
	gtrace("t_gsinskf::fsiam");
	_mutex.lock();

	sins.imu.Update(_wm, _vm, _shm);
	//if (double_eq(pos.norm(), 0.0))return 0;
	Eigen::Vector3d blh = Cart2Geod(pos, false);
	Eigen::Vector3d vn = t_gbase::Cen(blh).transpose()*vel;
	sins.set_posvel(blh, vn);
	sins.eth.Update(blh, vn);
	MOTION_TYPE state = motion_state();
	state = m_static;
	bool ok = false;
	if (state == m_static)
	{
		ok = align_coarse(_wm, _vm, _shm);
	}
	else if (state == m_straight)
	{
		if (SQRT(SQR(vn(0)) + SQR(vn(1))) > 2)ok = align_vva(vn, _shm);
		else ok = align_pva(pos, _shm);
	}
	else
	{
		_first_align = true;
	}

	_mutex.unlock();
	return ok;
}

bool gsins::t_gsinskf::align_moving_antenna(const Eigen::Vector3d& pos,
	const vector<Eigen::Vector3d>& lb,
	const vector<Eigen::Vector3d>& le,
	const t_scheme& scm)
{
	gtrace("t_gsinkf::align_moving_antenna");
	_mutex.lock();

	sins.imu.Update(_wm, _vm, _shm);
	if (lb.size() != le.size() || lb.size() == 0) return false;

	// *TEMPORARILY* : using only the first dl data

	gtrace("t_gsinskf::align_coarse");
	sins.t = scm.t;
	for (int j = 0; j < scm.nSamples; j++) {
		wmm = wmm + _wm[j];
		vmm = vmm + _vm[j];
		_align_count++;
	}
	if (le.front().norm()>0.05 && abs(_align_count - scm.align_time * scm.freq) >= 0)
	{
		wmm = wmm / _align_count; vmm = vmm / _align_count;

		sins.eth.Update(pos, Vector3d::Zero());

		Vector3d vb2 = lb.front();
		Vector3d vr2 = sins.eth.Cne * le.front();
		Vector3d vb1 = vmm.normalized();
		Vector3d vr1 = Vector3d(0, 0, 1);

		Matrix3d Mr = Matrix3d::Zero();
		Mr.block(0, 0, 1, 3) = vr1.transpose();
		Mr.block(1, 0, 1, 3) = (vr1.cross(vr2)).normalized().transpose();
		Mr.block(2, 0, 1, 3) = ((vr1.cross(vr2)).cross(vr1)).normalized().transpose();

		Matrix3d Mb = Matrix3d::Zero();
		Mb.block(0, 0, 1, 3) = vb1.transpose();
		Mb.block(1, 0, 1, 3) = (vb1.cross(vb2)).normalized().transpose();
		Mb.block(2, 0, 1, 3) = ((vb1.cross(vb2)).cross(vb1)).normalized().transpose();

		Matrix3d Cnb = Mr.inverse()*Mb;
		Vector3d att = t_gbase::m2att(Cnb);

		sins.qnb = t_gbase::a2qua(t_gbase::m2att(Cnb));

#ifdef  DEBUG
		std::cerr << lb.front() << std::endl;
		std::cerr << le.front() << std::endl;
		std::cerr << vb1.transpose() << std::endl;
		std::cerr << vr1.transpose() << std::endl;
		std::cerr << vb2.transpose() << std::endl;
		std::cerr << vr2.transpose() << std::endl;
		std::cerr << Mr.determinant() << std::endl;
		std::cerr << Mb.determinant() << std::endl;
		std::cerr << Mr << std::endl;
		std::cerr << Mb << std::endl;
		std::cerr << "att: " << (att*57.3).transpose() << std::endl;
		std::cerr << Cnb.determinant() << std::endl;
#endif //  DEBUG

		_mutex.unlock();
		return true;
	}
	_mutex.unlock();
	return false;
}


/************************************************************
** Function: Integration Continuous time state vector coefficient matrix
** Ft = [ Faa    Fav    Fap    -Cnb	   O33
**        Fva    Fvv    Fvp     O33    Cnb
**        O33    Fpv    Fpp     O33    O33
**        zeros(6,15) ];
** See also P159 - P161
***********************************************************/
void gsins::t_gsinskf::set_Ft(void)
{
	sins.err_mat();
	Ft.block(0, 0, 3, 3) = sins.Faa; Ft.block(0, 3, 3, 3) = sins.Fav; Ft.block(0, 6, 3, 3) = sins.Fap; Ft.block(0, 9, 3, 3) = -sins.Cnb;
	Ft.block(3, 0, 3, 3) = sins.Fva; Ft.block(3, 3, 3, 3) = sins.Fvv; Ft.block(3, 6, 3, 3) = sins.Fvp; Ft.block(3, 12, 3, 3) = sins.Cnb;
	Ft.block(6, 3, 3, 3) = sins.Fpv; Ft.block(6, 6, 3, 3) = sins.Fpp;
	Ft.block(9, 9, 3, 3) = sins._tauG.array().matrix().asDiagonal();
	Ft.block(12, 12, 3, 3) = sins._tauA.array().matrix().asDiagonal();
}

void gsins::t_gsinskf::set_Ft1()
{
	sins.err_mat();
	Ft.block(0, 0, 3, 3) = -t_gbase::askew(sins.eth.weie); Ft.block(0, 6, 3, 3) = Eigen::Matrix3d::Zero(); Ft.block(0, 9, 3, 3) = -sins.Ceb;
	Ft.block(3, 0, 3, 3) = t_gbase::askew(sins.eth.Cen*sins.fn); Ft.block(3, 3, 3, 3) = -2 * t_gbase::askew(sins.eth.weie); Ft.block(3, 12, 3, 3) = sins.Ceb;;
	Ft.block(6, 3, 3, 3) = Eigen::Matrix3d::Identity();
	Ft.block(9, 9, 3, 3) = sins._tauG.array().matrix().asDiagonal();
	Ft.block(12, 12, 3, 3) = sins._tauA.array().matrix().asDiagonal();
}

void gsins::t_gsinskf::set_meas_ZUPT(double t)
{
	MeasVel = sins.ve;
	tmeas = t;
	Flag = ZUPT_MEAS;
	_Cov_MeasVn << 0.0005, 0.0005, 0.0005;
	Rt << 0.5, 0.5, 0.5;
	Rt = Rt.array().abs2();
}

void gsins::t_gsinskf::set_meas_ZUPT_POS(double t)
{

	MeasPos = sins.pos_ecef;
	tmeas = t;
	Flag = ZUPT_POS_MEAS;
	_Cov_MeasPos << 0.000005, 0.000005, 0.000005;
	Rt << 0.5, 0.5, 0.5;
	Rt = Rt.array().abs2();
}

void gsins::t_gsinskf::set_meas_NHC(double t)
{
	MeasVel = sins.vb;
	tmeas = t;
	Flag = NHC_MEAS;
	_Cov_MeasVn << 0.1, 0.1, 0.1;
	Rt << 0.5, 0.5, 0.5, 1, 1, 1;
	Rt = Rt.array().abs2();
}

void gsins::t_gsinskf::set_meas_ATT(double t)
{
	tmeas = t;
	Flag = ATT_MEAS;
	_Cov_MeasAtt << 1, 1, 1;
	Rt << 0.5, 0.5, 0.5;// , 1, 1, 1, 1, 1, 1;  ///rn大小？
	Rt = Rt.array().abs2();

	Eigen::Vector3d mmm = Eigen::Vector3d::Zero(3), vmm = Eigen::Vector3d::Zero(3);
	double nts = _shm.nSamples * _shm.ts;
	for (auto mm : _mm)mmm = mmm + mm;
	for (auto vm : _vm)vmm = vmm + vm;
	Eigen::Vector3d m = mmm / nts;
	Eigen::Vector3d f = vmm / nts;
	double fai, theta, a;  /// roll pitch yaw
	fai = atan2(f[1], f[2]);   //   atan(fy/fz)
	theta = atan2(f[0], sqrt(SQR(f[1]) + SQR(f[2])  )  );
	m[0] = m[0] * cos(theta) + m[1] * sin(fai) * sin(theta) + m[2] * cos(fai) * sin(theta);
	m[1] = m[1] * cos(fai) - m[2] * sin(fai);
	a = atan2(m[0], m[1]);
	MeasAtt = Eigen::Vector3d(a, theta, fai);


}

void gsins::t_gsinskf::set_yaw(const double & yaw, double t)
{
	MeasYaw = yaw;
	tmeas = t;
	Flag = YAW_MEAS;
}

MOTION_TYPE gsins::t_gsinskf::motion_state()
{
	Eigen::Vector3d wmm = Eigen::Vector3d::Zero(3), vmm = Eigen::Vector3d::Zero(3);
	double nts = _shm.nSamples*_shm.ts;
	for (auto wm : _wm)wmm = wmm + wm;
	for (auto vm : _vm)vmm = vmm + vm;
	Eigen::Vector3d wbib = wmm / nts;
	Eigen::Vector3d fbib = vmm / nts;
	Eigen::VectorXd r(5);
	r << sins.an.norm(), sins.vn.norm(), sins.wnb.norm(), fbib.norm(), wbib.norm();
	_avar.update(r);
	if (!_aligned)
	{
		if (wbib.norm() < 0.05*glv.dps /*&& fbib.norm() - sins.eth.gn.norm() < 0.1 */&& _avar(3) < 0.5&& _avar(4) < 0.5*glv.dps)
			return m_static;
		else if (wbib.norm() < 0.1*glv.dps&&_avar(4) < 0.5*glv.dps)
			return m_straight;
		else
			return m_default;
	}
	else
	{
		if (_avar(2) < 2 * glv.dps && sins.wnb.norm() < 1 * glv.dps && abs(sins.vb(1)) > 2)
			return m_straight;
		else if (sins.an.norm() < 0.2 && sins.wnb.norm() < 0.23 * glv.dps && sins.vn.norm() < 0.02)
			return m_static;
		else
			return m_default;
	}

}

MEAS_TYPE gsins::t_gsinskf::meas_state()
{
	MOTION_TYPE motion = motion_state();
	switch (motion)
	{
	case m_static:
		return ZUPT_POS_MEAS; break;
		//return ZUPT_MEAS; break;
	case m_straight:
		return NHC_MEAS; break;
	default:
		return NO_MEAS;break;
	}
}

void gsins::t_gsinskf::feedback(double kfts)
{
	try {
		t_gkalman::feedback(kfts);
		sins.qnb = sins.qnb - FBXk.block(0, 0, 3, 1);
		sins.vn = sins.vn - FBXk.block(3, 0, 3, 1);
		sins.pos = sins.pos - FBXk.block(6, 0, 3, 1);
		sins.eb = sins.eb + FBXk.block(9, 0, 3, 1);
		sins.db = sins.db + FBXk.block(12, 0, 3, 1);
	}
	catch (...)
	{
		cout << "t_gsinskf::feedback" + to_string(kftk) + ": " + "FeedBack Failed!" << endl;
	}
}

void gsins::t_gsinskf::write()
{
	ostringstream os;
	_prt_ins_kml();
	sins.prt_sins(os);
	_fins->write(os.str().c_str(), os.str().size());
	os.str("");

	if (_fkf)
	{
		prt_Xk(kftk, Xk, os, true);
		_fkf->write(os.str().c_str(), os.str().size());
		os.str("");
	}
	if (_fpk)
	{
		prt_Pk(kftk, Pk, os);
		_fpk->write(os.str().c_str(), os.str().size());
		os.str("");
	}
}

void gsins::t_gsinskf::write_tmp(bool _beg_end)
{
	// todo: 
	if (sins.vec_qnb.size() < 2)
	{
		sins.vec_qnb[sins.t] = sins.qnb;
		return;
	}

	auto p = sins.vec_qnb.begin();
	if (_beg_end)
	{
		/// maintain slide window
		//p = sins.vec_qnb.begin();
		sins.vec_qnb.erase(p);
		///find last qnb
		p = sins.vec_qnb.end();
		p--;
	}
	else
	{
		p = sins.vec_qnb.end();
		p--;
		sins.vec_qnb.erase(p);

		p = sins.vec_qnb.begin();
	}

	sins.vec_qnb[sins.t] = sins.qnb;

	double t1 = p->first;
	double t2 = sins.t;
	double out_nts = 1.0 / sins.output_freq;

	/// it's time to print
	if (_beg_end && fmod(t1, out_nts) > fmod(t2, out_nts) || !_beg_end && fmod(t1, out_nts) < fmod(t2, out_nts))
	{
		/// int_sec output
		if (sins.output_intsec)
		{
			auto tmp_sins = sins;
			tmp_sins.int_sec(_beg_end);
			ostringstream os;
			tmp_sins.prt_sins(os);
			_prt_ins_kml();
			_fins->write(os.str().c_str(), os.str().size());
			os.str("");

			if (_fkf)
			{
				prt_Xk(sins.t, Xk, os, true);
				_fkf->write(os.str().c_str(), os.str().size());
				os.str("");
			}
			if (_fpk)
			{
				prt_Pk(sins.t, Pk, os);
				_fpk->write(os.str().c_str(), os.str().size());
				os.str("");
			}
		}
		else
		{
			ostringstream os;
			sins.prt_sins(os);
			_fins->write(os.str().c_str(), os.str().size());
			os.str("");

			if (_fkf)
			{
				prt_Xk(sins.t, Xk, os, true);
				_fkf->write(os.str().c_str(), os.str().size());
				os.str("");
			}
			if (_fpk)
			{
				prt_Pk(sins.t, Pk, os);
				_fpk->write(os.str().c_str(), os.str().size());
				os.str("");
			}

		}
	}
}

void gsins::t_gsinskf::prt_Xk(double t, const Eigen::VectorXd & Xk, ostringstream & os, bool IsError)
{
	if (IsError)
	{
		Eigen::VectorXd xk = Xk;
		for (int j = 0; j < 3; j++)
		{
			xk(j) = xk(j) / glv.deg * 60;
			xk(9 + j) = xk(9 + j) / glv.dph;
			xk(12 + j) = xk(12 + j) / glv.ug;
		}
		xk(6) = xk(6) / glv.deg*glv.hur;
		xk(7) = xk(7) / glv.deg*glv.hur;
		os << fixed << setprecision(6) << t;

		for (int j = 0; j < Xk.rows(); j++)
		{
			os << fixed << setprecision(8) << " " << xk(j) << " ";
		}
		os << endl;
	}
	else
	{
		Vector3d att1 = Xk.block(0, 0, 3, 1) / glv.deg;
		att1(2) = -att1(2);
		/*if (att1(2) > 0 && att1(2) < 180)att1(2) = 360 - att1(2);
		if (att1(2) > -180 && att1(2) < 0)att1(2) = -att1(2);*/
		Vector3d eb1 = Xk.block(9, 0, 3, 1) / glv.deg*glv.hur;
		Vector3d db1 = Xk.block(12, 0, 3, 1) / glv.mg;
		//Vector3d Car_pos = Geod2Cart(Xk.block(6, 0, 3, 1), false);
		Vector3d Car_pos = Xk.block(6, 0, 3, 1);
		//Eigen::Vector3d vn = Xk.block(3, 0, 3, 1);	
		//Eigen::Vector3d Ve = t_gbase::Cen(Xk.block(6, 0, 3, 1))*vn;
		Eigen::Vector3d Ve = Xk.block(3, 0, 3, 1);
		//os << fixed << setprecision(4) << t;
		//os << fixed << setprecision(8) <<
		//	" " << setw(13) << Car_pos(0) << " "
		//	" " << setw(13) << Car_pos(1) << " "
		//	" " << setw(13) << Car_pos(2) << " ";
		//os << fixed << setprecision(4) <<
		//	" " << setw(8) << Ve(0) << " "
		//	" " << setw(8) << Ve(1) << " "
		//	" " << setw(8) << Ve(2) << " ";
		//os << fixed << setprecision(4) <<
		//	" " << setw(8) << att1(0) << " "
		//	" " << setw(8) << att1(1) << " "
		//	" " << setw(10) << att1(2) << " ";
		//os << fixed << setprecision(8) <<
		//	" " << setw(13) << eb1(0) << " "
		//	" " << setw(13) << eb1(1) << " "
		//	" " << setw(13) << eb1(2) << " ";
		//os << fixed << setprecision(8) <<
		//	" " << setw(13) << db1(0) << " "
		//	" " << setw(13) << db1(1) << " "
		//	" " << setw(13) << db1(2) << " ";
		//os << endl;

		os << fixed << setprecision(6) << setw(18) << t;
		os << fixed << setprecision(3) <<
			setw(18) << Car_pos(0) <<
			setw(18) << Car_pos(1) <<
			setw(18) << Car_pos(2);
		os << fixed << setprecision(3) <<
			setw(10) << Ve(0) <<
			setw(10) << Ve(1) <<
			setw(10) << Ve(2);
		os << fixed << setprecision(4) <<
			setw(10) << att1(0) <<
			setw(10) << att1(1) <<
			setw(10) << att1(2);
		os << fixed << setprecision(4) <<
			setw(15) << eb1(0) <<
			setw(15) << eb1(1) <<
			setw(15) << eb1(2);
		os << fixed << setprecision(4) <<
			setw(15) << db1(0) <<
			setw(15) << db1(1) <<
			setw(15) << db1(2);
		os << endl;

		//Eigen::Vector3d att = Xk.block(0, 0, 3, 1) / glv.deg;
		//if (att(2) > 0 && att(2) < 180)att(2) = 360 - att(2);
		//if (att(2) > -180 && att(2) < 0)att(2) = -att(2);
		//Eigen::Vector3d eb = Xk.block(9, 0, 3, 1) / glv.deg*glv.hur;
		//Eigen::Vector3d db = Xk.block(12, 0, 3, 1) / glv.mg;
		//Eigen::Vector3d Car_pos = Geod2Cart(Xk.block(6, 0, 3, 1), false);
		//Eigen::Vector3d vn = Xk.block(3, 0, 3, 1);
		//Eigen::Vector3d Ve = t_gbase::Cen(Xk.block(6, 0, 3, 1))*vn;
		//os << fixed << setprecision(4) << t;
		//os << fixed << setprecision(8) <<
		//	" " << Car_pos(0) << " "
		//	" " << Car_pos(1) << " "
		//	" " << Car_pos(2) << " ";
		//os << fixed << setprecision(4) <<
		//	" " << Ve(0) << " "
		//	" " << Ve(1) << " "
		//	" " << Ve(2) << " ";
		//os << fixed << setprecision(4) <<
		//	" " << att(0) << " "
		//	" " << att(1) << " "
		//	" " << att(2) << " ";
		//os << fixed << setprecision(8) <<
		//	" " << eb(0) << " "
		//	" " << eb(1) << " "
		//	" " << eb(2) << " ";
		//os << fixed << setprecision(8) <<
		//	" " << db(0) << " "
		//	" " << db(1) << " "
		//	" " << db(2) << " ";
		//os << endl;
	}
}

void gsins::t_gsinskf::prt_Pk(double t, const Eigen::MatrixXd & Pk, ostringstream & os)
{
	os << fixed << setprecision(6) << t;
	if (Pk.rows() == Pk.cols())
	{
		for (int j = 0; j < Pk.rows(); j++)
		{
			os << scientific << setprecision(8) << " " << Pk(j, j) << " ";
		}
	}
	else
	{
		for (int j = 0; j < Pk.rows(); j++)
		{
			os << scientific << setprecision(8) << " " << Pk(j) << " ";
		}
	}
	os << endl;
}

void gsins::t_gsinskf::fusion(const Eigen::VectorXd & Xk1, const Eigen::VectorXd & Pk1, const Eigen::VectorXd & Xk2, const Eigen::VectorXd & Pk2, Eigen::VectorXd & Xk_Optimal, Eigen::VectorXd & Pk_Optimal)
{
	Pk_Optimal = Pk1 + Pk2;
	Xk_Optimal = (Pk2.cwiseProduct(Xk1) + Pk1.cwiseProduct(Xk2)).cwiseQuotient(Pk_Optimal);
	Pk_Optimal = Pk1.cwiseProduct(Pk2).cwiseQuotient(Pk_Optimal);
}

Eigen::Vector3d gsins::t_gsinskf::att_fusion(const Eigen::Vector3d & att1, const Eigen::Vector3d & att2, const Eigen::Vector3d & p1, const Eigen::Vector3d & p2)
{
	Eigen::Vector3d pf = p1.cwiseQuotient(p1 + p2);
	Eigen::Vector3d phi = t_gbase::a2qua(att1) - t_gbase::a2qua(att2);
	Eigen::Vector3d att = t_gbase::q2att(t_gbase::a2qua(att1) + phi.cwiseProduct(pf));
	return att;
}

void gsins::t_gsinskf::smoother()
{
	map<double, t_gsins>::iterator it1 = sins_forward.begin();
	map<double, t_gsins>::iterator it2 = sins_backward.begin();
	for (int i = 0; i < Xk_forward.size(); i++)
	{
		it2 = sins_backward.lower_bound(it1->first);
		//cout << it1->first << "  " << it2->first << endl;
		if (it2 == sins_backward.end()|| abs(it1->first - it2->first) > 1e-7)
		{
			it1++;
			if (it1 == sins_forward.end())break;
			continue;
		}
		Eigen::VectorXd Xk1 = Eigen::VectorXd::Zero(nq), Xk2 = Eigen::VectorXd::Zero(nq);
		Eigen::VectorXd Pk1 = Eigen::VectorXd::Zero(nq), Pk2 = Eigen::VectorXd::Zero(nq);
		Xk1 << it1->second.att, it1->second.ve, it1->second.pos_ecef, it1->second.eb, it1->second.db;
		Xk2 << it2->second.att, it2->second.ve, it2->second.pos_ecef, it2->second.eb, it2->second.db;
		//cerr << Xk_forward[it1->second.t] << endl;
		Xk1.block(0, 0, 3, 1) = t_gbase::q2att(t_gbase::a2qua(Xk1.block(0, 0, 3, 1)) - Xk_forward[it1->second.t].block(0, 0, 3, 1));
		Xk1.block(3, 0, 6, 1) = Xk1.block(3, 0, 6, 1) - Xk_forward[it1->second.t].block(3, 0, 6, 1);
		Xk1.block(9, 0, 6, 1) = Xk1.block(9, 0, 6, 1) + Xk_forward[it1->second.t].block(9, 0, 6, 1);
		Xk2.block(0, 0, 3, 1) = t_gbase::q2att(t_gbase::a2qua(Xk2.block(0, 0, 3, 1)) - Xk_backward[it2->second.t].block(0, 0, 3, 1));
		Xk2.block(3, 0, 6, 1) = Xk2.block(3, 0, 6, 1) - Xk_backward[it2->second.t].block(3, 0, 6, 1);
		Xk2.block(9, 0, 6, 1) = Xk2.block(9, 0, 6, 1) + Xk_backward[it2->second.t].block(9, 0, 6, 1);
		Xk2.block(3, 0, 3, 1) = -Xk2.block(3, 0, 3, 1); Xk2.block(9, 0, 3, 1) = -Xk2.block(9, 0, 3, 1);
		Pk1 = Pk_forward[it1->first].diagonal(); Pk2 = Pk_backward[it2->first].diagonal();
		Eigen::VectorXd Xk_Optimal, Pk_Optimal;
		fusion(Xk1, Pk1, Xk2, Pk2, Xk_Optimal, Pk_Optimal);
		Xk_Optimal.block(0, 0, 3, 1) = att_fusion(Xk1.block(0, 0, 3, 1), Xk2.block(0, 0, 3, 1), Pk1.block(0, 0, 3, 1), Pk2.block(0, 0, 3, 1));

		if (fmod(it1->first, 1.0 / sins.output_freq) < _shm.ts + (1e-7))
		{
			ostringstream os;
			prt_Xk(it1->first, Xk_Optimal, os, false);
			_fins->write(os.str().c_str(), os.str().size());
			os.str("");
		}

		//prt_Pk(it1->first, Pk_Optimal, os);
		//_fpk->write(os.str().c_str(), os.str().size());
		//os.str("");
		it1++;
	}
}

void gsins::t_gsinskf::write_tmp2(double t,const Eigen::VectorXd & Xk, const Eigen::VectorXd & Pk)
{
	// todo: 
	if (Xk_optimal.size() < 2)
	{
		Xk_optimal[t] = Xk;
		Pk_optimal[t] = Pk;
		return;
	}

	auto p = Xk_optimal.begin();
	Xk_optimal.erase(p);
	p = Pk_optimal.begin();
	Pk_optimal.erase(p);
	///find last qnb
	p = Xk_optimal.end();
	p--;
	Xk_optimal[t] = Xk;
	Pk_optimal[t] = Pk;

	double t1 = p->first;
	double t2 = t;
	double out_nts = 1.0 / sins.output_freq;

	/// it's time to print
	if (fmod(t1, out_nts) > fmod(t2, out_nts))
	{
		/// int_sec output
		if (sins.output_intsec)
		{
			auto tmp_sins = sins;
			tmp_sins.int_sec();
			ostringstream os;
			tmp_sins.prt_sins(os);
			_prt_ins_kml();
			_fins->write(os.str().c_str(), os.str().size());
			os.str("");

			if (_fkf)
			{
				prt_Xk(sins.t, Xk, os, true);
				_fkf->write(os.str().c_str(), os.str().size());
				os.str("");
			}
			if (_fpk)
			{
				prt_Pk(sins.t, Pk, os);
				_fpk->write(os.str().c_str(), os.str().size());
				os.str("");
			}
		}
		else
		{
			ostringstream os;
			sins.prt_sins(os);
			_fins->write(os.str().c_str(), os.str().size());
			os.str("");

			if (_fkf)
			{
				prt_Xk(sins.t, Xk, os, true);
				_fkf->write(os.str().c_str(), os.str().size());
				os.str("");
			}
			if (_fpk)
			{
				prt_Pk(sins.t, Pk, os);
				_fpk->write(os.str().c_str(), os.str().size());
				os.str("");
			}

		}
	}
}

void gsins::t_gsinskf::set_out()
{
	string tmp;
	tmp = dynamic_cast<t_gsetout*>(_setkf)->outputs("ins");
	if (tmp.empty())
	{
		tmp = _ofilename+"result.ins";
	}
	_fins = new t_giof;
	_fins->tsys(t_gtime::GPS);
	if (_ofilename != "") {
		substitute(tmp, "$(rec)", _ofilename, false);
	}
	_fins->mask(tmp);
	_fins->append(dynamic_cast<t_gsetout*>(_setkf)->append());

	tmp = dynamic_cast<t_gsetout*>(_setkf)->outputs("inskf");
	if (!tmp.empty())
	{
		_fkf = new t_giof;
		_fkf->tsys(t_gtime::GPS);
		if (_ofilename != "") {
			substitute(tmp, "$(rec)", _ofilename, false);
		}
		_fkf->mask(tmp);
		_fkf->append(dynamic_cast<t_gsetout*>(_setkf)->append());
	}else
	{
		_fkf = nullptr;
	}

	tmp = dynamic_cast<t_gsetout*>(_setkf)->outputs("inskfpk");
	if (!tmp.empty())
	{
		_fpk = new t_giof;
		_fpk->tsys(t_gtime::GPS);
		if (_ofilename != "") {
			substitute(tmp, "$(rec)", _ofilename, false);
		}
		_fpk->mask(tmp);
		_fpk->append(dynamic_cast<t_gsetout*>(_setkf)->append());
	}
	else
	{
		_fpk = nullptr;
	}


	ostringstream os;
	sins.prt_header(os);
	_fins->write(os.str().c_str(), os.str().size());

}

void gsins::t_gsinskf::_restore(const Eigen::VectorXd & X_store, const Eigen::MatrixXd & P_store)
{
	_Xk_store = X_store; _Pk_store = P_store;
}

void gsins::t_gsinskf::_recover(Eigen::VectorXd & X_recover, Eigen::MatrixXd & P_recover)
{
	X_recover = _Xk_store; P_recover = _Pk_store;
}

/******************************************* t_gmove_align ***********************************************/
t_gmove_align::t_gmove_align(const t_gsins& sins0, const Vector3d& lever0, double T0)
{
	sins = sins0; lever = lever0; T = T0;
	nt = 0;
	pos0 = sins.pos;
	dvk = LastPos = LastVn = Vector3d::Zero();
	LastVn = sins.vn;
	K = Matrix4d::Zero();
	Cne = Cee0 = Cnn0 = Cn0b0 = Cb0b = Matrix3d::Identity();
	Ce0n0 << 0, -sin(pos0(0)), cos(pos0(0)), 1, 0, 0, 0, cos(pos0(0)), sin(pos0(0));
}

void t_gmove_align::set_Meas(const Vector3d& vn, const Vector3d& pos, double t)
{
	MeasVel = vn; MeasPos = pos; tmeas = t;
}

/************************************************************
** Function: SINS Initial Alignment on moving base aided by External position and speed (GNSS)
** See also 联合单台GPS测速的捷联惯导单历元动基座粗对准
***********************************************************/
void t_gmove_align::align(const vector<Vector3d>& wm, const vector<Vector3d>& vm, const t_scheme& scm)
{
	//double nts = scm.ts, dt = sins.t - tmeas;
	double nts = scm.nSamples*scm.ts, dt = sins.t - tmeas;
	nt = nt + nts;
	sins.imu.Update(wm, vm, scm);
	sins.imu.phim = sins.Kg*sins.imu.phim - sins.eb*nts;
	sins.imu.dvbm = sins.Ka*sins.imu.dvbm - sins.db*nts;
	sins.wib = sins.imu.phim / nts;
	sins.fb = sins.imu.dvbm / nts;

	Vector3d IMUVn = MeasVel, IMUPos = MeasPos;
	if (nt>5.1) {
		IMUVn = MeasVel - sins.Cnb*t_gbase::askew(sins.web)*lever + sins.an*dt;
		IMUPos = MeasPos - sins.Fpv*sins.Cnb*lever - sins.eth.v2dp(sins.vn, dt);
	}

	sins.eth.Update(IMUPos, IMUVn);

	double Omega = glv.wie*nt;
	double dL = IMUPos(1) - pos0(1);

	Cee0 << cos(Omega), sin(Omega), 0, -sin(Omega), cos(Omega), 0, 0, 0, 1;
	Cne << -sin(dL), cos(dL), 0,
		-sin(IMUPos(0))*cos(dL), -sin(IMUPos(0))*sin(dL), cos(IMUPos(0)),
		cos(IMUPos(0))*cos(dL), cos(IMUPos(0))*sin(dL), sin(IMUPos(0));
	Cnn0 = Cne*Cee0*Ce0n0;

	Cb0b = Cb0b*t_gbase::rv2m(sins.imu.phim);
	dvk += sins.imu.dvbm;

	T_Alpha = Cb0b*dvk;
	T_Beta = (Cnn0.transpose())*(IMUVn - LastVn + T*(2 * sins.eth.wnie + sins.eth.wnen).cross(IMUVn) - T*sins.eth.gn);
	Matrix4d dK = ((t_gbase::m2m4(T_Beta) - t_gbase::m2m4_(T_Alpha)).transpose())*(t_gbase::m2m4(T_Beta) - t_gbase::m2m4_(T_Alpha))*(1 / T);
	K = K + dK;
	LastVn = IMUVn;
	dvk = Vector3d::Zero();

	EigenSolver<Matrix4d> es(K);
	Matrix4d val = es.pseudoEigenvalueMatrix();
	Matrix4d vec = es.pseudoEigenvectors();
	int minIndex = 0;
	for (int i = 0; i < 4; i++)
		if (val(minIndex, minIndex) > val(i, i))
		{
			minIndex = i;
		}
	qb0n0 = t_gquat(vec.block(0, minIndex, 4, 1));
	t_gquat::normlize(qb0n0);
	Cn0b0 = (t_gbase::q2mat(qb0n0)).transpose();

	sins.Cnb = Cnn0*Cn0b0*Cb0b;
	sins.qnb = t_gbase::m2qua(sins.Cnb);
	sins.fn = sins.qnb*sins.fb;
	sins.an = t_gbase::rv2q(-sins.eth.wnin*nts / 2)*sins.fn + sins.eth.gcc;
	Vector3d vn1 = sins.vn + sins.an*nts;  
	sins.pos = sins.pos + sins.eth.v2dp(sins.vn + vn1, nts / 2);    
	sins.vn = vn1;

	Vector3d angle = t_gbase::q2att(qb0n0);
	angle = angle / glv.deg;
	cout << fixed << setprecision(5) << (  nt) << ":";
	cout << angle(0) << " " << angle(1) << " " << angle(2) << endl;
}
