#include "gsimulation.h"
#include <iostream>
#include <random>
using namespace Eigen;
using namespace gsins;

#ifdef WIN32
#include <io.h>
#include <direct.h> 
#else
#include <unistd.h>
#include <sys/stat.h>
#endif

#ifdef WIN32
#define ACCESS(fileName,accessMode) _access(fileName,accessMode)
#define MKDIR(path) _mkdir(path)
#else
#define ACCESS(fileName,accessMode) access(fileName,accessMode)
#define MKDIR(path) mkdir(path,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif
t_gsim::t_gsim(string mark, t_gsetbase* set)
{
	_set  = set;
	_site = mark;
	Vector3d att = dynamic_cast<t_gsetins*>(set)->att();
	Vector3d vel = dynamic_cast<t_gsetins*>(set)->vel();
	Vector3d cart = dynamic_cast<t_gsetins*>(set)->pos();
	Vector3d pos = Cart2Geod(cart, false);
	double ts = 1.0 / dynamic_cast<t_gsetins*>(set)->freq();

	t_path path0 = { 0.0, att, vel, pos };
	_path.push_back(path0);
	double sb = sin(pos(0)), sb2 = sb * sb, sb4 = sb2 * sb2;
	_g = glv.g0*(1 + 5.27094e-3*sb2 + 2.32718e-5*sb4) - 3.086e-6*pos(2);
	_Vy = 0.0, _ts = ts, _T = 0;
	_MODE.insert(make_pair("m_keep", m_keep));
	_MODE.insert(make_pair("m_accelerate", m_accelerate));
	_MODE.insert(make_pair("m_decelerate", m_decelerate));
	_MODE.insert(make_pair("m_yawleft", m_yawleft));
	_MODE.insert(make_pair("m_yawright", m_yawright));
	_MODE.insert(make_pair("m_pitchup", m_pitchup));
	_MODE.insert(make_pair("m_pitchdown", m_pitchdown));
	_MODE.insert(make_pair("m_rollleft", m_rollleft));
	_MODE.insert(make_pair("m_rollright", m_rollright));
	_MODE.insert(make_pair("m_turnleft", m_turnleft));
	_MODE.insert(make_pair("m_turnright", m_turnright));
	_MODE.insert(make_pair("m_climb", m_climb));
	_MODE.insert(make_pair("m_descent", m_descent));
	_MODE.insert(make_pair("m_s", m_s));
	_MODE.insert(make_pair("m_8", m_8));

	set_out();

}
t_gsim::~t_gsim()
{
	if (_fins)delete _fins; _fins = nullptr;
}

void t_gsim::set_out()
{
	string tmp;
	tmp = dynamic_cast<t_gsetout*>(_set)->outputs("ins");
	if (tmp.empty())
	{
		tmp = "result.ins";
	}
	_fins = new t_giof;
	_fins->tsys(t_gtime::GPS);
	_fins->mask(tmp);
	_fins->append(dynamic_cast<t_gsetout*>(_set)->append());
}

t_gsim::t_gsim(const Vector3d& att, const Vector3d& vel, const Vector3d& pos, double ts)
{
	t_path path0 = { 0.0, att, vel, pos };
	_path.push_back(path0);
	double sb = sin(pos(0)), sb2 = sb * sb, sb4 = sb2 * sb2;
	_g = glv.g0*(1 + 5.27094e-3*sb2 + 2.32718e-5*sb4) - 3.086e-6*pos(2);
	_Vy = 0.0, _ts = ts, _T = 0;
	_MODE.insert(make_pair("m_keep", m_keep));
	_MODE.insert(make_pair("m_accelerate", m_accelerate));
	_MODE.insert(make_pair("m_decelerate", m_decelerate));
	_MODE.insert(make_pair("m_yawleft", m_yawleft));
	_MODE.insert(make_pair("m_yawright", m_yawright));
	_MODE.insert(make_pair("m_pitchup", m_pitchup));
	_MODE.insert(make_pair("m_pitchdown", m_pitchdown));
	_MODE.insert(make_pair("m_rollleft", m_rollleft));
	_MODE.insert(make_pair("m_rollright", m_rollright));
	_MODE.insert(make_pair("m_turnleft", m_turnleft));
	_MODE.insert(make_pair("m_turnright", m_turnright));
	_MODE.insert(make_pair("m_climb", m_climb));
	_MODE.insert(make_pair("m_descent", m_descent));
	_MODE.insert(make_pair("m_s", m_s));
	_MODE.insert(make_pair("m_8", m_8));
}

void t_gsim::trj_config(const string& filename)
{
	ifstream infile(filename);
	if (!infile.is_open())
	{
		cout << "Opening Trajectory Configure File is Failed" << endl;
	}
	string line, tmp1, tmp2, tmp3, tmp4, tmp5;
	string mode;
	stringstream buffer;
	double nT, w, a, T1;
	while (!infile.eof())
	{
		getline(infile, line);
		buffer << line;
		buffer >> mode; buffer >> nT; buffer >> w; buffer >> a; buffer >> T1;
		gen_mot(_MODE[mode], nT, w, a, T1);
		buffer.clear();
		buffer.str("");
	}
	for (int i = 0; i < _mot.size(); i++)_T += _mot[i].T;
}

void t_gsim::trj_config(const vector<string>& motion_strvec)
{
	try {
		string mode;
		double nT, w, a, T1;
		auto iter = motion_strvec.begin();
		while (iter != motion_strvec.end())
		{
			mode = *(iter++);
			transform(mode.begin(), mode.end(), mode.begin(), ::tolower);
			std::cout << mode << std::endl;
			nT = stod(*(iter++));
			w = stod(*(iter++));
			a = stod(*(iter++));
			T1 = stod(*(iter++));
			gen_mot(_MODE[mode], nT, w, a, T1);
		}
	}
	catch (...) {}

	for (int i = 0; i < _mot.size(); i++)_T += _mot[i].T;
}


void t_gsim::gen_mot(const MOTION_TYPE& mode, const double& T, const double& w, const double& a, const double& T1)
{
	double cf = w*glv.dps*_Vy;
	switch (mode)
	{
	case m_keep:
	{
		t_state mot_tmp;
		mot_tmp.T = T; mot_tmp.AngVel = Vector3d::Zero(); mot_tmp.Accel = Vector3d::Zero(); mot_tmp.Vy = _Vy;
		_mot.push_back(mot_tmp);
		break;
	}
	case m_accelerate:
	{
		t_state mot_tmp;
		mot_tmp.T = T; mot_tmp.AngVel = Vector3d::Zero(); mot_tmp.Accel = Vector3d(0, a, 0); mot_tmp.Vy = _Vy;
		_mot.push_back(mot_tmp);
		mot_tmp.Vy = _Vy + a*T; _Vy = mot_tmp.Vy;
		break;
	}
	case m_decelerate:
	{
		t_state mot_tmp;
		mot_tmp.T = T; mot_tmp.AngVel = Vector3d::Zero(); mot_tmp.Accel = Vector3d(0, -a, 0); mot_tmp.Vy = _Vy;
		_mot.push_back(mot_tmp);
		mot_tmp.Vy = _Vy - a*T; _Vy = mot_tmp.Vy;
		break;
	}
	case m_pitchup:
	{
		t_state mot_tmp;
		mot_tmp.T = T; mot_tmp.AngVel = Vector3d(w*glv.dps, 0.0, 0.0); mot_tmp.Accel = Vector3d(0, 0, cf); mot_tmp.Vy = _Vy;
		_mot.push_back(mot_tmp);
		break;
	}
	case m_pitchdown:
	{
		t_state mot_tmp;
		mot_tmp.T = T; mot_tmp.AngVel = Vector3d(-w*glv.dps, 0.0, 0.0); mot_tmp.Accel = Vector3d(0, 0, -cf); mot_tmp.Vy = _Vy;
		_mot.push_back(mot_tmp);
		break;
	}
	case m_rollleft:
	{
		t_state mot_tmp;
		mot_tmp.T = T; mot_tmp.AngVel = Vector3d(0.0, -w*glv.dps, 0.0); mot_tmp.Accel = Vector3d::Zero(); mot_tmp.Vy = _Vy;
		_mot.push_back(mot_tmp);
		break;
	}
	case m_rollright:
	{
		t_state mot_tmp;
		mot_tmp.T = T; mot_tmp.AngVel = Vector3d(0.0, w*glv.dps, 0.0); mot_tmp.Accel = Vector3d::Zero(); mot_tmp.Vy = _Vy;
		_mot.push_back(mot_tmp);
		break;
	}
	case m_yawleft:
	{
		t_state mot_tmp;
		mot_tmp.T = T; mot_tmp.AngVel = Vector3d(0.0, 0.0, w*glv.dps); mot_tmp.Accel = Vector3d(-cf, 0.0, 0.0); mot_tmp.Vy = _Vy;
		_mot.push_back(mot_tmp);
		break;
	}
	case m_yawright:
	{
		t_state mot_tmp;
		mot_tmp.T = T; mot_tmp.AngVel = Vector3d(0.0, 0.0, -w*glv.dps); mot_tmp.Accel = Vector3d(cf, 0.0, 0.0); mot_tmp.Vy = _Vy;
		_mot.push_back(mot_tmp);
		break;
	}
	case m_turnleft:
	{
		double rollw = atan2(cf, _g) / glv.dps / T1;
		gen_mot(m_rollleft, T1, rollw, 0, 0);
		gen_mot(m_yawleft, T, w, 0, 0);
		gen_mot(m_rollright, T1, rollw, 0, 0);
		break;
	}
	case m_turnright:
	{
		double rollw = atan2(cf, _g) / glv.dps / T1;
		gen_mot(m_rollright, T1, rollw, 0, 0);
		gen_mot(m_yawright, T, w, 0, 0);
		gen_mot(m_rollleft, T1, rollw, 0, 0);
		break;
	}
	case m_climb:
	{
		gen_mot(m_pitchup, T1, w, 0, 0);
		gen_mot(m_keep, T, 0, 0, 0);
		gen_mot(m_pitchdown, T1, w, 0, 0);
		break;
	}
	case m_descent:
	{
		gen_mot(m_pitchdown, T1, w, 0, 0);
		gen_mot(m_keep, T, 0, 0, 0);
		gen_mot(m_pitchup, T1, w, 0, 0);
		break;
	}
	case m_8:
	{
		double TT = 360 / w;
		gen_mot(m_turnleft, TT, w, 0, T1);
		gen_mot(m_turnright, TT, w, 0, T1);
		break;
	}
	default:
		break;
	}
}

void t_gsim::trj_imu()
{

	Vector3d att = _path[0].att, vn = _path[0].vn, pos = _path[0].pos;
	Vector3d phim, dvbm, wm, vm, wm_1 = Vector3d::Zero(), vm_1 = Vector3d::Zero();
	Vector3d vnr;
	Matrix3d Cbn = t_gbase::a2mat(att).transpose(), Cnb;
	double damping = 1 - exp(-_ts / 5.0), ts_2 = _ts / 2, t = 0;
	_path.erase(_path.begin());
	_eth.Update(pos, vn);
	int n = _T / _ts;
	for (int i = 0; i < _mot.size(); i++)
	{
		int nn = _mot[i].T / _ts, ufflag = 0;
		Vector3d wt = _mot[i].AngVel, at = _mot[i].Accel;
		if (wt.norm() == 0 && at.norm() == 0)
		{
			ufflag = 1;
			vnr = t_gbase::a2mat(att)*Vector3d(0, _mot[i].Vy, 0);
			ufflag = 0;
		}

		for (int j = 0; j < nn; j++)
		{
			Vector3d vn1, vn01, an;
			double sp = sin(att(0)), cp = cos(att(0)), sy = sin(att(2)), cy = cos(att(2));
			Matrix3d Cnt;
			Cnt << cy, -cp * sy, sp*sy,
				sy, cp*cy, -sp * cy,
				0, sp, cp;
			att = att + wt * _ts;
			Cnb = t_gbase::a2mat(att);
			att = t_gbase::m2att(Cnb); // ensure att is within the range
			if (ufflag == 1)
			{
				double na = (vn - vnr).norm() / _ts, maxa = 0.1;
				if (na < maxa)na = maxa;
				vn1 = vn - damping / na * (vn - vnr);
				an = (vn1 - vn) / _ts;
				vn01 = (vn + vn1) / 2; vn = vn1;
			}
			else
			{
				an = Cnt * at; vn1 = vn + an * _ts;
				vn01 = (vn + vn1) / 2; vn = vn1;
			}
			Vector3d dpos01 = Vector3d(vn01(1) / _eth.RMh, vn01(0) / _eth.cbRNh, vn01(2)) * ts_2;
			_eth.Update(pos + dpos01, vn01);
			pos = pos + dpos01 * 2;
			phim = t_gbase::q2rv(t_gbase::m2qua(Cbn*Cnb)) + (Cbn + Cnb.transpose())*(_eth.wnin*ts_2);
			wm = (Matrix3d::Identity() + t_gbase::askew(wm_1 / 12.0)).inverse()*phim;
			dvbm = Cbn * (t_gbase::rv2m(_eth.wnin*ts_2)*(an - _eth.gcc)*_ts);
			vm = (Matrix3d::Identity() + t_gbase::askew(wm / 2.0) + t_gbase::askew(wm_1 / 12.0)).inverse()*(dvbm - vm_1.cross(wm) / 12.0);
			wm_1 = wm; vm_1 = vm; Cbn = Cnb.transpose();
			t += _ts;

			t_imu_data data_tmp = { t, wm, vm };
			_data.push_back(data_tmp);
			t_path path_tmp = { t, att, vn, pos, Vector3d::Zero(), Vector3d::Zero()};
			_path.push_back(path_tmp);
		}
	}

	// set up IMU error/bias
	Vector3d eb, db, web, wdb;
	eb = dynamic_cast<t_gsetins*>(_set)->gyro_bias()*glv.deg / glv.hur;
	db = dynamic_cast<t_gsetins*>(_set)->acce_bias()*glv.mg;

	web = dynamic_cast<t_gsetign*>(_set)->misalignment_psd() * glv.dpsh / sqrt(_ts); /// rad/sqrt(s)
	wdb = dynamic_cast<t_gsetign*>(_set)->vel_psd() * glv.mgpsHz / sqrt(_ts);  /// m/s/sqrt(s)

	int len = _data.size();
	vector<vector<double>> whitenoise(6, vector<double>(len, 0.0));
	vector<vector<double>> markov(6, vector<double>(len, 0.0));
	for (int i = 0; i < 3; i++)
	{
		gen_whitenoise(0, web(i), len, whitenoise[  i]);
		gen_whitenoise(0, wdb(i), len, whitenoise[i+3]);
		gen_markov1(0.00000000 * sqrt(100.0 / 2), 100.0, _ts, len, markov[i  ]);
		gen_markov1(0.00000 * sqrt(100.0 / 2), 100.0, _ts, len, markov[i+3]);
	}
	
	for (int k = 0; k < len; k++)
	{
		for (int i = 0; i < 3; i++) /// add bias
		{
			_data[k].wm(i) += (eb(i) + markov[i    ][k] + whitenoise[i    ][k])* _ts; /// measurement
			_data[k].vm(i) += (db(i) + markov[i + 3][k] + whitenoise[i + 3][k])* _ts; /// measurement
			_path[k].eb(i) = eb(i) + markov[    i][k]; /// ground-truth
			_path[k].db(i) = db(i) + markov[i + 3][k]; /// ground-truth
		}
	}
}

void gsins::t_gsim::gen_path(t_gimudata * imu, t_gposdata * pos)
{
	vector<string> motion_strvec = dynamic_cast<t_gsetins*>(_set)->motion();
	trj_config(motion_strvec);
	trj_imu();
	t_gtime beg = dynamic_cast<t_gsetgen*>(_set)->beg();
	int out_freq = dynamic_cast<t_gsetins*>(_set)->out_freq();
	double t_beg = beg.sow();
	
	ostringstream os;

	Vector3d lever = dynamic_cast<t_gsetign*>(_set)->lever();
	for (int i = 0; i < _data.size(); i++)
	{
		double t = t_beg + _data[i].t;
		if (imu)
		{
			imu->add_IMU(t, _data[i].wm, _data[i].vm);
		}
		if (pos)
		{
			//double l_right = 0.5 * cos(2 * G_PI / 10.0*t);
			//double l_forward = 0.5 * sin(2 * G_PI / 10.0*t);
			//lever(0) = l_right;
			//lever(1) = l_forward;

			_eth.Update(_path[i].pos, _path[i].vn);
			Vector3d pos_gnss = Cart2Geod(Geod2Cart(_path[i].pos, false) + _eth.Cen*t_gbase::a2mat(_path[i].att)*lever, false);
			Vector3d vn_gnss = _path[i].vn; /// currently vn is not used in GNSS data simulation, so lever is (temporarily) not considered.
			pos->add_pos(t, pos_gnss, vn_gnss);
			if (_fins)
			{
				if (fmod(t - floor(t) + 1e-6, 1.0 / out_freq) > 1e-4) continue;
				os.clear(); os.str("");
				Vector3d att1 = _path[i].att / glv.deg;
				Vector3d eb1 = _path[i].eb / glv.deg*glv.hur;
				Vector3d db1 = _path[i].db / glv.mg;
				Vector3d Car_pos = Geod2Cart(_path[i].pos, false);
				Vector3d Geo_pos = Vector3d(_path[i].pos(0) / glv.deg, _path[i].pos(1) / glv.deg, _path[i].pos(2));
				Vector3d Ve = t_gbase::Cen(_path[i].pos)*(_path[i].vn);
				Vector3d vn = _path[i].vn;

				Vector3d pos_out, vel_out, att_out;

				auto order = dynamic_cast<t_gsetins*>(_set)->out_order();
				switch (order)
				{
					/// the att_out is positive if defination is from north to west[-Pi,Pi]
				case XYZ_XYZ_PRY_NW:
					pos_out = Car_pos;
					vel_out = Ve;
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

				os << fixed << setprecision(6) << setw(18) << t;
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

				_fins->write(os.str().c_str(), os.str().size());

			}
		}
	}

}

double gsins::t_gsim::gen_whitenoise(double mean, double stdvar)
{
	static std::default_random_engine engine;
	std::normal_distribution<double> distribution(mean, stdvar);
	double noise = distribution(engine);
	return noise;
}

void gsins::t_gsim::gen_whitenoise(double mean, double stdvar, int len, vector<double>& data)
{
	data.clear();
	static std::default_random_engine engine;
	std::normal_distribution<double> distribution(mean, stdvar);
	for (int i = 0; i < len; i++) data.push_back(distribution(engine));
}

void gsins::t_gsim::gen_markov1(double sigma, double tau, double ts, int len, vector<double>& data)
{
	double beta = 1.0 / tau;
	double q = 2 * SQR(sigma)*beta;
	double sQk = sqrt(q * ts);
	double a1 = exp(-beta * ts);
	int Ntao = ceil(tau / ts);
	int len1 = Ntao + len;
	vector<double> data_tmp=vector<double>(len1, 0.0);

	vector<double> noise;
	gen_whitenoise(0, 1.0, len1, noise);
	data_tmp[0] = sQk * noise[0];
	for (int i = 1; i < len1; i++)
	{
		data_tmp[i] = a1 * data_tmp[i - 1] + sQk * noise[i];
	}

	data.assign(data_tmp.begin() + Ntao, data_tmp.end());
}


void gsins::t_gsim::gen_path_from_groundtruth_test(const string& fgt, const vector<string>& fothers, t_gimudata *imudata, t_gposdata *posdata, int fs_gt, int fs_imu)
{
	struct ros_odom
	{
		long long stamp;
		Vector3d p;
		Quaterniond q;
		Vector3d v;
	};

	double ts_gt = 1.0 / fs_gt;
	double ts_imu = 1.0 / fs_imu;
	int nn = fs_gt / fs_imu;

	_data.clear();
	ifstream ifs_gt(fgt);
	vector<ros_odom> all_gt;

	string line;
	stringstream ss;
	while (ifs_gt.good())
	{
		getline(ifs_gt, line);
		ss.clear();
		ss.str("");
		ss << line;
		long long stamp;
		double px, py, pz;
		double qw, qx, qy, qz;
		double vx, vy, vz; // linear velocity
		double ax, ay, az; // angular velocity

		ss >> stamp >> px >> py >> pz >> qw >> qx >> qy >> qz
			>> vx >> vy >> vz >> ax >> ay >> az;
		while (all_gt.size() > 0 && (stamp - all_gt.back().stamp) / 1e9 > ts_gt + 1e-6)
		{
			double scale = ts_gt / ((stamp - all_gt.back().stamp) / 1e9);
			Vector3d p_inter;
			Quaterniond q_inter;
			Vector3d v_inter;
			p_inter = all_gt.back().p + (Vector3d(px, py, pz) - all_gt.back().p)*scale;
			v_inter = all_gt.back().v + (Vector3d(vx, vy, vz) - all_gt.back().v)*scale;
			q_inter = Quaterniond(all_gt.back().q.toRotationMatrix()
				* t_gbase::rv2m(t_gbase::q2rv(t_gbase::m2qua((all_gt.back().q.inverse()*Quaterniond(qw, qx, qy, qz)).toRotationMatrix()))*scale));
			all_gt.push_back({ stamp + (long long)(ts_gt * 1e9), p_inter, q_inter, v_inter });
		}
		all_gt.push_back({ stamp,Vector3d(px,py,pz),Quaterniond(qw,qx,qy,qz),Vector3d(vx,vy,vz) });
	}
	ifs_gt.close();

	if (all_gt.empty()) return;

	cerr << all_gt[0].stamp << std::endl;
	Vector3d att = t_gbase::m2att(all_gt[0].q.toRotationMatrix());
	Vector3d vn = all_gt[0].v;
	Vector3d pos = _path[0].pos;
	Vector3d phim, dvbm, wm, vm, wm_1 = Vector3d::Zero(), vm_1 = Vector3d::Zero();
	Vector3d vnr;
	Matrix3d Cbn = t_gbase::a2mat(att).transpose(), Cnb;
	// double damping = 1 - exp(-_ts / 5.0), ts_2 = _ts / 2, t = 0;
	double t = 0;

	_path.erase(_path.begin());
	_eth.Update(pos, vn);

	std::cerr << pos.transpose() << std::endl;
	// ATTENTION!: all_gt is refered to "n0"
	Matrix3d Cen0 = _eth.Cen;
	Vector3d pos0 = pos;

	int count = 0;
	Vector3d wm_accu = Vector3d::Zero();
	Vector3d vm_accu = Vector3d::Zero();

	for (int i = 1; i < all_gt.size(); i++)
	{
		count++;
		Vector3d an;
		double sp = sin(att(0)), cp = cos(att(0)), sy = sin(att(2)), cy = cos(att(2));
		Matrix3d Cnt;
		Cnt << cy, -cp * sy, sp*sy,
			sy, cp*cy, -sp * cy,
			0, sp, cp;

		pos = Cart2Geod(Geod2Cart(pos0, false) + Cen0 * all_gt[i].p, false);
		_eth.Update(pos, vn);

		vn = _eth.Cne * Cen0 * all_gt[i].v;
		an = _eth.Cne * Cen0 *(all_gt[i].v - all_gt[i - 1].v) / ts_gt;
		att = _eth.Cne*Cen0*t_gbase::m2att(all_gt[i].q.toRotationMatrix());
		Cnb = t_gbase::a2mat(att);

		phim = t_gbase::q2rv(t_gbase::m2qua(Cbn*Cnb)) + (Cbn + Cnb.transpose())*(_eth.wnin*(ts_gt / 2));
		wm = (Matrix3d::Identity() + t_gbase::askew(wm_1 / 12.0)).inverse()*phim;
		// std::cerr << _eth.gcc.transpose() << std::endl;
		dvbm = Cbn * (t_gbase::rv2m(_eth.wnin*(ts_gt / 2))*(an - _eth.gcc)*ts_gt);
		vm = (Matrix3d::Identity() + t_gbase::askew(wm / 2.0) + t_gbase::askew(wm_1 / 12.0)).inverse()*(dvbm - vm_1.cross(wm) / 12.0);
		wm_1 = wm; vm_1 = vm; Cbn = Cnb.transpose();
		t += ts_gt;

		wm_accu += wm;
		vm_accu += vm;
		if (count == nn)
		{
			t_imu_data data_tmp = { t, wm_accu, vm_accu };
			_data.push_back(data_tmp);
			t_path path_tmp = { t, att, vn, pos, Vector3d::Zero(), Vector3d::Zero() };
			_path.push_back(path_tmp);

			wm_accu = Vector3d::Zero();
			vm_accu = Vector3d::Zero();
			count = 0;
		}
	}


	// set up IMU error/bias
	Vector3d eb, db, web, wdb;
	eb = dynamic_cast<t_gsetins*>(_set)->gyro_bias()*glv.deg / glv.hur;
	db = dynamic_cast<t_gsetins*>(_set)->acce_bias()*glv.mg;

	web = dynamic_cast<t_gsetign*>(_set)->misalignment_psd() * glv.dpsh / sqrt(ts_imu); /// rad/sqrt(s)
	wdb = dynamic_cast<t_gsetign*>(_set)->vel_psd() * glv.mgpsHz / sqrt(ts_imu);  /// m/s/sqrt(s)

	int len = _data.size();
	vector<vector<double>> whitenoise(6, vector<double>(len, 0.0));
	vector<vector<double>> markov(6, vector<double>(len, 0.0));
	for (int i = 0; i < 3; i++)
	{
		gen_whitenoise(0, web(i), len, whitenoise[i]);
		gen_whitenoise(0, wdb(i), len, whitenoise[i + 3]);
		gen_markov1(0.00000000 * sqrt(100.0 / 2), 100.0, ts_imu, len, markov[i]);
		gen_markov1(0.00000 * sqrt(100.0 / 2), 100.0, ts_imu, len, markov[i + 3]);
	}

	for (int k = 0; k < len; k++)
	{
		for (int i = 0; i < 3; i++) /// add bias
		{
			_data[k].wm(i) += (eb(i) + markov[i][k] + whitenoise[i][k])* ts_imu; /// measurement
			_data[k].vm(i) += (db(i) + markov[i + 3][k] + whitenoise[i + 3][k])* ts_imu; /// measurement
			_path[k].eb(i) = eb(i) + markov[i][k]; /// ground-truth
			_path[k].db(i) = db(i) + markov[i + 3][k]; /// ground-truth
		}
	}



	t_gtime beg = dynamic_cast<t_gsetgen*>(_set)->beg();
	int out_freq = dynamic_cast<t_gsetins*>(_set)->out_freq();
	double t_beg = beg.sow();

	ostringstream os;

	Vector3d lever = dynamic_cast<t_gsetign*>(_set)->lever();
	for (int i = 0; i < _data.size(); i++)
	{
		double t = t_beg + _data[i].t;
		if (imudata)
		{
			imudata->add_IMU(t, _data[i].wm, _data[i].vm);
		}
		if (posdata)
		{
			_eth.Update(_path[i].pos, _path[i].vn);
			Vector3d pos_gnss = Cart2Geod(Geod2Cart(_path[i].pos, false) + _eth.Cen*t_gbase::a2mat(_path[i].att)*lever, false);
			Vector3d vn_gnss = _path[i].vn; /// currently vn is not used in GNSS data simulation, so lever is (temporarily) not considered.
			posdata->add_pos(t, pos_gnss, vn_gnss);
			if (_fins)
			{
				if (fmod(t - floor(t) + 1e-6, 1.0 / out_freq) > 1e-4) continue;
				os.clear(); os.str("");
				Vector3d att1 = _path[i].att / glv.deg;
				Vector3d eb1 = _path[i].eb / glv.deg*glv.hur;
				Vector3d db1 = _path[i].db / glv.mg;
				Vector3d Car_pos = Geod2Cart(_path[i].pos, false);
				Vector3d Geo_pos = Vector3d(_path[i].pos(0) / glv.deg, _path[i].pos(1) / glv.deg, _path[i].pos(2));
				Vector3d Ve = t_gbase::Cen(_path[i].pos)*(_path[i].vn);
				Vector3d vn = _path[i].vn;

				Vector3d pos_out, vel_out, att_out;

				auto order = dynamic_cast<t_gsetins*>(_set)->out_order();
				switch (order)
				{
					/// the att_out is positive if defination is from north to west[-Pi,Pi]
				case XYZ_XYZ_PRY_NW:
					pos_out = Car_pos;
					vel_out = Ve;
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

				os << fixed << setprecision(6) << setw(18) << t;
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

				_fins->write(os.str().c_str(), os.str().size());

			}
		}
	}

	for (const auto path : fothers)
	{
		ifstream ifs_other(path);
		ofstream ofs_other(path + "gpst");

		string line;
		while (ifs_other.good())
		{
			getline(ifs_other, line);
			if (line.size() < 2) continue;
			if (line[0] == '%' || line[0] == '#') continue;
			for (int i = 0; i < line.size(); i++) if (line[i] == ',') line[i] = ' ';
			if (line.find(' ') != string::npos)
			{
				long long stamp = stoll(line.substr(0, line.find(' ')));
				ofs_other << setprecision(4) << setiosflags(ios::fixed) << t_beg + (stamp - all_gt[0].stamp) / 1e9 << line.substr(line.find(' '), line.size() - 1) << std::endl;

			}
		}

		ifs_other.close();
		ofs_other.close();
	}
}