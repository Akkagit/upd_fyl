
#include"gsetvis.h"
using namespace gvis;

gvis::t_gsetvis::t_gsetvis() : t_gsetbase()
{
	_set.insert(XMLKEY_VIS);
}


gvis::t_gsetvis::~t_gsetvis()
{
}

void gvis::t_gsetvis::check()
{
	_gmutex.lock();
	_gmutex.unlock();
}

void gvis::t_gsetvis::help()
{
	_gmutex.lock();
	_gmutex.unlock();
}


double gvis::t_gsetvis::ts()
{
	_gmutex.lock();

	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("freq").as_double();
	res = 1 / res;
	_gmutex.unlock();
	return res;
}

double gvis::t_gsetvis::imu_ts()
{
	_gmutex.lock();

	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("imu_freq").as_double();
	res = 1 / res;
	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::freq()
{
	_gmutex.lock();

	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("freq").as_int();

	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::imu_freq()
{
	_gmutex.lock();

	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("imu_freq").as_int();

	_gmutex.unlock();
	return res;
}

gvis::Distortion_Type gvis::t_gsetvis::cam0_distortion_model()
{
	_gmutex.lock();

	string distortion_type = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("cam0_distortion_model").value();

	_gmutex.unlock();

	return str2distortion(distortion_type);
}

gvis::Distortion_Type gvis::t_gsetvis::cam1_distortion_model()
{
	_gmutex.lock();

	string distortion_type = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("cam1_distortion_model").value();
	_gmutex.unlock();

	return str2distortion(distortion_type);
}


Eigen::Vector2d gvis::t_gsetvis::cam0_resolution()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("cam0_resolution").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	double X, Y;
	ss >> X >> Y;
	_gmutex.unlock();
	return Eigen::Vector2d(X, Y);
}

Eigen::Vector2d gvis::t_gsetvis::cam1_resolution()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("cam1_resolution").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	double X, Y;
	ss >> X >> Y;
	_gmutex.unlock();
	return Eigen::Vector2d(X, Y);
}

Eigen::Vector4d gvis::t_gsetvis::cam0_intrinsics()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("cam0_intrinsics").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	double fx, fy, cx, cy;
	ss >> fx >> fy >> cx >> cy;
	_gmutex.unlock();
	return Eigen::Vector4d(fx, fy, cx, cy);
}

Eigen::Vector4d gvis::t_gsetvis::cam1_intrinsics()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("cam1_intrinsics").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	double fx, fy, cx, cy;
	ss >> fx >> fy >> cx >> cy;
	_gmutex.unlock();
	return Eigen::Vector4d(fx, fy, cx, cy);
}

Eigen::Vector4d gvis::t_gsetvis::cam0_distortion_coeffs()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("cam0_distortion_coeffs").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	double k1, k2, p1, p2;
	ss >> k1 >> k2 >> p1 >> p2;
	_gmutex.unlock();
	return Eigen::Vector4d(k1, k2, p1, p2);
}

Eigen::Vector4d gvis::t_gsetvis::cam1_distortion_coeffs()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("cam1_distortion_coeffs").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	double k1, k2, p1, p2;
	ss >> k1 >> k2 >> p1 >> p2;
	_gmutex.unlock();
	return Eigen::Vector4d(k1, k2, p1, p2);
}


Eigen::Matrix3d gvis::t_gsetvis::R_cam0_cam1()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("R_cam0_cam1").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	Eigen::Matrix3d res;
	double R11, R12, R13, R21, R22, R23, R31, R32, R33;
	ss >> R11 >> R12 >> R13 >> R21 >> R22 >> R23 >> R31 >> R32 >> R33;
	res << R11, R12, R13, R21, R22, R23, R31, R32, R33;
	_gmutex.unlock();
	return res;
}

Eigen::Vector3d gvis::t_gsetvis::t_cam0_cam1()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("t_cam0_cam1").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	Eigen::Matrix3d res;
	double X, Y, Z;
	ss >> X >> Y >> Z;
	_gmutex.unlock();
	return Eigen::Vector3d(X, Y, Z);
}

Eigen::Isometry3d gvis::t_gsetvis::T_cam0_cam1()
{
	_gmutex.lock();
	_gmutex.unlock(); Eigen::Matrix3d tmp_R_cam0_cam1 = R_cam0_cam1(); _gmutex.lock();
	_gmutex.unlock(); Eigen::Vector3d tmp_t_cam0_cam1 = t_cam0_cam1(); _gmutex.lock();
	
	Eigen::Isometry3d tmp_T_cam0_cam1=Eigen::Isometry3d::Identity();
	tmp_T_cam0_cam1.linear() = tmp_R_cam0_cam1;
	tmp_T_cam0_cam1.translation() = tmp_t_cam0_cam1;
	_gmutex.unlock();
	return tmp_T_cam0_cam1;
}

Eigen::Matrix3d gvis::t_gsetvis::R_cam0_imu()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("R_cam0_imu").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	Eigen::Matrix3d res;
	double R11, R12, R13, R21, R22, R23, R31, R32, R33;
	ss >> R11 >> R12 >> R13 >> R21 >> R22 >> R23 >> R31 >> R32 >> R33;
	res << R11, R12, R13, R21, R22, R23, R31, R32, R33;
	_gmutex.unlock();
	return res;
}

Eigen::Vector3d gvis::t_gsetvis::t_cam0_imu()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("t_cam0_imu").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	Eigen::Matrix3d res;
	double X, Y, Z;
	ss >> X >> Y >> Z;
	_gmutex.unlock();
	return Eigen::Vector3d(X, Y, Z);
}

Eigen::Isometry3d gvis::t_gsetvis::T_cam0_imu()
{
    _gmutex.lock();
	_gmutex.unlock(); Eigen::Matrix3d tmp_R_cam0_imu = R_cam0_imu(); _gmutex.lock();
	_gmutex.unlock(); Eigen::Vector3d tmp_t_cam0_imu = t_cam0_imu(); _gmutex.lock();
	
	Eigen::Isometry3d tmp_T_cam0_imu=Eigen::Isometry3d::Identity();
	tmp_T_cam0_imu.linear() = tmp_R_cam0_imu;
	tmp_T_cam0_imu.translation() = tmp_t_cam0_imu;
	_gmutex.unlock();
	return tmp_T_cam0_imu;
}



int gvis::t_gsetvis::num_of_cam()
{
	_gmutex.lock();

	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("num_of_cam").as_int();

	_gmutex.unlock();

	return res;
}

int gvis::t_gsetvis::grid_row()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("grid_row").as_int();
	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::grid_col()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("grid_col").as_int();
	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::grid_min_feature_num()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("grid_min_feature_num").as_int();
	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::grid_max_feature_num()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("grid_max_feature_num").as_int();
	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::pyramid_levels()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("pyramid_levels").as_int();
	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::patch_size()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("patch_size").as_int();
	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::fast_threshold()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("fast_threshold").as_int();
	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::fscalefactor()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("fscalefactor").as_int();
	_gmutex.unlock();
	return res;
}

double gvis::t_gsetvis::ransac_threshold()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("ransac_threshold").as_double();
	_gmutex.unlock();
	return res;
}

double gvis::t_gsetvis::stereo_threshold()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("stereo_threshold").as_double();
	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::max_iteration()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("max_iteration").as_int();
	_gmutex.unlock();
	return res;
}

double gvis::t_gsetvis::track_precision()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("track_precision").as_double();
	_gmutex.unlock();
	return res;
}

double gvis::t_gsetvis::mbf()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("mbf").as_double();
	_gmutex.unlock();
	return res;
}

bool gvis::t_gsetvis::usingstereorecify()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("usingstereorectify").as_int();
	bool b_res;
	//cout << res << endl;
	if (res == 1) b_res = true;
	else b_res = false;
	_gmutex.unlock();
	return b_res;
}

Eigen::Matrix3d gvis::t_gsetvis::K_L()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("K_L").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	Eigen::Matrix3d res;
	double R11, R12, R13, R21, R22, R23, R31, R32, R33;
	ss >> R11 >> R12 >> R13 >> R21 >> R22 >> R23 >> R31 >> R32 >> R33;
	res << R11, R12, R13, R21, R22, R23, R31, R32, R33;
	_gmutex.unlock();
	return res;
}

Eigen::MatrixXd gvis::t_gsetvis::D_L()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("D_L").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	Eigen::MatrixXd res;
	res = Eigen::MatrixXd::Zero(1, 5);
	double R1, R2, R3, R4, R5;
	ss >> R1 >> R2 >> R3 >> R4 >> R5;
	res << R1, R2, R3, R4, R5;
	_gmutex.unlock();

	return res;
}

Eigen::Matrix3d gvis::t_gsetvis::R_L()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("R_L").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	Eigen::Matrix3d res;
	double R11, R12, R13, R21, R22, R23, R31, R32, R33;
	ss >> R11 >> R12 >> R13 >> R21 >> R22 >> R23 >> R31 >> R32 >> R33;
	res << R11, R12, R13, R21, R22, R23, R31, R32, R33;
	_gmutex.unlock();
	return res;
}

Eigen::MatrixXd gvis::t_gsetvis::P_L()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("P_L").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	Eigen::MatrixXd res;
	res = Eigen::MatrixXd::Zero(3, 4);
	double R11, R12, R13, R14, R21, R22, R23, R24, R31, R32, R33, R34;
	ss >> R11 >> R12 >> R13 >> R14 >> R21 >> R22 >> R23 >> R24 >> R31 >> R32 >> R33 >> R34;
	res << R11, R12, R13, R14, R21, R22, R23, R24, R31, R32, R33, R34;
	_gmutex.unlock();
	return res;
}

Eigen::Matrix3d gvis::t_gsetvis::K_R()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("K_R").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	Eigen::Matrix3d res;
	double R11, R12, R13, R21, R22, R23, R31, R32, R33;
	ss >> R11 >> R12 >> R13 >> R21 >> R22 >> R23 >> R31 >> R32 >> R33;
	res << R11, R12, R13, R21, R22, R23, R31, R32, R33;
	_gmutex.unlock();
	return res;
}

Eigen::MatrixXd gvis::t_gsetvis::D_R()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("D_R").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	Eigen::MatrixXd res;
	res = Eigen::MatrixXd::Zero(1, 5);
	double R1, R2, R3, R4, R5;
	ss >> R1 >> R2 >> R3 >> R4 >> R5;
	res << R1, R2, R3, R4, R5;
	_gmutex.unlock();
	return res;
}

Eigen::Matrix3d gvis::t_gsetvis::R_R()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("R_R").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	Eigen::Matrix3d res;
	double R11, R12, R13, R21, R22, R23, R31, R32, R33;
	ss >> R11 >> R12 >> R13 >> R21 >> R22 >> R23 >> R31 >> R32 >> R33;
	res << R11, R12, R13, R21, R22, R23, R31, R32, R33;
	_gmutex.unlock();
	return res;
}

Eigen::MatrixXd gvis::t_gsetvis::P_R()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("P_R").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	Eigen::MatrixXd res;
	res = Eigen::MatrixXd::Zero(3, 4);
	double R11, R12, R13, R14, R21, R22, R23, R24, R31, R32, R33, R34;
	ss >> R11 >> R12 >> R13 >> R14 >> R21 >> R22 >> R23 >> R24 >> R31 >> R32 >> R33 >> R34;
	res << R11, R12, R13, R14, R21, R22, R23, R24, R31, R32, R33, R34;
	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::nfeatures()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("nfeature").as_int();
	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::max_fast()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("fmaxthfast").as_int();
	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::min_fast()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("fminthfast").as_int();
	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::max_cnt()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("maxcnt").as_int();
	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::max_cam_state_size()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("max_cam_state_size").as_int();
	_gmutex.unlock();
	return res;
}

double gvis::t_gsetvis::position_std_threshold()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("position_std_threshold").as_double();
	_gmutex.unlock();
	return res;
}

double gvis::t_gsetvis::rotation_threshold()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("rotation_threshold").as_double();
	_gmutex.unlock();
	return res;
}

double gvis::t_gsetvis::translation_threshold()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("translation_threshold").as_double();
	_gmutex.unlock();
	return res;
}

double gvis::t_gsetvis::feature_observation_noise()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("feature_observation_noise").as_double();
	_gmutex.unlock();
	return res;
}

double gvis::t_gsetvis::huber_epsilon()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("huber_epsilon").as_double();
	_gmutex.unlock();
	return res;
}

double gvis::t_gsetvis::estimation_precision()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("estimation_precision").as_double();
	_gmutex.unlock();
	return res;
}

double gvis::t_gsetvis::initial_damping()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("initial_damping").as_double();
	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::outler_loop_max_iteration()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("outler_loop_max_iteration").as_int();
	_gmutex.unlock();
	return res;
}

int gvis::t_gsetvis::inner_loop_max_iteration()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("inner_loop_max_iteration").as_int();
	_gmutex.unlock();
	return res;
}

bool gvis::t_gsetvis::estimate_extrinsic()
{
	_gmutex.lock();
	bool res = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("estimate_extrinsic").as_bool();
	_gmutex.unlock();
	return res;
}

Eigen::Vector3d gvis::t_gsetvis::initial_extrinsic_rotation_cov()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("initial_extrinsic_rotation_cov").value();
	if (str.empty()) return Eigen::Vector3d(3.0462e-4, 3.0462e-4, 3.0462e-4);
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	Eigen::Vector3d res;
	double s1,s2,s3;
	ss >> s1 >> s2 >> s3;
	res << s1, s2, s3;
	_gmutex.unlock();
	return res;
}

Eigen::Vector3d gvis::t_gsetvis::initial_extrinsic_translation_cov()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_VIS).attribute("initial_extrinsic_translation_cov").value();
	if (str.empty()) return Eigen::Vector3d(1e-4, 1e-4, 1e-4);
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	Eigen::Vector3d res;
	double s1, s2, s3;
	ss >> s1 >> s2 >> s3;
	res << s1, s2, s3;
	_gmutex.unlock();
	return res;
}

