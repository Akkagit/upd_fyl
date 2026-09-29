#include"gsetlidar.h"
using namespace glidar;

glidar::t_gsetlidar::t_gsetlidar() : t_gsetbase()
{
	_set.insert(XMLKEY_LIDAR);
}


glidar::t_gsetlidar::~t_gsetlidar()
{
}

void glidar::t_gsetlidar::check()
{
	_gmutex.lock();
	_gmutex.unlock();
}

void glidar::t_gsetlidar::help()
{
	_gmutex.lock();
	_gmutex.unlock();
}


double glidar::t_gsetlidar::ts()
{
	_gmutex.lock();

	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("freq").as_double();
	res = 1 / res;
	_gmutex.unlock();
	return res;
}

double glidar::t_gsetlidar::imu_ts()
{
	_gmutex.lock();

	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("imu_freq").as_double();
	res = 1 / res;
	_gmutex.unlock();
	return res;
}

int glidar::t_gsetlidar::freq()
{
	_gmutex.lock();

	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("freq").as_int();

	_gmutex.unlock();
	return res;
}

int glidar::t_gsetlidar::imu_freq()
{
	_gmutex.lock();

	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("imu_freq").as_int();

	_gmutex.unlock();
	return res;
}

Eigen::Matrix3d glidar::t_gsetlidar::R_lidar_imu()
{
	//_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("R_lidar_imu").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	Eigen::Matrix3d res;
	double R11, R12, R13, R21, R22, R23, R31, R32, R33;
	ss >> R11 >> R12 >> R13 >> R21 >> R22 >> R23 >> R31 >> R32 >> R33;
	res << R11, R12, R13, R21, R22, R23, R31, R32, R33;
	//_gmutex.unlock();
	return res;
}

Eigen::Vector3d glidar::t_gsetlidar::t_lidar_imu()
{
	//_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("t_lidar_imu").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	Eigen::Matrix3d res;
	double X, Y, Z;
	ss >> X >> Y >> Z;
	//_gmutex.unlock();
	return Eigen::Vector3d(X, Y, Z);
}

Eigen::Isometry3d glidar::t_gsetlidar::T_lidar_imu()
{
	_gmutex.lock();
	Eigen::Matrix3d tmp_R_lidar_imu = R_lidar_imu();
	Eigen::Vector3d tmp_t_lidar_imu = t_lidar_imu();
	Eigen::Isometry3d tmp_T_lidar_imu=Eigen::Isometry3d::Identity();
	tmp_T_lidar_imu.linear() = tmp_R_lidar_imu;
	tmp_T_lidar_imu.translation() = tmp_t_lidar_imu;
	_gmutex.unlock();
	return tmp_T_lidar_imu;
}


int glidar::t_gsetlidar::n_scans()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("n_scans").as_int();
	_gmutex.unlock();
	return res;
}

int glidar::t_gsetlidar::map_gap()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("map_gap").as_int();
	_gmutex.unlock();
	return res;
}

bool glidar::t_gsetlidar::use_scan()
{
	_gmutex.lock();
	bool res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("use_scan").as_bool();
	_gmutex.unlock();
	return res;
}

bool glidar::t_gsetlidar::use_map()
{
	_gmutex.lock();
	bool res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("use_map").as_bool();
	_gmutex.unlock();
	return res;
}

double glidar::t_gsetlidar::scan_period()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("scan_period").as_double();
	_gmutex.unlock();
	return res;
}

double glidar::t_gsetlidar::map_resolution_sharp()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("map_resolution_sharp").as_double();
	_gmutex.unlock();
	return res;
}

double glidar::t_gsetlidar::map_resolution_surf()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("map_resolution_surf").as_double();
	_gmutex.unlock();
	return res;
}

int glidar::t_gsetlidar::max_iteration()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("max_iteration").as_int();
	_gmutex.unlock();
	return res;
}

int glidar::t_gsetlidar::max_lidar_state_size()
{
	_gmutex.lock();
	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("max_lidar_state_size").as_int();
	_gmutex.unlock();
	return res;
}


double glidar::t_gsetlidar::screenout_threshold()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("screenout_threshold").as_double();
	_gmutex.unlock();
	return res;
}

double glidar::t_gsetlidar::lidar_observation_noise()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("lidar_observation_noise").as_double();
	_gmutex.unlock();
	return res;
}

double glidar::t_gsetlidar::surf_dis_threshold()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("surf_dis_threshold").as_double();
	_gmutex.unlock();
	return res;
}

bool glidar::t_gsetlidar::proj_lastmap()
{
	_gmutex.lock();
	bool res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("proj_lastmap").as_bool();
	_gmutex.unlock();
	return res;
}

bool glidar::t_gsetlidar::distortion()
{
	_gmutex.lock();
	bool res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("distortion").as_bool();
	_gmutex.unlock();
	return res;
}

double glidar::t_gsetlidar::close_threshold()
{
	_gmutex.lock();

	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("close_threshold").as_double();
	res = res;
	_gmutex.unlock();
	return res;
}

double glidar::t_gsetlidar::sharp_dis_threshold()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("sharp_dis_threshold").as_double();
	_gmutex.unlock();
	return res;
}

double glidar::t_gsetlidar::huber_parameter()
{
	_gmutex.lock();
	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("huber_parameter").as_double();
	_gmutex.unlock();
	return res;
}

string glidar::t_gsetlidar::output_path()
{
	_gmutex.lock();
	string res = _doc.child(XMLKEY_ROOT).child(XMLKEY_LIDAR).attribute("output_path").as_string();
	_gmutex.unlock();
	return res;
}





