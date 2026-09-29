
#include "glidarproc.h"
#include "gset/gsetlidar.h"
using namespace Eigen;
using namespace glidar;

glidar::t_glidarproc::t_glidarproc(t_gsetbase* _set)
{
	///< basic patameters
	ts = dynamic_cast<t_gsetlidar*>(_set)->ts();
	freq = dynamic_cast<t_gsetlidar*>(_set)->freq();
	imu_ts = dynamic_cast<t_gsetlidar*>(_set)->imu_ts();

	///< initial the parameters
	R_lidar_imu = dynamic_cast<t_gsetlidar*>(_set)->R_lidar_imu();
	t_lidar_imu = dynamic_cast<t_gsetlidar*>(_set)->t_lidar_imu();
	T_lidar_imu = dynamic_cast<t_gsetlidar*>(_set)->T_lidar_imu();

	use_scan = dynamic_cast<t_gsetlidar*>(_set)->use_scan();
	use_map=dynamic_cast<t_gsetlidar*>(_set)->use_map();
	map_gap=dynamic_cast<t_gsetlidar*>(_set)->map_gap();
	proj_lastmap = dynamic_cast<t_gsetlidar*>(_set)->proj_lastmap();
}

void glidar::t_glidarproc::load_imuobs(const double &t,const vector<Eigen::Vector3d> & wm,const vector<Eigen::Vector3d> & vm)
{
	
	assert(wm.size() == vm.size());
	/* only one sample can be processed successfully*/
	assert(wm.size() == 1);
	for (int i = 0; i < wm.size(); i++)
	{
		imumsg tmp_imu;
		tmp_imu.t = t;
		tmp_imu.angular_velocity = wm[i]/imu_ts;
		tmp_imu.linear_acceleration = vm[i]/imu_ts;
		_vecimu.push_back(tmp_imu);
	}
	
}

void glidar::t_glidarproc::load_lidarobs(const double &t, const lidarPath &lidar_path)
{
	cur_lidar_path = lidar_path;
}

LidarFrame glidar::t_glidarproc::ProcessBatch(const double& t)
{
	LidarFrame res;
	return res;
}
