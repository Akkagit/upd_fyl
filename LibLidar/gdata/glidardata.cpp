#include "glidardata.h"

glidar::t_glidardata::t_glidardata() :t_gdata()
{
	gtrace("t_glidardata::t_glidardata");
	id_type(LIDARDATA);
}

int glidar::t_glidardata::add_LIDAR(const double& t, const string& lidar_path)
{
	gtrace("t_glidardata::add_LIDAR");

#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	lidarPath tmp = { t,lidar_path};
	_veclidar.push_back(tmp);

	_gmutex.unlock();
	return 0;
}


void glidar::t_glidardata::load(double& imu_t, double& lidar_t, lidarPath& lidar_path)
{
	//cout << "imut_t" << setprecision(12) << imu_t << endl;
	for (int i = 0; i < _veclidar.size(); i++)
	{
		double dt = _veclidar.at(i).t - imu_t;
		if (dt > 0.005)
		{
			if (i != 0)
			{
				lidar_t = _veclidar.at(i - 1).t;
				lidar_path = _veclidar.at(i - 1);
			}
			else
			{
				lidar_t = 0.0;
			}
			break;
		}
		else
			lidar_t = 0.0;
	}
	//cout << "lidar_t:" << setprecision(12) << lidar_t << endl;
}

int glidar::t_glidardata::size()
{
	return _veclidar.size();
}


