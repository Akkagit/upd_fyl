
#include "gimgdata.h"

gvis::t_gimgdata::t_gimgdata() :t_gdata()
{
	gtrace("t_gimgdata::t_gimgdata");
	id_type(CAMDATA);
}

int gvis::t_gimgdata::add_IMG(const double& t, const string& img0_path, const string& img1_path)
{
	gtrace("t_gimgdata::add_CAM");

#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	imgPath tmp = { t,img0_path,img1_path };
	_vecimg.push_back(tmp);

	_gmutex.unlock();
	return 0;
}


void gvis::t_gimgdata::load(double& imu_t, double& img_t, imgPath& img_path)
{
	//cout << "imut_t" << setprecision(12) << imu_t << endl;
	for (int i = 0; i < _vecimg.size(); i++)
	{
		double dt = _vecimg.at(i).t - imu_t;
		if (dt > 1e-3)
		{
			if (i != 0)
			{
				img_t = _vecimg.at(i - 1).t;
				img_path = _vecimg.at(i - 1);
			}
			else
			{
				img_t = 0.0;
			}
			break;
		}
		else
			img_t = 0.0;
	}
	//cout << "img_t:" << setprecision(12) << img_t << endl;
}

int gvis::t_gimgdata::size()
{
	return _vecimg.size();
}


