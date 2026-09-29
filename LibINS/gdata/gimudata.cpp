#include "gimudata.h"
using namespace Eigen;

gsins::t_gimudata::t_gimudata():t_gdata()
{
	gtrace("t_gimudata::t_gimudata");
	first = end = 1;
	id_type(IMUDATA);
}

int gsins::t_gimudata::add_IMU(const double& t, const Vector3d& wm, const Vector3d& vm)
{
	gtrace("t_gimudata::add_IMU");

#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	dataIMU tmp = { t,wm,vm };
	_imu_forward.push_back(tmp);
	if (_FBS)_imu_back.push_back(tmp);

	_gmutex.unlock();
	return 0;
}

int gsins::t_gimudata::add_IMU(const double& t, const Vector3d& wm, const Vector3d& vm, const Vector3d& mm)
{
	gtrace("t_gimudata::add_IMU");

#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	dataIMU tmp = { t,wm,vm,mm };
	_imu_forward.push_back(tmp);
	if (_FBS)_imu_back.push_back(tmp);
	_gmutex.unlock();
	return 0;
}

void gsins::t_gimudata::set_ts(double ts)
{
	_ts = ts;
}

bool gsins::t_gimudata::load(vector<Eigen::Vector3d>& wm, vector<Eigen::Vector3d>& vm, double & t, double & ts, int nSamples,bool _beg_end)
{
	gtrace("t_gimudata::load");

#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();
	wm.clear(); vm.clear();

	if ((_beg_end&&_imu_forward.size() < nSamples) || (!_beg_end&&_imu_back.size() < nSamples)) {
		t = -1; _gmutex.unlock(); return false;
	}
	_pre_time = t;
	for (int j = 0; j < nSamples; j++)
	{
		if (_beg_end)
		{
			t = _imu_forward.front().t;
			wm.push_back(_imu_forward.front().wm);
			vm.push_back(_imu_forward.front().vm);
			_imu_forward.pop_front();
		}
		else
		{
			t = _imu_back.back().t;
			wm.push_back(-_imu_back.back().wm);
			vm.push_back(_imu_back.back().vm);
			_imu_back.pop_back();
		}

	}
	if (first == 1) {
		ts = nSamples * _ts;
		first = 0;
	}
	else ts = t - _pre_time;

	_gmutex.unlock();
	return true;
	}

bool gsins::t_gimudata::load(vector<Eigen::Vector3d>& wm, vector<Eigen::Vector3d>& vm, vector<Eigen::Vector3d>& mm, double & t, double & ts, int nSamples, bool _beg_end)
{
	gtrace("t_gimudata::load");

#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();
	wm.clear(); vm.clear(); mm.clear();

	if ((_beg_end&&_imu_forward.size() < nSamples)|| (!_beg_end&&_imu_back.size() < nSamples)){
		t = 0; _gmutex.unlock(); return false;
	}
	_pre_time = t;
	for (int j = 0; j < nSamples; j++)
	{
		if (_beg_end)
		{
			t = _imu_forward.front().t;
			wm.push_back(_imu_forward.front().wm);
			vm.push_back(_imu_forward.front().vm);
			mm.push_back(_imu_forward.front().mm);
			_imu_forward.pop_front();
		}
		else
		{
			t = _imu_back.back().t;
			wm.push_back(-_imu_back.back().wm);
			vm.push_back(_imu_back.back().vm);
			mm.push_back(_imu_back.back().mm);
			_imu_back.pop_back();
		}
	}
	if (first == 1) {
		ts = nSamples * _ts;
		first = 0;
	}
	else ts = t - _pre_time;

	_gmutex.unlock();
	return true;
}

t_gtime gsins::t_gimudata::erase_bef(t_gtime t, bool _beg_end)
{
	double t_ins;
	double t_gnss = t.sow() + t.dsec();
	if (_beg_end)
	{
		t_ins = _imu_forward.front().t;
		while (t_ins < t_gnss)
		{
			_imu_forward.pop_front();
			t_ins = _imu_forward.front().t;
		}
	}	
	else
	{
		t_ins = _imu_back.back().t;
		while (t_ins > t_gnss)
		{
			_imu_back.pop_back();
			t_ins = _imu_back.back().t;
		}
	}
	return t_gtime(t.gwk(), t_ins);
}

int gsins::t_gimudata::size(bool _beg_end)
{
	if(_beg_end) return _imu_forward.size();
	else return _imu_back.size();
}

bool gsins::t_gimudata::available(const t_gtime & now, bool _beg_end)
{
	try {
		double t = now.sow() + now.dsec();
		if ((_beg_end && t < _imu_forward.back().t) || (!_beg_end && t > _imu_back.front().t))
			return true;
		else
			return false;
	}

	catch(...)
	{ 
		return false;
		//throw("no enough imu data!");
	}
}

double gsins::t_gimudata::beg_obs(bool _beg_end)
{
	if (_beg_end)
	{
		if (_imu_forward.size() == 0)return 0;
		return _imu_forward.front().t;
	}
	else
	{
		if (_imu_back.size() == 0)return 0;
		return _imu_back.front().t;
	}
}

double gsins::t_gimudata::end_obs(bool _beg_end)
{
	if (_beg_end)
	{
		if (_imu_forward.size() == 0)return 0;
		return _imu_forward.back().t;
	}
	else
	{
		if (_imu_back.size() == 0)return 0;
		return _imu_back.back().t;
	}
}
