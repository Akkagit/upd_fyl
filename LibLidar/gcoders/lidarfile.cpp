#include "lidarfile.h"
#include "gdata/glidardata.h"
#include "gset/gsetlidar.h"

using namespace Eigen;

glidar::t_lidarfile::t_lidarfile(t_gsetbase * s, string version, int sz) :t_gcoder(s, version, sz),
_ts(0.1)
{
	gtrace("t_lidarfile::constructor");
	_ts = dynamic_cast<t_gsetlidar*>(s)->ts();
}

int glidar::t_lidarfile::decode_head(char * buff, int sz, vector<string>& errmsg)
{
#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_mutex.lock();

	// no header expected, but fill the buffer
	t_gcoder::_add2buffer(buff, sz);
	_mutex.unlock(); return -1;
}

//这一步做了什么?
int glidar::t_lidarfile::decode_data(char * buff, int sz, int & cnt, vector<string>& errmsg)
{
	gtrace("t_lidarfile::decode_data");
	_mutex.lock();

	if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

	int tmpsize = 0;
	string lidar_path = "data\\";
	string line;
	while ((tmpsize = t_gcoder::_getline(line, 0))>= 0)
	{
		for (int i = 0; i < line.size(); i++)
		{
			if (line[i] == ',' || line[i] == '*' || line[i] == ';')line[i] = ' ';
		}

		stringstream ss(line);
		string path;
		double time;
		ss >> time;
		if(time>1e11)
		time = time / 1e9;
		ss >> path;
		int idx = path.find("bin");
		if (idx == string::npos) path = path + ".bin";
	
		
		map<string, t_gdata*>::iterator it = _data.begin();
		while (it != _data.end())
		{
			if (it->second->id_type() == t_gdata::LIDARDATA)
			((t_glidardata*)it->second)->add_LIDAR(time, lidar_path+path);
			it++;
		}
		if (ss.fail()) {
			if (_log) _log->comment(1, "lidarfile", "warning: incorrect LIDAR data record: " + ss.str());
			else       cerr << "warning: incorrect IMU data record: " << ss.str() << endl;
			t_gcoder::_consume(tmpsize);
			_mutex.unlock(); return -1;
		}
		t_gcoder::_consume(tmpsize);
		cnt++;

	}
	_mutex.unlock();
	return 0;
}

