
#include "imgfile.h"
#include "gdata/gimgdata.h"
#include "gset/gsetvis.h"
using namespace Eigen;

gvis::t_imgfile::t_imgfile(t_gsetbase * s, string version, int sz) :t_gcoder(s, version, sz),
_ts(0.1)
{
	gtrace("t_imgfile::constructor");
	_ts = dynamic_cast<t_gsetvis*>(s)->ts();
}

int gvis::t_imgfile::decode_head(char * buff, int sz, vector<string>& errmsg)
{
#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_mutex.lock();

	// no header expected, but fill the buffer
	t_gcoder::_add2buffer(buff, sz);
	_mutex.unlock(); return -1;
}

int gvis::t_imgfile::decode_data(char * buff, int sz, int & cnt, vector<string>& errmsg)
{
	gtrace("t_imgfile::decode_data");
	_mutex.lock();

	if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

	int tmpsize = 0;
	string img0_path = "cam0\\";
	string img1_path = "cam1\\";
	string line;
	while ((tmpsize = t_gcoder::_getline(line, 0))>= 0)
	{
		for (int i = 0; i < line.size(); i++)
		{
			if (line[i] == ',' || line[i] == '*' || line[i] == ';')line[i] = ' ';
		}

		stringstream ss(line);
		string path0, path1, path;
		double time;
		ss >> time;
		if(time>1e11)
		time = time / 1e9;
		ss >> path;
		int idx = path.find("jpg");
		if (idx == string::npos) path = path + ".jpg";
	
		
		map<string, t_gdata*>::iterator it = _data.begin();
		while (it != _data.end())
		{
			if (it->second->id_type() == t_gdata::CAMDATA)
			((t_gimgdata*)it->second)->add_IMG(time, img0_path+path, img1_path+path);
			it++;
		}
		if (ss.fail()) {
			if (_log) _log->comment(1, "imgfile", "warning: incorrect IMG data record: " + ss.str());
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

