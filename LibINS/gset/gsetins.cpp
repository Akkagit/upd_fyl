#include "gsetins.h"
using namespace gsins;

gsins::t_gsetins::t_gsetins()
{
	_set.insert(XMLKEY_INS);
	_freq = 200;
	_ts = 1.0 / _freq;
	_pos = _vel = _att = Eigen::Vector3d::Zero();
	_align_time = 5 * 60;
	_GyroUnit = "DPS";
	_AcceUnit = "MPS2";
	_cps = "POLY";
	_subsample = 1;
	_order = "garfu";
	_out_order = "XYZ_XYZ_PRY_NW";
	_out_freq = 1;
	_int_sec = true;
	_align_type = "SA";
	_start = 0.0;
	_end = _start + 7 * 24 * 60 * 60;
	_acce_bias = Eigen::Vector3d::Zero();
	_gyro_bias = Eigen::Vector3d::Zero();
}

gsins::t_gsetins::~t_gsetins()
{
}

void gsins::t_gsetins::check()
{
	_gmutex.lock();
	xml_node parent = _doc.child(XMLKEY_ROOT);
	xml_node node = _default_node(parent, XMLKEY_INS);


	string str = to_string(_pos(0)) + ' ' + to_string(_pos(1)) + ' ' + to_string(_pos(2));
	_default_attr(node, "pos", str);

	str.clear(); str = "";
	str = to_string(_vel(0)) + ' ' + to_string(_vel(1)) + ' ' + to_string(_vel(2));
	_default_attr(node, "vel", str);

	str.clear(); str = "";
	str = to_string(_att(0)) + ' ' + to_string(_att(1)) + ' ' + to_string(_att(2));
	_default_attr(node, "att", str);

	_default_attr(node, "freq", 1 / _ts);
	_default_attr(node, "gyro_unit", _GyroUnit);
	_default_attr(node, "acce_unit", _AcceUnit);
	_default_attr(node, "align_type", _align_type);
	_default_attr(node, "align_time", _align_time);
	_default_attr(node, "compensation", _cps);
	_default_attr(node, "subsample", _subsample);
	_default_attr(node, "order", _order);
	_default_attr(node, "out_freq", _out_freq);
	_default_attr(node, "int_sec", _int_sec);
	_default_attr(node, "out_order", _out_order);
	_default_attr(node, "start", _start);
	_default_attr(node, "end", _end);

	_default_node(node, "motion", "");
	
	str.clear(); str = "";
	str = to_string(_gyro_bias(0)) + ' ' + to_string(_gyro_bias(1)) + ' ' + to_string(_gyro_bias(2));
	_default_attr(node, "gyro_bias", str);

	str.clear(); str = "";
	str = to_string(_acce_bias(0)) + ' ' + to_string(_acce_bias(1)) + ' ' + to_string(_acce_bias(2));
	_default_attr(node, "acce_bias", str);

	_gmutex.unlock();
}

void gsins::t_gsetins::help()
{
	_gmutex.lock();
	cerr << " <ins \n"
		<< "   ts=\"" << _ts << "\" \n"
		<< "   gyro_unit=\"" << _GyroUnit << "\" \n"
		<< "   acce_unit=\"" << _AcceUnit << "\" \n"

		<< "  />\n";

	_gmutex.unlock();
}

double gsins::t_gsetins::ts()
{
	_gmutex.lock();

	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("freq").as_double();
	_freq = res;
	res = 1 / res;

	_gmutex.unlock();
	return res;
}

int gsins::t_gsetins::freq()
{
	_gmutex.lock();

	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("freq").as_int();

	_gmutex.unlock();
	return res;
}


gsins::UNIT_TYPE gsins::t_gsetins::GyroUnit()
{
	_gmutex.lock();

	string gyrounit = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("gyro_unit").value();

	_gmutex.unlock();

	return str2Unit(gyrounit);
}

gsins::UNIT_TYPE gsins::t_gsetins::AcceUnit()
{
	_gmutex.lock();

	string acceunit = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("acce_unit").value();

	_gmutex.unlock();

	return str2Unit(acceunit);
}

gsins::UNIT_TYPE gsins::t_gsetins::MagUnit()
{
	_gmutex.lock();

	string acceunit = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("mag_unit").value();

	_gmutex.unlock();

	return str2Unit(acceunit);
}

string gsins::t_gsetins::order()
{
	_gmutex.lock();

	string res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("order").value();
	transform(res.begin(), res.end(), res.begin(), ::tolower);

	_gmutex.unlock();

	return res;
}

ODR_TYPE gsins::t_gsetins::out_order()
{
	_gmutex.lock();

	string res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("out_order").value();

	_gmutex.unlock();

	return str2odr(res);
}

Eigen::Vector3d gsins::t_gsetins::pos()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("pos").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	double X, Y, Z;
	ss >> X >> Y >> Z;
	_gmutex.unlock();
	return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d gsins::t_gsetins::vel()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("vel").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	double X, Y, Z;
	ss >> X >> Y >> Z;
	_gmutex.unlock();
	return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d gsins::t_gsetins::att()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("att").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	double X, Y, Z;
	ss >> X >> Y >> Z;
	_gmutex.unlock();
	return Eigen::Vector3d(X, Y, Z)*glv.deg;
}

bool gsins::t_gsetins::align()
{
	_gmutex.lock();
	bool x = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("align").as_bool();
	_gmutex.unlock();
	return x;
}


int gsins::t_gsetins::cps()
{
	_gmutex.lock();
	string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("compensation").as_string();
	_gmutex.unlock();
	return str2cps(x);
}

int gsins::t_gsetins::subsample()
{
	_gmutex.lock();
	int x = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("subsample").as_int();
	_gmutex.unlock();
	return x;
}

double gsins::t_gsetins::start()
{
	_gmutex.lock();

	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("start").as_double();

	_gmutex.unlock();
	return res;
}

double gsins::t_gsetins::end()
{
	_gmutex.lock();

	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("end").as_double();

	_gmutex.unlock();
	return res;
}

double gsins::t_gsetins::align_time()
{
	_gmutex.lock();

	double res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("align_time").as_double();

	_gmutex.unlock();
	return res;
}

ALIGN_TYPE gsins::t_gsetins::align_type()
{
	_gmutex.lock();

	string res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("align_type").value();

	_gmutex.unlock();

	return str2align(res);
}

bool gsins::t_gsetins::out_intsec()
{
	_gmutex.lock();
	bool x = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("int_sec").as_bool();
	_gmutex.unlock();
	return x;
}

int gsins::t_gsetins::out_freq()
{
	_gmutex.lock();

	int res = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("out_freq").as_double();

	_gmutex.unlock();
	return res;
}

vector<string> gsins::t_gsetins::motion()
{
	_gmutex.lock();

	vector<string> motion_strvec = t_gsetbase::_vecval(XMLKEY_INS, "motion");

	_gmutex.unlock();
	return motion_strvec;

}

Eigen::Vector3d gsins::t_gsetins::acce_bias()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("acce_bias").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	double X, Y, Z;
	ss >> X >> Y >> Z;
	_gmutex.unlock();
	return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d gsins::t_gsetins::gyro_bias()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("gyro_bias").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	double X, Y, Z;
	ss >> X >> Y >> Z;
	_gmutex.unlock();
	return Eigen::Vector3d(X, Y, Z);
}
