#include "gsetign.h"
using namespace Eigen;
using namespace gsins;

gsins::t_gsetign::t_gsetign() : t_gsetbase()
{
	_set.insert(XMLKEY_IGN);

	_initial_misalignment_std = Vector3d(1, 1, 10);
	_initial_vel_std = 0.1*Vector3d::Ones();
	_initial_pos_std = 3 * Vector3d::Ones();
	_initial_gyro_std = 10 * Vector3d::Ones();
	_initial_acce_std = 5 * Vector3d::Ones();

	_min_misalignment_std = 0.1*Vector3d::Ones();
	_min_vel_std = 0.005*Vector3d::Ones();
	_min_pos_std = 0.01*Vector3d::Ones();
	_min_gyro_std = 0.1*Vector3d::Ones();
	_min_acce_std = 0.1*Vector3d::Ones();

	_misalignment_psd = 1 * Vector3d::Ones();
	_vel_psd = 0.1*Vector3d::Ones();
	_pos_psd = 0 * Vector3d::Ones();
	_gyro_psd = 0 * Vector3d::Ones();
	_acce_psd = 0 * Vector3d::Ones();

	_meas_vel_noise = Vector3d(0.5, 0.5, 0.5);
	_meas_pos_noise = Vector3d(0.5, 0.5, 0.5);

	_nq = 15; _nr = 6;
	_lever = Vector3d::Zero();
	_max_pdop = 8.0; 
	_min_sat = 5;
	_order = "XYZXYZ";
	_TS = 1.0;
	_delay_t = 0.0;
	_gnssmode = LCI_TYPE::PPP;
	_fltmode = FORWARD;

}

gsins::t_gsetign::~t_gsetign()
{
}

void gsins::t_gsetign::check()
{
	_gmutex.lock();
	xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN);

	string str = to_string(_initial_misalignment_std(0)) + ' ' + to_string(_initial_misalignment_std(1)) + ' ' + to_string(_initial_misalignment_std(2));
	_default_attr(node, "initial_misalignment_std", str);

	str.clear(); str = "";
	str = to_string(_initial_vel_std(0)) + ' ' + to_string(_initial_vel_std(1)) + ' ' + to_string(_initial_vel_std(2));
	_default_attr(node, "initial_vel_std", str);

	str.clear(); str = "";
	str = to_string(_initial_pos_std(0)) + ' ' + to_string(_initial_pos_std(1)) + ' ' + to_string(_initial_pos_std(2));
	_default_attr(node, "initial_pos_std", str);

	_default_attr(node, "initial_gyro_std", _initial_gyro_std(0));
	_default_attr(node, "initial_acce_std", _initial_acce_std(0));

	_default_attr(node, "min_misalignment_std", _min_misalignment_std(0));
	_default_attr(node, "min_vel_std", _min_vel_std(0));
	_default_attr(node, "min_pos_std", _min_pos_std(0));
	_default_attr(node, "min_gyro_std", _min_gyro_std(0));
	_default_attr(node, "min_acce_std", _min_acce_std(0));

	_default_attr(node, "misalignment_psd", _misalignment_psd(0));
	_default_attr(node, "vel_psd", _vel_psd(0));
	_default_attr(node, "pos_psd", _pos_psd(0));
	_default_attr(node, "gyro_psd", _gyro_psd(0));
	_default_attr(node, "acce_psd", _acce_psd(0));

	_default_attr(node, "hor_vel_noise", _meas_vel_noise(0));
	_default_attr(node, "ver_vel_noise", _meas_vel_noise(2));
	_default_attr(node, "hor_pos_noise", _meas_pos_noise(0));
	_default_attr(node, "ver_pos_noise", _meas_pos_noise(2));

	str.clear(); str = "";
	str = to_string(_lever(0)) + ' ' + to_string(_lever(1)) + ' ' + to_string(_lever(2));
	_default_attr(node, "lever", str);

	_default_attr(node, "nq", _nq);
	_default_attr(node, "nr", _nr);
	_default_attr(node, "order", _order);
	_default_attr(node, "TS", _TS);
	_default_attr(node, "delay_t", _delay_t);
	_default_attr(node, "filtermode", _fltmode);
	_default_attr(node, "gnssmode", _gnssmode);
	_default_attr(node, "max_pdop", _max_pdop);
	_default_attr(node, "min_satnum", _min_sat);

	_gmutex.unlock();
}

void gsins::t_gsetign::help()
{
	//todo 不需要加锁
	_gmutex.lock();

	// no Implement

	_gmutex.unlock();
}

//todo 取值不需要加锁
Eigen::Vector3d gsins::t_gsetign::initial_misalignment_std()
{

	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("initial_misalignment_std").value();
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

Eigen::Vector3d gsins::t_gsetign::initial_vel_std()
{

	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("initial_vel_std").value();
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

Eigen::Vector3d gsins::t_gsetign::initial_pos_std()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("initial_pos_std").value();
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

Eigen::Vector3d gsins::t_gsetign::initial_gyro_std()
{
	_gmutex.lock();
	double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("initial_gyro_std").as_double();
	_gmutex.unlock();
	return Eigen::Vector3d::Ones()*x;
}

Eigen::Vector3d gsins::t_gsetign::initial_acce_std()
{
	_gmutex.lock();
	double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("initial_acce_std").as_double();
	_gmutex.unlock();
	return Eigen::Vector3d::Ones()*x;
}

Eigen::Vector3d gsins::t_gsetign::min_misalignment_std()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("min_misalignment_std").value();
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

Eigen::Vector3d gsins::t_gsetign::min_vel_std()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("min_vel_std").value();
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

Eigen::Vector3d gsins::t_gsetign::min_pos_std()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("min_pos_std").value();
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

Eigen::Vector3d gsins::t_gsetign::min_gyro_std()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("min_gyro_std").value();
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

Eigen::Vector3d gsins::t_gsetign::min_acce_std()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("min_acce_std").value();
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

Eigen::Vector3d gsins::t_gsetign::misalignment_psd()
{
	_gmutex.lock();
	double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("misalignment_psd").as_double();
	_gmutex.unlock();
	return Eigen::Vector3d::Ones()*x;
}

Eigen::Vector3d gsins::t_gsetign::vel_psd()
{
	_gmutex.lock();
	double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("vel_psd").as_double();
	_gmutex.unlock();
	return Eigen::Vector3d::Ones()*x;
}

Eigen::Vector3d gsins::t_gsetign::pos_psd()
{
	_gmutex.lock();
	double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("pos_psd").as_double();
	_gmutex.unlock();
	return Eigen::Vector3d::Ones()*x;
}

Eigen::Vector3d gsins::t_gsetign::gyro_psd()
{
	_gmutex.lock();
	double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("gyro_psd").as_double();
	_gmutex.unlock();
	return Eigen::Vector3d::Ones()*x;
}

Eigen::Vector3d gsins::t_gsetign::acce_psd()
{
	_gmutex.lock();
	double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("acce_psd").as_double();
	_gmutex.unlock();
	return Eigen::Vector3d::Ones()*x;
}

Eigen::Vector3d gsins::t_gsetign::meas_vel_noise()
{
	_gmutex.lock();
	double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("hor_vel_noise").as_double();
	double y = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("ver_vel_noise").as_double();
	_gmutex.unlock();
	return Eigen::Vector3d(x, x, y);
}

Eigen::Vector3d gsins::t_gsetign::meas_pos_noise()
{
	_gmutex.lock();
	double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("hor_pos_noise").as_double();
	double y = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("ver_pos_noise").as_double();
	_gmutex.unlock();
	return Eigen::Vector3d(x, x, y);
}

string gsins::t_gsetign::order()
{
	_gmutex.lock();
	string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("order").value();
	_gmutex.unlock();
	transform(x.begin(), x.end(), x.begin(), ::toupper);
	return x;
}

int gsins::t_gsetign::nq()
{
	_gmutex.lock();
	double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("nq").as_double();
	_gmutex.unlock();
	return x;
}

int gsins::t_gsetign::nr()
{
	_gmutex.lock();
	double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("nr").as_double();
	_gmutex.unlock();
	return x;
}

double gsins::t_gsetign::TS()
{
	_gmutex.lock();
	double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("TS").as_double();
	_gmutex.unlock();
	return x;
}

double gsins::t_gsetign::delay_t()
{
	_gmutex.lock();
	double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("delay_t").as_double();
	_gmutex.unlock();
	return x;
}

double gsins::t_gsetign::max_pdop()
{
	_gmutex.lock();
	double x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("max_pdop").as_double();
	_gmutex.unlock();
	return x;
}

int gsins::t_gsetign::min_sat()
{
	_gmutex.lock();
	int x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("min_satnum").as_int();
	_gmutex.unlock();
	return x;
}

Eigen::Vector3d gsins::t_gsetign::lever()
{
	_gmutex.lock();
	string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("lever").value();
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')str[i] = ' ';
	}
	stringstream ss(str);
	double lx, ly, lz;
	ss >> lx >> ly >> lz;
	_gmutex.unlock();
	return Vector3d(lx, ly, lz);
}

IGN_TYPE gsins::t_gsetign::ign_type()
{
	_gmutex.lock();
	string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("ign_type").value();
	_gmutex.unlock();
	return gsins::str2ign(x);
}

gsins::FLT_TYPE gsins::t_gsetign::fltmode()
{
	_gmutex.lock();
	string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("filtermode").value();
	_gmutex.unlock();
	return gsins::str2Flt(x);
}

gsins::LCI_TYPE gsins::t_gsetign::gnssmode()
{
	_gmutex.lock();
	string x = _doc.child(XMLKEY_ROOT).child(XMLKEY_IGN).attribute("gnssmode").value();
	_gmutex.unlock();
	return gsins::str2lci(x);
}
