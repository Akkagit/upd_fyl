/**
*
* The format :
* @verbatim
	History
	 -1.0 jdhuang	 2019-04-07 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gsetpar.cpp
* @brief		set pars from XML
*
* @author       jdhuang, Wuhan University
* @version		1.0.0
* @date			2019-04-07
*
*/

#include "gset/gsetxml.h"
#include "gutils/gstring.h"
//#define DEBUG
using namespace gnut;
namespace great
{
	t_gsetxml::t_gsetxml() : t_gsetbase()
	{
		_set.insert(XMLKEY_XML);
	}

	t_gsetxml::~t_gsetxml()
	{

	}

	void t_gsetxml::help()
	{
	}


	string t_gsetxml::xml_eg()
	{
		string tmp = gnut::trim(_doc.child(XMLKEY_ROOT).child(XMLKEY_XML).child_value("xml_eg"));
		if (tmp.empty()) return "pod.xml";
		else return tmp;
	}

	string t_gsetxml::app_type()
	{
		string tmp = gnut::trim(_doc.child(XMLKEY_ROOT).child(XMLKEY_XML).child_value("app_type"));
		if (tmp.empty()) return "PPP_IF";
		else return tmp;
	}

	string t_gsetxml::sp3_type()
	{
		string tmp = gnut::trim(_doc.child(XMLKEY_ROOT).child(XMLKEY_XML).child_value("sp3_type"));
		if (tmp.empty()) return "igs";
		else return tmp;
	}

	string t_gsetxml::rinexc_type()
	{
		string tmp = gnut::trim(_doc.child(XMLKEY_ROOT).child(XMLKEY_XML).child_value("rinexc_type"));
		if (tmp.empty()) return "igs";
		else return tmp;
	}

	string t_gsetxml::sp3_root()
	{
		string tmp = gnut::trim(_doc.child(XMLKEY_ROOT).child(XMLKEY_XML).child_value("sp3_root"));
		if (tmp.empty()) return project_root();
		else return tmp;
	}

	string t_gsetxml::brd_root()
	{
		string tmp = gnut::trim(_doc.child(XMLKEY_ROOT).child(XMLKEY_XML).child_value("brd_root"));
		if (tmp.empty()) return project_root();
		else return tmp;
	}

	string t_gsetxml::dcb_root()
	{
		string tmp = gnut::trim(_doc.child(XMLKEY_ROOT).child(XMLKEY_XML).child_value("dcb_root"));
		if (tmp.empty()) return project_root();
		else return tmp;
	}

	string t_gsetxml::sinex_root()
	{
		string tmp = gnut::trim(_doc.child(XMLKEY_ROOT).child(XMLKEY_XML).child_value("sinex_root"));
		if (tmp.empty()) return project_root();
		else return tmp;		
	}

	string t_gsetxml::rinexc_root()
	{
		string tmp = gnut::trim(_doc.child(XMLKEY_ROOT).child(XMLKEY_XML).child_value("rinexc_root"));
		if (tmp.empty()) return project_root();
		else return tmp;
	}

	string t_gsetxml::rinexo_root()
	{
		string tmp = gnut::trim(_doc.child(XMLKEY_ROOT).child(XMLKEY_XML).child_value("rinexo_root"));
		if (tmp.empty()) return project_root();
		else return tmp;
	}

	string t_gsetxml::sysfile_root()
	{
		string tmp = gnut::trim(_doc.child(XMLKEY_ROOT).child(XMLKEY_XML).child_value("sysfile_root"));
		if (tmp.empty()) return project_root();
		else return tmp;
	}

	string t_gsetxml::project_root()
	{
		string tmp = gnut::trim(_doc.child(XMLKEY_ROOT).child(XMLKEY_XML).child_value("project_root"));
		if (tmp.empty()) return "./";
		else return tmp;
	}

	double t_gsetxml::time_interval()
	{
		string tmp = gnut::trim(_doc.child(XMLKEY_ROOT).child(XMLKEY_XML).child_value("time_interval"));
		return str2dbl(tmp);
	}

	bool t_gsetxml::scan_sysfile()
	{
		string tmp = gnut::trim(_doc.child(XMLKEY_ROOT).child(XMLKEY_XML).child_value("scan_sysfile"));
		transform(tmp.begin(), tmp.end(), tmp.begin(), (int(*)(int))toupper);
		if (tmp == "TRUE") return true;
		else return false;
	}

	bool t_gsetxml::scan_rinexo()
	{
		string tmp = gnut::trim(_doc.child(XMLKEY_ROOT).child(XMLKEY_XML).child_value("scan_rinexo"));
		transform(tmp.begin(), tmp.end(), tmp.begin(), (int(*)(int))toupper);
		if (tmp == "TRUE") return true;
		else return false;
	}

}



