/**
*
* @verbatim
	History
	 -1.0 Yujie Qin		2019-12-02 creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gsetsp3orb.cpp
* @brief		sp3 to orb parameters' set from XML
*
* @author       Yujie Qin	, Wuhan University
* @version		1.0.0
* @date		    2019-12-02
*
*/

#include "gset/gsetsp3orb.h"
#include <algorithm>
#include <sstream>

//#define DEBUG

namespace gnut {

	t_gsetsp3orb::t_gsetsp3orb()
	{

	}

	t_gsetsp3orb::~t_gsetsp3orb()
	{
	}

	

	set<string> t_gsetsp3orb::sat()
	{
		set<string> csat;
		string sat;
		istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_SP3ORB).child_value(XMLKEY_SP3ORB_SAT));
		while (is >> sat) {
			transform(sat.begin(), sat.end(), sat.begin(), ::toupper);
			csat.insert(sat);
		}
		return csat;
	}


	/*int t_gsetsp3orb::drag()
	{
		xml_node sp3orb_model = _doc.child(XMLKEY_ROOT).child(XMLKEY_SP3ORB);
		int drag_min = 0;
		for (xml_node tmp_set = sp3orb_model.first_child(); tmp_set; tmp_set = tmp_set.next_sibling())
		{
			if ((strcmp("drag", tmp_set.name())) == 0)
			{
				drag_min = str2int(tmp_set.attribute("time").value());
				break;
			}				
		}
		if (drag_min == 0)
		{
			drag_min = 360;
			cout << "ERROR: t_gsetsp3orb::unknown drag min" << endl;
		}
		
		return drag_min;
	}

	map<string, int> t_gsetsp3orb::drag_parameter()
	{
		xml_node sp3orb_model = _doc.child(XMLKEY_ROOT).child(XMLKEY_SP3ORB);
		map<string, int> drag_para;
		set<string> drag;
		for (xml_node tmp_set = sp3orb_model.first_child(); tmp_set; tmp_set = tmp_set.next_sibling())
		{
			if ((strcmp("drag", tmp_set.name())) == 0)
			{
				drag_para.clear();
				string para = tmp_set.first_child().value();
				string word;
				istringstream is(para);
				while (is >> word) {
					drag.insert(word);
				}
				break;
			}
		}
		set<string>::iterator ipara = drag.begin();
		for (ipara; ipara != drag.end(); ipara++)
		{
			size_t m = (*ipara).find_first_of(":");
			size_t n = (*ipara).size();
			string name = (*ipara).substr(0, m);
			int min = str2int((*ipara).substr(m+1, n));
			drag_para[name] = min;
		}
		if (drag_para.size() == 0)
		{
			drag_para["DRAG_c"] = 360;
			cout << "ERROR: t_gsetsp3orb::unknown drag para" << endl;
		}
		return drag_para;
	}

	
	map<string, int> t_gsetsp3orb::cust_parameter()
	{
		xml_node sp3orb_model = _doc.child(XMLKEY_ROOT).child(XMLKEY_SP3ORB);
		map<string, int> cust_para;
		set<string> cust;
		for (xml_node tmp_set = sp3orb_model.first_child(); tmp_set; tmp_set = tmp_set.next_sibling())
		{
			if ((strcmp("customer", tmp_set.name())) == 0)
			{
				cust_para.clear();
				string para = tmp_set.first_child().value();
				string word;
				istringstream is(para);
				while (is >> word) {
					cust.insert(word);
				}
				break;
			}
		}
		set<string>::iterator ipara = cust.begin();
		for (ipara; ipara != cust.end(); ipara++)
		{
			size_t m = (*ipara).find_first_of(":");
			size_t n = (*ipara).size();
			string name = (*ipara).substr(0, m);
			int min = str2int((*ipara).substr(m+1, n));
			cust_para[name] = min;
		}
		if (cust_para.size() == 0)
		{
			cust_para["EMP_Sc1"] = 90;
			cust_para["EMP_Cc1"] = 90;
			cust_para["EMP_Sa1"] = 90;
			cust_para["EMP_Ca1"] = 90;
			cust_para["EMP_Sr1"] = 90;
			cust_para["EMP_Cr1"] = 90;
			cout << "ERROR: t_gsetsp3orb::unknown customer para" << endl;
		}
		return cust_para;
	}*/

	string t_gsetsp3orb::frame()
	{
		set<string> frame;
		string sframe;
		istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_SP3ORB).child_value(XMLKEY_SP3ORB_FRAME));
		while (is >> sframe) {
			transform(sframe.begin(), sframe.end(), sframe.begin(), ::tolower);
			frame.insert(sframe);
		}
		if (frame.size() == 0)
		{
			return "crs";
		}
		return *frame.begin();
	}

	void t_gsetsp3orb::check()
	{
	}

	void t_gsetsp3orb::help()
	{
		cerr << "<sp3orb>" << endl
			<< "  <sat>    G01    </sat> " << endl
			<< "<frame>    crs    </frame> " << endl
			<< "</sp3orb>" << endl;

	}
}