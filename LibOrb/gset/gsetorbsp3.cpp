/**
*
* @verbatim
	History
	 -1.0 Yujie Qin		2019-12-02 creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gsetorbsp3.cpp
* @brief		orb to sp3 parameters' set from XML
*
* @author       Yujie Qin	, Wuhan University
* @version		1.0.0
* @date		    2019-12-02
*
*/

#include "gset/gsetorbsp3.h"
#include <algorithm>
#include <sstream>

#define DEBUG

namespace gnut {

	t_gsetorbsp3::t_gsetorbsp3()
	{

	}

	t_gsetorbsp3::~t_gsetorbsp3()
	{
	}



	vector<string> t_gsetorbsp3::sat()
	{
		vector<string> csat;
		string sat;
		istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_ORBSP3).child_value(XMLKEY_ORBSP3_SAT));
		while (is >> sat) {
			transform(sat.begin(), sat.end(), sat.begin(), ::toupper);
			csat.push_back(sat);
		}
		return csat;
	}


	double t_gsetorbsp3::intv()
	{
		double intv = str2dbl(_doc.child(XMLKEY_ROOT).child(XMLKEY_ORBSP3).child_value(XMLKEY_ORBSP3_INTV));
		if(intv == 0.0)
		{
			intv = 300;
		}
		return intv;
	}

	string t_gsetorbsp3::frame()
	{
		string frame;
		string sframe;
		istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_ORBSP3).child_value(XMLKEY_ORBSP3_FRAME));
		while (is >> sframe) {
			transform(sframe.begin(), sframe.end(), sframe.begin(), ::tolower);
			frame=sframe;
		}
		if (frame.size() == 0)
		{
			return "trs";
		}
		if (frame != "crs"&&frame != "trs")
		{
			frame = "trs";

		}
		return frame;
	}


	void t_gsetorbsp3::check()
	{
	}

	void t_gsetorbsp3::help()
	{
		cerr << "<orbsp3>" << endl
			 << "<sat>      G01   </sat> " << endl
			 << "<intv>     30    </intv> " << endl
			<< "<frame>    trs    </frame> " << endl
			 << "</orbsp3>" << endl;

	}
}