/**
*
* @verbatim
	History
	 -1.0 Yujie Qin		2019-12-02 creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gsetsp3orb.h
* @brief		sp3 to orb parameters' set from XML
*
* @author       Yujie Qin	, Wuhan University
* @version		1.0.0
* @date		    2019-12-02
*
*/

#ifndef GSETSP3ORB_H
#define GSETSP3ORB_H

#define XMLKEY_SP3ORB 	           "sp3orb"
#define XMLKEY_SP3ORB_SAT          "sat"
#define XMLKEY_SP3ORB_FRAME        "frame"


#include "gset/gsetbase.h"
#include "gutils/ExportOI.h"

namespace gnut
{
	class OI_LIBRARY_EXPORT t_gsetsp3orb : public virtual t_gsetbase
	{
	public:
		/** @brief constructor. */
		t_gsetsp3orb();

		/** @brief destructor. */
		virtual ~t_gsetsp3orb();

 	    /**
		* @brief get the satellites that need calculation.
		* @return the sat
		*/
		set<string> sat();
		string frame();
		/*int drag();
		int customer();*/
		/*map<string, int>  drag_parameter();
		map<string, int>  cust_parameter();*/


		/** @brief settings check. */
		void check();

		/** @brief settings help. */
		void help();


	};
}
#endif