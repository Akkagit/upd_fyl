/**
*
* @verbatim
	History
	 -1.0 Yujie Qin		2019-12-02 creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gsetorbsp3.h
* @brief		orb to sp3 parameters' set from XML
*
* @author       Yujie Qin	, Wuhan University
* @version		1.0.0
* @date		    2019-12-02
*
*/

#ifndef GSETORBSP3_H
#define GSETORBSP3_H

#define XMLKEY_ORBSP3 	           "orbsp3"
#define XMLKEY_ORBSP3_SAT          "sat"
#define XMLKEY_ORBSP3_INTV         "intv"
#define XMLKEY_ORBSP3_FRAME        "frame"

#include "gset/gsetbase.h"
#include "gutils/ExportOI.h"

namespace gnut
{
	class OI_LIBRARY_EXPORT t_gsetorbsp3 : public virtual t_gsetbase
	{
	public:
		/** @brief constructor. */
		t_gsetorbsp3();

		/** @brief destructor. */
		virtual ~t_gsetorbsp3();

 	    /**
		* @brief get the satellites that need calculation.
		* @return the sat
		*/
		vector<string> sat();
		double intv();
		string frame();

		/** @brief settings check. */
		void check();

		/** @brief settings help. */
		void help();


	};
}
#endif