/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-12-24	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gallpreproc.h
 * @brief		The base class used to store the preprocessing result of slr observations.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-12-24
 */

#ifndef GALLPREPROC_H
#define GALLPREPROC_H

#include <iostream>
#include <string.h>
#include "gdata/gdata.h"
#include "gdata/gpreprocslr.h"
#include <vector>

using namespace std;
using namespace gnut;

namespace great
{
	class BASE_LIBRARY_EXPORT t_gallpreproc : public t_gdata
	{
	public:
		t_gallpreproc();
		virtual ~t_gallpreproc();

		void setSLROBS(string sat, t_gtime& t, t_gpreprocslr& onedata);

	    t_gpreprocslr getSLROBS(string sat, t_gtime t);

		map<string,map<t_gtime, t_gpreprocslr>>* getallSLROBS();

		map<int,string>& getStaCode();

		map<string, vector<string>>& getSatSta();

		//void setCode(int code);

	protected:
		map<string,map<t_gtime, t_gpreprocslr>> _slrobs;
		map<int, string> sta_code;
		map<string, vector<string>> sat_sta;
		//int  _code;
	};
}

#endif