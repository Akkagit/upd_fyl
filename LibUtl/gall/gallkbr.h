/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2020-03-20	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gallkbr.h
 * @brief		The base class used to store the information of kbr observations.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2020-03-20
 */

#ifndef GALLKBR_H
#define GALLKBR_H

#include <iostream>
#include <string.h>
#include "gdata/gdata.h"
#include "gutils/ExportUtl.h"
#include <vector>

using namespace std;
using namespace gnut;

namespace great
{
	class t_kbr_data
	{
	public:
		double bias_range;
		double range_rate;
		double iono;
		double ltime_cor;
		double ltime_rate;
		double ant_cor;
		double ant_rate;
		int flag;
	};

	class UTL_LIBRARY_EXPORT t_gallkbr : public t_gdata
	{
	public:
		t_gallkbr();
		virtual ~t_gallkbr();

		void setHead(t_gtime& beg, t_gtime& end, string sat1, string sat2, double intv);
		void setData(t_gtime& t, t_kbr_data& one_data);

		t_gtime&  getBeg();
		t_gtime&  getEnd();
		vector<string>& getSat();
		double   getIntv();
		t_kbr_data*   getkbr(t_gtime& t);
		void setOMC(t_gtime t,double omc);
		map<t_gtime, double>& getOMC();

	protected:
		t_gtime   _start;
		t_gtime   _end;
		string    _sat_1;
		string    _sat_2;
		vector<string> sat;
		double    _interval;
		map<t_gtime, t_kbr_data>  _kbr_data;
		map<t_gtime, double>  _kbr_omc;
	};
}

#endif