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

#include "gallkbr.h"
#include "gdata/gdata.h"

using namespace std;

namespace great
{

	/** default constructor */
	t_gallkbr::t_gallkbr()
	{
		id_type(t_gdata::ALLKBR);
	}

	/** default destructor */
	t_gallkbr::~t_gallkbr()
	{

	}

	void t_gallkbr::setHead(t_gtime& beg, t_gtime& end, string sat1, string sat2, double intv)
	{
		_start = beg;
		_end = end;
		_sat_1 = sat1;
		_sat_2 = sat2;
		_interval = intv;
	}

	void t_gallkbr::setData(t_gtime& t, t_kbr_data& one_data)
	{
		_kbr_data[t] = one_data;
	}

	void t_gallkbr::setOMC(t_gtime t, double omc)
	{
		_kbr_omc[t] = omc;
	}

	map<t_gtime,double>& t_gallkbr::getOMC()
	{
		return _kbr_omc;
	}

	t_gtime& t_gallkbr::getBeg()
	{
		return _start;
	}

	t_gtime& t_gallkbr::getEnd()
	{
		return _end;
	}

	vector<string>& t_gallkbr::getSat()
	{
		sat.push_back(_sat_1);
		sat.push_back(_sat_2);
		return sat;
	}

	t_kbr_data* t_gallkbr::getkbr(t_gtime& t)
	{
		return &_kbr_data[t];
	}

	double t_gallkbr::getIntv()
	{
		return _interval;
	}
}