/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2020-01-02	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gpreproc.cpp
 * @brief		The base class used to store the preprocessing result of slr observations.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2020-01-02
 */

#ifndef GPREPROCSLR_H
#define GPREPROCSLR_H

#include <iostream>
#include <string.h>
#include "gdata/gdata.h"
#include <vector>

using namespace std;
using namespace gnut;

namespace great
{

	class BASE_LIBRARY_EXPORT t_gpreprocslr
	{
	public:

		t_gpreprocslr();
		virtual ~t_gpreprocslr();

		//void setTime(t_gtime& time);
		void setCode(double code);
		void setTwowayTime(double sec);
		void setAtmPressure(double pressure);
		void setDryTem(double drytem);
		void setWave(double wave);
		void setAtmHumidity(double humidity);
		void setSta(string name);
		void setSat(string sat);

		//t_gtime _now;
		int _code;
		double _two_way_time;
		double _atm_pressure;
		double _dry_tem;
		double _wave;
		double _atm_humidity;
		string _sta_name;
		string _sat_name;

	protected:


	};






}

#endif