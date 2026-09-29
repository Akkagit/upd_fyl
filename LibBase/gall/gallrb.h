/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2020-03-06	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gallrb.h
 * @brief		The base class used to store range biase information.
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2020-03-06
 */

#ifndef GALLRB_H
#define GALLRB_H

#include <iostream>
#include <string.h>
#include "gdata/gdata.h"

using namespace std;
using namespace gnut;

namespace great
{
	class BASE_LIBRARY_EXPORT t_gallrb : public t_gdata
	{
	public:
		t_gallrb();
		virtual ~t_gallrb();

		void addrb(string name, string sat, double RB);
		void addstd(string name, string sat, double stdx);

		map<string, map<string, double>>& getrb();

		map<string, map<string, double>>* getstd();

	private:
		map<string, map<string, double>> _Range_Bias;
		map<string, map<string, double>> _Range_Bias_stdx;
	};
}

#endif