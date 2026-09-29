/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2020-05-14	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		galleop.h
 * @brief		The base class used to store eop information.
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2020-05-14
 */

#ifndef GALLEOP_H
#define GALLEOP_H

#include <iostream>
#include <string.h>
#include "gdata/gdata.h"

using namespace std;
using namespace gnut;

namespace great
{
	class BASE_LIBRARY_EXPORT eop
	{
	public:
		double xpole;
		double ypole;
		double dxpole;
		double dypole;
		double ut1;
		double dut1;
	
		void clear_value();
	};

	class BASE_LIBRARY_EXPORT t_galleop : public t_gdata
	{
	public:
		t_galleop();
		virtual ~t_galleop();

		void addeop(t_gtime t, eop eop_one);

		map<t_gtime, eop>* geteop();

	private:
		map<t_gtime, eop> eop_all;

	};
}
#endif