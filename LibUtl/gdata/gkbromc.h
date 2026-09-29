/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2020-03-20	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gkbromc.h
 * @brief		The base class used to store kbr caliboration result.
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2020-03-20
 */

#ifndef GKBROMC_H
#define GKBROMC_H

#include <iostream>
#include <string.h>
#include "gutils/ExportUtl.h"
#include "gutils/gmutex.h"
#include "gdata/gdata.h"
#include <vector>

using namespace std;
using namespace gnut;

namespace great
{

	class UTL_LIBRARY_EXPORT t_gkbromc: public t_gdata
	{
	public:

		t_gkbromc();
		virtual ~t_gkbromc();

		void setOMC(t_gtime t, double omc);
		map<t_gtime, double>& getOMC();

		map<t_gtime, double> _kbr_omc;

	private:
		

	};






}

#endif