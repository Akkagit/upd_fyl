/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2020-03-20	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gkbromc.cpp
 * @brief		The base class used to store kbr caliboration result.
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2020-03-20
 */
#include "gdata/gkbromc.h"
using namespace std;

namespace great
{
	t_gkbromc::t_gkbromc()
	{
		id_type(t_gdata::KBROMC);
	}
	t_gkbromc::~t_gkbromc()
	{

	}

	void t_gkbromc::setOMC(t_gtime t,double omc)
	{
		_kbr_omc[t] = omc;
	}

	map<t_gtime, double>& t_gkbromc::getOMC()
	{
		return _kbr_omc;
	}

}
