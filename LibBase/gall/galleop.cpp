/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2020-05-14	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		galleop.h
 * @brief		The base class used to store range biase information.
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2020-05-14
 */

#include "gall/galleop.h"

using namespace std;

namespace great
{
	t_galleop::t_galleop()
	{
		id_type(t_gdata::ALLEOP);
	}

	t_galleop::~t_galleop()
	{

	}

	void t_galleop::addeop(t_gtime t, eop eop_one)
	{
		eop_all[t] = eop_one;
	}

	map<t_gtime, eop>* t_galleop::geteop()
	{
		return &eop_all;
	}

	void eop::clear_value()
	{
		this->xpole = 0.0;
		this->ypole = 0.0;
		this->dxpole = 0.0;
		this->dypole = 0.0;
		this->ut1 = 0.0;
		this->dut1 = 0.0;
	}
}