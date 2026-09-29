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

#include "gall/gallrb.h"

using namespace std;

namespace great
{
	t_gallrb::t_gallrb()
	{
		
		id_type(t_gdata::ALLRB);
		
	}

	t_gallrb::~t_gallrb()
	{

	}

	void t_gallrb::addrb(string name,string sat, double RB)
	{
		_Range_Bias[sat][name] = RB;
	}

	void t_gallrb::addstd(string name, string sat, double stdx)
	{
		_Range_Bias_stdx[sat][name] = stdx;
	}

	map<string,map<string, double>>& t_gallrb:: getrb()
	{
		return _Range_Bias;
	}

	map<string, map<string, double>>* t_gallrb::getstd()
	{
		return &_Range_Bias_stdx;
	}

}