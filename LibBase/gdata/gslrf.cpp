/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-10-27	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gslrf.cpp
 * @brief		The base class used to store slrf2008 file information.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-10-27
 */

#include "gslrf.h"

using namespace std;

namespace great
{
	/** default constructor */
	t_gslrf::t_gslrf()
	{
		id_type(t_gdata::SLRF);
	}

	/** default destructor */
	t_gslrf::~t_gslrf()
	{

	}

	t_gslrf_one& t_gslrf::getOnegslrf(int i)
	{
		return _gslrf[i];
	}

	map<int,t_gslrf_one>& t_gslrf::getgslrf()
	{
		return _gslrf;
	}

	void t_gslrf::setGSLRF(t_gslrf& gslrf)
	{
		_gslrf = gslrf.getgslrf();
	}

	void t_gslrf::setCode(int i,int code)
	{
		_gslrf[i]._code = code;
	}

	void t_gslrf::setDomes(int i, string domes)
	{
		_gslrf[i]._DOmes = domes;
	}

	void t_gslrf::setNameGPS(int i, string name)
	{
		_gslrf[i]._name_gps = name;
	}

	int t_gslrf::getCode(int i)
	{
		return _gslrf[i]._code;
	}

	void t_gslrf::setXsta(int i, int num, double xv)
	{
		_gslrf[i]._xsta[num] = xv;
	}

	int t_gslrf::getSize()
	{
		return _gslrf.size();
	}

	void t_gslrf::setREF(int i,t_gtime now)
	{
		_gslrf[i].now = now;
	}
}