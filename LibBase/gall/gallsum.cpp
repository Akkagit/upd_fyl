/**
*
* @verbatim
	History
	 -1.0 Hongmin Zhang 2020-01-14 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gallsum.cpp
* @brief		    Storage the sum of residuals.
*
* @author       Hongmin Zhang, Wuhan University
* @version		1.0.0
* @date		    2020-01-14
*
*/

#include "gall/gallsum.h"

using namespace std;

namespace great
{
	/** default constructor */
	t_gallsum::t_gallsum()
	{
		id_type(t_gdata::ALLSUM);
	}

	/** default destructor */
	t_gallsum::~t_gallsum()
	{

	}

	void t_gallsum::setSat(vector<string>& satlist)
	{
		_sat_list = satlist;
	}

	void t_gallsum::setSite(vector<string>& sitelist)
	{
		_site_list = sitelist;
	}

	void t_gallsum::setRMS(map<t_gobscombtype, map<int, map<int, double>>>& rms)
	{
		_rmstable = rms;
	}

	void t_gallsum::setRMS_SLR(map<string, map<int, map<int, double>>>& rms)
	{
		_rmstable_slr = rms;
	}

	void t_gallsum::setNEPO(map<int, map<int, int>>& nepo)
	{
		_nepotbl = nepo;
	}

	void t_gallsum::setNEPO_SLR(map<int, map<int, int>>& nepo)
	{
		_nepotbl_slr = nepo;
	}

	vector<string>& t_gallsum::getSat()
	{
		return _sat_list;
	}

	vector<string>& t_gallsum::getSite()
	{
		return _site_list;
	}

	map<t_gobscombtype, map<int, map<int, double>>>& t_gallsum::getRMS()
	{
		return _rmstable;
	}

	map<int, map<int, int>>& t_gallsum::getNEPO()
	{
		return _nepotbl;
	}

	map<string, map<int, map<int, double>>>& t_gallsum::getRMS_SLR()
	{
		return _rmstable_slr;
	}

	map<int, map<int, int>>& t_gallsum::getNEPO_SLR()
	{
		return _nepotbl_slr;
	}

}
