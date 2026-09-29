/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-11-06	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gallresult.cpp
 * @brief		The base class used to store slr caliboration result.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-11-06
 */

#include "gallresult.h"
#include "gdata/gdata.h"

using namespace std;

namespace great
{

	/** default constructor */
	t_gresult::t_gresult()
	{
		id_type(t_gdata::ALLRESULT);
	}

	/** default destructor */
	t_gresult::~t_gresult()
	{

	}

	map<int, t_gresult_one>* t_gresult::getMap()
	{
		return &_map_result;
	}

	map<string,map<int, t_gresult_one>>* t_gresult::getMapSat()
	{
		return &_map_result_sat;
	}

	void t_gresult::clearmap()
	{
		_map_result.clear();
	}

	void t_gresult::setOneResult(int i, t_gresult_one* oneresult)
	{
		_map_result[i] = *oneresult;
	}

	void t_gresult::setOnesat(string sat, map<int, t_gresult_one>& mapresult)
	{
		_map_result_sat[sat] = mapresult;
	}

}