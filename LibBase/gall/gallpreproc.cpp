/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-12-24	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gallpreproc.cpp
 * @brief		The base class used to store the preprocessing result of slr observations.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-12-24
 */

#include "gallpreproc.h"
#include "gdata/gdata.h"

using namespace std;

namespace great
{

	/** default constructor */
	t_gallpreproc::t_gallpreproc()
	{
		id_type(t_gdata::ALLPREPROC);
	}

	/** default destructor */
	t_gallpreproc::~t_gallpreproc()
	{

	}

	void t_gallpreproc::setSLROBS(string sat,t_gtime& t, t_gpreprocslr& onedata)
	{

		_slrobs[sat][t] = onedata;
	}
	
	map<int,string>& t_gallpreproc::getStaCode()
	{
		int tab;
		
		for (auto iter = _slrobs.begin(); iter != _slrobs.end(); iter++)
		{
			for (auto itersat = iter->second.begin(); itersat != iter->second.end(); itersat++)
			{
				tab = 0;
				for (auto iter_sta = sta_code.begin(); iter_sta != sta_code.end(); iter_sta++)
				{
					if (iter_sta->first == itersat->second._code)
					{
						tab = 1;
					}
				}
				if (tab == 0)
				{
					sta_code[itersat->second._code] = itersat->second._sta_name;
				}
			}
			
		}

		return sta_code;
	}

	map<string,vector<string>>& t_gallpreproc::getSatSta()
	{
		
		for (auto iter = _slrobs.begin(); iter != _slrobs.end(); iter++)
		{
			for (auto itersat = iter->second.begin(); itersat != iter->second.end(); itersat++)
			{
				if (!sat_sta.empty())
				{
					if (sat_sta.find( iter->first) == sat_sta.end())
					{
						sat_sta[iter->first].push_back(itersat->second._sta_name);
					}
					if(find(sat_sta[iter->first].begin(), sat_sta[iter->first].end(), itersat->second._sta_name) == sat_sta[iter->first].end())
					{
						sat_sta[iter->first].push_back(itersat->second._sta_name);
					}
				}
				else
				{
					sat_sta[iter->first].push_back(itersat->second._sta_name);
				}
			}
		}

		return sat_sta;

	}

	t_gpreprocslr t_gallpreproc::getSLROBS(string sat,t_gtime t)
	{
		return _slrobs[sat][t];
	}

	map<string,map<t_gtime, t_gpreprocslr>>* t_gallpreproc ::getallSLROBS()
	{
		return &_slrobs;
	}

}