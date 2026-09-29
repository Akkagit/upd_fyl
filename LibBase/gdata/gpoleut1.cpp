/**
 *
 * @verbatim
	 History
	  -1.0	Yujie Qin	2019-01-06	creat the file.
	  -1.1  Jie Li		2019-04-01	add code in doxygen style.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		    gpoleut1.cpp
 * @brief		The base class used to storage the poleut1 data.
 *
 * @author      Yujie Qin, Wuhan University
 * @version		1.0.0
 * @date		    2019-01-06
 */


#include "gdata/gpoleut1.h"

using namespace std;
namespace great {

	/** @brief constructor. */
	t_gpoleut1::t_gpoleut1()
		:t_gdata()
	{
		id_type(t_gdata::ALLPOLEUT1);
	}

	/**
	* @brief set begin and end time for poleut1 data.
	* @param[in]   beg		begin time of the poleut1 data
	* @param[in]   end		end time of the poleut1 data
	*/
	void t_gpoleut1::setBegEndTime(int beg, int end)
	{
		_beg_time = beg;
		_end_end = end;
	}

	/**
	* @brief add poleut1 data.
	* @details add one poleut1 data record to _poleut1data and set timetype and interval.
	*
	* @param[in]   mjdtime		(mjd)time
	* @param[in]   data			one poleut1 data record(xpole,ypole and so on)
	* @param[in]   mode			type of UT1(such as UT1R)
	* @param[in]   intv			interal(unit:day)
	*/
	void t_gpoleut1::setEopData(t_gtime mjdtime, map<string, double> data, string mode, double intv)
	{
		_poleut1_data[mjdtime] = data;
		_UT1_mode = mode;
		_intv = intv;
	}

	/**
	* @brief get corresponding poleut1 data by time.
	* @param[in]   mjdtime		(mjd)time
	* @param[out]   data			one poleut1 data record(xpole,ypole and so on)
	*/
	void t_gpoleut1::getEopData(t_gtime mjdtime, map<string, double>& data)
	{
		if (_poleut1_data.find(mjdtime) != _poleut1_data.end())
		{
			data = _poleut1_data[mjdtime];
		}
		// interpolation
		else
		{
			auto iter_right = _poleut1_data.upper_bound(mjdtime);
			if (iter_right == _poleut1_data.begin() || iter_right == _poleut1_data.end())
			{
				cout << "interpolation poleut1 wrong!!" << endl;
				throw runtime_error("interpolation poleut1 wrong!!");
			}
			auto iter_left = iter_right--;

			double alpha = (mjdtime - iter_left->first) / (iter_right->first - iter_left->first);
			data.clear();
			for (auto iter = iter_left->second.begin(); iter != iter_left->second.end(); iter++)
			{
				double temp = iter_left->second[iter->first] + alpha * (iter_right->second[iter->first] - iter_left->second[iter->first]);
				data.insert(make_pair(iter->first, temp));
			}
		}
	}


	/** @brief whether the poleut1 data is empty.
	* @return  bool
	*	@retval   true   poleut1 data is empty
	*   @retval   false  poleut1 data is existent
	*/
	bool t_gpoleut1::isEmpty()
	{
		if (_poleut1_data.size() == 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}