/**
*
* @verbatim
    History
     -1.0 Hongmin Zhang 2019-01-06 creat the file.
     -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
	 -1.2 Huang Jiande  2020-09-12 Standard code style
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gallics_leo.cpp
* @brief		    Storage the LEO ICS data(more than one satellite)
*
* @author       Hongmin Zhang, Wuhan University
* @version		1.0.0
* @date		    2019-03-27
*
*/

#include "gall/gallics_leo.h"

using namespace std;

namespace great
{
	/** default constructor */
	t_gallics_leo::t_gallics_leo()
	{
		id_type(t_gdata::ALLICSLEO);
	}

	/** default destructor */
	t_gallics_leo::~t_gallics_leo()
	{
		
	}

	/* -------------
	* @brief get icsleo for one satellite
	* @param [in]  sat_name     LEO satellite name
	* @return  ICS data for one LEO satellite
	*/
	t_gics_leo* t_gallics_leo::get_leo_ics(string sat_name)
	{
		// jdhuang
		if (_map_ics_leo.count(sat_name) == 0)
		{
			if (_log)
			{
				_log->comment(1, "WARNING : t_gallics_leo::ics ==> Can not find the ics of " + sat_name);
			}
			else
			{
				cout << "WARNING : t_gallics_leo::ics ==> Can not find the attitude of " + sat_name << endl;
			}
			return NULL;
		}

		return &_map_ics_leo[sat_name];
		
	}

	/* ----------------------------
	* @brief whether the ics data is empty.
	* @return  bool
	*	@retval   true    ics data is empty
	*   @retval   false   ics data is existent
	*/
	bool t_gallics_leo::isEmpty()
	{
		if (_map_ics_leo.size() == 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	/* ----------------------------
	* @brief push one satellite with ics into the map
	* @param[in]  sat_name   LEO satellite name
	* @param[in]  gicsleo    ICS data for one LEO satellite
	*/
	void t_gallics_leo::add_leo_ics(string sat_name, t_gics_leo& ics_leo)
	{
		_map_ics_leo[sat_name] = ics_leo;
	}

	/* ----------------------------
	* @brief set begin and end time of the ICS data.
	* @param[in]  begtime   begin time
	* @param[in]  endtime   endtime
	*/
	void t_gallics_leo::set_time(const t_gtime& beg_time, const t_gtime& end_time)
	{
		_beg_time = beg_time;
		_end_time = end_time;
	}

	void t_gallics_leo::set_beg_time(const t_gtime & beg_time)
	{
		_beg_time = beg_time;
	}

	void t_gallics_leo::set_end_time(const t_gtime & end_time)
	{
		_end_time = end_time;
	}


	/* ----------------------------
	* @brief return begin time of the ICS data.
	* @return    begin time
	*/
	t_gtime t_gallics_leo::get_beg_time()
	{
		return _beg_time;
	}

	/* ----------------------------
	* @brief return end time of the ICS data.
	* @return    end time
	*/
	t_gtime t_gallics_leo::get_end_time()
	{
		return _end_time;
	}

	/* ----------------------------
	* @brief return LEO satellite name list.
	* @return    name list
	*/
	set<string>	t_gallics_leo::get_leo_sat_list()
	{
		std::set<std::string> leo_list;
		for (auto iter = _map_ics_leo.begin(); iter != _map_ics_leo.end(); iter++)
		{
			leo_list.insert(iter->first);
		}
		return leo_list;
	}


	/* ----------------------------
	* @brief return number of satellites.
	* @return    satellite number
	*/
	int t_gallics_leo::get_leo_sat_num()
	{
		return _map_ics_leo.size();
	}
}