/**
*
* @verbatim
	History
	 -1.0 Hongmin Zhang		2019-01-05 creat the file.
	 -1.1 Jie Li			2019-04-01 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gallics_nav.cpp
* @brief		Storage the GNS ICS data(more than one satellite)
*
* @author       Hongmin Zhang, Wuhan University
* @version		1.0.0
* @date		    2019-01-05
*
*/

#include "gall/gallics_nav.h"

using namespace std;

namespace great
{
	/** default constructor */
	t_gallics_nav::t_gallics_nav()
	{
		id_type(t_gdata::ALLICS);
	}

	/** @brief whether the ics data is empty.
	* @return  bool
	*	@retval   true    ics data is empty
	*   @retval   false   ics data is existent
	*/
	bool t_gallics_nav::isEmpty()
	{
		if (_map_nav_ics.size() == 0)
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
	* @param[in]  sat_name   GND satellite name
	* @param[in]  gics		 ICS data for one GNS satellite
	*/
	void t_gallics_nav::add_nav_ics(string sat_name, t_gics_nav gics)
	{
		//_map_ics[sat_name] = gics;

		_map_nav_ics.insert(pair<string, t_gics_nav>(sat_name, gics));
	}
	/* ----------------------------
	* @brief delete ICS data for one GNS satellite in the data map.
	* @param[in]  sat_name   GND satellite name
	*/
	void t_gallics_nav::del_nav_ics(string sat)
	{
		if (_map_nav_ics.count(sat) > 0)
		{
			_map_nav_ics.erase(sat);
		}
	}

	
	/* -------------
	* @brief get ics for one satellite
	* @param [in]  sat_name     GNS satellite name
	* @return  ICS data for one GNS satellite
	*/
	t_gics_nav* t_gallics_nav::get_nav_ics(string sat_name)
	{
		if (_map_nav_ics.count(sat_name) > 0)
		{
			return &_map_nav_ics[sat_name];
		}
		else
		{
			return NULL;
		}
	}
}
