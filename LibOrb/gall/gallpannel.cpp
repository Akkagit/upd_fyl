/**
*
* @verbatim
    History
     -1.0 Yujie Qin   2019-01-06 creat the file.
     -1.1 Wei Zhang   2019-03-27 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		   gallpannel.cpp
* @brief	   bStorage the LEO pannel data(more than one satellite)
*
* @author      Yujie Qin created, Wuhan University
* @version	   1.0.0
* @date		   2019-03-27
*
*/

#include "gall/gallpannel.h"

using namespace std;
namespace great
{
	/** @brief default constructor. */
	t_gallpannel::t_gallpannel()
	{
		id_type(t_gdata::ALLPANNEL);
	}

	/** @brief default destructor. */
	t_gallpannel::~t_gallpannel()
	{}

	/**
	* @brief return one sat pannels.
	* @param[in]   sat         LEO satellite name.
	* @return      pannel data for one LEO.
	*/
	t_gpannel * t_gallpannel::get_sat_pannel(string sat)
	{
		return _sat_pannels.count(sat) > 0 ? &_sat_pannels[sat] : NULL;
	}

	/**
	* @brief push pannel data for one LEO into the data map.
	* @param[in]   sat         LEO satellite name.
	* @param[in]  pannels      pannel data for one LEO.
	*/
	void t_gallpannel::add_sat_pannel(string sat, t_gpannel pannels)
	{
		_sat_pannels[sat] = pannels;
	}

	/** @brief whether the pannel data is empty.
	* @return  bool
	*	@retval   true    pannel data is empty
	*   @retval   false   pannel data is existent
	*/
	bool t_gallpannel::isEmpty()
	{
		if (_sat_pannels.size() == 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	/**
	* @brief return all LEO sat's pannels.
	* @return   all LEO sat's pannel data
	*/
	map<string, t_gpannel>& t_gallpannel::get_sat_pannels()
	{
		return _sat_pannels;
	}
}