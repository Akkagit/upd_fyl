/**
*
* @verbatim
    History
     -1.0 Yujie Qin     2019-01-06 creat the file.
     -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		        gallattitude.cpp
* @brief		    Storage the LEO attitude data(more than one satellite)
*
* @author           Yujie Qin created, Wuhan University
* @version		    1.0.0
* @date		        2019-03-27
*
*/

#include "gall/gallattitude.h"

using namespace std;

namespace great
{
	/** default constructor */
	t_gallatt::t_gallatt()
	{
		id_type(t_gdata::AllATTITUDE);
	}

	/** default destructor */
	t_gallatt:: ~t_gallatt()
	{

	}

	/**
	* @brief push attitude data for one LEO into the allattitude data map.
	* @param [in]  sat         LEO satellite name.
	* @param [in]  attitude    attitude data for one LEO.
	*/
	void t_gallatt::add_attitude(string sat, t_gattitude& attitude)
	{
		if (_attitudes.find(sat) != _attitudes.end())
		{
			t_gattitude& old = _attitudes[sat];
			map<t_gtime, t_attitude_record>* newatt = attitude.attitude();
			t_gtime begold = old.beg();
			t_gtime endold = old.end();
			t_gtime begnew = attitude.beg();
			t_gtime endnew = attitude.end();
			if (begold > begnew)
			{
				begold = begnew;
			}
			if (endold < endnew)
			{
				endold = endnew;
			}
			old.set_time(begold, endold);
			auto oneatt = newatt->begin();
			auto lastatt = newatt->end();
			while (oneatt != lastatt)
			{
				old.set_attitude_record(oneatt->first, oneatt->second);
				oneatt++;
			}
		}
		else
		{
			_attitudes[sat] = attitude;
		}
		
	}

	/** @brief whether the attitude data is empty.
	* @return  bool
	*	@retval   true    attitude data is empty
	*   @retval   false   attitude data is existent
	*/
	bool t_gallatt::isEmpty()
	{
		if (_attitudes.size() == 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	/**
	* @brief return attitude data map for all LEO satellites.
	* @return    attitude data map for all LEO satellites.
	*/
	map<string,t_gattitude>* t_gallatt::attitudes()
	{
		return &_attitudes;
	}

	/**
	* @brief return attitude data for one LEO satellites.
	* @return    attitude data for one LEO satellite.
	*/
	t_gattitude * t_gallatt::attitude(string sat)
	{
		if (_attitudes.count(sat) == 0)
		{
			if (_log)
			{
				_log->comment(1, "WARNING : t_gallattitude::attitude ==> Can not find the attitude of " + sat);
			}
			else
			{
				cout << "WARNING : t_gallattitude::attitude ==> Can not find the attitude of " + sat << endl;
			}
			return NULL;
		}
		return &_attitudes[sat];
	}
}