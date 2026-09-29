/**
*
* @verbatim
    History
     -1.0 Yujie Qin   2019-01-06 creat the file.
     -1.1 Wei Zhang   2019-03-27 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gsetleo.cpp
* @brief        part of LEO settings from XML file
*
* @author       Yujie Qin, Wuhan University
* @version		1.0.0
* @date		    2019-03-27
*
*/


#include "gset/gsetleo.h"
#include "gutils/gleo.h"
#include "gutils/gmutex.h"

using namespace std;

namespace great
{
	/** @brief constructor. */
	t_gsetleo::t_gsetleo()
	: t_gsetbase()
	{
		_set.insert(XMLKEY_LEO);
	}
	
	/** @brief destructor. */
	t_gsetleo::~t_gsetleo()
	{
		
	}
	
	/** @brief settings check. */
	void t_gsetleo::check()
	{
		
	}
	
	/** @brief settings help. */
	void t_gsetleo::help()
	{
		
	}
	
	/**
	* @brief return all LEO satellites in config file(string)
	* @return  all LEO satellites in config file(string)
	*/
	set<string> t_gsetleo::sat()
	{
		try
		{
			_gmutex.lock();

			set<string> allsat;
			set<string>::const_iterator it;
			set<string> sats = t_gsetbase::_setval("LEO", "sat");
			for (it = sats.begin(); it != sats.end(); ++it)
			{
				if (IS_SUPPORTED_LEO(*it))
				{
					allsat.insert(*it);
				}
			}

			_gmutex.unlock();
			return allsat;
		}
		catch (...)
		{
			cout << "ERROR: t_gsetleo::sat  unknown mistake" << endl;
			_gmutex.unlock();
			return set<string>();
		}
	}
	
	/**
	* @brief return all LEO satellites in config file(GLEOLONG enum)
	* @return  all LEO satellites in config file(GLEOLONG enum)
	*/
	vector<GLEOLONG> t_gsetleo::esat()
	{
		try
		{
			_gmutex.lock();

			vector<GLEOLONG> allsat;
			set<string>::const_iterator it;
			set<string> sats = t_gsetbase::_setval("LEO", "sat");
			for (it = sats.begin(); it != sats.end(); ++it)
			{
				if (IS_SUPPORTED_LEO(*it))
				{
					allsat.push_back(str2long(*it));
				}
			}

			_gmutex.unlock();
			return allsat;
		}
		catch (...)
		{
			cout << "ERROR: t_gsetleo::esat  unknown mistake" << endl;
			_gmutex.unlock();
			return vector<GLEOLONG>();
		}
	}
	
} // namespace