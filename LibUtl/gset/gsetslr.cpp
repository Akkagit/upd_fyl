/**
*
* @verbatim
	History
	 -1.0 Hongmin Zhang   2019-11-20 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gsetslr.cpp
* @brief        part of SLR settings from XML file
*
* @author       Hongmin Zhang, Wuhan University
* @version		1.0.0
* @date		    2019-11-20
*
*/


#include "gset/gsetslr.h"
#include "gutils/gleo.h"
#include "gutils/gmutex.h"

using namespace std;

namespace great
{
	/** @brief constructor. */
	t_gsetslr::t_gsetslr()
		: t_gsetbase()
	{
		_set.insert(XMLKEY_SLR);
	}

	/** @brief destructor. */
	t_gsetslr::~t_gsetslr()
	{

	}

	/** @brief settings check. */
	void t_gsetslr::check()
	{

	}

	/** @brief settings help. */
	void t_gsetslr::help()
	{

	}

	/**
	* @brief return all LEO satellites in config file(string)
	* @return  all LEO satellites in config file(string)
	*/
	string t_gsetslr::getSLRXmlSet()
	{
		try
		{
			_gmutex.lock();

			string sats = t_gsetbase::_strval("SLR", "type");
			_gmutex.unlock();
			return sats;
		}
		catch (...)
		{
			cout << "ERROR: t_gsetleo::sat  unknown mistake" << endl;
			_gmutex.unlock();
			return string();
		}
	}


} // namespace