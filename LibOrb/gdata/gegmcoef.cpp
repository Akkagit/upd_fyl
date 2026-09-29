/**
 *
 * @verbatim
	 History
	  -1.0	Jie Li	2019-01-06	creat the file.
	  -1.1  Jie Li	2019-04-01	add code in doxygen style.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gegmcoef.cpp
 * @brief		The base class used to storage EGM file information.
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 */

#include "gdata/gegmcoef.h"

using namespace std;

namespace great
{
	/** @brief constructor. */
	t_gegmcoef::t_gegmcoef()
	{
		gtrace("t_gegmcoef::constructor");
		id_type(t_gdata::EGM);
		_name = "EGM2008";
	}

	/** @brief destructor. */
	t_gegmcoef:: ~t_gegmcoef(){}


	/**
	* @brief add time GM and radius data.
	* @param [in]  t       time of EGM data
	* @param [in]  GM      GM of the Earth
	* @param [in]  radius  radius of the Earth
	* @return	if add data successfully
	*/
	int t_gegmcoef::add_head(t_gtime t, double GM, double radius)
	{
		gtrace("t_gegmcoef::add_head");
		_gmutex.lock();

		_t = t;
		_GM = GM;
		_radius = radius;

		_gmutex.unlock(); return 0;
	}

	/**
	* @brief add EGM order data.
	* @param [in]  order    EGM order
	* @return	if add data successfully
	*/
	int t_gegmcoef::add_order(int order)
	{
		gtrace("t_gegmcoef::add_order");
		_gmutex.lock();

		_order = order;

		_gmutex.unlock(); return 0;
	}

	/**
	* @brief add EGM coefficient data.
	* @param [in]  n    index
	* @param [in]  m    index
	* @param [in]  c    c value
	* @param [in]  s    s value
	* @param [in]  sigma_c    sigma_c value
	* @param [in]  sigma_s    sigma_s value
	* @return	index n
	*/
	int t_gegmcoef::add_egmcoef(int n, int m, double c, double s, double sigma_c, double sigma_s)
	{
		gtrace("t_gegmcoef::add_egmcoef");
		_gmutex.lock();

		t_gravity_coef tmp_coef;
		tmp_coef.c = c;
		tmp_coef.s = s;
		tmp_coef.sigma_c = sigma_c;
		tmp_coef.sigma_s = sigma_s;

		_map_egmcoef[n][m] = tmp_coef;

		_gmutex.unlock(); return n;
	}

	/** @brief whether the EGM coefficien data is empty.
	* @return  bool
	*	@retval   true    EGM coefficien data is empty
	*   @retval   false   EGM coefficien data is existent
	*/
	bool t_gegmcoef::is_empty()
	{
		if (_map_egmcoef.size() == 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	/**
	* @brief get corresponding c value.
	* @param [in]  n    index
	* @param [in]  m    index
	* @return  c value
	*/
	double t_gegmcoef::get_c(int n, int m)
	{
		gtrace("t_gegmcoef::get_c");
		_gmutex.lock();
		double c = 0.0;
		if (_map_egmcoef.count(n) == 0 || _map_egmcoef[n].count(m) == 0)
		{
			_log->comment(0, "ERROR : gravity degreee is out of range! max gravity degree in EGM file is 150");
			_gmutex.unlock();
			return 0.0;
		}
		c = _map_egmcoef[n][m].c;
		_gmutex.unlock();
		return c;
	}

	/**
	* @brief get corresponding s value.
	* @param [in]  n    index
	* @param [in]  m    index
	* @return  s value
	*/
	double t_gegmcoef::get_s(int n, int m)
	{
		gtrace("t_gegmcoef::get_s");
		_gmutex.lock();

		double s = 0.0;
		if (_map_egmcoef.count(n) == 0 || _map_egmcoef[n].count(m) == 0)
		{
			_log->comment(0, "ERROR : gravity degreee is out of range! max gravity degree in EGM file is 150");
			_gmutex.unlock();
			return 0.0;
		}
		s = _map_egmcoef[n][m].s;
		_gmutex.unlock();
		return s;
	}

	/**
	* @brief get reference time of EGM file.
	* @return  mjd of reference time.
	*/
	double t_gegmcoef::mjd_ref()
	{
		return _t.dmjd();
	}

}