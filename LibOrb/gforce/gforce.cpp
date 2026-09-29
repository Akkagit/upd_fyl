/**
 *
 * @verbatim
	 History
	  -1.0 Jiande Huang  2019-01-14 creat the file.
	  -1.1 Jie Li       2019-04-01 add code in doxygen style.
	  -1.2 Jiande Huang  2019-04-01 Adding exception throwing mechanism for all functions
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gforce.cpp
 * @brief		The base class for all the force models.
 *
 * @author      Jiande Huang, Wuhan University
 * @version		1.0.0
 * @date		2019-01-14
 */
#include "gforce/gforce.h"
#include "gutils/ginfolog.h"

namespace great
{
	/** @brief Constructor */
	t_gforce::t_gforce()
	{
		_acc.resize(3);
	}

	/** @brief Constructor */
	t_gforce::t_gforce(t_glog * log)
	{
		if (!log)
		{
			write_log(log, 0, "t_gforce::t_gforce",
				"ERROR : the log point is empty");
		}
		else
		{
			_glog = log;
			_acc.resize(3);
		}
	}

	/** @brief Destructor */
	t_gforce::~t_gforce()
	{

	}

	/**
	* @brief an empty function for subclass to overrider for force calculation.
	* @param[in]   argument     force calculate argument
	* @return	bool
	*/
	bool t_gforce::calculate(t_gforce_calarg * argument)
	{
		return false;
	}

	/**
	* @brief get the acceleration matrix.
	* @return	acceleration matrix
	*/
	ColumnVector& t_gforce::acc()
	{
		return _acc;
	}

	FORCE_MODEL t_gforce::force_model()
	{
		return _force_model;
	}

	/**
	* @brief get partial matrix.
	* @return   partial matrix
	*/
	Matrix& t_gforce::partial()
	{
		return _partial;
	}
}

