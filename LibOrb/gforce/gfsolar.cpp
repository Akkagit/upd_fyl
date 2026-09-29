/**
 *
 * @verbatim
	 History
	  -1.0 Jiande Huang		2019-01-14 creat the file.
	  -1.1 Jie Li			2019-04-01 add code in doxygen style.
   @endverbatim
 *
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gfsolar.cpp
 * @brief		The base class for all the solar force models.
 *
 * @author      Jiande Huang, Wuhan University
 * @version		1.0.0
 * @date		2019-01-14
 */

#include "gforce/gfsolar.h"
#include "gutils/ginfolog.h"


namespace great
{
	/** @brief Constructor */
	t_gfsolar::t_gfsolar(t_glog * log) : t_gforce(log)
	{
		if (!log)
		{
			write_log(log, 0, "t_gfsolar::t_gfsolar",
				"ERROR : the log point is empty");
		}
	}

	/** @brief	Destructor */
	t_gfsolar::~t_gfsolar()
	{

	}

}

