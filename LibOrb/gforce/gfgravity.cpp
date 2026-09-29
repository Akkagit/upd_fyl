/**
 *
 * @verbatim
	 History
	  -1.0 Jie Li  2019-01-06 creat the file.
	  -1.1 Jie Li  2019-03-31 add code in doxygen style
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gfgravity.cpp
 * @brief		The base class for gravity.
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 */
#include "gforce/gfgravity.h"
#include "gutils/ginfolog.h"

namespace great
{
	/** @brief constructor. */
	t_gfgravity::t_gfgravity(t_glog * log) : t_gforce(log)
	{
		if (!log)
		{
			write_log(log, 4, "t_gfgravity::t_gfgravity",
				"ERROR : the log point is empty");
		}
	}

}