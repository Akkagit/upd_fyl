/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-10-27	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gallnpt.cpp
 * @brief		The base class used to store npt file information.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-10-27
 */

#include "gall/gallnpt.h"

using namespace std;

namespace great
{
	/** default constructor */
	t_gallnpt::t_gallnpt()
	{
		id_type(t_gdata::ALLNPT);
	}

	/** default destructor */
	t_gallnpt::~t_gallnpt()
	{

	}

	void t_gallnpt::setAllNpt(string sat,int& now, t_gnpt& gnpt)
	{
		_gallnpt[sat][now] = gnpt;
	}

	t_gnpt& t_gallnpt::getAllNpt(string sat,int& now)
	{
		return _gallnpt[sat][now];
	}

	int t_gallnpt::getSize(string sat)
	{
		return _gallnpt[sat].size();
	}
}