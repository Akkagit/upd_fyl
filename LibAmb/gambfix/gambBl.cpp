/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambBl.cpp
* @brief	function for ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#include "gambfix/gambBlRec.h"
#include "gset/gsetamb.h"
#include "gset/gsetrec.h"
#include <algorithm>
#include "gambBl.h"

namespace great
{
	t_gamb_baseline::t_gamb_baseline()
	{
	}

	t_gamb_baseline::t_gamb_baseline(string obj1, string obj2, double distance = 0.0, int nsys = 1)
	{
		this->obj1 = obj1;
		this->obj2 = obj2;
		this->nsys = nsys;
		this->distance = distance;
	}

	t_gamb_baseline::~t_gamb_baseline()
	{
	}

	pair<string, string> t_gamb_baseline::bl_pair() const
	{
		return make_pair(this->obj1, this->obj2);
	}

	bool t_gamb_baseline::isSame(const string& obj_1st, const string& obj_2nd)
	{
		if (obj_1st == obj1 && obj_2nd == obj2) return true;
		return false;
	}

}

