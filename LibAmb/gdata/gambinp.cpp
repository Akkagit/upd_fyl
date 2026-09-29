/**
*
* @verbatim
History
	-1.0 Jie Li  2019-10-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambinp.cpp
* @brief	Storage the ambinp files' data include par_amb and neq
*
* @author   Jie Li, Wuhan University
* @version	1.0.0
* @date		2019-10-25
*
*/

#include "gdata/gambinp.h"
#include "gutils/gtime.h"

using namespace std;
namespace great 
{

	t_gambinp_head::t_gambinp_head()
	{
		num_amb  = 0;
		interval = 0.0;
	}

	t_gambinp_head::~t_gambinp_head()
	{
	}



	t_gambinp::t_gambinp()
	{
		gtrace("t_gambinp::constructor");
		id_type(t_gdata::AMBINP);
	}

	t_gambinp::~t_gambinp()
	{
	}

	void t_gambinp::add_hdinfo(t_gambinp_head * hd)
	{
		_ambinp_hd = hd;
	}

	void t_gambinp::add_ambpar(t_gpar ambpar)
	{
		_amb.push_back(ambpar);
	}

	void t_gambinp::add_neq(string parname, vector<double>& q)
	{
		_par_list.push_back(parname);
		_neq.push_back(q);
	}

	void t_gambinp::get_hdinfo(t_gambinp_head * hd)
	{
		hd = _ambinp_hd;
	}

	void t_gambinp::get_allambpar(vector<t_gpar>& ambpar)
	{
		ambpar = _amb;
	}

	void t_gambinp::get_neq(vector<string>& parlist, vector<vector<double>>& neq)
	{
		parlist = _par_list;
		neq = _neq;
	}

}//namespace