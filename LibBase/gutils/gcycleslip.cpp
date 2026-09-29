/**
*
* @verbatim
	History
	 -1.0 GREAT	    2019-01-04 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			  gcycleslip.h
* @brief		  CRS coordinate transform to ACR coordinate
*
* @author         GREAT, Wuhan University
* @version		  1.0.0
* @date			  2019-01-04
*
*/


#include "gcycleslip.h"
using namespace gnut;
great::t_gcycleslip::t_gcycleslip()
{
}

great::t_gcycleslip::t_gcycleslip(t_gsetbase * set, t_glog * log)
{
	_slip_model = dynamic_cast<t_gsetproc*>(set)->slip_model();
	_set = set;
	_log = log;
}

great::t_gcycleslip::~t_gcycleslip()
{
}
