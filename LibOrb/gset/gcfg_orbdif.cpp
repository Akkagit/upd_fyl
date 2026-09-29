/**
*
* @verbatim
	History
	 -1.0 Jie Li    2019-09-05  creat the file
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gcfg_orbdif.cpp
* @brief		orbit differ control set from XML
*
* @author       Jie Li , Wuhan University
* @version		1.0.0
* @date		    2019-09-05
*
*/


#include "gcfg_orbdif.h"

namespace great {

	t_gcfg_orbdif::t_gcfg_orbdif()
		: t_gsetbase(),
		t_gsetgen(),
		t_gsetinp(),
		t_gsetout(),
		t_gsetgnss(),
		t_gsetrec(),
		t_gsetorbdif()
	{
		_IFMT_supported.insert(SP3_INP);
		_IFMT_supported.insert(ORB_INP);
		_IFMT_supported.insert(POLEUT1_INP);

		_OFMT_supported.insert(ORBDIF_OUT);
	}

	t_gcfg_orbdif::~t_gcfg_orbdif()
	{

	}

	void t_gcfg_orbdif::check()
	{
		t_gsetgen::check();
		t_gsetinp::check();
		t_gsetout::check();
		t_gsetgnss::check();
		t_gsetrec::check();
		t_gsetorbdif::check();
	}

	void t_gcfg_orbdif::help()
	{
		t_gsetbase::help_header();
		t_gsetgen::help();
		t_gsetinp::help();
		t_gsetout::help();
		t_gsetgnss::help();
		t_gsetrec::help();
		t_gsetorbdif::help();
		t_gsetbase::help_footer();
	}

}

