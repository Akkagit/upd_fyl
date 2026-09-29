/**
*
* @verbatim
History
-1.0 Hongmin Zhang  2020-03-20 creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		       gcfg_kbr.h
* @brief	       KBR
* @author          Hongmin Zhang, Wuhan University
* @version	       1.0.0
* @date		       2020-03-20
*
*/

#include <iomanip>
#include <sstream>

#include "gset/gcfg_kbr.h"

using namespace std;
using namespace pugi;
using namespace gnut;

namespace great {

	// Constructor
	// ----------
	t_gcfg_kbr::t_gcfg_kbr()
		: t_gsetbase(),
		t_gsetgen(),
		t_gsetinp(),
		t_gsetout(),
		t_gsetgnss()
	{
		
		_IFMT_supported.insert(ORB_INP);
		_IFMT_supported.insert(KBR_INP);


		_OFMT_supported.insert(LOG_OUT);
		_OFMT_supported.insert(PPP_OUT);
		_OFMT_supported.insert(KBROMC_OUT);
		_OFMT_supported.insert(RECCLK_OUT);
		_OFMT_supported.insert(SATCLK_OUT);
		_OFMT_supported.insert(PRE_OUT);
		_OFMT_supported.insert(RESULT_OUT);
	}


	// Destructor
	// ----------
	t_gcfg_kbr::~t_gcfg_kbr()
	{}


	// settings check
	// ----------
	void t_gcfg_kbr::check()
	{
		t_gsetgen::check();
		t_gsetinp::check();
		t_gsetout::check();
		t_gsetgnss::check();
	}

	// settings help
	// ----------
	void t_gcfg_kbr::help()
	{
		t_gsetbase::help_header();
		t_gsetgen::help();
		t_gsetinp::help();
		t_gsetout::help();
		t_gsetgnss::help();
	}

} // namespace
