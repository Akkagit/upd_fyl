/**
*
* @verbatim
    History
     -1.0 Yujie Qin  2019-10-10 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		       gcfg_sp3orb.cpp
* @brief	       sp3 to orb set from XML
*
* @author          Yujie Qin, Wuhan University
* @version	       1.0.0
* @date		       2019-10-10
*
*/

#include <iomanip>
#include <sstream>

#include "gset/gcfg_sp3orb.h"

using namespace std;
using namespace pugi;
using namespace gnut;

namespace great {

	// Constructor
	// ----------
	t_gcfg_sp3orb::t_gcfg_sp3orb()
		: t_gsetbase(),
		t_gsetgen(),
		t_gsetinp(),
		t_gsetrec(),
		t_gsetout(),
		t_gsetgnss(),
		t_gsetsp3orb(),
		t_gsetforce()
	{
		
		_IFMT_supported.insert(SP3_INP);
		_IFMT_supported.insert(POLEUT1_INP);

		_OFMT_supported.insert(ORB_OUT);
		_OFMT_supported.insert(ICS_OUT);
		_OFMT_supported.insert(ICSLEO_OUT);
	}


	// Destructor
	// ----------
	t_gcfg_sp3orb::~t_gcfg_sp3orb()
	{}


	// settings check
	// ----------
	void t_gcfg_sp3orb::check()
	{
		t_gsetgen::check();
		t_gsetinp::check();
		t_gsetout::check();
		t_gsetrec::check();
		t_gsetgnss::check();
		t_gsetforce::check();
		t_gsetsp3orb::check();
	}

	// settings help
	// ----------
	void t_gcfg_sp3orb::help()
	{
		t_gsetbase::help_header();
		t_gsetgen::help();
		t_gsetinp::help();
		t_gsetout::help();
		t_gsetrec::help();
		t_gsetgnss::help();
		t_gsetforce::help();
		t_gsetsp3orb::help();
	}

} // namespace
