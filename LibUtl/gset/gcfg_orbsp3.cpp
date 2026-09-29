/**
*
* @verbatim
    History
     -1.0 Yujie Qin  2019-12-02 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		       gcfg_orbsp3.cpp
* @brief	       orb to sp3 set from XML
*
* @author          Yujie Qin, Wuhan University
* @version	       1.0.0
* @date		       2019-12-02
*
*/

#include <iomanip>
#include <sstream>

#include "gset/gcfg_orbsp3.h"

using namespace std;
using namespace pugi;
using namespace gnut;

namespace great {

	// Constructor
	// ----------
	t_gcfg_orbsp3::t_gcfg_orbsp3()
		: t_gsetbase(),
		t_gsetgen(),
		t_gsetinp(),
		t_gsetrec(),
		t_gsetout(),
		t_gsetgnss(),
		t_gsetorbsp3()
	{

		_IFMT_supported.insert(ORB_INP);	
		_IFMT_supported.insert(POLEUT1_INP);

		_OFMT_supported.insert(SP3_OUT);
	}


	// Destructor
	// ----------
	t_gcfg_orbsp3::~t_gcfg_orbsp3()
	{}


	// settings check
	// ----------
	void t_gcfg_orbsp3::check()
	{
		t_gsetgen::check();
		t_gsetinp::check();
		t_gsetout::check();
		t_gsetrec::check();
		t_gsetgnss::check();
		t_gsetorbsp3::check();
	}

	// settings help
	// ----------
	void t_gcfg_orbsp3::help()
	{
		t_gsetbase::help_header();
		t_gsetgen::help();
		t_gsetinp::help();
		t_gsetout::help();
		t_gsetrec::help();
		t_gsetgnss::help();
		t_gsetorbsp3::help();
	}

} // namespace
