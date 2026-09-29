/**
*
* @verbatim
History
-1.0 Hongmin Zhang  2019-10-24 creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		       gcfg_slr.h
* @brief	       satellite laser ranging
*
* @author          Hongmin Zhang, Wuhan University
* @version	       1.0.0
* @date		       2019-10-24
*
*/

#include <iomanip>
#include <sstream>

#include "gset/gcfg_slr.h"

using namespace std;
using namespace pugi;
using namespace gnut;

namespace great {

	// Constructor
	// ----------
	t_gcfg_slr::t_gcfg_slr()
		: t_gsetbase(),
		t_gsetgen(),
		t_gsetinp(),
		t_gsetrec(),
		t_gsetout(),
		t_gsetgnss(),
		t_gsetsp3orb(),
		t_gsetforce()
	{
		_IFMT_supported.insert(RINEXO_INP);
		_IFMT_supported.insert(RINEXC_INP);
		_IFMT_supported.insert(RINEXN_INP);
		_IFMT_supported.insert(ATX_INP);
		_IFMT_supported.insert(BLQ_INP);
		_IFMT_supported.insert(SP3_INP);
		_IFMT_supported.insert(BIABERN_INP);
		_IFMT_supported.insert(SINEX_INP);
		_IFMT_supported.insert(ORB_INP);
		_IFMT_supported.insert(ICS_INP);
		_IFMT_supported.insert(OCEANTIDE_INP);
		_IFMT_supported.insert(POLEUT1_INP);
		_IFMT_supported.insert(SATPARS_INP);
		_IFMT_supported.insert(LEAPSECOND_INP);
		_IFMT_supported.insert(DE_INP);
		_IFMT_supported.insert(LRA_INP);
		_IFMT_supported.insert(NPT_INP);
		_IFMT_supported.insert(SLRF_INP);
		_IFMT_supported.insert(ATTITUDE_INP);
		_IFMT_supported.insert(BLQ_INP);
		_IFMT_supported.insert(RB_INP);
		_IFMT_supported.insert(STA_INP);


		_OFMT_supported.insert(LOG_OUT);
		_OFMT_supported.insert(PPP_OUT);
		_OFMT_supported.insert(RECCLK_OUT);
		_OFMT_supported.insert(SATCLK_OUT);
		_OFMT_supported.insert(PRE_OUT);
		_OFMT_supported.insert(RESULT_OUT);
	}


	// Destructor
	// ----------
	t_gcfg_slr::~t_gcfg_slr()
	{}


	// settings check
	// ----------
	void t_gcfg_slr::check()
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
	void t_gcfg_slr::help()
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
