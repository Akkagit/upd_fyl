/**
*
* @verbatim
History
-1.0 glFeng   2019-02-28 creat the file.
-1.1 BoWong   2019-04-10 Adding Doxygen Style Code Remarks
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		   gcfg_upd.cpp
* @brief	   control set from XML
*
* @author      glFeng, Wuhan University
* @version	   1.0.0
* @date		   2019-04-10
*
*/

#include "gset/gcfg_upd.h"

using namespace std;
using namespace pugi;
using namespace gnut;
namespace great
{
	/** @brief default constructor. */
	t_gcfg_upd::t_gcfg_upd()
		: t_gsetbase(),
		t_gsetgen(),
		t_gsetinp(),
		t_gsetout(),
		t_gsetgnss(),
		t_gsetamb(),
		t_gsetproc()
	{
		_IFMT_supported.insert(AMBFLAG12_INP);
		_IFMT_supported.insert(AMBFLAG13_INP);
		_IFMT_supported.insert(AMBFLAG14_INP);
		_IFMT_supported.insert(AMBFLAG15_INP);
		_IFMT_supported.insert(AMBUPD_INP);
		_IFMT_supported.insert(UPD_INP);
		_IFMT_supported.insert(RINEXO_INP);
		_IFMT_supported.insert(BIABERN_INP);
		_IFMT_supported.insert(RINEXN_INP);
		_IFMT_supported.insert(IFCB_INP);

		_OFMT_supported.insert(UPD_OUT);
		_OFMT_supported.insert(LOG_OUT);
		_OFMT_supported.insert(AMBFLAG_DIR_OUT);
	}

	/** @brief default destructor. */
	t_gcfg_upd::~t_gcfg_upd()
	{

	}

	/** @brief settings check. */
	void t_gcfg_upd::check()
	{
		t_gsetgen::check();
		t_gsetinp::check();
		t_gsetout::check();
		t_gsetgnss::check();
		t_gsetamb::check();
	}

	/** @brief settings help. */
	void t_gcfg_upd::help()
	{
		t_gsetbase::help_header();
		t_gsetgen::help();
		t_gsetinp::help();
		t_gsetout::help();
		t_gsetgnss::help();
		t_gsetamb::help();
		t_gsetbase::help_footer();
	}



}