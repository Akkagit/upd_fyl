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

#include "gcfg_editres.h"

using namespace std;
using namespace pugi;

namespace great
{
	t_gcfg_editres::t_gcfg_editres():
		t_gsetbase(),
		t_gsetinp(),
		t_gsetres()
	{
		_IFMT_supported.insert(AMBFLAG12_INP);
		_IFMT_supported.insert(AMBFLAG13_INP);
		_IFMT_supported.insert(RECOVER_INP);
		_IFMT_supported.insert(PREOBS_INP);
		_IFMT_supported.insert(PANDA_RECOVER_INP);
	}

	t_gcfg_editres::~t_gcfg_editres()
	{
	}

	void t_gcfg_editres::check()
	{
		t_gsetinp::check();
		t_gsetres::check();
	}

	void t_gcfg_editres::help()
	{
		t_gsetbase::help_header();
		t_gsetinp::help();
		t_gsetres::help();
		t_gsetbase::help_footer();
	}

}