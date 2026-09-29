/**
*
* @verbatim
	History
	 -1.0 liyuhao	 2020-07-28 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gcfg_clkdif.cpp
* @brief		clock differ control set from XML
*
* @author		liyuhao, Wuhan University
* @version		1.0.0
* @date			2020-07-28
*
*/

#include "gcfg_clkdif.h"

using namespace std;
using namespace pugi;

namespace great
{
	t_gcfg_clkdif::t_gcfg_clkdif() :
		t_gsetbase(),
		t_gsetinp()
	{
		_IFMT_supported.insert(RINEXC_PRD);
		_IFMT_supported.insert(RINEXC_REF);
	}

	t_gcfg_clkdif::~t_gcfg_clkdif()
	{
	}
	void t_gcfg_clkdif::check()
	{
		t_gsetinp::check();
	}

	void t_gcfg_clkdif::help()
	{
		t_gsetbase::help_header();
		t_gsetinp::help();
		t_gsetbase::help_footer();
	}
}