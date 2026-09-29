/* ----------------------------------------------------------------------
* G-Nut - GNSS software development library
*
(c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

(c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
Research Institute of Geodesy, Topography and Cartography
Ondrejov 244, 251 65, Czech Republic

This file is part of the G-Nut C++ library.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 3 of
the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses>.

-*/

#include <iomanip>
#include <sstream>

#include "gset/gcfg_orb_com.h"

using namespace std;
using namespace pugi;
using namespace gnut;
namespace great {

	// Constructor
	// ----------
	t_gcfg_orb_com::t_gcfg_orb_com()
		: t_gsetbase(),
		t_gsetgen(),
		t_gsetinp(),
		t_gsetout(),
		t_gsetgnss(),
		t_gsetproc(),
		t_gsetrec(),
		t_gsetleo(),
		t_gsetpar(),
		t_gsetamb()
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
		_IFMT_supported.insert(ICSLEO_INP);
		_IFMT_supported.insert(OCEANTIDE_INP);
		_IFMT_supported.insert(POLEUT1_INP);
		_IFMT_supported.insert(SATPARS_INP);
		_IFMT_supported.insert(LEAPSECOND_INP);
		_IFMT_supported.insert(DE_INP);
		_IFMT_supported.insert(ATTITUDE_INP);
		_IFMT_supported.insert(AMBCON_INP);

		_OFMT_supported.insert(LOG_OUT);
		_OFMT_supported.insert(PPP_OUT);
		_OFMT_supported.insert(FLT_OUT);
		_OFMT_supported.insert(RECOVER_OUT);
		_OFMT_supported.insert(RECCLK_OUT);
		_OFMT_supported.insert(SATCLK_OUT);
		_OFMT_supported.insert(ICS_OUT);
		_OFMT_supported.insert(ICSLEO_OUT);
	}


	// Destructor
	// ----------
	t_gcfg_orb_com::~t_gcfg_orb_com()
	{}


	// settings check
	// ----------
	void t_gcfg_orb_com::check()
	{
		t_gsetgen::check();
		t_gsetinp::check();
		t_gsetout::check();
		t_gsetrec::check();
		t_gsetproc::check();
		t_gsetgnss::check();
		t_gsetleo::check();
	}

	// settings help
	// ----------
	void t_gcfg_orb_com::help()
	{
		t_gsetbase::help_header();
		t_gsetgen::help();
		t_gsetinp::help();
		t_gsetout::help();
		t_gsetrec::help();
		t_gsetproc::help();
		t_gsetgnss::help();
		t_gsetleo::help();
		t_gsetbase::help_footer();
	}

	bool t_gcfg_orb_com::check_inp()
	{
		string out = "check your xml file with input node, there is no : ";


		// CHECK INPUTS
		bool valid = true;

		auto tmp = this->inputs_all();
		if (tmp.empty())
		{
			cerr << out + "file" << endl;
			valid = false;
		}

		if (this->input_size("orb") == 0)
		{
			cerr << out + "orb" << endl;
			valid = false;
		}

		if (this->input_size("rinexc") == 0 && this->input_size("rinexn") == 0)
		{
			cerr << out + "rinexc or rinexn" << endl;
			valid = false;
		}

		if (this->input_size("rinexo") == 0)
		{
			cerr << out + "rinexo" << endl;
			valid = false;
		}

		if (this->input_size("ics") == 0)
		{
			cerr << out + "ics" << endl;
			valid = false;
		}

		if (this->input_size("poleut1") == 0)
		{
			cerr << out + "poleut1" << endl;
			valid = false;
		}

		if (this->input_size("satpars") == 0)
		{
			cerr << out + "satpars" << endl;
			valid = false;
		}

		if (this->input_size("icsleo") == 0)
		{
			cerr << out + "icsleo" << endl;
			valid = false;
		}

		if (this->input_size("atx") == 0)
		{
			cerr << out + "atx" << endl;
			valid = false;
		}

		/*if (this->input_size("attitude") == 0)
		{
			cerr << out + "attitude" << endl;
			valid = false;
		}*/

		return valid;
	}

} // namespace
