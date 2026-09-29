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

#include "gset/gcfg_podlsq.h"

using namespace std;
using namespace pugi;

namespace great {

	// Constructor
	// ----------
	t_gcfg_orb_leo::t_gcfg_orb_leo()
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
		_IFMT_supported.insert(IFCB_INP);
		_IFMT_supported.insert(SINEX_INP);
		_IFMT_supported.insert(ORB_INP);
		_IFMT_supported.insert(ORBLEO_INP);
		_IFMT_supported.insert(ICS_INP);
		_IFMT_supported.insert(ICSLEO_INP);
		_IFMT_supported.insert(AMBCON_INP);
		_IFMT_supported.insert(AMBCONLEO_INP);
		_IFMT_supported.insert(OCEANTIDE_INP);
		_IFMT_supported.insert(POLEUT1_INP);
		_IFMT_supported.insert(SATPARS_INP);
		_IFMT_supported.insert(LEAPSECOND_INP);
		_IFMT_supported.insert(PREOBS_INP);
		_IFMT_supported.insert(SLRF_INP);
		_IFMT_supported.insert(LRA_INP);
		_IFMT_supported.insert(DE_INP);
		_IFMT_supported.insert(RB_INP);
		_IFMT_supported.insert(ATTITUDE_INP);
		_OFMT_supported.insert(LOG_OUT);
		_OFMT_supported.insert(PPP_OUT);
		_OFMT_supported.insert(FLT_OUT);
		_OFMT_supported.insert(RECOVER_OUT);
		_OFMT_supported.insert(RECCLK_OUT);
		_OFMT_supported.insert(SATCLK_OUT);
		_OFMT_supported.insert(ICS_OUT);
		_OFMT_supported.insert(ICSLEO_OUT);
		_OFMT_supported.insert(RB_OUT);
		_OFMT_supported.insert(STA_OUT);
	}


	// Destructor
	// ----------
	t_gcfg_orb_leo::~t_gcfg_orb_leo()
	{}


	// settings check
	// ----------
	void t_gcfg_orb_leo::check()
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
	void t_gcfg_orb_leo::help()
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

	bool t_gcfg_orb_leo::check_inp()
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
		if (this->input_size("atx") == 0)
		{
			cerr << out + "atx" << endl;
			valid = false;
		}
		return valid;
	}

	void t_gcfg_orb_leo::arg(int argc, char* argv[], bool add, bool thin)
	{
		string conf("");
		string save("");
		bool stdi = false; // ENFORCE STDINP
		bool thrd = false;

		//cerr << "command-line opt: " << argc << endl;  for( int i=1; i<argc; ++i ) cerr << argv[i] << " "; cerr << endl;

		if (argc == 1)
		{
			cerr << "This is a app used for pod of gns satellite.      " << endl
				<< "your can use  --help get how to use the cmd.      " << endl
				<< "         use  --X    get how to set the xml file. " << endl
				<< "         use  --x    to  use the xml file and cal." << endl
				<< "         use  --v    get the verb information.    " << endl;
			return;
		}

		for (int i = 1; i < argc; ++i)
		{
			if (!strcmp(argv[i], "--help") && i < argc) { usage(); exit(0); }
			else if (!strcmp(argv[i], "-h") && i < argc) { usage(); exit(0); }
			else if (!strcmp(argv[i], "-V") && i < argc) { cout << app() << endl; exit(0); }
			else if (!strcmp(argv[i], "-X") && i < argc) { help(); exit(0); }
			else if (!strcmp(argv[i], "-m") && i < argc) { thrd = true; }
			else if (!strcmp(argv[i], "-brdm") && i < argc) { isUseBrdm = true; }
			else if (!strcmp(argv[i], "-ambfix") && i < argc) { isAmbfix = true; }
			else if (!strcmp(argv[i], "-ion") && i < argc) { isUseIon = true; }
			else if (!strcmp(argv[i], "-x") && i + 1 < argc) { conf = argv[++i]; }
			else if (!strcmp(argv[i], "-Z") && i + 1 < argc) { save = argv[++i]; }
			else if (!strcmp(argv[i], "-v") && i + 1 < argc) { if (_log) { _log->verb(atoi(argv[++i])); } } // ZDE PRILIS BRZY!!
			else if (!strcmp(argv[i], "-l") && i + 1 < argc) { if (_log) { _log->mask(argv[++i]); } } // ZDE PRILIS BRZY!!
			//warning if no other commandline option is supported
			else if (!add) fprintf(stderr, "Unknown option: %s\n", argv[i]);
			else continue;
		}

		if (!conf.empty())
		{
			if (read(conf)) { exit(1); } // file XML
		}

		if (thin)
		{
			xml_node node = _doc.child(XMLKEY_ROOT);
			_default_attr(node, "thin", "true", true);
		}

		this->check();

		if (!conf.empty()) { if (_log) _log->comment(0, "gsetbase", "xconfig: read from file " + _name); }
		else { if (_log) _log->comment(0, "gsetbase", "xconfig: read from istream"); }

		if (!save.empty()) {
			if (_log) _log->comment(0, "gsetbase", "xconfig: save to file " + _name);
			write(save);
		}

		// COMMENTS
		if (_log && !thin) {
			_log->comment(0, "gsetbase", "program: " + app());
			_log->comment(0, "gsetbase", "==========");
			_log->comment(1, "gsetbase", "setting: " + conf);
			_log->comment(1, "gsetbase", "logfile: " + _log->mask());
			_log->comment(1, "gsetbase", "verbose: " + int2str(_log->verb()));
			_log->comment(1, "gsetbase", "==========");
		}

		return;

	}

} // namespace
