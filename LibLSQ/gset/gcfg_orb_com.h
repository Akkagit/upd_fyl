
#ifndef GCFG_ORB_H
#define GCFG_ORB_H

/* ----------------------------------------------------------------------
(c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

(c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
Research Institute of Geodesy, Topography and Cartography
Ondrejov 244, 251 65, Czech Republic

Purpose: implements settings for testing Tefnut PP lite
Version: $ Rev: $

2012-10-23 /JD: created

-*/

#include <string>
#include <iostream>
#include <signal.h>

#include "gutils/ExportLSQ.h"
#include "gall/gallprec.h"
#include "gall/gallpcv.h"
#include "gall/gallobs.h"
#include "gall/gallotl.h"
#include "gall/gallbias.h"
#include "gproc/gpppflt.h"
#include "gproc/gpreproc.h"
#include "gio/gio.h"
#include "gio/glog.h"
#include "gio/gfile.h"
#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"
#include "gdata/gdata.h"
#include "gcoders/gcoder.h"
#include "gcoders/rinexo.h"
#include "gcoders/rinexc.h"
#include "gcoders/rinexn.h"
#include "gcoders/biasinex.h"
#include "gcoders/biabernese.h"
#include "gcoders/sp3.h"
#include "gcoders/atx.h"
#include "gcoders/blq.h"
#include "gmodels/gpcv.h"
#include "gmodels/gbancroft.h"

#include "gset/gsetgen.h"
#include "gset/gsetinp.h"
#include "gset/gsetout.h"
#include "gset/gsetproc.h"
#include "gset/gsetgnss.h"
#include "gset/gsetflt.h"
#include "gset/gsetrec.h"
#include "gset/gsetleo.h"
#include "gset/gsetpar.h"
#include "gset/gsetamb.h"

using namespace std;
using namespace pugi;
using namespace gnut;

namespace great {

	class LSQ_LIBRARY_EXPORT t_gcfg_orb_com : public t_gsetgen,
		public t_gsetinp,
		public t_gsetout,
		public t_gsetgnss,
		public t_gsetproc,
		public t_gsetrec,
		public t_gsetleo,
		public t_gsetpar,
		public t_gsetamb
	{

	public:
		t_gcfg_orb_com();
		~t_gcfg_orb_com();

		void check();                                 // settings check
		void help();                                  // settings help

		bool check_inp();

	protected:

	private:

	};

} // namespace

#endif