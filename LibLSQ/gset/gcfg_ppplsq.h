
#ifndef GCFG_PPPLSQ_H
#define GCFG_PPPLSQ_H

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

#include "gset/gsetgen.h"
#include "gset/gsetinp.h"
#include "gset/gsetout.h"
#include "gset/gsetproc.h"
#include "gset/gsetgnss.h"
#include "gset/gsetpar.h"
#include "gset/gsetrec.h"
#include "gset/gsetamb.h"
#include "gset/gsetturboedit.h"
#include "gset/gsetflt.h"

#include "gall/gallproc.h"
#include "gall/gallprec.h"
#include "gmodels/glsqproc.h"

#include "gproc/gpreproc.h"
#include "gproc/gpreproc.h"
#include "gproc/gppplsq.h"

#include "gproc/gexeturboedit.h"

#include "gcoders/gcoder.h"
#include "gcoders/rinexo.h"
#include "gcoders/rinexo_leo.h"
//#include "gcoders/rinexc.h"
#include "gcoders/rinexc_leo.h"
#include "gcoders/rinexn.h"
#include "gcoders/biasinex.h"
#include "gcoders/biabernese.h"
#include "gcoders/sp3.h"
#include "gcoders/sp3leo.h"
#include "gcoders/atx.h"
#include "gcoders/blq.h"
#include "gcoders/upd.h"
#include "gcoders/aug.h"
#include "gcoders/ambupd.h"

#include "gcoders/dvpteph405.h"
#include "gcoders/poleut1.h"
#include "gio/gfile.h"
#include "gdata/gleapsecond.h"
//#include "gdata/gupd_epoch.h"
//#include "gcoders/upd_epoch.h"
#include "gdata/gifcb.h"
#include "gcoders/ifcb.h"

using namespace std;
using namespace pugi;

namespace great
{

	class LSQ_LIBRARY_EXPORT t_gcfg_ppplsq : public t_gsetgen,
		public t_gsetinp,
		public t_gsetout,
		public t_gsetgnss,
		public t_gsetproc,
		public t_gsetrec,
		public t_gsetamb,
		public t_gsetpar,
		public t_gsetflt,
		public t_gsetturboedit
	{

	public:
		t_gcfg_ppplsq();
		~t_gcfg_ppplsq();

		void check();                                 // settings check
		void help();                                  // settings help

	protected:

	private:

	};

} // namespace

#endif