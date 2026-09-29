
#ifndef GCFG_PPP_H
#define GCFG_PPP_H

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

#include "gio/gfile.h"
#include "gio/gtcp.h"
#include "gio/gserial.h"

#include "gset/gsetgen.h"
#include "gset/gsetinp.h"
#include "gset/gsetout.h"
#include "gset/gsetproc.h"
#include "gset/gsetgnss.h"
#include "gset/gsetflt.h"
#include "gset/gsetrec.h"
#include "gset/gsetpar.h"//lvhb 202007
#include "gset/gsetamb.h"
#include "gall/gallprec.h"

#include "gproc/gpreproc.h"
#include "gproc/gexeturboedit.h"

#include "gcoders/bncobs.h"
#include "gcoders/bnccorr.h"
#include "gcoders/gcoder.h"
#include "gcoders/rinexo.h"
#include "gcoders/rinexc.h"
#include "gcoders/rinexn.h"
#include "gcoders/biasinex.h"
#include "gcoders/biabernese.h"
#include "gcoders/sp3.h"
#include "gcoders/atx.h"
#include "gcoders/blq.h"
#include "gcoders/upd.h"
#include "gcoders/aug.h"//lvhb 202007
#include "gcoders/ionex.h" // Fgl 202007

#include "gcoders/dvpteph405.h"
#include "gcoders/poleut1.h"
#include "gdata/gleapsecond.h"
#include "gdata/gifcb.h"
#include "gdata/gnavde.h"
#include "gcoders/ifcb.h"

#include "gutils/ExportINS.h"

using namespace std;
using namespace pugi;
using namespace great;

namespace gsins
{

	class INS_LIBRARY_EXPORT t_gcfg_ppp : virtual public t_gsetgen,
		virtual public t_gsetinp,
		virtual public t_gsetout,
		virtual public t_gsetgnss,
		virtual public t_gsetproc,
		virtual public t_gsetflt,
		virtual public t_gsetrec,
		virtual public t_gsetamb,
		virtual public t_gsetturboedit,
		virtual public t_gsetpar
	{

	public:
		t_gcfg_ppp();
		~t_gcfg_ppp();

		void check();                                 // settings check
		void help();                                  // settings help

	protected:

	private:

	};

} // namespace

#endif