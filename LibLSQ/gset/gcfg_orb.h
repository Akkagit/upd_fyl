/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2019-04-02  created
	 -1.1 ZhengHJ  2019-04-18  add Doxygen format comment
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gcfg_orb.h
* @brief		main setting class for orblsq application
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-18
*
*/


#ifndef GCFG_ORB_H
#define GCFG_ORB_H

#include <string>
#include <iostream>
#include <signal.h>

#include "gutils/ExportLSQ.h"
#include "gall/gallprec.h"
#include "gall/gallpcv.h"
#include "gall/gallobs.h"
#include "gall/gallotl.h"
#include "gall/gallbias.h"

#include "gproc/gpreproc.h"
#include "gproc/gpodlsqIF.h"
#include "gproc/gpodlsqALL.h"

#include "gmodels/glsqproc.h"
#include "gio/gio.h"
#include "gio/glog.h"
#include "gio/gfile.h"
#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"

#include "gdata/gdata.h"
#include "gdata/gnavde.h"
#include "gcoders/gcoder.h"
#include "gcoders/rinexo.h"
#include "gcoders/rinexc.h"
#include "gcoders/rinexn.h"
#include "gcoders/biasinex.h"
#include "gcoders/biabernese.h"
#include "gcoders/sp3.h"
#include "gcoders/atx.h"
#include "gcoders/blq.h"
#include "gcoders/icsnav.h"
#include "gcoders/orb.h"
#include "gcoders/poleut1.h"
#include "gcoders/satparameters.h"
#include "gcoders/dvpteph405.h"
#include "gmodels/gpcv.h"
#include "gmodels/gbancroft.h"

#include "gset/gsetgen.h"
#include "gset/gsetinp.h"
#include "gset/gsetout.h"
#include "gset/gsetproc.h"
#include "gset/gsetgnss.h"
#include "gset/gsetflt.h"
#include "gset/gsetrec.h"
#include "gset/gsetpar.h"

#include "gdata/gpoleut1.h"
#include "gall/gallsatparam.h"
#include "gall/gallproc.h"
#include "gdata/gleapsecond.h"

using namespace std;
using namespace pugi;

namespace great
{

	class LSQ_LIBRARY_EXPORT t_gcfg_orb : public t_gsetgen,
		public t_gsetinp,
		public t_gsetout,
		public t_gsetgnss,
		public t_gsetproc,
		public t_gsetrec,
		public t_gsetpar,
		public t_gsetamb
	{

	public:
		t_gcfg_orb();
		~t_gcfg_orb();

		void check();                                 // settings check
		void help();                                  // settings help

	protected:

	private:

	};

} // namespace

#endif