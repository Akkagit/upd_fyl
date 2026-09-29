/**
*
* @verbatim
    History
     -1.0 Yujie Qin  2019-10-10 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		       gcfg_sp3orb.h
* @brief	       sp3 to orb set from XML
*
* @author          Yujie Qin, Wuhan University
* @version	       1.0.0
* @date		       2019-10-10
*
*/


#ifndef GCFG_ORB_H
#define GCFG_ORB_H

#include <string>
#include <iostream>
#include <signal.h>

#include "gutils/ExportUtl.h"
#include "gall/gallprec.h"
#include "gall/gallpcv.h"
#include "gall/gallobs.h"
#include "gall/gallotl.h"
#include "gall/gallbias.h"

#include "gproc/gpreproc.h"
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
#include "gset/gsetsp3orb.h"
#include "gset/gsetforce.h"
//#include "gset/gsetforce_leo.h"
#include "gset/gsetrec.h"
#include "gdata/gpoleut1.h"
#include "gall/gallsatparam.h"
#include "gall/gallproc.h"
#include "gdata/gleapsecond.h"

using namespace std;
using namespace pugi;

using namespace gnut;

namespace great
{

	class UTL_LIBRARY_EXPORT t_gcfg_sp3orb : public t_gsetgen,
		public t_gsetinp,
		public t_gsetout,
		public t_gsetrec,
		public t_gsetgnss,
		public t_gsetsp3orb,
		public t_gsetforce
	{

	public:
		t_gcfg_sp3orb();
		~t_gcfg_sp3orb();

		void check();                                 // settings check
		void help();                                  // settings help

	protected:

	private:

	};

} // namespace

#endif