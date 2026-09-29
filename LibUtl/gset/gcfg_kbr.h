/**
*
* @verbatim
History
-1.0 Hongmin Zhang  2020-03-20 creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		       gcfg_kbr.h
* @brief	       KBR
*
* @author          Hongmin Zhang, Wuhan University
* @version	       1.0.0
* @date		       2020-03-20
*
*/


#ifndef GCFG_KBR_H
#define GCFG_KBR_H

#include <string>
#include <iostream>
#include <signal.h>

#include "gutils/ExportUtl.h"
#include "gall/gallkbr.h"
#include "gall/gallorb.h"

#include "gio/gio.h"
#include "gio/glog.h"
#include "gio/gfile.h"
#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"

#include "gdata/gdata.h"
#include "gcoders/orb.h"
#include "gcoders/kbr.h"

#include "gset/gsetgen.h"
#include "gset/gsetinp.h"
#include "gset/gsetout.h"

#include "gall/gallkbr.h"
#include "gall/gallorb.h"

using namespace std;
using namespace pugi;

using namespace gnut;

namespace great
{

	class UTL_LIBRARY_EXPORT t_gcfg_kbr : public t_gsetgen,
		public t_gsetinp,
		public t_gsetout,
		public t_gsetgnss
	{

	public:
		t_gcfg_kbr();
		~t_gcfg_kbr();

		void check();                                 // settings check
		void help();                                  // settings help

	protected:

	private:

	};

} // namespace

#endif