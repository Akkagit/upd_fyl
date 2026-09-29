
#ifndef GSYNC_H
#define GSYNC_H 

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: sync/filter epochs
  Version: $ Rev: $

  2015-12-17 /JD: created

-*/

#include <cmath>

#include "gio/glog.h"
#include "gutils/gtime.h"
#include "gutils/ExportBase.h"

using namespace std;

namespace gnut
{
	BASE_LIBRARY_EXPORT bool time_sync(const t_gtime& epo, double smp, double scl, t_glog* glog);
} // namespace

#endif // # GSYNC_H
