
#ifndef GBANCROFT_H
#define GBANCROFT_H
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  Purpose: statistical function (1D)
  Version: $ Rev: $

  2014-04-18 /PV: created
  2018-09-28 /JD: revised

-*/

#include "newmat/newmatap.h"

#include "gall/gallobs.h"
#include "gall/gallnav.h"
#include "gutils/gtriple.h"
#include "gutils/gmutex.h"
#include "gset/gsetbase.h"

namespace gnut {   

	BASE_LIBRARY_EXPORT int gbancroft(const Matrix& BBpass, ColumnVector& pos);

	BASE_LIBRARY_EXPORT inline double lorentz(const ColumnVector& aa, const ColumnVector& bb);

} // namespace

#endif

