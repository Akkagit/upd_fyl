
#ifndef GEPHPLAN_H
#define GEPHPLAN_H
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
 
  Purpose: planetary ephemerises
  Version: $ Rev: $

  2013-03-26 /PV: created

-*/

#include <vector>

#include "gutils/gconst.h"
#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gmodels/geop80.h"

using namespace std;

namespace gnut {   

// ----------
class t_gephplan {

 public:
  t_gephplan();
  virtual ~t_gephplan();

 t_gtriple sunPos(double mjd, bool itrf = true);
 t_gtriple moonPos(double mjd);
 double    gmst(double mjd);
   
 protected:
   double _frac(double x);

   t_geop80 _eop;

private:
    map<double, t_gtriple> _record_itrf_sunpos;
   
};

} // namespace

#endif
