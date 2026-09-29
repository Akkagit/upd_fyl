
#ifndef GEOP_H
#define GEOP_H
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
 
  Purpose: earth orientation parameters, base class
           Neither precession nor nutation apllied -> return identity matrix
  Version: $ Rev: $

  2013-03-27 /PV: created
  
-*/

#include "gdata/gdata.h"
#include "gutils/gconst.h"
#include "gutils/gtime.h"
#include "gutils/gtriple.h"

using namespace std;

namespace gnut {   

// ----------
class t_geop {

 public:
  t_geop();
  virtual ~t_geop();
 
  Matrix nutMatrix(double mjd);
  Matrix precMatrix(double mjd);
   
// protected:
   double _normangle(double x);   
};

} // namespace

#endif
