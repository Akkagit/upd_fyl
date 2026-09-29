
#ifndef GTIDE96_H
#define GTIDE96_H
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  Purpose: implements tides
  Version: $ Rev: $

  2012-11-05 /JD: created
 
 References:
    [1] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
-*/

#include <vector>

#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

#include "gutils/gconst.h"
#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gmodels/gephplan.h"
#include "gmodels/gtide.h"

using namespace std;

namespace gnut {   

// ----------
class BASE_LIBRARY_EXPORT t_gtide96 : public t_gtide {

 public:
  t_gtide96(t_glog* l);
  virtual ~t_gtide96();

  virtual t_gtriple tide_searth(const t_gtime& epo, t_gtriple& xyz);        // solid earth tides
  virtual t_gtriple tide_pole();                                            // pole tides
  virtual t_gtriple load_ocean(const t_gtime& epoch, const string& site, const t_gtriple& xRec);   // ocean tide loading
  virtual t_gtriple load_atmosph();                                         // atmospheric tide loading

 protected:

};

} // namespace

#endif
