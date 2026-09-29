
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com) 
  This file is part of the G-Nut C++ library.
 
-*/

#include <cmath>
#include <iomanip>

#include "gmodels/gtropoblind.h"

using namespace std;

namespace gnut {   

// ----------------------------- BLIND MODELS ------------------------------------

// ---------
// MOPS ZHD interface
// ---------
double t_blindmops::getZHD(const t_gtriple& ell, const t_gtime& epo) // RAD
{
  gtrace("t_blindmops::getZHD");
      
  double zhd, zwd;
  double ele = 90.0;
  double doy = double( epo.doy() );
  double lat = ell[0]*R2D;
  double hel = ell[2];

  _model->mops(lat, hel, ele, doy, zhd, zwd);  // DEG

#ifdef DEBUG
  cout << " MOPS: " << fixed << setprecision(3)
       << "  DOY: " << setw(12) << doy
       << "  ELE: " << setw(12) << ele
       << "  ZHD: " << setw(12) << zhd
       << endl;           
#endif

  return zhd;
}

// ---------
// MOPS ZWD interface
// ---------
double t_blindmops::getZWD(const t_gtriple& ell, const t_gtime& epo) // RAD
{
  gtrace("t_blindmops::getZWD");
      
  double zhd, zwd;
  double ele = 90.0;
  double doy = double( epo.doy() );
  double lat = ell[0]*R2D;
  double hel = ell[2];

  _model->mops(lat, hel, ele, doy, zhd, zwd); // DEG

#ifdef DEBUG
  cout << " MOPS: " << fixed << setprecision(3)
       << "  DOY: " << setw(12) << doy
       << "  ELE: " << setw(12) << ele
       << "  ZWD: " << setw(12) << zwd
       << endl;           
#endif

  return zwd;
}

} // namespace
