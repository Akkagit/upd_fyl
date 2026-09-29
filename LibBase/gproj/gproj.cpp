
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  This file is part of the G-Nut C++ library.
 
-*/

#include <math.h>
#include <iostream>

#include "gproj/gproj.h"
#include "gutils/gcommon.h"

using namespace std;

namespace gnut {


// Conversion from LatLon [deg] sphere to generic project projection, keeping the HEIGHT!
// ---------- 
t_gtriple t_gproj::ll2proj(const t_gtriple& latlon, int dec )
{
  gtrace("t_gproj::ll2proj(t_gtriple)");

  t_gpair ll( latlon.gpair() );
  t_gpair xx = ll2proj(ll, dec);
   
  t_gtriple xxh(xx[0], xx[1], latlon[2]);

#ifdef DEBUG
  cout << " XXH: " << xx[0] << " " << xx[1]
       << " LLH: " << ll[0] << " " << ll[1] << endl;
#endif

  return xxh;
}

// Conversion from generic project projection to LatLon [deg] sphere, keeping the HEIGHT!
// ---------- 
t_gtriple t_gproj::proj2ll(const t_gtriple& xy_idx)
{
  gtrace("t_gproj::proj2ll(t_gtriple)");

  t_gpair xy( xy_idx.gpair() );
  t_gpair xx = proj2ll(xy);

  t_gtriple llh(xx[0], xx[1], xy_idx[2]);

#ifdef DEBUG
  cout << " XX: " <<  xx[0] << " " <<  xx[1]
       << " LL: " << llh[0] << " " << llh[1] << " " << llh[2] << endl;
#endif

  return llh;
}

   
} // namespace
