
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  This file is part of the G-Nut C++ library.
 
-*/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "gmodels/gfunc.h"

using namespace std;

namespace gnut {   

// constructor
// ----------
t_gfunc::t_gfunc( const map<string,double>& coefficients )
{
  _coeff = coefficients;   

#ifdef DEBUG
  map<string,double>::iterator it;
  for( it = _coeff.begin(); it != _coeff.end(); ++it )
     cout << "COEFF:  " << fixed   << setprecision(3)
          << " " << setw(3) << it->first << " = " << setw(9) << it->second
          << endl;
#endif
}
   
// destructor
// ----------
t_gfunc::~t_gfunc()
{}
   
} // namespace
