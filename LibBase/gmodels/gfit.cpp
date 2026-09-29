
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  This file is part of the G-Nut C++ library.
 
-*/

#include "gmodels/gfit.h"

using namespace std;

namespace gnut {   

// constructor
// ----------
t_gfit::t_gfit( t_gfunc* model, int max_iter )
{
  _log      = 0;
  _func     = model;
  _lim_par  =  1e-2;   // convergence criterium limit for parameter // 1e-2 generally ok!

   _max_iter = max_iter;
}
   
} // namespace
