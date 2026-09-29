
#ifndef GFUNC_H
#define GFUNC_H
 
/* ----------------------------------------------------------------------
  Purpose: implements different functions (and derivatives) for 2D data fitting
  Version: $ Rev: $

  2014-03-21 /JD: created

-*/

#include <map>
#include <string>

#include "newmat/newmat.h"
#include "newmat/newmatio.h"

using namespace std;

namespace gnut {   

class t_gfunc 
{   
 public:
            t_gfunc( const map<string, double>& coefficients );
   virtual ~t_gfunc();
      
   virtual int value( const ColumnVector& dat, const ColumnVector& fit, ColumnVector& val ){ return 0; };
   virtual int deriv( const ColumnVector& dat, const ColumnVector& fit, Matrix& val ){ return 0; };

 protected:
   map<string,double> _coeff;

};

} // namespace

#endif