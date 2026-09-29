
#ifndef  GOCEANTIDE_H
#define  GOCEANTIDE_H

/* --------------------------------------------------------

-*/

#include "gdata/gdata.h"

using namespace std;

namespace gnut {
class t_gfocean_tide : public t_gdata {
	
  public:
    t_gfocean_tide();
	virtual ~t_gfocean_tide();
	
	int m();
	int n();
	int* dindex();
	double* coeff();
	void setdata(const int& m, const int& n, const int* dindex(), const double* coeff());
	
  private:
    int n;            // the index of spheric harmoics
    int m;            // the index of spheric harmoics
    int index[6];    // the index of Doodson
    double coeff[4];  // coefficients of the ocean tide(C_PROG, C_RETG, S_PROG, S_RETG)
};
} //namespace

#endif