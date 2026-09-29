
#ifndef STOCHASTIC_H
#define STOCHASTIC_H 
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: stochastic models
  Version: $ Rev: $

  2012-05-02 /PV: created

-*/

#include "gutils/gtime.h"
#include "newmat/newmat.h"

using namespace std;

namespace gnut {   

class BASE_LIBRARY_EXPORT t_stochastic
{
 public:
  t_stochastic();
  virtual ~t_stochastic() {};
  
  virtual double getQ()
  {
    return 0.0;
  }
  

 protected:
  
 private:

};

class BASE_LIBRARY_EXPORT t_randomwalk : public t_stochastic
{
 public:
  t_randomwalk();
  virtual ~t_randomwalk() {};
   
  virtual double getQ();
  void setTprev(const t_gtime&);
  void setTcurr(const t_gtime&);
  void updateTime(const t_gtime&);
  void setq(double q);
  double get_dt();
  
 protected:
 private:
  t_gtime _Tprev;
  t_gtime _Tcurr;
  double _dSig;
};

class BASE_LIBRARY_EXPORT t_whitenoise : public t_stochastic
{
 public:
  t_whitenoise(double);
  virtual ~t_whitenoise() {};
   
   virtual double getQ();
   void setVar(double);
   
 private:
   double _var;
};

class BASE_LIBRARY_EXPORT t_statemode
{
public:
	t_statemode();
	t_statemode(int order,double dt,double noise);
	virtual ~t_statemode();

	int order;
	// x(t+1) = M*x(t)  P
	Matrix M;
	SymmetricMatrix P;

private:
	static const double _coeff[6];
	
	
};

} // namespace

#endif // STOCHASTIC_H
