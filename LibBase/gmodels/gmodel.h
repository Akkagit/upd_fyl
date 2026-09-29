
#ifndef GMODEL_H
#define GMODEL_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
 
  Purpose: base abstract class for models
  Version: $ Rev: $

  2014-11-27 /PV: created

  2019-3-25  /ZHJ:add new data and compute
-*/

#include <string>
#include <map>
#include <cmath>

#include "gio/glog.h"
#include "gall/gallpar.h"
#include "gutils/gtime.h"
#include "gdata/gsatdata.h"
#include "gset/gsetbase.h"
#include "gmodels/gtropo.h"
#include "gall/gallbias.h"

using namespace std;

namespace gnut {   

class BASE_LIBRARY_EXPORT t_gmodel
{
 public:
  virtual ~t_gmodel();

   void setSite(string site);

   void   glog(t_glog* l){ _log = l; }         
   
   virtual double cmpObs(t_gtime& epoch, t_gallpar& param, t_gsatdata&, t_gobs& gobs, bool = false) = 0;
   virtual double cmpObsD(t_gtime& epoch, t_gallpar& param, t_gsatdata& gsatdata, t_gobs& gobs) = 0;
   virtual double windUp(t_gtime& epoch, const string, const ColumnVector&, const ColumnVector&){return 0.0;};

   virtual double tropoDelay(t_gtime& epoch, t_gallpar& param, t_gtriple ell, t_gsatdata& satdata){return 0.0;};   

   virtual double getZHD(const string& site , const t_gtime& epo) { return 0.0; };
   virtual double getZWD(const string& site , const t_gtime& epo) { return 0.0; };

   virtual void reset_observ(OBSCOMBIN obs) { };
   virtual void reset_tropmf(ZTDMPFUNC mf) {};	     //lvhb added for npp, 20200918
   virtual void reset_PCO(bool cal = true) {};		 //lvhb added for npp, 20200918

   virtual int outlierDetect(      vector<t_gsatdata>& data, SymmetricMatrix& Qx, const SymmetricMatrix&, const ColumnVector&) = 0;
   virtual int outlierDetect(vector<t_gsatdata>& data, SymmetricMatrix& Qx, const SymmetricMatrix&, vector<t_gsatdata>::iterator &) = 0;
   virtual int outlierDetect(      vector<t_gsatdata>& data, SymmetricMatrix& Qx, const SymmetricMatrix&) = 0;
   virtual int outlierDetect_chi(  vector<t_gsatdata>& data, SymmetricMatrix& Qx, const SymmetricMatrix&, const ColumnVector&) = 0;
   
   shared_ptr<t_gtropo> tropoModel(){ return _tropoModel; }
   
   void  setBIAS(t_gallbias* bia) {_gallbias = bia;}

 protected:
   shared_ptr<t_gtropo> _tropoModel;
   t_gsetbase*         _settings;
   string              _site;
   t_glog*             _log;
   bool                _phase;
   t_gallbias*         _gallbias;

};

} // namespace

#endif //  GMODEL_H
