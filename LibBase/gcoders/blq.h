
#ifndef BLQ_H
#define BLQ_H
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  Purpose: 
  Version: $Rev:$

  2013-03-29 /PV: created

-*/

#include <vector> 

#include "newmat/newmat.h"

#include "gcoders/gcoder.h"
#include "gmodels/gotl.h"
#include "gall/gallotl.h"
#include "gutils/gtime.h"

using namespace std;

namespace gnut {

class BASE_LIBRARY_EXPORT t_blq : public t_gcoder {

 public:
   t_blq( t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE );
 virtual ~t_blq(){};

  virtual  int decode_head(char* buff, int sz,           vector<string>& errmsg);
  virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

 protected:
   double _lat;
   double _lon;
   string _site;
   Matrix _blqdata;

 private:

};

} // namespace

#endif
