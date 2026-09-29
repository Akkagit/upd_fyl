
#ifndef BIASINEX_H
#define BIASINEX_H
 
/* ---------------------------------------------------------------------- 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  Purpose: bias SINEX
  Version: $Rev:$

  2018-03-09 /PV: created

-*/ 

#include "gcoders/gcoder.h"
#include "gcoders/sinex.h"
#include "gutils/gtime.h"
#include "gall/gallbias.h"

#define SNX_WINDOW_EXTENTION 6*3600    // [s] window for SNX temporal approximation

using namespace std;

namespace gnut {

class BASE_LIBRARY_EXPORT t_biasinex : public t_sinex {

 public:
   t_biasinex( t_gsetbase* s, string version, int sz = DEFAULT_BUFFER_SIZE, string id = "biasinex" );
  virtual ~t_biasinex(){};

 protected:
   virtual int _decode_comm();
   virtual int _decode_block();   

   virtual void _add_data(string id, t_gdata* pt_data);

   t_gallbias* _allbias;
};

} // namespace

#endif
