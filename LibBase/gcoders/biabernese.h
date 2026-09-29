
#ifndef BIABERNESE_H
#define BIABERNESE_H
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  Purpose: bias SINEX
  Version: $Rev:$

  2018-03-12 /PV: created
  2019-04-08 glfeng modified: add decode_data for BSX file from CAS & change the function combination

-*/ 

#include "gcoders/gcoder.h"
#include "gutils/gtime.h"
#include "gall/gallbias.h"

using namespace std;

namespace gnut {

class BASE_LIBRARY_EXPORT t_biabernese : public t_gcoder {

 public:
   t_biabernese( t_gsetbase* s, string version, int sz = DEFAULT_BUFFER_SIZE, string id = "biabernese" );
  virtual ~t_biabernese(){};

   virtual  int decode_head(char* buff, int bufLen,           vector<string>& errmsg);
   virtual  int decode_data(char* buff, int bufLen, int& cnt, vector<string>& errmsg);
   
 protected:
   int _decode_data_CODE(char* buff, int sz, int& cnt, vector<string>& errmsg);  // add by glfeng
   int _decode_data_sinex(char* buff, int sz, int& cnt, vector<string>& errmsg);  // add by glfeng
   int _decode_data_sinex_0(char* buff, int sz, int& cnt, vector<string>& errmsg);  // add by jqwu

   t_gallbias* _allbias;
   t_gtime _beg;     // add by glfeng
   t_gtime _end;     // add by glfeng

   string _ac;
   double _version; // 1.00 or 0.01
   bool _is_bias = false; // sinex bias file
   bool _is_absolute = false; // relative/absolute bias
};

} // namespace

#endif
