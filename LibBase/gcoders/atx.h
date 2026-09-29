
#ifndef ATX_H
#define ATX_H
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  Purpose: Antenna Exchange format (ATX)
  Version: $Rev:$

  2011-04-20 /JD: created

-*/

#include <vector> 

#include "gcoders/gcoder.h"

using namespace std;

namespace gnut {

class BASE_LIBRARY_EXPORT t_atx : public t_gcoder {

 public:
   t_atx( t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE );
  virtual ~t_atx(){};

  virtual  int decode_head(char* buff, int bufLen,           vector<string>& errmsg);
  virtual  int decode_data(char* buff, int bufLen, int& cnt, vector<string>& errmsg);

 protected:

 private:

};

} // namespace

#endif
