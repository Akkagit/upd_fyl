
#ifndef SP3_H
#define SP3_H
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  Purpose: 
  Version: $Rev:$

  2011-04-20 /JD: created

-*/

#include <vector> 

#include "gcoders/gcoder.h"
#include "gutils/gtime.h"
#include "gall/gallprec.h"

using namespace std;

namespace gnut {

class BASE_LIBRARY_EXPORT t_sp3 : public t_gcoder {

 public:
   t_sp3( t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE );
  virtual ~t_sp3(){};

  virtual  int decode_head(char* buff, int sz,           vector<string>& errmsg);
  virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);
  virtual  int encode_head(char* buff, int sz,           vector<string>& errmsg);
  virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

 protected:

 private:
  t_gtime        _start;
  t_gtime        _lastepo;
  long           _orbintv; // [sec]
  int            _nepochs;  
  int            _nrecmax;
  int            _nrecord;
  string         _orbrefs;
  string         _orbtype;
  string         _agency;
  int            _num_leo = 0;
  map<string,string> _leo;
  vector<string> _prn;  
  vector<string> _sat;
  vector<int>    _acc;    
  vector<string> _timesys;
  vector<int>    _accbase;
  int            _maxsats;
  string         _sattype;
  string         _data_type;
  t_map_prn      _mapsp3;
  t_map_dat      _sp3data;
};

} // namespace

#endif
