
#ifndef RINEXO_H
#define RINEXO_H
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: Clock RINEX encoder/decoder
  Version: $Rev:$

  2011-11-04 /JD: created
  2012-09-04 /JD: RINEX 3.01 implemented
  2012-09-24 /JD: RINEX 2.11 implemented
  2012-09-27 /JD: support of multi-GNSS (phase not converted to [m])
  2013-03-08 /JD: filtering via general settings
  2013-11-19 /JD: shared pointers
 
  Todo: header information completed 
        gobj getting information from rinexo header
        encoders implementation
        to finish obs-filtering via settings

-*/

#include <string> 
#include <vector> 

#include "gall/gallobs.h"
#include "gall/gallobj.h"
#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gutils/gsys.h"
#include "gutils/gobs.h"
#include "gutils/gtypeconv.h"
#include "gcoders/gcoder.h"
#include "gcoders/rinexo2.h"
#include "gcoders/rinexo3.h"
#include "gdata/grnxhdr.h"

using namespace std;

namespace gnut {

class BASE_LIBRARY_EXPORT t_rinexo : public t_rinexo3 {

 public:
   t_rinexo( t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE );
   t_rinexo(t_gtime beg, t_gtime end, t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);
  virtual ~t_rinexo(){};

  virtual  int decode_head(char* buff, int sz,           vector<string>& errmsg);
  virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);
  
  virtual  int encode_head(char* buff, int sz,           vector<string>& errmsg);
  virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);
  //READ_OFILE_MODE _read_ofile_mode = READ_OFILE_MODE::READALL;// Hwang Shih added;

  void setSite(string site) { _encode_site = site; };

 protected:
  virtual  int _decode_head();
  virtual  int _decode_data();

 private:

	 string  _encode_site;

};

} // namespace

#endif
