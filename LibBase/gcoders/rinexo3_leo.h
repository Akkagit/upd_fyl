/**
*
* @verbatim
History
-1.0 xjhan  2019-11-22  creat the file.

@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		rinexo3_leo.h
* @brief	read rinexo3 containing leo data.
*
* @author   xjhan , Wuhan University
* @version	1.0.0
* @date		2019-11-22
*
*/

#ifndef RINEXO3_LEO_H
#define RINEXO3_LEO_H


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
#include "gdata/grnxhdr.h"

using namespace std;

namespace gnut {

	class BASE_LIBRARY_EXPORT t_rinexo3_leo : public t_rinexo2 {

	public:
		t_rinexo3_leo(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);
		virtual ~t_rinexo3_leo() {};

		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg) = 0;
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg) = 0;

		virtual int encode_head(char* buff, int sz, vector<string>& errmsg) = 0;
		virtual int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg) = 0;

	protected:

		virtual int _decode_head();
		virtual int _decode_data();

		virtual int _check_head();                                           // fill header information
		virtual int _read_epoch();                                           // read epoch & number of satellites, return flag
		virtual int _read_obstypes(const string& sat, const string& sys);    // read single satellite observation types
		virtual int _fix_band(string sys, string& go);                       // fix band (BDS)

		t_rnxhdr::t_obstypes       _mapcyc;     // map of GOBS phase quater-cycle shifts
		t_rnxhdr::t_obstypes       _glofrq;     // map of GLONASS slot/frequency
		t_rnxhdr::t_vobstypes      _globia;     // vec of GLONASS obs code-phase biases


	private:

	};

} // namespace

#endif
