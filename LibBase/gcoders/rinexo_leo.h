
#ifndef RINEXO_LEO_H
#define RINEXO_LEO_H

#include <string> 
#include <vector> 

#include "gall/gallobs.h"
#include "gall/gallobj.h"
#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gutils/gsys.h"

#include "gutils/gobs.h"
#include "gutils/gtypeconv.h"
/**
*
* @verbatim
History
-1.0 xjhan  2019-11-22  creat the file.

@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		rinexo_leo.h
* @brief	read rinexo containing leo data.
*
* @author   xjhan , Wuhan University
* @version	1.0.0
* @date		2019-11-22
*
*/

#include "gcoders/gcoder.h"
#include "gcoders/rinexo2.h"
#include "gcoders/rinexo3_leo.h"
#include "gdata/grnxhdr.h"

using namespace std;

namespace gnut {

	class BASE_LIBRARY_EXPORT t_rinexo_leo : public t_rinexo3_leo {

	public:
		t_rinexo_leo(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);
		virtual ~t_rinexo_leo() {};

		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg);
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

		virtual  int encode_head(char* buff, int sz, vector<string>& errmsg);
		virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

	protected:
		virtual  int _decode_head();
		virtual  int _decode_data();

	private:

	};

} // namespace

#endif
