/**
*
* @verbatim
History
-1.0 xjhan  2019-11-22  creat the file.

@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		rinexc_leo.h
* @brief	read rinexc containing leo data.
*
* @author   xjhan , Wuhan University
* @version	1.0.0
* @date		2019-11-22
*
*/

#ifndef RINEXC_LEO_H
#define RINEXC_LEO_H

#include "gcoders/rinexc.h"

#define RINEXC_BUFFER_LEN 81

using namespace std;

namespace gnut {

	class BASE_LIBRARY_EXPORT t_rinexc_leo : public t_rinexc {

	public:
		t_rinexc_leo(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);
		virtual ~t_rinexc_leo();

		virtual int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg) override;

	};

} // namespace

#endif

