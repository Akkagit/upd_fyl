/**
*
* @verbatim
History
-1.0 xjhan  2019-11-22  creat the file.

@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		sp3leo.h
* @brief	read sp3 containing leo data.
*
* @author   xjhan , Wuhan University
* @version	1.0.0
* @date		2019-11-22
*
*/

#ifndef SP3LEO_H
#define SP3LEO_H

#include <vector> 

#include "gcoders/gcoder.h"
#include "gutils/gtime.h"

using namespace std;

namespace gnut {

	class BASE_LIBRARY_EXPORT t_sp3leo : public t_gcoder {

	public:
		t_sp3leo(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);
		virtual ~t_sp3leo() {};

		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg);
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

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
		vector<string> _prn;
		vector<int>    _acc;
		vector<string> _timesys;
		vector<int>    _accbase;
	};

} // namespace

#endif

