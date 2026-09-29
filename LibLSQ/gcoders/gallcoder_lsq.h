/**
*
* @verbatim
History
	-1.0 Jiande Huang  2019-10-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gallcoder_lsq.h
* @brief	main decode function for all files used in podlsq
*
*
* @author   Jiande Huang, Wuhan University
* @version	1.0.0
* @date		2019-10-25
*
*/


#ifndef GALLCODER_LSQ_H
#define GALLCODER_LSQ_H

#include "gdata/gdata.h"
#include "gcoders/gcoder.h"
#include "gset/gsetinp.h"
#include "gall/gallobj.h"
#include "gutils/ExportLSQ.h"
using namespace gnut;
using namespace std;

namespace great
{
	class LSQ_LIBRARY_EXPORT t_gallcoder_lsq : public virtual t_gallcoder
	{
	public:
		t_gallcoder_lsq();
		virtual ~t_gallcoder_lsq();

		bool support(const IFMT& ifmt) override;
		void coder(const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj) override;
		void creat(const IFMT& ifmt, map<t_gdata::ID_TYPE, t_gdata*>& map_gdata) override;
	private:
		std::set<IFMT> _lsq_support;
		string _class_id = "t_gallcoder_lsq";
		int    _data_id = 0;
	};
}

#endif