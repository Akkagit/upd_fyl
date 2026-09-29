/**
*
* @verbatim
History
	-1.0 Jiande Huang  2019-10-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gallcoder_orb.h
* @brief	main decode function for all files used in oi
*
*
* @author   Jiande Huang, Wuhan University
* @version	1.0.0
* @date		2019-10-25
*
*/

#ifndef GALLCODER_ORB_H
#define GALLCODER_ORB_H

#include "gdata/gdata.h"
#include "gcoders/gcoder.h"
#include "gset/gsetinp.h"
#include "gutils/ExportOI.h"
using namespace gnut;
using namespace std;

namespace great
{
	class OI_LIBRARY_EXPORT t_gallcoder_orb : public virtual t_gallcoder
	{
	public:
		t_gallcoder_orb();
		virtual ~t_gallcoder_orb();

		bool support(const IFMT& ifmt) override;
		void coder(const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj) override;
		void creat(const IFMT& ifmt, map<t_gdata::ID_TYPE, t_gdata*>& map_gdata) override;
	private:
		std::set<IFMT>    _orb_support;
		string _class_id = "t_gallcoder_orb";
		int    _data_id = 0;
	};
}

#endif