/**
*
* @verbatim
History
	-1.0 Jiande Huang  2019-10-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gallcoder_amb.h
* @brief	main decode function for all files used in ambfix
*         
*
* @author   Jiande Huang, Wuhan University
* @version	1.0.0
* @date		2019-10-25
*
*/

#ifndef GALLCODER_AMB_H
#define GALLCODER_AMB_H

#include "gdata/gdata.h"
#include "gcoders/gcoder.h"
#include "gcoders/gallcoder.h"
#include "gutils/ExportAMB.h"
using namespace gnut;
namespace great
{
	class AMB_LIBRARY_EXPORT t_gallcoder_amb : public virtual t_gallcoder
	{
	public:
		t_gallcoder_amb();
		virtual ~t_gallcoder_amb();

		bool support(const IFMT& ifmt) override;
		void coder(const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj) override;
		void creat(const IFMT& ifmt, map<t_gdata::ID_TYPE, t_gdata*>& map_gdata) override;
	protected:
		std::set<IFMT>   _base_support;
		string  _class_id = "t_gallcoder_amb";
		int  _data_id = 0;
	};
}

#endif