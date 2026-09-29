#ifndef GALLCODER_BASE_LEO_H
#define GALLCODER_BASE_LEO_H

#include "gdata/gdata.h"
#include "gcoders/gcoder.h"
#include "gset/gsetinp.h"
#include "gall/gallobj.h"
#include "gcoders/gallcoder.h"

using namespace gnut;
namespace great
{
	class BASE_LIBRARY_EXPORT t_gallcoder_base_leo : public virtual t_gallcoder
	{
	public:
		t_gallcoder_base_leo();
		virtual ~t_gallcoder_base_leo();

		bool support(const IFMT& ifmt) override;
		void coder(const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj) override;
		void creat(const IFMT& ifmt, map<t_gdata::ID_TYPE, t_gdata*>& map_gdata) override;

		//IFMT             data2ifmt(t_gdata::ID_TYPE type) override;
		// t_gdata::ID_TYPE ifmt2data(IFMT ifmt) override;

	protected:
		std::set<IFMT>   _base_support;
		string  _class_id = "t_gallcoder_base_leo";
		int  _data_id = 0;
	};
}

#endif