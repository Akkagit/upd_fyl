#ifndef GALLCODER_UTL_H
#define GALLCODER_UTL_H

#include "gdata/gdata.h"
#include "gcoders/gcoder.h"
#include "gset/gsetinp.h"
#include "gall/gallobj.h"
#include "gcoders/gallcoder.h"
#include "gutils/ExportUtl.h"

using namespace std;
using namespace gnut;
using namespace great;
namespace great
{
	class UTL_LIBRARY_EXPORT t_gallcoder_utl : public virtual t_gallcoder
	{
	public:
		t_gallcoder_utl();
		virtual ~t_gallcoder_utl();

		bool support(const IFMT& ifmt) override;
		void coder(const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj) override;
		void creat(const IFMT& ifmt, map<t_gdata::ID_TYPE, t_gdata*>& map_gdata) override;

	protected:
		std::set<IFMT>   _base_support;
		string           _class_id;
		int              _data_id;
	};
}

#endif