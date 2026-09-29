#ifndef GALLCODER_ORB_LEO_H
#define GALLCODER_ORB_LEO_H

#include "gdata/gdata.h"
#include "gcoders/gcoder.h"
#include "gset/gsetinp.h"
#include "gall/gallobj.h"
#include "gutils/ExportOI.h"
using namespace gnut;
using namespace std;

namespace great
{
	class OI_LIBRARY_EXPORT t_gallcoder_orb_leo : public virtual t_gallcoder
	{
	public:
		t_gallcoder_orb_leo();
		virtual ~t_gallcoder_orb_leo();

		bool support(const IFMT& ifmt) override;
		void coder(const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj) override;
		void creat(const IFMT& ifmt, map<t_gdata::ID_TYPE, t_gdata*>& map_gdata) override;
	private:
		std::set<IFMT> _orb_support;
		string _class_id = "t_gallcoder_orb_leo";
		int    _data_id = 0;
	};
}

#endif