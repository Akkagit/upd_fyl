#ifndef GALLNPT_H
#define GALLNPT_H

#include <iostream>
#include <string.h>
#include "gutils/ExportUtl.h"
#include "gutils/gmutex.h"
#include "gdata/gdata.h"
#include "gdata/gnpt.h"

using namespace std;
using namespace gnut;

namespace great
{
	class UTL_LIBRARY_EXPORT  t_gallnpt : public t_gdata
	{
	public:
		t_gallnpt();
		virtual ~t_gallnpt();

		void setAllNpt(string sat, int& now, t_gnpt& gnpt);

		t_gnpt& getAllNpt(string sat, int& now);

		int getSize(string sat);

	private:

		map<string ,map<int, t_gnpt>> _gallnpt;

	};
}

#endif