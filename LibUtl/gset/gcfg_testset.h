
#ifndef GCFG_TESTSET_H
#define GCFG_TESTSET_H

#include <string>
#include <iostream>
#include <signal.h>

#include "gutils/ExportUtl.h"
#include "gset/gsetpar.h"
#include "gset/gsetintgra.h"

using namespace std;
using namespace pugi;

namespace great {

	class UTL_LIBRARY_EXPORT t_gcfg_testset : 
		public t_gsetpar,
		public t_gsetintgra
	{

	public:
		t_gcfg_testset();
		~t_gcfg_testset();

		void check();                                 // settings check
		void help();                                  // settings help

	protected:

	private:

	};

} // namespace

#endif