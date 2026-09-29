

#include <iomanip>
#include <sstream>

#include "gset/gcfg_testset.h"

using namespace std;
using namespace pugi;

namespace great {

	// Constructor
	// ----------
	t_gcfg_testset::t_gcfg_testset()
		: t_gsetbase(),
		t_gsetpar()
	{

	}


	// Destructor
	// ----------
	t_gcfg_testset::~t_gcfg_testset()
	{}


	// settings check
	// ----------
	void t_gcfg_testset::check()
	{
		t_gsetpar::check();
		t_gsetintgra::check();
	}

	// settings help
	// ----------
	void t_gcfg_testset::help()
	{
		t_gsetbase::help_header();
		t_gsetpar::help();
		t_gsetintgra::help();
	}

} // namespace
