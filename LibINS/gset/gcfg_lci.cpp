#include "gcfg_lci.h"

gsins::t_gcfg_lci::t_gcfg_lci() : t_gsetbase(),
t_gsetgen(),
t_gsetinp(),
t_gsetgnss(),
t_gsetout(),
t_gsetins(), 
t_gsetign()
{
	_IFMT_supported.insert(IMU_INP);
	_IFMT_supported.insert(IMU_POS_INP);

	_OFMT_supported.insert(INS_OUT);
	_OFMT_supported.insert(INSKF_OUT);
	_OFMT_supported.insert(INSKFPK_OUT);
}

gsins::t_gcfg_lci::~t_gcfg_lci()
{
}

void gsins::t_gcfg_lci::check()
{
	t_gsetgen::check();
	t_gsetinp::check();
	t_gsetout::check();
	t_gsetgnss::check();
	t_gsetins::check();
	t_gsetign::check();
}

void gsins::t_gcfg_lci::help()
{
	t_gsetgen::help();
	t_gsetinp::help();
	t_gsetout::help();
	t_gsetgnss::help();
	t_gsetins::help();
	t_gsetign::help();
}
