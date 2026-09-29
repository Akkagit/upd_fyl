#include "gcfg_ign.h"

gsins::t_gcfg_ign::t_gcfg_ign() :
	t_gcfg_ppp(),
	t_gsetins(), 
	t_gsetign()
{
	_IFMT_supported.insert(IMU_INP);
	_IFMT_supported.insert(IMU_POS_INP);

	_OFMT_supported.insert(INS_OUT);
	_OFMT_supported.insert(INSKF_OUT);
	_OFMT_supported.insert(INSKFPK_OUT);
}

gsins::t_gcfg_ign::~t_gcfg_ign()
{
}

void gsins::t_gcfg_ign::check()
{
	t_gcfg_ppp::check();
	t_gsetins::check();
	t_gsetign::check();
}

void gsins::t_gcfg_ign::help()
{
	t_gcfg_ppp::help();
	t_gsetins::help();
	t_gsetign::help();
}
