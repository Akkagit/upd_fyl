#ifdef USE_VIS
#include "gcfg_multign.h"

gsins::t_gcfg_multign::t_gcfg_multign() :
	t_gcfg_ppp(),
	t_gsetins(),
	t_gsetign(),
	t_gsetvis()
{
	_IFMT_supported.insert(IMU_INP);
	_IFMT_supported.insert(IMU_POS_INP);
	_IFMT_supported.insert(IMAGE_INP);

	_OFMT_supported.insert(INS_OUT);
	_OFMT_supported.insert(INSKF_OUT);
	_OFMT_supported.insert(INSKFPK_OUT);
}

gsins::t_gcfg_multign::~t_gcfg_multign()
{
}

void gsins::t_gcfg_multign::check()
{
	t_gcfg_ppp::check();
	t_gsetins::check();
	t_gsetign::check();
	t_gsetvis::check();
}

void gsins::t_gcfg_multign::help()
{
	t_gcfg_ppp::help();
	t_gsetins::help();
	t_gsetign::help();
	t_gsetvis::help();
}

#endif // USE_VIS