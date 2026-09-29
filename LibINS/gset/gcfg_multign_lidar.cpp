#ifdef USE_LIDAR
#include "gcfg_multign_lidar.h"

gsins::t_gcfg_multign_lidar::t_gcfg_multign_lidar() :
	t_gcfg_ppp(),
	t_gsetins(),
	t_gsetign(),
	t_gsetlidar()
{
	_IFMT_supported.insert(IMU_INP);
	_IFMT_supported.insert(IMU_POS_INP);
	_IFMT_supported.insert(LIDAR_INP);

	_OFMT_supported.insert(INS_OUT);
	_OFMT_supported.insert(INSKF_OUT);
	_OFMT_supported.insert(INSKFPK_OUT);
}

gsins::t_gcfg_multign_lidar::~t_gcfg_multign_lidar()
{
}

void gsins::t_gcfg_multign_lidar::check()
{
	t_gcfg_ppp::check();
	t_gsetins::check();
	t_gsetign::check();
	t_gsetlidar::check();
}

void gsins::t_gcfg_multign_lidar::help()
{
	t_gcfg_ppp::help();
	t_gsetins::help();
	t_gsetign::help();
	t_gsetlidar::help();
}

#endif // USE_LIDAR