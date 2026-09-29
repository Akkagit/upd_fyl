#ifdef USE_LIDAR
#ifndef GCFG_MULTIGN_LIDAR_H
#define GCFG_MULTIGN_LIDAR_H

/**
* @file gcfg_multign_lidar.h
* @brief implements settings for testing mult ign
* @mainpage
* @author WangHuidan
* @email
* @version 1.0
* @date 2020-10-13
* @license
*/

#include "gcoders/imufile.h"
#include "gcoders/posfile.h"
#include "gcoders/lidarfile.h"
#include "gset/gsetins.h"
#include "gset/gsetign.h"
#include "gset/gcfg_ppp.h"
//#include "gset/gcfg_rtk.h"
#include "gset/gsetlidar.h"
#include "gutils/ExportLIDAR.h"


using namespace std;
using namespace gnut;
using namespace great;
using namespace glidar;

namespace gsins
{
	/**
	* @class t_gcfg_multign
	* @brief ppp/ins/lidar configure
	* @author WangHuidan
	*/
	class INS_LIBRARY_EXPORT t_gcfg_multign_lidar : public t_gcfg_ppp,
		public t_gsetins,
		public t_gsetign ,
		public t_gsetlidar
	{
	public:
		t_gcfg_multign_lidar();
		~t_gcfg_multign_lidar();

		void check();
		void help();

	protected:

	};

}



#endif
#endif // USE_LIDAR