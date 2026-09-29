#ifdef USE_VIS
#ifndef GCFG_MULTIGN_H
#define GCFG_MULTIGN_H

/**
* @file gcfg_multign.h
* @brief implements settings for testing mult ign
* @mainpage
* @author LiShengyu
* @email
* @version 1.0
* @date 2020-06-18
* @license
*/

#include "gcoders/imufile.h"
#include "gcoders/posfile.h"
#include "gcoders/imgfile.h"
#include "gset/gsetins.h"
#include "gset/gsetign.h"
#include "gset/gcfg_ppp.h"
//#include "gset/gcfg_rtk.h"
#include "gset/gsetvis.h"
#include "gutils/ExportINS.h"


using namespace std;
using namespace gnut;
using namespace great;
using namespace gvis;

namespace gsins
{
	/**
	* @class t_gcfg_multign
	* @brief ppp/ins/vis configure
	* @author LiShengyu
	*/
	class INS_LIBRARY_EXPORT t_gcfg_multign : public t_gcfg_ppp,
		public t_gsetins,
		public t_gsetign,
		public t_gsetvis
	{
	public:
		t_gcfg_multign();
		~t_gcfg_multign();

		void check();
		void help();

	protected:

	};

}



#endif
#endif // USE_VIS