#ifndef GCFG_SIMUINS_H
#define GCFG_SIMUINS_H

/**
* @file gcfg_simuins.h
* @brief implements settings for simulating GNSS/INS data
* @mainpage
* @author ZhouYuxuan
* @email
* @version 1.0
* @date 2020-07-03
* @license
*/

#include "gcoders/imufile.h"
#include "gcoders/posfile.h"
#include "gset/gsetins.h"
#include "gset/gsetign.h"
#include "gset/gcfg_ppp.h"
#include "gutils/ExportINS.h"


using namespace std;
using namespace gnut;
using namespace great;

namespace gsins
{
	/**
	* @class t_gcfg_simuins
	* @brief GNSS/INS simulation configure
	* @author ShenZhiheng
	*/
	class INS_LIBRARY_EXPORT t_gcfg_simuins : public t_gcfg_ppp,
		public t_gsetins,
		public t_gsetign
	{
	public:
		t_gcfg_simuins();
		~t_gcfg_simuins();

		void check();
		void help();
		int generateXML(stringstream & ss);

	protected:

	};

}



#endif