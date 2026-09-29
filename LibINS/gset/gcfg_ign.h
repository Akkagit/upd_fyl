#ifndef GCFG_IGN_H
#define GCFG_IGN_H

/**
* @file gcfg_ign.h
* @brief implements settings for testing lci
* @mainpage
* @author ShenZhiheng
* @email
* @version 1.0
* @date 2019-12-28
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
	* @class t_gcfg_ign
	* @brief lci configure
	* @author ShenZhiheng
	*/
	class INS_LIBRARY_EXPORT t_gcfg_ign : public virtual t_gcfg_ppp,
					  public virtual t_gsetins,
					  public virtual t_gsetign
	{
	public:
		t_gcfg_ign();
		~t_gcfg_ign();
		
		void check();
		void help();

	protected:
		
	};

}



#endif