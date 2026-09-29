#ifndef GCFG_LCI_H
#define GCFG_LCI_H

//todo  more description
/**
* @file gcfg_lci.h
* @brief implements settings for testing lci
* @mainpage
* @author ShenZhiheng
* @email
* @version 1.0
* @date 2018-12-25
* @license
*/


#include "gdata/gimudata.h"
#include "gcoders/imufile.h"
#include "gio/gtcp.h"
#include "gcoders/bncobs.h"
#include "gcoders/bnccorr.h"
#include "gdata/gposdata.h"
#include "gcoders/posfile.h"
#include "gset/gsetgen.h"
#include "gset/gsetinp.h"
#include "gset/gsetout.h"
#include "gset/gsetproc.h"
#include "gset/gsetgnss.h"
#include "gset/gsetflt.h"
#include "gset/gsetrec.h"
#include "gset/gsetins.h"
#include "gset/gsetign.h"
using namespace std;
using namespace gnut;

namespace gsins
{
	/**
	* @class t_gcfg_lci
	* @brief lci configure
	* @author ShenZhiheng
	*/
	class INS_LIBRARY_EXPORT t_gcfg_lci : public t_gsetgen,
		              public t_gsetinp,
		              public t_gsetout,
		              public t_gsetgnss,
					  public t_gsetins,
					  public t_gsetign
	{
	public:
		t_gcfg_lci();
		~t_gcfg_lci();
		
		void check();
		void help();

	protected:
		
	};

}



#endif