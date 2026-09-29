/**
*
* @verbatim
	History
	 -1.0 Jie Li		2019-09-06  creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gorbfit.h
* @brief		orbit fit processing for GNS satellite
* @author       Jie Li, Wuhan University
* @version		1.0.0
* @date		    2019-09-06
*
*/

#ifndef G_ORBFIT_H
#define G_ORBFIT_H

#include "gio/glog.h"
#include "gio/giof.h"

#include "gall/gallics_nav.h"
#include "gall/gallprec.h"
#include "gall/gallorb.h"

#include "gdata/gics_nav.h"
#include "gdata/gorb.h"
#include "gdata/gpoleut1.h"

#include "newmat/newmat.h"
#include "gproc/glsq.h"

#include "gproc/gorbdif.h"

#include "gutils/ExportUtl.h"

using namespace std;

namespace great
{
	/**
	*@brief	   Class for GNS orbit fit processing
	*/
	class UTL_LIBRARY_EXPORT t_gorbfit : public t_gorbdif
	{
	public:
		/** @brief constructor. */
		t_gorbfit();

		/**
		* @brief  constructor.
		* @param[in]	set    setting information
		* @param[in]	log    log file
		*/
		t_gorbfit(t_gsetbase* set, t_glog* log);

		/** @brief destructor. */
		virtual ~t_gorbfit();


	public:
		/**
		* @brief  batch processing for GNS orbit fit
		* @return
			@retval true	orbit fit correctly
			@retval false   orbit fit failed
		*/
		bool processBatch();

	protected:


#ifdef BMUTEX
		boost::mutex  _mutex;
#endif
		t_gmutex      _gmutex;

	protected:
		bool  _init_data();


		void _init_ERP_par();
		void _init_SAT_par();
		void _set_out();
		void _update_icsfile();
	};
}
#endif
