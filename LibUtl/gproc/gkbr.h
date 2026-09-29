/**
*
* @verbatim
	History
	 -1.0 Hongmin Zhang		2020-03-20  creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gkbr.h
* @brief		KBR
*
* @author       Hongmin Zhang, Wuhan University
* @version		1.0.0
* @date		    2020-03-20
*
*/

#ifndef G_KBR_H
#define G_KBR_H

#include "gio/glog.h"
#include "gio/giof.h"
#include "gio/gfile.h"

#include "gall/gallorb.h"
#include "gall/gallobj.h"
#include "gall/gallkbr.h"

#include "gdata/gorb.h"

#include "gproc/glsq.h"

#include "newmat/newmat.h"
#include "gutils/ExportUtl.h"
#include "gutils/gtrs2crs.h"
#include "gutils/gscf2crsleo.h"
#include "gutils/gcrs2acr.h"
#include "gutils/gscf2crs.h"

#include "gmodels/gtideIERS.h"

using namespace std;

namespace great
{
	/**
	*@brief	   Class for slr validation
	*/
	class UTL_LIBRARY_EXPORT t_gkbr
	{
	public:
		/** @brief constructor. */
		t_gkbr();

		/**
		* @brief  constructor.
		* @param[in]	set    setting information
		* @param[in]	log    log file
		*/
		t_gkbr(t_gsetbase* set, t_glog* log);

		void setORB(t_gallorb* gorb);
		void setKBR(t_gallkbr* gkbr);
		void setLEO(set<string>& leo);

		/** @brief destructor. */
		virtual ~t_gkbr();

		bool processBatch(t_gtime& beg, t_gtime& end);

		void GenerateProduct();
		
	protected:
		// For log
		t_glog*				_glog;             ///< Genereal   log output. 

		t_gtide*            _tide = nullptr;
		// For set
		t_gsetbase*			_set;              ///< Base setting. 

		t_gallorb*          _gorb;

		t_gallkbr*        _gkbr;

		//t_gkbromc*        _gkbromc;

		t_gtrs2crs*  _trs2crs_cord;              ///< class of TRS to CRS

		t_gscf2crsleo*	 _scf2crsleo_cord;       ///< class of SCF to CRS 

		t_gscf2crs*      _scf2crs_cord;

		t_gcrs2acr*     _crs2acr;

		t_gtrs2crs*  _trs2crs_time;              ///< time of TRS to CRS

		t_gtime             _start;

		t_gtime             _stop;

		t_gtrs2crs*    _trs2crs_2000 = nullptr;

		t_gpoleut1*     _allerp = nullptr;

		shared_ptr<t_gobj>	_crt_obj;

		set<string>     _leo;
	};
}

#endif