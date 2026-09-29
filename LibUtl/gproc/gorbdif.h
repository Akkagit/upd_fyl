/**
*
* @verbatim
	History
	 -1.0 Jie Li		2019-09-06  creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gorbdif.h
* @brief		orbit differ processing for GNS satellite
* @author       Jie Li, Wuhan University
* @version		1.0.0
* @date		    2019-09-06
*
*/

#ifndef G_ORBDIF_H
#define G_ORBDIF_H

#include "gio/glog.h"
#include "gio/giof.h"

#include "gall/gallics_nav.h"
#include "gall/gallprec.h"
#include "gall/gallorb.h"

#include "gdata/gics_nav.h"
#include "gdata/gorb.h"
#include "gdata/gpoleut1.h"

#include "gproc/glsq.h"

#include "newmat/newmat.h"
#include "gutils/ExportUtl.h"

using namespace std;

namespace great
{
	/**
	*@brief	   Class for GNS orbit differ processing
	*/
	class UTL_LIBRARY_EXPORT t_gorbdif
	{
	public:
		/** @brief constructor. */
		t_gorbdif();

		/**
		* @brief  constructor.
		* @param[in]	set    setting information
		* @param[in]	log    log file
		*/
		t_gorbdif(t_gsetbase* set, t_glog* log);

		/** @brief destructor. */
		virtual ~t_gorbdif();


	public:
		/**
		* @brief  batch processing for GNS orbit differ
		* @param[in]  beg         begin time from XML
		* @param[in]  end         end time from XML
		* @return
			@retval true	orbit differ correctly
			@retval false   orbit differ failed
		*/
		bool processBatch();

		/**
		* @brief set  ICS data
		* @param[in] gics    ICS data
		*/
		void setICS(t_gallics_nav*		gics);

		/**
		* @brief  set orb funct data
		* @param[in]  gorb   orb funct data
		*/
		void setORB(t_gallorb*			gorb);

		/**
		* @brief set sp3  data
		* @param[in]  gsp3   sp3 data
		*/
		void setSP3(t_gallnav*			gsp3);

		/**
		* @brief set sp3  data
		* @param[in]  gsp3   sp3 data
		*/
		void setPOLEUT1(t_gpoleut1*		gpoleut1);

	protected:
		// For log
		t_glog*				_glog;             ///< Genereal   log output. 

		t_giof*				_orbdif;		   ///< write orbit differ result
		// For set
		t_gsetbase*			_set;              ///< Base setting. 

		// data for calulate
		t_gallics_nav*		_gics;				///< all  ics data.
		t_gallorb*          _gorb;				///< orb data.
		t_gallnav*			_gsp3;				///< sp3 data
		t_gpoleut1*			_gpoleut1;			///< poleut1 data

		// estimator
		t_glsq*				_lsqEstimator;		///< least square estimator
		
		map<string, map<t_gtime, t_gtriple > >	_res;
		// save for orbdif settings
		t_gtime					_beg;			///< begin time from xml
		t_gtime					_end;			///< end time from xml
		t_gtime					_ref_orb;		///< reference time from ics
		double					_sampling;		///< sampling
		set<string>				_sat;			///< satellite prn.
		string					_trans;			///< reference frame
		set<string>				_excsys;
		set<string>				_excsat;
#ifdef BMUTEX
		boost::mutex  _mutex;
#endif
		t_gmutex      _gmutex;

	protected:

		/**
		* @brief   check the data needed.
		* @return
			@retval true	data is completed
			@retval false   some data is lack
		*/
		bool  _init_data();

		bool _get_setting();
		// added by yqyuan for orbfit
		bool _get_setting_orbfit();			   

		bool _check_time(string sat);

		void _init_par();

		void _set_out();

		void _recover_resid();

		void _read_resfile();

		void _rms_statistics(map<string, t_gtriple>& max_acr, map<string, t_gtriple>& min_acr,
			map<string, t_gtriple> &fit_acr, map<string, double>	&rms, double &rmsmean);
	};
}
#endif
