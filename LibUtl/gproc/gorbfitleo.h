/**
*
* @verbatim
	History
	 -1.0 Yujie Qin		2019-12-12 creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gorbfitleo.h
* @brief		orbit fit processing for leos
*
* @author       Yujie Qin, Wuhan University
* @version		1.0.0
* @date		    2019-12-12
*
*/

#ifndef G_ORBFITLEO_H
#define G_ORBFITLEO_H

#include "gio/glog.h"
#include "gio/giof.h"
#include "gall/gallics_leo.h"
#include "gall/gallprec.h"
#include "gall/gallorb.h"

#include "gdata/gics_leo.h"
#include "gdata/gorb.h"
#include "gdata/gpoleut1.h"

#include "newmat/newmat.h"
#include "gproc/glsq.h"


#include "gutils/ExportUtl.h"

using namespace std;

namespace great
{
	/**
	*@brief	   Class for LEO orbit fit processing
	*/
	class UTL_LIBRARY_EXPORT t_gorbfitleo 
	{
	public:
		/** @brief constructor. */
		t_gorbfitleo();

		/**
		* @brief  constructor.
		* @param[in]	set    setting information
		* @param[in]	log    log file
		*/
		t_gorbfitleo(t_gsetbase* set, t_glog* log);

		/** @brief destructor. */
		virtual ~t_gorbfitleo();


	public:
		/**
		* @brief  batch processing for LEO orbit fit
		* @return
			@retval true	orbit fit correctly
			@retval false   orbit fit failed
		*/
		bool processBatch();
        /**
		* @brief  set orb funct data
		* @param[in]  gorb   orb funct data
		*/
		void setORB(t_gallorb*			gorb);
		/**
		* @brief set sp3  data
		* @param[in]  gsp3   sp3 data
		*/
		void setSP3(t_gallprec*			gsp3);
		
		void setICS(t_gallics_leo*		gics);
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

		t_gallorb*          _gorb;				///< orb data.
		t_gallprec*			_gsp3;				///< sp3 data
		t_gpoleut1*			_gpoleut1;			///< poleut1 data
        t_gallics_leo*		_gicsleo;
		t_glsq*				_lsqEstimator;		///< least square estimator
		map<string, map<t_gtime, t_gtriple > >	_res;
		t_gtime					_beg;			///< begin time from xml
		t_gtime					_end;			///< end time from xml
		t_gtime					_ref_orb;		///< reference time from ics
		double					_sampling;		///< sampling
		set<string>				_sat;			///< satellite prn.
		string					_trans;			///< reference frame
		double                  _unit;
		bool                    _isfit;
		bool                    _isinter;
		int                     _icspar = 0;
#ifdef BMUTEX
		boost::mutex  _mutex;
#endif
		t_gmutex      _gmutex;

	protected:
		bool  _init_data();	
		bool _get_setting();
		void _read_resfile();
		bool _check_time(string sat);
		bool _init_leo();
		void _init_ERP_par();
		//void _init_SAT_par(string  satname);
		void _set_out();
		void _update_icsfile();
		void _num2sstream(ostringstream&, int, int);
		void _recover_resid_leo();
		void _rms_statistics(map<string, t_gtriple>& max_acr, map<string, t_gtriple>& min_acr,
			map<string, t_gtriple> &fit_acr, map<string, t_gtriple> &mean_acr, map<string, double>	&rms, double &rmsmean);


		
	};
}
#endif
