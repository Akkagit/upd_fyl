/**
*
* @verbatim
	History
	 -1.0 liyuhao		2020-07-26  creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gclkdif.h
* @brief		clock differ processing for GNSS satellite
* @author       liyuhao, Wuhan University
* @version		1.0.0
* @date		    2020-07-26
*
*/
#ifndef GCLKDIF_H
#define GCLKDIF_H

#include "gio/glog.h"
#include "gio/giof.h"

#include "gall/gallprec.h"
#include "gall/gallorb.h"

#include "gdata/gorb.h"

#include "gproc/glsq.h"

#include "newmat/newmat.h"
#include "gutils/ExportUtl.h"

using namespace std;

namespace great
{
	/**
	*@brief	   Class for GNSS clock differ processing
	*/
	class UTL_LIBRARY_EXPORT t_gclkdif
	{
	public:
		/** @brief constructor. */
		t_gclkdif();

		/**
		* @brief  constructor.
		* @param[in]	set    setting information
		* @param[in]	log    log file
		*/
		t_gclkdif(t_gsetbase* set, t_glog* log);

		/** @brief destructor. */
		virtual ~t_gclkdif();

	public:
		/**
		* @brief  batch processing for GNSS clock differ
		* @param[in]  beg         begin time from XML
		* @param[in]  end         end time from XML
		* @return
			@retval true	clock differ correctly
			@retval false   clock differ failed
		*/
		bool processBatch();

		/**
		* @brief  set rinexc data
		* @param[in]  gorb rinexc data
		*/
		void setCLKPRD(t_gallnav* gclkprd);
		void setCLKREF(t_gallnav* gclkref);

	protected:
		// For log
		t_glog* _glog;             ///< Genereal   log output. 

		t_giof* _clkdif;		   ///< write clock differ result
		// For set
		t_gsetbase* _set;              ///< Base setting. 

		// data for calulate
		t_gallnav* _gclkprd;				///< rinexc product data.
		t_gallnav* _gclkref;				///< rinexc reference data.

		t_gtime					_beg;			///< begin time from xml
		t_gtime					_end;			///< end time from xml
		double					_sampling;		///< sampling
		set<string>				_sat;			///< satellite prn.
		string					_refsat;		///<reference satallite

		int						_nepoch;		///< number of epochs
		int						_numsat;		///< number of satallites
		vector<vector<double>>  _res;			///< clock differ result
		vector<vector<double>>	_res_rms;		///< rms statistic result
		vector<vector<int>>     _flag;			///< whether the data is valid
		vector<int>				_res_del;		///< number of delete epochs
		double					_mean_rms;
		double					_mean_std;
		int						_ncount;		///< ncount for the position of refsat

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
		bool _get_setting();

		bool _set_out();

		bool _check_time(string sat);

		bool _get_sat();

		bool _init_par();

		bool _rms_statistics();

		bool _generate_result();
	};
}
#endif // !GCLKDIF_H#pragma once