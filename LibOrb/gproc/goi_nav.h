/**
*
* @verbatim
	History
	 -1.0 Yujie Qin		2019-01-10  creat
	 -1.1 Yujie Qin		2019-02-18  optim
	 -1.2 Jiande Huang	2019-02-22  final 1.0 version
	 -1.2 Jie Li		2019-03-31  add Doxygen style code
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    goi_nav.h
* @brief		orbit integration processing for GNS satellite
*
* @author       Yujie Qin, Wuhan University
* @version		1.0.0
* @date		    2019-01-10
*
*/


#ifndef GOI_NAV_H
#define GOI_NAV_H

#include "gutils/gtime.h"
#include "gutils/gmutex.h"
#include "gutils/gnss.h"
#include "gutils/gtrs2crs.h"
#include "gutils/gscf2crs.h"
#include "gutils/gcrs2acr.h"

#include "gio/glog.h"
#include "gio/giof.h"

#include "gall/gallprod.h"
#include "gall/gallpannel.h"
#include "gall/gallics_nav.h"
#include "gall/gallsatparam.h"
#include "gdata/gnavde.h"
#include "gdata/gegmcoef.h"
#include "gdata/gics_nav.h"
#include "gdata/gorb.h"
#include "gdata/gleapsecond.h"
#include "gdata/gsatparam.h"
#include "gdata/gotdata.h"
#include "gforce/gforcenav.h"
#include "gutils/gadams.h"
#include "gutils/grkf.h"
#include "gutils/grkfnav.h"
#include "gset/gsetinp.h"
#include "newmat/newmat.h"
#include "gproc/goi.h"
#include "gutils/ExportOI.h"

using namespace std;

namespace great
{
	typedef map<int, vector<double> > t_map_vector;

	/**
	*@brief	   Class for GNS orbit integration processing
	*/
	class OI_LIBRARY_EXPORT t_goi_nav :public t_goi
	{
	public:
		/** @brief constructor. */
		t_goi_nav() {};

		/** 
		* @brief  constructor.
		* @param[in]	mark   satellite prn
		* @param[in]	set    setting information
		*/
		t_goi_nav(string mark, t_gsetbase* set);

		/** @brief destructor. */
		virtual ~t_goi_nav();

	public:
		/**
		* @brief  batch processing for GNS OI
		* @param[in]  beg         begin time from XML
		* @param[in]  end         end time from XML
		* @return
			@retval true	OI correctly
			@retval false   OI failed
		*/
		bool processBatch(const t_gtime& beg, const t_gtime& end); 

		/**
		* @brief set  ICS data
		* @param[in] gics    ICS data
		*/
		void set_data_used(t_gallics_nav*		gics)  override;

	protected:
		//// data for calulate
		t_gallics_nav*		_gics        = nullptr;		    ///< all  ics data.
		t_gics_nav*			_sat_ics     = nullptr;			///< one satellite ics data.
		
		// data for result
		t_gforce_nav        _gforce;			            ///< force model for GNS.
		GSYS			    _sys;			                ///< satellite system.

		//// for trs to crs(j2000) and 
		t_gscf2crs*	 _scf2crs_cord = nullptr;				///< class of SCF to CRS.
		t_gcrs2acr*  _crs2acr_cord = nullptr;				///< class of CRS to ACR.
		t_gtrs2crs*  _trs2crs_cord = nullptr;				///< class of TRS to CRS.	

		//// satellite mass
		double          _gnsmass = 0.0;               ///< satellite mass for physical model e.g. SRP, ERP and Antenna thrust

	protected:
		/**
		* @brief   Get settings from XML file and set local variables.
		* @return
			@retval true	get settings correctly
			@retval false   get settings failed
		*/
		bool  _get_settings();

		/**
		* @brief   check the data needed.
		* @return
			@retval true	data is completed
			@retval false   some data is lack
		*/
		bool  _check_data();

		/** @brief Set output file (empty). */
		void  _setOut();

		/** @brief Set processing log file (<ppp>). */
		void  _setLog();

		

		/**
		* @brief   initial the sat infomation
		* @param[in] sat_name   satellite name
		* @return
			@retval true	initial sat infomation successfully
			@retval false   initial sat infomation unsuccessfully
		*/
		bool  _init_sat_info(string sat_name);

		/**
		* @brief      initial the state matrix
		* @param[in]  neq_size   number of OI funct params(£¨np+6£©*6+6)
		* @param[out] state      state data of OI processing
		* @return
			@retval true	initial state matrix successfully
			@retval false   initial state matrix unsuccessfully
		*/
		bool  _init_state_matrix(int neq_size, vector<double>& state);

		
		/**
		* @brief  calculate the acc and partial in epoch t
		* @param[in]  t			current epoch
		* @param[in]  state     state at current epoch
		* @param[out] partial   partial of the state 
		* @return
			@retval true	calculate successfully
			@retval false   calculate unsuccessfully
		*/
		bool _calculate(const t_gtime& t,const vector<double>& state, vector<double>& partial) override;
		
		/**
		* @brief  creat initial information of force model
		* @param[in]  force_ID			force type
		* @param[in]  force_contain     force model initial information
		* @return
			@retval true	creat force model successfully
			@retval false   creat force model unsuccessfully
		*/
		bool _creat_force_model(FORCE force_ID, t_gforce_nav* force_contain);
	};
}
#endif // !GOI_H
