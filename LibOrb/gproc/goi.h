/**
*
* @verbatim
    History
     -1.0 Yujie Qin  2019-05-23  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		        goi_leo.h
* @brief		    orbit integration processing for LEO satellite
*
* @author           Yujie Qin, Wuhan University
* @version		    1.0.0
* @date		        2019-05-23
*
*/

#ifndef GOI_H
#define GOI_H

#include "gutils/gtime.h"
#include "gutils/gmutex.h"
#include "gutils/gnss.h"

#include "gio/glog.h"
#include "gio/giof.h"

#include "gall/gallprod.h"
#include "gdata/gnavde.h"
#include "gall/gallpannel.h"
#include "gall/gallobj.h"
#include "gall/gallics_leo.h"
#include "gall/gallics_nav.h"
#include "gall/gallsatparam.h"
#include "gall/gallattitude.h"

#include "gdata/gegmcoef.h"
#include "gdata/gics_leo.h"
#include "gdata/gics_nav.h"
#include "gdata/gorb.h"
#include "gdata/gleapsecond.h"
#include "gdata/gsatparam.h"
#include "gdata/gotdata.h"
#include "gdata/gsolargeomag.h"
#include "gdata/gdeadata.h"
#include "gdata/gpoleut1.h"
#include "gdata/galbedo.h"
#include "gutils/grkf.h"
#include "gutils/grkfnav.h"
#include "gutils/gadams.h"
#include "gset/gsetinp.h"

#include "newmat/newmat.h"

#include <vector>



using namespace std;

namespace great
{
	typedef map<int, vector<double> > t_map_vector;
	/**
	*@brief	   Class for all sats orbit integration processing
	*/
	class OI_LIBRARY_EXPORT t_goi
	{
	public:
		/** @brief constructor. */
		t_goi();

		/**
	    * @brief constructor.
	    * @param[in]	mark   leo name
	    * @param[in]	set    setting information
	    */
		t_goi(string mark, t_gsetbase* set);

		/** @brief destructor. */
		virtual ~t_goi();

	public:
		/**
	    * @brief  batch processing for OI
	    * @param[in]  beg   begin time from XML
	    * @param[in]  end   end time from XML
	    * @return           calculate correctly or not
	    */
		virtual bool processBatch(const t_gtime& beg, const t_gtime& end);
		

		/**
		* @brief set  ICS data
		* @param[in] gics    ICS data
		*/
		virtual void set_data_used(t_gallics_nav* gics);

		/**
		* @brief set LEO ICS data
		* @param[in] gics   LEO ICS data
		*/
		virtual void set_data_used(t_gallics_leo* gics);

		/**
	    * @brief  set EGM data
	    * @param[in]  gegmcoef   EGM data
	    */
		void set_data_used(t_gegmcoef* gegmcoef);

		/**
	    * @brief set ocean load data
	    * @param[in]  gotl   ocean load data
	    */
		void set_data_used(t_gallotl*			    gotl);

		virtual void set_data_used(t_gdeadata*	 gotl);
		/**
	    * @brief   set leap second data
	    * @param[in]  leapsecond   leap second data
	    */
		void set_data_used(t_gleapsecond*		    leapsecond);

		/** @brief set ocean_tide data */
		void set_data_used(t_gocean_tide*		    gotdata);

		/**
	    * @brief   set ocean_tide data
	    * @param[in]  gotdata   ocean tide data
	    */
		void set_data_used(t_gnavde*        gplanet);

		/**
	    * @brief   set LEO pannel data
	    * @param[in] gpannel  pannel data for one LEO satellite
	    */
		void set_data_used(t_gpannel*        gpannel);

		/**
	    * @brief set poleut1 data
	    * @param[in]  poleut1   poleut1 data
	    */
		void set_data_used(t_gpoleut1*         poleut1);

		/**
	    * @brief   set sat_parameters_new data
	    * @param[in] gallsat   data for all satellites
	    */
		void set_data_used(t_gallsatinfo*	    gallsat);

		/**
		* @brief   set solar_flux and geomag_kp data
		* @param[in] solargeomag   solar_flux and geomag_kp data
		*/
		virtual void set_data_used(t_gsolargeomag* solargeomag);

		/**
		* @brief   set albedo emit and refl data
		* @param[in] albedodata   emit and refl data
		*/
		virtual void set_data_used(t_galbedo* albedodata);

		/**
		* @brief    set LEO attitude data
		* @param[in] attitude  attitude data for one LEO satellite
		*/
		virtual void set_data_used(t_gattitude*       attitude);


		void set_data_used(t_gorb * gorb);
		
		/** @brief set general log file. */
		void glog(t_glog* glog) { _glog = glog; };

		
		
	protected:
		// For log
		t_glog*				_glog               = nullptr;           ///< Genereal   log output. 

		// For set
		t_gsetbase*			_set                = nullptr;          ///< Base setting. 

		string				_sat;                                   ///< satellite name
		
		// data 													// data for calulate
		t_gnavde*           _gplanet            = nullptr;          ///< data for all planets
		t_gpoleut1*         _gpoleut1           = nullptr;          ///< poleut1 data
		t_gocean_tide*      _gotdata            = nullptr;          ///< ocean tide data
		t_gleapsecond*      _gleapsecond        = nullptr;          ///< leap second data
		t_gegmcoef*			_gegmcoef           = nullptr;		    ///< EGM data
		t_gallsatinfo*      _gallsat            = nullptr;          ///< satellite data
		t_gallotl*			_gotl               = nullptr;          ///< ocean load data
		
		t_gsolargeomag*     _gsolargeomag       = nullptr;          ///< solar_flux and geomag_kp data
		t_galbedo*          _galbedo            = nullptr;          ///< emit and refl data
		t_gpannel*          _gpannel            = nullptr;          ///< LEO pannel data
		t_gattitude*        _gattitude          = nullptr;          ///< LEO attitude data
		t_gdeadata*         _gdeadata           = nullptr;
		
		// data for result
		t_gsatparam 		_gsat;			   ///< one satelite data
		t_gorb*             _gorb               = nullptr;          ///< orb data (output)
		t_gforce_calarg*    _arg_base           = nullptr;          ///< force argument
		
	    // save for oi settings
		int                     _step_inp=0;       ///< OI input step(unit: sec)
		int                     _step_out=0;       ///< OI output step(unit: sec)
		int                     _integra_degree=0; ///< OI degree
		int                     _gravity_degree=0; ///< degree of gravity model

		t_gtime                 _beg;            ///< begin time for OI
		t_gtime                 _end;            ///< end time for OI
		t_gtime                 _ref;            ///< reference time of OI
		set<FORCE>				_forces;		 ///< List of considered forces
		map<FORCE, bool>		_force_valid;   ///< map of forces(true or false)
		map<FORCE, string>		_force_model;	 ///< List of models for considered forces

		// for intergra info
		t_gtime                 _epoch;          ///< current epoch
		t_gtime					_extend_beg;     ///< extended begin time
		t_gtime					_extend_end;     ///< extended end time

		double                  _epoch_mjd      = 0.0; ///< modified julian day of current epoch
		double                  _extend_beg_mjd = 0.0; ///< extended begin time(MJD)
		double                  _extend_end_mjd = 0.0; ///< extended end time(MJD)

		int                     _neq_size=0;       ///< number of OI funct params
		int					    _par_num =0;       ///< number of params
		double                  _sat_mass = 0.0;   ///< mass of satellite
		double                  _period   = 0.0;   ///< period of satellite


	    //save the state for writing orb
		t_map_vector          _state_back;             ///< state data of back integration
		t_map_vector          _state_frwd;             ///< state data of forward integration

		map<ICS_PAR, double>     _xics;                ///< par values from ICS data
		
		
		
		/**
	    * @brief numerical integration for OI
	    * @param[in]  beg          extended begin time
	    * @param[in]  end          extended end time
	    * @param[in]  state_init   initial state data
	    * @param[out] state_result the state data after the numerical integration
	    * @return bool
			@retval true     numerical integrate correctly
			@retval false    numerical integrate incorrectly
	   */
		bool  _numerical_integration(const t_gtime& beg,const t_gtime& end,
			const vector<double>& state_init, t_map_vector& state_result);
			
	    /**
	    * @brief OI processsing for one epoch.
	    * @param[in]  epoch          current epoch
	    * @param[in]  adams          adams integrator
	    * @param[out] state_result   state result after calculation
	    * @param[out] partial_result partial result after calculation
	    * @return bool
			@retval true   calculate correctly
			@retval false  calculate incorrectly
	    */
		bool    _processEpoch(t_gtime epoch, t_adams* adams,
			vector<list<double> >& state_result, vector<list<double> >& partial_result);
	
		
		/**
	    * @brief   extend the time period
	    * @param[in]  beg         begin time according to the XML file
	    * @param[in]  end         end time according to the XML file
	    * @param[in]  step_out    OI output step
	    * @param[out] extend_beg  extended begin time
	    * @param[out] extend_end  extended end time
	    * @return bool
			@retval true     extend time correctly
			@retval false    extend time incorrectly
	    */
		bool  _extend_time_period(const t_gtime& beg, const t_gtime& end, const int& step_out, t_gtime& extend_beg, t_gtime& extend_end);

		/**
	    * @brief   prepare adams with rkf
	    * @param[in]  adams_integra_degree      order of adams integration
	    * @param[in]  epoch     begin time of integration
	    * @param[in]  step_h    step of rkf
	    * @param[in]  rkf_init_state    initial state of rkf integration
	    * @param[out] adams_state		initial state of adams integration
	    * @param[out] adams_partial		initial partial of adams integration
	    * @return bool
			@retval true    prepare adams with rkf correctly
			@retval false   prepare adams with rkf incorrectly
	    */
		bool  _prepare_adams_with_rkf(const int& adams_integra_degree,const t_gtime& epoch,const double& step_h, const vector<double>& rkf_init_state, vector<list<double> >& state_adams, vector<list<double> >& admas_partial);

		/**
		* @brief  calculate the acc and partial in epoch t
		* @param[in]  t        current epoch for calculation
		* @param[out] state    state data after calculation
		* @param[out] partial  partial data after calculation
		* @return bool
			@retval true     calculate correctly
			@retval false    calculate incorrectly
		*/
		virtual bool _calculate(const t_gtime& t, const vector<double>& state, vector<double>& partial) = 0;

	};
}
#endif // !GOI_H
