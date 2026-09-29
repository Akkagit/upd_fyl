/**
*
* @verbatim
    History
     -1.0 Yujie Qin  2019-01-10  created
     -1.1 Yujie Qin  2019-02-18  optim
     -1.2 Wei Zhang  2019-03-27  Adding Doxygen Style Code Remarks
	 -1.3 Yujie Qin  2019-04-15  debug and optim
	 -1.4 Yujie Qin  2019-07-03  Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		        goi_leo.h
* @brief		    orbit integration processing for LEO satellite
*
* @author           Yujie Qin, Wuhan University
* @version		    1.0.0
* @date		        2019-03-27
*
*/

#ifndef GOI_LEO_H
#define GOI_LEO_H

#include "gutils/gtime.h"
#include "gutils/gmutex.h"
#include "gutils/gnss.h"
#include "gutils/gtrs2crs.h"
#include "gutils/gscf2crsleo.h"
#include "gutils/gcrs2acr.h"

#include "gio/glog.h"
#include "gio/giof.h"

#include "gall/gallprod.h"
#include "gdata/gnavde.h"
#include "gall/gallpannel.h"
#include "gall/gallobj.h"
#include "gall/gallics_leo.h"
#include "gall/gallsatparam.h"
#include "gall/gallattitude.h"

#include "gdata/gegmcoef.h"
#include "gdata/gics_leo.h"
#include "gdata/gorb.h"
#include "gdata/gleapsecond.h"
#include "gdata/gsatparam.h"
#include "gdata/gotdata.h"
#include "gdata/gsolargeomag.h"
#include "gdata/gdeadata.h"
#include "gforce/gforceleo.h"

#include "gutils/grkf.h"
#include "gutils/grkfnav.h"
#include "gutils/gadams.h"
#include "gproc/goi.h"
#include "gset/gsetinp.h"

#include "newmat/newmat.h"

#include <vector>

//#define NEXT_POINT_NUM 5
//#define RKF_DEGREE     7
//#define RKF_COEF_I     9
//#define RKF_POINT_NUM  30
//#define OCEAN_DEGREE   30
//#define GRAVITY_DEGREE 120

using namespace std;

namespace great
{
	typedef map<int, vector<double> > t_map_vector;
	/**
	*@brief	   Class for LEO orbit integration processing
	*/
	class OI_LIBRARY_EXPORT t_goi_leo :public t_goi
	{
	public:
		/** @brief constructor. */
		t_goi_leo();

		/**
	    * @brief constructor.
	    * @param[in]	mark   leo name
	    * @param[in]	set    setting information
	    */
		t_goi_leo(string mark, t_gsetbase* set);

		/** @brief destructor. */
		virtual ~t_goi_leo();

	public:
		/**
	    * @brief  batch processing for LEO OI
	    * @param[in]  beg   begin time from XML
	    * @param[in]  end   end time from XML
	    * @return           calculate correctly or not
	    */
		bool processBatch(const t_gtime& beg, const t_gtime& end);

		//todo 命名使用重载 set_data
		/**
	    * @brief set LEO ICS data
	    * @param[in] gics   LEO ICS data
	    */
		void set_data_used(t_gallics_leo* gics) override;

			
		/**
	    * @brief   get orb funct according to the satellite name
	    * @param[in]  sat   satellite name
	    * @param[in]  gorb  orb data
	    */
		void getOrbFunct(string sat, t_gorb* gorb);

		/**
	    * @brief  get orb header data according to the satellite name
	    * @param[in]  sys   system of the satelite
	    * @param[in]  sats  all satellite names
	    * @param[in]  gorb  orb data
	    */
		void getOrbHead(set<string>& sys, set<string>& sats, t_gorb* gorb);
	protected:

		int             _ocean_degree = 30;
		int             _pole_degree = 30;

		// data for calulate
		t_gics_leo*			_sat_ics = nullptr;            ///< one LEO ics data
		t_gallics_leo*		_gics    = nullptr;			   ///< all LEO ics data
		
														   ///< data for result
		t_gforce_leo        _gforce;                       ///< force model for LEO

		map<string, int>       _cus_par;
		map<string, int>       _drag_par;
		map<string, int>       _srp_par;

		bool                    _isPoletide;
		bool                    _isPannel;       ///< use pannel data or not for leo

		// for trs to crs(j2000) and 
		t_gscf2crsleo*	 _scf2crsleo_cord;       ///< class of SCF to CRS 
		t_gtrs2crs*  _trs2crs_cord;              ///< class of TRS to CRS
		t_gtrs2crs*  _trs2crs_time;              ///< time of TRS to CRS

	protected:
		/**
	    * @brief   Get settings from XML file and set local variables.
	    * @return  get settings correctly or not
		* @return bool 
			@retval true     get settings normally
			@retval false    have some error
	    */
	 	bool   _get_settings();
		
		/**
	    * @brief   check all file data.
	    * @return  get filedata correctly or not
		* @return bool
			@retval true     all data exists
			@retval false    lack of data
	    */

		bool _check_par();
		bool  _checkDataUsed();

		/** @brief Set output file (empty). */
		void  _setOut();
		
		/** @brief Set processing log file (<ppp>). */
		void  _setLog();
		
	
		/**
	    * @brief   init the sat infomation
	    * @param[in] sat_name   satellite name
	    * @return bool
			@retval true    initial leo correctly
			@retval false   initial leo incorrectly
	    */
		bool   _initSatInfo(string sat_name);
		
		/**
	    * @brief      init the state matrix
	    * @param[in]  neq_size   number of OI funct params
	    * @param[out] state      state data of OI processing
	    * @return bool
			@retval true    initial state correctly
			@retval false   initial state incorrectly
	    */
		bool   _initStateMat(int neq_size, vector<double>& state);
		
		/**
	    * @brief   init the partial matrix
	    * @param[in]  par_num   param number
	    * @param[out] partial   partial matrix of OI
	    * @return bool
			@retval true    initial leo partial correctly
			@retval false   initial leo partial incorrectly
	    */
		bool   _initPartialMat(int par_num, Matrix& partial);

		/**
	    * @brief  calculate the acc and partial in epoch t
	    * @param[in]  t        current epoch for calculation
	    * @param[out] state    state data after calculation
	    * @param[out] partial  partial data after calculation
	    * @return bool
			@retval true     calculate correctly
			@retval false    calculate incorrectly
	    */
		bool _calculate(const t_gtime& t, const vector<double>& state, vector<double>& partial) override;

		/**
	    * @brief  creat the force model for leo
	    * @param[in]  force_ID   the force to creat
	    * @param[in]  force_contain         end time from XML
	    * @return bool
			@retval true      creat force model correctly
			@retval false     creat force model incorrectly
	    */
		bool _creatForceModel(FORCE force_ID, t_gforce_leo* force_contain);
	};
}
#endif // !GOI_H
