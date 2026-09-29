/**
 *
 * @verbatim
     History
      -1.0 Yujie Qin and  Bo Wong  2019-01-06 Creat the file.
      -1.1 Hongmin Zhang           2019-03-04 update.
      -1.2 Yujie Qin               2019-03-05 edit for final version.
      -1.3 Yujie Qin               2019-03-15 Code comments are added in Doxygen format.
	  -1.4 Yujie Qin               2019-07-03 Adding Doxygen Style Code Remarks
   @endverbatim
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		    gfatmosphere.h
 * @brief		    The class for atmospheric drag mode.
*
 * @author          Yujie Qin, Wuhan University
 * @version		    1.0.0
 * @date		    2019-01-06
 *
 */

#ifndef GFATMOSPHERE_H
#define GFATMOSPHERE_H

#include "gutils/gconst.h"
#include "gforce/gforce.h"
#include "gforce/gfarg_init.h"
#include "gio/glog.h"
#include "gutils/gconst_oi.h"
#include "gutils/ExportOI.h"
#include "gforce/gdtm94.h"
#include "gforce/gnrlmsise00.h"
#include "gforce/gdtm2013.h"
#include "newmat/newmat.h"


#define COEFFICIENT_DRAG 2.2   ///< the coefficient of drag(Cd)

using namespace std;

namespace great
{
	/** @brief class for atmospheric drag model. */
	class OI_LIBRARY_EXPORT t_gfatmosphere : public t_gforce
	{
	public:
		
		/** @brief constructor */
		t_gfatmosphere(t_init_arg* arg, t_glog* log);

		/** @brief destructor */
		virtual ~t_gfatmosphere();

		/** 
		* @brief calculate the accleration and the corresponding partial matrix.
	    * @param[in] argument Parameters to be used in calculation
		* @return
			@retval true   calculate normally
			@retval false  have some error
	    */
		bool calculate(t_gforce_calarg* argument);


	protected:

		ColumnVector               _sat_pos;         ///< the position of the satellite in the TRS.
		ColumnVector               _sat_vel;         ///< the velocity of the satellite in the TRS.
		ColumnVector               _sun_pos;         ///< the position of the sun in the TRS.
		double                     _epoch;           ///< the current gps time.
		string                     _model;           ///< the model of calculation of density
		t_gsolargeomag*            _gsolmag;         ///< the data of kp and solar.
		t_gpannel*                 _gpannel;         ///< the pannel data of sat.
		double                     _mass;            ///< the mass of sat.
		double                     _drag;            ///< the parameter of atmos.
		double                     _kp;              ///< the kp index of geomagnetism for DTM94 density computation.
		vector<double>             _ap;              ///< the ap index of geomagnetism for MSISE00 density computation.
		double                     _akp[4];             ///< the akp index of geomagnetism for DTM2013 density computation.
		vector<double>             _solar;           ///< the solar flux value for density computation.
		double                     _density;         ///< dendity of the atmosphere.

		/**
		* @brief calculate kp and solar of reference time.
	    * @param[in] time the current gps_t
		* @return
			@retval true   calculate normally
			@retval false  have some error
	    */
		bool _calculate_sol_geo(double time);

		/** 
		* @brief liner interpolation.
	    * @param[in] start_value the start value
	    * @param[in] end_value the end value
	    * @param[in] start_time the start time
	    * @param[in] end_time the end time
	    * @param[in] ref_time the reference time
		* @return the value interpolated
	    */
		double _linear_interpolation(double start_value, double end_value, double start_time, double end_time, double ref_time);

		/**
		* @brief liner interpolation.
		* @param[in] ap_value the ap value of 4 days
		* @param[in] leftvalue the index of ap value on current epoch in ap_value
		*/
		double _linear_interpolation(vector<double> ap_value, int leftvalue);

		/**
		* @brief liner interpolation.
		* @param[in] start_value the start value
		* @param[in] end_value the end value
		* @param[in] start_time the start time
		* @param[in] end_time the end time
		* @param[in] ref_time the reference time
		* @param[in] kp_value the kp value of 24h from current epoch on 
		*/
		double _linear_interpolation(double start_value, double end_value, double start_time, double end_time, double ref_time, vector<double> kp_value);

	};
} //namespace

#endif // GATMOSPHERE_H
