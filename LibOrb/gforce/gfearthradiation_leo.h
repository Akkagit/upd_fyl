#ifndef GFEARTHRADIATION_H
#define GFEARTHRADIATION_H

#include "gutils/gconst.h"
#include "gforce/gforce.h"
#include "gforce/gfarg_init.h"
#include "gio/glog.h"
#include "gutils/gconst_oi.h"
/**
*
* @verbatim
History
	-1.0 Qian Zhang      2019-11-10 Creat the file.
	-1.1 Jiande Huang    2020-10-01 Standard code style and Optimize code logic
@endverbatim
* Copyright (c) 2018, Wuhan University.All rights reserved.
*
* @file		      gfearthradiation.h
* @brief		    The class for earth radiation model.
*
* @author          Qian Zhang, Wuhan University
* @version		    1.0.0
* @date		    2019-11-10
*
*/

#include "gutils/ExportOI.h"
#include "gdata/galbedo.h"
#include "newmat/newmat.h"

using namespace std;

namespace great
{
	/** @brief class for earth radiation model. */
	class OI_LIBRARY_EXPORT t_gfearthradiation_leo : public t_gforce
	{
		public:

			/** @brief constructor */
			t_gfearthradiation_leo(t_init_arg* arg, t_glog * log);

			/** @brief destructor */
			virtual ~t_gfearthradiation_leo();

			/**
			* @brief value latkmx/lonkmx.
			* @param[in] grd the grid
			* @return
			@retval true   calculate normally
			@retval false  have some error
			*/
			bool cal_max(int grd);

			/**
			* @brief calculate the reflection/emission properties value.
			* @param[in] grd          the grid
			* @param[in] ceresdata_re the reflection properties value
			* @param[in] ceresdata_em the emission properties value
			* @return
			@retval true   calculate normally
			@retval false  have some error
			*/
			bool cal_cerg_re_em(int grd, Matrix& ceresdata_re, Matrix& ceresdata_em);

			/**
			* @brief calculate the accleration caused by earth radiation.
			* @param[in] argument Parameters to be used in calculation
			* @return
			@retval true   calculate normally
			@retval false  have some error
			*/
			bool calculate(t_gforce_calarg* argument);

		protected:
			string     _model;      ///< the model of earth radiation

			int _latkmx = 0;
			int _lonkmx = 0;
			int _grdcer = 0;

			double  _toa = 0.0;
			double  _pi  = 0.0;
			double  _sat_mas = 0.0;       ///< the mass of sat.
			double  _grd_num= 0.0;        ///< the number of the grid

			ColumnVector sat_crs;      ///< satellite position (unit:m) and velocity (unit:m/s) in crs
			ColumnVector sun_crs;       ///< sun position vector (unit:m)

			Matrix _v_ns_all1;       ///< the earth properties.
			Matrix _v_ns_all2;       ///< the earth properties.
			Matrix _v_ns_all3;       ///< the earth properties.
			Matrix _cergre;          ///< the reflection data.
			Matrix _cergem;          ///< the emission data.
			Matrix _d_area_all;      ///< the earth properties.

			t_gpannel*  _gpannel = nullptr;    ///< the pannel data of sat.
			t_galbedo*  _galbedo = nullptr;    ///< the ceres data.

	};
}

#endif