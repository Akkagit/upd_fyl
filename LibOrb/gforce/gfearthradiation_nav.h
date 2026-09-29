#ifndef GFEARTHRADIATION2_H
#define GFEARTHRADIATION2_H

#include "gutils/gconst.h"
#include "gforce/gforce.h"
#include "gforce/gfarg_init.h"
#include "gio/glog.h"
#include "gutils/gconst_oi.h"
/**
*
* @verbatim
History
	-1.0     Liu Chengbo     2020-03-03 Creat the file.
	-1.1     Jiande Huang    2020-10-01 Standard code style and Optimize code logic
@endverbatim
* Copyright (c) 2018, Wuhan University.All rights reserved.
*
* @file		      gfearthradiation2.h
* @brief		    The class for earth radiation model.
*
* @author          Liu Chengbo, Wuhan University
* @version		    2.0.0
* @date		    2020-03-03
*
*/

#include "gutils/ExportOI.h"
#include "gdata/galbedo.h"
#include "newmat/newmat.h"

using namespace std;

namespace great
{
	/** @brief class for earth radiation model. */
	class OI_LIBRARY_EXPORT t_gfearthradiation_nav : public t_gforce
	{
	public:

		/** @brief constructor */
		t_gfearthradiation_nav(t_init_arg* arg, t_glog* log);

		/** @brief destructor */
		virtual ~t_gfearthradiation_nav();

		/**
		* @brief value latkmx/lonkmx.
		* @param[in] grd the grid
		* @return
		@retval true   calculate normally
		@retval false  have some error
		*/
		bool cal_max(int grd);

		bool cal_sat_panel(string block);

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
		double toa;
		double pi;
		Matrix v_ns_all1;       ///< the earth properties.
		Matrix v_ns_all2;       ///< the earth properties.
		Matrix v_ns_all3;       ///< the earth properties.
		Matrix cergre;          ///< the reflection data.
		Matrix cergem;          ///< the emission data.
		Matrix d_area_all;      ///< the earth properties.
		t_galbedo* _galbedo;    ///< the ceres data.
		double     _mass;       ///< the mass of sat.
		string     _model;      ///< the model of earth radiation
		ColumnVector ysat;      ///< satellite position (unit:m) and velocity (unit:m/s) in crs
		ColumnVector sun;       ///< sun position vector (unit:m)
		ColumnVector vsat;      ///< satellite velocity from ics
		ColumnVector Xatt;      ///< nominal satellite attitude
		ColumnVector Yatt;      ///< nominal satellite attitude
		ColumnVector Zatt;      ///< nominal satellite attitude
		ColumnVector Datt;      ///< nominal satellite attitude
		string blockname;       ///< sat block name
		int latkmx;
		int lonkmx;
		double grdnum;          ///< the number of the grid
		int grdcer;
		double area[4][2];
		double difu[4][2];
		double refl[4][2];
		double area2[4][2];
		double refl2[4][2];
		double difu2[4][2];
		double difuir[4][2];
		double reflir[4][2];
		double difuir2[4][2];
		double reflir2[4][2];

	};
}

#endif