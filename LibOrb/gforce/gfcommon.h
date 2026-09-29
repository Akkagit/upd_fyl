/**
 * This file contains some common defination and functions for OI,
 *		include translation bewteen different data formate.
 *
 * @verbatim
	 History
	  -1.0 Jiande Huang  2019-01-06 creat the file.
	  -1.1 Jie Li	  	 2019-03-31 add code in doxygen style
	  -1.2 Yujie Qin     2019-07-04 Adding Doxygen Style Code Remarks
	  -1.3 Jiande Huang  2020-10-01 Standard code style and Optimize code logic
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gfcommon.h
 * @brief		Some common defination and functions for orbit integration process.
 *
 * @author      Jiande Huang, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 *
 */

#ifndef GFCOMMON_H
#define GFCOMMON_H

#include "newmat/newmat.h"
#include "gutils/gtime.h"
#include "gutils/gnss.h"
#include "gutils/ExportOI.h"
#include <math.h>
using namespace gnut;

namespace great
{	
	/** @brief Force considered */
	enum FORCE
	{
		GRAVITY = 0,     ///< gravity.
		RELATIVITY = 1,  ///< relativity.
		POINT_MASS = 2,  ///< point mass(N body).
		SOLAR_RADI = 3,  ///< solar radiation.
		EARTH_RADI = 4,  ///< earth radiation.
		CUSTOMER = 5,    ///< customer.
		ATMOSPHERE = 6   ///< atmosphere.
	};

	/** @brief Force model considered */
	enum FORCE_MODEL
	{
		F_GRAVI_DEF,      ///< default gravity model (EGM08 12).
		F_RELAV_DEF,      ///< default relativity model.
		F_SOLAR_DEF,      ///< default solar radiation model(ECOM 5).
		F_NBODY_DEF,      ///< default point mass model.

		F_GRAVI_EGM08,    ///< EGM 2008.
		F_SOLAR_ECOM_3,   ///< ECOM with 3 parameters(only have constant term)
        F_SOLAR_ECOM_5,	  ///< ECOM with 5 parameters.
		F_SOLAR_ECOM_9,   ///< ECOM with 9 parameters.
		F_SOLAR_ECOM2,    ///< ECOM2 with 9 parameters.
		F_SOLAR_PBW_5,    ///< The a Priori Box-wing model with ECOM-5 parameters.
		F_SOLAR_PBW_9,    ///< The a Priori Box-wing model with ECOM-9 parameters.
        F_SOLAR_BOXWING,  ///< box wing.
		F_CUSTOMER,       ///< customer.
		F_ATMOSPHERE ,    ///< atmosphere.
		F_EARTHRADIATION, ///< earth radiation
		
		NO_DEF_F_MODEL
	};

	/** @brief all kinds of ics par */
	enum class ICS_PAR
	{
		PXSAT,		///< coordinate of satellite.
		PYSAT,      ///< coordinate of satellite.
		PZSAT,      ///< coordinate of satellite.
		VXSAT,		///< velocity of satellite.
		VYSAT,      ///< velocity of satellite.
		VZSAT,      ///< velocity of satellite.

		ECOM_D0,	///< parameters of ECOM model.
		ECOM_Dc,
		ECOM_Ds,
		ECOM_Y0,
		ECOM_Yc,
		ECOM_Ys,
		ECOM_X0,
		ECOM_Xc,
		ECOM_Xs,

		ECOM2_D0,	///< parameters of ECOM2 model.
		ECOM2_D2c,
		ECOM2_D2s,
		ECOM2_D4c,
		ECOM2_D4s,
		ECOM2_Y0,
		ECOM2_X0,
		ECOM2_X1c,
		ECOM2_X1s,

		SR_scale,   ///< parameters of solar radiation pressure scale.
		DRAG_c,     ///< parameters of atmosphere drag.

		EMP_Ac,     ///< parameters customer.
		EMP_Bc,
		EMP_Aa,
		EMP_Ba,
		EMP_Ar,
		EMP_Br,

		EMP_Sc2, 
		EMP_Cc2,
		EMP_Sa2, 
		EMP_Ca2,
		EMP_Sr2, 
		EMP_Cr2,

		EMP_Sc1,
		EMP_Cc1,
		EMP_Sa1,
		EMP_Ca1,
		EMP_Sr1,
		EMP_Cr1,

		NO_DEF
	};

	enum CHANGE_TYPE
	{
		GPS2TDT,
		TAI2TDT
	};

	/**
	* @brief calculate the shadow factor.
	* @param[in]   sat_pos      position of satellite
	* @param[in]   sun_pos		position of sun
	* @param[in]   moon_pos		position of moon
	* @return		shadow factor
	*/
	double get_shadow_factor(const ColumnVector& sat_pos,const ColumnVector& sun_pos,const ColumnVector& moon_pos);

	/**
    * @brief calculate the shadow area.
    * @param[in]   a_radius     radius of object a
    * @param[in]   b_radius		radius of object b
    * @param[in]   ab_distance	distance between a and b
    * @return		shadow area
    */
	double get_shadow_area(double a_radius, double b_radius, double ab_distance);

	/**
     * @brief change time system.
     * @param[in]   type			time change type(enum CHANGE_TYPE)
     * @param[in]   inp_time		time need to be changed(mjd)
     * @param[out]  out_time		time have been changed(mjd)
     */
	void change_time_sys(CHANGE_TYPE type ,double inp_time, double& out_time);

	/**
    * @brief change time system.
    * @param[in]   type			time change type(enum CHANGE_TYPE)
    * @param[in]   inp_time		time need to be changed
    * @param[out]  out_time		time have been changed
    */
	void OI_LIBRARY_EXPORT change_time_sys(CHANGE_TYPE type, t_gtime inp_time, t_gtime& out_time);

	/**
    * @brief string type to ICS_PAR
    * @param[in]   par_name		name of parameters
    * @return	enum ICS_PAR
    */
	ICS_PAR OI_LIBRARY_EXPORT str2ics(string par_name);

	/**
	* @brief string type to ICS_PAR for leos
	* @param[in]   par_name		name of parameters
	* @return	enum ICS_PAR
	*/
	ICS_PAR OI_LIBRARY_EXPORT str2icsleo(string par_name);
	
	/**
	* @brief FORCE type to string
	* @param[in]   force	force model name in FORCE
	* @return	force model name in string type
	*/
	string  OI_LIBRARY_EXPORT force2str(FORCE force);
	
	/**
	* @brief ICS_PAR type to string
	* @param[in]   par_name		parameter name in ICS_PAR
	* @return	parameter namee in string type
	*/
	string  OI_LIBRARY_EXPORT ics2str(ICS_PAR par_name);
	
	/**
	* @brief string type to bool
	* @param[in]   s	bool expressed by string
	* @return	bool
	*/
	bool    str2bool(const string& s);
	
	/**
	* @brief GSYS to string type
	* @param[in]   sys	system expressed by GSYS
	* @return	system expressed by string
	*/
	string  OI_LIBRARY_EXPORT gsys2str(GSYS sys);
	
	/**
	* @brief string type to PLANET
	* @param[in]   planet	planet name expressed by string
	* @return	planet name expressed by PLANET(enum)
	*/
	int     str2planet(const string& planet);

    /**
	 * @brief calculate mean pole by time
	 * @param[in]   mjd		mjd time
	 * @param[out]  xpm		mean xpole
	 * @param[out]  ypm		mean ypole
	 */
	void OI_LIBRARY_EXPORT  get_mean_pole(double mjd, double& xpm, double& ypm);
}

#endif 
// !GFCOMMON_H
