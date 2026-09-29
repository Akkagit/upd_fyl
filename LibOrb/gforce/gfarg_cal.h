/**
 * This file contains all calculating arguments of diffetent force models used in OI.
 *
 * @verbatim
	 History
	  -1.0 Jiande Huang	 2019-01-14 creat the file.
	  -1.1 Jiande Huang  2019-04-01 Adding abnormal throwing mechanism of functions
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gfarg_cal.h
 * @brief		The argument used in the caculation of force model.
 *
 * @author      Jiande Huang, Wuhan University
 * @version		1.0.0
 * @date		2019-01-14
 */

#ifndef GFARG_CAL_H
#define GFARG_CAL_H

#include <map>
#include "gutils/ExportOI.h"
#include "gutils/gconst_oi.h"
#include "gutils/gtime.h"
#include "gforce/gfcommon.h"

#include "newmat/newmat.h"

namespace great
{
	/** @brief argument used for init force model*/
	enum ARG_TYPE
	{	
		ARG_POINT_MASS,	   ///< N body attraction.
		ARG_SOLAR_ECOM,	   ///< ECOM model for solar radiation.
		ARG_SOLAR_ECOM2,   ///< ECOM2 model for solar radiation.
		ARG_SOLAR_PBW,     ///< Priori Box wing model for solar radiation.	
		ARG_GRAVI_EGM,	   ///< EGM 2008 gravity model.
		ARG_RELATIVITY,	   ///< Default relativity force model.
        ARG_SOLAR_BOXWING, ///< Box wing model for solar radiation.	
		ARG_CUSTOMER,      ///< Default customer force model.
		ARG_ATMOSPHERE,     ///< Default atmosphere force model.
		ARG_EARTHRADIATION,
		ARG_UNDEF
	};

	/** @brief base calss for caculate argument */
	class t_gforce_calarg
	{
	public:
		/** @brief default constructor */
		t_gforce_calarg() { arg_type = ARG_TYPE::ARG_UNDEF; };
		
		/** @brief default destructor */
		virtual ~t_gforce_calarg() {};

		double*			time = nullptr;		///< current epoch, should be TT or TAT.
		ARG_TYPE		arg_type;	        ///< force model type.
		ColumnVector*	sat_pos = nullptr;	///< sat position in J2000.
		ColumnVector*   sat_vel = nullptr;	///< sat velocity in J2000.
	};

	/** @brief point mass argument */
	class t_gfcalarg_point_mass : public t_gforce_calarg
	{
	public:
	    /** @brief default constructor */
		t_gfcalarg_point_mass() { arg_type = ARG_POINT_MASS; };
		
		/** @brief default destructor */
		virtual ~t_gfcalarg_point_mass() {};
	};

	/** @brief ECOM model argument */
	class t_gfcalarg_solar_ecom : public t_gforce_calarg
	{
	public:
		/** @brief default constructor */
		t_gfcalarg_solar_ecom() { arg_type = ARG_SOLAR_ECOM; }
		
		/** @brief default destructor */
		virtual ~t_gfcalarg_solar_ecom() {};

		
		ColumnVector* sun_pos  = nullptr;	///< sun pos in J2000.
		ColumnVector* moon_pos = nullptr;	///< moon pos in J2000.
		double*  shadow_factor = nullptr;	///< shadow factor, which is relative with the positions of satellite, moon and sun.
		string  sat_type;		///< type of satellite, specially for BLOCK IIR.

	};

	/** @brief ECOM2 model argument */
	class t_gfcalarg_solar_ecom2 : public t_gforce_calarg
	{
	public:
		/** @brief default constructor */
		t_gfcalarg_solar_ecom2() { arg_type = ARG_SOLAR_ECOM2; }

		/** @brief default destructor */
		virtual ~t_gfcalarg_solar_ecom2() {};


		ColumnVector* sun_pos = nullptr;	///< sun pos in J2000.
		ColumnVector* moon_pos= nullptr;	///< moon pos in J2000.
		double*  shadow_factor= nullptr;	///< shadow factor, which is relative with the positions of satellite, moon and sun.
		string  sat_type;		///< type of satellite, specially for BLOCK IIR.

	};

	/** @brief BOXWING model argument */
	class t_gfcalarg_solar_boxwing : public t_gforce_calarg
	{
	public:
		/** @brief default constructor */
		t_gfcalarg_solar_boxwing() { arg_type = ARG_SOLAR_BOXWING; }

		/** @brief default destructor */
		virtual ~t_gfcalarg_solar_boxwing() {};

		map<ICS_PAR, double>* solar_par = nullptr;     ///< solar par in ics.
		ColumnVector*         sun_pos   = nullptr;       ///< sun pos in J2000.
		ColumnVector*         moon_pos  = nullptr;	///< moon pos in J2000.
		double*               shadow_factor = nullptr; ///< shadow factor, which is relative with the positions of satellite, moon and sun.
		string  sat_type;		///< type of satellite, specially for BLOCK IIR.
		
	};


	/** @brief PRIORI BOX WING model argument */
	class t_gfcalarg_priori_solar_boxwing : public t_gforce_calarg
	{
	public:
		/** @brief default constructor */
		t_gfcalarg_priori_solar_boxwing() { arg_type = ARG_SOLAR_PBW; }

		/** @brief default destructor */
		virtual ~t_gfcalarg_priori_solar_boxwing() {};

		map<ICS_PAR, double>* solar_par     = nullptr;     ///< solar par in ics.
		ColumnVector*         sun_pos       = nullptr;     ///< sun pos in J2000.
		ColumnVector*         moon_pos      = nullptr;	   ///< moon pos in J2000.
		double*               shadow_factor = nullptr;     ///< shadow factor, which is relative with the positions of satellite, moon and sun.
		string  sat_type;		        ///< type of satellite, specially for BLOCK IIR.
		double  sat_mass  = 0.0;		///< sat mass
	};


	/** @brief customer argument */
	class t_gfcalarg_customer : public t_gforce_calarg
	{
	public:
		/** @brief default constructor */
		t_gfcalarg_customer() { arg_type = ARG_CUSTOMER; }

		/** @brief default destructor */
		virtual ~t_gfcalarg_customer() {};
		double*              period          = nullptr;         ///< period of sat.
		Matrix*              scf2crs         = nullptr;         ///< rotmat of scf2crs.
		std::map<ICS_PAR, double>* cust_pars = nullptr;         ///< customer par in ics.
	};

	/** @brief atmosphere argument */
	class t_gfcalarg_atmosphere : public t_gforce_calarg
	{
	public:
		/** @brief default constructor */
		t_gfcalarg_atmosphere() { arg_type = ARG_ATMOSPHERE; }

		/** @brief default destructor */
		virtual ~t_gfcalarg_atmosphere() {};
		ColumnVector*	     sat_pos_trs = nullptr;  ///< sat pos in TRS.
		ColumnVector*        sat_vel_trs = nullptr;  ///< sat vel in TRS.
		ColumnVector*        sun_pos     = nullptr;      ///< sun pos in J2000.
		Matrix*              e2j         = nullptr;          ///< rotmat of trs2crs.
		double*              atmos_par   = nullptr;    ///< atmosphere par in ics.
		int                  year = 0;
		int                   mon = 0;
		int                   doy = 0;
		double               xics_x = 0.0;
		double               xics_y = 0.0;
		double               xics_z = 0.0;
	};

	class t_gfcalarg_earthradiation : public t_gforce_calarg
	{
	public:
		/** @brief default constructor */
		t_gfcalarg_earthradiation() { arg_type = ARG_EARTHRADIATION; }

		/** @brief default destructor */
		virtual ~t_gfcalarg_earthradiation() {};
//		ColumnVector*	     sat_pos_crs;  ///< sat pos in CRS.
//		ColumnVector*        sat_vel_crs;  ///< sat vel in CRS.
		ColumnVector*        sun_pos_crs = nullptr;      ///< sun pos in J2000.
		Matrix*              e2j = nullptr;          ///< rotmat of trs2crs.
		double sat_mass = 0.0;
		string sat_blocknum;
	};

	/** @brief gravity argument */
	class t_gfcalarg_gravity_egm : public t_gforce_calarg
	{
	public:
		/** @brief default constructor */
		t_gfcalarg_gravity_egm() { arg_type = ARG_GRAVI_EGM; };
		
		/** @brief default destructor */
		virtual ~t_gfcalarg_gravity_egm() {};

		
		ColumnVector *sat_pos_trs = nullptr;			///< sat pos in TRS.
		map<string, ColumnVector> *allXE = nullptr;	///< planet positions in TRS.
		double* xpole = nullptr;						///< X coordinates for pole point. 
		double* ypole = nullptr;						///< Y coordinates for pole point. 
		double* gst   = nullptr;						///< Greenwich Time.
		Matrix* e2j = nullptr;						///< Rot matrix from TRS to J2000.
	};

	
	/** @brief relativity argument */
	class t_gfcalarg_relativity : public t_gforce_calarg
	{
	public:
		/** @brief default constructor */
		t_gfcalarg_relativity() { arg_type = ARG_RELATIVITY; };
		
		/** @brief default destructor */
		virtual ~t_gfcalarg_relativity() {};

		
		double  gm = EARTH_GM;				///< GM for earth.
	};
}
#endif // !GFARGUMENT_H
