/**
 *
 *	Copyright (c) 2018, Wuhan University
 *	All rights reserved.
 *
 * @file		gforcenav.cpp
 * @brief		The functions in gforcenav.h are implemented here.\n
 *
 * @author      Jiande Huang, Wuhan University
 * @version		1.0.0
 * @date		2019-03-01
 *
 * @verbatim
	 History
	  -1.0 Jiande Huang 2019-02-25 creat the file.
   @endverbatim
 *
 */

#include "gforce/gforcenav.h"
#include "gforce/gfsolarecom.h"
#include "gforce/gfsolarecom2.h"
#include "gforce/gfsolarprioriboxwing.h"
#include "gforce/gfgravityegm.h"
#include "gforce/gfearthradiation_nav.h"

namespace great
{
	/** default destructor, all the point used in the class will be deleted */
	t_gforce_nav::~t_gforce_nav()
	{
		if (point_mass)
		{
			delete	point_mass; 
			point_mass = NULL;
		}
		if (solar_radi)
		{
			delete	solar_radi;
			solar_radi = NULL;
		}
		if (gravity)
		{
			delete	gravity;
			gravity = NULL;
		}
		if (relativity)
		{
			delete	relativity; 
			relativity = NULL;
		}
		if (earthradi)
		{
			delete earthradi;
			earthradi = NULL;
		}
	}

	/**
	 *
	 @verbatim
	  Currently,we consider solar radiation, gravity, point mass and relativity.
	  - solar radiation : ECOM with 5 parameters, ECOM with 9 parameters.
	  - gravity         : 12th-order EGM08.
	  - point mass      : defalut model.
	  - relativity      : defalut model.
	 @endverbatim
	 *
	 * @param [in] arg : the init argument for init one force model.
	 * @param [in] log : the log file for write log.
	 * -# t_init_solar_ecom    : init ecom force model.
	 * -# t_init_pointmass     : init pointmass force model.
	 * -# t_init_gravity_egm08 : init EGM08 force model.
	 * -# t_init_relativity    : init relativity force model.
	 *
	 * @note t_init_arg is a base class for the force init argument.
	 */
	bool t_gforce_nav::creatForceModel(t_init_arg * arg, t_glog* log)
	{
		switch (arg->force_model)
		{
		    // default solar radiation model is ECOM 5.
		    case F_SOLAR_DEF:
		    {
		    	solar_radi = new t_gfsolar_ecom(arg, log);
		    	break;
		    }
		    // ECOM with 5 parameters 
		    case F_SOLAR_ECOM_5:
		    {
		    	solar_radi = new t_gfsolar_ecom(arg, log);
		    	break;
		    }
		    // ECOM with 9 parameters
		    case F_SOLAR_ECOM_9:
		    {
		    	solar_radi = new t_gfsolar_ecom(arg, log);
		    	break;
		    }
			// ECOM2 with 9 parameters
			case F_SOLAR_ECOM2:
			{
				solar_radi = new t_gfsolar_ecom2(arg, log);
				break;
			}
			// PBW with 5 ECOM parameters
			case F_SOLAR_PBW_5:
			{
				solar_radi = new t_gfsolar_pbw(arg, log);
				break;
			}
			// PBW with 5 ECOM parameters
			case F_SOLAR_PBW_9:
			{
				solar_radi = new t_gfsolar_pbw(arg, log);
				break;
			}
		    // default relativity model
		    case F_RELAV_DEF:
		    {
		    	relativity = new t_gfrelativity(arg, log);
		    	break;
		    }
		    // default N body attraction model
		    case F_NBODY_DEF:
		    {
		    	point_mass = new t_gfpoint_mass(arg, log);
		    	break;
		    }
		    // default gravity model is 12th-order EGM08
		    case F_GRAVI_DEF:
		    {
		    	gravity = new t_gfgravity_egm08(arg, log);
		    	break;
		    }
		    // 12th-order EGM08 model
		    case F_GRAVI_EGM08:
		    {
		    	gravity = new t_gfgravity_egm08(arg, log);
		    	break;
		    }
			case F_EARTHRADIATION:
			{
				earthradi = new t_gfearthradiation_nav(arg, log);
				break;
			}
		    default:
		    {
		    	return false;
		    	break;
		    }
		}
		return true;
	}
}

