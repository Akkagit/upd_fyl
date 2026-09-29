/**
 *
 * @verbatim
     History
      -1.0 Yujie Qin 2019-03-05 creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		    gforceleo.h
 * @brief		    The base class for all the force-leo models.
*
 * @author          Yujie Qin, Wuhan University
 * @version		    1.0.0
 * @date		    2019-03-05
 */
#include "gforceleo.h"
#include "gforce/gfsolarboxwing.h"
#include "gforce/gfcustomer.h"
#include "gforce/gfgravityegm.h"
#include "gforce/gfatmosphere.h"
#include "gforce/gfearthradiation_leo.h"

namespace great
{
	/** default destructor, all the point used in the class will be deleted */
	t_gforce_leo::~t_gforce_leo()
	{
		delete	point_mass; point_mass = NULL;
		delete	solar_radi; solar_radi = NULL;
		delete	gravity;	gravity = NULL;
		delete	relativity; relativity = NULL;
		delete	atmosphere; atmosphere = NULL;
		/*if (point_mass){ delete	point_mass; point_mass = NULL; }
		if (solar_radi) { delete	solar_radi; solar_radi = NULL; }
		if (gravity) { delete	gravity;	gravity = NULL; }
		if (relativity) { delete	relativity; relativity = NULL; }
		if (atmosphere) { delete	atmosphere; atmosphere = NULL; }
		delete	earthradi; earthradi = NULL;
		if (customer) { delete	customer;	customer = NULL; }*/
		delete	customer;	customer = NULL;
	}
	/**
	*
	@verbatim
	 Currently,we consider gravity, N body attraction, relativity, solar radiation, atmosphere and customer.
	- gravity         : 120th-order EGM08.
	- point mass      : defalut model.
	- relativity      : defalut model.
	- solar radiation : BOXWING model.
	- atmosphere      : dtm94 model for density, one parameter.
	- customer        : six parameters.
	@endverbatim
	*
	* @param [in] arg : the init argument for init one force model.
	
	* -# t_init_gravity_egm08 : init EGM08 force model.
	* -# t_init_pointmass     : init pointmass force model.
	* -# t_init_relativity    : init relativity force model.	
	* -# t_init_solar_boxwing : init box wing force model.
	* -# t_init_atmosphere    : init atmosphere force model.
	* -# t_init_customer      : init customer force model.
	*
	* @note t_init_arg is a base class for the force init argument.
	*/
	bool t_gforce_leo::creatForceModel(t_init_arg * arg, t_glog* log)
	{
		switch (arg->force_model)
		{
		case F_RELAV_DEF:
			relativity = new t_gfrelativity(arg, log);
			break;
		case F_NBODY_DEF:
			point_mass = new t_gfpoint_mass(arg, log);
			break;
		case F_GRAVI_EGM08:
			gravity = new t_gfgravity_egm08(arg, log);
			break;
		case F_SOLAR_BOXWING:
			solar_radi = new t_gfsolarbow(arg, log);
			break;
		case F_CUSTOMER:
			customer = new t_gfcustomer(arg, log);
			break;
		case F_ATMOSPHERE:
			atmosphere = new t_gfatmosphere(arg, log);
			break;
		case F_EARTHRADIATION:
			earthradi = new t_gfearthradiation_leo(arg, log);
			break;
		default:
			break;
		}
		return true;
	}
}

