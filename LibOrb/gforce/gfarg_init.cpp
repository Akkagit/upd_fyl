/**
 *
 * @verbatim
	 History
	  -1.0	Jiande Huang	  2019-01-11	 creat the file.
	  -1.1  Jie Li	      2019-04-01	 add code in doxygen style.
	  -1.2  Yujie Qin     2019-04-15 change sth about leo
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		    gfarg_init.cpp
 * @brief		    The base class used to initialize the force model parameters。
*
 * @author          Jiande Huang, Wuhan University
 * @version		    1.0.0
 * @date		    2019-01-06
 */

#include "gforce/gfarg_init.h"

namespace great
{
	/**
	* @brief constructor, initialize gravity EGM2008 model.
	* @param[in]   model      force model enum
	* @param[in]   order      EGM model order
	* @param[in]   egmcoef    EGM coefficient data from EGM file
	* @param[in]   navde      jpleph_de405 data
	* @param[in]   otl        oceanload data
	* @param[in]   ot         oceantide data
	*/
	t_init_gravity_egm08::t_init_gravity_egm08(FORCE_MODEL model, int order, t_gegmcoef * egmcoef,
		t_gnavde * navde, t_gallotl * otl, t_gocean_tide * ot) : t_init_arg(model)
	{
		gravity_order = order;		///< gravity order.
		gracoef = egmcoef;	///< EGM-file data.
		gplanet = navde;	///< DE行星星历.
		gotl = otl;		    ///< ocean load data.
		gotdata = ot;	    ///< ocean tide data.
	}

	t_init_gravity_egm08::t_init_gravity_egm08(FORCE_MODEL model, int order, t_gegmcoef * egmcoef,
		t_gnavde * navde, t_gallotl * otl, t_gocean_tide * ot, t_gdeadata* de, int oceandegree, bool poletide,int poledegree) : t_init_arg(model)
	{
		gravity_order = order;		///< gravity order.
		gracoef = egmcoef;	///< EGM-file data.
		gplanet = navde;	///< DE行星星历.
		gotl = otl;		    ///< ocean load data.
		gotdata = ot;	    ///< ocean tide data.
		dedata = de;
		ocean_degree = oceandegree;
		ispoletide = poletide;
		pole_degree = poledegree;
	}
	/**
	 * @brief constructor, base initialize force model.
	 * @param[in]   model      force model enum
	 */
	t_init_arg::t_init_arg(FORCE_MODEL model)
	{
		force_model = model;
	}

	/**
	 * @brief constructor, initialize relativity model.
	 * @param[in]   model      force model enum
	 */
	t_init_relativity::t_init_relativity(FORCE_MODEL model) : t_init_arg(model)
	{
	}



	/**
	 * @brief constructor, initialize point_mass(N-body) model.
	 * @param[in]   model      force model enum
	 * @param[in]   navde	   JPL data
	 */
	t_init_point_mass::t_init_point_mass(FORCE_MODEL model,
		t_gnavde * navde) : t_init_arg(model)
	{
		gplanet = navde;
	}

	/**
	 * @brief constructor, initialize customer model.
	 * @param[in]   model      force model enum
	 * @param[in]   pars_map   parameters list of ecom model
	 */
	t_init_solar_ecom::t_init_solar_ecom(FORCE_MODEL model,
		std::map<ICS_PAR, double> pars_map) : t_init_arg(model)
	{
		solar_pars = pars_map;
	}

	/**
	* @brief constructor, initialize customer model.
	* @param[in]   model      force model enum
	* @param[in]   pars_map   parameters list of ecom model
	*/
	t_init_solar_ecom2::t_init_solar_ecom2(FORCE_MODEL model,
		std::map<ICS_PAR, double> pars_map) : t_init_arg(model)
	{
		solar_pars = pars_map;
	}

	/**
	 * @brief constructor, initialize customer model.
	 * @param[in]   model      force model enum
	 * @param[in]   partial    calculate the partial or not
	 * @param[in]   panneldata pannel data of leo
	 * @param[in]   massage	   mass of leo

	 */
	t_init_solar_boxwing::t_init_solar_boxwing(FORCE_MODEL model, bool partial, t_gpannel* panneldata, double massage) : t_init_arg(model)
	{
		ispartial = partial;
		pannel = panneldata;
		mass = massage;
	}

	/**
	* @brief constructor, initialize customer model.
	* @param[in]   model      force model enum
	* @param[in]   partial    calculate the partial or not
	* @param[in]   panneldata pannel data of leo
	* @param[in]   massage	   mass of leo

	*/
	t_init_solar_priori_boxwing::t_init_solar_priori_boxwing(FORCE_MODEL model, t_gpannel* panneldata, std::map<ICS_PAR, double> pars_map) : t_init_arg(model)
	{
		//ispartial = partial;
		pannel = panneldata;
		solar_pars = pars_map;
	}

	/**
	 * @brief constructor, initialize customer model.
	 * @param[in]   model      force model enum
	 * @param[in]   mass_	   mass of leo
	 * @param[in]   partial_   calculate the partial or not
	 * @param[in]   period_	   period of leo
	 */
	t_init_customer::t_init_customer(FORCE_MODEL model, double mass_, bool partial_, double period_) : t_init_arg(model)
	{
		mass = mass_;          ///< mass of sat.
		ispartial = partial_;  ///< calculate the partial or not.
		period = period_;      ///< period of sat.
	}

	/**
	 * @brief constructor, initialize customer model.
	 * @param[in]   model        force model enum
	 * @param[in]   massage	     mass of leo
	 * @param[in]   atm_model    model of atmosphere
	 * @param[in]   gsolmagdata	 solar data of leo
	 * @param[in]   gpanneldata	 pannel data of leo
	 */
	t_init_atmosphere::t_init_atmosphere(FORCE_MODEL model_, double massage, string atm_model, t_gsolargeomag*  gsolmagdata, t_gpannel*  gpanneldata) : t_init_arg(model_)
	{
		mass = massage;
		model = atm_model;
		gsolmag = gsolmagdata;
		gpannel = gpanneldata;
	}

	t_init_earthradiation::t_init_earthradiation(FORCE_MODEL model_, double massage, string ear_model, t_galbedo*  galbedodata, t_gpannel*  gpanneldata) : t_init_arg(model_)
	{
		mass = massage;
		model = ear_model;
		galbedo = galbedodata;
		gpannel = gpanneldata;
	}

	t_init_earthradiation2::t_init_earthradiation2(FORCE_MODEL model_, string ear_model, t_galbedo* galbedodata) : t_init_arg(model_)
	{
		model = ear_model;
		galbedo = galbedodata;
	}
}



