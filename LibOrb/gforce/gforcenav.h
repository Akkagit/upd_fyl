/**
 *
 * @verbatim
	 History
	  -1.0 Jiande Huang  2019-02-25 creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gforcenav.h
 * @brief		A container of all the force models for navigation satellite. 
 *              We consider solar radiation, gravity and so on.\n
 *
 * @author      Jiande Huang, Wuhan University
 * @version		1.0.0
 * @date		2019-03-01
 *
 */

#ifndef GFORCENAV_H
#define GFORCENAV_H

#include "gforce/gforcemodel.h"
#include "gutils/ExportOI.h"
#include "gutils/gconst_oi.h"

namespace great
{
	/**
	* @brief	Class for all the force model used in GNS satellites orbit integral.
	*			Currently, we only consider gravity, relativity, N body attraction, solar radiation.
	*/
	class OI_LIBRARY_EXPORT t_gforce_nav : public t_gforce_model
	{
	public:
		/** 
		 * @brief default constructor 
		 */
		t_gforce_nav() {};
		/** 
		 * @brief default destructor, all the point used in the class will be deleted 
		 */
		virtual ~t_gforce_nav();

		/** 
		 * @brief init force model for force, for example : ECOM 5 for solar radition.
		 * @param [in] arg   the init argument for init one force model.
		 * @param [in] log   the log file for write log.
		 */
		bool creatForceModel(t_init_arg* arg, t_glog* log);

		t_gforce*	point_mass = NULL;	///< the point for N body force.
		t_gforce*	solar_radi = NULL;	///< the point for solar radiation.
		t_gforce*	gravity    = NULL;    ///< the point for gravity.
		t_gforce*	relativity = NULL; ///< the point for relativity.
		t_gforce*   earthradi = NULL;  ///<the point for earthradiation.>
	};
}

#endif // !GFORCENAV_H
