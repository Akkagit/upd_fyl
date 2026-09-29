/**
 *
 * @verbatim
	 History
	  -1.0 Jiande Huang 2019-01-14 creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gforcemodel.h
 * @brief		The container of all the force models used in space craft orbit intergration.
 *
 * @author      Jiande Huang, Wuhan University
 * @version		1.0.0
 * @date		2019-01-14
 */

#ifndef GFORCEMODEL_H
#define GFORCEMODEL_H

#include "gforce/gfsolar.h"
#include "gforce/gfgravity.h"
#include "gforce/gfpointmass.h"
#include "gforce/gfrelativity.h"

#include "gdata/gics_nav.h"
#include "gdata/gegmcoef.h"
namespace great
{
	/** @brief The container of all the force models. */
	class t_gforce_model
	{
	public:
		/** @brief default constructor. */
		t_gforce_model() {};

		/** @brief default destructor. */
		virtual ~t_gforce_model() {};

		/** 
		* @brief init force model for force. 
		* @param [in] arg   the init argument for init one force model.
		* @param [in] log   the log file for write log.
		*/
		virtual bool creatForceModel(t_init_arg* arg, t_glog* log) = 0;
		
		/** 
		* @brief clear the force model and delete the point.
		* @param[in]  ptr	 force pointer need to be delete
		*/
		void destoryForceModel(t_gforce *ptr);

	};
}
#endif // !GFORCEMODEL_H
