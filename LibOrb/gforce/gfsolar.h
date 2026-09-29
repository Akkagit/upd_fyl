/**
 * 
 * @verbatim
	 History
	  -1.0 Jiande Huang 2019-01-14 creat the file.
   @endverbatim
 *
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gfsolar.h
 * @brief		The base class for all the solar force models.
 *
 * @author      Jiande Huang, Wuhan University
 * @version		1.0.0
 * @date		2019-01-14
 */

#ifndef GFSOLAR_H
#define GFSOLAR_H

#include "gforce/gforce.h"
#include <map>

using namespace std;
namespace great
{
	/** The base class for the solar force model.*/
	class OI_LIBRARY_EXPORT t_gfsolar : public t_gforce
	{
	public:
		/**
		 * @brief constructor.
		 * @param[in]  log        the log file for write log.
		 */
		t_gfsolar(t_glog* log);

		/** @brief default destructor. */
		virtual ~t_gfsolar();

		/** 
		 * @brief Show the detail information of the force model.
		 */
		virtual void show() = 0; 

		/** 
		 * @brief Check the information of the force model 
		 */
		virtual bool check() = 0;

		/** An empty function for the subclass to inherit. 
		* @param [in] argument The argument used in the calculation.
		*/
		//bool calculate(t_gforce_calarg* argument);
	};
}
#endif // !GSOLAR_H
