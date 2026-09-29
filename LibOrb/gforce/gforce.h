/**
 *
 * @verbatim
	 History
	  -1.0 Jiande Huang  2019-01-14 creat the file.
	  -1.1 Jiande Huang  2019-03-29 Add the value log and the function write_log.
      -1.2 Yujie Qin    2019-07-04 Adding Doxygen Style Code Remarks
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gforce.h
 * @brief		The base class for all the force models.
 *
 * @author      Jiande Huang, Wuhan University
 * @version		1.0.0
 * @date		2019-01-14
 */

#ifndef GFORCE_H
#define GFORCE_H

#include "newmat/newmat.h"

#include "gutils/gtime.h"
#include "gforce/gfarg_cal.h"
#include "gio/glog.h"

#include "gutils/ExportOI.h"

namespace great
{
	/** The base class for all the force model. */
	class OI_LIBRARY_EXPORT t_gforce
	{
	public:
		/** @brief Constructor.*/
		t_gforce();

		/** @brief Constructor.
		* @param[in] log   the log file for write log.
		*/
		t_gforce(t_glog* log);

		/** @brief Destructor.*/
		virtual ~t_gforce();

		/**
	    * @brief an empty function for subclass to overrider for force calculation.
	    * @param[in]   argument     force calculate argument
	    * @return	bool
		*	 @retval  true  calculate successfuly
		*	 @retval  false calculate failed
	    */
		virtual bool calculate(t_gforce_calarg* argument);

		/**
	    * @brief get partial matrix.
	    * @return   partial matrix
	    */
		Matrix&			partial();
		
		/**
	    * @brief get the acceleration matrix.
	    * @return	acceleration matrix
	    */
		ColumnVector&	acc();

		/**
		* @brief get the forcemodel.
		* @return	force model 
		*/
		FORCE_MODEL	force_model();

	protected:
		Matrix		 _partial;	   ///< partial matrix.
		ColumnVector _acc;         ///< accelertion matrix.

		t_glog* 	 _glog;        ///< log file to save some important information.
		FORCE_MODEL  _force_model; ///< the Force model
		
	};
}
#endif // !GFORCE_H
