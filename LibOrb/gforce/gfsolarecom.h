/**
 *
 * @verbatim
	 History
	  -1.0 Bo Wong			 2018-12-23 creat the file.
	  -1.1 Jiande Huang      2019-02-20 edit for final version
	  -1.2 Jiande Huang      2019-03-14 Adding Doxygen Style Code Remarks
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gfsolarecom.h
 * @brief		ECOM force model, which used for getting solar radiation.
 *
 * @author      Bo Wong, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 *
 */
 
#ifndef GFSOLARECOM_H
#define GFSOLARECOM_H

#include "gio/glog.h"
#include "gforce/gfsolar.h"
#include "gforce/gfarg_init.h"
#include "gutils/gconst_oi.h"

#include "newmat/newmat.h"

namespace great
{
	/**
	 * @brief The class for ECOM force model.
	 */
	class OI_LIBRARY_EXPORT t_gfsolar_ecom : public t_gfsolar
	{
	public:
		/// constructor
		//t_gfsolar_ecom() {};

		/** 
		 * @brief constructor
		 * @param [in] log	The log file point
		 */
		t_gfsolar_ecom(t_glog* log);

		/**
		 * @brief constructor
		 * @param [in] arg The init argument for the force model
		 * @param [in] log The log file point
		 */
		t_gfsolar_ecom(t_init_arg* arg, t_glog* log);

		/**
		 * @brief destructor
		 */
		virtual ~t_gfsolar_ecom() {};

		/**
		 * @brief show the type of model
		 */
		void show() override;

		/**
		 * @brief check some important pars in the model
		 * @return  whether the model is valid
		 *  @retval false The model is not valid
		 *  @retval true  The model is valid
		 */	
		bool check();

		/**
		* @brief calculate the acc and partial
		* Calculate the solar radiation with ECOM force model, The shadow factor is
		* considered here.
		* @param [in] argument the calculate argument for ECOM models.
		* @return the calculate result
		*  @retval false calculate failed
		*  @retval true  calculate succeed
		*/
		bool calculate(t_gforce_calarg* argument) override;

	protected:
		
		std::map<ICS_PAR,double>	_ecom_pars;	      ///< 9 solar pars of ECOM
		const double _aunit = 1.4959787061e8; ///< Intermediate variables
	};

}
#endif // !GFSOLARECOM_H
