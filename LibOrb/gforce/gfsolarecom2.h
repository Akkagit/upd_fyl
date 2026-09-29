/**
 *
 * @verbatim
	 History
	  -1.0 Yongqinag Yuan			 2019-11-02 creat the file.
   @endverbatim
 * Copyright (c) 2019, Wuhan University. All rights reserved.
 *
 * @file		gfsolarecom2.h
 * @brief		ECOM2 force model, which used for getting solar radiation.
 *
 * @author      Yongqiang Yuan, Wuhan University
 * @version		1.0.0
 * @date		2019-11-02
 *
 */
 
#ifndef GFSOLARECOM2_H
#define GFSOLARECOM2_H

#include "gio/glog.h"
#include "gforce/gfsolar.h"
#include "gforce/gfarg_init.h"
#include "gutils/gconst_oi.h"

#include "newmat/newmat.h"

namespace great
{
	/**
	 * @brief The class for ECOM2 force model.
	 */
	class t_gfsolar_ecom2 : public t_gfsolar
	{
	public:
		/// constructor
		//t_gfsolar_ecom2() {};

		/** 
		 * @brief constructor
		 * @param [in] log	The log file point
		 */
		t_gfsolar_ecom2(t_glog* log);

		/**
		 * @brief constructor
		 * @param [in] arg The init argument for the force model
		 * @param [in] log The log file point
		 */
		t_gfsolar_ecom2(t_init_arg* arg, t_glog* log);

		/**
		 * @brief destructor
		 */
		virtual ~t_gfsolar_ecom2() {};

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
		
		std::map<ICS_PAR,double> _ecom2_pars;	      ///< 9 solar pars of ECOM2
		const double _aunit = 1.4959787061e8; ///< Intermediate variables
	};

}
#endif // !GFSOLARECOM2_H
