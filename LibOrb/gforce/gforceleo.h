/**
 *
 * @verbatim
     History
      -1.0 Yujie Qin   2019-03-05 creat the file.
	  -1.1 Yujie Qin   2019-07-03 Adding Doxygen Style Code Remarks
   @endverbatim
 *
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		    gforceleo.h
 * @brief		    The base class for all the force-leo models.
*
 * @author          Yujie Qin, Wuhan University
 * @version		    1.0.0
 * @date		    2019-03-05
 */
#ifndef GFORCELEO_H
#define GFORCELEO_H

#include "gforce/gforcemodel.h"
#include "gutils/ExportOI.h"
#include "gutils/gconst_oi.h"

namespace great
{
	/**
	*@brief	   Class for all the force model used in LEO satellites orbit integral.
	* 
	* Currently, we only consider gravity, N body attraction, relativity, 
	*	solar radiation, atmosphere and customer.
	*/
	class OI_LIBRARY_EXPORT t_gforce_leo : public t_gforce_model
	{
	public:
		/** 
		 * @brief constructor. 
		 */
		t_gforce_leo() {};

		/** 
		 * @brief destructor. 
		 */
		virtual ~t_gforce_leo();
		/**
	    *
		@verbatim
			Currently, we consider gravity, N body attraction, relativity, solar radiation, atmosphere and customer.
			- gravity         : 120th - order EGM08.
			- point mass : defalut model.
			- relativity : defalut model.
			- solar radiation : BOXWING model.
			- atmosphere : dtm94 model for density, one parameter.
			- customer : six parameters.
	    @endverbatim
	    *
	    * @param[in] arg : the init argument for init one force model.
	    
	    * -# t_init_gravity_egm08 : init EGM08 force model.
	    * -# t_init_pointmass : init pointmass force model.
	    * -# t_init_relativity : init relativity force model.
	    * -# t_init_solar_boxwing : init box wing force model.
	    * -# t_init_atmosphere : init atmosphere force model.
	    * -# t_init_customer : init customer force model.
	    *
	    * @note t_init_arg is a base class for the force init argument.
		*/
		bool creatForceModel(t_init_arg* arg, t_glog* log);

		//todo 变量的初始化
		t_gforce*	point_mass;	///< the point for N body force.
		t_gforce*	solar_radi;	///< the point for solar radiation.
		t_gforce*	gravity;    ///< the point for gravity.
		t_gforce*	relativity; ///< the point for relativity.
		t_gforce*   atmosphere; ///< the point for atmosphere.
		t_gforce*   earthradi; ///< the point for earthradiation.
		t_gforce*	customer;   ///< the point for customer.
	};
}

#endif // !GFORCELEO_H
