/**
 *
 * @verbatim
	 History	
	  -1.0 Hongmin Zhang    2019-01-06 creat the file.
	  -1.1 Yujie Qin        2019-02-18 optim
	  -1.2 Jiande Huang     2019-02-20 edit for final version
	  -1.3 Jiande Huang     2019-03-14 Adding Doxygen Style Code Remarks
	  -1.4 Yujie Qin        2019-03-14 Adding Doxygen Style Code Remarks
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gfpointmass.h
 * @brief		Point mass force model, which used for getting the N body attrations.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 *
 */
 
 
#ifndef GFPOINTMASS_H
#define GFPOINTMASS_H

#include "gdata/gnavde.h"
#include "gutils/gtriple.h"
#include "gall/gallics_nav.h"
#include "gforce/gforce.h"
#include "gforce/gfarg_init.h"

//todo 三方库移出源文件
#include "newmat/newmat.h"
#include "newmat/newmatio.h"

using namespace std;
namespace great
{
	/** @brief The class for point mass force model.*/
	class OI_LIBRARY_EXPORT t_gfpoint_mass : public t_gforce
	{
	public:
		/** @brief constructor. */
		t_gfpoint_mass() {};

		/**
	    * @brief constructor,init the point mass force model by the init argument.
	    * @param[in] arg init argument for point mass.
	    */
		t_gfpoint_mass(t_init_arg* arg , t_glog* log);
		
		/** @brief destructor. */
		virtual ~t_gfpoint_mass() {};
		
		/**
	    * @brief calculate the acc and partial for all planet.
	    * @param[in] argument calculate argument for point mass
		* @return bool
		*	 @retval  true  calculate successfuly
		*	 @retval  false calculate failed
	    */
		bool calculate(t_gforce_calarg* argument);

	protected:
		/**
	    * @brief calculate the acc for one planet, the planet should not be earth.
	    * @param[in]  planet_name planet name.
	    * @param[in]  epoch       epoch time in TT.
	    * @param[in]  sat_pos     sat pos in J2000.
	    * @param[out] planet_acc  the planet acceleration.
		* @return bool
		*	 @retval  true  calculate successfuly
		*	 @retval  false calculate failed
	    */
		bool _cal_acc_mat(const string& planet_name, const double& epoch,const ColumnVector& sat_pos, ColumnVector& planet_acc);		
		
		/**
	    * @brief calculate the partial mat for one planet. The planet should not be earth.
	    * @param[in]  planet_name planet name.
	    * @param[in]  epoch       epoch time in TT.
	    * @param[in]  sat_pos     sat pos in J2000.
	    * @param[out] partial  the planet acceleration partial matrix.
	    */
		bool _cal_partial_mat(const string& planet_name, const double& epoch,const ColumnVector& sat_pos, Matrix& partial);									

		t_gnavde*      _gplanet;	 ///< all planets
		set<string>    _planet_list; ///< the planet type for this point mass calculation
	};
} //namespace

#endif // !GPOINTMASS_H
