/**
 *
 * @verbatim
	 History
	  -1.0 Yujie Qin        2019-01-06  creat the file.
	  -1.1 Yujie Qin        2019-02-18  optim
	  -1.2 Jiande Huang      2019-02-20  edit for final version
	  -1.3 Jiande Huang      2019-03-14  Adding Doxygen Style Code Remarks
	  -1.4 Yujie Qin        2019-07-03  Adding Doxygen Style Code Remarks
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gfrelativity.h
 * @brief		Relativity model, which used for getting force caused by relativistic effect.
 *
 * @author      Yujie Qin, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 *
 */
 
#ifndef GFRELATIVITY_H
#define GFRELATIVITY_H

#include <string.h>
#include <stdio.h>
#include <iostream>

#include "gutils/gconst.h"
#include "gutils/gtriple.h"
#include "gforce/gforce.h"
#include "gforce/gfarg_init.h"

#include "newmat/newmat.h"

using namespace std;

namespace great
{
	/** @brief class for relativity force */
	class t_gfrelativity : public t_gforce
	{
	public:
		/** @brief constructor. */
		t_gfrelativity() {};

		/** @brief constructor. */
		t_gfrelativity(t_init_arg* arg, t_glog* log);

		/** @brief destructor. */
		virtual ~t_gfrelativity() {};
		
		/**
	    * @brief calculate the acceleration derived by the relativity.
	    * @param[in]  argument  calculated parameters
	    * @return  whether the calculate valid
	    *	 @retval  true  calculate successfuly
	    *	 @retval  false calculate failed
	    */
		bool calculate(t_gforce_calarg* argument);                                               

	};
} // namespace

#endif // !GRELATIVITY_H
