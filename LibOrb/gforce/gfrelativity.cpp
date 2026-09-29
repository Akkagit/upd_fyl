/**
 * Purpose: get the acceleration caused by relativity.
 * @verbatim
	 History
	  -1.0 Yujie Qin       2019-01-06  creat the file.
	  -1.1 Yujie Qin       2019-02-18  optim
	  -1.2 Jiande Huang     2019-02-20  edit for final version
	  -1.3 Jiande Huang     2019-03-14  Adding Doxygen Style Code Remarks
	  -1.4 Jiande Huang     2019-03-29  Add the function of write_log to report.
	                       2019-04-01  Adding exception throwing mechanism for all functions
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gfrelativity.cpp
 *
 * @author      Yujie Qin, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 *
 */


#include "gutils/ginfolog.h"
#include "gforce/gfrelativity.h"

using namespace std;

namespace great
{
	/** @brief constructor. */
	t_gfrelativity::t_gfrelativity(t_init_arg* arg, t_glog* log) : t_gforce(log)
	{
		if (!arg)
		{
			write_log(_glog, 0, "t_gfrelativity::t_gfrelativity", "init failed");
			write_log(_glog, 0, "t_gfrelativity::t_gfrelativity", "no init argument!");
		}
		else
		{
			_force_model = arg->force_model;
		}
	}

	/**
	* @brief calculate the acceleration derived by the relativity.
	* @param[in]  argument  calculated parameters
	* @return  whether the calculate valid
	*		@retval  true  calculate successfuly
	*		@retval  false calculate failed
	*/
	bool t_gfrelativity::calculate(t_gforce_calarg* argument)
	{
		string function_name = "t_gfrelativity::calculate";
		
		// change the type of argument
		t_gfcalarg_relativity* arg = NULL;
		arg = dynamic_cast<t_gfcalarg_relativity*>(argument);

		// if not argument, calculate failed
		if (!arg)
		{
			write_log(_glog, 0, "ERROR : " + function_name + " ==> no relativity argument!");
			return false;
		}
		// if wrong argument
		if (arg->arg_type != ARG_RELATIVITY)
		{
			write_log(_glog, 0, "ERROR : " + function_name + " ==> wrong relativity argument!");
			return false;
		}

		try
		{
			// reset for safe
			_acc = 0.0;
			//double vright = 299792458E-3;
			// 卫星1D速度和位置
			double r = arg->sat_pos->norm_Frobenius();
			double v = arg->sat_vel->norm_Frobenius();
			double rv = dotproduct(arg->sat_pos->t(), *(arg->sat_vel));

			// 计算加速度，具体公式可以参考博士论文或者相关书籍
			double c1 = (arg->gm) / r;
			double c2 = c1 / r / r / CLIGHT_2 / 1e-6;
			//double c2 = c1 / r / r / vright / vright;
			double alpha = (4.0 * c1 - v * v) * c2;
			double beta  = 4.0 * rv * c2;

			_acc = *(arg->sat_pos) * alpha + *(arg->sat_vel) * beta;

			return true;
#ifdef DEBUG
			cout << "The relativity model : " << endl
				<< "Sat pos is :" << setw(30) << setprecision(15) << scientific << *(arg->satpos) << endl
				<< "Sat vel is :" << setw(30) << setprecision(15) << scientific << *(arg->satvel) << endl
				<< "Acc     is :" << setw(30) << setprecision(15) << scientific << _acc << endl;
#endif // DEBUG

		}
		catch (...)
		{
			write_log(_glog, 0, "ERROR : " + function_name + " ==> throw exception");
			return false;
		}
	}
}
