/**
 *
 * @verbatim
     History
      -1.0 Wei Zhang    2019-01-20 Creat the file.
      -1.1 Yujie Qin    2019-03-04 edit for final version.
      -1.2 Yujie Qin    2019-03-15 Code comments are added in Doxygen format.
	  -1.3 Jiande Huang  2019-04-01 Adding exception throwing mechanism for all functions

   @endverbatim
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		      gfsolarboxwing.cpp
 * @brief	          general solar radiation model, calculate as sum of accelerations all pannels.
 *
 * @author            Wei Zhang, Wuhan University
 * @version	          1.0.0
 * @date		      2019-01-20
 *
 */
#include "gforce/gfsolarboxwing.h"
#include "gutils/ginfolog.h"

using namespace std;

namespace great
{

	/** constructor. */
	t_gfsolarbow::t_gfsolarbow(t_init_arg * arg,t_glog* log) : t_gfsolar(log)
	{
		string function_name = "t_gfsolarboxwing::t_gfsolarboxwing";

		t_init_solar_boxwing *box_wing = NULL;
		box_wing = dynamic_cast<t_init_solar_boxwing*>(arg);

		if (!box_wing)
		{
			write_log(log, 0, function_name, "ERROR : init failed");
			write_log(log, 0, function_name, "ERROR : no init argument");

		}
		else
		{
			if (!box_wing->pannel)
			{
			  write_log(log, 0, function_name, "ERROR : init failed");
			  write_log(log, 0, function_name, "ERROR : no pannel data");
			}
			else
			{
				_pannel    = box_wing->pannel;
				_ispartial = box_wing->ispartial;
				_mass      = box_wing->mass;
			}
		}
	}
	
	/** deconstructor. */
	t_gfsolarbow::~t_gfsolarbow()
	{
		
	}

	/** show the detail of the BOXWING model. */
	void t_gfsolarbow::show()
	{
		cout << "solar model : boxwing" << endl;
	}

	/** check the important parameters and information in the model */
	bool t_gfsolarbow::check()
	{
		cout << "check some important pars in model : boxwing" << endl;
		return true;
		// TODO
	}

	/** calculate the acceleration and the corresponding partial matrix.
	*@param [in] argument Parameters to be used in calculation.
	*/
	bool t_gfsolarbow::calculate(t_gforce_calarg* argument)
	{
		string function_name = "t_gfsolarboxwing::calculate";
		t_gfcalarg_solar_boxwing* arg = NULL;
		arg = dynamic_cast<t_gfcalarg_solar_boxwing*>(argument);
		if (!arg)
		{
			write_log(_glog, 0, function_name, "ERROR : no calculate argument");
			return false;
		}
		try
		{

			_acc = 0.0;
			_boxw_par = *(arg->solar_par);
			static bool isfirst = true;
			int index = 0;
			static int parindex = 0;
			map<ICS_PAR, double>::iterator it;
			double parvalue = 1.0;
			bool lpar = true;
			if (_boxw_par.size() == 0)
			{
				lpar = false;
			}
			// For one satellite, the index of solar par in the parameter list will not change
			// So it's not necessary to judge the index every time, just for the first time
			if (lpar)
			{
				if (isfirst)
				{
					isfirst = false;
					for (it = _boxw_par.begin(); it != _boxw_par.end(); it++)
					{
						index++;
						if (it->first == _parname)
						{
							parindex = index;
							parvalue += it->second;
							break;
						}
					}
				}
			}

			double lambda = *(arg->shadow_factor);

			// unit vector of sat
			ColumnVector unit_sat = *(arg->sat_pos) / arg->sat_pos->norm_Frobenius();

			// vector of satellite to sun
			ColumnVector sat2sun = *(arg->sun_pos) - *(arg->sat_pos);

			// unit vector of satellite to sun
			ColumnVector unit_sat2sun = sat2sun / sat2sun.norm_frobenius();


			double dist_factor = pow(_aunit / sat2sun.NormFrobenius(), 2);
			int n = _pannel->get_pannels()->size();

			// acceleration from all panels
			ColumnVector normj;
			normj.resize(3);
			for (int i = 1; i <= n; i++)
			{
				normj = (*_pannel->get_pannels())[i].crs_norm_dir;
				double cosa = dotproduct(normj, normj);
				if (cosa < 1e-10)
				{
					normj = unit_sat2sun;
				}
				cosa = dotproduct(unit_sat2sun, normj);
				if (cosa <= 0)
				{
					continue;
				}
				_acc -= (*_pannel->get_pannels())[i].area*cosa*
					(2.0*((*_pannel->get_pannels())[i].dr_coef / 3.0 + (*_pannel->get_pannels())[i].sr_coef*cosa)*normj +
					(1.0 - (*_pannel->get_pannels())[i].sr_coef)*unit_sat2sun);
			}
			if (_mass != 0.0)
			{
				_acc *= parvalue * lambda*dist_factor*_ps*1.0e-3 / _mass;
			}

			// jdhuang
			//if (_ispartial)
			//{
			//	if (lpar)
			//	{
			//		if (parindex != 0)
			//		{
			//			_partial = _acc / parvalue;
			//		}
			//	}
			//}
			// partial matrix
			if (_ispartial && lpar && parindex != 0)
			{
				_partial = _acc / parvalue;
			}
			return true;
		}
		catch (...)
		{
			write_log(_glog, 0, "ERROR : " + function_name + " throw exception while calculating Density");
			return false;
		}
	}

	
} // namespace