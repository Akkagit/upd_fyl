/**
 *
 * @verbatim
	 History
	  -1.0 Hongmin Zhang   2019-01-06 creat the file.
	  -1.1 Yujie Qin       2019-02-18 optim
	  -1.2 Jiande Huang     2019-02-20 edit for final version
	  -1.3 Jiande Huang     2019-03-14 Adding Doxygen Style Code Remarks
	  -1.4 Jiande Huang     2019-03-29 Adding the write_log function to Reporting procedural anomalies.
	                       2019-04-02 Adding exception throwing mechanism for all functions
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gfpointmass.cpp
 * @brief		The functions in gfpointmass.h are implemented in this file.
 *
 * @author      Hongmin Zhang , Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 *
 */

#include "gforce/gfpointmass.h"
#include "gutils/gtriple.h"
#include "gutils/ginfolog.h"

using namespace std;

namespace great
{

	/**
	* init the point mass force model by the init argument.
	* @param [in] arg init argument for point mass.
	*/
	t_gfpoint_mass::t_gfpoint_mass(t_init_arg * arg, t_glog* log) : t_gforce(log)
	{
		t_init_point_mass* point_mass = NULL;
		point_mass = dynamic_cast<t_init_point_mass*>(arg);
		if (!point_mass)
		{
			write_log(log, 0, "t_gfpoint_mass::t_gfpoint_mass", "ERROR : init failed");
			write_log(log, 0, "t_gfpoint_mass::t_gfpoint_mass", "ERROR : no init arg!");
		}
		else
		{
			if (!point_mass->gplanet)
			{
				write_log(log, 0, "t_gfpoint_mass::t_gfpoint_mass", "ERROR : init failed");
				write_log(log, 0, "ERROR : t_gfpoint_mass::t_gfpoint_mass", "no DE data !");
			}
			else
			{
				_force_model = arg->force_model;

				_gplanet = point_mass->gplanet;
				_partial.resize(3, 3);
				_acc.resize(3);

				// 插入卫星列表，此列表为N体引力考虑的行星的默认列表
				_planet_list.insert("EART");
				_planet_list.insert("SUN");
				_planet_list.insert("MOON");
				_planet_list.insert("MERC");
				_planet_list.insert("VENU");
				_planet_list.insert("MARS");
				_planet_list.insert("JUPI");
				_planet_list.insert("SATU");
				_planet_list.insert("URAN");
				_planet_list.insert("NEPT");
				_planet_list.insert("PLUT");

			}
		}	
	}
	 
	/**
	* calculate the acc for one planet, the planet should not be earth.
	* @param [in]  planet_name planet name.
	* @param [in]  epoch       epoch time in TT.
	* @param [in]  sat_pos     sat pos in J2000.
	* @param [out] planet_acc  the planet acceleration.
	*/
	bool t_gfpoint_mass::_cal_acc_mat(const string& planet_name, const double& epoch,const ColumnVector& sat_pos, ColumnVector& planet_acc)
	{
		t_gtriple point_pos(0.0,0.0,0.0);
		t_gtriple earth_pos(0.0,0.0,0.0);
		t_gtriple sat_xyz(sat_pos);
		t_gtriple sat_acc(0.0,0.0,0.0) ;
		
		// get the planet information from de405
		int idx	= str2planet(planet_name);
		double GM = _gplanet->get_GM(planet_name);

		//map<int, t_gplaneteph> all_planet = _gplanet->get_allplanets();
		if (GM == 0)
		{
			write_log(_glog, 0, "t_gfpoint_mass::calculate_acc", "ERROR : " + 
				planet_name + "doesn't have enough information in the function of pointmass");
			return false;
		}

		// get the planet and earth pos in j2000
		_gplanet->get_pos(epoch, planet_name, point_pos);
		_gplanet->get_pos(epoch, "EARTH", earth_pos);

		try
		{
			//calculate the distance between satellite and planet
			double dist1 = (point_pos - sat_xyz).norm();
			double dTmp1 = -GM / (dist1 * dist1 * dist1);
			sat_acc = (sat_xyz - point_pos) * dTmp1;

			// 计算一颗行星的加速度，具体公式可参考刘林或者Oliver的定轨书
			int idx_earth = str2planet("EARTH");
			if (idx != idx_earth)
			{
				double dist2 = (point_pos - earth_pos).norm();
				double dTmp2 = -GM / (dist2 * dist2 * dist2);
				sat_acc = sat_acc + point_pos * dTmp2;
			}

			// change triple to colum
			planet_acc = sat_acc.crd_cvect();
			return true;
		}
		catch (...)
		{
			write_log(_glog, 0, "t_gfpoint_mass::calculate_acc", "ERROR : throw exception");
			return false;
		}
	}

	/**
	* calculate the partial mat for one planet. The planet should not be earth.
	* @param [in]  planet_name planet name.
	* @param [in]  epoch       epoch time in TT.
	* @param [in]  sat_pos     sat pos in J2000.
	* @param [out] partial  the planet acceleration partial matrix.
	*/
	bool t_gfpoint_mass::_cal_partial_mat(const string& planet_name, const double& epoch,const ColumnVector& sat_pos, Matrix& partial)
	{
		t_gtriple point_pos(0.0,0.0,0.0);
		t_gtriple sat_xyz(sat_pos);

		//map<int, t_gplaneteph> all_planet = _gplanet->get_allplanets();
		_gplanet->get_pos(epoch, planet_name, point_pos);

		// judge for the planet information
		double GM = _gplanet->get_GM(planet_name);
		if (GM == 0)
		{
			write_log(_glog, 0, "t_gfpoint_mass::_calculate_partial","ERROR : " 
				+ planet_name + "doesn't have enough information in the function of pointmass");
			return false;
		}

		try
		{
			double dist1 = (point_pos - sat_xyz).norm();
			double dTmp1 = -GM / (dist1 * dist1 * dist1);
			double dTmp2 = -3 / (dist1 * dist1);

			//_partial_mat = 0;
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k <= j; k++)
				{
					double dDummy = dTmp2 * (sat_xyz[j] - point_pos[j]) * (sat_xyz[k] - point_pos[k]);
					if (j == k)
					{
						dDummy += 1;
					}
					dDummy = dDummy * dTmp1;
					partial(j + 1, k + 1) = partial(j + 1, k + 1) + dDummy;
					if (j != k)
					{
						partial(k + 1, j + 1) = partial(k + 1, j + 1) + dDummy;
					}
				}
			}
			return true;
		}
		catch (...)
		{
			write_log(_glog, 0, "t_gfpoint_mass::_calculate_partial", "ERROR : throw exception");
			return false;
		}
	}

	/** 
	* calculate the acc and partial for all planet.
	* @param [in] argument calculate argument for point mass.
	*/
	bool t_gfpoint_mass::calculate(t_gforce_calarg* argument)
	{
		string function_name = "t_gfpoint_mass::calculate";
		t_gfcalarg_point_mass* arg = NULL;
		arg = dynamic_cast<t_gfcalarg_point_mass*>(argument);
		if (!arg)
		{
			write_log(_glog, 0, "t_gfpoint_mass::calculate", "no pointmass calculate argument!");
			return false;
		}

		// 判断传入参数是否正确
		if (arg->arg_type != ARG_POINT_MASS)
		{
			write_log(_glog, 0, "t_gfpoint_mass::calculate", "wrong type calculate argument, not point mass!");
			return false;
		}

		// reset for zeros for safe
		_partial = 0.0;
		_acc 	 = 0.0;

		// the acc that one planet cause
		ColumnVector acc_one_planet;

		// planet name
		string point_param;

		set<string>::iterator iter;
		for (iter = _planet_list.begin(); iter != _planet_list.end(); ++iter)
		{
			point_param = *iter;
			// reset for zeros
			acc_one_planet = 0.0;
			// one planet
			if (!_cal_acc_mat(point_param, *(arg->time), *(arg->sat_pos), acc_one_planet))
			{
				write_log(_glog, 0, function_name, "ERROR : can not calculate acc of " + point_param);
				return false;
			};
			if (!_cal_partial_mat(point_param, *(arg->time), *(arg->sat_pos), _partial))
			{
				write_log(_glog, 0, function_name, "ERROR : can not calculate partial of " + point_param);
				return false;
			};
			// add for all planet
			_acc += acc_one_planet;
		}

		//write_log(_glog, 4, "t_gfpoint_mass::calculate", "normal end calculate!");
		return true;
	}
}