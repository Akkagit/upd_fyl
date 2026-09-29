/**
 *
 * @verbatim
     History
      -1.0 Yongqiang Yuan    2019-11-02 Creat the file.

   @endverbatim
 * Copyright (c) 2019, Wuhan University.All rights reserved.
 *
 * @file		      gfsolarboxwing.cpp
 * @brief	          general solar radiation model, calculate as sum of accelerations all pannels.
 *
 * @author            Yongqiang Yuan, Wuhan University
 * @version	          1.0.0
 * @date		      2019-11-02
 *
 */
#include <iostream>
#include "gutils/gconst.h"
#include "gutils/ginfolog.h"
#include "gforce/gfarg_cal.h"
#include "gforce/gfcommon.h"
#include "gforce/gfsolarprioriboxwing.h"

using namespace std;

namespace great
{

	/** constructor. */
	t_gfsolar_pbw::t_gfsolar_pbw(t_init_arg * arg,t_glog* log) : t_gfsolar(log)
	{
		string function_name = "t_gfsolarprioribowwing::t_gfsolarprioribowwing";

		t_init_solar_priori_boxwing *priori_box_wing = NULL;
		priori_box_wing = dynamic_cast<t_init_solar_priori_boxwing*>(arg);

		if (!priori_box_wing)
		{
			write_log(log, 0, function_name, "ERROR : init failed");
			write_log(log, 0, function_name, "ERROR : no init argument");

		}
		else
		{
			_force_model = priori_box_wing->force_model;
			if (!priori_box_wing->pannel)
			{
			  write_log(log, 0, function_name, "ERROR : init failed");
			  write_log(log, 0, function_name, "ERROR : no pannel data");
			}
			else
			{
				_pannel    = priori_box_wing->pannel;
				//_ispartial = priori_box_wing->ispartial;
				
				_pbw_pars[ICS_PAR::ECOM_D0] = (priori_box_wing->solar_pars.count(ICS_PAR::ECOM_D0) > 0) ? priori_box_wing->solar_pars[ICS_PAR::ECOM_D0] : 0;
				_pbw_pars[ICS_PAR::ECOM_Y0] = (priori_box_wing->solar_pars.count(ICS_PAR::ECOM_Y0) > 0) ? priori_box_wing->solar_pars[ICS_PAR::ECOM_Y0] : 0;
				_pbw_pars[ICS_PAR::ECOM_X0] = (priori_box_wing->solar_pars.count(ICS_PAR::ECOM_X0) > 0) ? priori_box_wing->solar_pars[ICS_PAR::ECOM_X0] : 0;
				_pbw_pars[ICS_PAR::ECOM_Dc] = (priori_box_wing->solar_pars.count(ICS_PAR::ECOM_Dc) > 0) ? priori_box_wing->solar_pars[ICS_PAR::ECOM_Dc] : 0;
				_pbw_pars[ICS_PAR::ECOM_Ds] = (priori_box_wing->solar_pars.count(ICS_PAR::ECOM_Ds) > 0) ? priori_box_wing->solar_pars[ICS_PAR::ECOM_Ds] : 0;
				_pbw_pars[ICS_PAR::ECOM_Yc] = (priori_box_wing->solar_pars.count(ICS_PAR::ECOM_Yc) > 0) ? priori_box_wing->solar_pars[ICS_PAR::ECOM_Yc] : 0;
				_pbw_pars[ICS_PAR::ECOM_Ys] = (priori_box_wing->solar_pars.count(ICS_PAR::ECOM_Ys) > 0) ? priori_box_wing->solar_pars[ICS_PAR::ECOM_Ys] : 0;
				_pbw_pars[ICS_PAR::ECOM_Xc] = (priori_box_wing->solar_pars.count(ICS_PAR::ECOM_Xc) > 0) ? priori_box_wing->solar_pars[ICS_PAR::ECOM_Xc] : 0;
				_pbw_pars[ICS_PAR::ECOM_Xs] = (priori_box_wing->solar_pars.count(ICS_PAR::ECOM_Xs) > 0) ? priori_box_wing->solar_pars[ICS_PAR::ECOM_Xs] : 0;

				_partial.resize(9, 3);
			}
		}
	}
	
	/** deconstructor. */
	t_gfsolar_pbw::~t_gfsolar_pbw()
	{
		
	}

	/** show the detail of the PRIORI BOXWING model. */
	void t_gfsolar_pbw::show()
	{
		cout << "solar model : priori box wing and ECOM" << endl;
	}

	/** check the important parameters and information in the model */
	bool t_gfsolar_pbw::check()
	{
		cout << "check some important pars in model : priori box wing and ECOM" << endl;
		return true;
		// TODO
	}

	/** calculate the acceleration and the corresponding partial matrix.
	*@param [in] argument Parameters to be used in calculation.
	*/
	bool t_gfsolar_pbw::calculate(t_gforce_calarg* argument)
	{
		string function_name = "t_gfsolarprioribowwing::calculate";
		t_gfcalarg_priori_solar_boxwing* arg = NULL;
		arg = dynamic_cast<t_gfcalarg_priori_solar_boxwing*>(argument);
		if (!arg)
		{
			write_log(_glog, 0, function_name, "ERROR : no calculate argument");
			return false;
		}
		try
		{

			// unit vector z=(0,0,1.0)
			ColumnVector unit_Z(0.0, 3);
			unit_Z(3) = 1.0;

			// zunit is defined from satellite to the earth
			ColumnVector unit_sat = *(arg->sat_pos) / arg->sat_pos->norm_Frobenius();

			// vector sun to satellite
			// unit vector of sun2sat 
			ColumnVector sun2sat = *(arg->sat_pos) - *(arg->sun_pos);

			// unit vector of sun to satellite
			ColumnVector unit_sun2sat = sun2sat / sun2sat.norm_frobenius();

			// the sqrt of sum of 3 elements in es2s vector
			double length = sun2sat.norm_Frobenius();
			double dist_factor = pow(_aunit / length, 2);

			// Unit vector form satellite to solar direction 
			ColumnVector unit_D = -unit_sun2sat;

			// unit vector of orbit plane
			ColumnVector unit_Nop = crossproduct(*(arg->sat_pos), *(arg->sat_vel));
			unit_Nop = unit_Nop / unit_Nop.norm_Frobenius();

			// unit vector of _sun_pos
			ColumnVector unit_sun = *(arg->sun_pos) / arg->sun_pos->norm_Frobenius();

			//beta angel between orbit plane and sun (not used)
			//dBeta = asin(dotproduct(_eNop, esun) / sqrt(dotproduct(esun,esun)));

			//u angel from node to satellite along orbit plane
			ColumnVector unit_Node = crossproduct(unit_Z, unit_Nop);
			unit_Node = unit_Node / unit_Node.norm_Frobenius();
			double fx = dotproduct(unit_Node, unit_sat);
			ColumnVector unit_tmp = crossproduct(unit_Node, unit_sat);
			double fz = dotproduct(unit_tmp, unit_Nop);
			double u = atan2(fz, fx);

			//u angel from midnight to satellite along orbit plane  (u_midnight,'orbit angle')
			ColumnVector unit_Terminator = crossproduct(unit_sun, unit_Nop);
			unit_Terminator = unit_Terminator / unit_Terminator.norm_Frobenius();
			ColumnVector unit_midnight = crossproduct(unit_Terminator, unit_Nop);
			unit_midnight = unit_midnight / unit_midnight.norm_Frobenius();
			fx = dotproduct(unit_midnight, unit_sat);
			unit_tmp = crossproduct(unit_midnight, unit_sat);
			fz = dotproduct(unit_tmp, unit_Nop);
			double u_midnight = atan2(fz, fx);

			// update the shadow factor
			double lambda = *(arg->shadow_factor);

			//eY : Unit vector from satellite to geocentric direction
			ColumnVector unit_Y(0.0, 3);
			if (arg->sat_type.substr(0, 8) == "BLOCKIIR" && lambda < 1.0)
			{
				unit_Y = crossproduct(-unit_sat, *(arg->sat_vel));
			}
			else
			{
				unit_Y = crossproduct(*(arg->sat_pos), unit_sun2sat);
			}
			lambda = (lambda < 0.999 ? 0.0 : lambda);
			unit_Y = unit_Y / unit_Y.norm_Frobenius();

			// Unit vector of satellite solar panel axial 
			ColumnVector unit_X = crossproduct(unit_Y, unit_sun2sat);
			unit_X = unit_X / unit_X.norm_Frobenius();

			ColumnVector DYXu(0.0, 3);
			ColumnVector DYXf(0.0, 3);

			double cos_u = cos(u_midnight);
			double sin_u = sin(u_midnight);

			// calcluate acc in sat fix
			DYXu(1) = _pbw_pars[ICS_PAR::ECOM_D0] + (_pbw_pars[ICS_PAR::ECOM_Dc] * cos_u + _pbw_pars[ICS_PAR::ECOM_Ds] * sin_u) * lambda;
			DYXu(2) = _pbw_pars[ICS_PAR::ECOM_Y0] + (_pbw_pars[ICS_PAR::ECOM_Yc] * cos_u + _pbw_pars[ICS_PAR::ECOM_Ys] * sin_u) * lambda;
			DYXu(3) = _pbw_pars[ICS_PAR::ECOM_X0] + (_pbw_pars[ICS_PAR::ECOM_Xc] * cos_u + _pbw_pars[ICS_PAR::ECOM_Xs] * sin_u) * lambda;

			// acc in crs
			DYXf = (unit_D * DYXu(1) + unit_Y * DYXu(2) - unit_X * DYXu(3)) * (*(arg->shadow_factor));

			// calculate acc, km/s^2
			_acc = DYXf * dist_factor * 1e-9;

			// Calculate the partial derivative of DYX, parameter sequence: d y b 
			// reset zeros for safe
			_partial = 0.0;
			for (int i = 1; i <= 3; i++)
			{
				if (_force_model == F_SOLAR_PBW_5)
				{
					_partial(1, i) = unit_D(i) * (*(arg->shadow_factor));
					_partial(2, i) = unit_Y(i) * (*(arg->shadow_factor));
					_partial(3, i) = -unit_X(i) * (*(arg->shadow_factor));
					_partial(4, i) = -unit_X(i) *  cos_u * lambda;
					_partial(5, i) = -unit_X(i) *  sin_u * lambda;
				}
				if (_force_model == F_SOLAR_PBW_9)
				{
					_partial(1, i) = unit_D(i) * (*(arg->shadow_factor));	
					_partial(2, i) = unit_D(i) *  cos_u * lambda;
					_partial(3, i) = unit_D(i) *  sin_u * lambda;

					_partial(4, i) = unit_Y(i) * (*(arg->shadow_factor));
					_partial(5, i) = unit_Y(i) *  cos_u * lambda;
					_partial(6, i) = unit_Y(i) *  sin_u * lambda;

					_partial(7, i) = -unit_X(i) * (*(arg->shadow_factor));
					_partial(8, i) = -unit_X(i) *  cos_u * lambda;
					_partial(9, i) = -unit_X(i) *  sin_u * lambda;
				}
			}
			_partial = _partial * dist_factor * 1e-9;



			ColumnVector unit_sat2sun = (-1) * unit_sun2sat;
			ColumnVector f_box_wing(0.0, 3);

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

				// Note: a provisional practice ot use <idr of LEO> for <borw ids for GNS>
				//       a provisional practice ot use <sr of LEO> for <abso ids for GNS>
				//       a provisional practice ot use <dr of LEO> for <spec ids for GNS>
				//       a provisional practice ot use <drag of LEO> for <diff ids for GNS>
				double borw = (*_pannel->get_pannels())[i].idr_coef;
				double abso = (*_pannel->get_pannels())[i].sr_coef;
				double spec = (*_pannel->get_pannels())[i].dr_coef;
				double diff = (*_pannel->get_pannels())[i].drag_coef;
				double area = (*_pannel->get_pannels())[i].area;
				//double ps = 4.5604e-6;
				// for solar panels
				if (borw > -2.0 && borw < 0.0) {
					f_box_wing -= area*cosa * ((1.0-spec)*unit_D+2.0*(1.0/3.0*diff+spec*cosa)*(*_pannel->get_pannels())[i].crs_norm_dir);
				}
				// for satellite bus
				else {
					f_box_wing -= area*cosa * ((abso + diff)*(unit_D + 2.0 / 3.0*(*_pannel->get_pannels())[i].crs_norm_dir) + 2.0*spec*cosa*(*_pannel->get_pannels())[i].crs_norm_dir);
				}
			}
			if (arg->sat_mass != 0.0)
			{
				f_box_wing = f_box_wing * lambda * dist_factor * _ps * 1.0e-3 / (arg->sat_mass);
				_acc = _acc + f_box_wing;
			}
			else {
				f_box_wing = f_box_wing * 0.0;
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