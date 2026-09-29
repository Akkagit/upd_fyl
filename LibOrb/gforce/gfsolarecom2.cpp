/**
 *
 * @verbatim
	 History
	  -1.0 Yongqiang Yuan			 2019-11-02 creat the file.
   @endverbatim
 * Copyright (c) 2019, Wuhan University. All rights reserved.
 *
 * @file		gfsolarecom2.cpp
 * @brief		ECOM2 force model, which used for getting solar radiation.
 *
 * @author      Yongqiang Yuan, Wuhan University
 * @version		1.0.0
 * @date		2019-11-02
 *
 */

#include <iostream>
#include "gutils/gconst.h"
#include "gutils/ginfolog.h"
#include "gforce/gfarg_cal.h"
#include "gforce/gfcommon.h"
#include "gforce/gfsolarecom2.h"


namespace great
{
	 t_gfsolar_ecom2::t_gfsolar_ecom2(t_glog * log) : t_gfsolar(log)
	 {
		if (!log)
		{
			write_log(log, 0, "t_gfsolar_ecom2::t_gfsolar_ecom2",
				"ERROR : the log point is empty");
		}
	}

	 t_gfsolar_ecom2::t_gfsolar_ecom2(t_init_arg * arg, t_glog * log) : t_gfsolar(log)
	 {
		 string function_name = "t_gfsolar_ecom2::t_gfsolar_ecom2";
		 try
		 {
			 t_init_solar_ecom2* ecom2 = NULL;
			 ecom2 = dynamic_cast<t_init_solar_ecom2*>(arg);
			 if (!ecom2)
			 {
				 write_log(log, 0, "ERROR : " + function_name + " ==> init failed");
				 write_log(log, 0, "ERROR : " + function_name + " ==> no init argument");
			 }
			 else
			 {
				 if (ecom2->solar_pars.size() == 9 && arg->force_model == F_SOLAR_ECOM2)
				 {
					 _force_model = arg->force_model;
					 // get params' value from ics file according to epoch
					 _ecom2_pars[ICS_PAR::ECOM2_D0]  = (ecom2->solar_pars.count(ICS_PAR::ECOM2_D0) > 0)  ? ecom2->solar_pars[ICS_PAR::ECOM2_D0]  : 0;
					 _ecom2_pars[ICS_PAR::ECOM2_D2c] = (ecom2->solar_pars.count(ICS_PAR::ECOM2_D2c) > 0) ? ecom2->solar_pars[ICS_PAR::ECOM2_D2c] : 0;
					 _ecom2_pars[ICS_PAR::ECOM2_D2s] = (ecom2->solar_pars.count(ICS_PAR::ECOM2_D2s) > 0) ? ecom2->solar_pars[ICS_PAR::ECOM2_D2s] : 0;
					 _ecom2_pars[ICS_PAR::ECOM2_D4c] = (ecom2->solar_pars.count(ICS_PAR::ECOM2_D4c) > 0) ? ecom2->solar_pars[ICS_PAR::ECOM2_D4c] : 0;
					 _ecom2_pars[ICS_PAR::ECOM2_D4s] = (ecom2->solar_pars.count(ICS_PAR::ECOM2_D4s) > 0) ? ecom2->solar_pars[ICS_PAR::ECOM2_D4s] : 0;
					 _ecom2_pars[ICS_PAR::ECOM2_Y0]  = (ecom2->solar_pars.count(ICS_PAR::ECOM2_Y0) > 0)  ? ecom2->solar_pars[ICS_PAR::ECOM2_Y0]  : 0;
					 _ecom2_pars[ICS_PAR::ECOM2_X0]  = (ecom2->solar_pars.count(ICS_PAR::ECOM2_X0) > 0)  ? ecom2->solar_pars[ICS_PAR::ECOM2_X0]  : 0;
					 _ecom2_pars[ICS_PAR::ECOM2_X1c] = (ecom2->solar_pars.count(ICS_PAR::ECOM2_X1c) > 0) ? ecom2->solar_pars[ICS_PAR::ECOM2_X1c] : 0;
					 _ecom2_pars[ICS_PAR::ECOM2_X1s] = (ecom2->solar_pars.count(ICS_PAR::ECOM2_X1s) > 0) ? ecom2->solar_pars[ICS_PAR::ECOM2_X1s] : 0;

					 _partial.resize(9, 3);

				 }
				 else
				 {
					 write_log(log, 0, "ERROR : " + function_name + " ==> init failed");
					 write_log(log, 0, "ERROR : " + function_name + " ==> error data for ics init value");
				 }
			 }

#ifdef DEBUG
			 cout << "The solar force model is " << force2str(_force_model) << endl
				 << "The init value of solar pars is : " << endl
				 << " ECOM2_D0 " << setw(30) << setprecision(15) << scientific << _ecom2_pars[ECOM2_D0] << endl
				 << " ECOM2_D2c " << setw(30) << setprecision(15) << scientific << _ecom2_pars[ECOM2_D2c] << endl
				 << " ECOM2_D2s " << setw(30) << setprecision(15) << scientific << _ecom2_pars[ECOM2_D2s] << endl
				 << " ECOM2_D4c " << setw(30) << setprecision(15) << scientific << _ecom2_pars[ECOM2_D4c] << endl
				 << " ECOM2_D4s " << setw(30) << setprecision(15) << scientific << _ecom2_pars[ECOM2_D4s] << endl
				 << " ECOM2_Y0 " << setw(30) << setprecision(15) << scientific << _ecom2_pars[ECOM2_Y0] << endl
				 << " ECOM2_X0 " << setw(30) << setprecision(15) << scientific << _ecom2_pars[ECOM2_X0] << endl
				 << " ECOM2_X1c " << setw(30) << setprecision(15) << scientific << _ecom2_pars[ECOM2_X1c] << endl
				 << " ECOM2_X1s " << setw(30) << setprecision(15) << scientific << _ecom2_pars[ECOM2_X1s] << endl;
#endif // DEBUG

		 }
		 catch (...)
		 {
			 write_log(log, 0, "ERROR : " + function_name + " ==> throw exception while constructing ");
		 }
	 }

	/** @todo Functional functions still need to be improved */
	void t_gfsolar_ecom2::show()
	{
		cout << "solar model : ECOM2" << endl;
	}

	/** Functional functions still need to be improved */
	bool t_gfsolar_ecom2::check()
	{
		write_log(_glog, 4, "check some important pars in model : ECOM2");
		if (_ecom2_pars.size() == 0)
		{
			write_log(_glog, 0, "ERROR : no ics data in ECOM model!");
			return false;
		}
		return true;
		// TODO
	}


	bool t_gfsolar_ecom2::calculate(t_gforce_calarg* argument)
	{
		string function_name = "t_gfsolar_ecom2::calculate";
		// 判断计算的充分条件，并获取参数
		t_gfcalarg_solar_ecom2* arg = NULL;
		arg = dynamic_cast<t_gfcalarg_solar_ecom2*>(argument);
		if (!arg)
		{
			write_log(_glog, 0, "ERROR : " + function_name + " ==> no calculate argument");
			return false;
		}

		try
		{
			// unit vector z=(0,0,1.0)
			ColumnVector unit_Z(0.0, 3);
			unit_Z(3) = 1.0;

			// unit_sat is defined from satellite to the earth
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

			//u angel from the ascending node to satellite along orbit plane
			ColumnVector unit_Node_ascending = crossproduct(unit_Z, unit_Nop);
			unit_Node_ascending = unit_Node_ascending / unit_Node_ascending.norm_Frobenius();
			double fx = dotproduct(unit_Node_ascending, unit_sat);
			ColumnVector unit_tmp = crossproduct(unit_Node_ascending, unit_sat);
			double fz = dotproduct(unit_tmp, unit_Nop);
			double u = atan2(fz, fx);

			//u angel from the ascending node to the Sun's projection along orbit plane
			// unit_terminator: unit vector of the twilight terminator
			ColumnVector unit_terminator = crossproduct(unit_sun, unit_Nop);
			unit_terminator = unit_terminator / unit_terminator.norm_Frobenius();
			// unit_projected_sun: unit vector of Sun's projection on the orbital plane
			ColumnVector unit_projected_sun = crossproduct(unit_Nop, unit_terminator);
			unit_projected_sun = unit_projected_sun / unit_projected_sun.norm_Frobenius();
			fx = dotproduct(unit_Node_ascending, unit_projected_sun);
			unit_tmp = crossproduct(unit_Node_ascending, unit_projected_sun);
			fz = dotproduct(unit_tmp, unit_Nop);
			double u_s = atan2(fz, fx);

			// du in ECOM2 model
			double du = u - u_s;

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

			double cos_du = cos(du);
			double sin_du = sin(du);
			double cos_2du = cos(2 * du);
			double sin_2du = sin(2 * du);
			double cos_4du = cos(4*du);
			double sin_4du = sin(4*du);

			// calcluate acc in DYB frame
			DYXu(1) = _ecom2_pars[ICS_PAR::ECOM2_D0] + (_ecom2_pars[ICS_PAR::ECOM2_D2c] * cos_2du + _ecom2_pars[ICS_PAR::ECOM2_D2s] * sin_2du + _ecom2_pars[ICS_PAR::ECOM2_D4c] * cos_4du + _ecom2_pars[ICS_PAR::ECOM2_D4s] * sin_4du) * lambda;
			DYXu(2) = _ecom2_pars[ICS_PAR::ECOM2_Y0] + (0.0) * lambda;
			DYXu(3) = _ecom2_pars[ICS_PAR::ECOM2_X0] + (_ecom2_pars[ICS_PAR::ECOM2_X1c] * cos_du + _ecom2_pars[ICS_PAR::ECOM2_X1s] * sin_du) * lambda;

			// acc in crs
			DYXf = (unit_D * DYXu(1) + unit_Y * DYXu(2) - unit_X * DYXu(3)) * (*(arg->shadow_factor));

			// calculate acc, km/s^2
			_acc = DYXf * dist_factor * 1e-9;

			// Calculate the partial derivative of DYX, parameter sequence: d y b 
			// reset zeros for safe
			_partial = 0.0;
			for (int i = 1; i <= 3; i++)
			{
				if (_force_model == F_SOLAR_ECOM2)
				{
					_partial(1, i) = unit_D(i) * (*(arg->shadow_factor));
					_partial(2, i) = unit_D(i) * cos_2du * lambda;
					_partial(3, i) = unit_D(i) * sin_2du * lambda;
					_partial(4, i) = unit_D(i) * cos_4du * lambda;
					_partial(5, i) = unit_D(i) * sin_4du * lambda;

					_partial(6, i) = unit_Y(i) * (*(arg->shadow_factor));
					_partial(7, i) = -unit_X(i) * (*(arg->shadow_factor));
					_partial(8, i) = -unit_X(i) *  cos_du * lambda;
					_partial(9, i) = -unit_X(i) *  sin_du * lambda;
				}
			}
			_partial = _partial * dist_factor * 1e-9;


#ifdef DEBUG
			cout << "==> The force model : ECOM " << endl
				<< "==> And Acc in DYB is : " << endl
				<< "==> D : " << setw(30) << setprecision(15) << scientific << DYXu(1) << endl
				<< "==> Y : " << setw(30) << setprecision(15) << scientific << DYXu(2) << endl
				<< "==> B : " << setw(30) << setprecision(15) << scientific << DYXu(3) << endl
				<< "==> And Acc in CRS is : " << endl
				<< "==> X : " << setw(30) << setprecision(15) << scientific << _acc(1) << endl
				<< "==> Y : " << setw(30) << setprecision(15) << scientific << _acc(2) << endl
				<< "==> Z : " << setw(30) << setprecision(15) << scientific << _acc(3) << endl
				<< "==> And parital is :" << endl;
			if (_force_model == F_SOLAR_ECOM2)
			{
				cout << "==> Partial " << setw(30) << setprecision(15) << scientific << _partial(1, 1) << _partial(1, 2) << _partial(1, 3) << endl
				     << "==> Partial " << setw(30) << setprecision(15) << scientific << _partial(2, 1) << _partial(2, 2) << _partial(2, 3) << endl
				     << "==> Partial " << setw(30) << setprecision(15) << scientific << _partial(3, 1) << _partial(3, 2) << _partial(3, 3) << endl
				     << "==> Partial " << setw(30) << setprecision(15) << scientific << _partial(4, 1) << _partial(4, 2) << _partial(4, 3) << endl
				     << "==> Partial " << setw(30) << setprecision(15) << scientific << _partial(5, 1) << _partial(5, 2) << _partial(5, 3) << endl
				     << "==> Partial " << setw(30) << setprecision(15) << scientific << _partial(6, 1) << _partial(6, 2) << _partial(6, 3) << endl
				     << "==> Partial " << setw(30) << setprecision(15) << scientific << _partial(7, 1) << _partial(7, 2) << _partial(7, 3) << endl
				     << "==> Partial " << setw(30) << setprecision(15) << scientific << _partial(8, 1) << _partial(8, 2) << _partial(8, 3) << endl
				     << "==> Partial " << setw(30) << setprecision(15) << scientific << _partial(9, 1) << _partial(9, 2) << _partial(9, 3) << endl;
			}
			else
			{
				cout << "==> Partial " << setw(30) << setprecision(15) << scientific << _partial(1, 1) << _partial(1, 2) << _partial(1, 3) << endl
				     << "==> Partial " << setw(30) << setprecision(15) << scientific << _partial(2, 1) << _partial(2, 2) << _partial(2, 3) << endl
				     << "==> Partial " << setw(30) << setprecision(15) << scientific << _partial(3, 1) << _partial(3, 2) << _partial(3, 3) << endl
				     << "==> Partial " << setw(30) << setprecision(15) << scientific << _partial(4, 1) << _partial(4, 2) << _partial(4, 3) << endl
				     << "==> Partial " << setw(30) << setprecision(15) << scientific << _partial(5, 1) << _partial(5, 2) << _partial(5, 3) << endl;
			}
#endif // DEBUG

			return true;
		}
		catch (...)
		{
			write_log(_glog, 0, "ERROR : " + function_name + " ==> Throw while calculating");
			return false;
		}
	}

}

