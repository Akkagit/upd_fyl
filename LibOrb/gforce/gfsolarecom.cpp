/**
 *
 * @verbatim
	 History
	  -1.0 Bo Wong			 2018-12-23 creat the file.
	  -1.1 Jiande Huang      2019-02-20 edit for final version
	  -1.2 Jiande Huang      2019-03-14 Adding Doxygen Style Code Remarks
	  -1.3 Jiande Huang      2019-04-01 Adding exception throwing mechanism for all functions
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gfsolarecom.cpp
 * @brief		ECOM force model, which used for getting solar radiation.
 *
 * @author      Bo Wong, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 *
 */

#include <iostream>
#include "gutils/gconst.h"
#include "gutils/ginfolog.h"
#include "gforce/gfarg_cal.h"
#include "gforce/gfcommon.h"
#include "gforce/gfsolarecom.h"


namespace great
{
	 t_gfsolar_ecom::t_gfsolar_ecom(t_glog * log) : t_gfsolar(log)
	 {
		if (!log)
		{
			write_log(log, 0, "t_gfsolar_ecom::t_gfsolar_ecom",
				"ERROR : the log point is empty");
		}
	}

	 t_gfsolar_ecom::t_gfsolar_ecom(t_init_arg * arg, t_glog * log) : t_gfsolar(log)
	 {
		 string function_name = "t_gfsolar_ecom::t_gfsolar_ecom";
		 try
		 {
			 t_init_solar_ecom* ecom = NULL;
			 ecom = dynamic_cast<t_init_solar_ecom*>(arg);
			 if (!ecom)
			 {
				 write_log(log, 0, "ERROR : " + function_name + " ==> init failed");
				 write_log(log, 0, "ERROR : " + function_name + " ==> no init argument");
			 }
			 else
			 {
				 if (ecom->solar_pars.size() == 5 && (arg->force_model == F_SOLAR_ECOM_5 || arg->force_model == F_SOLAR_DEF))
				 {
					 _force_model = arg->force_model;
					 // get params' value from ics file according to epoch
					 _ecom_pars[ICS_PAR::ECOM_D0] = (ecom->solar_pars.count(ICS_PAR::ECOM_D0) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_D0] : 0;
					 _ecom_pars[ICS_PAR::ECOM_Y0] = (ecom->solar_pars.count(ICS_PAR::ECOM_Y0) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_Y0] : 0;
					 _ecom_pars[ICS_PAR::ECOM_X0] = (ecom->solar_pars.count(ICS_PAR::ECOM_X0) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_X0] : 0;
					 _ecom_pars[ICS_PAR::ECOM_Xc] = (ecom->solar_pars.count(ICS_PAR::ECOM_Xc) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_Xc] : 0;
					 _ecom_pars[ICS_PAR::ECOM_Xs] = (ecom->solar_pars.count(ICS_PAR::ECOM_Xs) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_Xs] : 0;

					 _ecom_pars[ICS_PAR::ECOM_Dc] = 0;
					 _ecom_pars[ICS_PAR::ECOM_Ds] = 0;
					 _ecom_pars[ICS_PAR::ECOM_Yc] = 0;
					 _ecom_pars[ICS_PAR::ECOM_Ys] = 0;

					 _partial.resize(5, 3);

				 }
				 else if (ecom->solar_pars.size() == 9 && arg->force_model == F_SOLAR_ECOM_9)
				 {
					 _force_model = arg->force_model;
					 // get params' value from ics file according to epoch
					 _ecom_pars[ICS_PAR::ECOM_D0] = (ecom->solar_pars.count(ICS_PAR::ECOM_D0) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_D0] : 0;
					 _ecom_pars[ICS_PAR::ECOM_Dc] = (ecom->solar_pars.count(ICS_PAR::ECOM_Dc) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_Dc] : 0;
					 _ecom_pars[ICS_PAR::ECOM_Ds] = (ecom->solar_pars.count(ICS_PAR::ECOM_Ds) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_Ds] : 0;
					 _ecom_pars[ICS_PAR::ECOM_Y0] = (ecom->solar_pars.count(ICS_PAR::ECOM_Y0) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_Y0] : 0;
					 _ecom_pars[ICS_PAR::ECOM_Yc] = (ecom->solar_pars.count(ICS_PAR::ECOM_Yc) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_Yc] : 0;
					 _ecom_pars[ICS_PAR::ECOM_Ys] = (ecom->solar_pars.count(ICS_PAR::ECOM_Ys) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_Ys] : 0;
					 _ecom_pars[ICS_PAR::ECOM_X0] = (ecom->solar_pars.count(ICS_PAR::ECOM_X0) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_X0] : 0;
					 _ecom_pars[ICS_PAR::ECOM_Xc] = (ecom->solar_pars.count(ICS_PAR::ECOM_Xc) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_Xc] : 0;
					 _ecom_pars[ICS_PAR::ECOM_Xs] = (ecom->solar_pars.count(ICS_PAR::ECOM_Xs) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_Xs] : 0;
					 _partial.resize(9, 3);
				 }
				 else if (ecom->solar_pars.size() == 3 && arg->force_model == F_SOLAR_ECOM_3)
				 {
					 _force_model = arg->force_model;
					 // get params' value from ics file according to epoch
					 _ecom_pars[ICS_PAR::ECOM_D0] = (ecom->solar_pars.count(ICS_PAR::ECOM_D0) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_D0] : 0;
					 _ecom_pars[ICS_PAR::ECOM_Y0] = (ecom->solar_pars.count(ICS_PAR::ECOM_Y0) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_Y0] : 0;
					 _ecom_pars[ICS_PAR::ECOM_X0] = (ecom->solar_pars.count(ICS_PAR::ECOM_X0) > 0) ? ecom->solar_pars[ICS_PAR::ECOM_X0] : 0;
					 
					 _ecom_pars[ICS_PAR::ECOM_Dc] = 0;
					 _ecom_pars[ICS_PAR::ECOM_Ds] = 0;
					 _ecom_pars[ICS_PAR::ECOM_Yc] = 0;
					 _ecom_pars[ICS_PAR::ECOM_Ys] = 0;
					 _ecom_pars[ICS_PAR::ECOM_Xc] = 0;
					 _ecom_pars[ICS_PAR::ECOM_Xs] = 0;
					 
					 _partial.resize(3, 3);
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
				 << " ECOM_D0 " << setw(30) << setprecision(15) << scientific << _ecom_pars[ECOM_D0] << endl
				 << " ECOM_Dc " << setw(30) << setprecision(15) << scientific << _ecom_pars[ECOM_Dc] << endl
				 << " ECOM_Ds " << setw(30) << setprecision(15) << scientific << _ecom_pars[ECOM_Ds] << endl
				 << " ECOM_Y0 " << setw(30) << setprecision(15) << scientific << _ecom_pars[ECOM_Y0] << endl
				 << " ECOM_Yc " << setw(30) << setprecision(15) << scientific << _ecom_pars[ECOM_Yc] << endl
				 << " ECOM_Ys " << setw(30) << setprecision(15) << scientific << _ecom_pars[ECOM_Ys] << endl
				 << " ECOM_X0 " << setw(30) << setprecision(15) << scientific << _ecom_pars[ECOM_X0] << endl
				 << " ECOM_Xc " << setw(30) << setprecision(15) << scientific << _ecom_pars[ECOM_Xc] << endl
				 << " ECOM_Xs " << setw(30) << setprecision(15) << scientific << _ecom_pars[ECOM_Xs] << endl;
#endif // DEBUG

		 }
		 catch (...)
		 {
			 write_log(log, 0, "ERROR : " + function_name + " ==> throw exception while constructing ");
		 }
	 }

	/** @todo Functional functions still need to be improved */
	void t_gfsolar_ecom::show()
	{
		cout << "solar model : ECOM" << endl;
	}

	/** Functional functions still need to be improved */
	bool t_gfsolar_ecom::check()
	{
		write_log(_glog, 4, "check some important pars in model : ECOM");
		if (_ecom_pars.size() == 0)
		{
			write_log(_glog, 0, "ERROR : no ics data in ECOM model!");
			return false;
		}
		return true;
		// TODO
	}


	bool t_gfsolar_ecom::calculate(t_gforce_calarg* argument)
	{
		string function_name = "t_gfsolar_ecom::calculate";
		// 判断计算的充分条件，并获取参数
		t_gfcalarg_solar_ecom* arg = NULL;
		arg = dynamic_cast<t_gfcalarg_solar_ecom*>(argument);
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

			double cos_u = cos(u);
			double sin_u = sin(u);

			// calcluate acc in sat fix
			DYXu(1) = _ecom_pars[ICS_PAR::ECOM_D0] + (_ecom_pars[ICS_PAR::ECOM_Dc] * cos_u + _ecom_pars[ICS_PAR::ECOM_Ds] * sin_u) * lambda;
			DYXu(2) = _ecom_pars[ICS_PAR::ECOM_Y0] + (_ecom_pars[ICS_PAR::ECOM_Yc] * cos_u + _ecom_pars[ICS_PAR::ECOM_Ys] * sin_u) * lambda;
			DYXu(3) = _ecom_pars[ICS_PAR::ECOM_X0] + (_ecom_pars[ICS_PAR::ECOM_Xc] * cos_u + _ecom_pars[ICS_PAR::ECOM_Xs] * sin_u) * lambda;

			// acc in crs
			DYXf = (unit_D * DYXu(1) + unit_Y * DYXu(2) - unit_X * DYXu(3)) * (*(arg->shadow_factor));

			// calculate acc, km/s^2
			_acc = DYXf * dist_factor * 1e-9;

			// Calculate the partial derivative of DYX, parameter sequence: d y b 
			// reset zeros for safe
			_partial = 0.0;
			for (int i = 1; i <= 3; i++)
			{
				if (_force_model == F_SOLAR_ECOM_5)
				{
					_partial(1, i) = unit_D(i) * (*(arg->shadow_factor));
					_partial(2, i) = unit_Y(i) * (*(arg->shadow_factor));
					_partial(3, i) = -unit_X(i) * (*(arg->shadow_factor));
					_partial(4, i) = -unit_X(i) *  cos_u * lambda;
					_partial(5, i) = -unit_X(i) *  sin_u * lambda;
				}
				else if (_force_model == F_SOLAR_ECOM_9)
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
				else if (_force_model == F_SOLAR_ECOM_3)
				{
					_partial(1, i) = unit_D(i) * (*(arg->shadow_factor));
					_partial(2, i) = unit_Y(i) * (*(arg->shadow_factor));
					_partial(3, i) = -unit_X(i) * (*(arg->shadow_factor));
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
			if (_force_model == F_SOLAR_ECOM_9)
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

