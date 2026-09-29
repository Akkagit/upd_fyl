/**
 *
 * @verbatim
	 History
	  -1.0 Jiande Huang   2019-01-06 creat the file.
	  -1.1 Jie Li        2019-03-31 add code in doxygen style
	  -1.2 Jiande Huang   2019-04-02 Adding abnormal throwing mechanism of functions
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gfcommon.cpp
 * @brief		Some common defination and functions for orbit integration process.
 *
 * @author      Jiande Huang, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 */

#include "gforce/gfcommon.h"
#include "gutils/gconst.h"
#include "gutils/gconst_oi.h"
#include "gdata/gnavde.h"
#include <algorithm>

namespace great
{

	double get_shadow_factor(const ColumnVector& sat_pos, const ColumnVector& sun_pos, const ColumnVector& moon_pos)
	{
		// reset lambda(shadow factor) to zeros
		double lambda = 0.0;

		try
		{
			// position vector refer to satellite
			ColumnVector sat2sun = sun_pos - sat_pos;
			ColumnVector sat2moon = moon_pos - sat_pos;
			ColumnVector sat2earth = -sat_pos;

			double rsat2sun = sqrt(dotproduct(sat2sun, sat2sun));
			double rsat2moon = sqrt(dotproduct(sat2moon, sat2moon));
			double rsat2earth = sqrt(dotproduct(sat2earth, sat2earth));

			// sight distances between solar and lunar, solar and earth
			// see from satellite
			double d01 = acos(dotproduct(sat2sun, sat2moon) / rsat2sun / rsat2moon);
			double d02 = acos(dotproduct(sat2sun, sat2earth) / rsat2sun / rsat2earth);

			// sight radius of solar, lunar and earth see from satellite
			double r0 = SUN_R / rsat2sun;
			double r1 = MOON_R / rsat2moon;
			double r2 = EARTH_R / rsat2earth;

			// check the eclipse( d,s02,s123 are in unit PI )
			double s01 = get_shadow_area(r0, r1, d01);
			double s02 = get_shadow_area(r0, r2, d02);

			// eclipse factor
			lambda = 1.0 - (s01 + s02) / r0 / r0 / G_PI;
			//lambda = (lambda < 0.0) ? 0.0 : lambda;
			//if (lambda < 0)
			//{
			//	lambda = 0;
			//}
			return (lambda < 0.0) ? 0.0 : lambda;
		}
		catch (...)
		{
			cout << "ERROR : throw exception while calculating shadow_factor, we return 0.0 as default value" << endl;
			cout << "NOTE : This may bring some error" << endl;
			return 0.0;
		}
	}


	double get_shadow_area(double a_radius, double b_radius, double ab_distance)
	{
		double lam = 0.0;
		double alpha_a = 0.0;
		double alpha_b = 0.0;

		try
		{
			if (ab_distance >= (a_radius + b_radius))
			{
				lam = 0;
			}
			else if (ab_distance < fabs(a_radius - b_radius))
			{
				if (a_radius >= b_radius)
				{
					lam = b_radius * b_radius* G_PI;
				}
				else
				{
					lam = a_radius * a_radius * G_PI;
				}
			}
			else
			{
				alpha_a = acos((a_radius*a_radius + ab_distance * ab_distance - b_radius * b_radius) / 2 / a_radius / ab_distance);
				alpha_b = acos((b_radius*b_radius + ab_distance * ab_distance - a_radius * a_radius) / 2 / b_radius / ab_distance);
				lam = a_radius * a_radius*(alpha_a - sin(alpha_a)*cos(alpha_a)) + b_radius * b_radius*(alpha_b - sin(alpha_b)*cos(alpha_b));
			}
			return lam;
		}
		catch (...)
		{
			cout << "ERROR : throw exception while calculating shadow_area, we return 0.0 as default value" << endl;
			cout << "NOTE : This may bring some error" << endl;
			return 0.0;
		}
	}

	void change_time_sys(CHANGE_TYPE type, double inp_time, double& out_time)
	{
		switch (type)
		{
		case GPS2TDT:
			out_time = inp_time + 51.184 / 86400.0;
			break;
		default:
			out_time = inp_time;
		}
		return;
	}


	void change_time_sys(CHANGE_TYPE type, t_gtime inp_time, t_gtime & out_time)
	{
		switch (type)
		{
		case GPS2TDT:
			out_time = inp_time + 51.184;
			out_time.tsys(out_time.str2tsys("TT"));
			break;
		default:
			out_time = inp_time;
		}
	}


	ICS_PAR str2ics(string par_name)
	{
		if (par_name == "PXSAT") return ICS_PAR::PXSAT;
		if (par_name == "PYSAT") return ICS_PAR::PYSAT;
		if (par_name == "PZSAT") return ICS_PAR::PZSAT;
		if (par_name == "VXSAT") return ICS_PAR::VXSAT;
		if (par_name == "VYSAT") return ICS_PAR::VYSAT;
		if (par_name == "VZSAT") return ICS_PAR::VZSAT;

		if (par_name == "D0" || par_name == "ECOM_D0" || par_name == "Kd_BERN")   return ICS_PAR::ECOM_D0;
		if (par_name == "Dc" || par_name == "ECOM_Dc" || par_name == "Kdc1_BERN") return ICS_PAR::ECOM_Dc;
		if (par_name == "Ds" || par_name == "ECOM_Ds" || par_name == "Kds1_BERN") return ICS_PAR::ECOM_Ds;
		if (par_name == "Y0" || par_name == "ECOM_Y0" || par_name == "Ky_BERN")   return ICS_PAR::ECOM_Y0;
		if (par_name == "Yc" || par_name == "ECOM_Yc" || par_name == "Kyc1_BERN") return ICS_PAR::ECOM_Yc;
		if (par_name == "Ys" || par_name == "ECOM_Ys" || par_name == "Kys1_BERN") return ICS_PAR::ECOM_Ys;
		if (par_name == "X0" || par_name == "ECOM_X0" || par_name == "Kb_BERN")   return ICS_PAR::ECOM_X0;
		if (par_name == "Xc" || par_name == "ECOM_Xc" || par_name == "Kbc1_BERN") return ICS_PAR::ECOM_Xc;
		if (par_name == "Xs" || par_name == "ECOM_Xs" || par_name == "Kbs1_BERN") return ICS_PAR::ECOM_Xs;

		if (par_name == "ECOM2_D0")   return ICS_PAR::ECOM2_D0;
		if (par_name == "ECOM2_D2c" ) return ICS_PAR::ECOM2_D2c;
		if (par_name == "ECOM2_D2s" ) return ICS_PAR::ECOM2_D2s;
		if (par_name == "ECOM2_D4c")  return ICS_PAR::ECOM2_D4c;
		if (par_name == "ECOM2_D4s")  return ICS_PAR::ECOM2_D4s;
		if (par_name == "ECOM2_Y0")   return ICS_PAR::ECOM2_Y0;
		if (par_name == "ECOM2_X0")   return ICS_PAR::ECOM2_X0;
		if (par_name == "ECOM2_X1c" ) return ICS_PAR::ECOM2_X1c;
		if (par_name == "ECOM2_X1s" ) return ICS_PAR::ECOM2_X1s;

		return ICS_PAR::NO_DEF;
	}


	ICS_PAR str2icsleo(string par_name)
	{
		if (par_name == "PXSAT") return ICS_PAR::PXSAT;
		if (par_name == "PYSAT") return ICS_PAR::PYSAT;
		if (par_name == "PZSAT") return ICS_PAR::PZSAT;
		if (par_name == "VXSAT") return ICS_PAR::VXSAT;
		if (par_name == "VYSAT") return ICS_PAR::VYSAT;
		if (par_name == "VZSAT") return ICS_PAR::VZSAT;

		if (par_name == "EMP_Ac") return ICS_PAR::EMP_Ac;
		if (par_name == "EMP_Aa") return ICS_PAR::EMP_Aa;
		if (par_name == "EMP_Ar") return ICS_PAR::EMP_Ar;
		if (par_name == "EMP_Bc") return ICS_PAR::EMP_Bc;
		if (par_name == "EMP_Ba") return ICS_PAR::EMP_Ba;
		if (par_name == "EMP_Br") return ICS_PAR::EMP_Br;

		if (par_name == "EMP_Sc1") return ICS_PAR::EMP_Sc1;
		if (par_name == "EMP_Cc1") return ICS_PAR::EMP_Cc1;
		if (par_name == "EMP_Sa1") return ICS_PAR::EMP_Sa1;
		if (par_name == "EMP_Ca1") return ICS_PAR::EMP_Ca1;
		if (par_name == "EMP_Sr1") return ICS_PAR::EMP_Sr1;
		if (par_name == "EMP_Cr1") return ICS_PAR::EMP_Cr1;

		if (par_name == "EMP_Sc2") return ICS_PAR::EMP_Sc2;
		if (par_name == "EMP_Cc2") return ICS_PAR::EMP_Cc2;
		if (par_name == "EMP_Sa2") return ICS_PAR::EMP_Sa2;
		if (par_name == "EMP_Ca2") return ICS_PAR::EMP_Ca2;
		if (par_name == "EMP_Sr2") return ICS_PAR::EMP_Sr2;
		if (par_name == "EMP_Cr2") return ICS_PAR::EMP_Cr2;

		if (par_name == "SR_scale") return ICS_PAR::SR_scale;
		if (par_name == "DRAG_c") return ICS_PAR::DRAG_c;
		
		return ICS_PAR::NO_DEF;
	}

	/**
	* @brief FORCE type to string
	* @param[in]   force	force model name in FORCE
	* @return	force model name in string type
	*/
	string force2str(FORCE force)
	{
		switch (force)
		{
		case GRAVITY:
		{
			return "gravity";
		}
		case RELATIVITY:
		{
			return "relativity";
		}
		case POINT_MASS:
		{
			return "point_mass";
		}
		case SOLAR_RADI:
		{
			return "solar_radiation";
		}
		case EARTH_RADI:
		{
			return "earth_radiation";
		}
		case CUSTOMER:
		{
			return "customer";
		}
		case ATMOSPHERE:
		{
			return "atmospher";
		}
		default:
		{
			return "No_defination";
		}
		}
	}


	/**
	* @brief ICS_PAR type to string
	* @param[in]   par_name		parameter name in ICS_PAR
	* @return	parameter namee in string type
	*/
	string ics2str(ICS_PAR par_name)
	{
		switch (par_name)
		{
		case ICS_PAR::PXSAT: return "PXSAT";
		case ICS_PAR::PYSAT: return "PYSAT";
		case ICS_PAR::PZSAT: return "PZSAT";
		case ICS_PAR::VXSAT: return "VXSAT";
		case ICS_PAR::VYSAT: return "VYSAT";
		case ICS_PAR::VZSAT: return "VZSAT";
		case ICS_PAR::ECOM_D0: return "D0";
		case ICS_PAR::ECOM_Dc: return "Dc";
		case ICS_PAR::ECOM_Ds: return "Ds";
		case ICS_PAR::ECOM_Y0: return "Y0";
		case ICS_PAR::ECOM_Yc: return "Yc";
		case ICS_PAR::ECOM_Ys: return "Ys";
		case ICS_PAR::ECOM_X0: return "X0";
		case ICS_PAR::ECOM_Xc: return "Xc";
		case ICS_PAR::ECOM_Xs: return "Xs";
		case ICS_PAR::ECOM2_D0:  return "ECOM2_D0";
		case ICS_PAR::ECOM2_D2c: return "ECOM2_D2c";
		case ICS_PAR::ECOM2_D2s: return "ECOM2_D2s";
		case ICS_PAR::ECOM2_D4c: return "ECOM2_D4c";
		case ICS_PAR::ECOM2_D4s: return "ECOM2_D4s";
		case ICS_PAR::ECOM2_Y0:  return "ECOM2_Y0";
		case ICS_PAR::ECOM2_X0:  return "ECOM2_X0";
		case ICS_PAR::ECOM2_X1c: return "ECOM2_X1c";
		case ICS_PAR::ECOM2_X1s: return "ECOM2_X1s";
		case ICS_PAR::EMP_Sc1: return "EMP_Sc1";
		case ICS_PAR::EMP_Cc1: return "EMP_Cc1";
		case ICS_PAR::EMP_Sa1: return "EMP_Sa1";
		case ICS_PAR::EMP_Ca1: return "EMP_Ca1";
		case ICS_PAR::EMP_Sr1: return "EMP_Sr1";
		case ICS_PAR::EMP_Cr1: return "EMP_Cr1";
					
		case ICS_PAR::EMP_Sc2: return "EMP_Sc2";
		case ICS_PAR::EMP_Cc2: return "EMP_Cc2";
		case ICS_PAR::EMP_Sa2: return "EMP_Sa2";
		case ICS_PAR::EMP_Ca2: return "EMP_Ca2";
		case ICS_PAR::EMP_Sr2: return "EMP_Sr2";
		case ICS_PAR::EMP_Cr2: return "EMP_Cr2";

		case ICS_PAR::EMP_Ac:return "EMP_Ac";
		case ICS_PAR::EMP_Bc:return "EMP_Bc";
		case ICS_PAR::EMP_Aa:return "EMP_Aa";
		case ICS_PAR::EMP_Ba:return "EMP_Ba";
		case ICS_PAR::EMP_Ar:return "EMP_Ar";
		case ICS_PAR::EMP_Br:return "EMP_Br";
	
		case ICS_PAR::SR_scale: return "SR_scale";
		case ICS_PAR::DRAG_c:   return "DRAG_c";
		
		default:
			return "NO_DEF";
		}
	}

	/**
	* @brief string type to bool
	* @param[in]   s	bool expressed by string
	* @return	bool
	*/
	bool str2bool(const string& s)
	{
		string tmp = s;
		transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
		if (tmp == "true" || tmp == "yes")
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	/**
	* @brief GSYS to string type
	* @param[in]   sys	system expressed by GSYS
	* @return	system expressed by string
	*/
	string gsys2str(GSYS sys)
	{

		gtrace("t_gsys::gsys2str");

		// PROBLEM WITH STATIC FUNCTION
		//  boost::mutex::scoped_lock lock(_mutex);

		switch (sys) {
		case GPS: return "GPS";
		case GLO: return "GLO";
		case GAL: return "GAL";
		case BDS: return "BDS";
		case SBS: return "SBS";
		case QZS: return "QZS";
		case IRN: return "IRN";
		case LEO: return "LEO";
		case GNS: return "GNS";

		default: { cout << "*** warning: unknown GNSS system!\n"; cout.flush(); }

				 // cannot be while static function !!!
				 // if( _log ) _log->comment(0,"gobsgnss","warning: unknown frequency code ");
				 //     else 
		}

		return "GNS";
	}

	/**
	* @brief string type to PLANET
	* @param[in]   planet	planet name expressed by string
	* @return	planet name expressed by PLANET(enum)
	*/
	int str2planet(const string& planet)
	{
		string tmp = planet;
		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
		     if (tmp == "MERCURY" || tmp == "MERC") return PLANET_MERC;
		else if (tmp == "VENUS"   || tmp == "VENU") return PLANET_VENU;
		else if (tmp == "EARTH"   || tmp == "EART") return PLANET_EART;
		else if (tmp == "MARS")                     return PLANET_MARS;
		else if (tmp == "JUPITER" || tmp == "JUPI") return PLANET_JUPI;
		else if (tmp == "SATURN"  || tmp == "SATU") return PLANET_SATU;
		else if (tmp == "URANUS"  || tmp == "URAN") return PLANET_URAN;
		else if (tmp == "NEPTUNE" || tmp == "NEPT") return PLANET_NEPT;
		else if (tmp == "PLUTO"   || tmp == "PLUT") return PLANET_PLUT;
		else if (tmp == "MOON")                     return PLANET_MOON;
		else if (tmp == "SUN")                      return PLANET_SUN;
		else if (tmp == "SOLAR_SYSTEM")             return SOLAR_SYSTEM;
		else if (tmp == "EARTH_MOON")               return EARTH_MOON;
		else { cout << "unknow planet" << endl; return -1; }
	}


	void get_mean_pole(double mjd, double & xpm, double & ypm)
	{
		// 51544 is the mjd2000
		// 55197 is the mjd2010
		double dt = (mjd - 51544) / 365.25;

		if (mjd > 55197)
		{
			xpm = 23.513 + 7.6141 * dt;
			ypm = 358.891 - 0.6287 * dt;
		}
		else
		{
			xpm = 55.974 + (1.8243 + (0.18413 + 0.007024 * dt) * dt) * dt;
			ypm = 346.346 + (1.7896 - (0.10729 - 0.000908 * dt) * dt) * dt;
		}

		//to seconds
		xpm = xpm * (1e-3);
		ypm = ypm * (1e-3);
	}
}

