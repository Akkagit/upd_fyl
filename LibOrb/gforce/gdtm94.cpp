/**
 *
 * @verbatim
	 History
	  -1.0 Yujie Qin and Bo Wong  2019-01-06 Creat the file.
	  -1.1 Zhang Hongmin          2019-03-04 update.
	  -1.2 Yujie Qin              2019-03-05 edit for final version.
	  -1.3 Yujie Qin              2019-03-15 Code comments are added in Doxygen format.
	  -1.4 Jiande Huang            2019-04-01 Adding abnormal throwing mechanism of functions
   @endverbatim
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		    gdtm94.cpp
 * @brief		    The class for dtm94 model, serving for atmospheric drag model by providing the density.
*
 * @author          Yujie Qin, Wuhan University
 * @version		    1.0.0
 * @date		    2019-01-06
 *
 */

#include "gforce/gdtm94.h"
#include "gutils/ginfolog.h"

using namespace std;
namespace great
{
	/** constructor. */
	t_gfatmos_dtm94::t_gfatmos_dtm94(t_glog* log) : t_gforce(log)
	{
		string fuction_name = "t_gfatmos_dtm94::t_gfatmos_dtm94";
		if (!log)
		{
			write_log(log, 0, fuction_name,
				"ERROR : the log point is empty");
		}
	}

	/** destructor. */
	t_gfatmos_dtm94::~t_gfatmos_dtm94()
	{

	}

	/** Calculate Legendre polynomial. */
	bool t_gfatmos_dtm94::cal_Pnm()
	{
		string function_name = "t_gfatmos_dtm94::calPnm";
		try
		{
			double sinlat = sin(_sat_geo_crd(1));
			double sinlat2 = pow(sinlat, 2);
			double sinlat4 = pow(sinlat2, 2);
			double coslat = cos(_sat_geo_crd(1));
			double coslat2 = pow(coslat, 2);
			_pnm[0] = sinlat;
			_pnm[1] = 0.5*(3.0*sinlat2 - 1.0);
			_pnm[2] = 0.5*(5.0*sinlat2 - 3.0)*sinlat;
			_pnm[3] = 0.125*(35.0*sinlat4 - 30.0*sinlat2 + 3.0);
			_pnm[4] = 0.125*(63.0*sinlat4 - 70.0*sinlat2 + 15.0)*sinlat;

			_pnm[5] = coslat;
			_pnm[6] = 3.0*sinlat*coslat;
			_pnm[7] = 0.5*(15.0*sinlat2 - 3.0)*coslat;
			_pnm[8] = 0.125*(315.0*sinlat4 - 210.0*sinlat2 + 15.0)*coslat;

			_pnm[9] = 3.0*coslat2;
			_pnm[10] = 15.0*sinlat*coslat2;
			_pnm[11] = 15.0*coslat2*coslat;
			return true;
		}
		catch (...)
		{
			write_log(_glog, 0, function_name + "ERROR : throw exception while calculating Pnm");
			write_log(_glog, 0, "Maybe you should check the value : _geosat(1)");
			return false;
		}
	}


	/** calculate the Gl, Gl:Empirical Function.
	* @param[in] gas the gas to calculate.
	* @param[in] isTemp use index counting or not.
	* @param[in] isHorhe use solar_activity or not.
	*/
	bool t_gfatmos_dtm94::cal_GL(const double* gas, bool isTemp, bool isHorhe)
	{
		string function_name = "t_gfatmos_dtm94::calGl";
		try
		{
			t_gtime tdoy;
			tdoy.from_dmjd(_epoch);
			double doy = tdoy.doy() + tdoy.dmjd() - tdoy.mjd();

			double scale = 0.0;
			double omega = 3.14159265 / 12.0;
			double omega_U = 2.0* 3.14159265 / 365.0;
			double t = fmod((_sat_geo_crd(2) - _sun_geo_crd(2)) / 3.14159265 * 12.0 + 12.0, 24.0);

			t = (t < 0) ? t + 24.0 : t;
			
			double solar_activity = gas[3] * (_solar[0] - _solar[1]) +
				gas[4] * pow((_solar[0] - _solar[1]), 2) +
				gas[5] * (_solar[1] - 150.0) +
				gas[37] * pow((_solar[1] - 150.0), 2);

			scale = (isHorhe) ? 1.0 : 1.0 + solar_activity;
			

			_GL = gas[1] * _pnm[1] + gas[2] * _pnm[3]
				+ gas[36] * _pnm[0] + solar_activity
				+ (gas[6] + gas[7] * _pnm[1])*_kp + scale * ((gas[8] + gas[9] * _pnm[1]) * cos(omega_U * (doy - gas[10]))
					+ (gas[11] + gas[12] * _pnm[1])*cos(2 * omega_U * (doy - gas[13]))
					+ (gas[14] * _pnm[0] + gas[15] * _pnm[2] + gas[16] * _pnm[4]) * cos(omega_U * (doy - gas[17])) + gas[18] * _pnm[0] * cos(2 * omega_U* (doy - gas[19]))
					+ (gas[20] * _pnm[5] + gas[21] * _pnm[7] + gas[22] * _pnm[8] + (gas[23] * _pnm[5] + gas[24] * _pnm[6]) * cos(omega_U* (doy - gas[17]))) * cos(omega * t)
					+ (gas[25] * _pnm[5] + gas[26] * _pnm[7] + gas[27] * _pnm[8] + (gas[28] * _pnm[5] + gas[29] * _pnm[6]) * cos(omega_U * (doy - gas[17]))) * sin(omega * t)
					+ (gas[30] * _pnm[9] + gas[31] * _pnm[10] * cos(omega_U * (doy - gas[17]))) * cos(2 * omega * t)
					+ (gas[32] * _pnm[9] + gas[33] * _pnm[10] * cos(omega_U * (doy - gas[17]))) * sin(2 * omega * t) + gas[34] * _pnm[11] * cos(3 * omega * t) + gas[35] * _pnm[11] * sin(3 * omega * t));

			_GL = (isTemp) ? _GL + gas[38] * exp(_kp) : _GL + gas[38] * _kp * _kp;
			

			return true;
		}
		catch (...)
		{
			write_log(_glog, 0, "ERROR : " + function_name + " throw exception while calculating G1");
			return false;
		}
	}

	/** calculate the dendity.
	* @param [in] sat_pos the sat pos in trs.
	* @param [in] sun_pos the sun pos in trs.
	* @param [in] epoch the dynamic time.
	* @param [in] kp the kp data interpolated from geomag_kp file.
	* @param [in] solar the solar data interpolated from solar_flux file.
	* @return  density.
	*/
	double t_gfatmos_dtm94::cal_density(ColumnVector &sat_pos, ColumnVector &sun_pos, double &epoch, double &kp, vector<double> &solar)// result is a little different , need more check
	{

		double sigma = 0.0;
		double Kusi = 0.0;
		double T_unlimit = 0.0;

		double Tz = 0.0;
		double gamma = 0.0;
		double Rho = 0.0;
		double Fi = 0.0;

		double density = 0.0;

		string function_name = "t_gfatmos_dtm94::calDensity";
		try
		{
			ColumnVector temp_sat, temp_sun;
			temp_sat << sat_pos;
			temp_sun << sun_pos;
			_epoch = epoch;
			_kp = kp;
			_solar = solar;
		
			temp_sat *= 1000;
			temp_sun *= 1000;
			
			xyz2blh(temp_sat, _sat_geo_crd);
			xyz2blh(temp_sun, _sun_geo_crd);

			cal_Pnm();
			cal_GL(_atemp, true, false);
			double dh = _sat_geo_crd(3) / 1000.0;

			T_unlimit = _atemp[0] * (1.0 + _GL);                          //the temperature on the top of Thermosphere
			sigma = _tgrd / (T_unlimit - _t120);                          //Relative vertical temperature gradient 
			Kusi = (_radius + _h120) / (_radius + dh)*(dh - _h120);        //geopotential height
			Tz = T_unlimit - (T_unlimit - _t120)*exp(-sigma * Kusi);       //the temperature at the height of z(km)

			bool horhe;
			for (int i = 0; i < 5; i++)
			{
				// jdhuang
				horhe = (i > 1) ? false : true;
				//horhe = true;
				//if (i > 1)
				//{
				//	horhe = false;
				//}

				gamma = _mi[i] * _uamc*_g120 / sigma / _boltman / T_unlimit * 1e3;          //Acceleration of gravity at the height of 120km
				Fi = (pow((_t120 / Tz), (1 + _alpha[i] + gamma)))*exp(-sigma * gamma*Kusi);

				//the first and the second kind atmosphere have the scale 1.0
				if (i < 4)
				{
					switch (i)
					{
					case 0:
					{
						cal_GL(_ah, false, horhe);
						_GL = _ah[0] * exp(_GL);
						break;
					}
					case 1:
					{
						cal_GL(_ahd, false, horhe);
						_GL = _ahd[0] * exp(_GL);
						break;
					}
					case 2:
					{
						cal_GL(_ao, false, horhe);
						_GL = _ao[0] * exp(_GL);
						break;
					}
					case 3:
					{
						cal_GL(_an2, false, horhe);
						_GL = _an2[0] * exp(_GL);
						break;
					}
					}
				}
				else
				{
					_GL = _o2;
				}
				Rho += _mi[i] * _GL * Fi;
			}
			density = 1.0* Rho*_uamc*1e6;
			return density;
		}
		catch (...)
		{
			write_log(_glog, 0, "ERROR : " + function_name + " throw exception while calculating Density");
			write_log(_glog, 0, " NOTE : The function will return the default value 0.0, which may bring the error but you may ignore");
			return 0.0;
		}
	}

	bool t_gfatmos_dtm94::check()
	{
		return true;
	}

	/** coordinate transforming, get the BLH according to the position in trs.
	* @param [in] crd the position in trs.
	* @param [out] ell the BLH.
	*/
	bool t_gfatmos_dtm94::xyz2blh(ColumnVector &crd, ColumnVector &ell)
	{
		string function_name = "t_gfatmos_dtm94::xyz2blh";
		try
		{
			// 这个是不是需要改写啊？？？写的感觉不是很合适
			ell.resize(3);
			const double a = Aell;
			double b = Finv;
			double e2;
			double s, rl, zps, rb, rh, n, hp, rbi, rbd, rhd;
			int i;
			// jdhuang 
			b = (b <= 6000000.0) ? a - a / b : b;
			//if (b <= 6000000.0)
			//{
			//	b = a - a / b;
			//}
			e2 = (a*a - b * b) / (a*a);
			s = sqrt(crd(1) * crd(1) + crd(2)  * crd(2));
			rl = atan2(crd(2), crd(1));
			// jdhaung
			rl = (rl < 0.0) ? 2.0*G_PI + rl : rl;
			//if (rl < 0.0)
			//	rl = 2.0*G_PI + rl;
			zps = crd(3) / s;

			rh = sqrt(crd(1)  * crd(1) + crd(2)  * crd(2) + crd(3)  * crd(3)) - a;
			rb = atan(zps / (1 - e2 * a / (a + rh)));
			/*rb = atan(zps / (1.0 - e2));
			rh = s / cos(rb) - a / sqrt(1 - e2*sin(rb)*sin(rb));*/
			for (i = 1; i <= 100; i++)
			{
				n = a / sqrt(1.0 - e2 * pow(sin(rb), 2));
				hp = rh;
				rbi = rb;
				rh = s / cos(rb) - n;
				rb = atan(zps / (1.0 - e2 * n / (n + rh)));
				rbd = abs(rbi - rb);
				rhd = abs(hp - rh);
				if (rbd <= 1.0e-11 && rhd <= 1.0e-5)
				{
					break;
				}
			}
			ell(1) = rb;
			ell(2) = rl;
			ell(3) = rh;

			return true;
		}
		catch (...)
		{
			write_log(_glog, 0, function_name, "ERROR : throw exception while change xyz to blh");
			return false;
		}

	}

}