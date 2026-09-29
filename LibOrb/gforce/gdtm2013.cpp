/**
*
* @verbatim
History
-1.0 Qian Zhang  2019-11-15 Creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University.All rights reserved.
*
* @file		    gdtm94.cpp
* @brief		    The class for dtm2013 model, serving for atmospheric drag model by providing the density.
*
* @author          Qian Zhang, Wuhan University
* @version		    1.0.0
* @date		    2019-01-06
*
*/
#include "gforce/gdtm2013.h"
#include "gutils/ginfolog.h"
#include "gforce/DTM2013.h"

using namespace std;
namespace great
{
	/** constructor. */
	t_gfatmos_dtm2013::t_gfatmos_dtm2013(t_glog* log) : t_gforce(log)
	{
		string fuction_name = "t_gfatmos_dtm2013::t_gfatmos_dtm2013";
		if (!log)
		{
			write_log(log, 0, fuction_name,
				"ERROR : the log point is empty");
		}
	}

	/** destructor. */
	t_gfatmos_dtm2013::~t_gfatmos_dtm2013()
	{

	}

	/** coordinate transforming, get the BLH according to the position in trs.
	* @param [in] crd the position in trs.
	* @param [out] ell the BLH.
	*/
	bool t_gfatmos_dtm2013::xyz2blh(ColumnVector &crd, ColumnVector &ell)
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
			e2 = (a*a - b * b) / (a*a);
			s = sqrt(crd(1) * crd(1) + crd(2)  * crd(2));
			rl = atan2(crd(2), crd(1));
			// jdhaung
			rl = (rl < 0.0) ? 2.0*G_PI + rl : rl;
			zps = crd(3) / s;

			rh = sqrt(crd(1)  * crd(1) + crd(2)  * crd(2) + crd(3)  * crd(3)) - a;
			rb = atan(zps / (1 - e2 * a / (a + rh)));

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

	double t_gfatmos_dtm2013::calDensity(float day, float f, float fbar, double * akp, float alt, float stl, float lat, float lon)   
	{
		string function_name = "t_gfatmos_dtm2013::calDensity";
		try
		{
			float tz, tinf, tp120, wmm, ro, density, ro_unc;
			float d[6] = {0.0};
			tz = tinf = tp120 = wmm = ro = density = ro_unc = 0.0;
			/*char configfile[200];
			struct dtm_date today;
			strncpy(configfile, "config.cfg", 200);
			today.type_flag = 2;
			today.day = day;
			today.month = mon;
			today.year = year;
			today.hour = hr;
			today.minute = min;
			today.second = sec;*/
			
			
			//load_config(configfile);
			//dtm_wrapper(&today, &alt, &lat, &lon, &tz, &tinf, &tp120, &ro, &ro_unc, d, &wmm);
			/*P_ReadDTM13();
			dtm2013(&day, &f, &fbar, akp, &alt, &stl, &lat, &lon, &tz, &tinf, &tp120, &ro, d, &wmm);*/
			density = ro*1000.0;
			return density;
		}
		catch (...)
		{
			write_log(_glog, 0, function_name, "ERROR : throw exception while calculating atmosphere density using dtm2013 model");
			return -1;
		}
	}
}