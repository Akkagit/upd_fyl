/**
*
* @verbatim
History
-1.0 Qian Zhang  2019-11-10 Creat the file.

@endverbatim
* Copyright (c) 2018, Wuhan University.All rights reserved.
*
* @file		    gnrlmsise00.cpp
* @brief		    The class for nrlmsise00 model, serving for atmospheric drag model by providing the density.
*
* @author          Qian Zhang, Wuhan University
* @version		    1.0.0
* @date		    2019-11-10
*
*/

#include "gforce/gnrlmsise00.h"
#include "gutils/ginfolog.h"
#include "gutils/gconst_oi.h"

using namespace std;
namespace great
{
	/** constructor. */
	t_gfatmos_nrlmsise00::t_gfatmos_nrlmsise00(t_glog* log) : t_gforce(log)
	{
		string fuction_name = "t_gfatmos_nrlmsise00::t_gfatmos_nrlmsise00";
		if (!log)
		{
			write_log(log, 0, fuction_name,
				"ERROR : the log point is empty");
		}

		for (int i = 1; i <= 24; i++)
		{
			switches.push_back(1);
		}
	}

	/** destructor. */
	t_gfatmos_nrlmsise00::~t_gfatmos_nrlmsise00()
	{

	}

	bool t_gfatmos_nrlmsise00::gtd7(int year, int doy, double sec, double alt, double g_lat, double g_long, double lst, vector<double> &d, vector<double> &t)
	{
		string function_name = "t_gfatmos_nrlmsise00::gtd7";
		try
		{		
			double xlat;
			double xmm;
			int mn3 = 5;
			double zn3[5] = { 32.5,20.0,15.0,10.0,0.0 };
			int mn2 = 4;
			double zn2[4] = { 72.5,55.0,45.0,32.5 };
			double altt;
			double zmix = 62.5;
			double tmp;
			double dm28m;
			double tz;
			double dmc;
			double dmr;
			double dz28;
			vector<double> sd;
			vector<double> st;
			int i;

			tselec();

			xlat = g_lat;
			if (sw[2] == 0)
				xlat = 45.0;
			glatf(xlat);
			xmm = pdm[2][4];
			if (alt > zn2[0])
				altt = alt;
			else
				altt = zn2[0];

			tmp = alt;
			alt = altt;
			if(!gts7(year, doy, sec, alt, g_lat, g_long, lst, sd, st))
				write_log(_glog, 0, function_name, "ERROR : when call the gts7 function in nrlmsise density calculation");
			altt = alt;
			alt = tmp;
			if (sw[0])
				dm28m = dm28*1.0e6;
			else
				dm28m = dm28;
			t.push_back(st[0]);
			t.push_back(st[1]);
			if (alt >= zn2[0])
			{
				for (i = 0; i < 9; i++)
					d.push_back(sd[i]);
				return true;
			}
			meso_tgn2[0] = meso_tgn1[1];
			meso_tn2[0] = meso_tn1[4];
			meso_tn2[1] = pma[0][0] * pavgm[0] / (1.0 - sw[20] * glob7s(pma[0], year, doy, sec, alt, g_lat, g_long, lst));
			meso_tn2[2] = pma[1][0] * pavgm[1] / (1.0 - sw[20] * glob7s(pma[1], year, doy, sec, alt, g_lat, g_long, lst));
			meso_tn2[3] = pma[2][0] * pavgm[2] / (1.0 - sw[20] * sw[22] * glob7s(pma[2], year, doy, sec, alt, g_lat, g_long, lst));
			meso_tgn2[1] = pavgm[8] * pma[9][0] * (1.0 + sw[20] * sw[22] * glob7s(pma[9], year, doy, sec, alt, g_lat, g_long, lst))*meso_tn2[3] * meso_tn2[3] / (pow((pma[2][0] * pavgm[2]), 2.0));
			meso_tn3[0] = meso_tn2[3];
			if (alt <= zn3[0])
			{
				meso_tgn3[0] = meso_tgn2[1];
				meso_tn3[1] = pma[3][0] * pavgm[3] / (1.0 - sw[22] * glob7s(pma[3], year, doy, sec, alt, g_lat, g_long, lst));
				meso_tn3[2] = pma[4][0] * pavgm[4] / (1.0 - sw[22] * glob7s(pma[4], year, doy, sec, alt, g_lat, g_long, lst));
				meso_tn3[3] = pma[5][0] * pavgm[5] / (1.0 - sw[22] * glob7s(pma[5], year, doy, sec, alt, g_lat, g_long, lst));
				meso_tn3[4] = pma[6][0] * pavgm[6] / (1.0 - sw[22] * glob7s(pma[6], year, doy, sec, alt, g_lat, g_long, lst));
				meso_tgn3[1] = pma[7][0] * pavgm[7] * (1.0 + sw[22] * glob7s(pma[7], year, doy, sec, alt, g_lat, g_long, lst)) *meso_tn3[4] * meso_tn3[4] / (pow((pma[6][0] * pavgm[6]), 2.0));
			}

			dmc = 0;
			if(alt>zmix)
				dmc = 1.0 - (zn2[0] - alt) / (zn2[0] - zmix);
			dz28 = sd[2];

			/**** N2 density ****/
			dmr = sd[2] / dm28m - 1.0;
			d[2] = densm(alt, dm28m, xmm, &tz, mn3, zn3, meso_tn3, meso_tgn3, mn2, zn2);
			d[2] = d[2] * (1.0 + dmr*dmc);

			/**** HE density ****/
			dmr = sd[0] / (dz28 * pdm[0][1]) - 1.0;
			d[0] = d[2] * pdm[0][1] * (1.0 + dmr*dmc);

			/**** O density ****/
			d[1] = 0;
			d[8] = 0;

			/**** O2 density ****/
			dmr = sd[3] / (dz28 * pdm[3][1]) - 1.0;
			d[3] = d[2] * pdm[3][1] * (1.0 + dmr*dmc);

			/**** AR density ***/
			dmr = sd[4] / (dz28 * pdm[4][1]) - 1.0;
			d[4] = d[2] * pdm[4][1] * (1.0 + dmr*dmc);

			/**** Hydrogen density ****/
			d[6] = 0;

			/**** Atomic nitrogen density ****/
			d[7] = 0;

			/**** Total mass density */
			d[5] = 1.66E-24 * (4.0 * d[0] + 16.0 * d[1] + 28.0 * d[2] + 32.0 * d[3] + 40.0 * d[4] + d[6] + 14.0 * d[7]);

			if (sw[0])
				d[5] = d[5] / 1000;

			/**** temperature at altitude ****/
			dd = densm(alt, 1.0, 0, &tz, mn3, zn3, meso_tn3, meso_tgn3, mn2, zn2);
			t[1] = tz;

			return true;
		}
		catch (...)
		{
			write_log(_glog, 0, function_name, "ERROR : in the gtd7 function in nrlmsise density calculation");
			return false;
		}
	}

	bool t_gfatmos_nrlmsise00::gtd7d(int year, int doy, double sec, double alt, double g_lat, double g_long, double lst, double f107a, double f107, vector<double> &ap, vector<double> &d, vector<double> &t)
	{
		string function_name = "t_gfatmos_nrlmsise00::gtd7d";
		try
		{
			_f107a = f107a;
			_f107 = f107;
			_ap = ap;
			_dap = ap[0];
			if (!gtd7(year, doy, sec, alt, g_lat, g_long, lst, d, t))
			{
				write_log(_glog, 0, function_name + "ERROR : in gtd7d while call gtd7");
			}
			d[5] = 1.66E-24 * (4.0 * d[0] + 16.0 * d[1] + 28.0 * d[2] + 32.0 * d[3] + 40.0 * d[4] + d[6] + 14.0 * d[7] + 16.0 *d[8]);
			if (sw[0])
				d[5] = d[5] / 1000;
			return true;
		}
		catch (...)
		{
			write_log(_glog, 0, function_name + "ERROR : throw exception while calculating density in gtd7d");
			return false;
		}
	}

	bool t_gfatmos_nrlmsise00::ghp7(int year, int doy, double sec, double alt, double g_lat, double g_long, double lst, vector<double> &d, vector<double> &t, double press)
	{
		string function_name = "t_gfatmos_nrlmsise00::gtd7";
		try
		{
			double bm = 1.3806E-19;
			double rgas = 831.4;
			double test = 0.00043;
			double ltest = 12;
			double pl, p;
			double zi;
			double z;
			double cl, cl2;
			double ca, cd;
			double xn, xm, diff;
			double g, sh;
			int l;
			pl = log10(press);

			if (pl >= -5.0) 
			{
				if (pl>2.5)
					zi = 18.06 * (3.00 - pl);
				else if ((pl>0.075) && (pl <= 2.5))
					zi = 14.98 * (3.08 - pl);
				else if ((pl>-1) && (pl <= 0.075))
					zi = 17.80 * (2.72 - pl);
				else if ((pl>-2) && (pl <= -1))
					zi = 14.28 * (3.64 - pl);
				else if ((pl>-4) && (pl <= -2))
					zi = 12.72 * (4.32 - pl);
				else
					zi = 25.3 * (0.11 - pl);
				cl = g_lat / 90.0;
				cl2 = cl*cl;
				if (doy<182)
					cd = (1.0 - (double)doy) / 91.25;
				else
					cd = ((double)doy) / 91.25 - 3.0;
				ca = 0;
				if ((pl > -1.11) && (pl <= -0.23))
					ca = 1.0;
				if (pl > -0.23)
					ca = (2.79 - pl) / (2.79 + 0.23);
				if ((pl <= -1.11) && (pl>-3))
					ca = (-2.93 - pl) / (-2.93 + 1.11);
				z = zi - 4.87 * cl * cd * ca - 1.64 * cl2 * ca + 0.31 * ca * cl;
			}
			else
				z = 22.0 * pow((pl + 4.0), 2.0) + 110.0;

			/* iteration  loop */
			l = 0;
			do 
			{
				l++;
				alt = z;
				if(!gtd7(year,doy,sec,alt,g_lat,g_long,lst, d,t))
					write_log(_glog, 0, function_name + "ERROR : when call the gtd7 function in t_gfatmos_nrlmsise00::ghp7");
				z = alt;
				xn = d[0] + d[1] + d[2] + d[3] + d[4] + d[6] + d[7];
				p = bm * xn * t[1];
				if (sw[0])
					p = p*1.0E-6;
				diff = pl - log10(p);
				if (sqrt(diff*diff)<test)
					return false;
				if (l == ltest) {
					write_log(_glog, 0, function_name + "ERROR: ghp7 not converging for press " + to_string(press) + ", diff " + to_string(diff));
					return false;
				}
				xm = d[5] / xn / 1.66E-24;
				if (sw[0])
					xm = xm * 1.0E3;
				g = gsurf / (pow((1.0 + z / re), 2.0));
				sh = rgas * t[1] / (xm * g);

				/* new altitude estimate using scale height */
				if (l <  6)
					z = z - sh * diff * 2.302;
				else
					z = z - sh * diff;
			} while (1 == 1);
			return true;
		}
		catch (...)
		{
			write_log(_glog, 0, function_name + "ERROR : throw exception while calculating density in ghp7");
			return false;
		}
	}

	/*      Calculate Temperature and Density Profiles for lower atmos.  */
	double t_gfatmos_nrlmsise00::densm(double alt, double d0, double xm, double *tz, int mn3, double *zn3, double *tn3, double *tgn3, int mn2, double *zn2)
	{
		string function_name = "t_gfatmos_nrlmsise00::densm";
		try
		{
			double xs[10], ys[10], y2out[10];
			double rgas = 831.4;
			double z, z1, z2, t1, t2, zg, zgdif;
			double yd1, yd2;
			double x, y, yi;
			double expl, gamm, glb;
			double densm_tmp;
			int mn;
			int k;
			densm_tmp = d0;
			if (alt>zn2[0]) 
			{
				if (xm == 0.0)
					return *tz;
				else
					return d0;
			}

			/* STRATOSPHERE/MESOSPHERE TEMPERATURE */
			if (alt>zn2[mn2 - 1])
				z = alt;
			else
				z = zn2[mn2 - 1];
			mn = mn2;
			z1 = zn2[0];
			z2 = zn2[mn - 1];
			t1 = meso_tn2[0];
			t2 = meso_tn2[mn - 1];
			zg = zeta(z, z1);
			zgdif = zeta(z2, z1);

			/* set up spline nodes */
			for (k = 0; k<mn; k++) {
				xs[k] = zeta(zn2[k], z1) / zgdif;
				ys[k] = 1.0 / meso_tn2[k];
			}
			yd1 = -meso_tgn2[0] / (t1*t1) * zgdif;
			yd2 = -meso_tgn2[1] / (t2*t2) * zgdif * (pow(((re + z2) / (re + z1)), 2.0));

			/* calculate spline coefficients */
			spline(xs, ys, mn, yd1, yd2, y2out);
			x = zg / zgdif;
			splint(xs, ys, y2out, mn, x, &y);

			/* temperature at altitude */
			*tz = 1.0 / y;
			if (xm != 0.0) 
			{
				/* calaculate stratosphere / mesospehere density */
				glb = gsurf / (pow((1.0 + z1 / re), 2.0));
				gamm = xm * glb * zgdif / rgas;

				/* Integrate temperature profile */
				splini(xs, ys, y2out, mn, x, &yi);
				expl = gamm*yi;
				if (expl>50.0)
					expl = 50.0;

				/* Density at altitude */
				densm_tmp = densm_tmp * (t1 / *tz) * exp(-expl);
			}

			if (alt>zn3[0]) 
			{
				if (xm == 0.0)
					return *tz;
				else
					return densm_tmp;
			}

			/* troposhere / stratosphere temperature */
			z = alt;
			mn = mn3;
			z1 = zn3[0];
			z2 = zn3[mn - 1];
			t1 = tn3[0];
			t2 = tn3[mn - 1];
			zg = zeta(z, z1);
			zgdif = zeta(z2, z1);

			/* set up spline nodes */
			for (k = 0; k<mn; k++) {
				xs[k] = zeta(zn3[k], z1) / zgdif;
				ys[k] = 1.0 / tn3[k];
			}
			yd1 = -tgn3[0] / (t1*t1) * zgdif;
			yd2 = -tgn3[1] / (t2*t2) * zgdif * (pow(((re + z2) / (re + z1)), 2.0));

			/* calculate spline coefficients */
			spline(xs, ys, mn, yd1, yd2, y2out);
			x = zg / zgdif;
			splint(xs, ys, y2out, mn, x, &y);

			/* temperature at altitude */
			*tz = 1.0 / y;
			if (xm != 0.0) 
			{
				/* calaculate tropospheric / stratosphere density */
				glb = gsurf / (pow((1.0 + z1 / re), 2.0));
				gamm = xm * glb * zgdif / rgas;

				/* Integrate temperature profile */
				splini(xs, ys, y2out, mn, x, &yi);
				expl = gamm*yi;
				if (expl>50.0)
					expl = 50.0;

				/* Density at altitude */
				densm_tmp = densm_tmp * (t1 / *tz) * exp(-expl);
			}
			if (xm == 0.0)
				return *tz;
			else
				return densm_tmp;
		}
		catch (...)
		{
			write_log(_glog, 0, function_name + "ERROR : throw exception while calculating density in densm");
			return 0.0;
		}
	}

	void t_gfatmos_nrlmsise00::tselec()
	{
		int i;
		for (i = 0; i<24; i++) 
		{
			if (i != 9) {
				if (switches[i] == 1)
					sw.push_back(1);
				else
					sw.push_back(0);
				if (switches[i]>0)
					swc.push_back(1);
				else
					swc.push_back(0);
			}
			else {
				sw.push_back(switches[i]);
				swc.push_back(switches[i]);
			}
		}
	}

	void t_gfatmos_nrlmsise00::glatf(double xlat)
	{
		double c2;
		double dgtr = 1.74533e-2;
		c2 = cos(2.0*dgtr*xlat);
		gsurf = 980.616*(1.0 - 0.0026373*c2);
		re = 2.0*gsurf*1.0e-5 / (3.085462e-6 + 2.27e-9*c2);
	}

	//int t_gfatmos_nrlmsise00::vtst7(int iyd, double sec, double glat, double glong, double stl, double f107a, double f107, vector<double> &ap, int ic)
	//{
	//	int iydl[2] = { -999, -999 };
	//	double secl[2] = { -999.0, -999.0 };
	//	double glatl[2] = { -999.0, -999.0 };
	//	double gll[2] = { -999.0, -999.0 };
	//	double stll[2] = { -999.0, -999.0 };
	//	double fal[2] = { -999.0, -999.0 };
	//	double fl[2] = { -999.0, -999.0 };
	//	double apl[14] = { -999.0};
	//	double swl[50] = { -999.0};  //如何整体赋值？
	//	double swcl[50] = { -999.0};
	//	int _vtst7 = 0;
	//	int    flag = 0;
	//	if (iyd == iydl[ic - 1] && sec == secl[ic - 1] && glat == glatl[ic - 1] && glong == gll[ic - 1] && stl == stll[ic - 1] && f107a == fal[ic - 1] && f107 == fl[ic - 1])
	//	{
	//		for (int i = 1; i <= 7; i++)
	//		{
	//			if (ap[i - 1] != apl[(ic - 1) * 7 + i - 1])
	//			{
	//				flag = 1;
	//				break;
	//			}
	//		}
	//		if (flag == 0)
	//		{
	//			for (int i = 1; i <= 25; i++)
	//			{
	//				if (sw[i - 1] != swl[(ic - 1) * 25 + i - 1] || swc[i - 1] != swcl[(ic - 1) * 25 + i - 1])
	//				{
	//					flag = 1;
	//					break;
	//				}
	//			}
	//		}
	//	}
	//	else
	//	{
	//		flag = 1;
	//	}
	//		if (flag == 1)
	//		{
	//			_vtst7 = 1;
	//			iydl[ic - 1] = iyd;
	//			secl[ic - 1] = sec;
	//			glatl[ic - 1] = glat;
	//			gll[ic - 1] = glong;
	//			stll[ic - 1] = stl;
	//			fal[ic - 1] = f107a;
	//			fl[ic - 1] = f107;
	//			for (int i = 1; i <= 7; i++)
	//			{
	//				apl[(ic - 1) * 7 + i - 1] = ap[i - 1];
	//			}
	//			for (int i = 1; i <= 25; i++)
	//			{
	//				swl[i - 1] = sw[(ic - 1) * 25 + i - 1];
	//				swcl[i - 1] = swc[(ic - 1) * 25 + i - 1];
	//			}
	//		}
	//		return _vtst7;
	//}

	bool t_gfatmos_nrlmsise00::gts7(int year, int doy, double sec, double alt, double g_lat, double g_long, double lst, vector<double> &d, vector<double> &t)
	{
		string function_name = "t_gfatmos_nrlmsise00::gts7";
		try
		{
			double za;
			int i, j;
			double ddum, z;
			double zn1[5] = { 120.0, 110.0, 100.0, 90.0, 72.5 };
			double tinf;
			int mn1 = 5;
			double g0;
			double tlb;
			double s;
			double db01, db04, db14, db16, db28, db32, db40;
			double zh28, zh04, zh16, zh32, zh40, zh01, zh14;
			double zhm28, zhm04, zhm16, zhm32, zhm40, zhm01, zhm14;
			double xmd;
			double b28, b04, b16, b32, b40, b01, b14;
			double tz;
			double g28, g4, g16, g32, g40, g1, g14;
			double zhf, xmm;
			double zc04, zc16, zc32, zc40, zc01, zc14;
			double hc04, hc16, hc32, hc40, hc01, hc14;
			double hcc16, hcc32, hcc01, hcc14;
			double zcc16, zcc32, zcc01, zcc14;
			double rc16, rc32, rc01, rc14;
			double rl;
			double g16h, db16h, tho, zsht, zmho, zsho;
			double dgtr = 1.74533e-2;
			double dr = 1.72142e-2;
			double alpha[9] = { -0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0 };
			double altl[8] = { 200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0 };
			double dd;
			double hc216, hcc232;

			za = pdl[1][15];
			zn1[0] = za;
			for (j = 0; j < 9; j++)
				d.push_back(0);

			/* TINF VARIATIONS NOT IMPORTANT BELOW ZA OR ZN1(1) */
			if (alt > zn1[0])
				tinf = ptm[0] * pt[0] * \
				(1.0 + sw[16] * globe7(pt, year, doy, sec, alt, g_lat, g_long, lst));
			else
				tinf = ptm[0] * pt[0];
			t.push_back(tinf);

			/*  GRADIENT VARIATIONS NOT IMPORTANT BELOW ZN1(5) */
			if (alt > zn1[4])
				g0 = ptm[3] * ps[0] * \
				(1.0 + sw[19] * globe7(ps, year, doy, sec, alt, g_lat, g_long, lst));
			else
				g0 = ptm[3] * ps[0];
			tlb = ptm[1] * (1.0 + sw[17] * globe7(pd[3], year, doy, sec, alt, g_lat, g_long, lst))*pd[3][0];
			s = g0 / (tinf - tlb);

			/*      Lower thermosphere temp variations not significant for
			*       density above 300 km */
			if (alt < 300.0)
			{
				meso_tn1[1] = ptm[6] * ptl[0][0] / (1.0 - sw[18] * glob7s(ptl[0], year, doy, sec, alt, g_lat, g_long, lst));
				meso_tn1[2] = ptm[2] * ptl[1][0] / (1.0 - sw[18] * glob7s(ptl[1], year, doy, sec, alt, g_lat, g_long, lst));
				meso_tn1[3] = ptm[7] * ptl[2][0] / (1.0 - sw[18] * glob7s(ptl[2], year, doy, sec, alt, g_lat, g_long, lst));
				meso_tn1[4] = ptm[4] * ptl[3][0] / (1.0 - sw[18] * sw[20] * glob7s(ptl[3], year, doy, sec, alt, g_lat, g_long, lst));
				meso_tgn1[1] = ptm[8] * pma[8][0] * (1.0 + sw[18] * sw[20] * glob7s(pma[8], year, doy, sec, alt, g_lat, g_long, lst))*meso_tn1[4] * meso_tn1[4] / (pow((ptm[4] * ptl[3][0]), 2.0));
			}
			else {
				meso_tn1[1] = ptm[6] * ptl[0][0];
				meso_tn1[2] = ptm[2] * ptl[1][0];
				meso_tn1[3] = ptm[7] * ptl[2][0];
				meso_tn1[4] = ptm[4] * ptl[3][0];
				meso_tgn1[1] = ptm[8] * pma[8][0] * meso_tn1[4] * meso_tn1[4] / (pow((ptm[4] * ptl[3][0]), 2.0));
			}

			/* N2 variation factor at Zlb */
			g28 = sw[21] * globe7(pd[2], year, doy, sec, alt, g_lat, g_long, lst);

			/* VARIATION OF TURBOPAUSE HEIGHT */
			zhf = pdl[1][24] * (1.0 + sw[5] * pdl[0][24] * sin(dgtr*g_lat)*cos(dr*(doy - pt[13])));
			t[0] = tinf;
			t.push_back(0.0);
			xmm = pdm[2][4];
			z = alt;

			/**** N2 DENSITY ****/
			/* Diffusive density at Zlb */
			db28 = pdm[2][0] * exp(g28)*pd[2][0];
			/* Diffusive density at Alt */
			d[2] = densu(z, db28, tinf, tlb, 28.0, alpha[2], &t[1], ptm[5], s, mn1, zn1);
			dd = d[2];
			/* Turbopause */
			zh28 = pdm[2][2] * zhf;
			zhm28 = pdm[2][3] * pdl[1][5];
			xmd = 28.0 - xmm;
			/* Mixed density at Zlb */
			b28 = densu(zh28, db28, tinf, tlb, xmd, (alpha[2] - 1.0), &tz, ptm[5], s, mn1, zn1);
			if ((sw[15]) && (z <= altl[2]))
			{
				/*  Mixed density at Alt */
				dm28 = densu(z, b28, tinf, tlb, xmm, alpha[2], &tz, ptm[5], s, mn1, zn1);
				/*  Net density at Alt */
				d[2] = dnet(d[2], dm28, zhm28, xmm, 28.0);
			}

			/**** HE DENSITY ****/
			/*   Density variation factor at Zlb */
			g4 = sw[21] * globe7(pd[0], year, doy, sec, alt, g_lat, g_long, lst);
			/*  Diffusive density at Zlb */
			db04 = pdm[0][0] * exp(g4)*pd[0][0];
			/*  Diffusive density at Alt */
			d[0] = densu(z, db04, tinf, tlb, 4., alpha[0], &t[1], ptm[5], s, mn1, zn1);
			dd = d[0];
			if ((sw[15]) && (z < altl[0]))
			{
				/*  Turbopause */
				zh04 = pdm[0][2];
				/*  Mixed density at Zlb */
				b04 = densu(zh04, db04, tinf, tlb, 4. - xmm, alpha[0] - 1., &t[1], ptm[5], s, mn1, zn1);
				/*  Mixed density at Alt */
				dm04 = densu(z, b04, tinf, tlb, xmm, 0., &t[1], ptm[5], s, mn1, zn1);
				zhm04 = zhm28;
				/*  Net density at Alt */
				d[0] = dnet(d[0], dm04, zhm04, xmm, 4.);
				/*  Correction to specified mixing ratio at ground */
				rl = log(b28*pdm[0][1] / b04);
				zc04 = pdm[0][4] * pdl[1][0];
				hc04 = pdm[0][5] * pdl[1][1];
				/*  Net density corrected at Alt */
				d[0] = d[0] * ccor(z, rl, hc04, zc04);
			}

			/**** O DENSITY ****/
			/*  Density variation factor at Zlb */
			g16 = sw[21] * globe7(pd[1], year, doy, sec, alt, g_lat, g_long, lst);
			/*  Diffusive density at Zlb */
			db16 = pdm[1][0] * exp(g16)*pd[1][0];
			/*   Diffusive density at Alt */
			d[1] = densu(z, db16, tinf, tlb, 16., alpha[1], &t[1], ptm[5], s, mn1, zn1);
			dd = d[1];
			if ((sw[15]) && (z <= altl[1]))
			{
				/*   Turbopause */
				zh16 = pdm[1][2];
				/*  Mixed density at Zlb */
				b16 = densu(zh16, db16, tinf, tlb, 16.0 - xmm, (alpha[1] - 1.0), &t[1], ptm[5], s, mn1, zn1);
				/*  Mixed density at Alt */
				dm16 = densu(z, b16, tinf, tlb, xmm, 0., &t[1], ptm[5], s, mn1, zn1);
				zhm16 = zhm28;
				/*  Net density at Alt */
				d[1] = dnet(d[1], dm16, zhm16, xmm, 16.);
				rl = pdm[1][1] * pdl[1][16] * (1.0 + sw[1] * pdl[0][23] * (_f107a - 150.0));
				hc16 = pdm[1][5] * pdl[1][3];
				zc16 = pdm[1][4] * pdl[1][2];
				hc216 = pdm[1][5] * pdl[1][4];
				d[1] = d[1] * ccor2(z, rl, hc16, zc16, hc216);
				/*   Chemistry correction */
				hcc16 = pdm[1][7] * pdl[1][13];
				zcc16 = pdm[1][6] * pdl[1][12];
				rc16 = pdm[1][3] * pdl[1][14];
				/*  Net density corrected at Alt */
				d[1] = d[1] * ccor(z, rc16, hcc16, zcc16);
			}

			/**** O2 DENSITY ****/
			/*   Density variation factor at Zlb */
			g32 = sw[21] * globe7(pd[4], year, doy, sec, alt, g_lat, g_long, lst);
			/*  Diffusive density at Zlb */
			db32 = pdm[3][0] * exp(g32)*pd[4][0];
			/*   Diffusive density at Alt */
			d[3] = densu(z, db32, tinf, tlb, 32., alpha[3], &t[1], ptm[5], s, mn1, zn1);
			dd = d[3];
			if (sw[15])
			{
				if (z <= altl[3]) {
					/*   Turbopause */
					zh32 = pdm[3][2];
					/*  Mixed density at Zlb */
					b32 = densu(zh32, db32, tinf, tlb, 32. - xmm, alpha[3] - 1., &t[1], ptm[5], s, mn1, zn1);
					/*  Mixed density at Alt */
					dm32 = densu(z, b32, tinf, tlb, xmm, 0., &t[1], ptm[5], s, mn1, zn1);
					zhm32 = zhm28;
					/*  Net density at Alt */
					d[3] = dnet(d[3], dm32, zhm32, xmm, 32.);
					/*   Correction to specified mixing ratio at ground */
					rl = log(b28*pdm[3][1] / b32);
					hc32 = pdm[3][5] * pdl[1][7];
					zc32 = pdm[3][4] * pdl[1][6];
					d[3] = d[3] * ccor(z, rl, hc32, zc32);
				}
				/*  Correction for general departure from diffusive equilibrium above Zlb */
				hcc32 = pdm[3][7] * pdl[1][22];
				hcc232 = pdm[3][7] * pdl[0][22];
				zcc32 = pdm[3][6] * pdl[1][21];
				rc32 = pdm[3][3] * pdl[1][23] * (1. + sw[1] * pdl[0][23] * (_f107a - 150.));
				/*  Net density corrected at Alt */
				d[3] = d[3] * ccor2(z, rc32, hcc32, zcc32, hcc232);
			}

			/**** AR DENSITY ****/
			/*   Density variation factor at Zlb */
			g40 = sw[21] * globe7(pd[5], year, doy, sec, alt, g_lat, g_long, lst);
			/*  Diffusive density at Zlb */
			db40 = pdm[4][0] * exp(g40)*pd[5][0];
			/*   Diffusive density at Alt */
			d[4] = densu(z, db40, tinf, tlb, 40., alpha[4], &t[1], ptm[5], s, mn1, zn1);
			dd = d[4];
			if ((sw[15]) && (z <= altl[4]))
			{
				/*   Turbopause */
				zh40 = pdm[4][2];
				/*  Mixed density at Zlb */
				b40 = densu(zh40, db40, tinf, tlb, 40. - xmm, alpha[4] - 1., &t[1], ptm[5], s, mn1, zn1);
				/*  Mixed density at Alt */
				dm40 = densu(z, b40, tinf, tlb, xmm, 0., &t[1], ptm[5], s, mn1, zn1);
				zhm40 = zhm28;
				/*  Net density at Alt */
				d[4] = dnet(d[4], dm40, zhm40, xmm, 40.);
				/*   Correction to specified mixing ratio at ground */
				rl = log(b28*pdm[4][1] / b40);
				hc40 = pdm[4][5] * pdl[1][9];
				zc40 = pdm[4][4] * pdl[1][8];
				/*  Net density corrected at Alt */
				d[4] = d[4] * ccor(z, rl, hc40, zc40);
			}

			/**** HYDROGEN DENSITY ****/
			/*   Density variation factor at Zlb */
			g1 = sw[21] * globe7(pd[6], year, doy, sec, alt, g_lat, g_long, lst);
			/*  Diffusive density at Zlb */
			db01 = pdm[5][0] * exp(g1)*pd[6][0];
			/*   Diffusive density at Alt */
			d[6] = densu(z, db01, tinf, tlb, 1., alpha[6], &t[1], ptm[5], s, mn1, zn1);
			dd = d[6];
			if ((sw[15]) && (z <= altl[6]))
			{
				/*   Turbopause */
				zh01 = pdm[5][2];
				/*  Mixed density at Zlb */
				b01 = densu(zh01, db01, tinf, tlb, 1. - xmm, alpha[6] - 1., &t[1], ptm[5], s, mn1, zn1);
				/*  Mixed density at Alt */
				dm01 = densu(z, b01, tinf, tlb, xmm, 0., &t[1], ptm[5], s, mn1, zn1);
				zhm01 = zhm28;
				/*  Net density at Alt */
				d[6] = dnet(d[6], dm01, zhm01, xmm, 1.);
				/*   Correction to specified mixing ratio at ground */
				rl = log(b28*pdm[5][1] * sqrt(pdl[1][17] * pdl[1][17]) / b01);
				hc01 = pdm[5][5] * pdl[1][11];
				zc01 = pdm[5][4] * pdl[1][10];
				d[6] = d[6] * ccor(z, rl, hc01, zc01);
				/*   Chemistry correction */
				hcc01 = pdm[5][7] * pdl[1][19];
				zcc01 = pdm[5][6] * pdl[1][18];
				rc01 = pdm[5][3] * pdl[1][20];
				/*  Net density corrected at Alt */
				d[6] = d[6] * ccor(z, rc01, hcc01, zcc01);
			}

			/**** ATOMIC NITROGEN DENSITY ****/
			/*   Density variation factor at Zlb */
			g14 = sw[21] * globe7(pd[7], year, doy, sec, alt, g_lat, g_long, lst);
			/*  Diffusive density at Zlb */
			db14 = pdm[6][0] * exp(g14)*pd[7][0];
			/*   Diffusive density at Alt */
			d[7] = densu(z, db14, tinf, tlb, 14., alpha[7], &t[1], ptm[5], s, mn1, zn1);
			dd = d[7];
			if ((sw[15]) && (z <= altl[7]))
			{
				/*   Turbopause */
				zh14 = pdm[6][2];
				/*  Mixed density at Zlb */
				b14 = densu(zh14, db14, tinf, tlb, 14. - xmm, alpha[7] - 1., &t[1], ptm[5], s, mn1, zn1);
				/*  Mixed density at Alt */
				dm14 = densu(z, b14, tinf, tlb, xmm, 0., &t[1], ptm[5], s, mn1, zn1);
				zhm14 = zhm28;
				/*  Net density at Alt */
				d[7] = dnet(d[7], dm14, zhm14, xmm, 14.);
				/*   Correction to specified mixing ratio at ground */
				rl = log(b28*pdm[6][1] * sqrt(pdl[0][2] * pdl[0][2]) / b14);
				hc14 = pdm[6][5] * pdl[0][1];
				zc14 = pdm[6][4] * pdl[0][0];
				d[7] = d[7] * ccor(z, rl, hc14, zc14);
				/*   Chemistry correction */
				hcc14 = pdm[6][7] * pdl[0][4];
				zcc14 = pdm[6][6] * pdl[0][3];
				rc14 = pdm[6][3] * pdl[0][5];
				/*  Net density corrected at Alt */
				d[7] = d[7] * ccor(z, rc14, hcc14, zcc14);
			}

			/**** Anomalous OXYGEN DENSITY ****/
			g16h = sw[21] * globe7(pd[8], year, doy, sec, alt, g_lat, g_long, lst);
			db16h = pdm[7][0] * exp(g16h)*pd[8][0];
			tho = pdm[7][9] * pdl[0][6];
			dd = densu(z, db16h, tho, tho, 16., alpha[8], &t[1], ptm[5], s, mn1, zn1);
			zsht = pdm[7][5];
			zmho = pdm[7][4];
			double rgas = 831.4;
			zsho = gsurf / (pow((1.0 + zmho / re), 2.0));
			zsho = rgas * tho / (zsho * 16.0);
//			zsho = scalh(zmho, 16.0, tho);
			d[8] = dd*exp(-zsht / zsho*(exp(-(z - zmho) / zsht) - 1.));

			/* total mass density */
			d[5] = 1.66E-24*(4.0*d[0] + 16.0*d[1] + 28.0*d[2] + 32.0*d[3] + 40.0*d[4] + d[6] + 14.0*d[7]);

			/* temperature */
			z = sqrt(alt*alt);
			ddum = densu(z, 1.0, tinf, tlb, 0.0, 0.0, &t[1], ptm[5], s, mn1, zn1);
			(void)ddum; /* silence gcc */
			if (sw[0]) {
				for (i = 0; i < 9; i++)
					d[i] = d[i] * 1.0E6;
				d[5] = d[5] / 1000;
			}
			return true;
		}
		catch (...)
		{
			write_log(_glog, 0, function_name, "ERROR : throw exception in gts7 function while calculate the nrlmsise00 density model");
			return false;
		}
	}

	double t_gfatmos_nrlmsise00::glob7s(double * par, int year, int doy, double sec, double alt, double g_lat, double g_long, double lst)   //不会调用到而且fortran版本看不懂
	{
		string function_name = "t_gfatmos_nrlmsise00::glob7s";
		try
		{
			double pset = 2.0;
			double t[14];
			double tt;
			double cd32, cd18, cd14, cd39;
			int i, j;
			double dr = 1.72142E-2;
			double dgtr = 1.74533E-2;

			/* confirm parameter set */
			if (par[99] == 0)
				par[99] = pset;
			if (par[99] != pset) 
			{
				write_log(_glog, 0, function_name, "ERROR : Wrong parameter set for glob7s");
				throw;
			}
			for (j = 0; j<14; j++)
				t[j] = 0.0;
			cd32 = cos(dr*(doy - par[31]));
			cd18 = cos(2.0*dr*(doy - par[17]));
			cd14 = cos(dr*(doy - par[13]));
			cd39 = cos(2.0*dr*(doy - par[38]));

			/* F10.7 */
			t[0] = par[21] * dfa;

			/* time independent */
			t[1] = par[1] * plg[0][2] + par[2] * plg[0][4] + par[22] * plg[0][6] + par[26] * plg[0][1] + par[14] * plg[0][3] + par[59] * plg[0][5];

			/* SYMMETRICAL ANNUAL */
			t[2] = (par[18] + par[47] * plg[0][2] + par[29] * plg[0][4])*cd32;

			/* SYMMETRICAL SEMIANNUAL */
			t[3] = (par[15] + par[16] * plg[0][2] + par[30] * plg[0][4])*cd18;

			/* ASYMMETRICAL ANNUAL */
			t[4] = (par[9] * plg[0][1] + par[10] * plg[0][3] + par[20] * plg[0][5])*cd14;

			/* ASYMMETRICAL SEMIANNUAL */
			t[5] = (par[37] * plg[0][1])*cd39;

			/* DIURNAL */
			if (sw[7]) 
			{
				double t71, t72;
				t71 = par[11] * plg[1][2] * cd14*swc[5];
				t72 = par[12] * plg[1][2] * cd14*swc[5];
				t[6] = ((par[3] * plg[1][1] + par[4] * plg[1][3] + t71) * ctloc + (par[6] * plg[1][1] + par[7] * plg[1][3] + t72) * stloc);
			}

			/* SEMIDIURNAL */
			if (sw[8]) {
				double t81, t82;
				t81 = (par[23] * plg[2][3] + par[35] * plg[2][5])*cd14*swc[5];
				t82 = (par[33] * plg[2][3] + par[36] * plg[2][5])*cd14*swc[5];
				t[7] = ((par[5] * plg[2][2] + par[41] * plg[2][4] + t81) * c2tloc + (par[8] * plg[2][2] + par[42] * plg[2][4] + t82) * s2tloc);
			}

			/* TERDIURNAL */
			if (sw[14]) {
				t[13] = par[39] * plg[3][3] * s3tloc + par[40] * plg[3][3] * c3tloc;
			}

			/* MAGNETIC ACTIVITY */
			if (sw[9]) {
				if (sw[9] == 1)
					t[8] = apdf * (par[32] + par[45] * plg[0][2] * swc[2]);
				if (sw[9] == -1)
					t[8] = (par[50] * apt[0] + par[96] * plg[0][2] * apt[0] * swc[2]);
			}

			/* LONGITUDINAL */
			if (!((sw[10] == 0) || (sw[11] == 0) || (g_long <= -1000.0))) 
			{
				t[10] = (1.0 + plg[0][1] * (par[80] * swc[5] * cos(dr*(doy - par[81]))\
					+ par[85] * swc[6] * cos(2.0*dr*(doy - par[86])))\
					+ par[83] * swc[3] * cos(dr*(doy - par[84]))\
					+ par[87] * swc[4] * cos(2.0*dr*(doy - par[88])))\
					*((par[64] * plg[1][2] + par[65] * plg[1][4] + par[66] * plg[1][6]\
						+ par[74] * plg[1][1] + par[75] * plg[1][3] + par[76] * plg[1][5]\
						)*cos(dgtr*g_long)\
						+ (par[90] * plg[1][2] + par[91] * plg[1][4] + par[92] * plg[1][6]\
							+ par[77] * plg[1][1] + par[78] * plg[1][3] + par[79] * plg[1][5]\
							)*sin(dgtr*g_long));
			}
			tt = 0;
			for (i = 0; i<14; i++)
				tt += fabs(sw[(int64_t)i + 1])*t[i];
			return tt;
		}
		catch (...)
		{
			write_log(_glog, 0, function_name, "ERROR : throw exception in glob7s function while calculate the nrlmsise00 density model");
			return 0.0;
		}
	}

	/*      Calculate Temperature and Density Profiles for MSIS models
	*      New lower thermo polynomial
	*/
	double t_gfatmos_nrlmsise00::densu(double alt, double dlb, double tinf, double _tlb, double xm, double alpha, double *tz, double zlb, double s2, int mn1, double * zn1)
	{
		string function_name = "t_gfatmos_nrlmsise00::densu";
		try
		{
			double yd2, yd1, x = 0, y;
			double rgas = 831.4;
			double densu_temp = 1.0;
			double za, z, zg2, tt, ta;
			double dta, z1 = 0, z2, t1 = 0, t2, zg, zgdif = 0;
			int mn = 0;
			int k;
			double glb;
			double expl;
			double yi;
			double densa;
			double gamma, gamm;
			double xs[5], ys[5], y2out[5];
			/* joining altitudes of Bates and spline */
			za = zn1[0];
			if (alt>za)
				z = alt;
			else
				z = za;

			/* geopotential altitude difference from ZLB */
			zg2 = zeta(z, zlb);

			/* Bates temperature */
			tt = tinf - (tinf - _tlb) * exp(-s2*zg2);
			ta = tt;
			*tz = tt;
			densu_temp = *tz;

			if (alt < za)  /* calculate temperature below ZA  temperature gradient at ZA from Bates profile */
			{
				dta = (tinf - ta)*s2*((re + zlb) / (re + za))*((re + zlb) / (re + za));
				meso_tgn1[0] = dta;
				meso_tn1[0] = ta;
				z = (alt > zn1[mn1 - 1]) ? alt : zn1[mn1 - 1];
				mn = mn1;
				z1 = zn1[0];
				z2 = zn1[mn - 1];
				t1 = meso_tn1[0];
				t2 = meso_tn1[mn - 1];
				/* geopotental difference from z1 */
				zg = zeta(z, z1);
				zgdif = zeta(z2, z1);
				/* set up spline nodes */
				for (k = 0; k < mn; k++)
				{
					xs[k] = zeta(zn1[k], z1) / zgdif;
					ys[k] = 1.0 / meso_tn1[k];
				}
				/* end node derivatives */
				yd1 = -meso_tgn1[0] / (t1*t1) * zgdif;
				yd2 = -meso_tgn1[1] / (t2*t2) * zgdif * pow(((re + z2) / (re + z1)), 2.0);
				/* calculate spline coefficients */
				spline(xs, ys, mn, yd1, yd2, y2out);
				x = zg / zgdif;
				splint(xs, ys, y2out, mn, x, &y);
				/* temperature at altitude */
				*tz = 1.0 / y;
				densu_temp = *tz;
			}
			if (xm == 0)
				return densu_temp;

			/* calculate density above za */
			glb = gsurf / pow((1.0 + zlb / re), 2.0);
			gamma = xm * glb / (s2 * rgas * tinf);
			expl = exp(-s2 * gamma * zg2);
			if (expl>50.0)
				expl = 50.0;
			if (tt <= 0)
				expl = 50.0;

			/* density at altitude */
			densa = dlb * pow((_tlb / tt), ((1.0 + alpha + gamma))) * expl;
			densu_temp = densa;
			if (alt >= za)
				return densu_temp;

			/* calculate density below za */
			glb = gsurf / pow((1.0 + z1 / re), 2.0);
			gamm = xm * glb * zgdif / rgas;
				
			/* integrate spline temperatures */
			splini(xs, ys, y2out, mn, x, &yi);
			expl = gamm * yi;
			if (expl>50.0)
				expl = 50.0;
			if (*tz <= 0)
				expl = 50.0;

			/* density at altitude */
			densu_temp = densu_temp * pow((t1 / *tz), (1.0 + alpha)) * exp(-expl);
			return densu_temp;
		}
		catch (...)
		{
			write_log(_glog, 0, function_name, "ERROR : throw exception in densu function while calculate the nrlmsise00 density model");
			return 0.0;
		}
	}

	double t_gfatmos_nrlmsise00::dnet(double _dd, double dm, double zhm, double xmm, double xm)
	{
		string function_name = "t_gfatmos_nrlmsise00::dnet";
		try
		{
			double a, ylog;
			a = zhm / (xmm - xm);
			if (dm <= 0.0 || _dd <= 0.0) // goto 5
			{
				write_log(_glog, 0, function_name, "DNET LOG ERROR");
				if (_dd == 0.0&&dm == 0.0)
					_dd = 1.0;
				if (dm == 0)
					return _dd;
				if (_dd == 0)
					return dm;
			}
			ylog = a * log(dm / _dd);
			if (ylog<-10)
				return _dd;
			if (ylog>10)
				return dm;
			a = _dd*pow((1.0 + exp(ylog)), (1.0 / a));
			return a;
		}
		catch (...)
		{
			write_log(_glog, 0, function_name, "ERROR : throw exception in dnet function while calculate the nrlmsise00 density model");
			return 0.0;
		}
	}

	double t_gfatmos_nrlmsise00::globe7(double *p, int year, int doy, double sec, double alt, double g_lat, double g_long, double lst)
	{
		string function_name = "t_gfatmos_nrlmsise00::globe7";
		try
		{
			double t[15];
			int i, j;
			double apd;
			double tloc;
			double c, s, c2, c4, s2;
			double sr = 7.2722E-5;
			double dgtr = 1.74533E-2;
			double dr = 1.72142E-2;
			double hr = 0.2618;
			double cd32, cd18, cd14, cd39;
			double df;
			double f1, f2;
			double tinf;
			double tmp_ap[7];

			tloc = lst;
			for (j = 0; j<14; j++)
				t[j] = 0;
			/* calculate legendre polynomials */
			c = sin(g_lat * dgtr);
			s = cos(g_lat * dgtr);
			c2 = c*c;
			c4 = c2*c2;
			s2 = s*s;

			plg[0][1] = c;
			plg[0][2] = 0.5*(3.0*c2 - 1.0);
			plg[0][3] = 0.5*(5.0*c*c2 - 3.0*c);
			plg[0][4] = (35.0*c4 - 30.0*c2 + 3.0) / 8.0;
			plg[0][5] = (63.0*c2*c2*c - 70.0*c2*c + 15.0*c) / 8.0;
			plg[0][6] = (11.0*c*plg[0][5] - 5.0*plg[0][4]) / 6.0;
			/*      plg[0][7] = (13.0*c*plg[0][6] - 6.0*plg[0][5])/7.0; */
			plg[1][1] = s;
			plg[1][2] = 3.0*c*s;
			plg[1][3] = 1.5*(5.0*c2 - 1.0)*s;
			plg[1][4] = 2.5*(7.0*c2*c - 3.0*c)*s;
			plg[1][5] = 1.875*(21.0*c4 - 14.0*c2 + 1.0)*s;
			plg[1][6] = (11.0*c*plg[1][5] - 6.0*plg[1][4]) / 5.0;
			/*      plg[1][7] = (13.0*c*plg[1][6]-7.0*plg[1][5])/6.0; */
			/*      plg[1][8] = (15.0*c*plg[1][7]-8.0*plg[1][6])/7.0; */
			plg[2][2] = 3.0*s2;
			plg[2][3] = 15.0*s2*c;
			plg[2][4] = 7.5*(7.0*c2 - 1.0)*s2;
			plg[2][5] = 3.0*c*plg[2][4] - 2.0*plg[2][3];
			plg[2][6] = (11.0*c*plg[2][5] - 7.0*plg[2][4]) / 4.0;
			plg[2][7] = (13.0*c*plg[2][6] - 8.0*plg[2][5]) / 5.0;
			plg[3][3] = 15.0*s2*s;
			plg[3][4] = 105.0*s2*s*c;
			plg[3][5] = (9.0*c*plg[3][4] - 7.*plg[3][3]) / 2.0;
			plg[3][6] = (11.0*c*plg[3][5] - 8.*plg[3][4]) / 3.0;

			if (!(((sw[7] == 0) && (sw[8] == 0)) && (sw[14] == 0))) {
				stloc = sin(hr*tloc);
				ctloc = cos(hr*tloc);
				s2tloc = sin(2.0*hr*tloc);
				c2tloc = cos(2.0*hr*tloc);
				s3tloc = sin(3.0*hr*tloc);
				c3tloc = cos(3.0*hr*tloc);
			}

			cd32 = cos(dr*(doy - p[31]));
			cd18 = cos(2.0*dr*(doy - p[17]));
			cd14 = cos(dr*(doy - p[13]));
			cd39 = cos(2.0*dr*(doy - p[38]));

			/* F10.7 EFFECT */
			df = _f107 - _f107a;
			dfa = _f107a - 150.0;
			t[0] = p[19] * df*(1.0 + p[59] * dfa) + p[20] * df*df + p[21] * dfa + p[29] * pow(dfa, 2.0);
			f1 = 1.0 + (p[47] * dfa + p[19] * df + p[20] * df*df)*swc[1];
			f2 = 1.0 + (p[49] * dfa + p[19] * df + p[20] * df*df)*swc[1];

			/*  TIME INDEPENDENT */
			t[1] = (p[1] * plg[0][2] + p[2] * plg[0][4] + p[22] * plg[0][6]) + \
				(p[14] * plg[0][2])*dfa*swc[1] + p[26] * plg[0][1];

			/*  SYMMETRICAL ANNUAL */
			t[2] = p[18] * cd32;

			/*  SYMMETRICAL SEMIANNUAL */
			t[3] = (p[15] + p[16] * plg[0][2])*cd18;

			/*  ASYMMETRICAL ANNUAL */
			t[4] = f1*(p[9] * plg[0][1] + p[10] * plg[0][3])*cd14;

			/*  ASYMMETRICAL SEMIANNUAL */
			t[5] = p[37] * plg[0][1] * cd39;

			/* DIURNAL */
			if (sw[7]) 
			{
				double t71, t72;
				t71 = (p[11] * plg[1][2])*cd14*swc[5];
				t72 = (p[12] * plg[1][2])*cd14*swc[5];
				t[6] = f2*((p[3] * plg[1][1] + p[4] * plg[1][3] + p[27] * plg[1][5] + t71) * \
					ctloc + (p[6] * plg[1][1] + p[7] * plg[1][3] + p[28] * plg[1][5] \
						+ t72)*stloc);
			}

			/* SEMIDIURNAL */
			if (sw[8]) 
			{
				double t81, t82;
				t81 = (p[23] * plg[2][3] + p[35] * plg[2][5])*cd14*swc[5];
				t82 = (p[33] * plg[2][3] + p[36] * plg[2][5])*cd14*swc[5];
				t[7] = f2*((p[5] * plg[2][2] + p[41] * plg[2][4] + t81)*c2tloc + (p[8] * plg[2][2] + p[42] * plg[2][4] + t82)*s2tloc);
			}

			/* TERDIURNAL */
			if (sw[14]) 
			{
				t[13] = f2 * ((p[39] * plg[3][3] + (p[93] * plg[3][4] + p[46] * plg[3][6])*cd14*swc[5])* s3tloc + (p[40] * plg[3][3] + (p[94] * plg[3][4] + p[48] * plg[3][6])*cd14*swc[5])* c3tloc);
			}

			/* magnetic activity based on daily ap */
			if (sw[9] == -1) 
			{
				for (int k = 0; k < 7; k++)
					tmp_ap[k] = _ap[k];
				if (p[51] != 0) 
				{
					double exp1;
					exp1 = exp(-10800.0*sqrt(p[51] * p[51]) / (1.0 + p[138] * (45.0 - sqrt(g_lat*g_lat))));
					if (exp1>0.99999)
						exp1 = 0.99999;
					if (p[24]<1.0E-4)
						p[24] = 1.0E-4;
					apt[0] = sg0(exp1, p, tmp_ap);
					/* apt[1]=sg2(exp1,p,ap->a);
					apt[2]=sg0(exp2,p,ap->a);
					apt[3]=sg2(exp2,p,ap->a);
					*/
					if (sw[9]) 
					{
						t[8] = apt[0] * (p[50] + p[96] * plg[0][2] + p[54] * plg[0][4] + \
							(p[125] * plg[0][1] + p[126] * plg[0][3] + p[127] * plg[0][5])*cd14*swc[5] + \
							(p[128] * plg[1][1] + p[129] * plg[1][3] + p[130] * plg[1][5])*swc[7] * \
							cos(hr*(tloc - p[131])));
					}
				}
			}
			else 
			{
				double p44, p45;
				apd = _dap - 4.0;
				p44 = p[43];
				p45 = p[44];
				if (p44<0)
					p44 = 1.0E-5;
				apdf = apd + (p45 - 1.0)*(apd + (exp(-p44 * apd) - 1.0) / p44);
				if (sw[9]) {
					t[8] = apdf*(p[32] + p[45] * plg[0][2] + p[34] * plg[0][4] + \
						(p[100] * plg[0][1] + p[101] * plg[0][3] + p[102] * plg[0][5])*cd14*swc[5] +
						(p[121] * plg[1][1] + p[122] * plg[1][3] + p[123] * plg[1][5])*swc[7] *
						cos(hr*(tloc - p[124])));
				}
			}

			if ((sw[10]) && (g_long > -1000.0))
			{
				/* longitudinal */
				if (sw[11])
				{
					t[10] = (1.0 + p[80] * dfa*swc[1])* \
						((p[64] * plg[1][2] + p[65] * plg[1][4] + p[66] * plg[1][6]\
							+ p[103] * plg[1][1] + p[104] * plg[1][3] + p[105] * plg[1][5]\
							+ swc[5] * (p[109] * plg[1][1] + p[110] * plg[1][3] + p[111] * plg[1][5])*cd14)* \
							cos(dgtr*g_long) \
							+ (p[90] * plg[1][2] + p[91] * plg[1][4] + p[92] * plg[1][6]\
								+ p[106] * plg[1][1] + p[107] * plg[1][3] + p[108] * plg[1][5]\
								+ swc[5] * (p[112] * plg[1][1] + p[113] * plg[1][3] + p[114] * plg[1][5])*cd14)* \
							sin(dgtr*g_long));
				}

				/* ut and mixed ut, longitude */
				if (sw[12])
				{
					t[11] = (1.0 + p[95] * plg[0][1])*(1.0 + p[81] * dfa*swc[1])*\
						(1.0 + p[119] * plg[0][1] * swc[5] * cd14)*\
						((p[68] * plg[0][1] + p[69] * plg[0][3] + p[70] * plg[0][5])*\
							cos(sr*(sec - p[71])));
					t[11] += swc[11] * \
						(p[76] * plg[2][3] + p[77] * plg[2][5] + p[78] * plg[2][7])*\
						cos(sr*(sec - p[79]) + 2.0*dgtr*g_long)*(1.0 + p[137] * dfa*swc[1]);
				}

				/* ut, longitude magnetic activity */
				if (sw[13])
				{
					if (sw[9] == -1) {
						if (p[51]) {
							t[12] = apt[0] * swc[11] * (1. + p[132] * plg[0][1])*\
								((p[52] * plg[1][2] + p[98] * plg[1][4] + p[67] * plg[1][6])*\
									cos(dgtr*(g_long - p[97])))\
								+ apt[0] * swc[11] * swc[5] * \
								(p[133] * plg[1][1] + p[134] * plg[1][3] + p[135] * plg[1][5])*\
								cd14*cos(dgtr*(g_long - p[136])) \
								+ apt[0] * swc[12] * \
								(p[55] * plg[0][1] + p[56] * plg[0][3] + p[57] * plg[0][5])*\
								cos(sr*(sec - p[58]));
						}
					}
					else
					{
						t[12] = apdf*swc[11] * (1.0 + p[120] * plg[0][1])*\
							((p[60] * plg[1][2] + p[61] * plg[1][4] + p[62] * plg[1][6])*\
								cos(dgtr*(g_long - p[63])))\
							+ apdf*swc[11] * swc[5] * \
							(p[115] * plg[1][1] + p[116] * plg[1][3] + p[117] * plg[1][5])* \
							cd14*cos(dgtr*(g_long - p[118])) \
							+ apdf*swc[12] * \
							(p[83] * plg[0][1] + p[84] * plg[0][3] + p[85] * plg[0][5])* \
							cos(sr*(sec - p[75]));
					}
				}
			}

			/* parms not used: 82, 89, 99, 139-149 */
			tinf = p[30];
			for (i = 0; i<14; i++)
				tinf = tinf + fabs(sw[(int64_t)i + 1])*t[i];
			return tinf;
		}
		catch (...)
		{
			write_log(_glog, 0, function_name, "ERROR : throw exception in globe7 function while calculate the nrlmsise00 density model");
		}
	}

	bool t_gfatmos_nrlmsise00::spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
	{
		string function_name = "t_gfatmos_nrlmsise00::spline";
		try
		{
			vector<double> u;
			double sig, p, qn, un;
			int i, k;
			
			if (yp1>0.99e30) 
			{
				y2[0] = 0;
				u.push_back(0);
			}
			else 
			{
				y2[0] = -0.5;
				u.push_back((3.0 / (x[1] - x[0]))*((y[1] - y[0]) / (x[1] - x[0]) - yp1));
			}
			for (i = 1; i<(n - 1); i++) 
			{
				sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
				p = sig * y2[i - 1] + 2.0;
				y2[i] = (sig - 1.0) / p;
				u.push_back((6.0 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1])) / (x[i + 1] - x[i - 1]) - sig * u[(int64_t)i - 1]) / p);
			}
			if (ypn>0.99e30) 
			{
				qn = 0;
				un = 0;
			}
			else 
			{
				qn = 0.5;
				un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
			}
			y2[n - 1] = (un - qn * u[(int64_t)n - 2]) / (qn * y2[n - 2] + 1.0);
			for (k = n - 2; k >= 0; k--)
				y2[k] = y2[k] * y2[k + 1] + u[k];

			return true;
		}
		catch (...)
		{
			write_log(_glog, 0, function_name, "ERROR : throw exception in spline function while calculate the nrlmsise00 density model");
			return false;
		}
	}

	bool t_gfatmos_nrlmsise00::splint(double *xa, double *ya, double *y2a, int n, double x, double *y)
	{
		string function_name = "t_gfatmos_nrlmsise00::splint";
		try
		{
			int klo = 0;
			int khi = n - 1;
			int k;
			double h;
			double a, b, yi;
			while ((khi - klo)>1) 
			{
				k = (khi + klo) / 2;
				if (xa[k]>x)
					khi = k;
				else
					klo = k;
			}
			h = xa[khi] - xa[klo];
			if (h == 0.0)
				write_log(_glog, 0, function_name, "bad XA input to splint");
			a = (xa[khi] - x) / h;
			b = (x - xa[klo]) / h;
			yi = a * ya[klo] + b * ya[khi] + ((a*a*a - a) * y2a[klo] + (b*b*b - b) * y2a[khi]) * h * h / 6.0;
			*y = yi;
			return true;
		}
		catch (...)
		{
			write_log(_glog, 0, function_name, "ERROR : throw exception in splint function while calculate the nrlmsise00 density model");
			return false;
		}
	}

	bool t_gfatmos_nrlmsise00::splini(double *xa, double *ya, double *y2a, int n, double x, double *y)
	{
		string function_name = "t_gfatmos_nrlmsise00::splini";
		try
		{
			double yi = 0;
			int klo = 0;
			int khi = 1;
			double xx, h, a, b, a2, b2;
			while ((x>xa[klo]) && (khi<n)) 
			{
				xx = x;
				if (khi<(n - 1)) {
					if (x<xa[khi])
						xx = x;
					else
						xx = xa[khi];
				}
				h = xa[khi] - xa[klo];
				a = (xa[khi] - xx) / h;
				b = (xx - xa[klo]) / h;
				a2 = a*a;
				b2 = b*b;
				yi += ((1.0 - a2) * ya[klo] / 2.0 + b2 * ya[khi] / 2.0 + ((-(1.0 + a2*a2) / 4.0 + a2 / 2.0) * y2a[klo] + (b2*b2 / 4.0 - b2 / 2.0) * y2a[khi]) * h * h / 6.0) * h;
				klo++;
				khi++;
			}
			*y = yi;
			return true;
		}
		catch (...)
		{
			write_log(_glog, 0, function_name, "ERROR : throw exception in splini function while calculate the nrlmsise00 density model");
			return false;
		}
	}

	bool t_gfatmos_nrlmsise00::xyz2blh(ColumnVector &crd, ColumnVector &ell)
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

	double t_gfatmos_nrlmsise00::ccor(double alt, double r, double h1, double zh)
	{
		double e;
		double ex;
		e = (alt - zh) / h1;
		if (e>70)
			return exp(0);
		if (e<-70)
			return exp(r);
		ex = exp(e);
		e = r / (1.0 + ex);
		return exp(e);
	}

	double t_gfatmos_nrlmsise00::ccor2(double alt, double r, double h1, double zh, double h2)
	{
		double e1, e2;
		double ex1, ex2;
		double ccor2v;
		e1 = (alt - zh) / h1;
		e2 = (alt - zh) / h2;
		if ((e1 > 70) || (e2 > 70))
			return exp(0);
		if ((e1 < -70) && (e2 < -70))
			return exp(r);
		ex1 = exp(e1);
		ex2 = exp(e2);
		ccor2v = r / (1.0 + 0.5 * (ex1 + ex2));
		return exp(ccor2v);
	}

	double t_gfatmos_nrlmsise00::sg0(double ex, double *p, double *t_ap)
	{
		double sumex;
		sumex = 1.0 + (1.0 - pow(ex, 19.0)) / (1.0 - ex) * pow(ex, 0.5);
		return (g0(t_ap[1], p) + (g0(t_ap[2], p)*ex + g0(t_ap[3], p)*ex*ex + \
			g0(t_ap[4], p)*pow(ex, 3.0) + (g0(t_ap[5], p)*pow(ex, 4.0) + \
				g0(t_ap[6], p)*pow(ex, 12.0))*(1.0 - pow(ex, 8.0)) / (1.0 - ex))) / sumex;
	}

	double t_gfatmos_nrlmsise00::g0(double a, double *p)
	{
		return (a - 4.0 + (p[25] - 1.0) * (a - 4.0 + (exp(-sqrt(p[24] * p[24]) * (a - 4.0)) - 1.0) / sqrt(p[24] * p[24])));
	}

	double t_gfatmos_nrlmsise00::zeta(double zz, double zl)
	{
		return ((zz - zl)*(re + zl) / (re + zz));
	}
}