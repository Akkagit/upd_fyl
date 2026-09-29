/**
*
* @verbatim
History
-2.0         Liu Chengbo              2020-03-03 Creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University.All rights reserved.
*
* @file		      gfearthradiation2.h
* @brief		    The class for earth radiation model.
*
* @author          Liu Chengbo, Wuhan University
* @version		    2.0.0
* @date		    2020-03-03
*
*/

#include "gforce/gfearthradiation_nav.h"
#include "gutils/ginfolog.h"

using namespace std;

namespace great
{
	/** constructor. */
	t_gfearthradiation_nav::t_gfearthradiation_nav(t_init_arg* arg, t_glog* log) : t_gforce(log)
	{
		string function_name = "t_gfearthradiation2::t_gfearthradiation2";
		if (!log)
		{
			write_log(log, 0, function_name,
				"ERROR : the log point is empty");
		}

		t_init_earthradiation2* earad = NULL;
		earad = dynamic_cast<t_init_earthradiation2*>(arg);

		if (!earad)
		{
			write_log(log, 4, function_name, "ERROR : init failed ");
			if (!earad)
			{
				write_log(log, 4, function_name, "ERROR : no init argument ");
			}
			if (!earad->galbedo)
			{
				write_log(log, 4, function_name, "ERROR : no albedo data ");
			}
			/*if (!earad->gpannel)
			{
				write_log(log, 4, function_name, "ERROR : no pannel data ");
			}*/
		}
		else
		{
			// set the init value
			_galbedo = earad->galbedo;
			_model = earad->model;
			v_ns_all1.resize(72, 144);
			v_ns_all2.resize(72, 144);
			v_ns_all3.resize(72, 144);
			d_area_all.resize(72, 144);
			cergem.resize(72, 144);
			cergre.resize(72, 144);
			toa = 6371000.0 + 30000.0;
			pi = 4.0 * atan(1.0);
			Xatt.resize(3, 1);
			Yatt.resize(3, 1);
			Zatt.resize(3, 1);
			Datt.resize(3, 1);
		}
	}

	/**
	* defalut deconstructor.
	*/
	t_gfearthradiation_nav::~t_gfearthradiation_nav()
	{

	}
	/**value satellite panel.
	* @param[in] block satellite block name*/
	
	bool t_gfearthradiation_nav::cal_sat_panel(string block)
	{
		string function_name = "t_gfearthradiation::calsatpanel";
		try
		{
			const char* blockc = block.c_str();
			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < 2; k++)
				{
					area[j][k] = 0;
					area2[j][k] = 0;
					difu[j][k] = 0;
					difu2[j][k] = 0;
					refl[j][k] = 0;
					refl2[j][k] = 0;
					difuir[j][k] = 0;
					reflir[j][k] = 0;
					reflir2[j][k] = 0;
					difuir2[j][k] = 0;
				}
			}
			if (strcmp(blockc,"BLOCKIIF")==0)
			{
				area[0][0] = 5.40; area[0][1] = 5.40;
				area[1][0] = 7.01; area[1][1] = 7.01;
				area[2][0] = 5.72; area[2][1] = 5.72;
				area[3][0] = 22.25; area[3][1] = 22.25;
				difu[1][0] = 0.448; difu[1][1] = 0.448;
				difu[2][0] = 0.448; difu[2][1] = 0.448;
				difu[3][0] = 0.0345; difu[3][1] = 0.055;
				refl[1][0] = 0.112; refl[1][1] = 0.112;
				refl[2][0] = 0.112; refl[2][1] = 0.112;
				refl[3][0] = 0.1955; refl[3][1] = 0.055;
				for (int j = 0; j < 3; j++)
				{
					reflir[j][0] = 0.1;
					difuir[j][0] = 0.1;
					reflir[j][1] = 0.1;
					difuir[j][1] = 0.1;
				}
			}
			else if (strcmp(blockc, "BLOCKIIR-M") == 0|| strcmp(blockc, "BLOCKIIR-A") == 0|| strcmp(blockc, "BLOCKIIR-B") == 0)
			{
				area[0][0] = 4.25; area[0][1] = 4.25;
				area[1][0] = 3.40; area[1][0] = 3.40;
				area[2][0] = 3.05; area[2][1] = 3.05;
				area[3][0] = 13.92; area[3][1] = 13.92;
				refl[3][0] = 0.2491; refl[3][1] = 0.055;
				difu[0][0] = 0.06; difu[0][1] = 0.06;
				difu[1][0] = 0.06; difu[1][1] = 0.06;
				difu[2][0] = 0.06; difu[2][1] = 0.06;
				difu[3][0] = 0.04397; difu[3][1] = 0.055;
				area2[1][0] = 1.06; area2[1][1] = 1.06;
				area2[2][0] = 1.06; area2[2][1] = 1.06;
				difu2[1][0] = 0.06; difu2[1][1] = 0.06;
				difu2[2][0] = 0.06; difu2[2][1] = 0.06;
				for (int j = 0; j < 3; j++)
				{
					reflir[j][0] = 0.1;
					difuir[j][0] = 0.1;
					reflir[j][1] = 0.1;
					difuir[j][1] = 0.1;
				}
			}
			else if (strcmp(blockc, "BEIDOU-2I") == 0 || strcmp(blockc, "BEIDOU-2M") == 0) //BDS-2 MEO IGSO
			{
				area[0][0] = 3.440; area[0][1] = 3.440;
				area[1][0] = 4.400; area[1][1] = 7.01;
				area[2][0] = 3.784; area[2][1] = 3.784;
				area[3][0] = 22.680; area[3][1] = 22.680;
				difu[0][0] = 0.679; difu[0][1] = 0.640;
				difu[1][0] = 0.317; difu[1][1] = 0.317;
				difu[2][0] = 0.640; difu[2][1] = 0.640;
				difu[3][0] = 0.040; difu[3][1] = 0.040;
				refl[1][0] = 0.439; refl[1][1] = 0.431;
				refl[3][0] = 0.240; refl[3][1] = 0.240;
				for (int j = 0; j < 3; j++)
				{
					reflir[j][0] = 0.1;
					difuir[j][0] = 0.1;
					reflir[j][1] = 0.1;
					difuir[j][1] = 0.1;
				}
			}
			else if (strcmp(blockc, "GALILEO-1") == 0)  //Galileo IOV
			{
				area[0][0] = 3.000; area[0][1] = 3.000;
				area[1][0] = 3.000; area[1][1] = 3.000;
				area[2][0] = 1.320; area[2][1] = 1.320;
				area[3][0] = 10.820; area[3][1] = 10.820;
				difu[0][0] = 0.124; difu[0][1] = 0.060;
				difu[1][0] = 0.139; difu[1][1] = 0.140;
				difu[2][0] = 0.060; difu[2][1] = 0.131;
				refl[0][0] = 0.094; 
				refl[1][0] = 0.473; refl[1][1] = 0.480;
				refl[2][0] = 0.000; refl[2][1] = 0.425;
				refl[3][0] = 0.086; refl[3][1] = 0.086;
				for (int j = 0; j < 3; j++)
				{
					reflir[j][0] = 0.1;
					difuir[j][0] = 0.1;
					reflir[j][1] = 0.1;
					difuir[j][1] = 0.1;
				}
			}
			else if (strcmp(blockc, "GALILEO-2") == 0)  //Galileo FOC
			{
				area[0][0] = 3.036; area[0][1] = 3.036;
				area[1][0] = 2.783; area[1][1] = 2.783;
				area[2][0] = 1.320; area[2][1] = 1.320;
				area[3][0] = 10.820; area[3][1] = 10.820;
				difu[0][0] = 0.162; difu[0][1] = 0.108;
				difu[1][0] = 0.136; difu[1][1] = 0.141;
				difu[2][0] = 0.070; difu[2][1] = 0.150;
				refl[0][0] = 0.143; refl[0][1] = 0.231;
				refl[1][0] = 0.404; refl[1][1] = 0.434;
				refl[2][0] = 0.000; refl[2][1] = 0.487;
				refl[3][0] = 0.086; refl[3][1] = 0.086;
				for (int j = 0; j < 3; j++)
				{
					reflir[j][0] = 0.1;
					difuir[j][0] = 0.1;
					reflir[j][1] = 0.1;
					difuir[j][1] = 0.1;
				}
			}
			else
			{
				write_log(_glog, 4, function_name, "ERROR : invalid satellite type");
				return false;
			}
			return true;
		}
		catch (...)
		{
			write_log(_glog, 4, function_name, "ERROR : throw exception while calculating erp satpanel");
			return false;
		}
	}
	/**value latkmx/lonkmx.
	* @param[in] grd the grid
	*/
	bool t_gfearthradiation_nav::cal_max(int grd)
	{
		string function_name = "t_gfearthradiation::calmax";
		try
		{
			if (grd == 1)
			{
				grdnum = 1.0;
				latkmx = 72;
				lonkmx = 144;
			}
			else if (grd == 2)
			{
				grdnum = 2.0;
				grdcer = 2;
				latkmx = 36;
				lonkmx = 72;
			}
			else if (grd == 3)
			{
				grdnum = 4.0;
				grdcer = 4;
				latkmx = 18;
				lonkmx = 36;
			}
			return true;
		}
		catch (...)
		{
			write_log(_glog, 4, function_name, "ERROR : throw exception while calculating atmosphere acc");
			return false;
		}
	}
	/**calculate the reflection/emission properties value.
	* @param[in] grd          the grid
	* @param[in] ceresdata_re the reflection properties value
	* @param[in] ceresdata_em the emission properties value
	*/
	bool t_gfearthradiation_nav::cal_cerg_re_em(int grd, Matrix& ceresdata_re, Matrix& ceresdata_em)
	{
		double grdang = 2.5;
		double latind = 36.5;
		double lonind = 72.5;
		double d_ang, lat_in, lon_in, coslat;
		int latk, lonk;

		string function_name = "t_gfearthradiation::calcerg_re_em";
		try
		{
			grdang = grdang * grdnum;
			latind = (latind - 0.5) / grdnum + 0.5;
			lonind = (lonind - 0.5) / grdnum + 0.5;
			d_ang = (pi * grdang / 180.0) * (pi * grdang / 180.0);
			for (latk = 1; latk <= latkmx; latk++)
			{
				for (lonk = 1; lonk <= lonkmx; lonk++)
				{
					lat_in = (latk - latind) * grdang * (pi / 180.0);
					lon_in = (lonk - lonind) * grdang * (pi / 180.0);
					coslat = cos(lat_in);
					d_area_all(latk, lonk) = toa * toa * coslat * d_ang;
					v_ns_all1(latk, lonk) = coslat * cos(lon_in);
					v_ns_all2(latk, lonk) = coslat * sin(lon_in);
					v_ns_all3(latk, lonk) = sin(lat_in);
					cergem(latk, lonk) = 0.0;
					cergre(latk, lonk) = 0.0;
					if (grd == 1)
					{
						cergre(latk, lonk) = ceresdata_re(latk, lonk);
						cergem(latk, lonk) = ceresdata_em(latk, lonk);
					}
					else if (grd == 2 || grd == 3)
					{
						for (int ii = 0; ii < grdcer; ii++)
						{
							for (int jj = 0; jj < grdcer; jj++)
							{
								cergre(latk, lonk) = cergre(latk, lonk) + ceresdata_re(grdcer * latk - ii, grdcer * lonk - jj);
								cergem(latk, lonk) = cergem(latk, lonk) + ceresdata_em(grdcer * latk - ii, grdcer * lonk - jj);
							}
						}
						cergre(latk, lonk) = cergre(latk, lonk) / (grdnum * grdnum);
						cergem(latk, lonk) = cergem(latk, lonk) / (grdnum * grdnum);
					}
				}
			}
			return true;
		}
		catch (...)
		{
			write_log(_glog, 4, function_name, "ERROR : throw exception while calculating cerg_re_em");
			return false;
		}
	}

	/**calculate the accleration caused by earth radiation.
	* @param[in] argument Parameters to be used in calculation
	*/
	bool t_gfearthradiation_nav::calculate(t_gforce_calarg* argument)
	{
		string function_name = "t_gfearthradiation2::calcearthradmod";
		try
		{
			t_gfcalarg_earthradiation* arg = NULL;
			arg = dynamic_cast<t_gfcalarg_earthradiation*>(argument);
			if (!arg)
			{
				write_log(_glog, 4, function_name, "ERROR : no calculate argument");
				return false;
			}
			if (arg->e2j->is_zero() || arg->e2j->nrows() != 3 || arg->e2j->ncols() != 3)
			{
				write_log(_glog, 4, function_name, "ERROR : no rot matrix for TRS to CRS");
				return false;
			}
			double abspos, abssun, rsun, absy0v, psidot, psi, albfac, phasevi, phaseir, absncfvi, absncfir;
			double facabdi, facrefl;
			double eau = 149597870691.0;
			double s0 = 1367.0;
			double alb = 0.3;
			double c = 299792458.0;
			double cosa = 0.0;
			double isr = 0.10;
			double idr = 0.10;
			int erm, grd, reff;
			ColumnVector e_sun;
			Matrix kappa;
			ColumnVector radvec, ncfvec, radvel;
			double crs_vec[3], alg_vec[3];
			ColumnVector side(1, 4);
			ColumnVector cosang(0.0, 4);
			ColumnVector forcevis(0.0, 3);
			ColumnVector forceir(0.0, 3);
			ColumnVector force(0.0, 3);
			ColumnVector xtmp(0.0, 3);
			ColumnVector forceACR(0.0, 3);
			double ATTSURF[3][4];
			double norvec[3][4];
			double tabsp[3], tabsp2[3], trefl[3], trefl2[3], tdifu[3], tdifu2[3];
			double prtflt[3][3], prtcyl[3][3], prtfltir[3][3], prtcylir[3][6];
			double y0sat[3];
			double abscrs;
			double absvel;
			ysat = *(arg->sat_pos);
			vsat = *(arg->sat_vel);
			blockname = arg->sat_blocknum;
			cal_sat_panel(blockname);
			sun = *(arg->sun_pos_crs);
			sun = sun * 1000.0;
			const char* cmodel = _model.c_str();
			if (strcmp(cmodel, "ERM_A") == 0)
			{
				erm = 1;
			}
			else 
			{
				erm = 2;
			}
			//printf("%s,%d\n", cmodel, erm);
			grd = 1;
			reff = 0;
			kappa = (*arg->e2j).t();
			_mass = (arg->sat_mass);
			ysat = ysat * 1000.0;
			vsat = vsat * 1000.0;
			absvel = sqrt(vsat[0] * vsat[0] + vsat[1] * vsat[1] + vsat[2] * vsat[2]);		
			radvel = vsat / absvel;
			abspos = sqrt(ysat[0] * ysat[0] + ysat[1] * ysat[1] + ysat[2] * ysat[2]);
			radvec = ysat / abspos;
			crs_vec[0] = ysat[1] * vsat[2] - ysat[2] * vsat[1];
			crs_vec[1] = ysat[2] * vsat[0] - ysat[0] * vsat[2];
			crs_vec[2] = ysat[0] * vsat[1] - ysat[1] * vsat[0];
			abscrs = sqrt(crs_vec[0] * crs_vec[0] + crs_vec[1] * crs_vec[1] + crs_vec[2] * crs_vec[2]);
			for (int k = 0; k < 3; k++)
			{
				crs_vec[k] = crs_vec[k] / abscrs;
			}
			alg_vec[0] = crs_vec[1] * radvec[2] - crs_vec[2] * radvec[1];
			alg_vec[1] = crs_vec[2] * radvec[0] - crs_vec[0] * radvec[2];
			alg_vec[2] = crs_vec[0] * radvec[1] - crs_vec[1] * radvec[0];
			rsun = sqrt(pow(ysat[0] - sun[0], 2) + pow(ysat[1] - sun[1], 2) + pow(ysat[2] - sun[2], 2));
			if (erm > 0)
			{
				abssun = sqrt(sun[0] * sun[0] + sun[1] * sun[1] + sun[2] * sun[2]);
				e_sun = sun / abssun;
				psidot = e_sun[0] * radvec[0] + e_sun[1] * radvec[1] + e_sun[2] * radvec[2];
				if (fabs(psidot) > 1 - 1.0e-6)
				{
					psi = 0.0;
				}
				else
				{
					psi = acos(psidot);
				}
				for (int k = 0; k < 3; k++)
				{
					Datt[k] = (sun[k] - ysat[k]) / rsun;
					Zatt[k] = -ysat[k] / abspos;
				}
				y0sat[0] = Zatt[1] * Datt[2] - Zatt[2] * Datt[1];
				y0sat[1] = Zatt[2] * Datt[0] - Zatt[0] * Datt[2];
				y0sat[2] = Zatt[0] * Datt[1] - Zatt[1] * Datt[0];
				absy0v = sqrt(y0sat[0] * y0sat[0] + y0sat[1] * y0sat[1] + y0sat[2] * y0sat[2]);
				for (int k = 0; k < 3; k++)
				{
					Yatt[k] = y0sat[k] / absy0v;
				}
				Xatt[0] = Yatt[1] * Zatt[2] - Yatt[2] * Zatt[1];
				Xatt[1] = Yatt[2] * Zatt[0] - Yatt[0] * Zatt[2];
				Xatt[2] = Yatt[0] * Zatt[1] - Yatt[1] * Zatt[0];
				s0 = s0 * (eau / abssun) * (eau / abssun);
				for (int k = 0; k < 3; k++)
				{
					ATTSURF[k][0] = Zatt[k];
					ATTSURF[k][1] = Yatt[k];
					ATTSURF[k][2] = Xatt[k];
					ATTSURF[k][3] = Datt[k];
				}
			}
			/*ColumnVector Enop;
			Enop = crossproduct(radvec, radvel);
			double absEnop = sqrt(Enop[0] * Enop[0] + Enop[1] * Enop[1] + Enop[2] * Enop[2]);
			Enop = Enop / absEnop;
			double sinb = dotproduct(Enop, e_sun);
			double beita = 90 - acos(sinb) * 180 / pi;
			double cosb = sqrt(1 - pow(sinb, 2));
			ColumnVector Enod;
			Enod = crossproduct(e_sun, Enop);
			double absEnod = sqrt(Enod[0] * Enod[0] + Enod[1] * Enod[1] + Enod[2] * Enod[2]);
			Enod = Enod / absEnod;
			ColumnVector Emid;
			Emid = crossproduct(Enod, Enop);*/
			if (erm == 1)
			{
				ncfvec = radvec;
				albfac = (pi * toa * toa) * (s0 / c) / (abspos * abspos);
				phasevi = (2 * alb / (3 * pi * pi)) * ((pi - psi) * cos(psi) + sin(psi));
				phaseir = (1 - alb) / (4 * pi);
				absncfvi = albfac * phasevi;
				absncfir = albfac * phaseir;

				for (int k = 0; k < 4; k++) 
				{
					cosang[k] = ATTSURF[0][k] * ncfvec[0] + ATTSURF[1][k] * ncfvec[1] + ATTSURF[2][k] * ncfvec[2];
					//printf("%f", cosang[k]);
					if (cosang[k] > 0)
					{
						side[k] = 2;
						norvec[0][k] = ATTSURF[0][k];
						norvec[1][k] = ATTSURF[1][k];
						norvec[2][k] = ATTSURF[2][k];
					}
					else if (cosang[k] <= 0)
					{
						side[k] = 1;
						norvec[0][k] = -ATTSURF[0][k];
						norvec[1][k] = -ATTSURF[1][k];
						norvec[2][k] = -ATTSURF[2][k];
					}
				}
				for (int k = 0; k < 4; k++)
				{
					double fviflt[3], firflt[3], fvi[3], fir[3], fvicyl[3], fircyl[3];
					facabdi = fabs(cosang[k]);
					facrefl = pow(cosang[k], 2);
					for (int j = 0; j < 3; j++)
					{
						tabsp[j] = facabdi * ncfvec[j];
						//printf("tabsp %.15f\n", tabsp[j]);
						tabsp2[j] = tabsp[j];
						trefl[j] = 2.0 * facrefl * norvec[j][k];
						//printf("trefl %.15f\n", trefl[j]);
						trefl2[j] = 4.0 / 3.0 * facrefl * norvec[j][k];
						//printf("trefl2 %.15f\n", trefl2[j]);
						tdifu[j] = facabdi * (ncfvec[j] + 2.0 / 3.0 * norvec[j][k]);
						//printf("tdifu %.15f\n", tdifu[j]);
						tdifu2[j] = facabdi * (ncfvec[j] + (pi / 6.0) * norvec[j][k]);
						//printf("tdifu2 %.15f\n\n", tdifu2[j]);
						if (k < 3)
						{
							int flag = side[k] - 1;
							//printf("%d\n", flag);
							prtflt[0][j] = tdifu[j] * (1 - refl[k][flag]);
							prtflt[1][j] = trefl[j] * refl[k][flag];
							prtflt[2][j] = 0.0;

							prtcyl[0][j] = tdifu2[j] * (1 - refl2[k][flag]);
							prtcyl[1][j] = trefl2[j] * refl2[k][flag];
							prtcyl[2][j] = 0.0;

							prtfltir[0][j] = tdifu[j] * (1 - isr);
							prtfltir[1][j] = trefl[j] * isr;
							prtfltir[2][j] = 0.0;

							prtcylir[0][j] = tdifu2[j] * (1 - isr);
							prtcylir[1][j] = trefl2[j] * isr;
							prtcylir[2][j] = 0.0;
							/*printf("%.15f %.15f %.15f\n", prtflt[0][j], prtflt[1][j], prtflt[2][j]);
							printf("%.15f %.15f %.15f\n\n", prtfltir[0][j], prtfltir[1][j], prtfltir[2][j]);*/
						}
						else  //else if(k==3)
						{
							int flag = side[k] - 1;
							prtflt[0][j] = tabsp[j] * (1 - refl[k][flag] - difu[k][flag]);
							prtflt[1][j] = trefl[j] * refl[k][flag];
							prtflt[2][j] = tdifu[j] * difu[k][flag];

							prtcyl[0][j] = tabsp2[j] * (1 - refl2[k][flag] - difu2[k][flag]);
							prtcyl[1][j] = trefl2[j] * refl2[k][flag];
							prtcyl[2][j] = tdifu2[j] * difu2[k][flag];
							
							prtfltir[0][j] = tabsp[j] * (1 - isr - idr);
							prtfltir[1][j] = trefl[j] * isr;
							prtfltir[2][j] = tdifu[j] * idr;

							prtcylir[0][j] = tabsp2[j] * (1 - isr - idr);
							prtcylir[1][j] = trefl2[j] * isr;
							prtcylir[2][j] = tdifu2[j] * idr;
							//printf("%.15f %.15f %.15f\n", prtflt[0][j], prtflt[1][j], prtflt[2][j]);
							//printf("%.15f %.15f %.15f\n\n", prtfltir[0][j], prtfltir[1][j], prtfltir[2][j]);
						}
						
						fviflt[j] = prtflt[0][j] + prtflt[1][j] + prtflt[2][j];
						firflt[j] = prtfltir[0][j] + prtfltir[1][j] + prtfltir[2][j];
						fvicyl[j] = prtcyl[0][j] + prtcyl[1][j] + prtcyl[2][j];
						fircyl[j] = prtcylir[0][j] + prtcylir[1][j] + prtcylir[2][j];

						fviflt[j] = fviflt[j] * area[k][1];
						firflt[j] = firflt[j] * area[k][1];
						fvicyl[j] = fvicyl[j] * area2[k][1];
						fircyl[j] = fircyl[j] * area2[k][1];
						fvi[j] = (fviflt[j] + fvicyl[j]) * absncfvi;
						fir[j] = (firflt[j] + fircyl[j]) * absncfir;
						force[j] = force[j] + fvi[j] + fir[j];
					}
				}	
			}
			else if (erm == 2)
			{
				double d_area, dist2, absdist, cos_in, cos_re;
				double refl_cf, emit_cf, e_refl, e_emit;
				ColumnVector v_ns(0.0, 3);
				ColumnVector v_ins(0.0, 3);
				ColumnVector v_dist(0.0, 3);
				ColumnVector v_sat(0.0, 3);
				if (!cal_max(grd))
				{
					write_log(_glog, 4, function_name, "ERROR : in t_gfearthradiation::calmax");
				}			
				cal_cerg_re_em(grd, _galbedo->get_ceresr_data(), _galbedo->get_cerese_data());
				for (int latk = 1; latk <= latkmx; latk++)
				{
					for (int lonk = 1; lonk <= lonkmx; lonk++)
					{
						d_area = d_area_all(latk, lonk);
						v_ns(1) = v_ns_all1(latk, lonk);
						v_ns(2) = v_ns_all2(latk, lonk);
						v_ns(3) = v_ns_all3(latk, lonk);
						for (int ii = 1; ii <= 3; ii++)
						{
							v_ins(ii) = 0.0;
							for (int jj = 1; jj <= 3; jj++)
							{
								v_ins(ii) = v_ins(ii) + kappa(jj, ii) * v_ns(jj);
							}
						}
						//v_ins += kappa*v_ns;
						v_dist = ysat - toa * v_ins;
						dist2 = v_dist(1) * v_dist(1) + v_dist(2) * v_dist(2) + v_dist(3) * v_dist(3);
						absdist = sqrt(dist2);
						v_sat = v_dist / absdist;
						cos_in = e_sun[0] * v_ins(1) + e_sun[1] * v_ins(2) + e_sun[2] * v_ins(3);
						cos_re = v_sat(1) * v_ins(1) + v_sat(2) * v_ins(2) + v_sat(3) * v_ins(3);
						if (cos_re >= 0)
						{
							refl_cf = cergre(latk, lonk);
							emit_cf = cergem(latk, lonk);

							if (cos_in >= 0)
							{
								e_refl = (refl_cf / (pi * dist2)) * cos_re * cos_in * s0 * d_area;
							}
							else
							{
								e_refl = 0.0;
							}
							e_emit = (emit_cf / (4 * pi * dist2)) * cos_re * s0 * d_area;
							absncfvi = e_refl / c;
							absncfir = e_emit / c;
							ncfvec = v_sat;

							for (int k = 0; k < 4; k++)
							{
								cosang[k] = ATTSURF[0][k] * ncfvec[0] + ATTSURF[1][k] * ncfvec[1] + ATTSURF[2][k] * ncfvec[2];
								//printf("%f", cosang[k]);
								if (cosang[k] > 0)
								{
									side[k] = 2;
									norvec[0][k] = ATTSURF[0][k];
									norvec[1][k] = ATTSURF[1][k];
									norvec[2][k] = ATTSURF[2][k];
								}
								else if (cosang[k] <= 0)
								{
									side[k] = 1;
									norvec[0][k] = -ATTSURF[0][k];
									norvec[1][k] = -ATTSURF[1][k];
									norvec[2][k] = -ATTSURF[2][k];
								}
							}
							////
							for (int k = 0; k < 4; k++)
							{
								double fviflt[3], firflt[3], fvi[3], fir[3], fvicyl[3], fircyl[3];
								facabdi = fabs(cosang[k]);
								facrefl = pow(cosang[k], 2);
								for (int j = 0; j < 3; j++)
								{
									tabsp[j] = facabdi * ncfvec[j];
									//printf("tabsp %.15f\n", tabsp[j]);
									tabsp2[j] = tabsp[j];
									trefl[j] = 2.0 * facrefl * norvec[j][k];
									//printf("trefl %.15f\n", trefl[j]);
									trefl2[j] = 4.0 / 3.0 * facrefl * norvec[j][k];
									//printf("trefl2 %.15f\n", trefl2[j]);
									tdifu[j] = facabdi * (ncfvec[j] + 2.0 / 3.0 * norvec[j][k]);
									//printf("tdifu %.15f\n", tdifu[j]);
									tdifu2[j] = facabdi * (ncfvec[j] + (pi / 6.0) * norvec[j][k]);
									//printf("tdifu2 %.15f\n\n", tdifu2[j]);
									if (k < 3)
									{
										int flag = side[k] - 1;
										//printf("%d\n", flag);
										prtflt[0][j] = tdifu[j] * (1 - refl[k][flag]);
										prtflt[1][j] = trefl[j] * refl[k][flag];
										prtflt[2][j] = 0.0;

										prtcyl[0][j] = tdifu2[j] * (1 - refl2[k][flag]);
										prtcyl[1][j] = trefl2[j] * refl2[k][flag];
										prtcyl[2][j] = 0.0;

										prtfltir[0][j] = tdifu[j] * (1 - isr);
										prtfltir[1][j] = trefl[j] * isr;
										prtfltir[2][j] = 0.0;

										prtcylir[0][j] = tdifu2[j] * (1 - isr);
										prtcylir[1][j] = trefl2[j] * isr;
										prtcylir[2][j] = 0.0;
										/*printf("%.15f %.15f %.15f\n", prtflt[0][j], prtflt[1][j], prtflt[2][j]);
										printf("%.15f %.15f %.15f\n\n", prtfltir[0][j], prtfltir[1][j], prtfltir[2][j]);*/
									}
									else  //else if(k==3)
									{
										int flag = side[k] - 1;
										prtflt[0][j] = tabsp[j] * (1 - refl[k][flag] - difu[k][flag]);
										prtflt[1][j] = trefl[j] * refl[k][flag];
										prtflt[2][j] = tdifu[j] * difu[k][flag];

										prtcyl[0][j] = tabsp2[j] * (1 - refl2[k][flag] - difu2[k][flag]);
										prtcyl[1][j] = trefl2[j] * refl2[k][flag];
										prtcyl[2][j] = tdifu2[j] * difu2[k][flag];

										prtfltir[0][j] = tabsp[j] * (1 - isr - idr);
										prtfltir[1][j] = trefl[j] * isr;
										prtfltir[2][j] = tdifu[j] * idr;

										prtcylir[0][j] = tabsp2[j] * (1 - isr - idr);
										prtcylir[1][j] = trefl2[j] * isr;
										prtcylir[2][j] = tdifu2[j] * idr;
										//printf("%.15f %.15f %.15f\n", prtflt[0][j], prtflt[1][j], prtflt[2][j]);
										//printf("%.15f %.15f %.15f\n\n", prtfltir[0][j], prtfltir[1][j], prtfltir[2][j]);
									}

									fviflt[j] = prtflt[0][j] + prtflt[1][j] + prtflt[2][j];
									firflt[j] = prtfltir[0][j] + prtfltir[1][j] + prtfltir[2][j];
									fvicyl[j] = prtcyl[0][j] + prtcyl[1][j] + prtcyl[2][j];
									fircyl[j] = prtcylir[0][j] + prtcylir[1][j] + prtcylir[2][j];

									fviflt[j] = fviflt[j] * area[k][1];
									firflt[j] = firflt[j] * area[k][1];
									fvicyl[j] = fvicyl[j] * area2[k][1];
									fircyl[j] = fircyl[j] * area2[k][1];
									fvi[j] = (fviflt[j] + fvicyl[j]) * absncfvi;
									fir[j] = (firflt[j] + fircyl[j]) * absncfir;
									force[j] = force[j] + fvi[j] + fir[j];
								}
							}
						}
					}
				}
			}
			/*double cosu = cos(psi) / cosb;
			double deltau;
			double at = dotproduct(radvel, Emid);
			if (at >= 0)
			{
				deltau = acos(cosu) * 180 / pi;
			}
			else
			{
				deltau = 360 - acos(cosu) * 180 / pi;
			}*/
		
			for (int k = 0; k < 3; k++)
			{
				forceACR[0] = forceACR[0] + force[k] * radvec[k];
				forceACR[1] = forceACR[1] + force[k] * alg_vec[k];
				forceACR[2] = forceACR[2] + force[k] * crs_vec[k];
			}
			forceACR = forceACR / _mass;

			_acc = force / _mass;
			
			/*printf("%.2f ", beita);*/
			//printf("%.10f %.10f %.10f \n", forceACR[0] * 1e9, forceACR[1] * 1e9, forceACR[2] * 1e9);
			/*printf("%f\n", deltau);*/
			return true;
		}
		catch (...)
		{
			write_log(_glog, 4, function_name, "ERROR : throw exception while calculating atmosphere acc");
			return false;
		}
	}
}