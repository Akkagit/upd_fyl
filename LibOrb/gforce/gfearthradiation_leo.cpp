/**
*
* @verbatim
History
-1.0 Qian Zhang     2019-11-10 Creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University.All rights reserved.
*
* @file		    gfearthradiation.cpp
* @brief		    The class for earth radiation model.
*
* @author          Qian Zhang, Wuhan University
* @version		    1.0.0
* @date		    2019-11-10
*
*/

#include "gforce/gfearthradiation_leo.h"
#include "gutils/ginfolog.h"

using namespace std;

namespace great
{
	/** constructor. */
	t_gfearthradiation_leo::t_gfearthradiation_leo(t_init_arg* arg, t_glog* log) : t_gforce(log)
	{
		string function_name = "t_gfearthradiation::t_gfearthradiation";
		if (!log)
		{
			write_log(log, 0, function_name,
				"ERROR : the log point is empty");
		}
		
		t_init_earthradiation *earad = NULL;
		earad = dynamic_cast<t_init_earthradiation*>(arg);

		if (!earad || !earad->galbedo || !earad->gpannel)
		{
			write_log(log, 4, function_name, "ERROR : init failed ");
			if (!earad)
			{
				write_log(log, 4, function_name, "ERROR : no init argument ");
			}
			// if no data in argument
			if (!earad->galbedo)
			{
				write_log(log, 4, function_name, "ERROR : no solar data ");
			}
			if (!earad->gpannel)
			{
				write_log(log, 4, function_name, "ERROR : no pannel data ");
			}
		}
		else
		{
			// set the init value
			_galbedo = earad->galbedo;
			_gpannel = earad->gpannel;

			_sat_mas = earad->mass;
			_model = earad->model;

			_v_ns_all1.resize(72, 144);
			_v_ns_all2.resize(72, 144);
			_v_ns_all3.resize(72, 144);
			_d_area_all.resize(72, 144);
			_cergem.resize(72, 144);
			_cergre.resize(72, 144);
			_toa = 6371000.0 + 30000.0;
			_pi = 4.0*atan(1.0);
		}
	}

	/**
	* defalut deconstructor.
	*/
	t_gfearthradiation_leo::~t_gfearthradiation_leo()
	{

	}

	/**value latkmx/lonkmx.
	* @param[in] grd the grid
	*/
	bool t_gfearthradiation_leo::cal_max(int grd)
	{
		string function_name = "t_gfearthradiation::calmax";
		try
		{
			if (grd == 1)
			{
				_grd_num = 1.0;
				_latkmx = 72;
				_lonkmx = 144;
			}
			else if (grd == 2)
			{
				_grd_num = 2.0;
				_grdcer = 2;
				_latkmx = 36;
				_lonkmx = 72;
			}
			else if (grd == 3)
			{
				_grd_num = 4.0;
				_grdcer = 4;
				_latkmx = 18;
				_lonkmx = 36;
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
	bool t_gfearthradiation_leo::cal_cerg_re_em(int grd, Matrix& ceresdata_re, Matrix& ceresdata_em)
	{
		double grdang = 2.5;
		double latind = 36.5;
		double lonind = 72.5;
		double d_ang, lat_in, lon_in, coslat;
		int latk, lonk;
		
		string function_name = "t_gfearthradiation::calcerg_re_em";
		try
		{
			grdang = grdang*_grd_num;
			latind = (latind - 0.5) / _grd_num + 0.5;
			lonind = (lonind - 0.5) / _grd_num + 0.5;
			d_ang = (_pi*grdang / 180.0)*(_pi*grdang / 180.0);
			for (latk = 1; latk <= _latkmx; latk++)
			{
				for (lonk = 1; lonk <= _lonkmx; lonk++)
				{
					lat_in = (latk - latind)*grdang*(_pi / 180.0);
					lon_in = (lonk - lonind)*grdang*(_pi / 180.0);
					coslat = cos(lat_in);
					_d_area_all(latk, lonk) = _toa*_toa*coslat*d_ang;
					_v_ns_all1(latk, lonk) = coslat*cos(lon_in);
					_v_ns_all2(latk, lonk) = coslat*sin(lon_in);
					_v_ns_all3(latk, lonk) = sin(lat_in);
					_cergem(latk, lonk) = 0.0;
					_cergre(latk, lonk) = 0.0;
					if (grd == 1)
					{
						_cergre(latk, lonk) = ceresdata_re(latk, lonk);
						_cergem(latk, lonk) = ceresdata_em(latk, lonk);
					}
					else if (grd == 2 || grd == 3)
					{
						for (int ii = 0; ii < _grdcer; ii++)
						{
							for (int jj = 0; jj < _grdcer; jj++)
							{
								_cergre(latk, lonk) = _cergre(latk, lonk) + ceresdata_re(_grdcer*latk - ii, _grdcer*lonk - jj);
								_cergem(latk, lonk) = _cergem(latk, lonk) + ceresdata_em(_grdcer*latk - ii, _grdcer*lonk - jj);
							}
						}
						_cergre(latk, lonk) = _cergre(latk, lonk) / (_grd_num*_grd_num);
						_cergem(latk, lonk) = _cergem(latk, lonk) / (_grd_num*_grd_num);
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
	bool t_gfearthradiation_leo::calculate(t_gforce_calarg* argument)
	{
		string function_name = "t_gfearthradiation::calcearthradmod";
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
		    double abspos, abssun, psidot, psi, albfac, phasevi, phaseir, absncfvi, absncfir;
		    double eau = 149597870691.0;
		    double s0 = 1367.0;
		    double alb = 0.3;
		    double c = 299792458.0;
			double cosa = 0.0;
			int erm, grd, reff;
			vector<double> esun;
			Matrix kappa;
			ColumnVector radvec, ncfvec;
			ColumnVector forcevis(0.0, 3);
			ColumnVector forceir(0.0, 3);
			ColumnVector force(0.0, 3);
			ColumnVector xtmp(0.0, 3);
			sat_crs = *(arg->sat_pos);
			sun_crs = *(arg->sun_pos_crs);
			erm = 2;
			grd = 1;
			reff = 0;
			kappa = (*arg->e2j).t();
		    abspos = sqrt(sat_crs[0] * sat_crs[0] + sat_crs[1] * sat_crs[1] + sat_crs[2] * sat_crs[2])*1000.0;
			radvec = sat_crs*1000.0 / abspos;		   
		    if (erm>0)
		    {
			    abssun = sqrt(sun_crs[0] * sun_crs[0] + sun_crs[1] * sun_crs[1] + sun_crs[2] * sun_crs[2])*1000.0;
			    for (int k = 0; k < 3; k++)
			    {
				    esun.push_back(sun_crs[k]*1000.0 / abssun);
			    }
			    psidot = esun[0] * radvec[0] + esun[1] * radvec[1] + esun[2] * radvec[2];
				if (fabs(psidot)>1 - 1.0e-6)
			    {
				    psi = 0.0;
			    }
			    else
			    {
				    psi = acos(psidot);
			    }
			    s0 = s0*(eau / abssun)*(eau / abssun);
		    }
		    if (erm == 1)
		    {
				ncfvec = radvec;
			    albfac = (_pi*_toa*_toa)*(s0 / c) / (abspos*abspos);
			    phasevi = (2 * alb / (3 * _pi*_pi))*((_pi - psi)*cos(psi) + sin(psi));
			    phaseir = (1 - alb) / (4 * _pi);
			    absncfvi = albfac*phasevi;
			    absncfir = albfac*phaseir;				
				
				map<int, t_oneside>::iterator it_pannel = _gpannel->get_pannels()->begin();
				for (it_pannel; it_pannel != _gpannel->get_pannels()->end(); it_pannel++)
				{
					xtmp = it_pannel->second.crs_norm_dir;
					cosa = dotproduct(xtmp, xtmp);
					if (cosa < 1.0e-10)
					{
						xtmp = ncfvec;
					}
					cosa = dotproduct(ncfvec, xtmp);
					if (cosa <= 0)
					{
						continue;
					}
					for (int j = 1; j <= 3; j++)
					{
						forcevis(j) += it_pannel->second.area*cosa*(2.0*(it_pannel->second.dr_coef / 3.0 + it_pannel->second.sr_coef*cosa)*xtmp(j) + (1.0 - it_pannel->second.sr_coef)*ncfvec(j));
						forceir(j) += it_pannel->second.area*cosa*(2.0*(it_pannel->second.idr_coef / 3.0 + it_pannel->second.isr_coef*cosa)*xtmp(j) + (1.0 - it_pannel->second.isr_coef)*ncfvec(j));
					}
				}
				forcevis = forcevis*absncfvi;
				forceir = forceir*absncfir;
				force += forcevis + forceir;
		    }
			else if (erm == 2)
			{
				double d_area, dist2, absdist, cos_in, cos_re;
				double refl_cf, emit_cf, e_refl, e_emit;
				ColumnVector v_ns(0.0, 3);
				ColumnVector v_ins(0.0, 3);
				ColumnVector v_dist(0.0, 3);
				ColumnVector v_sat(0.0, 3);
				abssun = sqrt(sun_crs[0] * sun_crs[0] + sun_crs[1] * sun_crs[1] + sun_crs[2] * sun_crs[2])*1000.0;
				for (int k = 0; k < 3; k++)
				{
					esun.push_back(sun_crs[k]*1000.0 / abssun);
				}	
				if (!cal_max(grd))
				{
					write_log(_glog, 4, function_name, "ERROR : in t_gfearthradiation::calmax");
				}
				static bool is_ref_emi;
				static bool lfirst = true;
				if (lfirst)
				{
					is_ref_emi = cal_cerg_re_em(grd, _galbedo->get_ceresr_data(), _galbedo->get_cerese_data());
					lfirst = false;
				}			
				for (int latk = 1; latk <= _latkmx; latk++)
				{
					for (int lonk = 1; lonk <= _lonkmx; lonk++)
					{
						d_area =  _d_area_all(latk, lonk);
						v_ns(1) = _v_ns_all1(latk, lonk);
						v_ns(2) = _v_ns_all2(latk, lonk);
						v_ns(3) = _v_ns_all3(latk, lonk);
						for (int ii = 1; ii <= 3; ii++)
						{
							v_ins(ii) = 0.0;
							for (int jj = 1; jj <= 3; jj++)
							{
								v_ins(ii) = v_ins(ii) + kappa(jj, ii)*v_ns(jj);
							}
						}
						//v_ins += kappa*v_ns;
						v_dist = sat_crs*1000.0 - _toa*v_ins;
						dist2 = v_dist(1)*v_dist(1) + v_dist(2)*v_dist(2) + v_dist(3)*v_dist(3);
						absdist = sqrt(dist2);
						v_sat = v_dist / absdist;
						cos_in = esun[0] * v_ins(1) + esun[1] * v_ins(2) + esun[2] * v_ins(3);
						cos_re = v_sat(1)*v_ins(1) + v_sat(2)*v_ins(2) + v_sat(3)*v_ins(3);
						if (cos_re >= 0)
						{
							if (is_ref_emi) ////????????
							{
								refl_cf = _cergre(latk, lonk);
								emit_cf = _cergem(latk, lonk);
							}
							   
							if (cos_in >= 0)
							{
								e_refl = (refl_cf / (_pi*dist2))*cos_re*cos_in*s0*d_area;
							}
							else
							{
								e_refl = 0.0;
							}
							e_emit = (emit_cf / (4 * _pi*dist2))*cos_re*s0*d_area;
							absncfvi = e_refl / c;
							absncfir = e_emit / c;
							ncfvec = v_sat;
							forcevis(1) = forcevis(2) = forcevis(3) = 0.0;
							forceir(1) = forceir(2) = forceir(3) = 0.0;
							map<int, t_oneside>::iterator it_pannel = _gpannel->get_pannels()->begin();
							for (it_pannel; it_pannel != _gpannel->get_pannels()->end(); it_pannel++)
							{
								xtmp = it_pannel->second.crs_norm_dir;
								cosa = dotproduct(xtmp, xtmp);
								if (cosa < 1.0e-10)
								{
									xtmp = ncfvec;
								}
								cosa = dotproduct(ncfvec, xtmp);
								if (cosa >= 0)
								{
									continue;
								}
								for (int j = 1; j <= 3; j++)
								{
									forcevis(j) += it_pannel->second.area*cosa*(2.0*(it_pannel->second.dr_coef / 3.0 + it_pannel->second.sr_coef*cosa)*xtmp(j) + (1.0 - it_pannel->second.sr_coef)*ncfvec(j));
									forceir(j) += it_pannel->second.area*cosa*(2.0*(it_pannel->second.idr_coef / 3.0 + it_pannel->second.isr_coef*cosa)*xtmp(j) + (1.0 - it_pannel->second.isr_coef)*ncfvec(j));
								}
							}
							forcevis = forcevis*absncfvi;
							forceir = forceir*absncfir;
							force += forcevis + forceir;
						}
					}
				}
			}
			_acc = force / _sat_mas;
			return true;
		}
		catch (...)
		{
			write_log(_glog, 4, function_name, "ERROR : throw exception while calculating atmosphere acc");
			return false;
		}
	}
}