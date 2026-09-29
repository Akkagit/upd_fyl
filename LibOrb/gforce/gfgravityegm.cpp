/**
 *
 * @verbatim
	 History
	  -1.0 Jie Li&Yujie Qin  2019-02-21 creat the file.
	  -1.1 Jiande Huang       2019-03-28 modfied t_gfgravity_egm08::_cal_pnm
	  -1.2 Jiande Huang       2019-04-01 Adding exception throwing mechanism for all functions
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		    gfgravityegm.cpp
 * @brief		    The class for EGM08 gravity model, including 12th-order EGM08, 120th-order EGM08 and so on.\n
*
 * @author          Jie Li, Wuhan University
 * @version		    1.0.0
 * @date		    2019-02-21
 *
 */

#include "gutils/gconst.h"
#include "gutils/gconst_oi.h"
#include "gutils/ginfolog.h"
#include "gforce/gfgravityegm.h"
#include "gforce/gfarg_cal.h"

using namespace std;

namespace great
{
	t_gfgravity_egm08::t_gfgravity_egm08(t_init_arg* arg, t_glog* log) : t_gforce(log)
	{
		string function_name = "t_gfgravity_egm08_12::t_gfgravity_egm08_12";

		t_init_gravity_egm08* egm08 = nullptr;
		egm08 = dynamic_cast<t_init_gravity_egm08*> (arg);
		if (!egm08 || !egm08->gracoef || !egm08->gotl || !egm08->gotdata || !egm08->gplanet)
		{
			write_log(log, 0, "ERROR : " + function_name + " ==> init failed");
			if (!egm08)
			{
				write_log(log, 0, "ERROR : " + function_name + " ==> no init argument");
			}

			// put the EGM data
			if (!egm08->gracoef)
			{
				write_log(log, 0, "ERROR : " + function_name + " ==> no EGM data, check your EGM file");
			}

			// put the oceanload data
			if (!egm08->gotl)
			{
				write_log(log, 0, "ERROR : " + function_name + " ==> no ocean load data, check your oceanload file");
			}

			// put the oceantide data
			if (!egm08->gotdata)
			{
				write_log_info(log, 0, "ERROR : " + function_name + " ==> no ocean ocean tide data, check your oceantide file");
			}

			// put the de nav data
			if (!egm08->gplanet)
			{
				write_log_info(log, 0, "ERROR : " + function_name + " ==> no ocean de nav data, check your de file");
			}

			// put the order of model
			if (egm08->gravity_order == 0)
			{
				write_log_info(log, 0, "ERROR : " + function_name + " ==> no degree order, check your xml file");
			}
		}
		else
		{
			_force_model = egm08->force_model;

			_gracoef = egm08->gracoef;
			_gotl = egm08->gotl;
			_gotdata = egm08->gotdata;
			_gplanets = egm08->gplanet;
			_gravi_degree = egm08->gravity_order;
			_ocean_degree = egm08->ocean_degree;
			_poletide = egm08->ispoletide;
			_ocean_tide_Cnm.resize(_max_ocean_degree, _max_ocean_degree + 1);
			_ocean_tide_Snm.resize(_max_ocean_degree, _max_ocean_degree + 1);
			_solid_tide_Cnm.resize(3, 4);
			_solid_tide_Snm.resize(3, 4);
			if (_poletide)
			{
				_pole_degree = egm08->pole_degree;
				_gdeadata = egm08->dedata;
				_load = _gdeadata->load();
				_AR = _gdeadata->Anm();
				_BR = _gdeadata->Bnm();
				_AI = _gdeadata->AInm();
				_BI = _gdeadata->BInm();
				_solid_tide_Cnm.resize(4, 4);
				_solid_tide_Snm.resize(4, 4);
				_ocean_pole_tide_Cnm.resize(30, 31);
				_ocean_pole_tide_Snm.resize(30, 31);
			}

			// force model for solide tide calculation(list: sun moom...)
			_planet_list.clear();
			_planet_list.insert("SUN");
			_planet_list.insert("MOON");

			// init the Matrix list
			_Cnm_save.resize(_gravi_degree + 1, _gravi_degree + 1);
			_Snm_save.resize(_gravi_degree + 1, _gravi_degree + 1);

			_Anm.resize(_gravi_degree + 3, _gravi_degree + 3);
			_SUM_F.resize(4, 1);

			// get the init Cnm and Snm
			if (!_isGetCnm())
			{
				write_log_info(_glog, 0, "ERROR : " + function_name + " ==> we set the Cnm zero");
			};

			if (!_isGetSnm())
			{
				write_log_info(_glog, 0, "ERROR : " + function_name + " ==> we set the Snm zero");
			};
		}
	}

	/** show the detail of the gravity model. */
	void t_gfgravity_egm08::show()
	{
		cout << "gravity model : EGM " << endl;
	}

	/** check the important parameters and information in the model */
	bool t_gfgravity_egm08::check()
	{
		cout << "check some important pars in model : EGM " << endl;
		return true;
	}


	/** The acceleration and partial matrix will be calculated here.
	* @param [in] argument Parameters to be used in calculation.
	*/
	bool t_gfgravity_egm08::calculate(t_gforce_calarg* argument)
	{
		string function_name = "t_gfgravity_egm08::calculate";

		// 判断计算的充分条件
		t_gfcalarg_gravity_egm* arg = NULL;
		arg = dynamic_cast<t_gfcalarg_gravity_egm*>(argument);

		if (!arg)
		{
			write_log_info(_glog, 0, "ERROR : " + function_name + " ==> no gravity_egm argument");
			return false;
		}

		if (arg->arg_type != ARG_GRAVI_EGM)
		{
			write_log_info(_glog, 0, "ERROR : " + function_name + " ==> error gravity_egm argument");
			return false;
		}

		try
		{
			bool isValid;
			// get the moon and sun pos
			_planet_crd_crs = *(arg->allXE);

			//get solitide and oceantide Cnm Snm correct for gravity
			isValid = _calculate_solid_tide(*(arg->time), *(arg->gst));
			if (!isValid)
			{
				write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::calculate ==> can not calculate solid tide!");
				return false;
			};

			isValid = _calculate_ocean_tide(*(arg->time), *(arg->gst));
			if (!isValid)
			{
				write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::calculate ==> can not calculate ocean tide!");
				return false;
			};

			// get the distance and normalize sat position
			_distance = arg->sat_pos_trs->NormFrobenius();
			_sat_pos_norm = *(arg->sat_pos_trs) / _distance;

			// get the coef from solid tide and ocean tide
			double mjd = *(arg->time);

			// c20 rate corrections
			// constant values
			// 计算初始化系数
			if (_isfirst)
			{
				_isfirst = false;
				// mjd ref is got from EGM file
				double mjd_dif = mjd - _gracoef->mjd_ref();

				const double _c20rate = 1.162755e-11;
				const double _c30rate = 0.490000e-11;
				const double _c40rate = 0.470000e-11;
				const double _c21rate = -0.33700e-11;
				const double _s21rate = 1.60600e-11;

				_Cnm_save(3, 1) = _Cnm_save(3, 1) + _c20rate * mjd_dif / 365.25;
				_Cnm_save(4, 1) = _Cnm_save(4, 1) + _c30rate * mjd_dif / 365.25;
				_Cnm_save(5, 1) = _Cnm_save(5, 1) + _c40rate * mjd_dif / 365.25;
				_Cnm_save(3, 2) = _Cnm_save(3, 2) + _c21rate * (mjd - 51545) / 365.25;
				_Snm_save(3, 2) = _Snm_save(3, 2) + _s21rate * (mjd - 51545) / 365.25;
			}

			// 保存初始化系数。注意，每次计算都是使用初始化系数重新开始计算的
			_Cnm << _Cnm_save;
			_Snm << _Snm_save;

			//mean pole according to IERS 2010
			double xpm = 0.0;
			double ypm = 0.0;
			_get_mean_pole(mjd, xpm, ypm);
			double m1 = *(arg->xpole) - xpm;
			double m2 = -(*(arg->ypole) - ypm);

			//Solid earth pole tide
			_Cnm(3, 2) = _Cnm(3, 2) - 1.333e-9 * (m1 - 0.0115 * m2);
			_Snm(3, 2) = _Snm(3, 2) - 1.333e-9 * (m2 + 0.0115 * m1);

			//ocean pole tide 
			
			if (_poletide)
			{
				isValid = _calculate_ocean_pole_tide(m1, m2, 30);
				if (!isValid)
				{
					write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::calculate ==> can not calculate ocean pole tide!");
					return false;
				};
				for (int n = 1; n <= 30; n++)
				{
					for (int m = 1; m <= n; m++)
					{
						_Cnm(n + 1, m + 1) = _Cnm(n + 1, m + 1) + _ocean_pole_tide_Cnm(n, m + 1);
						_Snm(n + 1, m + 1) = _Snm(n + 1, m + 1) + _ocean_pole_tide_Snm(n, m + 1);
					}
				}
			}
			else
			{
				_Cnm(3, 2) = _Cnm(3, 2) - 2.1778e-10 * (m1 - 0.01724 * m2);
				_Snm(3, 2) = _Snm(3, 2) - 1.7232e-10 * (m2 - 0.03365 * m1);
			}
			

			//Cnm Snm correction due to solid earth tide and/or ocean tide 
			int max_degree = max(3, _ocean_degree);
			for (int n = 2; n < max_degree + 1; n++)
			{
				for (int m = 0; m < max_degree + 1; m++)
				{
					if (n > _gravi_degree || m > _gravi_degree)
					{
						break;
					}

					_Cnm(n + 1, m + 1) = _Cnm(n + 1, m + 1) + _ocean_tide_Cnm(n, m + 1);
					_Snm(n + 1, m + 1) = _Snm(n + 1, m + 1) + _ocean_tide_Snm(n, m + 1);
					
					if (_poletide)
					{
						if (n <= 4 && m <= 3)
						{
							_Cnm(n + 1, m + 1) = _Cnm(n + 1, m + 1) + _solid_tide_Cnm(n, m + 1);
							_Snm(n + 1, m + 1) = _Snm(n + 1, m + 1) + _solid_tide_Snm(n, m + 1);
						}
					}
					else 
					{
						if (n <= 3 && m <= 3)
						{
							_Cnm(n + 1, m + 1) = _Cnm(n + 1, m + 1) + _solid_tide_Cnm(n, m + 1);
							_Snm(n + 1, m + 1) = _Snm(n + 1, m + 1) + _solid_tide_Snm(n, m + 1);
						}
					}
				}
			}

			isValid = _cal_Anm_coeff();
			if (!isValid)
			{
				write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::calculate ==> can not calculate Anm!");
				return false;
			}

			isValid = _cal_acc_mat(*(arg->e2j));
			if (!isValid)
			{
				write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::calculate ==> can not calculate Acc!");
				return false;
			}

			isValid = _cal_Amat_coeff(*(arg->e2j));
			if (!isValid)
			{
				write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::calculate ==> can not calculate Amat!");
				return false;
			}

			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 0, "t_gfgravity_egm08::calculate", "ERROR : throw exception");
			return false;
		}
	}

	/** Calculating and Evaluating Pole Shift Point Location.Note that the default coordinate system is TRS
	* @param [in] mjd Dynamic time.
	* @param [out] xpm The x pos of pole point.
	* @param [out] ypm The y pos of pole point.
	*/
	void t_gfgravity_egm08::_get_mean_pole(const double& mjd, double& xpm, double& ypm)
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


	bool t_gfgravity_egm08::_isGetCnm()
	{
		try
		{
			_Cnm_save = 0.0;
			//matrix start from 1
			for (int m = 0; m <= _gravi_degree; m++)
			{
				for (int n = m; n <= _gravi_degree; n++)
				{
					_Cnm_save(n + 1, m + 1) = _gracoef->get_c(n, m);
				}
			}

			//unit
			if (fabs(_Cnm_save(3, 1)) > 1.0)
			{
				for (int m = 0; m <= _gravi_degree; m++)
				{
					for (int n = m; n <= _gravi_degree; n++)
					{
						_Cnm_save(n + 1, m + 1) = _Cnm_save(n + 1, m + 1) * 1.0e-6;
					}
				}
			}
			return true;
		}
		catch (...)
		{
			_Cnm_save = 0.0;
			write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::_get_Cnm ==> throw exception");
			write_log_info(_glog, 0, "NOTE  : t_gfgravity_egm08::_get_Cnm ==> We reset the Cnm to zeros for safety");
			return false;
		}
	}


	bool t_gfgravity_egm08::_isGetSnm()
	{
		_Snm_save = 0.0;
		try
		{
			for (int m = 0; m <= _gravi_degree; m++)
			{
				for (int n = m; n <= _gravi_degree; n++)
				{
					_Snm_save(n + 1, m + 1) = _gracoef->get_s(n, m);
				}
			}

			//unit
			if (fabs(_Cnm_save(3, 1)) > 1.0)
			{
				for (int m = 0; m <= _gravi_degree; m++)
				{
					for (int n = m; n <= _gravi_degree; n++)
					{
						_Snm_save(n + 1, m + 1) = _Snm_save(n + 1, m + 1) * 1.0e-6;
					}
				}
			}
			return true;
		}
		catch (...)
		{
			_Snm_save = 0.0;
			write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::_get_Cnm ==> throw exception");
			write_log_info(_glog, 0, "NOTE  : t_gfgravity_egm08::_get_Cnm ==> We reset the Snm to zeros for safety");
			return false;
		}

	}

	/** caculate the Anm coefficents by using the EGM08 coefficents read from file
	* @remarks Specific formulas can be referred to Geodetic books.
	*/
	bool t_gfgravity_egm08::_cal_Anm_coeff()
	{
		try
		{
			_Anm = 0.0;
			_Anm(1, 1) = 1.0;
			_Anm(2, 2) = sqrt(3.0);

			for (int n = 2; n < _gravi_degree + 3; n++)
			{
				_Anm(n + 1, n + 1) = _Anm(n, n) * sqrt(1.0 + 0.5 / n);
			}

			for (int j = 1; j < _gravi_degree + 3; j++)
			{
				double rat1 = sqrt((1.0 * ((int64_t)j + j + 1) * ((int64_t)j + 1)) / (1.0 * j * ((int64_t)4 * j + 2)));
				double rat2 = sqrt((1.0 * ((int64_t)j + j + 1) * ((int64_t)j - 1)) / (1.0 * j * ((int64_t)4 * j - 2)));
				_Anm(j + 1, 1) = rat1 * _sat_pos_norm(3) * _Anm(j + 1, 2) - rat2 * _Anm(j, 2);

				if ((j + 1) < (_gravi_degree + 3))
				{
					rat1 = sqrt((1.0 * ((int64_t)j + j + 3) * j) / (1.0 * ((int64_t)j + j + 1) * ((int64_t)j + 2)));
					rat2 = sqrt((1.0 * ((int64_t)4 * j + 6) * ((int64_t)j + 1)) / (1.0 * ((int64_t)j + j + 1) * ((int64_t)j + 2)));
					_Anm(j + 2, 2) = rat2 * _Anm(j + 1, 1) + rat1 * _sat_pos_norm(3) * _Anm(j + 1, 2);
				}

				for (int n = j + 2; n < _gravi_degree + 3; n++)
				{
					int m = n - j;
					rat1 = sqrt((1.0 * ((int64_t)n + n + 1) * ((int64_t)n - m)) / (1.0 * ((int64_t)n + n - 1) * ((int64_t)n + m)));
					rat2 = rat1 * sqrt((1.0 * ((int64_t)n + m - 1)) / (1.0 * ((int64_t)n - m)));
					_Anm(n + 1, m + 1) = rat2 * _Anm(n, m) + rat1 * _sat_pos_norm(3) * _Anm(n, m + 1);
				}
			}

			return true;
		}
		catch (...)
		{
			_Anm = 0.0;
			write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::_calculate_Anm ==> throw exception");
			write_log_info(_glog, 0, "NOTE  : t_gfgravity_egm08::_calculate_Anm ==> We reset the Anm to zeros for safety");
			return false;
		}
	}


	/** caculate gravity acceleration in J2000 by using the coefficients.
	* @remarks Specific formulas can be referred to Geodetic books.And the solid tide and
	* ocean tide were considered.
	* @param [in] ej2 The rotation matrix from TRS to J2000.
	*/
	bool t_gfgravity_egm08::_cal_acc_mat(const Matrix& e2j)
	{
		try
		{
			vector<double> rm, im, ro;
			for (int i = 0; i < _gravi_degree + 3; i++)
			{
				if (i < 3)
				{
					rm.push_back(1.0);
					im.push_back(0.0);
					continue;
				}
				double tmp = _sat_pos_norm(1) * rm[(int64_t)i - 1] - _sat_pos_norm(2) * im[(int64_t)i - 1];
				rm.push_back(tmp);

				tmp = _sat_pos_norm(1) * im[(int64_t)i - 1] + _sat_pos_norm(2) * rm[(int64_t)i - 1];
				im.push_back(tmp);
			}

			double gm = _gracoef->GM();
			double radius = _gracoef->radius();

			double dummy = radius / _distance;
			ro.push_back(gm / _distance / radius);
			for (int i = 1; i < _gravi_degree + 3; i++)
			{
				ro.push_back(dummy * ro[(int64_t)i - 1]);
			}

			_SUM_F = 0.0;
			double SUM_Fn[4] = { 0.0 };
			double ratio01 = 0.0;
			double ratio11 = 0.0;
			double Dnm = 0.0;
			double Enm = 0.0;
			double Fnm = 0.0;

			// loop over n
			for (int n = 2; n <= _gravi_degree; n++)
			{
				for (int i = 0; i < 4; i++)
				{
					SUM_Fn[i] = 0.0;
				}

				for (int m = 0; m <= n; m++)
				{
					//(de)normalization
					ratio01 = sqrt(1.0 * (((int64_t)n + m + 1) * ((int64_t)n - m)));
					ratio11 = sqrt(1.0 * (((int64_t)2 * n + 1) * ((int64_t)n + m + 1) * ((int64_t)n + m + 2)) / ((int64_t)2 * n + 3));

					if (m == 0)
					{
						ratio01 = ratio01 * sqrt(0.5);
						ratio11 = ratio11 * sqrt(0.5);
					}
					Dnm = _Cnm(n + 1, m + 1) * rm[(int64_t)m + 2] + _Snm(n + 1, m + 1) * im[(int64_t)m + 2];
					Enm = _Cnm(n + 1, m + 1) * rm[(int64_t)m + 1] + _Snm(n + 1, m + 1) * im[(int64_t)m + 1];
					Fnm = _Snm(n + 1, m + 1) * rm[(int64_t)m + 1] - _Cnm(n + 1, m + 1) * im[(int64_t)m + 1];

					SUM_Fn[0] = SUM_Fn[0] + _Anm(n + 1, m + 1) * m * Enm;
					SUM_Fn[1] = SUM_Fn[1] + _Anm(n + 1, m + 1) * m * Fnm;
					SUM_Fn[2] = SUM_Fn[2] + _Anm(n + 1, m + 2) * Dnm * ratio01;
					SUM_Fn[3] = SUM_Fn[3] - _Anm(n + 2, m + 2) * Dnm * ratio11;
				}

				for (int i = 0; i < 4; i++)
				{
					_SUM_F(i + 1, 1) += ro[(int64_t)n + 1] * SUM_Fn[i];
				}

			}

			_SUM_F(1, 1) += _SUM_F(4, 1) * _sat_pos_norm(1);
			_SUM_F(2, 1) += _SUM_F(4, 1) * _sat_pos_norm(2);
			_SUM_F(3, 1) += _SUM_F(4, 1) * _sat_pos_norm(3);

			//acceleration to inertial system
			_acc = e2j * _SUM_F.submatrix(1, 3, 1, 1);
			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::_calculate_Acc ==> throw exception");
			return false;
		}
	}

	/** caculate gravity acceleration partial matrix in J2000 by using the coefficients.
	* @remarks Specific formulas can be referred to Geodetic books.And the solid tide and
	* ocean tide were considered.
	* @param [in] ej2 The rotation matrix from TRS to J2000.
	*/
	bool t_gfgravity_egm08::_cal_Amat_coeff(const Matrix& e2j)
	{
		try
		{
			Matrix SUM_P(4, 4), SUM_Pn(4, 4);
			SUM_P = 0.0;
			double tmp = 0.0;
			vector<double> rm, im, ro;

			double gm = _gracoef->GM();
			double radius = _gracoef->radius();
			double dummy = radius / _distance;
			ro.push_back(gm / _distance / radius);

			for (int i = 0; i < _gravi_degree + 3; i++)
			{
				if (i < 3)
				{
					rm.push_back(1.0);
					im.push_back(0.0);
					continue;
				}
				tmp = _sat_pos_norm(1) * rm[(int64_t)i - 1] - _sat_pos_norm(2) * im[(int64_t)i - 1];
				rm.push_back(tmp);

				tmp = _sat_pos_norm(1) * im[(int64_t)i - 1] + _sat_pos_norm(2) * rm[(int64_t)i - 1];
				im.push_back(tmp);
			}

			for (int i = 1; i < _gravi_degree + 3; i++)
			{
				ro.push_back(dummy * ro[(int64_t)i - 1]);
			}

			for (int n = 2; n <= _gravi_degree; n++)
			{
				SUM_Pn = 0.0;
				for (int m = 0; m <= n; m++)
				{
					double ratio01 = sqrt(1.0 * (((int64_t)n + m + 1) * ((int64_t)n - m)));
					double ratio11 = sqrt(1.0 * (((int64_t)2 * n + 1) * ((int64_t)n + m + 1) * ((int64_t)n + m + 2)) / ((int64_t)2 * n + 3));

					if (m == 0)
					{
						ratio01 = ratio01 * sqrt(0.5);
						ratio11 = ratio11 * sqrt(0.5);
					}

					double ratio12 = ratio11 * sqrt(1.0 * (((int64_t)n - m) * ((int64_t)n + m + 3)));
					double ratio22 = ratio11 * sqrt(1.0 * (((int64_t)2 * n + 3) * ((int64_t)n + m + 3) * ((int64_t)n + m + 4)) / ((int64_t)2 * n + 5));

					double Dnm = _Cnm(n + 1, m + 1) * rm[(int64_t)m + 2] + _Snm(n + 1, m + 1) * im[(int64_t)m + 2];
					double Enm = _Cnm(n + 1, m + 1) * rm[(int64_t)m + 1] + _Snm(n + 1, m + 1) * im[(int64_t)m + 1];
					double Fnm = _Snm(n + 1, m + 1) * rm[(int64_t)m + 1] - _Cnm(n + 1, m + 1) * im[(int64_t)m + 1];
					double Gnm = _Cnm(n + 1, m + 1) * rm[m] + _Snm(n + 1, m + 1) * im[m];
					double Hnm = _Snm(n + 1, m + 1) * rm[m] - _Cnm(n + 1, m + 1) * im[m];

					SUM_Pn(1, 1) = SUM_Pn(1, 1) + _Anm(n + 1, m + 1) * m * (m - 1.0) * Gnm;
					SUM_Pn(1, 2) = SUM_Pn(1, 2) + _Anm(n + 1, m + 1) * m * (m - 1.0) * Hnm;
					SUM_Pn(1, 3) = SUM_Pn(1, 3) + _Anm(n + 1, m + 2) * m * Enm * ratio01;
					SUM_Pn(1, 4) = SUM_Pn(1, 4) - _Anm(n + 2, m + 2) * m * Enm * ratio11;
					SUM_Pn(2, 3) = SUM_Pn(2, 3) + _Anm(n + 1, m + 2) * m * Fnm * ratio01;
					SUM_Pn(2, 4) = SUM_Pn(2, 4) - _Anm(n + 2, m + 2) * m * Fnm * ratio11;
					SUM_Pn(3, 4) = SUM_Pn(3, 4) - _Anm(n + 2, m + 3) * Dnm * ratio12;
					SUM_Pn(4, 4) = SUM_Pn(4, 4) + _Anm(n + 3, m + 3) * Dnm * ratio22;

				}
				for (int i = 1; i < 5; i++)
				{
					for (int j = 1; j < 5; j++)
					{
						SUM_P(i, j) = SUM_P(i, j) + ro[(int64_t)n + 2] / radius * SUM_Pn(i, j);
					}
				}
			}

			SUM_Pn(1, 1) = SUM_P(1, 1) + _sat_pos_norm(1) * _sat_pos_norm(1) * SUM_P(4, 4) + _SUM_F(4, 1) / _distance + 2.0 * _sat_pos_norm(1) * SUM_P(1, 4);
			SUM_Pn(1, 2) = SUM_P(1, 2) + _sat_pos_norm(1) * _sat_pos_norm(2) * SUM_P(4, 4) + _sat_pos_norm(1) * SUM_P(2, 4) + _sat_pos_norm(2) * SUM_P(1, 4);
			SUM_Pn(1, 3) = SUM_P(1, 3) + _sat_pos_norm(1) * _sat_pos_norm(3) * SUM_P(4, 4) + _sat_pos_norm(1) * SUM_P(3, 4) + _sat_pos_norm(3) * SUM_P(1, 4);
			SUM_Pn(2, 2) = -SUM_P(1, 1) + _sat_pos_norm(2) * _sat_pos_norm(2) * SUM_P(4, 4) + _SUM_F(4, 1) / _distance + 2.0 * _sat_pos_norm(2) * SUM_P(2, 4);
			SUM_Pn(2, 3) = SUM_P(2, 3) + _sat_pos_norm(2) * _sat_pos_norm(3) * SUM_P(4, 4) + _sat_pos_norm(3) * SUM_P(2, 4) + _sat_pos_norm(2) * SUM_P(3, 4);

			SUM_Pn(3, 3) = -(SUM_Pn(1, 1) + SUM_Pn(2, 2));
			SUM_Pn(2, 1) = SUM_Pn(1, 2);
			SUM_Pn(3, 1) = SUM_Pn(1, 3);
			SUM_Pn(3, 2) = SUM_Pn(2, 3);
			SUM_Pn.resize_keep(3, 3);

			_partial = e2j * SUM_Pn * e2j.t();
			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::_calculate_Amat ==> throw exception.");
			return false;
		}
	}


	/** caculate solid tide coefficients.
	* @remarks Specific formulas can be referred to Geodetic books.
	* @param [in] tdt Dynamic time.
	* @param [in] gst Greenwich sidereal time.
	*/
	bool t_gfgravity_egm08::_calculate_solid_tide(const double& tdt, const double& gst)
	{
		double theta = 0.0;
		double sin_theta = 0.0;
		double cos_theta = 0.0;
		Matrix CF(0.0, 3, 4);
		Matrix SF(0.0, 3, 4);
		try
		{
		double love_2p[3] = { -0.00087, -0.00079, -0.00057 };
		if(_poletide)
		{
			CF.resize(4, 4);
			SF.resize(4, 4);
		}
		
			

			_solid_tide_Cnm = 0.0;
			_solid_tide_Snm = 0.0;

			Matrix norm;
			if (!_normalize_coeff("norm", 3, 3, 3, 3, norm))
			{
				write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::_calculate_solid_tide ==> can not get normalize coeff");
				return false;
			};

			for (int n = 2; n < 4; n++)
			{
				for (int m = 0; m <= 3; m++)
				{
					norm(n, m + 1) = norm(n, m + 1) / ((int64_t)2 * n + 1.0);
				}
			}

			for (auto i = _planet_list.begin(); i != _planet_list.end(); i++)
			{
				Matrix pnm;
				if (!_cal_pnm(*i, tdt, pnm))
				{
					write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::_calculate_solid_tide ==> can not calculate the Pnm");
					return false;
				};

				for (int n = 2; n < 4; n++)
				{
					for (int m = 0; m < n + 1; m++)
					{
						_solid_tide_Cnm(n, m + 1) += norm(n, m + 1) * _factor * pnm(n, m + 1) * cos(m * _lon) * pow(_rho, (n - 2));
						_solid_tide_Snm(n, m + 1) += norm(n, m + 1) * _factor * pnm(n, m + 1) * sin(m * _lon) * pow(_rho, (n - 2));
					}
				}
				if (_poletide)
				{
					for (int m = 0; m < 3; m++)
					{
						_solid_tide_Cnm(4, m + 1) += norm(2, m + 1) * _factor * pnm(2, m + 1) * cos(m * _lon) * love_2p[m];
						_solid_tide_Snm(4, m + 1) += norm(2, m + 1) * _factor * pnm(2, m + 1)  *sin(m * _lon) * love_2p[m];
					}
				}
			}

			for (int n = 2; n < 4; n++)
			{
				for (int m = 0; m < n + 1; m++)
				{
					_solid_tide_Cnm(n, m + 1) = LOVE_NO[n - 1][m] * _solid_tide_Cnm(n, m + 1);
					_solid_tide_Snm(n, m + 1) = LOVE_NO[n - 1][m] * _solid_tide_Snm(n, m + 1);
					CF(n, m + 1) = 0;
					SF(n, m + 1) = 0;
				}
			}
			ColumnVector beta;
			_cal_beta_arg(tdt, gst, beta);
			int index = 0;
			if(_poletide)
			{
				for (int i = 0; i < 72; i++)
				{
					theta += ETF[0][i] * beta(1);
					theta += ETF[1][i] * beta(2);
					theta += ETF[2][i] * beta(3);
					theta += ETF[3][i] * beta(4);
					theta += ETF[4][i] * beta(5);
					theta += ETF[5][i] * beta(6);

					sin_theta = sin(theta);
					cos_theta = cos(theta);
					index = int(ETF[0][i] + 0.01);
					if (index == 0)
					{
						CF(2, index + 1) = CF(2, index + 1) + (ETF[6][i] * cos_theta - ETF[7][i] * sin_theta) * 1E-12;
					}
					else
					{
						CF(2, index + 1) += (ETF[6][i] * sin_theta + ETF[7][i] * cos_theta) * 1E-12;
						//yjqin
						SF(2, index + 1) += (ETF[6][i] * cos_theta - ETF[7][i] * sin_theta) * 1E-12;
						//SF(2, index + 1) += (ETF[6][i] * sin_theta - ETF[7][i] * cos_theta) * 1E-12;
					}
					theta = 0.0;
				}
			}
			else
			{

				for (int i = 0; i < 47; i++)
				{
					theta += ETF[0][i] * beta(1);
					theta += ETF[1][i] * beta(2);
					theta += ETF[2][i] * beta(3);
					theta += ETF[3][i] * beta(4);
					theta += ETF[4][i] * beta(5);
					theta += ETF[5][i] * beta(6);

					sin_theta = sin(theta);
					cos_theta = cos(theta);
					index = int(ETF[0][i] + 0.01);
					CF(2, index + 1) += (ETF[6][i] * sin_theta + ETF[7][i] * cos_theta) * 1E-12;
					SF(2, index + 1) += (ETF[6][i] * sin_theta - ETF[7][i] * cos_theta) * 1E-12;
					theta = 0.0;
				}
			}

			for (int n = 2; n < 4; n++)
			{
				for (int m = 0; m < n + 1; m++)
				{
					_solid_tide_Cnm(n, m + 1) += CF(n, m + 1);
					_solid_tide_Snm(n, m + 1) += SF(n, m + 1);

				}
			}
			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::_calculate_solid_tide ==> throw exception");
			return false;
		}
	}

	/** caculate normalization coefficient.
	* @param [in] mode Type of planetary ephemeris.
	* @param [in] ndegree The ocean degree n that we used.
	* @param [in] morder  The ocean degree m taht we used.
	* @param [in] ndim    Max ocean degree n.
	* @param [in] ndim    Max ocean degree m.
	* @param [out] norm   Normalize coefficients.
	*/
	bool t_gfgravity_egm08::_normalize_coeff(const string& mode, const int& ndegree, const int& morder, const int& ndim, const int& mdim, Matrix& norm)
	{
		try
		{
			norm.ReSize(ndim, mdim + 1);
			norm = 0.0;

			for (int n = 1; n <= ndegree; n++)
			{
				norm(n, 1) = sqrt(2.0 * n + 1.0);
				norm(n, 2) = sqrt((2.0 * n + 1.0) / (0.5 * n * (n + 1.0)));
			}

			for (int m = 2; m <= morder; m++)
			{
				for (int n = m; n <= ndegree; n++)
				{
					norm(n, m + 1) = norm(n, m) / sqrt(1.0 * ((int64_t)n + m) * ((int64_t)n - m + 1));
				}
			}

			if (mode.substr(0, 2) == "de")
			{
				for (int m = 1; m <= ndegree; m++)
				{
					for (int n = 0; n <= morder; n++)
					{
						if (n <= m)
						{
							norm(m, n + 1) = 1.0 / norm(m, n + 1);
						}
					}
				}
			}
			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::_normalize_coeff ==> throw exception");
			return false;
		}
	}

	/** caculate beta angle.
	* @param [in] tdt Dynamic time.
	* @param [in] gst Greenwich sidereal time.
	* @param [out] beta angle.
	*/
	void t_gfgravity_egm08::_cal_beta_arg(const double& tdt, const double& gst, ColumnVector& beta)
	{
		double arg[6] = { 0.0 };
		//应使用从J2000起算的儒略世纪数，而非简化儒略日
		//double t = tdt;
		double t = (tdt - 51544.5) / 36525;
		arg[1] = fmod(485868.249036 +
			t * (1717915923.2178 +
				t * (31.8792 +
					t * (0.051635 +
						t * (-0.00024470)))), TURNAS) * SEC2RAD;

		arg[2] = fmod(1287104.79305 +
			t * (129596581.0481 +
				t * (-0.5532 +
					t * (0.000136 +
						t * (-0.00001149)))), TURNAS) * SEC2RAD;

		arg[3] = fmod(335779.526232 +
			t * (1739527262.8478 +
				t * (-12.7512 +
					t * (-0.001037 +
						t * (0.00000417)))), TURNAS) * SEC2RAD;

		arg[4] = fmod(1072260.70369 +
			t * (1602961601.2090 +
				t * (-6.3706 +
					t * (0.006593 +
						t * (-0.00003169)))), TURNAS) * SEC2RAD;

		arg[5] = fmod(450160.398036 +
			t * (-6962890.5431 +
				t * (7.4722 +
					t * (0.007702 +
						t * (-0.00005939)))), TURNAS) * SEC2RAD;

		beta.resize(6);
		beta = 0.0;
		beta(2) = arg[3] + arg[5];
		beta(3) = beta(2) - arg[4];
		beta(4) = beta(2) - arg[1];
		beta(5) = -arg[5];
		beta(6) = beta(2) - arg[4] - arg[2];
		beta(1) = gst + atan(1.0) * 4.0 - beta(2);
	}

	/** caculate Pnm coefficients.
	* @remarks Specific formulas can be referred to Geodetic books.
	* @param [in] tdt Dynamic time.
	* @param [in] name Planet name.
	* @param [out] Pnm coefficients.
	*/
	bool t_gfgravity_egm08::_cal_pnm(const string& name, const double& tdt, Matrix& pnm)
	{
		try
		{
			t_gtriple xj_model(0.0, 0.0, 0.0);
			t_gtriple xj_center(0.0, 0.0, 0.0);

			ColumnVector xe(0.0, 3);
			if (_planet_crd_crs.count(name) > 0)
			{
				xe = _planet_crd_crs[name];
			}
			else
			{
				write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::_cal_pnm ==> There is no data for " + name + ", can not get the pos in J2000");
				return false;
			}

			pnm.ReSize(3, 4);
			_gplanets->get_pos(tdt, name, xj_model);
			_gplanets->get_pos(tdt, "EARTH", xj_center);

			double dist = (xj_model - xj_center).norm();

			double sin_lat = xe(3) / dist;
			_lon = atan2(xe(2), xe(1));

			if (dist == 0.0)
			{
				write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::_cal_pnm ==> Can not get the sun/moon and earth pos in J2000, check your DE file");
				return false;
			}
			else
			{
				_rho = EARTH_R / dist;
			}


			// jdhuang
			double model_GM = _gplanets->get_GM(name);
			if (model_GM == 0.0)
			{
				write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::_cal_pnm ==> Can not get the GM of " + name);
				return false;
			}
			else
			{
				_factor = model_GM / EARTH_GM * _rho * _rho * _rho;
			}

			pnm = 0.0;
			double sin2_lat = sin_lat * sin_lat;
			pnm(2, 1) = 0.5 * (3.0 * sin2_lat - 1.0);
			pnm(2, 2) = 3.0 * sqrt(1.0 - sin2_lat) * sin_lat;
			pnm(2, 3) = 3.0 * (1.0 - sin2_lat);
			pnm(3, 1) = 0.5 * (5.0 * sin2_lat - 3.0) * sin_lat;
			pnm(3, 2) = 1.5 * sqrt(1.0 - sin2_lat) * (5.0 * sin2_lat - 1.0);
			pnm(3, 3) = 15.0 * (1.0 - sin2_lat) * sin_lat;
			pnm(3, 4) = 15.0 * sqrt(pow((1.0 - sin2_lat), 3));

			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 0, "ERROR : t_gfgravity_egm08::_cal_pnm ==> throw exception");
			return false;
		}
	}

	/** caculate ocean tide coefficients.
	* @remarks Specific formulas can be referred to Geodetic books.
	* @param [in] tdt Dynamic time.
	* @param [in] gst Greenwich sidereal time.
	*/
	bool t_gfgravity_egm08::_calculate_ocean_tide(const double& tdt, const double& gst)
	{
		double theta = 0.0;
		double sin_theta = 0.0;
		double cos_theta = 0.0;
		int n, m = 0;

		try
		{
			double fact = 16 * atan(1.0) * 1025.0 * 6.67259e-11 / 9.780327 * 0.01;
			Matrix fnm(0.0, _max_ocean_degree, _max_ocean_degree + 1);
			_ocean_tide_Cnm = 0.0;
			_ocean_tide_Snm = 0.0;

			Matrix norm;
			_normalize_coeff("denorm", _ocean_degree, _ocean_degree, _max_ocean_degree, _max_ocean_degree, norm);

			_ocean_load = _gotdata->load();
			_ocean_tide = _gotdata->oct();
			for (n = 2; n <= _ocean_degree; n++)
			{
				for (m = 0; m <= _ocean_degree; m++)
				{
					fnm(n, m + 1) = fact * norm(n, m + 1) * (1.0 + _ocean_load[(int64_t)n - 1]) / (1.0 * ((int64_t)2 * n + 1));
				}
			}

			ColumnVector beta;
			_cal_beta_arg(tdt, gst, beta);
			_max_constits = _ocean_tide.size();

			double tmp_coef = 0.0;

			for (int i = 0; i < _max_constits; i++)
			{
				theta += _ocean_tide[i].index[0] * beta(1);
				theta += _ocean_tide[i].index[1] * beta(2);
				theta += _ocean_tide[i].index[2] * beta(3);
				theta += _ocean_tide[i].index[3] * beta(4);
				theta += _ocean_tide[i].index[4] * beta(5);
				theta += _ocean_tide[i].index[5] * beta(6);

				sin_theta = sin(theta);
				cos_theta = cos(theta);

				n = (_ocean_tide)[i].n;
				m = (_ocean_tide)[i].m;

				tmp_coef = 0.0;
				tmp_coef += ((_ocean_tide)[i].coeff[0] + (_ocean_tide)[i].coeff[1]) * cos_theta;
				tmp_coef += ((_ocean_tide)[i].coeff[2] + (_ocean_tide)[i].coeff[3]) * sin_theta;
				_ocean_tide_Cnm(n, m + 1) += fnm(n, m + 1) * tmp_coef;

				if (m > 0)
				{
					tmp_coef = 0.0;
					tmp_coef += ((_ocean_tide)[i].coeff[2] - (_ocean_tide)[i].coeff[3]) * cos_theta;
					tmp_coef -= ((_ocean_tide)[i].coeff[0] - (_ocean_tide)[i].coeff[1]) * sin_theta;
					_ocean_tide_Snm(n, m + 1) += fnm(n, m + 1) * tmp_coef;
				}
				theta = 0.0;
			}
			return true;
		}
		catch (...)
		{
			write_log(_glog, 0, "ERROR : t_gfgravity_egm08::_calculate_ocean_tide ==> throw exception");
			return false;
		}
	}
	bool t_gfgravity_egm08::_calculate_ocean_pole_tide(const double& m1, const double& m2, const int& ndegree)
	{
		double a, b;
		_ocean_pole_tide_Cnm = 0.0;
		_ocean_pole_tide_Snm = 0.0;
		double ae = 6378136.600000;
		double ge = 9.78032780000000;
		double G = 6.6742800000e-11;
		double gme = 3.98600441500000e+14;
		double wearth = 7.29211514670000e-5;
		double fact1 = wearth * wearth*pow(ae, 4) / gme;
		double fact2 = 16 * atan(1.0)*G*1025.0 / ge;
		double das2r = 4.848136811095359935899141e-6;
		double rn[30];
		map<int, int> nm;
		for (int i = 0; i < 30; i++)
		{
			rn[i] = fact1 * fact2*(1 + _load[i]) / (2.0 * i + 3);
		}
		for (int n = 1; n <= 30; n++)
		{
			for (int m = 0; m <= n; m++)
			{
				a= rn[n - 1] * (_AR[n][m] * (0.6870*m1 + 0.0036*m2) + _AI[n][m] * (0.6870*m2 - 0.0036*m1))*das2r;
				b= rn[n - 1] * (_BR[n][m] * (0.6870*m1 + 0.0036*m2) + _BI[n][m] * (0.6870*m2 - 0.0036*m1))*das2r;
				_ocean_pole_tide_Cnm(n, m + 1) = a;
				_ocean_pole_tide_Snm(n, m + 1) = b;
			}
		}
		return true;
	}
}