/**
 *
 * @verbatim
	 History
	  -1.0 Yujie Qin and Bo Wong   2019-01-06 Creat the file.
	  -1.1 Hongmin Zhang           2019-03-04 update.
	  -1.2 Yujie Qin               2019-03-05 edit for final version.
	  -1.3 Yujie Qin               2019-03-15 Code comments are added in Doxygen format.
	  -1.4 Jiande Huang             2019-04-01 Adding abnormal throwing mechanism of functions

   @endverbatim
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		    gfatmosphere.cpp
 * @brief		    The class for atmospheric drag mode.
*
 * @author          Yujie Qin, Wuhan University
 * @version		    1.0.0
 * @date		    2019-01-06
 *
 */
#include "gforce/gfatmosphere.h"
#include "gutils/ginfolog.h"

namespace great
{
	/** constructor. */
	t_gfatmosphere::t_gfatmosphere(t_init_arg* arg, t_glog* log) : t_gforce(log)
	{
		// set the function name 
		string function_name = "t_gfatmosphere::t_gfatmosphere";

		// change the type of argument
		t_init_atmosphere *atmos = NULL;
		atmos = dynamic_cast<t_init_atmosphere*>(arg);

		// if no argument
		if (!atmos || !atmos->gsolmag ||!atmos->gpannel)
		{
			write_log(log, 4, function_name, "ERROR : init failed ");
			if (!atmos)
			{
				write_log(log, 4, function_name, "ERROR : no init argument ");
			}
			// if no data in argument
			if (!atmos->gsolmag)
			{
				write_log(log, 4, function_name, "ERROR : no solar data ");
			}
			if (!atmos->gpannel)
			{
				write_log(log, 4, function_name, "ERROR : no pannel data ");
			}
		}
		else
		{
			// set the init value
			_gsolmag = atmos->gsolmag;
			_gpannel = atmos->gpannel;

			_sat_pos.resize(3);
			_sat_vel.resize(3);
			_sun_pos.resize(3);

			_sat_pos = 0;
			_sat_vel = 0;
			_sun_pos = 0;

			_mass = atmos->mass;
			_model = atmos->model;
			_gsolmag = atmos->gsolmag;
			_gpannel = atmos->gpannel;

			for (int i = 0; i < 7; i++)
				_ap.push_back(0.0);

		}
	}



	/** 
	 * defalut deconstructor. 
	 */
	t_gfatmosphere::~t_gfatmosphere()
	{

	}

	/** calculate the accleration and the corresponding partial matrix.
	 *@param [in] argument Parameters to be used in calculation.
	 */
	bool t_gfatmosphere::calculate(t_gforce_calarg* argument)
	{
		// set the function name 
		string function_name = "t_gfatmosphere::calculate";
		try
		{
			// reset zero for safe
			_acc = 0.0;

			// if no argument
			t_gfcalarg_atmosphere* arg = NULL;
			arg = dynamic_cast<t_gfcalarg_atmosphere*>(argument);
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

			// get the data form argument
			_sat_pos = *(arg->sat_pos_trs);
			_sat_vel = *(arg->sat_vel_trs);
			_sun_pos = *(arg->sun_pos);
			_epoch   = *(arg->time);
			_drag    = *(arg->atmos_par);

			//calculate kp and solar 
			if (!_calculate_sol_geo(_epoch))
			{
				return false;
			}

			//get density
			if (_model == "DTM94")
			{
				t_gfatmos_dtm94* dtm94 = new t_gfatmos_dtm94(_glog);
				_density = dtm94->cal_density(_sat_pos, _sun_pos, _epoch, _kp, _solar);
				delete dtm94;
				dtm94 = NULL;
			}
			else if (_model == "MSISE00")
			{
				t_gfatmos_nrlmsise00* msise00 = new t_gfatmos_nrlmsise00(_glog);
				ColumnVector temp_xsc, gsc;
				Matrix xics;
				double lat, lon, stl, alt;
				vector<double> density, t;
				xics.resize(3, 1);
				xics(1, 1) = arg->xics_x;
				xics(2, 1) = arg->xics_y;
				xics(3, 1) = arg->xics_z;
				temp_xsc << _sat_pos;
				temp_xsc *= 1000.0;
				msise00->xyz2blh(temp_xsc, gsc);
				lat = gsc(1)*R2D;
				lon = gsc(2)*R2D;
				if (lon > 180)
					lon = lon - 360.0;
				else if (lon < -180)
					lon = lon + 360.0;
				stl = (_epoch - (int)_epoch) * 86400 / 3600 + lon / 15.0;
				alt = sqrt(dotproduct(xics, xics)) - 6.37813646e3;
				msise00->gtd7d(arg->year, arg->doy, (_epoch - (int)_epoch) * 86400, alt, lat, lon, stl, _solar[1], _solar[0], _ap, density, t);
				_density = density[5];
				delete msise00;
				msise00 = NULL;
			}
			else if (_model == "DTM2013")
			{
				t_gfatmos_dtm2013* dtm2013 = new t_gfatmos_dtm2013(_glog);
				//int hour, min;
				//double sec;
				float lat, lon, stl, alt;
				Matrix xics;
				ColumnVector temp_xsc, gsc;				
				/*hour = (int)((_epoch - (int)_epoch)*3600.0);
				min = (int)(((_epoch - (int)_epoch)*3600.0 - hour)*60.0);
				sec = _epoch - (int)_epoch - 3600.0*hour - 60.0*min;*/
				xics.resize(3, 1);
				xics(1, 1) = arg->xics_x;
				xics(2, 1) = arg->xics_y;
				xics(3, 1) = arg->xics_z;
				alt = sqrt(dotproduct(xics, xics)) - 6.37813646e3;				
				temp_xsc << _sat_pos;
				temp_xsc *= 1000.0;
				dtm2013->xyz2blh(temp_xsc, gsc);
				lat = gsc(1);
				lon = gsc(2);
				stl = (_epoch - (int)_epoch) * 24 * 15 * G_PI / 180.0 + lon;
				cout << setw(25) << setprecision(16) << arg->doy + _epoch - (int)_epoch << _solar[0] << _solar[1] << _akp[0] << _akp[2] << alt << stl << lat << lon << endl;
				_density = dtm2013->calDensity(arg->doy + _epoch - (int)_epoch, _solar[0], _solar[1], _akp, alt, stl, lat, lon); //arg->doy, arg->mon, arg->year, hour, min, sec,
				delete dtm2013;
				dtm2013 = NULL;
			}
			else
			{
				_density = 2.125e-11;
			}


			// wind velocity ( not used)
			ColumnVector wind_vel;
			wind_vel.resize(3);
			wind_vel = 0;

			// sat velocity in earth fixed system
			ColumnVector efs_sat_vel = *(arg->e2j) * _sat_vel;
			double vel_norm = efs_sat_vel.norm_Frobenius();
			efs_sat_vel = efs_sat_vel / vel_norm;				// unit vector

			// effect area (sv_norm must be in inertial system)
			ColumnVector xtmp(0.0,3);
			ColumnVector gsun(0.0,3);

			//take all the areas of satellite into consideration
			double cosa = 0.0;
			double area = 0.0;
			map<int, t_oneside>::iterator it_pannel = _gpannel->get_pannels()->begin();
			for (it_pannel; it_pannel != _gpannel->get_pannels()->end(); it_pannel++)
			{
				xtmp = it_pannel->second.crs_norm_dir;
				cosa = dotproduct(xtmp, xtmp);
				if (cosa < 1.0e-10)
				{
					xtmp = efs_sat_vel;
				}
				cosa = dotproduct(efs_sat_vel, xtmp);
				if (cosa < 0)
				{
					continue;
				}

				gsun = crossproduct(efs_sat_vel, xtmp);
				gsun = crossproduct(gsun, efs_sat_vel);
				area = it_pannel->second.area * it_pannel->second.drag_coef;
				
				_acc += area * cosa*efs_sat_vel;
			}
			
			_acc = (_mass == 0.0) ? _acc * 0.0 : _acc * -1100.0 * _density * vel_norm * vel_norm / _mass;

			_partial = _acc;
			_acc    *= _drag;

			return true;
		}
		catch (...)
		{
			write_log(_glog, 4, function_name, "ERROR : throw exception while calculating atmosphere acc");
			return false;
		}
	}

	/**calculate kp and solar of reference time.
	*@param [in] time the current gps_t.
	*/
	bool t_gfatmosphere::_calculate_sol_geo(double time)
	{
		// set the function name 
		string function_name = "t_gfatmosphere::calculate_sol_geo";
		
		try
		{
			_solar.clear();
			int flag_sol;
			int flag_mag;

			static bool lfirst = true;
			//the data in file solar_flux
			static map<double, double[2]>  solar_data = _gsolmag->get_solar_data();
			//the data in file geomag_kp
			static map<double, double>     magkp_data = _gsolmag->get_geomag_data();

			if (_model == "DTM94")
			{
				flag_sol = 1;
				flag_mag = 1;
			}
			else if (_model == "MSISE00")
			{
				flag_sol = 1;
				flag_mag = 2;
			}
			else if (_model == "DTM2013")
			{
				flag_sol = 1;
				flag_mag = 3;
			}
			
			if (flag_sol == 1)
			{
				double interval_sol = _gsolmag->get_intv("solar_flux");
				double t_solar;
				if (_model == "DTM94")
					t_solar = _epoch - interval_sol;
				else if (_model == "MSISE00" || _model == "DTM2013")
					t_solar = _epoch;
				
				map<double, double[2]>::iterator beg_solar, end_solar;     //time
				auto iter1 = solar_data.begin();
			auto iter2 = solar_data.end();
			--iter2;
			if (iter1->first > t_solar || iter2->first < t_solar)
			{
				write_log(_glog, 4, function_name, "ERROR : Input time beyond the time table of solar_flux file");
				return false;
			}
				beg_solar = solar_data.lower_bound(t_solar);
				end_solar = solar_data.upper_bound(t_solar);
				if (beg_solar == end_solar)
				{
					beg_solar--;
				}

				//the value at the start and end time that used for linear interpolation
				double solar1_beg = beg_solar->second[0];
				double solar2_beg = beg_solar->second[1];
				double solar1_end = end_solar->second[0];
				double solar2_end = end_solar->second[1];

				double dcalso1 = _linear_interpolation(solar1_beg, solar1_end, beg_solar->first, end_solar->first, t_solar);
				double dcalso2 = _linear_interpolation(solar2_beg, solar2_end, beg_solar->first, end_solar->first, t_solar);

				_solar.push_back(dcalso1);
				_solar.push_back(dcalso2);
			}
			//interval in solar_flux
			
			if (flag_mag == 1)
			{
				double interval_mag = _gsolmag->get_intv("geomag_kp");
				double t_kp = _epoch - interval_mag;
                auto iter3 = magkp_data.begin();
			auto iter4 = magkp_data.end();
			--iter4;
			if (iter3->first > t_kp || iter4->first < t_kp)
			{
				write_log(_glog, 4, function_name, "ERROR : Input time beyond the time table of geomag_kp file");
				return false;
			}
				map<double, double>::iterator beg_kp = magkp_data.lower_bound(t_kp);
				map<double, double>::iterator end_kp = magkp_data.upper_bound(t_kp);
				if (beg_kp == end_kp)
				{
					beg_kp--;
				}

				//the value at the start and end time that used for linear interpolation
				double kp_beg_value = beg_kp->second;
				double kp_end_value = end_kp->second;

				_kp = _linear_interpolation(kp_beg_value, kp_end_value, beg_kp->first, end_kp->first, t_kp);
			}
			else if (flag_mag == 2)
			{
				//double interval_mag = _gsolmag->interval("geomag_ap");
				double t_ap = _epoch;
				double t_ap0 = t_ap - 3;
				int    it_ap0 = (int)t_ap0;
				double t_ap1 = t_ap + 1;
				int    it_ap1 = (int)t_ap1;
				t_ap1 = (double)it_ap1; // - interval_mag
				int    lval;
				map<double, double>::iterator beg_ap = magkp_data.lower_bound((double)it_ap0); //?应该取一个32大小的数组出来插值
				map<double, double>::iterator end_ap = magkp_data.lower_bound(t_ap1);
				map<double, double>::iterator fb_ap = magkp_data.lower_bound(t_ap);
				map<double, double>::iterator fe_ap = magkp_data.upper_bound(t_ap);
				vector<double> tmp_ap;
				if (fb_ap == fe_ap)
				{
					fb_ap--;
				}
				map<double, double>::iterator it = beg_ap;
				for (it; it != end_ap; it++)
				{
					tmp_ap.push_back(it->second);
					if (it == fb_ap)
						lval = tmp_ap.size();
				}
				//the value at the start and end time that used for linear interpolation
				_linear_interpolation(tmp_ap, lval);
			}
			else if (flag_mag == 3)
			{
				double interval_mag = _gsolmag->get_intv("geomag_kp");
				double t_kp = _epoch-interval_mag;
				double t_kp0 = t_kp - 1;
				vector<double> tmp_kp;
				map<double, double>::iterator beg_kp = magkp_data.lower_bound(t_kp); //?应该取一个32大小的数组出来插值
				map<double, double>::iterator end_kp = magkp_data.upper_bound(t_kp);
				if (beg_kp == end_kp)
				{
					beg_kp--;
				}
				map<double, double>::iterator fb_kp = magkp_data.lower_bound(t_kp0);
				map<double, double>::iterator it = fb_kp;
				for (it; it != end_kp; it++)
				{
					tmp_kp.push_back(it->second);
				}
				_linear_interpolation(beg_kp->second, end_kp->second, beg_kp->first, end_kp->first, t_kp, tmp_kp);
			}
			return true;
		}
		catch (...)
		{
			write_log(_glog, 4, function_name, "ERROR : throw exception while calculate_sol_geo");
			return false;
		}
	}

	/** liner interpolation.
	*@param [in] start_value the start value.
	*@param [in] end_value the end value.
	*@param [in] start_time the start time.
	*@param [in] end_time the end time.
	*@param [in] ref_time the reference time.
	*/
	double t_gfatmosphere::_linear_interpolation(double start_value, double end_value, double start_time, double end_time, double ref_time)
	{
		try
		{
			double alpha = (ref_time - start_time) / (end_time - start_time);
			return start_value + alpha * (end_value - start_value);
		}
		catch(...)
		{
			write_log(_glog, 4, "t_gfatmosphere::LinearInterpolation", "ERROR : throw exception while Linear Interpolation, we return 0.0 as default value");
			write_log(_glog, 4, "NOTE : May be you should check your end time and start time ");
			return 0.0;
		}
	}

	double t_gfatmosphere::_linear_interpolation(vector<double> ap_value, int leftvalue)
	{
		try
		{	
			for (int j = 0; j < 7; j++)
				_ap[j] = 0.0;
			for (int j = 25; j <= 32; j++)
				_ap[0] += ap_value[j - 1] / 8.0;
			_ap[1] = ap_value[leftvalue - 1];
			_ap[2] = ap_value[leftvalue - 2];
			_ap[3] = ap_value[leftvalue - 3];
			_ap[4] = ap_value[leftvalue - 4];
			for (int j = leftvalue - 11; j <= leftvalue - 4; j++)
				_ap[5] += ap_value[j - 1] / 8.0;
			for (int j = leftvalue - 19; j <= leftvalue - 12; j++)
				_ap[6] += ap_value[j - 1] / 8.0;

			return 1;
		}
		catch (...)
		{
			write_log(_glog, 4, "t_gfatmosphere::LinearInterpolation", "ERROR : throw exception while Linear Interpolation, we return 0.0 as default value");
			write_log(_glog, 4, "NOTE : May be you should check your end time and start time ");
			return 0.0;
		}
	}

	double t_gfatmosphere::_linear_interpolation(double start_value, double end_value, double start_time, double end_time, double ref_time, vector<double> kp_value)
	{
		try
		{
			for (int j = 0; j < 4; j++)
				_akp[j] = 0.0;
			if (ref_time >= 57855.0)
				cout << "notice" << endl;
			double alpha = (ref_time - start_time) / (end_time - start_time);
			_akp[0]= start_value + alpha * (end_value - start_value);
			for (int j = 0; j < 8; j++)
				_akp[2] += kp_value[j] / 8;
			return 1;
		}
		catch (...)
		{
			write_log(_glog, 4, "t_gfatmosphere::LinearInterpolation", "ERROR : throw exception while Linear Interpolation, we return 0.0 as default value");
			write_log(_glog, 4, "NOTE : May be you should check your end time and start time ");
			return 0.0;
		}
	}
}