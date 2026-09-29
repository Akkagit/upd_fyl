//#define _DEBUG_LOG
/**
*
* @verbatim
    History
     -1.0 Yujie Qin  2019-01-10  created
     -1.1 Yujie Qin  2019-02-18  optim
     -1.2 Wei Zhang  2019-03-27  Adding Doxygen Style Code Remarks
	 -1.3 Yujie Qin  2019-04-15  debug and optim
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    goi_leo.cpp
* @brief		orbit integration processing for LEO satellite
*
* @author       Yujie Qin, Wuhan University
* @version		1.0.0
* @date		    2019-03-27
*
*/


#include <iostream>
#include <iomanip>
#include <algorithm>

#include "gproc/goi_leo.h"
#include "gforce/gforce.h"
#include "gforce/gfcommon.h"
#include "gset/gsetforce.h"
#include "gutils/ginfolog.h"
#include "newmat/newmatio.h"


namespace great
{
	/** @brief constructor. */
	t_goi_leo::t_goi_leo()
	{

	}

	/** 
	* @brief constructor. 
	* @param[in]	mark   leo name
	* @param[in]	set    setting information
	*/
	t_goi_leo::t_goi_leo(string mark, t_gsetbase* set)
	:t_goi(mark, set)
	{
		gtrace("t_goi_leo::t_goi_leo");

		_sat = mark;
		_epoch.tsys(t_gtime::GPS);
		_set = set;
		_isPannel = true;
	}

	/** @brief destructor. */
	t_goi_leo::~t_goi_leo()
	{

	}

	
	
	/**
	* @brief set LEO ICS data
	* @param[in] gics   LEO ICS data
	*/
	void t_goi_leo::set_data_used(t_gallics_leo* gics)
	{
		_gics = gics;
	}

	

	
	
	/**
	* @brief   get orb funct according to the satellite name 
	* @param[in]  sat   satellite name
	* @param[in]  gorb  orb data
	*/
	void t_goi_leo::getOrbFunct(string sat, t_gorb* gorb)
	{
		map<int, vector< t_gtriple> > funct_tmp;
		_gorb->get_funct(sat, funct_tmp);

		gorb->add_funct(sat, funct_tmp);
	}
	
	/**
	* @brief  get orb header data according to the satellite name
	* @param[in]  sys   system of the satelite
	* @param[in]  sats  all satellite names
	* @param[in]  gorb  orb data
	*/
	void t_goi_leo::getOrbHead(set<string>& sys, set<string>& sats, t_gorb* gorb)
	{
		t_orb_hd orb_hd_tmp;
		string function_name = "t_goi_leo::get_orb_head";
		try
		{
			_gorb->get_head(orb_hd_tmp);
			//sat group
			if (sys.size() > 1)
			{
				orb_hd_tmp.sat_group = "LEO";
			}
			else
			{
				orb_hd_tmp.sat_group = *(sys.begin());
			}
			//sat type
			set<string>::iterator itsys = sys.begin();
			for (itsys; itsys != sys.end(); ++itsys)
			{
				//char sys_tmp = t_gsys::str2char(*itsys);
				orb_hd_tmp.sat_type.insert('L');
			}
			//sat num
			orb_hd_tmp.sat_num = sats.size();
			//sat list
			orb_hd_tmp.sat_prn = sats;
			gorb->add_head(orb_hd_tmp);
		}
		catch (...)
		{
			writeLogInfo(_glog, 2, function_name, "ERROR : unknown mistake");
		}
	}

	/**
	* @brief   init the sat infomation
	* @param[in] sat_name   satellite name
	* @return    initial leo correctly or not
	*/
	bool t_goi_leo::_initSatInfo(string sat_name)
	{
		t_gtime beg;
		t_gtime end;
		string function_name = "t_goi_leo::_init_sat_info";
		try
		{
			beg = _extend_beg + NEXT_POINT_NUM * _step_out;
			end = _extend_end - NEXT_POINT_NUM * _step_out;

			if (!_gallsat->getSatinfoUsed(sat_name, beg, end, _gsat))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not get the sat information of "
					+ sat_name + ", check your sat parameter file or the read process of the file.");
				return false;
			}


			_sat_mass = _gsat.mass();

			double p = 0.0;
			double v = 0.0;

			p += pow(_sat_ics->get_par_value(ICS_PAR::PXSAT, _beg), 2);
			p += pow(_sat_ics->get_par_value(ICS_PAR::PYSAT, _beg), 2);
			p += pow(_sat_ics->get_par_value(ICS_PAR::PZSAT, _beg), 2);

			v += pow(_sat_ics->get_par_value(ICS_PAR::VXSAT, _beg), 2);
			v += pow(_sat_ics->get_par_value(ICS_PAR::VYSAT, _beg), 2);
			v += pow(_sat_ics->get_par_value(ICS_PAR::VZSAT, _beg), 2);

			_period = 2.0 / sqrt(p) - v / EARTH_GM;
			_period = _period * sqrt(_period * EARTH_GM);

			return true;
		}
		catch (...)
		{
			write_log_info(_glog, 3, function_name, "ERROR : throw while getting the sat information");
			return false;
		}

	}

	/**
	* @brief      init the state matrix
	* @param[in]  neq_size   number of OI funct params
	* @param[out] state      state data of OI processing
	* @return  initial state correctly or not
	*/
	bool t_goi_leo::_initStateMat(int neq_size, vector<double>& state)
	{
		string function_name = "t_goi_leo::_init_state_matrix";
		try
		{
			// x,y,z,vx,vy,vz and solar pars will be read from ICS
			// get initial value from ics-file
			state.push_back(_sat_ics->get_par_value(ICS_PAR::PXSAT, _beg));
			state.push_back(_sat_ics->get_par_value(ICS_PAR::PYSAT, _beg));
			state.push_back(_sat_ics->get_par_value(ICS_PAR::PZSAT, _beg));
			state.push_back(_sat_ics->get_par_value(ICS_PAR::VXSAT, _beg));
			state.push_back(_sat_ics->get_par_value(ICS_PAR::VYSAT, _beg));
			state.push_back(_sat_ics->get_par_value(ICS_PAR::VZSAT, _beg));


			// initial partial can be ensured directly
			// storaging with column vector is easy to integral
			for (int i = 1; i < neq_size / 6; ++i)
			{
				for (int j = 1; j < 7; j++)
				{
					if (i == j)
					{
						state.push_back(1.0);
					}
					else
					{
						state.push_back(0.0);
					}
				}
			}
		

#ifdef DEBUG
		cout << "The init state matrix is : " << endl;
		cout << "From 1 - 6 are x, y, z, vx, vy, vz" << endl;
		cout << "Form 7 - 12 or 7 - 16 are 5 or 9 parameters of solar model" << endl;
		cout << "===============================================================" << endl;
		for (int i = 0; i < state.size(); i++)
		{
			cout << "Num = " << i + 1 << setw(30) << setprecision(15) << state[i] << endl;
		}
		cout << "===============================================================" << endl;
#endif // DEBUG
		
		    return true;
		
		}
		catch (...)
		{
			writeLogInfo(_glog, 2, function_name, "ERROR : throw while init_state_matrix");
			return false;
		}
	}

	/**
	* @brief   init the partial matrix
	* @param[in]  par_num   param number
	* @param[out] partial   partial matrix of OI
	* @return  initial leo partial correctly or not
	*/
	bool t_goi_leo::_initPartialMat(int par_num, Matrix & partial)
	{
		partial.resize(par_num + 1, 6);
		partial = 0.0;
		return true;
	}

	/**
	* @brief   Get settings from XML file and set local variables. 
	* @return  get settings correctly or not
	*/
	bool t_goi_leo::_get_settings()
	{

		string function_name = "t_goi_leo::_get_settings";

		// beg and end time
		_beg = dynamic_cast<t_gsetgen*>(_set)->beg();
		_end = dynamic_cast<t_gsetgen*>(_set)->end();

		// settings for intergrator
		_step_inp = dynamic_cast<t_gsetforce*>(_set)->getInpStep(_sat);
		_step_out = dynamic_cast<t_gsetforce*>(_set)->getOutStep(_sat);
		_integra_degree = dynamic_cast<t_gsetforce*>(_set)->getIntegraDegree(_sat);

		// settings for forces
		// List of models for considered forces
		_forces = dynamic_cast<t_gsetforce*>(_set)->getForceUsedSet(_sat);
		_force_model = dynamic_cast<t_gsetforce*>(_set)->getForceModelMap(_sat);
		_force_valid = dynamic_cast<t_gsetforce*>(_set)->getForcesValidMap(_sat);
		_cus_par = dynamic_cast<t_gsetforce*>(_set)->getCustomerPars(_sat);
		_drag_par = dynamic_cast<t_gsetforce*>(_set)->getAtmospherePars(_sat);
		_srp_par = dynamic_cast<t_gsetforce*>(_set)->getSolarRadiationPars(_sat);
		dynamic_cast<t_gsetforce*>(_set)->getGravitySet(_sat, _ocean_degree, _isPoletide, _pole_degree);


#ifdef _DEBUG_LOG

		// Creat and set the log file : oi.log

		_glog->comment(0, _sat, "------------------------------------------------");
		_glog->comment(0, _sat, "beg time       : " + dbl2str(_beg.dmjd()));
		_glog->comment(0, _sat, "ned time       : " + dbl2str(_end.dmjd()));
		_glog->comment(0, _sat, "step_inp       : " + dbl2str(_step_inp));
		_glog->comment(0, _sat, "step_out       : " + dbl2str(_step_out));
		_glog->comment(0, _sat, "integra_degree : " + dbl2str(_integra_degree));
		_glog->comment(0, _sat, "------------------------------------------------");
#endif
		return true;
	}
	bool t_goi_leo::_check_par()
	{
		map<string, int> cuspar_ics;
		map<string, int> drag_ics;
		map<string, int> srp_ics;
		map<ICS_PAR, vector<t_icsleo_record>> * par_ics = _sat_ics->leoics();
		auto par_beg = par_ics->begin();
		auto par_end = par_ics->end();
		ICS_PAR partype;
		int partime;
		for (auto ipar = par_beg; ipar != par_end; ipar++)
		{
			partype = ipar->first;
			if (partype <= ICS_PAR::VZSAT) continue;
			if (partype == ICS_PAR::SR_scale) srp_ics[ics2str(partype)] = ipar->second.begin()->piecewise_time;
			if (partype == ICS_PAR::DRAG_c) drag_ics[ics2str(partype)] = ipar->second.begin()->piecewise_time;
			if (partype > ICS_PAR::DRAG_c) cuspar_ics[ics2str(partype)] = ipar->second.begin()->piecewise_time;
		}

		if (cuspar_ics.size() != _cus_par.size() || drag_ics.size() != _drag_par.size() || srp_ics.size() != _srp_par.size())
			return false;
		if (cuspar_ics == _cus_par && drag_ics == _drag_par && srp_ics == _srp_par)
			return true;
		else
			return false;
	}
	/**
	* @brief   check all file data.
	* @return  get filedata correctly or not
	*/
	bool t_goi_leo::_checkDataUsed()
	{
		if (_gplanet->is_empty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_leo::_check_data ==> no de data");
			return false;
		};			///< data for all planets.

		if (_gpoleut1->isEmpty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_leo::_check_data ==> no poleut1 data");
			return false;
		};			///< poleut1 data.

		if (_gotdata->is_empty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_leo::_check_data ==> no ocean tide data");
			return false;
		};			///< ocean tide data.

		if (_gdeadata->is_empty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_leo::_check_data ==> no ocean pole tide data");
			return false;
		};     ///< ocean pole tide data.
		if (_gleapsecond->is_empty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_leo::_check_data ==> no leapsecond data");
			return false;
		};		///< leap second data.

		if (_gics->isEmpty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_leo::_check_data ==> no ics data");
			return false;
		};				///< all  ics data.

		if (_gegmcoef->is_empty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_leo::_check_data ==> no EGM data");
			return false;
		};			///< EGM data.

		if (_gallsat->isEmpty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_leo::_check_data ==> no sat parameters data");
			return false;
		};			///< satellite parameters data.

		if (!_gotl)
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_leo::_check_data ==> no ocean load data");
			return false;
		};				///< ocean load data.

		if (!_gsolargeomag)
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_leo::_check_data ==> no solar data");
			return false;
		};            ///< solar data.

		if (!_gpannel)
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_leo::_check_data ==> no pannel data");
			return false;
		};            ///< pannel data.

		if (!_gattitude)
		{
			writeLogInfo(_glog, 1, "ERROR : t_goi_leo::_check_data ==> no attitude data");
			//return false; attitude data not neccesarry by zhangwei
		};          ///< pannel data.
		return true;
	}


	/** @brief  Set output file (empty).*/
	void t_goi_leo::_setOut()
	{
	}

	/** @brief  Set processing log file*/
	void t_goi_leo::_setLog()
	{
		
	}

	


	/**
	* @brief  calculate the acc and partial in epoch t
	* @param[in]  t        current epoch for calculation
	* @param[out] state    state data after calculation
	* @param[out] partial  partial data after calculation
	* @return  calculate correctly or not
	*/
	bool t_goi_leo::_calculate(const t_gtime& t, const vector<double>& state, vector<double>& partial)
	{
		string function_name = "t_goi_leo::_calculate";		
		try
		{
			//change time sys to tdt
			t_gtime tdt_gtime = t;
			tdt_gtime.tsys(t_gtime::TT);
			//ChangeTimeSys(GPS2TDT, t, tdt_gtime);
			double tdt = tdt_gtime.dmjd();

			double dt = (t.mjd() - _ref.mjd()) * 86400 + t.sod() + t.dsec() - _ref.sod() - _ref.dsec();
			double gps_t = t.dmjd();
			_trs2crs_cord->calcRotMat(tdt_gtime, true, true, true);

			// Matrix
			ColumnVector satXJ(0.0, 3);       // satellite position in J2000
			ColumnVector satVJ(0.0, 3);       // satellite velocity in J2000
			ColumnVector satXE(0.0, 3);       // satellite position in TRS
			ColumnVector satVE(0.0, 3);       // satellite velocity in TRS

			map<string, ColumnVector> allXE;  // planet position in TRS,used for ocean tide correcting etc.
			map<string, ColumnVector> allXJ;  // planet position in J2000,used for point_mass etc. 

			Matrix rot_scf2crsleo(0.0, 3, 3); //rotation matrix of SCF to CRS, just for leo!
			Matrix rot_trs2crs(0.0, 3, 3);   // rotation matrix of TRS to CRS

			// get information of coordinate transformation from TRS to CRS like rotation matrix of TRS to CRS 
			Matrix _Dgmat = _trs2crs_cord->getMatDu();
			Matrix _Dxmat = _trs2crs_cord->getMatDx();
			Matrix _Dymat = _trs2crs_cord->getMatDy();
			rot_trs2crs = _trs2crs_cord->getRotMat();
			double gst = _trs2crs_cord->getGmst();
			double xpole = _trs2crs_cord->getXpole();
			double ypole = _trs2crs_cord->getYpole();

			// for planet
			ColumnVector sun_pos_crs(0.0, 3);
			ColumnVector sun_vel_crs(0.0, 3);
			ColumnVector moon_pos_crs(0.0, 3);
			ColumnVector moon_vel_crs(0.0, 3);

			for (int i = 0; i < 3; i++)
			{
				satXJ(i + 1) = state[i];
				satVJ(i + 1) = state[i + 3];
			}


			ColumnVector x_earth = (satXJ.t() * _Dgmat).t();
			satXE = (satXJ.t() * rot_trs2crs).t();
			satVE = (satVJ.t() * rot_trs2crs).t() + x_earth / 43200.0 * G_PI;
			// get the sun and moon pos in crs
			_gplanet->get_pos(tdt, "SUN", sun_pos_crs);
			_gplanet->get_pos(tdt, "MOON", moon_pos_crs);
			allXE["SUN"] = (sun_pos_crs.t() * rot_trs2crs).t();
			allXJ["SUN"] = sun_pos_crs;
			allXE["MOON"] = (moon_pos_crs.t() * rot_trs2crs).t();
			allXJ["MOON"] = moon_pos_crs;

#ifdef _DEBUG_LOG_1
			cout << "xpole ypole" << endl;
			cout << setw(25) << scientific << setprecision(16) << xpole << endl;
			cout << setw(25) << scientific << setprecision(16) << ypole << endl;
			cout << "e2j" << endl;
			for (int i = 1; i <= rot_trs2crs.Nrows(); i++)
			{
				for (int j = 1; j <= rot_trs2crs.Ncols(); j++)
				{
					cout << setw(25) << scientific << setprecision(16) << rot_trs2crs(j, i);
				}
				cout << endl;
			}

			cout << "dgmat" << endl;
			for (int i = 1; i <= _Dgmat.Nrows(); i++)
			{
				for (int j = 1; j <= _Dgmat.Ncols(); j++)
				{
					cout << setw(25) << scientific << setprecision(16) << _Dgmat(j, i);
				}
				cout << endl;
			}

			cout << "satXJ" << endl;
			for (int i = 1; i <= 3; i++)
			{
				cout << scientific << setw(30) << setprecision(20) << satXJ(i) << "  ";
			}
			cout << endl;

			cout << "satVJ" << endl;
			for (int i = 1; i <= 3; i++)
			{
				cout << scientific << setw(30) << setprecision(20) << satVJ(i) << "  ";
			}
			cout << endl;
			cout << "satXE" << endl;
			for (int i = 1; i <= 3; i++)
			{
				cout << scientific << setw(30) << setprecision(20) << satXE(i) << "  ";
			}
			cout << endl;

			cout << "satVE" << endl;
			for (int i = 1; i <= 3; i++)
			{
				cout << scientific << setw(30) << setprecision(20) << satVE(i) << "  ";
			}
			cout << endl;
			cout << "sunXE" << endl;
			for (int i = 1; i <= 3; i++)
			{
				cout << scientific << setw(30) << setprecision(20) << allXE["SUN"](i) << "  ";
			}
			cout << endl;
#endif

			//prepare for solar
			double shadow_factor_boxwing = get_shadow_factor(satXJ, sun_pos_crs, moon_pos_crs);

			//get the param from ics-file according the time
			map<ICS_PAR, double> solar_par;
			if (_srp_par.size() > 0) 
				solar_par[ICS_PAR::SR_scale] = _sat_ics->get_par_value(ICS_PAR::SR_scale, t);
			map<ICS_PAR, double> cust_par = _sat_ics->get_pars_value("EMP", t);
			double               xics_x = _sat_ics->get_par_value(ICS_PAR::PXSAT, t);
			double               xics_y = _sat_ics->get_par_value(ICS_PAR::PYSAT, t);
			double               xics_z = _sat_ics->get_par_value(ICS_PAR::PZSAT, t);

			double atmos_par = _sat_ics->get_par_value(ICS_PAR::DRAG_c, t);

			// get scf2crs
			_scf2crsleo_cord->calcRotMat(t, _gattitude, satXJ, satVJ, rot_trs2crs, rot_scf2crsleo);

#ifdef _DEBUG_LOG_1
			cout << "scf2crs is" << endl;
			for (int i = 1; i <= rot_scf2crsleo.Nrows(); i++) {
				for (int j = 1; j <= rot_scf2crsleo.Ncols(); j++) {
					cout << setw(25) << scientific << setprecision(16) << rot_scf2crsleo(i, j);
				}
				cout << endl;
			}
#endif
			ColumnVector sat2sun = sun_pos_crs - satXJ;
			ColumnVector unit_sat2sun = sat2sun / sat2sun.norm_frobenius();
			if (_isPannel)
			{
				_gpannel->get_rotmat(rot_scf2crsleo, unit_sat2sun);
			}
			// calculate acc
			ColumnVector _acc_point_mass(0.0, 3);
			ColumnVector _acc_gravity(0.0, 3); 
			ColumnVector _acc_relativity(0.0, 3);
			ColumnVector _acc_solar_radi(0.0, 3); 
			ColumnVector _acc_customer(0.0, 3);
			ColumnVector _acc_atmosphere(0.0, 3); 
			ColumnVector _acc_earthradiation(0.0, 3);
			if (_force_valid[POINT_MASS] == true)
			{
				t_gfcalarg_point_mass* arg_point_mass = new t_gfcalarg_point_mass();
				arg_point_mass->arg_type = ARG_POINT_MASS;
				arg_point_mass->time = &tdt;
				arg_point_mass->sat_pos = &satXJ;
				arg_point_mass->sat_vel = &satVJ;
				if (!_gforce.point_mass->calculate(arg_point_mass))
				{
					writeLogInfo(_glog, 4, function_name, "ERROR : can not calculate point mass");
					delete arg_point_mass; arg_point_mass = NULL;
					return false;
				}
				delete arg_point_mass; arg_point_mass = NULL;
				_acc_point_mass = _gforce.point_mass->acc();
			}


			if (_force_valid[RELATIVITY] == true)
			{
				t_gfcalarg_relativity* arg_relativity = new t_gfcalarg_relativity();
				arg_relativity->arg_type = ARG_RELATIVITY;
				arg_relativity->time = &tdt;
				arg_relativity->sat_pos = &satXJ;
				arg_relativity->sat_vel = &satVJ;
				if (!_gforce.relativity->calculate(arg_relativity))
				{
					writeLogInfo(_glog, 4, function_name, "ERROR : can not calculate relativity");
					delete arg_relativity; arg_relativity = NULL;
					return false;
				};
				delete arg_relativity; arg_relativity = NULL;
				_acc_relativity = _gforce.relativity->acc();
			}



			if (_force_valid[SOLAR_RADI] == true)
			{
				t_gfcalarg_solar_boxwing* arg_solar_boxwing = new t_gfcalarg_solar_boxwing();
				arg_solar_boxwing->arg_type = ARG_SOLAR_BOXWING;
				arg_solar_boxwing->time = &tdt;
				arg_solar_boxwing->sun_pos = &sun_pos_crs;
				arg_solar_boxwing->sat_pos = &satXJ;
				arg_solar_boxwing->sat_vel = &satVJ;
				arg_solar_boxwing->solar_par = &solar_par;
				arg_solar_boxwing->shadow_factor = &shadow_factor_boxwing;
				if (!_gforce.solar_radi->calculate(arg_solar_boxwing))
				{
					writeLogInfo(_glog, 4, function_name, "ERROR : can not calculate relativity");
					delete arg_solar_boxwing; arg_solar_boxwing = NULL;
					return false;
				};
				delete arg_solar_boxwing; arg_solar_boxwing = NULL;
				_acc_solar_radi = _gforce.solar_radi->acc();
			}



			if (_force_valid[GRAVITY] == true)
			{
				t_gfcalarg_gravity_egm* arg_gravity = new t_gfcalarg_gravity_egm();
				arg_gravity->arg_type = ARG_GRAVI_EGM;
				arg_gravity->xpole = &xpole;
				arg_gravity->ypole = &ypole;
				arg_gravity->sat_pos_trs = &satXE;
				arg_gravity->e2j = &rot_trs2crs;
				arg_gravity->time = &tdt;
				arg_gravity->sat_pos = &satXJ;
				arg_gravity->sat_vel = &satVJ;
				arg_gravity->gst = &gst;
				arg_gravity->allXE = &allXE;
				if (!_gforce.gravity->calculate(arg_gravity))
				{
					writeLogInfo(_glog, 4, function_name, "ERROR : can not calculate solar_radi");
					delete arg_gravity; arg_gravity = NULL;
					return false;
				}
				delete arg_gravity; arg_gravity = NULL;
				_acc_gravity = _gforce.gravity->acc();
			}



			if (_force_valid[ATMOSPHERE] == true)
			{
				t_gfcalarg_atmosphere* arg_atmosphere = new t_gfcalarg_atmosphere();
				arg_atmosphere->arg_type = ARG_ATMOSPHERE;
				arg_atmosphere->sat_pos = &satXJ;
				arg_atmosphere->sat_vel = &satVJ;
				arg_atmosphere->sat_pos_trs = &satXE;
				arg_atmosphere->sat_vel_trs = &satVE;
				arg_atmosphere->time = &gps_t;
				arg_atmosphere->year = t.year();
				arg_atmosphere->mon = t.mon();
				arg_atmosphere->doy = t.doy();
				arg_atmosphere->xics_x = xics_x;
				arg_atmosphere->xics_y = xics_y;
				arg_atmosphere->xics_z = xics_z;
				arg_atmosphere->sun_pos = &allXE["SUN"];
				arg_atmosphere->e2j = &rot_trs2crs;
				arg_atmosphere->atmos_par = &atmos_par;
				if (!_gforce.atmosphere->calculate(arg_atmosphere))
				{
					writeLogInfo(_glog, 4, function_name, "ERROR : can not calculate atmosphere");
					delete arg_atmosphere; arg_atmosphere = NULL;
					return false;
				}
				delete arg_atmosphere; arg_atmosphere = NULL;
				_acc_atmosphere = _gforce.atmosphere->acc();
			}


			if (_force_valid[EARTH_RADI] == true)
			{
				t_gfcalarg_earthradiation* arg_earthradiation = new t_gfcalarg_earthradiation();
				arg_earthradiation->arg_type = ARG_EARTHRADIATION;
				arg_earthradiation->sat_pos = &satXJ;
				arg_earthradiation->sat_vel = &satVJ;
				arg_earthradiation->sun_pos_crs = &sun_pos_crs;
				arg_earthradiation->e2j = &rot_trs2crs;
				if (!_gforce.earthradi->calculate(arg_earthradiation))
				{
					writeLogInfo(_glog, 4, function_name, "ERROR : can not calculate earthradiation");
					delete arg_earthradiation; arg_earthradiation = NULL;
					return false;
				}
				delete arg_earthradiation; arg_earthradiation = NULL;
				_acc_earthradiation = _gforce.earthradi->acc()*1e-3;
			}

			/*cout << "dt" << endl;
			cout << setw(25) << scientific << setprecision(16) << dt;*/

			if (_force_valid[CUSTOMER] == true)
			{
				t_gfcalarg_customer* arg_customer = new t_gfcalarg_customer();
				arg_customer->arg_type = ARG_CUSTOMER;
				arg_customer->time = &dt;
				arg_customer->sat_pos = &satXJ;
				arg_customer->sat_vel = &satVJ;
				arg_customer->period = &_period;
				arg_customer->cust_pars = &cust_par;
				arg_customer->scf2crs = &rot_scf2crsleo;
				if (!_gforce.customer->calculate(arg_customer))
				{
					writeLogInfo(_glog, 4, function_name, "ERROR : can not calculate customer");
					delete arg_customer; arg_customer = NULL;
					return false;
				}
				delete arg_customer; arg_customer = NULL;
				_acc_customer = _gforce.customer->acc();
			}

		/*	_gforce.point_mass->calculate(arg_point_mass);
			_gforce.gravity->calculate(arg_gravity);
			_gforce.relativity->calculate(arg_relativity);
			_gforce.solar_radi->calculate(arg_solar_boxwing);
			_gforce.atmosphere->calculate(arg_atmosphere);
			_gforce.customer->calculate(arg_customer);*/

			/*delete arg_gravity; arg_gravity = NULL;
			delete arg_point_mass; arg_point_mass = NULL;
			delete arg_solar_boxwing; arg_solar_boxwing = NULL;
			delete arg_atmosphere; arg_atmosphere = NULL;
			delete arg_customer; arg_customer = NULL;
			delete arg_relativity; arg_relativity = NULL;*/

		
			
			ColumnVector _acc_all =
				(_acc_point_mass
				+ _acc_gravity
				+ _acc_solar_radi
				+ _acc_customer
				+ _acc_relativity
				+ _acc_atmosphere
				+ _acc_earthradiation);

#ifdef _DEBUG_LOG_1
			cout << "gravity acc:" << endl;
			for (int i = 1; i <= 3; i++)
			{
				cout << setw(25) << scientific << setprecision(16) << _acc_gravity(i) << "  ";
			}
			cout << endl;

			cout << "pointmass acc:" << endl;
			for (int i = 1; i <= 3; i++)
			{
				cout << setw(25) << scientific << setprecision(16) << _acc_point_mass(i) << "  ";
			}
			cout << endl;

			cout << "relativity acc:" << endl;
			for (int i = 1; i <= 3; i++)
			{
				cout << setw(25) << scientific << setprecision(16) << _acc_relativity(i) << "  ";
			}
			cout << endl;
			cout << "solar acc:" << endl;
			for (int i = 1; i <= 3; i++)
			{
				cout << setw(25) << scientific << setprecision(16) << _acc_solar_radi(i) << "  ";
			}
			cout << endl;
			cout << "customer acc:" << endl;
			for (int i = 1; i <= 3; i++)
			{
				cout << setw(25) << scientific << setprecision(16) << _acc_customer(i) << "  ";
			}
			cout << endl;
			cout << "atmosphere acc:" << endl;
			for (int i = 1; i <= 3; i++)
			{
				cout << setw(25) << scientific << setprecision(16) << _acc_atmosphere(i) << "  ";
			}
			cout << endl;
			cout << "acc_all is:" << endl;
			for (int i = 1; i <= 3; i++)
			{
				cout << setw(25) << scientific << setprecision(16) << _acc_all(i) << "  ";
			}
			cout << endl;
#endif
			
			// calculate partial
			Matrix _partial_point_mass = _gforce.point_mass->partial();
			Matrix _partial_gravity = _gforce.gravity->partial();
			Matrix _partial_solar_radi = _gforce.solar_radi->partial();
			Matrix _partial_atmosphere = _gforce.atmosphere->partial();
			Matrix _partial_customer = _gforce.customer->partial();

			Matrix _partial_pos = _partial_gravity + _partial_point_mass;
			Matrix _partial_force = _partial_atmosphere.t()&_partial_customer;
			if (_srp_par.size() > 0)
				_partial_force = _partial_solar_radi.t()&_partial_force;
			Matrix _partial_vel(0.0, 3, 3);

#ifdef _DEBUG_LOG
			double dsod = t.sod() + t.dsec();
			cout << "EPO" << setw(8) << t.mjd() 
				 << setw(12) << fixed << setprecision(5) << dsod << endl;
			cout << " amat,cmat" << endl;
			for (int i = 1; i <= _partial_pos.Nrows(); i++) {
				for (int j = 1; j <= _partial_pos.Ncols(); j++)
					cout << scientific << uppercase << setw(24) << setprecision(14) << _partial_pos(j, i);
				cout << endl;
			}

			//cout << "vel is" << endl;
			//for (int i = 1; i <= _partial_vel.Nrows(); i++)
			//{
			//	for (int j = 1; j <= _partial_vel.Ncols(); j++)
			//	{
			//		cout << scientific << setw(25) << setprecision(16) << _partial_vel(j, i);
			//	}
			//	cout << endl;
			//}

			//cout << "force is" << endl;
			for (int i = 1; i <= _partial_force.Nrows(); i++) {
				for (int j = 1; j <= _partial_force.Ncols(); j++)
					cout << scientific << uppercase << setw(24) <<setprecision(14) << _partial_force(i, j);
				cout << endl;
			}
#endif

			// update partial 
			int index = 0;
			Matrix mat(0.0, 3, 1);
			mat = 0.0;
			vector<double> partial_tmp;
			for (int i = 1; i < 7; i++)
			{
				ColumnVector tmp(0.0, 3);
				index = 6 * i;
				tmp << state[index] << state[index + 1] << state[index + 2];
				mat = _partial_pos * tmp;
				tmp << state[index + 3] << state[index + 4] << state[index + 5];
				mat = mat + _partial_vel * tmp;
				for (int j = 1; j < 4; j++)
				{
					partial_tmp.push_back(mat(j, 1));
				}
			}

			mat = 0.0;
			for (int i = 1; i < _neq_size / 6 - 6; i++)
			{
				ColumnVector tmp(0.0, 3);
				index = 6 * i + 36;
				mat = _partial_force.Row(i).as_column();
				tmp << state[index] << state[index + 1] << state[index + 2];
				mat = mat + _partial_pos * tmp;
				tmp << state[index + 3] << state[index + 4] << state[index + 5];
				mat = mat + _partial_vel * tmp;
				for (int j = 1; j < 4; j++)
				{
					partial_tmp.push_back(mat(j, 1));
				}
			}
			//velocity
			partial[0] = state[3];
			partial[1] = state[4];
			partial[2] = state[5];

			//acc
			for (int i = 1; i <= _acc_all.size(); ++i)
			{
				partial[i + 2] = _acc_all(i);
			}
			//
			index = _acc_all.size() + 3;
			for (int i = 1; i < _neq_size / 6; i++)
			{
				partial[index] = state[i * 6 + 3];
				partial[index + 1] = state[i * 6 + 4];
				partial[index + 2] = state[i * 6 + 5];

				int j = (i - 1) * 3;
				partial[index + 3] = partial_tmp[j];
				partial[index + 4] = partial_tmp[j + 1];
				partial[index + 5] = partial_tmp[j + 2];

				index += 6;
			}
			
#ifdef _DEBUG_LOG
			cout << " funct" << endl;
			for (int i = 0; i < _neq_size / 6; ++i) {
				for (int j = 0; j < 6; ++j)
					cout << setw(24) << setprecision(14) << scientific << uppercase << right << state[i * 6 + j];
				cout << endl;
			}
			//cout << endl;
			//cout << "partial:" << endl;
			//for (int i = 0; i < _neq_size / 6; ++i) {
			//	for (int j = 0; j < 6; ++j)
			//		cout << setw(25) << setprecision(16) << scientific << right << partial[i * 6 + j];
			//	cout << endl;
			//}
			
#endif
			return true;
		}
		catch (...)
		{
			writeLogInfo(_glog, 2, function_name, "ERROR : unknown mistake");
			return false;
		}

	}

	/**
	* @brief  creat the force model for leo
	* @param[in]  force_ID   the force to creat
	* @param[in]  force_contain         end time from XML
	* @return  creat firce model correctly or not
	*/
	bool t_goi_leo::_creatForceModel(FORCE force_ID, t_gforce_leo * force_contain)
	{
		string function_name = "t_goi_leo::_creat_force_model";
		try
		{
			switch (force_ID)
			{
				// for gravity
			case great::GRAVITY:
			{
				if (_force_valid[force_ID] == true)
				{
				    t_init_gravity_egm08* gravity = NULL;
				
					if (_force_model[force_ID].find("EGM08_120") != string::npos)
					{
						gravity = new t_init_gravity_egm08(F_GRAVI_EGM08, 120, _gegmcoef, _gplanet, _gotl, _gotdata,_gdeadata,_ocean_degree,_isPoletide,_pole_degree);
						force_contain->creatForceModel(gravity, _glog);
					}
					else if (_force_model[force_ID].find("EGM08_12") != string::npos)
					{
						gravity = new t_init_gravity_egm08(F_GRAVI_EGM08, 12, _gegmcoef, _gplanet, _gotl, _gotdata);
						force_contain->creatForceModel(gravity, _glog);
					}
					else if (_force_model[force_ID].find("EGM08_150") != string::npos)
					{
						gravity = new t_init_gravity_egm08(F_GRAVI_EGM08, 150, _gegmcoef, _gplanet, _gotl, _gotdata);
						force_contain->creatForceModel(gravity, _glog);
					}
					else
					{
						writeLogInfo(_glog, 1, "WARNING : no setting for gravity model, we used EGM08 08th-order as default");
						gravity = new t_init_gravity_egm08(F_GRAVI_DEF, 8, _gegmcoef, _gplanet, _gotl, _gotdata);
						force_contain->creatForceModel(gravity, _glog);
					}

					delete gravity;
					gravity = NULL;
					break;
				}
				
				
			}
			// for relavity
			case great::RELATIVITY:
			{
				if (_force_valid[force_ID] == true)
				{
					t_init_relativity* relativity = NULL;
					if (_force_model[force_ID].find("Default") != string::npos)
					{
						relativity = new t_init_relativity(F_RELAV_DEF);
						force_contain->creatForceModel(relativity, _glog);
					}
					else
					{
						writeLogInfo(_glog, 1, "WARNING : no setting for relativity model, we used default model");
						relativity = new t_init_relativity(F_RELAV_DEF);
						force_contain->creatForceModel(relativity, _glog);
					}
					delete relativity;
					relativity = NULL;
					break;
				}
			}
			case great::POINT_MASS:
			{
				if (_force_valid[force_ID] == true)
				{
					t_init_point_mass* point_mass = NULL;
					if (_force_model[force_ID].find("Default") != string::npos)
					{
						point_mass = new t_init_point_mass(F_NBODY_DEF, _gplanet);
						force_contain->creatForceModel(point_mass, _glog);
					}
					else
					{
						writeLogInfo(_glog, 1, "WARNING : no setting for point_mass model, we used default model");
						point_mass = new t_init_point_mass(F_NBODY_DEF, _gplanet);
						force_contain->creatForceModel(point_mass, _glog);
					}
					delete point_mass;
					point_mass = NULL;
					break;
				}
			}
			case great::SOLAR_RADI:
			{
				if (_force_valid[force_ID] == true)
				{
					t_init_solar_boxwing* boxwing = NULL;
					if (_force_model[force_ID].find("BOX-WING") != string::npos)
					{
						boxwing = new t_init_solar_boxwing(F_SOLAR_BOXWING, _srp_par.size() != 0, _gpannel, _sat_mass);
						force_contain->creatForceModel(boxwing, _glog);
					}
					else
					{
						writeLogInfo(_glog, 1, "WARNING : no setting for solar radiation model, we used boxwing model");
						boxwing = new t_init_solar_boxwing(F_SOLAR_BOXWING, _srp_par.size() != 0, _gpannel, _sat_mass);
						force_contain->creatForceModel(boxwing, _glog);
					}
					delete boxwing;
					boxwing = NULL;
					break;
				}
			}
			
			case great::CUSTOMER:
			{
				if (_force_valid[force_ID] == true)
				{
					t_init_customer* customer = NULL;
					if (_force_model[force_ID].find("CUSTOMER") != string::npos)
					{
						customer = new t_init_customer(F_CUSTOMER, _sat_mass, true, _period);
						force_contain->creatForceModel(customer, _glog);
					}
					else
					{
						writeLogInfo(_glog, 1, "WARNING : no setting for customer model");
						customer = new t_init_customer(F_CUSTOMER, _sat_mass, true, _period);
						force_contain->creatForceModel(customer, _glog);
					}
					delete customer;
					customer = NULL;
					break;
				}
			}
			case great::ATMOSPHERE:
			{
				if (_force_valid[force_ID] == true)
				{
					t_init_atmosphere* atmosphere = NULL;
					if (_force_model[force_ID].find("DTM94") != string::npos)
					{
						atmosphere = new t_init_atmosphere(F_ATMOSPHERE, _sat_mass, "DTM94", _gsolargeomag, _gpannel);
						force_contain->creatForceModel(atmosphere, _glog);
					}
					else if (_force_model[force_ID].find("MSISE00") != string::npos)
					{
						atmosphere = new t_init_atmosphere(F_ATMOSPHERE, _sat_mass, "MSISE00", _gsolargeomag, _gpannel);
						force_contain->creatForceModel(atmosphere, _glog);
					}
					else if (_force_model[force_ID].find("DTM2013") != string::npos)
					{
						atmosphere = new t_init_atmosphere(F_ATMOSPHERE, _sat_mass, "DTM2013", _gsolargeomag, _gpannel);
						force_contain->creatForceModel(atmosphere, _glog);
					}
					else
					{
						writeLogInfo(_glog, 1, "WARNING : no setting for atmosphere model,we use atmosphere with dtm94");
						atmosphere = new t_init_atmosphere(F_ATMOSPHERE, _sat_mass, "DTM94", _gsolargeomag, _gpannel);
						force_contain->creatForceModel(atmosphere, _glog);
					}
					delete atmosphere;
					atmosphere = NULL;
					break;
				}
			}
			case great::EARTH_RADI:
			{
				if (_force_valid[force_ID] == true)
				{
					t_init_earthradiation* earthradi = NULL;

					if (_force_model[force_ID].find("CERES") != string::npos)
					{
						earthradi = new t_init_earthradiation(F_EARTHRADIATION, _sat_mass, "CERES", _galbedo, _gpannel);
						force_contain->creatForceModel(earthradi, _glog);
					}
					else
					{
						writeLogInfo(_glog, 1, "WARNING : no setting for earth radiation model,we use CERES");
						earthradi = new t_init_earthradiation(F_EARTHRADIATION, _sat_mass, "CERES", _galbedo, _gpannel);
						force_contain->creatForceModel(earthradi, _glog);
					}
					delete earthradi;
					earthradi = NULL;
					break;
				}
			}
			//case great::EARTH_RADI:
			//	break;
			default:
				break;
			}

			return true;
		}
		catch (...)
		{
			writeLogInfo(_glog, 0, function_name, "ERROR : throw exception while creating force model by the setting.");
			return false;
		}
	}
	/**
	* @brief  batch processing for LEO OI
	* @param[in]  beg         begin time from XML
	* @param[in]  end         end time from XML
	* @return  calculate correctly or not
	*/
	bool  t_goi_leo::processBatch(const t_gtime& beg, const t_gtime& end)
	{
		gtrace("t_gpppflt::processBatch");
		string function_name = "t_goi_leo::processBatch";
		_glog->comment(0, function_name, _sat + "processBatch started !");
		try
		{

			// check all data
			if (!_checkDataUsed())
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : no compelete data, check your xml file.");
				return false;
			}

			// get all the setting information
			if (!_get_settings())
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not get the settings form xml, check your xml file.");
				return false;
			}

			// prepare the trs2crs and crs2scf and so on;
			_trs2crs_cord = new t_gtrs2crs(false, _gpoleut1);
			_scf2crsleo_cord = new t_gscf2crsleo(_sat);

			// get the ics for the sat
			_sat_ics = _gics->get_leo_ics(_sat);
			
			if (!_check_par())
			{
				return false;
			}

			if (!_sat_ics)
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not get the ics data of " + _sat);
				return false;
			}

			// init sat info
			_initSatInfo(_sat);
			if (_gsat.prn().empty())
			{
				writeLogInfo(_glog, 1, function_name, "ERROR : no information for " + _sat + ", check your sat_paramters file.");
				return false;
			}

			// init all the force models
			// gravity model
			if (!_creatForceModel(GRAVITY, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(GRAVITY));
				return false;
			}
			// N-body model
			if (!_creatForceModel(POINT_MASS, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(POINT_MASS));
				return false;
			}
			// solar radiation model
			if (!_creatForceModel(SOLAR_RADI, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(SOLAR_RADI));
				return false;
			}
			//relativity
			if (!_creatForceModel(RELATIVITY, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(RELATIVITY));
				return false;
			}
			//customer
			if (!_creatForceModel(CUSTOMER, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(CUSTOMER));
				return false;
			};
			//atmosphere
			if (!_creatForceModel(ATMOSPHERE, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(ATMOSPHERE));
				return false;
			}
			//earth radiation
			if (!_creatForceModel(EARTH_RADI, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(EARTH_RADI));
				return false;
			}

			// extend the time with the settings
			if (!_extend_time_period(beg, end, _step_out, _extend_beg, _extend_end))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not extend time period ==>  beg : "
					+ dbl2str(beg.dmjd()) + " end : " + dbl2str(end.dmjd()));
				return false;
			}

			//check the time whether in the poleut1 data time range 
			if (_extend_beg.mjd() < _gpoleut1->getBegTime() || _extend_end.mjd() > _gpoleut1->getEndTime())
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : input time is beyond poleut1 data's time range");
				return false;
			}

			// get the number of equation and init the partial
			if (_sat_ics->get_par_num() == 0)
			{
				writeLogInfo(_glog, 1, function_name, "ERROR : no data in the " + _sat + " ics file, check your ics file.");
				return false;
			}
			else
			{
				_neq_size = 6 * (_sat_ics->get_par_num() + 1);
			}

			
			// current time must between the start and the end
			// reference time from ics
			t_gtime integra_ref_time = _gics->get_beg_time();
			_ref = integra_ref_time;
			if ((_extend_beg >= integra_ref_time) || (_extend_end <= integra_ref_time))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : There are some conflicts between the xml and ics file in beg and end time");
				return false;
			}

			// get init state and partial form ics
			vector<double> state_init;
			if (!_initStateMat(_neq_size, state_init))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not init the state matrix.");
				return false;
			}

			// get the integration direction and get the statr and end time of integration
			if (!_numerical_integration(integra_ref_time, _extend_beg, state_init, _state_frwd))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : something goes wrong while intergrating the time during "
					+ dbl2str(integra_ref_time.dmjd()) + " " + dbl2str(_extend_beg.dmjd()));
				return false;
			};
			if (!_numerical_integration(integra_ref_time, _extend_end, state_init, _state_back))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : something goes wrong while intergrating the time during "
					+ dbl2str(integra_ref_time.dmjd()) + " " + dbl2str(_extend_end.dmjd()));
				return false;
			};


			//write state info to orb file
			//init state 
			vector<double> init_funct(_neq_size / 2, 0.0);
			for (int i = 0; i < _neq_size / 6; i++)
			{
				init_funct[i * 3] = state_init[i * 6];
				init_funct[i * 3 + 1] = state_init[i * 6 + 1];
				init_funct[i * 3 + 2] = state_init[i * 6 + 2];

			}

			//save data into gorb
			_gorb->add_funct(0, _sat, init_funct);
			map<int, vector<double> >::iterator itState = _state_frwd.begin();
			for (itState; itState != _state_frwd.end(); ++itState)
			{
				int epo = itState->first;
				vector<double> one_epo_funct = itState->second;
				_gorb->add_funct(epo, _sat, one_epo_funct);
			}
			for (itState = _state_back.begin(); itState != _state_back.end(); ++itState)
			{
				int epo = itState->first;
				vector<double> one_epo_funct = itState->second;
				_gorb->add_funct(epo, _sat, one_epo_funct);
			}
			//set orb head info for only one time

			//set orb partial parlist
			map<ICS_PAR, vector<t_icsleo_record>>* ics_par = _sat_ics->leoics();
			vector<string > par_list;
			for (auto itpar=ics_par->begin(); itpar != ics_par->end(); itpar++)
			{
				
				par_list.push_back(ics2str(itpar->first));
				
				
			}
			_gorb->add_parlist(_sat, par_list);

			static bool isfirst_orbhd = true;
			if (isfirst_orbhd)
			{
				isfirst_orbhd = false;
				t_orb_hd orb_hd;
				orb_hd.adminp = _step_inp;
				orb_hd.admout = _step_out;
				orb_hd.admord = _integra_degree;
				orb_hd.nequ = _neq_size / 2;
				orb_hd.nepo = (_state_frwd.size() + _state_back.size() + 1)*_step_inp / _step_out;
				orb_hd.beg = _extend_beg;
				orb_hd.end = _extend_end;
				orb_hd.t_ref = integra_ref_time;
				orb_hd.partial = true;
				_gorb->add_head(orb_hd);
			}

			delete _trs2crs_cord; _trs2crs_cord = NULL;
			delete _scf2crsleo_cord; _scf2crsleo_cord = NULL;

			return true;

			//	// init all the force models
			//	t_init_gravity_egm08* gravity = new t_init_gravity_egm08();
			//	gravity->gracoef = _gegmcoef;
			//	gravity->gotdata = _gotdata;
			//	gravity->gotl = _gotl;
			//	gravity->gplanet = _gplanet;
			//	gravity->force_model = F_GRAVI_EGM08;
			//	gravity->gorder = 150;
			//	_gforce.creat_force_model(gravity, _glog);
			//	delete gravity;	      gravity = NULL;

			//	t_init_point_mass* point_mass = new t_init_point_mass();
			//	point_mass->force_model = F_NBODY_DEF;
			//	point_mass->gplanet = _gplanet;
			//	_gforce.creat_force_model(point_mass, _glog);
			//	delete point_mass;    point_mass = NULL;

			//	map<ICS_PAR, double> solarpars;
			//	t_init_solar_boxwing* boxwing = new t_init_solar_boxwing();
			//	boxwing->ispartial = false;
			//	boxwing->mass = _satmass;
			//	boxwing->pannel = _gpannel;
			//	boxwing->force_model = F_SOLAR_BOXWING;
			//	_gforce.creat_force_model(boxwing, _glog);
			//	delete boxwing;         boxwing = NULL;

			//	t_init_relativity* relativity = new t_init_relativity();
			//	relativity->force_model = F_RELAV_DEF;
			//	_gforce.creat_force_model(relativity, _glog);
			//	delete relativity;	  relativity = NULL;

			//	t_init_atmosphere* atmosphere = new t_init_atmosphere();
			//	atmosphere->mass = _satmass;
			//	atmosphere->model = "DTM94";
			//	atmosphere->gsolmag = _gsolargeomag;
			//	atmosphere->gpannel = _gpannel;
			//	atmosphere->force_model = F_ATMOSPHERE;
			//	_gforce.creat_force_model(atmosphere, _glog);
			//	delete atmosphere; atmosphere = NULL;

			//	t_init_customer* customer = new t_init_customer();
			//	customer->mass = _satmass;
			//	customer->ispartial = false;
			//	customer->period = _period;
			//	customer->force_model = F_CUSTOMER;
			//	_gforce.creat_force_model(customer, _glog);
			//	delete customer; customer = NULL;

			//	// reference time from ics
			//	t_gtime integra_ref_time = _gics->beg();
			//	_REF = integra_ref_time;

			//	// current time must between the start and the end
			//	bool ics_valid = (_extend_beg <= integra_ref_time)
			//		&& (_extend_end >= integra_ref_time);
			//	if (!ics_valid)
			//	{
			//		cout << "There are some conflicts between the xml and ics file in beg and end time" << endl;
			//		return 0;
			//	}

			//	// get init state and partial form ics
			//	vector<double> state_init;
			//	_init_state_matrix(_neq_size, state_init);

			//	// get the integration direction and get the statr and end time of integration
			//	_numerical_integration(integra_ref_time, _extend_beg, state_init, _state_for);
			//	_numerical_integration(integra_ref_time, _extend_end, state_init, _state_back);


			//	//write state info to orb file
			//	//init state 
			//	vector<double> init_funct(_neq_size / 2, 0.0);
			//	for (int i = 0; i < _neq_size / 6; i++)
			//	{
			//		init_funct[i * 3] = state_init[i * 6];
			//		init_funct[i * 3 + 1] = state_init[i * 6 + 1];
			//		init_funct[i * 3 + 2] = state_init[i * 6 + 2];

			//	}
			//	//save data into gorb
			//	_gorb->add_funct(0, _sat, init_funct);
			//	map<int, vector<double> >::iterator itState = _state_for.begin();
			//	for (itState; itState != _state_for.end(); ++itState)
			//	{
			//		int epo = itState->first;
			//		vector<double> one_epo_funct = itState->second;
			//		_gorb->add_funct(epo, _sat, one_epo_funct);
			//	}
			//	for (itState = _state_back.begin(); itState != _state_back.end(); ++itState)
			//	{
			//		int epo = itState->first;
			//		vector<double> one_epo_funct = itState->second;
			//		_gorb->add_funct(epo, _sat, one_epo_funct);
			//	}
			//	//set orb head info for only one time
			//	static bool isfirst_orbhd = true;
			//	if (isfirst_orbhd)
			//	{
			//		isfirst_orbhd = false;
			//		t_orb_hd orb_hd;
			//		orb_hd.adminp = _step_inp;
			//		orb_hd.admout = _step_out;
			//		orb_hd.admord = _integra_degree;
			//		orb_hd.nequ = _neq_size / 2;
			//		orb_hd.nepo = (_state_for.size() + _state_back.size() + 1)*_step_inp / _step_out;
			//		orb_hd.beg = _extend_beg;
			//		orb_hd.end = _extend_end;
			//		orb_hd.t_ref = integra_ref_time;
			//		orb_hd.partial = true;
			//		_gorb->add_head(orb_hd);
			//	}
			//	delete _trs2crs_cord; _trs2crs_cord = NULL;
			//	delete _scf2crsleo_cord; _scf2crsleo_cord = NULL;

			//	return 0;
			//}

		}
		catch (...)
		{
			if (_trs2crs_cord)
			{
				delete _trs2crs_cord; _trs2crs_cord = NULL;
			}
			if (_scf2crsleo_cord)
			{
				delete _scf2crsleo_cord; _scf2crsleo_cord = NULL;
			}

			writeLogInfo(_glog, 1, function_name, "ERROR : throw exception while doing procseeBatch of" + _sat);
			return false;
		}
	}
}