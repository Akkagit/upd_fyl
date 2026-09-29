/**
*
* @verbatim
	History
	 -1.0 Jiande Huang   2019-01-04  creat the file
	 -1.1 Jie Li		 2019-04-01  add code in doxygen style.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    goi_nav.cpp
* @brief		orbit integration processing for GNS satellite
*
* @author       Jiande Huang, Wuhan University
* @version		1.0.0
* @date		    2019-01-04
*
*/



#include "gproc/goi_nav.h"
#include "gforce/gforce.h"
#include "gforce/gfcommon.h"
#include "gset/gsetforce.h"
#include "gutils/ginfolog.h"
#include "newmat/newmatio.h"

namespace great
{

	/** @brief  constructor.
	* @param[in]	mark   satellite prn
	* @param[in]	set    setting information
	*/
	t_goi_nav::t_goi_nav(string mark, t_gsetbase* set)
	:t_goi(mark,set)
	{
		gtrace("t_goi_nav::t_goi_nav");

		_sat = mark;
		_set = set;

		// the last column of _partial_all and _state_all
		//_gorb = new t_gorb;
	}

	t_goi_nav::~t_goi_nav()
	{

	}

	

	/**
	* @brief set  ICS data
	* @param[in] gics    ICS data
	*/
	void t_goi_nav::set_data_used(t_gallics_nav* gics)
	{
		_gics = gics;
	}

	
	/**
	* @brief   get orb funct according to the satellite name
	* @param[in]  sat   satellite name
	* @param[out] gorb  orb data
	*/
	/*void t_goi_nav::get_orb_funct(string sat, t_gorb* gorb)
	{
		map<int, vector< t_gtriple> > funct_tmp;
		_gorb->get_funct(sat, funct_tmp);

		gorb->add_funct(sat, funct_tmp);
	}*/

	/**
	* @brief  get orb header data according to the satellite name
	* @param[in]  sys   system of the satelite
	* @param[in]  sats  all satellite names
	* @param[out] gorb  orb data
	*/
	//void t_goi_nav::get_orb_head(set<string>& sys, set<string>& sats, t_gorb* gorb)
	//{
	//	t_orb_hd orb_hd_tmp;
	//	_gorb->get_head(orb_hd_tmp);
	//	//sat group
	//	if (sys.size() > 1)
	//	{
	//		orb_hd_tmp.satgroup = "GNS";
	//	}
	//	else
	//	{
	//		orb_hd_tmp.satgroup = *(sys.begin());
	//	}
	//	//sat type
	//	set<string>::iterator itsys = sys.begin();
	//	for (itsys; itsys != sys.end(); ++itsys)
	//	{
	//		char sys_tmp = t_gsys::str2char(*itsys);
	//		orb_hd_tmp.sattype.insert(sys_tmp);
	//	}
	//	//sat num
	//	orb_hd_tmp.nsat = sats.size();
	//	//sat list
	//	orb_hd_tmp.satprn = sats;
	//	gorb->add_head(orb_hd_tmp);
	//}

	/**
	* @brief   init the sat infomation
	* @param[in] sat_name   satellite name
	*/
	bool t_goi_nav::_init_sat_info(string sat_name)
	{
		string function_name = "t_goi_nav::_init_sat_info";
		try
		{
			// 可能会有问题，是否需要扩展时间不确定
			t_gtime beg = _extend_beg + NEXT_POINT_NUM * _step_out;
			t_gtime end = _extend_end - NEXT_POINT_NUM * _step_out;
			if (!_gallsat->getSatinfoUsed(sat_name, beg, end, _gsat))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not get the sat information of "
					+ sat_name + ", check your sat parameter file or the read process of the file.");
				return false;
			}
			_gnsmass = _gsat.mass();
			return true;
		}
		catch (...)
		{
			writeLogInfo(_glog, 3, function_name, "ERROR : throw while getting the sat information");
			return false;
		}
	}

	/**
	* @brief      init the state matrix
	* @param[in]  neq_size   number of OI funct params(（np+6）*6+6)
	* @param[out] state      state data of OI processing
	*/
	bool t_goi_nav::_init_state_matrix(int neq_size, vector<double>& state)
	{
		string function_name = "t_goi_nav::_init_state_matrix";
		try
		{
			// x,y,z,vx,vy,vz and solar pars will be read from ICS
			// 从ics文件获取x，y，z和vx，vy，vz状态矩阵的初值
			state.push_back(_sat_ics->get_par_value(ICS_PAR::PXSAT));
			state.push_back(_sat_ics->get_par_value(ICS_PAR::PYSAT));
			state.push_back(_sat_ics->get_par_value(ICS_PAR::PZSAT));
			state.push_back(_sat_ics->get_par_value(ICS_PAR::VXSAT));
			state.push_back(_sat_ics->get_par_value(ICS_PAR::VYSAT));
			state.push_back(_sat_ics->get_par_value(ICS_PAR::VZSAT));

			// initial partial can be ensured derectly
			// storaging with column vector is easy to integral
			// 偏导状态可以直接确定
			// 因为将矩阵存储为了列向量形式（为了便于积分），所以采用neq_size / 6循环
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
			for (int i = 0; i < state.size; i++)
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

	//// 初始化偏导矩阵
	//void t_goi_nav::_init_partial_matrix(int par_num, Matrix & partial)
	//{
	//	partial.resize(par_num + 1, 6);
	//	partial = 0.0;
	//}

	/**
	* @brief   Get settings from XML file and set local variables.
	* @return  get settings correctly or not
	*/
	bool t_goi_nav::_get_settings()
	{
		string function_name = "t_goi_nav::_get_settings";
		try
		{
			// sys and sat name
			_sys = t_gsys::str2gsys(_sat.substr(0, 1));

			// settings for intergrator
			_step_inp = dynamic_cast<t_gsetforce*>(_set)->getInpStep(_sat);
			_step_out = dynamic_cast<t_gsetforce*>(_set)->getOutStep(_sat);
			_integra_degree = dynamic_cast<t_gsetforce*>(_set)->getIntegraDegree(_sat);

			// settings for forces
			// List of models for considered forces
			_forces = dynamic_cast<t_gsetforce*>(_set)->getForceUsedSet(_sat);
			_force_model = dynamic_cast<t_gsetforce*>(_set)->getForceModelMap(_sat);
			_force_valid = dynamic_cast<t_gsetforce*>(_set)->getForcesValidMap(_sat);
			/*t_gsatparam onesatinfo;
			t_gtime beg;
			t_gtime end;
			beg = _extend_beg + NEXT_POINT_NUM * _step_out;
			end = _extend_end - NEXT_POINT_NUM * _step_out;
			if (!_gallsat->getSatinfoUsed(_sat, beg, end, onesatinfo)) {
				writeLogInfo(_glog, 0, function_name, "ERROR : can not get the sat information of "
					+ _sat + ", check your sat parameter file or the read process of the file.");
				return false;
			}
			_gnsmass = onesatinfo.mass();*/


#ifdef DEBUG
			// Creat and set the log file : oi.log
			cout << _sat + "------------------------------------------------" << endl
				<< _sat + "step_inp       : " + dbl2str(_step_inp) << endl
				<< _sat + "step_out       : " + dbl2str(_step_out) << endl
				<< _sat + "integra_degree : " + dbl2str(_integra_degree) << endl
				<< _sat + "------------------------------------------------" << endl;

			set<int>::iterator it;
			for (it = _forces.begin(); it != _forces.end(); it++)
			{
				cout << _sat + " force model : " + _force_models[*it] + " is " << _forces_valid[*it] << endl;
			}
#endif
			return true;
		}
		catch (...)
		{
			writeLogInfo(_glog, 0, function_name, "ERROR : throw exception while get settings form xml file");
			return false;
		}

	}

	bool t_goi_nav::_check_data()
	{
		if (!_gplanet || _gplanet->is_empty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no de data");
			return false;
		};			///< data for all planets.

		if (!_gpoleut1 || _gpoleut1->isEmpty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no poleut1 data");
			return false;
		};			///< poleut1 data.

		if (!_gotdata || _gotdata->is_empty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no ocean tide data");
			return false;
		};			///< ocean tide data.
		
		if (!_gleapsecond || _gleapsecond->is_empty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no leapsecond data");
			return false;
		};		///< leap second data.
		
		if (!_gics || _gics->isEmpty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no ics data");
			return false;
		};				///< all  ics data.
		
		if (!_gegmcoef || _gegmcoef->is_empty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no EGM data");
			return false;
		};			///< EGM data.
		
		if (!_gallsat || _gallsat->isEmpty())
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no sat parameters data");
			return false;
		};			///< satellite parameters data.
		
		if (!_gotl)
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no ocean load data");
			return false;
		};				///< ocean load data.

		if (_force_model[(SOLAR_RADI)].find("PBW") != string::npos)
		{
			if (!_gpannel || _gpannel->isEmpty())
			{
				writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no pannel data");
					return false;
			};
		}
		///< _gpannel data when using PBW models.

		/*if (_force_models[(SOLAR_RADI)].find("PBW") != string::npos)
		{
			writeLogInfo(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no pannel data");
			return false;
		};		*/		///< _gnsmass when using PBW models.
		

		
		
		// _sat_ics data not be assigned
		/*if (_sat_ics->is_empty())
		{
			write_log(_glog, 0, "ERROR : t_goi_nav::_check_data ==> no ics data");
			return false;
		};		*/	///< one satellite ics data.
		return true;
	}

	/** @brief Set output file (empty). */
	void t_goi_nav::_setOut()
	{
	}

	/** @brief Set processing log file (<ppp>). */
	void t_goi_nav::_setLog()
	{
	}

	


	
	/**
	* @brief  calculate the acc and partial in epoch t
	* @param[in]  t			current epoch
	* @param[in]  state     state at current epoch
	* @param[out] partial   partial of the state
	* @return  calculate correctly or not
	*/
	bool t_goi_nav::_calculate(const t_gtime& t, const vector<double>& state, vector<double>& partial)
	{
		//change time sys to tdt
		t_gtime tdt_gtime = t;
		tdt_gtime.tsys(t_gtime::TT);
		//ChangeTimeSys(GPS2TDT, t, tdt_gtime);
		double tdt = tdt_gtime.dmjd();

		//double tdt;
		//change_time_sys(GPS2TDT, t, tdt);

		_trs2crs_cord->calcRotMat(tdt_gtime, true, true, true);

		// Matrix
		ColumnVector satXJ(0.0, 3);		// satellite position in J2000 卫星在J2000的坐标
		ColumnVector satVJ(0.0, 3);		// satellite VELOCITY in J2000 卫星在J2000的速度
		ColumnVector satXE(0.0, 3);		// satellite position in TRS   卫星在TRS的位置

		map<string, ColumnVector> allXE;	// planet position in TRS,used for ocean tide correcting 行星天体在TRS坐标下的位置，用于海潮改正等
		map<string, ColumnVector> allXJ;	// planet position in J2000,used for point_mass and so on 行星天体在J2000坐标下的位置，用于N体引力计算等

		Matrix rot_scf2crs(0.0, 3, 3);	// rotation matrix of SCF to CRS 星固系到惯性系的旋转矩阵
		Matrix rot_crs2acr(0.0, 3, 3);	// rotation matrix of CRS to ARC 惯性系到ARC坐标系的旋转矩阵
		Matrix rot_trs2crs(0.0, 3, 3);	// rotation matrix of TRS to CRS 地固系到惯性系的旋转矩阵

		// 获取地固系到惯性系的坐标旋转矩阵等坐标转换信息
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
			satVJ(i + 1) = state[(int64_t)i + 3];
		}


		ColumnVector x_earth = (satXJ.t() * _Dgmat).t();
		satXE = (satXJ.t() * rot_trs2crs).t();
		//satVE = (satVJ.t() * rot_trs2crs).t() + x_earth / 43200.0 * G_PI;

		// get the sun and moon pos in crs
		_gplanet->get_pos(tdt, "SUN", sun_pos_crs);
		_gplanet->get_pos(tdt, "MOON", moon_pos_crs);
		allXE["SUN"] = (sun_pos_crs.t() * rot_trs2crs).t();
		allXJ["SUN"] = sun_pos_crs;
		allXE["MOON"] = (moon_pos_crs.t() * rot_trs2crs).t();;

		double shadow_factor_ecom = get_shadow_factor(satXJ, sun_pos_crs, moon_pos_crs);

		// get scf2crs
		//_scf2crs_cord->calcRotMat(satXJ, satVJ, allXJ["SUN"], rot_scf2crs);
		//extend by yqyuan for attitude modeling
		_scf2crs_cord->calcRotMat(_gsat.blocktype(),_sat,satXJ, satVJ, allXJ["SUN"], rot_scf2crs);

		// If the box-wing model is used, the update the normdirj of _gpannel
		if (_force_model[(SOLAR_RADI)].find("PBW") != string::npos) {
			map<int, t_oneside>::iterator it_pannel = _gpannel->get_pannels()->begin();
			for (it_pannel; it_pannel != _gpannel->get_pannels()->end(); it_pannel++)
			{
				Matrix norm;
				norm = it_pannel->second.scf_norm_dir.crd_cvect();
				norm.resize_keep(3, 1);
				it_pannel->second.crs_norm_dir = rot_scf2crs*norm;
			}
		}

		// get crs2acr
		_crs2acr_cord->calcRotMat(satXJ, satVJ, rot_crs2acr);

		// calculate acc
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
			writeLogInfo(_glog, 4, "t_goi_nav::_calculate", "ERROR : can not calculate gravity");
			delete arg_gravity; arg_gravity = NULL;
			return false;
		};
		delete arg_gravity; arg_gravity = NULL;


		t_gfcalarg_point_mass* arg_point_mass = new t_gfcalarg_point_mass();
		arg_point_mass->arg_type = ARG_POINT_MASS;
		arg_point_mass->time = &tdt;
		arg_point_mass->sat_pos = &satXJ;
		arg_point_mass->sat_vel = &satVJ;
		if (!_gforce.point_mass->calculate(arg_point_mass))
		{
			writeLogInfo(_glog, 4, "t_goi_nav::_calculate", "ERROR : can not calculate point mass");
			delete arg_point_mass; arg_point_mass = NULL;
			return false;
		};
		delete arg_point_mass; arg_point_mass = NULL;


		t_gfcalarg_relativity* arg_relativity = new t_gfcalarg_relativity();
		arg_relativity->arg_type = ARG_RELATIVITY;
		arg_relativity->time = &tdt;
		arg_relativity->sat_pos = &satXJ;
		arg_relativity->sat_vel = &satVJ;
		if (!_gforce.relativity->calculate(arg_relativity))
		{
			writeLogInfo(_glog, 4, "t_goi_nav::_calculate", "ERROR : can not calculate relativity");
			delete arg_relativity; arg_relativity = NULL;
			return false;
		};
		delete arg_relativity; arg_relativity = NULL;
		if (_gforce.solar_radi->force_model() == F_SOLAR_ECOM_3 || 
			_gforce.solar_radi->force_model()== F_SOLAR_ECOM_5 ||
			_gforce.solar_radi->force_model() == F_SOLAR_ECOM_9) 
		{
			t_gfcalarg_solar_ecom* arg_solar_ecom = new t_gfcalarg_solar_ecom();
			arg_solar_ecom->arg_type = ARG_SOLAR_ECOM;
			arg_solar_ecom->sat_type = _gsat.blocktype();
			arg_solar_ecom->moon_pos = &moon_pos_crs;
			arg_solar_ecom->sun_pos = &sun_pos_crs;
			arg_solar_ecom->shadow_factor = &shadow_factor_ecom;
			arg_solar_ecom->time = &tdt;
			arg_solar_ecom->sat_pos = &satXJ;
			arg_solar_ecom->sat_vel = &satVJ;
			if (!_gforce.solar_radi->calculate(arg_solar_ecom))
			{
				writeLogInfo(_glog, 4, "t_goi_nav::_calculate", "ERROR : can not calculate solar_radi");
				delete arg_solar_ecom; arg_solar_ecom = NULL;
				return false;
			};
			delete arg_solar_ecom; arg_solar_ecom = NULL;
		}
		if (_gforce.solar_radi->force_model() == F_SOLAR_ECOM2)
		{
			t_gfcalarg_solar_ecom2* arg_solar_ecom2 = new t_gfcalarg_solar_ecom2();
			arg_solar_ecom2->arg_type = ARG_SOLAR_ECOM2;
			arg_solar_ecom2->sat_type = _gsat.blocktype();
			arg_solar_ecom2->moon_pos = &moon_pos_crs;
			arg_solar_ecom2->sun_pos = &sun_pos_crs;
			arg_solar_ecom2->shadow_factor = &shadow_factor_ecom;
			arg_solar_ecom2->time = &tdt;
			arg_solar_ecom2->sat_pos = &satXJ;
			arg_solar_ecom2->sat_vel = &satVJ;
			if (!_gforce.solar_radi->calculate(arg_solar_ecom2))
			{
				writeLogInfo(_glog, 4, "t_goi_nav::_calculate", "ERROR : can not calculate solar_radi");
				delete arg_solar_ecom2; arg_solar_ecom2 = NULL;
				return false;
			};
			delete arg_solar_ecom2; arg_solar_ecom2 = NULL;
		}
		if (_gforce.solar_radi->force_model() == F_SOLAR_PBW_5||
			_gforce.solar_radi->force_model() == F_SOLAR_PBW_9)
		{
			t_gfcalarg_priori_solar_boxwing* arg_solar_pbw = new t_gfcalarg_priori_solar_boxwing();
			arg_solar_pbw->arg_type = ARG_SOLAR_PBW;
			arg_solar_pbw->sat_type = _gsat.blocktype();
			arg_solar_pbw->sat_mass = _gsat.mass();
			arg_solar_pbw->moon_pos = &moon_pos_crs;
			arg_solar_pbw->sun_pos = &sun_pos_crs;
			arg_solar_pbw->shadow_factor = &shadow_factor_ecom;
			arg_solar_pbw->time = &tdt;
			arg_solar_pbw->sat_pos = &satXJ;
			arg_solar_pbw->sat_vel = &satVJ;
			if (!_gforce.solar_radi->calculate(arg_solar_pbw))
			{
				writeLogInfo(_glog, 4, "t_goi_nav::_calculate", "ERROR : can not calculate solar_radi");
				delete arg_solar_pbw; arg_solar_pbw = NULL;
				return false;
			};
			delete arg_solar_pbw; arg_solar_pbw = NULL;
		}
		ColumnVector _acc_earthradiation(0.0, 3);
		if (_force_valid[EARTH_RADI] == true)
		{
			t_gfcalarg_earthradiation* arg_earthradiation = new t_gfcalarg_earthradiation();
			arg_earthradiation->arg_type = ARG_EARTHRADIATION;
			arg_earthradiation->sat_mass = _gsat.mass();
			arg_earthradiation->sat_blocknum = _gsat.blocktype();
			arg_earthradiation->sat_pos = &satXJ;
			arg_earthradiation->sat_vel = &satVJ;
			arg_earthradiation->sun_pos_crs = &sun_pos_crs;
			arg_earthradiation->e2j = &rot_trs2crs;
			if (!_gforce.earthradi->calculate(arg_earthradiation))
			{
				writeLogInfo(_glog, 4, "t_goi_nav::_calculate", "ERROR : can not calculate earthradiation");
				delete arg_earthradiation; arg_earthradiation = NULL;
				return false;
			}
			_acc_earthradiation = _gforce.earthradi->acc() * 1e-3;
			delete arg_earthradiation; arg_earthradiation = NULL;
		}

		//todo &
		ColumnVector _acc_point_mass = _gforce.point_mass->acc();
		ColumnVector _acc_gravity = _gforce.gravity->acc();
		ColumnVector _acc_relativity = _gforce.relativity->acc();
		ColumnVector _acc_solar_radi = _gforce.solar_radi->acc();

		ColumnVector _acc_all =
			(_acc_point_mass
				+ _acc_gravity
				+ _acc_relativity
				+ _acc_solar_radi
				+ _acc_earthradiation);

		// calculate partial
		Matrix _partial_point_mass = _gforce.point_mass->partial();
		Matrix _partial_gravity = _gforce.gravity->partial();
		Matrix _partial_solar_radi = _gforce.solar_radi->partial();

		Matrix _partial_pos = _partial_gravity + _partial_point_mass;
		Matrix _partial_force = _partial_solar_radi;
		Matrix _partial_vel(0.0, 3, 3);

		// update partial 
		int index = 0;
		Matrix mat(0.0, 3, 1);
		mat = 0.0;
		vector<double> partial_tmp;
		for (int i = 1; i < 7; i++)
		{
			ColumnVector tmp(0.0, 3);
			index = 6 * i;
			tmp << state[(int64_t)index] << state[(int64_t)index + 1] << state[(int64_t)index + 2];
			mat = _partial_pos * tmp;
			tmp << state[(int64_t)index + 3] << state[(int64_t)index + 4] << state[(int64_t)index + 5];
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
			tmp << state[(int64_t)index] << state[(int64_t)index + 1] << state[(int64_t)index + 2];
			mat = mat + _partial_pos * tmp;
			tmp << state[(int64_t)index + 3] << state[(int64_t)index + 4] << state[(int64_t)index + 5];
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
			partial[index] = state[(uint64_t)i * 6 + 3];
			partial[(uint64_t)index + 1] = state[(uint64_t)i * 6 + 4];
			partial[(uint64_t)index + 2] = state[(uint64_t)i * 6 + 5];

			int j = (i - 1) * 3;
			partial[(uint64_t)index + 3] = partial_tmp[j];
			partial[(uint64_t)index + 4] = partial_tmp[(uint64_t)j + 1];
			partial[(uint64_t)index + 5] = partial_tmp[(uint64_t)j + 2];

			index += 6;
		}
		return true;
	}

	/**
	* @brief  creat initial information of force model
	* @param[in]  force_ID			force type
	* @param[in]  force_contain     force model initial information
	* @return  creat force model correctly or not
	*/
	bool t_goi_nav::_creat_force_model(FORCE force_ID, t_gforce_nav * force_contain)
	{
		string function_name = "t_goi_nav::_creat_force_model";
		try
		{
			switch (force_ID)
			{
				// for gravity
			case great::GRAVITY:
			{
				t_init_gravity_egm08* gravity = NULL;
				if (_force_model[force_ID].find("EGM08_12") != string::npos)
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
			// for relavity
			case great::RELATIVITY:
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
			case great::POINT_MASS:
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
			case great::SOLAR_RADI:
			{
				t_init_solar_ecom* ecom = NULL;
				t_init_solar_ecom2* ecom2 = NULL;
				t_init_solar_priori_boxwing* pbw = NULL;
				if (_force_model[force_ID].find("ECOM_5") != string::npos)
				{
					ecom = new t_init_solar_ecom(F_SOLAR_ECOM_5, _gics->get_nav_ics(_sat)->get_ecom5_solar_pars_values());
					force_contain->creatForceModel(ecom, _glog);
				}
				else if (_force_model[force_ID].find("ECOM_9") != string::npos)
				{
					ecom = new t_init_solar_ecom(F_SOLAR_ECOM_9, _gics->get_nav_ics(_sat)->get_ecom9_solar_pars_values());
					force_contain->creatForceModel(ecom, _glog);
				}
				else if (_force_model[force_ID].find("ECOM_3") != string::npos)
				{
					ecom = new t_init_solar_ecom(F_SOLAR_ECOM_3, _gics->get_nav_ics(_sat)->get_ecom3_solar_pars_values());
					force_contain->creatForceModel(ecom, _glog);
				}
				else if (_force_model[force_ID].find("ECOM2") != string::npos)
				{
					ecom2 = new t_init_solar_ecom2(F_SOLAR_ECOM2, _gics->get_nav_ics(_sat)->get_ecom2_solar_pars_values());
					force_contain->creatForceModel(ecom2, _glog);
				}
				else if (_force_model[force_ID].find("PBW_5") != string::npos)
				{
					pbw = new t_init_solar_priori_boxwing(F_SOLAR_PBW_5, _gpannel, _gics->get_nav_ics(_sat)->get_pbw_solar_pars_values());
					force_contain->creatForceModel(pbw, _glog);
				}
				else if (_force_model[force_ID].find("PBW_9") != string::npos)
				{
					pbw = new t_init_solar_priori_boxwing(F_SOLAR_PBW_9, _gpannel, _gics->get_nav_ics(_sat)->get_pbw_solar_pars_values());
					force_contain->creatForceModel(pbw, _glog);
				}
				else
				{
					writeLogInfo(_glog, 1, "WARNING : no setting for solar radiation model, we used ECOM with 5 pars as default model");
					ecom = new t_init_solar_ecom(F_SOLAR_ECOM_5, _gics->get_nav_ics(_sat)->get_ecom5_solar_pars_values());
					force_contain->creatForceModel(ecom, _glog);
				}
				delete ecom;
				delete ecom2;
				delete pbw;
				ecom = NULL;
				ecom2 = NULL;
				pbw = NULL;
				break;
			}
			case great::EARTH_RADI:
			{
				t_init_earthradiation2* earthradi = NULL;
				if (_force_valid[force_ID] == true)
				{
					if (_force_model[force_ID].find("CERES") != string::npos)
					{
						earthradi = new t_init_earthradiation2(F_EARTHRADIATION, "CERES", _galbedo);
						force_contain->creatForceModel(earthradi, _glog);
					}
					else
					{
						writeLogInfo(_glog, 1, "WARNING : no setting for earth radiation model,we use ERM-A");
						earthradi = new t_init_earthradiation2(F_EARTHRADIATION, "ERM_A", _galbedo);
						force_contain->creatForceModel(earthradi, _glog);
					}
				}
				delete earthradi;
				earthradi = NULL;
				break;
			}
			case great::CUSTOMER:
				break;
			case great::ATMOSPHERE:
				break;
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
	* @brief  batch processing for GNS OI
	* @param[in]  beg         begin time from XML
	* @param[in]  end         end time from XML
	* @return  OI correctly or not
	*/
	bool t_goi_nav::processBatch(const t_gtime& beg, const t_gtime& end)
	{
		gtrace("t_gpppflt::processBatch");
		string function_name = "t_goi_nav::processBatch";
		// jdhuang
		_glog->comment(1, function_name, _sat + " processBatch : started !");

		try
		{

			if (!_check_data())
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

			// get the ics for the sat
			_sat_ics = _gics->get_nav_ics(_sat);
			if (!_sat_ics)
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not get the ics data of " + _sat);
				return false;
			}

			// prepare the trs2crs and crs2scf and so on;
			_trs2crs_cord = new t_gtrs2crs(false, _gpoleut1);
			_scf2crs_cord = new t_gscf2crs();
			_crs2acr_cord = new t_gcrs2acr();

			// init all the force models
			// gravity model
			// 重力场模型
			if (!_creat_force_model(GRAVITY, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(GRAVITY));
				return false;
			};
			// N-body model
			// N体引力模型
			if (!_creat_force_model(POINT_MASS, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(POINT_MASS));
				return false;
			};
			// solar radiation model
			// 太阳光压模型
			if (!_creat_force_model(SOLAR_RADI, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(SOLAR_RADI));
				return false;
			};
			if (!_creat_force_model(EARTH_RADI, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(EARTH_RADI));
				return false;
			};
			//relativity
			// 相对论力
			if (!_creat_force_model(RELATIVITY, &_gforce))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not creat force model ==>" + force2str(RELATIVITY));
				return false;
			};

			// extend the time with the settings
			if (!_extend_time_period(beg, end, _step_out, _extend_beg, _extend_end))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not extend time period ==>  beg : "
					+ dbl2str(beg.dmjd()) + " end : " + dbl2str(end.dmjd()));
				return false;
			};

			//check the time whether in the poleut1 data time range 
			if (_extend_beg.mjd() < _gpoleut1->getBegTime() || _extend_end.mjd() > _gpoleut1->getEndTime())
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : input time is beyond poleut1 data's time range");
				_glog->comment(t_glog::LOG_LV::LOG_ERROR, "t_goi_nav", "processBatch", _extend_beg.str_mjdsod("extend beg"));
				_glog->comment(t_glog::LOG_LV::LOG_ERROR, "t_goi_nav", "processBatch", _extend_end.str_mjdsod("extend end"));
				_glog->comment(t_glog::LOG_LV::LOG_ERROR, "t_goi_nav", "processBatch", int2str(_gpoleut1->getBegTime()));
				_glog->comment(t_glog::LOG_LV::LOG_ERROR, "t_goi_nav", "processBatch", int2str(_gpoleut1->getEndTime()));
				return false;
			}

			// get the number of equation and init the partial
			//_par_num = _sat_ics->get_par_num();
			if (_sat_ics->get_par_num() == 0)
			{
				writeLogInfo(_glog, 1, function_name, "ERROR : no data in the " + _sat + " ics file, check your ics file.");
				return false;
			}
			else
			{
				_neq_size = 6 * (_sat_ics->get_par_num() + 1);
			}

			// init sat info
			_init_sat_info(_sat);
			if (_gsat.prn().empty())
			{
				writeLogInfo(_glog, 1, function_name, "ERROR : no information for " + _sat + ", check your sat_paramters file.");
				return false;
			}

			// current time must between the start and the end
			// reference time from ics
			t_gtime integra_ref_time = _gics->beg_time;
			if ((_extend_beg >= integra_ref_time) || (_extend_end <= integra_ref_time))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : There are some conflicts between the xml and ics file in beg and end time");
				return false;
			}

			// get init state and partial form ics
			vector<double> state_init;
			if (!_init_state_matrix(_neq_size, state_init))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not init the state matrix.");
				return false;
			};

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
			//_glog->comment(0, _sat, "processBatch : numberical intergration finshed !");

			//write state info to orb file
			//init state 
			vector<double> init_funct(_neq_size / 2, 0.0);
			for (int i = 0; i < _neq_size / 6; i++)
			{
				init_funct[(int64_t)i * 3]     = state_init[(int64_t)i * 6];
				init_funct[(int64_t)i * 3 + 1] = state_init[(int64_t)i * 6 + 1];
				init_funct[(int64_t)i * 3 + 2] = state_init[(int64_t)i * 6 + 2];

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
			//set orb partial parlist
			auto ics_par = _sat_ics->get_par_list();
			vector<string > par_list;
			for (auto itpar : ics_par)
			{
				par_list.push_back(ics2str(itpar.first));
			}
			_gorb->add_parlist(_sat, par_list);
			//set orb head info for only one time
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
			delete _scf2crs_cord; _scf2crs_cord = NULL;
			delete _crs2acr_cord; _crs2acr_cord = NULL;

			return true;
		}
		catch (...)
		{
			if (_trs2crs_cord)
			{
				delete _trs2crs_cord; _trs2crs_cord = NULL;
			}
			if (_scf2crs_cord)
			{
				delete _scf2crs_cord; _scf2crs_cord = NULL;
			}
			if (_crs2acr_cord)
			{
				delete _crs2acr_cord; _crs2acr_cord = NULL;
			}
			writeLogInfo(_glog, 1, function_name, "ERROR : throw exception while doing procseeBatch of" + _sat);
			return false;
		}
	}

}