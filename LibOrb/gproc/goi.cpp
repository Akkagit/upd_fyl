/**
*
* @verbatim
    History
     -1.0 Yujie Qin  2019-05-23  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    goi.cpp
* @brief		orbit integration processing for satellite
*
* @author       Yujie Qin, Wuhan University
* @version		1.0.0
* @date		    2019-05-23
*
*/


#include <iostream>
#include <iomanip>
#include <algorithm>

#include "gproc/goi.h"
#include "gforce/gforce.h"
#include "gforce/gfcommon.h"
#include "gutils/ginfolog.h"
#include "newmat/newmatio.h"


namespace great
{
	/** @brief constructor. */
	t_goi::t_goi()
	{

	}

	/** 
	* @brief constructor. 
	* @param[in]	mark   leo name
	* @param[in]	set    setting information
	*/
	t_goi::t_goi(string mark, t_gsetbase* set)
	{
		gtrace("t_goi::t_goi");
		_sat = mark;
		_set = set;
	}

	/** @brief destructor. */
	t_goi::~t_goi()
	{

	}

	/**
	* @brief set LEO ICS data
	* @param[in] gics   LEO ICS data
	*/
	void t_goi::set_data_used(t_gallics_leo* gics)
	{
		
	}

	void t_goi::set_data_used(t_gallics_nav* gics)
	{

	}

	/**
	* @brief set poleut1 data
	* @param[in]  poleut1   poleut1 data
	*/
	void t_goi::set_data_used(t_gpoleut1* poleut1)
	{
		_gpoleut1 = poleut1;
	}

	/**
	* @brief set ocean load data
	* @param[in]  gotl   ocean load data
	*/
	void t_goi::set_data_used(t_gallotl * gotl)
	{
		_gotl = gotl;
	}


	/**
	* @brief  set EGM data
	* @param[in]  gegmcoef   EGM data
	*/
	void t_goi::set_data_used(t_gegmcoef* gegmcoef)
	{
		_gegmcoef = gegmcoef;
	}

	/**
	* @brief   set planets data
	* @param[in] gplanet  data for all planets
	*/
	void t_goi::set_data_used(t_gnavde* gplanet)
	{
		_gplanet = gplanet;
	}

	/**
	* @brief   set LEO pannel data
	* @param[in] gpannel  pannel data for one LEO satellite
	*/
	void t_goi::set_data_used(t_gpannel*  gpannel)
	{
		_gpannel = gpannel;
	}

	/**
	* @brief   set ocean_tide data
	* @param[in]  gotdata   ocean tide data
	*/
	void t_goi::set_data_used(t_gocean_tide* gotdata)
	{
		_gotdata = gotdata;
	}


	/**
	* @brief   set leap second data
	* @param[in]  leapsecond   leap second data
	*/
	void t_goi::set_data_used(t_gleapsecond  *leapsecond)
	{
		_gleapsecond = leapsecond;
	}
	

	/**
	* @brief   set sat_parameters_new data
	* @param[in] gallsat   data for all satellites
	*/
	void t_goi::set_data_used(t_gallsatinfo *gallsat)
	{
		_gallsat = gallsat;
	}
	
	void t_goi::set_data_used(t_gorb * gorb)
	{
		_gorb = gorb;
	}


	/**
	* @brief   set solar_flux and geomag_kp data
	* @param[in] solargeomag   solar_flux and geomag_kp data
	*/
	void t_goi::set_data_used(t_gsolargeomag* solargeomag)
	{
		_gsolargeomag = solargeomag;
	}

	/**
	* @brief   set albedo emit and refl data
	* @param[in] solargeomag   emit and refl data
	*/
	void t_goi::set_data_used(t_galbedo* albedodata)
	{
		_galbedo = albedodata;
	}



	/**

	* @brief   set ocean_pole_tide data
	* @param[in]  gotdata   ocean tide data
	*/
	void t_goi::set_data_used(t_gdeadata* gdedata)
	{
		_gdeadata = gdedata;
	}



	/**
	* @brief    set LEO attitude data
	* @param[in] attitude  attitude data for one LEO satellite
	*/
	void t_goi::set_data_used(t_gattitude*    attitude)
	{
		_gattitude = attitude;
	}

	///**
	//* @brief   get orb funct according to the satellite name 
	//* @param[in]  sat   satellite name
	//* @param[in]  gorb  orb data
	//*/
	//void t_goi::getOrbFunct(string sat, t_gorb* gorb)
	//{
	//	map<int, vector< t_gtriple> > funct_tmp;
	//	_gorb->get_funct(sat, funct_tmp);

	//	gorb->add_funct(sat, funct_tmp);
	//}
	//
	///**
	//* @brief  get orb header data according to the satellite name
	//* @param[in]  sys   system of the satelite
	//* @param[in]  sats  all satellite names
	//* @param[in]  gorb  orb data
	//*/
	//void t_goi::getOrbHead(set<string>& sys, set<string>& sats, t_gorb* gorb)
	//{
	//	t_orb_hd orb_hd_tmp;
	//	string function_name = "t_goi::get_orb_head";
	//	try
	//	{
	//		_gorb->get_head(orb_hd_tmp);
	//		//sat group
	//		if (sys.size() > 1)
	//		{
	//			orb_hd_tmp.satgroup = "LEO";
	//		}
	//		else
	//		{
	//			orb_hd_tmp.satgroup = *(sys.begin());
	//		}
	//		//sat type
	//		set<string>::iterator itsys = sys.begin();
	//		for (itsys; itsys != sys.end(); ++itsys)
	//		{
	//			//char sys_tmp = t_gsys::str2char(*itsys);
	//			orb_hd_tmp.sattype.insert('L');
	//		}
	//		//sat num
	//		orb_hd_tmp.nsat = sats.size();
	//		//sat list
	//		orb_hd_tmp.satprn = sats;
	//		gorb->add_head(orb_hd_tmp);
	//	}
	//	catch (...)
	//	{
	//		writeLogInfo(_glog, 2, function_name, "ERROR : unknown mistake");
	//	}
	//}


	/**
	* @brief numerical integration for OI
	* @param[in]  beg          extended begin time
	* @param[in]  end          extended end time
	* @param[in]  state_init   initial state data
	* @param[out] state_result the state data after the numerical integration
	* @return  numerical integrate correctly or not
	*/
	bool t_goi::_numerical_integration(const t_gtime& beg, const t_gtime& end,
		const vector<double>& state_init, t_map_vector& state_result)
	{
		string function_name = "t_goi::_numerical_integration";

		// init rkf settings
		const t_gtime& beg_time = beg;
		const t_gtime& end_time = end;
		t_adams* adams = nullptr;

		try
		{
			bool forward = (beg_time >= end_time ? true : false);
			int  sign = (forward ? 1 : -1);
			double rkf_h = _step_inp * sign / 30.0;

			if (rkf_h == 0.0)
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : the _step in was zero, check your process or xml file settings");
				return false;
			}

			// init adams 
			vector<list<double> > admas_state;
			vector<list<double> > admas_partial;

			// get the adams' state and partial by using rkf
			if (!_prepare_adams_with_rkf(_integra_degree, beg_time, rkf_h, state_init, admas_state, admas_partial))
			{
				writeLogInfo(_glog, 0, function_name, "ERROR : can not get the adams points with the rkf");
				return false;
			}
	
			int num_of_all = 0;

			// start adams
			int adams_sign = (forward ? -1 : 1);
			double adams_step = _step_inp * 1.0 * adams_sign;
			t_gtime adams_ref_time = beg_time;
			if(!adams) adams = new t_adams(_integra_degree, adams_step);

			// try at first to get the real epoch num;
			int epoch_num = 0;
			bool isValid = true;		
			t_gtime epoch_try = beg_time;
			while (epoch_try.dmjd()* adams_sign < end_time.dmjd() * adams_sign)
			{
				epoch_num++;
				epoch_try = beg_time + epoch_num * adams_step;

				if (!_processEpoch(epoch_try, adams, admas_state, admas_partial))
				{
					isValid = false;
					break;
				}

				int epo = epoch_num * adams_sign;

				for (int i = 0; i < _neq_size / 6; i++)
				{
					state_result[epo].push_back(admas_state[(int64_t)i * 6].back());
					state_result[epo].push_back(admas_state[(int64_t)i * 6 + 1].back());
					state_result[epo].push_back(admas_state[(int64_t)i * 6 + 2].back());
				}
			}

			if (adams) { delete adams;  adams = nullptr; }
			if (isValid)
			{
				return true;
			}
			else
			{
				writeLogInfo(_glog, 0, function_name,
					"ERROR : There are some problems in process Epoch while calculating with adams");
				return false;
			}

		}
		catch (...)
		{
			if (adams)
			{
				delete adams;
				adams = NULL;
			}
			writeLogInfo(_glog, 3, function_name, "ERROR : throw exception while numerical_integration");
			return false;
		}
	}

	/**
	* @brief OI processsing for one epoch.
	* @param[in]  epoch          current epoch
	* @param[in]  adams          adams integrator
	* @param[out] state_result   state result after calculation
	* @param[out] partial_result partial result after calculation
	* @return  calculate correctly or not
	*/
	bool t_goi::_processEpoch(t_gtime epoch, t_adams* adams,
		vector<list<double> >& state_result, vector<list<double> >& partial_result)
	{
		string function_name = "t_goi::_processEpoch";

		try
		{
			for (int iequ = 0; iequ < _neq_size; iequ++)
			{
				if (!adams->update("PRE", state_result[iequ], partial_result[iequ]))
				{
					writeLogInfo(_glog, 0, function_name, "ERROR : unknown mistake ");
				}
			}

			// state matrix and partial matrix
			vector<double> state(_neq_size, 0.0);
			vector<double> partial(_neq_size, 0.0);
			for (int iequ = 0; iequ < _neq_size; iequ++)
			{
				state[iequ] = (state_result[iequ].front());
			}

			_calculate(epoch, state, partial);

			for (int iequ = 0; iequ < _neq_size; iequ++)
			{
				partial_result[iequ].push_back(partial[iequ]);
				partial_result[iequ].pop_front();
				adams->update("COR", state_result[iequ], partial_result[iequ]);
			}
			for (int iequ = 0; iequ < _neq_size; iequ++)
			{
				state[iequ] = (state_result[iequ].back());
			}

			_calculate(epoch, state, partial);
			for (int iequ = 0; iequ < _neq_size; iequ++)
			{
				partial_result[iequ].pop_back();
				partial_result[iequ].push_back(partial[iequ]);
			}

			return true;
		}
		catch (...)
		{
			writeLogInfo(_glog, 0, function_name, "ERROR : unknown mistake");
			return false;
		}
	}

	/**
	* @brief   extend the time period
	* @param[in]  beg         begin time according to the XML file
	* @param[in]  end         end time according to the XML file
	* @param[in]  step_out    OI output step
	* @param[out] extend_beg  extended begin time
	* @param[out] extend_end  extended end time
	* @return  extend time correctly or not
	*/
	bool t_goi::_extend_time_period(const t_gtime& beg, const t_gtime& end, const int& step_out, t_gtime& extend_beg, t_gtime& extend_end)
	{
		string function_name = "t_goi::_extend_time_period";
		if (step_out <= 0)
		{
			writeLogInfo(_glog, 0, function_name, "ERROR : the step_out should not be " + int2str(step_out));
			return false;
		}
		try
		{
			double sec_beg = round((beg.dsec() + beg.sod()) * 1.0 / (step_out * 1.0)) *  step_out;
			double sec_end = round((end.dsec() + end.sod()) * 1.0 / (step_out * 1.0)) *  step_out;
			extend_beg.from_mjd(beg.mjd(), sec_beg, sec_beg - floor(sec_beg) );
			extend_end.from_mjd(end.mjd(), sec_end, sec_end - floor(sec_end));

			sec_beg = extend_beg.sod() + extend_beg.dsec() - (int64_t)NEXT_POINT_NUM * _step_out;
			sec_end = extend_end.sod() + extend_end.dsec() + (int64_t)NEXT_POINT_NUM * _step_out;
			extend_beg.from_mjd(extend_beg.mjd(), sec_beg, sec_beg - floor(sec_beg));
			extend_end.from_mjd(extend_end.mjd(), sec_end, sec_end - floor(sec_end));
		
			return true;
		}
		catch (...)
		{
			writeLogInfo(_glog, 2, function_name, "ERROR : throw exception while extending time period!");
			return false;
		}
	}

	/**
	* @brief   prepare adams with rkf
	* @param[in]  adams_integra_degree      order of adams integration
	* @param[in]  epoch     begin time of integration
	* @param[in]  step_h    step of rkf
	* @param[in]  rkf_init_state    initial state of rkf integration
	* @param[out] adams_state		initial state of adams integration
	* @param[out] adams_partial		initial partial of adams integration
	* @return  prepare adams with rkf correctly or not
	*/
	bool t_goi::_prepare_adams_with_rkf(const int& adams_integra_degree, const t_gtime& epoch, const double& step_h,
		const vector<double>& rkf_init_state, vector<list<double> >& adams_state, vector<list<double> >& adams_partial)
	{
		//double tn = epoch.dmjd(false); debug,qyj
		string function_name = "t_goi::_prepare_adams_with_rkf";

		try
		{
			// init rkf
			t_grkf_nav* rkf_nav = nullptr;
			if(!rkf_nav) rkf_nav = new t_grkf_nav(_neq_size, RKF_DEGREE, step_h, epoch, rkf_init_state);

			// loop rkf to get the adams point
			bool isfirst = true;
			t_gtime ti = epoch;
			vector< vector<double> > adams_state_tmp;
			vector< vector<double> > adams_partial_tmp;
			vector<double> fi(_neq_size, 0.0);
			vector<double> Xi(_neq_size, 0.0);

			// cycle to get all adams state and partial
			Xi = rkf_init_state;
			adams_state_tmp.push_back(rkf_init_state);
			for (int i = 1; i < adams_integra_degree; ++i)
			{
				for (int j = 0; j < RKF_POINT_NUM; ++j)
				{
					for (int k = 0; k <= RKF_COEF_I; ++k)
					{
						_calculate(ti, Xi, fi);
						if (isfirst)
						{
							adams_partial_tmp.push_back(fi);
							isfirst = false;
						}
						rkf_nav->updatRkfPoint(k, fi, ti, Xi);
					}
					Xi = rkf_nav->getRkfResultVector();
					rkf_nav->resetRkf(Xi, ti);
				}
				// get the intergra result
				adams_state_tmp.push_back(Xi);
				_calculate(ti, Xi, fi);
				adams_partial_tmp.push_back(fi);
			}

#ifdef _DEBUG_LOG
			for (int i = 0; i < 11; ++i)
			{
			for (int j = 0; j < _neq_size/3; j++)
			{
			cout <<setw(25)<<setprecision(14)<<scientific<<right <<adams_state_tmp[i][j*3] << adams_state_tmp[i][j*3 + 1] << adams_state_tmp[i][j*3 + 2] << endl;
			}
			cout << endl;
			}
#endif
			for (int i = 0; i < _neq_size; ++i)
			{
				list<double> one_list_state;
				list<double> one_list_partial;
				for (int j = (adams_integra_degree - 1); j >= 0; --j)
				{
					one_list_state.push_back(adams_state_tmp[j][i]);
					one_list_partial.push_back(adams_partial_tmp[j][i]);
				}
				adams_state.push_back(one_list_state);
				adams_partial.push_back(one_list_partial);
			}
			if (rkf_nav) { delete rkf_nav; rkf_nav = nullptr; }
			return true;
		}
		catch (...)
		{
			writeLogInfo(_glog, 2, function_name, "ERROR: unknown mistake");
			return false;
		}

	}

	
	
	/**
	* @brief  batch processing for  OI
	* @param[in]  beg         begin time from XML
	* @param[in]  end         end time from XML
	* @return  calculate correctly or not
	*/
	bool  t_goi::processBatch(const t_gtime& beg, const t_gtime& end)
	{
		return true;
	}

	/**
    * @brief  calculate the acc and partial in epoch t
    * @param[in]  t        current epoch for calculation
    * @param[out] state    state data after calculation
    * @param[out] partial  partial data after calculation
    * @return bool
    	@retval true     calculate correctly
    	@retval false    calculate incorrectly
    */
	//bool  t_goi::_calculate(const t_gtime& t, const vector<double>& state, vector<double>& partial)
	//{
	//	return true;
	//}
}