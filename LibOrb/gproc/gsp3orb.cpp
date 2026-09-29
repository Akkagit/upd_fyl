/**
*
* @verbatim
	History
	 -1.0 Yujie Qin		2019-10-10 creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gsp3orb.cpp
* @brief		convert sp3 file or kin file to orb file
*
* @author       Yujie Qin	, Wuhan University
* @version		1.0.0
* @date		    2019-10-10
*
*/

#include <iostream>
#include <iomanip>
#include <algorithm>
#include "newmat/newmatio.h"
#include "gproc/gsp3orb.h"
using namespace gnut;
namespace great
{
	t_gsp3orb::t_gsp3orb()
	{

	}
	t_gsp3orb::t_gsp3orb(set<string> mark, set<string> sys, t_gsetbase* set, t_glog * log)
	{
		_flag = false;
		_sats = mark;
		_sys = sys;
		_set = set;
		_glog = log;
		_gorb = new t_gorb;
		_epoch.tsys(t_gtime::GPS);
		
	}

	t_gsp3orb::~t_gsp3orb()
	{
		delete _gorb; _gorb = NULL;
		
	}

	void t_gsp3orb::setDataUsed(t_gallprec* gsp3)
	{
		_sp3 = gsp3;
	}
	void t_gsp3orb::setDataUsed(t_gpoleut1* poleut1)
	{
		_poleut1 = poleut1;
	}
	bool t_gsp3orb::processBatch(const t_gtime& beg, const t_gtime& end)
	{
		_beg = beg;
		_end = end;
		if (!_initial_data())
		{
			_glog->comment(0, "sp3orb", "ERROR : sp3orb have no comeplete file !");
			return false;
		}
		if (!_get_settings())
		{
			_glog->comment(0, "sp3orb", "ERROR : sp3orb have something wrong in xml-file !");
			return false;
		}
		set<string> sat_sp3 = _sp3->satellites();
		set<string>::iterator itsat = _sats.begin();
		//_interval = int(dynamic_cast<t_gsetgen*>(_set)->sampling());
		/*if (!_interval)
		{*/
		_interval = _sp3->intv(*sat_sp3.begin());
			//if (_flag == false)
			//{
			//	_interval = 300;
			//	//to do
			//}
		//}
			
		_trs2crs = new t_gtrs2crs(false, _poleut1);
		if (!_extend_time_period(beg, end, _interval, _beg_extend, _end_extend))
		{
			return false;
		}
		_ref = _beg_extend + NEXT_POINT_NUM * _interval;
		_epoch = _beg_extend;
		int epo = -1*NEXT_POINT_NUM;
		Matrix rot_trs2crs(0.0, 3, 3);	// rotation matrix of TRS to CRS

		for (itsat; itsat != _sats.end(); itsat++)
		{
			if (sat_sp3.find(*itsat) == sat_sp3.end())
			{
				_glog->comment(1, "sp3orb", "WARNING : " + *itsat + "can not be found in sp3 file!");
				_sats.erase(itsat++);
				continue;
			}
			if (_beg_extend < _sp3->beg_data(*itsat) || _end_extend> _sp3->end_data(*itsat))
			{
				if (_flag&&_sats.size() == 1)
				{
					_beg_extend = _beg_extend < _sp3->beg_data(*itsat) ? _sp3->beg_data(*itsat) : _beg_extend;
					_end_extend = _end_extend > _sp3->end_data(*itsat) ? _sp3->end_data(*itsat) : _end_extend;
				}
				else
				{
					_glog->comment(0, "sp3orb", "WARNING : " + *itsat + "time of kin or sp3 file is wrong!");
					_sats.erase(itsat++);
					continue;
				}

			}
		}
		if (_sats.size() == 0)
		{
			return false;
		}
		while (_epoch <= _end_extend)
		{
			t_gtime tdt_gtime = _epoch;
			tdt_gtime.tsys(t_gtime::TT);
			//ChangeTimeSys(GPS2TDT, t, tdt_gtime);
			double tdt = tdt_gtime.dmjd();
			if (_frame == "crs")
			{
				_trs2crs->calcRotMat(tdt_gtime, true, true, true);
				// rotation matrix of TRS to CRS
				rot_trs2crs = _trs2crs->getRotMat();
			}
			
		
			ColumnVector trs_pos, trs_vel, crs_pos, crs_vel;
			ColumnVector target_pos;
			set<string>::iterator itsat = _sats.begin();
			for (itsat; itsat != _sats.end(); )
			{
				vector<double> one_epo_funct;
				one_epo_funct.resize(_neq_size[*itsat]);
				t_gtriple pos, vel;
				_sp3->lagrange_pos(*itsat, _epoch, pos);
				//_sp3->lagrange_vel(*itsat, _epoch, vel);
				trs_pos = pos.crd_cvect() / 1000;
				if (double_eq(trs_pos[0] * trs_pos[1] * trs_pos[2], 0.0) || trs_pos[0] >= 1E18) {
					// _sat.erase(iter++); // is it necessary?
					itsat++;
					continue;
				}
				//trs_vel = vel.crd_cvect()/1000 / _step_out;

				if (_frame == "crs")
				{
					crs_pos = rot_trs2crs * trs_pos;
					//crs_vel = rot_trs2crs * trs_vel;
					target_pos = crs_pos;
				}
				else
				{
					target_pos = trs_pos;
				}
				one_epo_funct[0] = target_pos(1);
				one_epo_funct[1] = target_pos(2);
				one_epo_funct[2] = target_pos(3);
				vector <string> list;
				for (auto i = _param[*itsat].begin(); i != _param[*itsat].end(); i++)
				{
					list.push_back(*i);
				}
				_gorb->add_parlist(*itsat, list);
				_gorb->add_funct(epo, *itsat, one_epo_funct);
				++itsat;
			}
			
			
			epo++;
			_epoch = _epoch + _interval;
		}
		if (_trs2crs)
		{
			delete _trs2crs;
			_trs2crs = nullptr;
		}
		return true;
	}
	bool t_gsp3orb::_initial_data()
	{
		if (!_poleut1 || !_sp3)
			return false;
		if (_sp3->get_sp3_size() == 0)
			return false;
		return true;
	}

	bool t_gsp3orb::_get_settings()
	{
		_frame = dynamic_cast<t_gsetsp3orb*>(_set)->frame();
		if (_frame == "")
		{
			_frame = "crs";
		}
		for (auto itsat = _sats.begin(); itsat != _sats.end(); itsat++)
		{
			_flag = IS_SUPPORTED_LEO(*itsat);//need stability	

			_param_size[*itsat] = 6;
			_param[*itsat].push_back("PXSAT");
			_param[*itsat].push_back("PYSAT");
			_param[*itsat].push_back("PZSAT");
			_param[*itsat].push_back("VXSAT");
			_param[*itsat].push_back("VYSAT");
			_param[*itsat].push_back("VZSAT");
			if (!_flag)
			{
				_step_inp[*itsat] = dynamic_cast<t_gsetforce*>(_set)->getInpStep(*itsat);
				_step_out[*itsat] = dynamic_cast<t_gsetforce*>(_set)->getOutStep(*itsat);
				_integra_degree[*itsat] = dynamic_cast<t_gsetforce*>(_set)->getIntegraDegree(*itsat);

				map<FORCE, string> force_models = dynamic_cast<t_gsetforce*>(_set)->getForceModelMap(*itsat);
				map<FORCE, bool> force_valid = dynamic_cast<t_gsetforce*>(_set)->getForcesValidMap(*itsat);

				if (force_valid[great::SOLAR_RADI] == true && force_models[great::SOLAR_RADI].find("ECOM_5") != string::npos)
				{
					_param_size[*itsat] = _param_size[*itsat] + 5;
					_param[*itsat].push_back("Kd_BERN");
					_param[*itsat].push_back("Ky_BERN");
					_param[*itsat].push_back("Kb_BERN");
					_param[*itsat].push_back("Kbc1_BERN");
					_param[*itsat].push_back("Kbs1_BERN");
				}
				else if (force_models[great::SOLAR_RADI].find("ECOM_9") != string::npos)
				{
					_param_size[*itsat] = _param_size[*itsat] + 9;
					//todo
				}
				else
				{
					_param_size[*itsat] = _param_size[*itsat] + 5;
					_param[*itsat].push_back("Kd_BERN");
					_param[*itsat].push_back("Ky_BERN");
					_param[*itsat].push_back("Kb_BERN");
					_param[*itsat].push_back("Kbc1_BERN");
					_param[*itsat].push_back("Kbs1_BERN");
				}
				_neq_size[*itsat] = (_param_size[*itsat] + 1) * 3;
			}
			else
			{
				_step_inp[*itsat] = dynamic_cast<t_gsetforce*>(_set)->getInpStep(*itsat);
				_step_out[*itsat] = dynamic_cast<t_gsetforce*>(_set)->getOutStep(*itsat);
				_integra_degree[*itsat] = dynamic_cast<t_gsetforce*>(_set)->getIntegraDegree(*itsat);

				_drag_para[*itsat] = dynamic_cast<t_gsetforce*>(_set)->getAtmospherePars(*itsat);
				_srp_para[*itsat]  = dynamic_cast<t_gsetforce*>(_set)->getSolarRadiationPars(*itsat);
				_cust_para[*itsat] = dynamic_cast<t_gsetforce*>(_set)->getCustomerPars(*itsat);

				for (auto ipa = _srp_para[*itsat].begin(); ipa != _srp_para[*itsat].end(); ipa++)
					_param[*itsat].push_back((*ipa).first);
				_param_size[*itsat] = _param_size[*itsat] + _srp_para[*itsat].size();
				
				for (auto ipa = _drag_para[*itsat].begin(); ipa != _drag_para[*itsat].end(); ipa++)
					_param[*itsat].push_back((*ipa).first);
				_param_size[*itsat] = _param_size[*itsat] + _drag_para[*itsat].size();
				
				for (auto ipa = _cust_para[*itsat].begin(); ipa != _cust_para[*itsat].end(); ipa++)
					_param[*itsat].push_back((*ipa).first);
				_param_size[*itsat] = _param_size[*itsat] + _cust_para[*itsat].size();
				//_param.insert(_drag_para.begin(), _drag_para.end());
				//_param.insert(_cust_drag.begin(), _cust_drag.end());

				_neq_size[*itsat] = (_param_size[*itsat] + 1) * 3;
			}
		}
		return true;
	}

	bool t_gsp3orb::_extend_time_period(t_gtime beg, t_gtime end, int step_out, t_gtime& extend_beg, t_gtime& extend_end)
	{
		try
		{
			t_gtime sp3_beg = _sp3->beg_data();
			t_gtime sp3_end = _sp3->end_data();
			double sec_beg = round((beg.dsec() + beg.sod()) * 1.0 / (step_out * 1.0)) *  step_out;
			double sec_end = round((end.dsec() + end.sod()) * 1.0 / (step_out * 1.0)) *  step_out;
			extend_beg.from_mjd(beg.mjd(), sec_beg, sec_beg - floor(sec_beg));
			extend_end.from_mjd(end.mjd(), sec_end, sec_end - floor(sec_end));

			sec_beg = extend_beg.sod() + extend_beg.dsec() - NEXT_POINT_NUM * step_out;
			sec_end = extend_end.sod() + extend_end.dsec() + NEXT_POINT_NUM * step_out;
			extend_beg.from_mjd(extend_beg.mjd(), sec_beg, sec_beg - floor(sec_beg));
			extend_end.from_mjd(extend_end.mjd(), sec_end, sec_end - floor(sec_end));
			return true;

			if (extend_beg < sp3_beg || extend_end > sp3_end)
			{
				_glog->comment(1, "sp3orb", "WARNING : time is too long!");
				return false;

			}
		}
		catch (...)
		{
			return false;
		}

	}


	void t_gsp3orb::_generate_icsfile()
	{

		set<string>::iterator isat = _sats.begin();
		string satname;
		bool isfirst = true;
		for (isat; isat != _sats.end(); isat++)
		{
			satname = *isat;

			vector<string>::iterator itparam = _param[*isat].begin();
			_ref_state[satname][*itparam] = 0.0;
			double x, y, z, vx, vy, vz;
			if (!_gorb->get_xyz(_ref, satname, x, y, z))
			{
				continue;
			}
			_gorb->get_vel(_ref, satname, vx, vy, vz);
			_ref_state[satname]["PXSAT"] = x / 1000;
			_ref_state[satname]["PYSAT"] = y / 1000;
			_ref_state[satname]["PZSAT"] = z / 1000;
			_ref_state[satname]["VXSAT"] = vx / 1000;
			_ref_state[satname]["VYSAT"] = vy / 1000;
			_ref_state[satname]["VZSAT"] = vz / 1000;
			if (abs(x + y + z) == 0)
			{
				_sats.erase(satname);
			}
			if (_flag == false)
			{
				if (isfirst == true)
				{
					_ics_nav = new t_gallics_nav();
					_ics_nav->beg_time = _beg;
					_ics_nav->end_time = _end;
					_ics_nav->sys_type = "GNS";
					_ics_nav->sys_list = _sys;
					_ics_nav->sat_list = _sats;
					_ics_nav->sat_num = _sats.size();
					isfirst = false;
				}
				
				t_gics_nav _gics;
				_gics.sat_name = satname;
				t_ics_record    _record;
				_record.beg_time = _beg;
				_record.end_time = _end;
				vector<string>::iterator itparam = _param[*isat].begin();
				for (itparam; itparam != _param[*isat].end(); ++itparam)
				{
					_record.str_name = *itparam;
					_record.enum_name = str2ics(*itparam);
				    _record.value = _ref_state[satname][*itparam];
					_record.sigma = 0.0;
					_gics.set_par_record(_record.str_name, _record);
				}
				_ics_nav->add_nav_ics(_gics.sat_name, _gics);

				
			}
			else
			{
				if (isfirst == true)
				{
					_ics_leo = new t_gallics_leo();
					_ics_leo->set_time(_beg, _end);
					isfirst = false;
				}

				t_gics_leo _gics;
				vector<t_icsleo_record> par_all;
				_gics.sat_name = satname;
				_gics.sat_type = satname;
				t_icsleo_record    record;

				vector<string>::iterator itparam = _param[*isat].begin();
				for (itparam; itparam != _param[*isat].end(); ++itparam)
				{
					string para_name = *itparam;
					if (para_name == "PXSAT" || para_name == "PYSAT" || para_name == "PZSAT"
						|| para_name == "VXSAT" || para_name == "VYSAT" || para_name == "VZSAT")
					{
						record.beg = _beg;
						record.end = _end;
						record.str_name = *itparam;
						record.enum_name = str2ics(para_name);
						record.value = _ref_state[satname][*itparam];
						record.sig = 0.0;
						par_all.push_back(record);
						_gics.set_par_record(par_all[0].str_name, par_all);
						par_all.clear();
					}
					else if (para_name == "DRAG_c" || para_name == "SR_scale")
					{
						int min, num;
						if (para_name == "DRAG_c")
							min = _drag_para[*isat][para_name];
						else
							min = _srp_para[*isat][para_name];
						if (min == 0) 
							num = 1;
						else 
							num = round((_end - _beg) / min / 60);
						record.str_name = para_name;
						record.enum_name = str2icsleo(para_name);
						for (int i = 0; i < num-1; i++)
						{
							record.beg = _beg + i * min * 60;
							record.end = record.beg + min * 60;
							record.value = 0;
							record.sig = 0.0;
							par_all.push_back(record);
						}
						record.beg = _beg + (num - 1) * min * 60;
						record.end = _end;
						record.value = 0;
						record.sig = 0.0;
						par_all.push_back(record);

						par_all[0].value = _ref_state[satname][para_name];
						_gics.set_par_record(par_all[0].str_name, par_all);
						par_all.clear();
					}
					else
					{

						int min = _cust_para[*isat][para_name];
						int num = round((_end - _beg) / min / 60);
						record.str_name = para_name;
						record.enum_name = str2icsleo(para_name);
						for (int i = 0; i < num-1; i++)
						{
							record.beg = _beg + i * min * 60;;
							record.end = record.beg + min * 60;
							record.value = 0;
							record.sig = 0.0;
							par_all.push_back(record);
						}
						record.beg = _beg + (num - 1) * min * 60;;
						record.end = _end;
						record.value = 0;
						record.sig = 0.0;

						par_all.push_back(record);
						par_all[0].value = _ref_state[satname][para_name];
						_gics.set_par_record(par_all[0].str_name, par_all);
						par_all.clear();
					}
				}
				_ics_leo->add_leo_ics(satname, _gics);
				
			}
		}
		
		


	}

	void t_gsp3orb::_generate_orbfile()
	{
		t_orb_hd orb_hd;
		orb_hd.adminp = _interval;//yjqin
		orb_hd.admout = _interval;
		orb_hd.admord = _integra_degree.begin()->second;
		orb_hd.nequ = _neq_size.begin()->second;
		orb_hd.nepo = (_end_extend - _beg_extend) / _interval + 1;
		orb_hd.beg = _beg_extend;
		orb_hd.end = _end_extend;
		orb_hd.t_ref = _ref;
		orb_hd.partial = true;

		if (_flag == false)
		{
			if (_sys.size() > 1)
			{
				orb_hd.sat_group = "GNS";
			}
			else
			{
				orb_hd.sat_group = *(_sys.begin());
			}
			set<string>::iterator itsys = _sys.begin();
			for (itsys; itsys != _sys.end(); ++itsys)
			{
				char sys_tmp = t_gsys::str2char(*itsys);
				orb_hd.sat_type.insert(sys_tmp);
			}
		}
		else
		{
			orb_hd.sat_group = "LEO";

			set<string>::iterator itsys = _sats.begin();
			for (itsys; itsys != _sats.end(); ++itsys)
			{
				char sys_tmp = 'L';
				orb_hd.sat_type.insert(sys_tmp);
			}
		}
		
		//sat num
		orb_hd.sat_num = _sats.size();
		//sat list
		orb_hd.sat_prn = _sats;
		_gorb->add_head(orb_hd);
	}

	void t_gsp3orb::GenerateProduct()
	{
		_generate_orbfile();
		_generate_icsfile();
		

		// Update ICS file
		if (_flag == false)
		{
			string ics_path = dynamic_cast<t_gsetout*>(_set)->outputs("ics");
			if (ics_path.empty())
			{
				ics_path = "ics_gnss";
			}
			t_gio* gout = new t_gfile; t_gcoder* gics_coder = new t_ics_nav(_set, "", 4096);
			gout->path(ics_path); gout->glog(_glog);
			gics_coder->clear(); gics_coder->glog(_glog); gics_coder->path(ics_path);
			gics_coder->add_data("ID1", _ics_nav);
			gout->coder(gics_coder);
			gout->run_write();
			
			string orb_path = dynamic_cast<t_gsetout*>(_set)->outputs("orb");
			if (orb_path.empty())
			{
				orb_path = "orb_gnss";
			}
			t_gio* gout_orb = new t_gfile;  t_gcoder* gorb_coder = new t_orb(_set, "", 4096);
			//string ics_path = (dynamic_cast<t_gsetinp*>(_set)->inputs("ics"))[0];
			gout_orb->path(orb_path); gout_orb->glog(_glog);
			gorb_coder->clear(); gorb_coder->glog(_glog); gorb_coder->path(orb_path);
			gorb_coder->add_data("ID1", _gorb);
			gout_orb->coder(gorb_coder);
			gout_orb->run_write();
			delete gout; gout = nullptr;
			delete gout_orb; gout_orb = nullptr;
			delete gics_coder; gics_coder = nullptr;
			delete gorb_coder; gorb_coder = nullptr;
			delete _ics_nav;_ics_nav = nullptr;
			
		}
		else
		{

			string ics_path = dynamic_cast<t_gsetout*>(_set)->outputs("icsleo");
			if (ics_path.empty())
			{
				ics_path = "ics_leo";
			}
			t_gio* gout = new t_gfile; 
			t_gcoder* gics_coder = new t_ics_leo(_set, "", 4096);
			gout->path(ics_path); gout->glog(_glog);
			gics_coder->clear(); gics_coder->glog(_glog); gics_coder->path(ics_path);
			gics_coder->add_data("ID1", _ics_leo);
			gout->coder(gics_coder);
			gout->run_write();

			string orb_path = dynamic_cast<t_gsetout*>(_set)->outputs("orb");
			if (orb_path.empty())
			{
				orb_path = "orb_leo";
			}
			t_gio* gout_orb = new t_gfile;  t_gcoder*gorb_coder = new t_orb(_set, "", 4096);
			//string ics_path = (dynamic_cast<t_gsetinp*>(_set)->inputs("ics"))[0];
			gout_orb->path(orb_path); gout_orb->glog(_glog);
			gorb_coder->clear(); gorb_coder->glog(_glog); gorb_coder->path(orb_path);
			gorb_coder->add_data("ID1", _gorb);
			gout_orb->coder(gorb_coder);
			gout_orb->run_write();
			
			delete gout; gout = nullptr;
			delete gout_orb; gout_orb = nullptr;
			delete gics_coder; gics_coder = nullptr;
			delete gorb_coder; gorb_coder = nullptr;
			delete _ics_leo;_ics_leo = nullptr;
			
		}
		


	}
}
