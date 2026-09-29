/**
*
* @verbatim
	History
	 -1.0 Jiande Huang  2019-09-20  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gpreedit.cpp
* @brief		main control class for estimate clk
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-09-20
*

* 功能设计：为了读取配置文件，准备lsq步骤前的所有文件
* 功能需求：
*	1. 生成ics文件
*	2. 生成tb 文件
*   3. 检查必要的数据文件是否存在
*
*/


#include "gproc/gpreedit.h"
#include "gset/gsetgen.h"
#include "gset/gsetgnss.h"
#include "gset/gsetxml.h"
#include "gset/gsetout.h"
#include "gio/gfile.h"
#include "gset/gsetforce.h"

#include "gutils/gtriple.h"
#include "gutils/gstring.h"


#include "gutils/ginfolog.h"
#include "gall/gallnav.h"

#include "gcoders/rinexc.h"

#ifdef WIN32
#include <io.h>
#include <direct.h> 
#else
#include <unistd.h>
#include <sys/stat.h>
#endif
#include <stdint.h>
#include <string>

#define MAX_PATH_LEN 256

#ifdef WIN32
#define ACCESS(fileName,accessMode) _access(fileName,accessMode)
#define MKDIR(path) _mkdir(path)
#else
#define ACCESS(fileName,accessMode) access(fileName,accessMode)
#define MKDIR(path) mkdir(path,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif

#include "gcoders/icsnav.h"
namespace great
{
	t_gpreedit::t_gpreedit()
	{
	}

	t_gpreedit::t_gpreedit(t_gsetbase * set, t_glog * log) :
		_gset(set),
		_log(log)
	{
		_beg_time = dynamic_cast<t_gsetgen*>(set)->beg();
		_end_time = dynamic_cast<t_gsetgen*>(set)->end();
		_rec_list = dynamic_cast<t_gsetgen*>(set)->recs();
		_sat_list = dynamic_cast<t_gsetgnss*>(set)->sat();
		_sys_list = dynamic_cast<t_gsetgen*>(set)->sys();
		_intv = dynamic_cast<t_gsetgen*>(set)->sampling();

		_def_ECOM_3.insert(ICS_PAR::PXSAT);
		_def_ECOM_3.insert(ICS_PAR::PYSAT);
		_def_ECOM_3.insert(ICS_PAR::PZSAT);
		_def_ECOM_3.insert(ICS_PAR::VXSAT);
		_def_ECOM_3.insert(ICS_PAR::VYSAT);
		_def_ECOM_3.insert(ICS_PAR::VZSAT);
		_def_ECOM_3.insert(ICS_PAR::ECOM_D0);
		_def_ECOM_3.insert(ICS_PAR::ECOM_Y0);
		_def_ECOM_3.insert(ICS_PAR::ECOM_X0);

		_def_ECOM_5 = _def_ECOM_3;
		/*_def_ECOM_5.insert(ICS_PAR::PXSAT);
		_def_ECOM_5.insert(ICS_PAR::PYSAT);
		_def_ECOM_5.insert(ICS_PAR::PZSAT);
		_def_ECOM_5.insert(ICS_PAR::VXSAT);
		_def_ECOM_5.insert(ICS_PAR::VYSAT);
		_def_ECOM_5.insert(ICS_PAR::VZSAT);
		_def_ECOM_5.insert(ICS_PAR::ECOM_D0);
		_def_ECOM_5.insert(ICS_PAR::ECOM_Y0);
		_def_ECOM_5.insert(ICS_PAR::ECOM_X0);*/
		_def_ECOM_5.insert(ICS_PAR::ECOM_Xc);
		_def_ECOM_5.insert(ICS_PAR::ECOM_Xs);

		_def_ECOM_9 = _def_ECOM_5;
		_def_ECOM_9.insert(ICS_PAR::ECOM_Dc);
		_def_ECOM_9.insert(ICS_PAR::ECOM_Ds);
		_def_ECOM_9.insert(ICS_PAR::ECOM_Yc);
		_def_ECOM_9.insert(ICS_PAR::ECOM_Ys);

		_def_ECOM2.insert(ICS_PAR::PXSAT);
		_def_ECOM2.insert(ICS_PAR::PYSAT);
		_def_ECOM2.insert(ICS_PAR::PZSAT);
		_def_ECOM2.insert(ICS_PAR::VXSAT);
		_def_ECOM2.insert(ICS_PAR::VYSAT);
		_def_ECOM2.insert(ICS_PAR::VZSAT);
		_def_ECOM2.insert(ICS_PAR::ECOM2_D0);
		_def_ECOM2.insert(ICS_PAR::ECOM2_Y0);
		_def_ECOM2.insert(ICS_PAR::ECOM2_X0);
		_def_ECOM2.insert(ICS_PAR::ECOM2_D2c);
		_def_ECOM2.insert(ICS_PAR::ECOM2_D2s);
		_def_ECOM2.insert(ICS_PAR::ECOM2_D4c);
		_def_ECOM2.insert(ICS_PAR::ECOM2_D4s);
		_def_ECOM2.insert(ICS_PAR::ECOM2_X1c);
		_def_ECOM2.insert(ICS_PAR::ECOM2_X1s);

		_def_solar_model = F_SOLAR_ECOM_5;
	}

	t_gpreedit::~t_gpreedit()
	{
		if (_trs2crs)
		{
			delete _trs2crs;
			_trs2crs = nullptr;
		}
	}

	bool t_gpreedit::ProcessBatch(const t_gtime& beg, const t_gtime& end)
	{
		if (beg < _beg_time) _beg_time = beg;
		if (end > _end_time) _end_time = end;

		_gnerate_ics_file(_beg_time, _end_time, _intv);
		//_gnerate_clk_file(_beg_time, _end_time, _intv);
		_gnerate_rec_crd_xml(_beg_time, _end_time);
		return true;
	}

	void t_gpreedit::set_nav(t_gallnav* nav)
	{
		_gnav = nav;
	}

	void t_gpreedit::set_erp(t_gpoleut1* erp)
	{
		_trs2crs = new t_gtrs2crs(false, erp);
	}

	void t_gpreedit::set_obj(t_gallobj* obj)
	{
		_gobj = obj;
	}

	void t_gpreedit::set_obs(t_gallobs* obs)
	{
		_gobs = obs;
	}

	bool t_gpreedit::_gnerate_ics_file(const t_gtime& beg, const t_gtime& end, double intv)
	{
		if (!_gnav) return false;

		_gics.beg_time = beg;
		_gics.end_time = end;
		_gics.sat_num = _sat_list.size();
		_gics.sat_list = _sat_list;
		_gics.sys_list = _sys_list;
		_gics.sys_type = "GNS";

		// update ics data
		for (auto sat_iter = _sat_list.begin(); sat_iter != _sat_list.end();)
		{
			string sat = *sat_iter;
			if (_sys_list.count(t_gsys::gsys2str(t_gsys::sat2gsys(sat))) == 0)
			{
				_log->comment(t_glog::LOG_LV::LOG_WARN, "preedit", "get_sat_xyz", "skip for no setting in xml gen sys : " + sat);
				_sat_list.erase(sat_iter++);
				_gics.sat_list = _sat_list;
				_gics.sat_num = _sat_list.size();
				continue;
			}

			t_gics_nav ics(sat);
			set<ICS_PAR> par_list;
			auto force = dynamic_cast<t_gsetforce*>(_gset)->getForceModelMap(sat);
			string model = force[SOLAR_RADI];
			if (model.find("5") != string::npos || model.find("Default") != string::npos) par_list = _def_ECOM_5;
			else if (model.find("ECOM2") != string::npos) par_list = _def_ECOM2;
			else if (model.find("ECOM_3") != string::npos) par_list = _def_ECOM_3;
			else par_list = _def_ECOM_9;

			// get the vel and pos of sat
			double xyz[3] = { 0.0 };
			double var[3] = { 0.0 };
			double vel[3] = { 0.0 };

			
			_gnav->pos(sat, beg, xyz, var, vel, true);
			if (double_eq(xyz[0], 0.0) && double_eq(xyz[1], 0.0) && double_eq(xyz[2], 0.0))
			{
				_gnav->t_gallnav::pos(sat, beg, xyz, var, vel, true);
			}

			if (double_eq(xyz[0], 0.0) && double_eq(xyz[1], 0.0) && double_eq(xyz[2], 0.0))
			{
				_log->comment(t_glog::LOG_LV::LOG_WARN, "preedit", "get_sat_xyz", "can not get the xyz from sp3 or rinexn file " + sat);
				_sat_list.erase(sat_iter++);
				_gics.sat_list = _sat_list;
				_gics.sat_num = _sat_list.size();
				continue;
			}
			t_gtime tdt_beg = beg;
			t_gtime tdt_end = end;
			tdt_beg.tsys(t_gtime::TT);
			tdt_end.tsys(t_gtime::TT);
			_trs2crs->calcRotMat(tdt_beg, true, true, true);

			ColumnVector trs_xyz(xyz, 3);
			ColumnVector trs_vel(vel, 3);
			Matrix trs2crs = _trs2crs->getRotMat();
			ColumnVector crs_xyz = trs2crs * trs_xyz;

			ColumnVector x_earth = (crs_xyz.t() * _trs2crs->getMatDu() / RAD2TSEC).t();
			ColumnVector crs_vel = trs2crs * (trs_vel + x_earth);
			
			//unit in ics is km
			crs_xyz = crs_xyz / 1E3;
			crs_vel = crs_vel / 1E3;
			for(auto par_iter = par_list.begin(); par_iter != par_list.end(); par_iter++)
			{		
				t_ics_record record;
				ICS_PAR par = *par_iter;
				record.beg_time = beg;
				record.end_time = end;
				record.enum_name = par;
				record.str_name = ics2str(par);

				switch (par) 
				{
					case ICS_PAR::PXSAT:
						record.value = crs_xyz(1);
						record.sigma = var[0];
						break;
					case ICS_PAR::PYSAT:
						record.value = crs_xyz(2);
						record.sigma = var[1];
						break;
					case ICS_PAR::PZSAT:
						record.value = crs_xyz(3);
						record.sigma = var[2];
						break;
					case ICS_PAR::VXSAT:
						record.value = crs_vel(1);
						record.sigma = 0.0;
						break;
					case ICS_PAR::VYSAT:
						record.value = crs_vel(2);
						record.sigma = 0.0;
						break;
					case ICS_PAR::VZSAT:
						record.value = crs_vel(3);
						record.sigma = 0.0;
						break;
					case ICS_PAR::ECOM_D0:
					case ICS_PAR::ECOM_Dc:
					case ICS_PAR::ECOM_Ds:
					case ICS_PAR::ECOM_Y0:
					case ICS_PAR::ECOM_Yc:
					case ICS_PAR::ECOM_Ys:
					case ICS_PAR::ECOM_X0:
					case ICS_PAR::ECOM_Xc:
					case ICS_PAR::ECOM_Xs:
					case ICS_PAR::ECOM2_D0:
					case ICS_PAR::ECOM2_Y0:
					case ICS_PAR::ECOM2_X0:
					case ICS_PAR::ECOM2_D2c:
					case ICS_PAR::ECOM2_D2s:
					case ICS_PAR::ECOM2_D4c:
					case ICS_PAR::ECOM2_D4s:
					case ICS_PAR::ECOM2_X1c:
					case ICS_PAR::ECOM2_X1s:
						record.value = 0.0;
						record.sigma = 0.0;
						break;
					default:
						break;
				}

				ics.set_par_record(record.str_name, record);
			}

			_gics.add_nav_ics(sat, ics);
			sat_iter++;
		}

		//output new icsfile
		//string ics_path = "./grt" + beg.str_gwkd(false) + ".ics";
		string ics_path = dynamic_cast<t_gsetout*>(_gset)->outputs("ics");
		if (ics_path.empty()) {
			ics_path = "grt$(date).ics";
		}

		if (ics_path.find("$(gwkd)") != string::npos)
		{
			substitute(ics_path, "$(gwkd)", beg.str_gwkd(), false);
		}
		else if (ics_path.find("$(date)") != string::npos)
		{
			substitute(ics_path, "$(date)", beg.str_yyyydoy(), false);
		}

		t_gfile gout;
		t_ics_nav gics_coder(_gset, "", 4096);
		gout.path(ics_path); 
		gout.glog(_log);
		gics_coder.clear(); 
		gics_coder.glog(_log); 
		gics_coder.path(ics_path);
		gics_coder.add_data("ID1", &_gics);
		gout.coder(&gics_coder);
		gout.run_write();

		return true;
	}

	bool t_gpreedit::_gnerate_clk_file(const t_gtime& beg, const t_gtime& end, double intv)
	{
		if (!_gnav) return false;

		//output new icsfile
		t_gfile gout;
		t_rinexc gclk(_gset, "", 4096);
		string ics_path = "./grt" + beg.str_gwkd() + ".clk";
		gout.path(ics_path);
		gout.glog(_log);
		gclk.clear();
		gclk.glog(_log);
		gclk.path(ics_path);
		gclk.add_data("ID1", dynamic_cast<t_gallprec*>(_gnav));
		gout.coder(&gclk);
		gout.run_write();
		return true;
	}

	bool t_gpreedit::_gnerate_tb_csh_file(const t_gtime& beg, const t_gtime& end, double intv)
	{

		return true;
	}

	bool t_gpreedit::_gnerate_rec_crd_xml(const t_gtime& beg, const t_gtime& end)
	{
		xml_document xdoc;

		auto root = xdoc.load_file(_gset->xml_name().c_str());
		auto config = xdoc.child("config");

		string obj = "SNX";
		// remove for case
		if (config.child("receiver")) config.remove_child("receiver");

		config.append_child("receiver");
		xml_node pars = config.child("receiver");

		pars.append_child(pugi::node_comment).set_value("The part for receiver settings. ");
		pars.append_child(pugi::node_comment).set_value("id    : for index, Rec name.     ");
		pars.append_child(pugi::node_comment).set_value("X     :  X pos.                   ");
		pars.append_child(pugi::node_comment).set_value("Y     :  Y pos.                   ");
		pars.append_child(pugi::node_comment).set_value("Z     :  Z pos.                   ");
		pars.append_child(pugi::node_comment).set_value("dX    : dX pos.                   ");
		pars.append_child(pugi::node_comment).set_value("dY    : dY pos.                   ");
		pars.append_child(pugi::node_comment).set_value("dZ    : dZ pos.                   ");

		// for STA ===============================================>
		set<string> rec_tmp = _gobj->objects_set(t_gdata::REC);
		int valid_sta_num = 0;
		for (auto tmp : _rec_list)
		{

			string rec = str2upper(tmp);
			shared_ptr<t_gobj> grec = _gobj->obj(rec);

			if (!grec) continue;

			t_gtriple rec_crd;
			t_gtriple rec_std;

			if (dynamic_cast<t_gsetproc*>(_gset)->crd_est() == CONSTRPAR::FIX)
			{
				grec->get_adjacent_crd(beg, 0.1, rec_crd, rec_std);
			}
			else
			{
				grec->get_recent_crd(beg, 500, rec_crd, rec_std);
			}
			obj = "SNX";
			
			if (rec_crd.zero())
			{
				_sppflt[rec] = new t_gsppflt(rec, _gset);
				_sppflt[rec]->setDAT(_gobs, _gnav);
				_sppflt[rec]->setOBJ(_gobj);

				// SPPFlt for kinematic station
				if (_sppflt[tmp]->processBatch(beg, end) < 0)
				{
					delete _sppflt[rec];
					_sppflt[rec] = nullptr;
					_log->comment(1, tmp, "Kinematic SPPFlt Fail !");
					continue;
				}

				rec_crd = _sppflt[rec]->getCrd(beg);

				delete _sppflt[rec];
				_sppflt[rec] = nullptr;

				rec_std = t_gtriple(100, 100, 100);
				obj = "SPP";
			}
			else
			{
				if (rec_std.crd(0) > 0.1)
				{
					rec_std = t_gtriple(100, 100, 100);
					obj = "RNO";
				}
			}

			if (rec_crd.zero())
			{
				continue;
			}

			// one station valid
			valid_sta_num++;

			xml_node par_sta = pars.append_child("rec");
			pugi::xml_attribute tmp_attri;
			tmp_attri = par_sta.append_attribute("id");
			tmp_attri.set_value(rec.c_str());

			string crd_x = format("%16.6lf", rec_crd[0]);
			tmp_attri = par_sta.append_attribute("X");
			tmp_attri.set_value(crd_x.c_str());

			string crd_y = format("%16.6lf", rec_crd[1]);
			tmp_attri = par_sta.append_attribute("Y");
			tmp_attri.set_value(crd_y.c_str());

			string crd_z = format("%16.6lf", rec_crd[2]);
			tmp_attri = par_sta.append_attribute("Z");
			tmp_attri.set_value(crd_z.c_str());

			crd_x = format("%16.6lf", rec_std[0]);
			tmp_attri = par_sta.append_attribute("dX");
			tmp_attri.set_value(crd_x.c_str());

			crd_y = format("%16.6lf", rec_std[1]);
			tmp_attri = par_sta.append_attribute("dY");
			tmp_attri.set_value(crd_y.c_str());

			crd_z = format("%16.6lf", rec_std[2]);
			tmp_attri = par_sta.append_attribute("dZ");
			tmp_attri.set_value(crd_z.c_str());

			tmp_attri = par_sta.append_attribute("obj");
			tmp_attri.set_value(obj.c_str());

		}

		xdoc.save_file(_gset->xml_name().c_str());

		return true;
	}

	bool t_gpreedit::_check_inp_files(set<string>& _inp_all)
	{
		return false;
	}
}


