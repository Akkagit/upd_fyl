/**
*
* @verbatim
History
-1.0 glfeng  2019-03-05  creat the file.
-1.2 BoWong  2019-04-09  Adding Doxygen Style Code Remarks.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambfixDd.cpp
* @brief	integer estimation with the LAMBDA method.
*
* @author   glfeng BoWong , Wuhan University
* @version	1.0.0
* @date		2019-04-09
*
*/

#include "gambfix/gambfixDd.h"
#include "gambfix/gallamb.h"
#include "gcoders/ambcon.h"
#include "gset/gsetout.h"
#include "gio/gfile.h"

#include "math.h"
#include <string>
#include <stdio.h>
using namespace std;

namespace great
{
	t_gambfixDd::t_gambfixDd()
	{

	}

	t_gambfixDd::t_gambfixDd(t_glog* log, t_gsetbase* set) :
		t_gambfix(log, set)
	{
		cerr << "================> constructor of t_gambfixDd" << endl;


		// add data to gall ambs
		this->_gambs.add_set(_gset);
		this->_gambs.add_log(_glog);
	}


	t_gambfixDd::~t_gambfixDd()
	{

	}

	bool t_gambfixDd::ProcessBatch()
	{
		string func_id = "ProcessBatch";
		// 判断是否有数据
		if (!this->_gpar || !this->_gset || !this->_gbias)
		{
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, this->_class_id, func_id, "check the gpar || gset || gbias");
			return false;
		}

		// add data to gambs
		// 设置格式，加入数据到all amb
		this->_gambs.add_data(_gpar, _gobs, _gbias, _gobj);

		// get baselines for data or setting
		t_gambfixDd::_screen_rec_and_sat(_gpar, _gset, _gobj);
		this->_gambs.add_WL_bdeci(this->_evaluator_WL);
		this->_gambs.add_NL_bdeci(this->_evaluator_NL);

		// 组单差模糊度
		DD_MODEL dd_mode = dynamic_cast<t_gsetamb*>(_gset)->dd_mode();
		for (const auto& ID : this->_gambs.IDs())
		{
			
			// 无论是L1+L2还是L1+L3,在UCUD中，都是L1和L2的基线完全相同的
			if (ID == t_gamb::AMB_ID::ALL_L1)    continue;
			//if (ID == t_gamb::AMB_ID::ALL_L1_L3) continue;

			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "get Sd ambs");
			this->_gambs.cmb_amb(t_gallambs::AMB_COMB_OPTION::CMB_SAT_SD, ID, this->_gSatBl.baselines(), this->_min_common_time);

			// 组双差模糊度
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "get Dd ambs");
			this->_gambs.cmb_amb(t_gallambs::AMB_COMB_OPTION::CMB_REC_DD, ID, this->_gRecBl.baselines(), this->_min_common_time);
		}
		// 判断双差模糊度是否固定
		if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "check Dd ambs");
		
#ifdef DEBUG_AMBFIX
		this->_gambs->cout_ambs();
#endif // DEBUG_AMBFIX

		// 获取独立星间基线

#ifdef DEBUG_AMBFIX
		// 以测站为标准，输出所有双差模糊度
		auto baselines = this->_rec_baselines->baselines();
		for (auto iter : baselines)
		{
			this->_gambs->cout_DD_ambs(iter);
		}
#endif // DEBUG_AMBFIX

		// write ambcon
		if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "write Dd ambs ambcon");
		this->_write_ambcon();

		return true;
	}

	bool t_gambfixDd::ProcessBatch_ALL()
	{
		string func_id = "ProcessBatch";
		// 判断是否有数据
		if (!this->_gpar || !this->_gset || !this->_gobj)
		{
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, this->_class_id, func_id, "check the gpar || gset");
			return false;
		}

		// add data to gambs
		// 设置格式，加入数据到all amb
		this->_gambs.add_data(_gpar, _gobs, nullptr);

		// get baselines for data or setting
		t_gambfixDd::_screen_rec_and_sat(_gpar, _gset, _gobj);

		// 组单差模糊度
		for (const auto& ID : this->_gambs.IDs())
		{

			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "get Sd ambs");
			this->_gambs.cmb_amb(t_gallambs::AMB_COMB_OPTION::CMB_SAT_SD, ID, this->_gSatBl.baselines(), this->_min_common_time);

			// 组双差模糊度
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "get Dd ambs");
			this->_gambs.cmb_amb(t_gallambs::AMB_COMB_OPTION::CMB_REC_DD, ID, this->_gRecBl.baselines(), this->_min_common_time);
		}

		// write ambcon
		if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "write Dd ambs ambcon");
		this->_write_ambcon();

		return true;

	}



	bool t_gambfixDd::_screen_rec_and_sat(t_gallpar* pars, t_gsetbase* set, t_gallobj* obj)
	{
		if (!pars && !set && !obj) return false;
		// 获取交集
		try
		{
			auto par_rec_list = pars->amb_recs();
			auto par_sat_list = pars->amb_prns();
			auto xml_rec_list = dynamic_cast<t_gsetgen*>(set)->recs();
			auto xml_sat_list = dynamic_cast<t_gsetgnss*>(set)->sat();

			set_intersection(par_rec_list.begin(), par_rec_list.end(), xml_rec_list.begin(), xml_rec_list.end(), inserter(_rec_list, _rec_list.begin()));//交集
			set_intersection(par_sat_list.begin(), par_sat_list.end(), xml_sat_list.begin(), xml_sat_list.end(), inserter(_sat_list, _sat_list.begin()));//交集
		}
		catch (...)
		{
			return false;
		}

		// 获取模糊度的测站和卫星列表
		const auto max = pars->parNumber();
		for (unsigned int i = 0; i < max; i++)
		{
			_gRecBl.add_rec((*pars)[i].site, (*pars)[i].prn);
			_gSatBl.add_sat((*pars)[i].prn);
		}

		// 再判断下是不是测站，卫星，系统列表都有
		if (_rec_list.empty() || _sat_list.empty())
		{
			return false;
		}

		// 放入数据
		_gRecBl.add_cut(_beg, _max_baseline_length);
		_gRecBl.add_gobj(obj);
		_gRecBl.add_gset(set);
		_gRecBl.add_glog(this->_glog);

		return true;

	}

	void t_gambfixDd::_write_ambcon()
	{
		//输出约束文件
		t_gambcon _gambcon;
		auto ddambs = this->_gambs.get_Dds();
		
		int num_dds = 0;
		for (auto& sys : ddambs)
		{
			for (auto& bl : sys.second)
			{
				for (auto& amb : bl.second)
				{
					num_dds++;
					_gambcon.add_ddamb(amb);
				}
			}
		}

		_gambcon.set_header(_beg, _intv, num_dds, _sat_list, _rec_list);
		_gambcon.all_NL_num = this->_gambs.all_NL_num;
		_gambcon.all_WL_num = this->_gambs.all_WL_num;
		_gambcon.fix_NL_num = this->_gambs.fix_NL_num;
		_gambcon.fix_WL_num = this->_gambs.fix_WL_num;

		t_gfile gout;
		t_ambcon gamb_coder(_gset, "", 4096);
		string path = (dynamic_cast<t_gsetout*>(_gset)->outputs("ambcon"));
		if (path.empty())
		{
			path = "grt$(date).con";
		}
		substitute(path, "$(date)", _beg.str_yyyydoy(), false);
		gout.path(path); 
		gout.glog(this->_glog);
		gamb_coder.clear(); 
		gamb_coder.glog(this->_glog); 
		gamb_coder.path(path);
		gamb_coder.add_data("ID1", &_gambcon);
		gamb_coder.add_data("OBJ", this->_gobj);
		gout.coder(&gamb_coder);
		gout.run_write();
	}
}



