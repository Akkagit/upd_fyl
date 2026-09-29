/**
*
* @verbatim
History
-1.0 ZhengHJ  2019-09-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		  gallrecover.cpp
* @brief	  Storage all recover info
*
* @author     ZhengHJ, Wuhan University
* @version	  1.0.0
* @date		  2019-09-25
*
*/

#include "gallrecover.h"

namespace great
{
	t_gallrecover::t_gallrecover() :
		t_gdata()
	{
		_type = t_gdata::ALLRECOVER;
	}


	t_gallrecover::~t_gallrecover()
	{
		for(t_grecover_data* data : _recover_data)
		{
			delete data;
			data = nullptr;
		}
		_time_equmap.clear();
		_time_parmap.clear();
		_site_sat_time_equmap.clear();
		_type_parmap.clear();
	}

	void t_gallrecover::add_recover_equation(const t_grecover_equation& recover_equ)
	{
		t_grecover_equation* equ_data = new t_grecover_equation(recover_equ);

		_add_common_data(equ_data);

		_recover_head.site_list.insert(recover_equ.site_name);                              
		_recover_head.sat_list.insert(recover_equ.sat_name);


		// default set && add
		if (_time_equmap.find(recover_equ.time) == _time_equmap.end()) {
			_time_equmap[recover_equ.time] = vector<t_grecover_equation*>();
		}
		// add equ map
		_time_equmap[recover_equ.time].push_back(equ_data);

		if (_site_sat_time_equmap.find(recover_equ.site_name) == _site_sat_time_equmap.end()) {
			_site_sat_time_equmap[recover_equ.site_name] = t_map_sat_equ();
		}
		t_map_sat_equ& temp_sat_equ = _site_sat_time_equmap[recover_equ.site_name];
		if (temp_sat_equ.find(recover_equ.sat_name) == temp_sat_equ.end()) {
			temp_sat_equ[recover_equ.sat_name] = t_map_time_equ();
		}
		t_map_time_equ& temp_time_equ = temp_sat_equ[recover_equ.sat_name];
		if (temp_time_equ.find(recover_equ.time) == temp_time_equ.end()) {
			temp_time_equ[recover_equ.time] = vector<t_grecover_equation*>();
		}
		temp_time_equ[recover_equ.time].push_back(equ_data);

	}

	void t_gallrecover::add_recover_par(const t_grecover_par & recover_par)
	{
		t_grecover_par* par_data = new t_grecover_par(recover_par);

		_add_common_data(par_data);

		if (_time_parmap.find(recover_par.get_recover_time()) == _time_parmap.end())
		{
			_time_parmap[recover_par.get_recover_time()] = vector<t_grecover_par*>();
			_time_parmap_pce[recover_par.get_recover_time()] = vector<t_grecover_par>();
		}
		_time_parmap[recover_par.get_recover_time()].push_back(par_data);
		_time_parmap_pce[recover_par.get_recover_time()].push_back(*par_data);
		
		if (_type_parmap.find(recover_par.par.parType) == _type_parmap.end())
		{
			_type_parmap[recover_par.par.parType]= vector<t_grecover_par*>();
		}
		_type_parmap[recover_par.par.parType].push_back(par_data);
	}

	void t_gallrecover::get_clkdata(t_gallprec& clkdata, t_gallprec::clk_type type)
	{
		set<t_gallprec::clk_type> type_list;
		if (type != t_gallprec::UNDEF) {
			type_list.insert(type);
		}
		else {
			type_list.insert(t_gallprec::AR);
			type_list.insert(t_gallprec::AS);
		}
		
		for (auto iter = _time_parmap.begin(); iter != _time_parmap.end(); iter++)
		{
			t_gtime epoch = iter->first;
			for (auto par : iter->second) {
				string obj = "";
				if (par->par.parType == t_gpar::CLK && type_list.find(t_gallprec::AR) != type_list.end() ){
					obj = par->par.site;
				}
				else if (par->par.parType == t_gpar::CLK_SAT && type_list.find(t_gallprec::AS) != type_list.end() ) {
					obj = par->par.prn;
				}
				else {
					continue;
				}

				// jdhuang
				if (par->correct_value == 0.0) {
					continue;
				}

				double clk[3] = { 0.0,0.0,0.0 };
				double var[3] = { 0.0,0.0,0.0 };

				clk[0] = (par->par.value() + par->correct_value) / CLIGHT;
				clkdata.addclk(obj, epoch, clk, var);
			}
		}
	}

	void t_gallrecover::get_rbdata(t_gallrb& rbdata)
	{

		for (auto iter = _time_parmap.begin(); iter != _time_parmap.end(); iter++)
		{
			t_gtime epoch = iter->first;
			for (auto par : iter->second) {
				string obj = "";
				string obj_sat = "";
				if (par->par.parType == t_gpar::RB) {
					obj = par->par.site;
					obj_sat = par->par.prn;
				}
				else {
					continue;
				}

				double rb = 0.0;

				rb = par->par.value() + par->correct_value;

				rbdata.addrb(obj, obj_sat, rb);
			}
		}
	}

	void t_gallrecover::get_eopdata(t_galleop& eopdata)
	{
		map<t_gtime, eop> eop_all;
		eop eop_one;
		int tab = 0;

		for (auto iter = _time_parmap.begin(); iter != _time_parmap.end(); iter++)
		{
			t_gtime epoch = iter->first;
			for (auto par : iter->second) {
				if (par->par.parType == t_gpar::XPOLE) {
					//eop_one.xpole = par->par.value() + par->correct_value;
					eop_one.xpole = par->correct_value;
					tab = 1;
				}
				else if (par->par.parType == t_gpar::YPOLE) {
					//eop_one.ypole = par->par.value() + par->correct_value;
					eop_one.ypole = par->correct_value;
				}
				else if (par->par.parType == t_gpar::DXPOLE) {
					//eop_one.dxpole = par->par.value() + par->correct_value;
					eop_one.dxpole = par->correct_value;
				}
				else if (par->par.parType == t_gpar::DYPOLE) {
					//eop_one.dypole = par->par.value() + par->correct_value;
					eop_one.dypole = par->correct_value;
				}
				else if (par->par.parType == t_gpar::UT1) {
					//eop_one.ut1 = par->par.value() + par->correct_value;
					eop_one.ut1 = par->correct_value;
				}
				else if (par->par.parType == t_gpar::DUT1) {
					//eop_one.dut1 = par->par.value() + par->correct_value;
					eop_one.dut1 = par->correct_value;
				}
				else {
					continue;
				}

			}
			if (tab == 1)
			{
				eop_all[iter->first] = eop_one;
				eop_one.clear_value();
				tab = 0;
			}

		}

		for (auto iter = eop_all.begin(); iter != eop_all.end(); iter++)
		{
			eopdata.addeop(iter->first, iter->second);
		}

	}

	void t_gallrecover::get_poledata(t_gpoleut1& poledata)
	{
		map<string, double> pole_one;
		int tab = 0;

		for (auto iter = _time_parmap.begin(); iter != _time_parmap.end(); iter++)
		{
			t_gtime epoch = iter->first;
			for (auto par : iter->second) {
				if (par->par.parType == t_gpar::XPOLE) {
					//eop_one.xpole = par->par.value() + par->correct_value;
					pole_one["XPOLE"] = par->correct_value + par->par.value();
					tab = 1;
				}
				else if (par->par.parType == t_gpar::YPOLE) {
					//eop_one.ypole = par->par.value() + par->correct_value;
					pole_one["YPOLE"] = par->correct_value+ par->par.value();
				}
				else if (par->par.parType == t_gpar::DXPOLE) {
					//eop_one.dxpole = par->par.value() + par->correct_value;
					pole_one["DXPOLE"] = par->correct_value+ par->par.value();
				}
				else if (par->par.parType == t_gpar::DYPOLE) {
					//eop_one.dypole = par->par.value() + par->correct_value;
					pole_one["DYPOLE"] = par->correct_value+ par->par.value();
				}
				else if (par->par.parType == t_gpar::UT1) {
					//eop_one.ut1 = par->par.value() + par->correct_value;
					pole_one["UT1"] = par->correct_value+ par->par.value();
				}
				else if (par->par.parType == t_gpar::DUT1) {
					//eop_one.dut1 = par->par.value() + par->correct_value;
					pole_one["DUT1"] = par->correct_value+ par->par.value();
				}
				else {
					continue;
				}

			}
			if (tab == 1)
			{
				poledata.setEopData(iter->first, pole_one, "UT1R", 1.0);
				pole_one.clear();
				tab = 0;
			}

		}

	}

	void t_gallrecover::get_stadata(t_gallsta& stadata)
	{
		map<string, double> dx;
		map<string, double> dy;
		map<string, double> dz;

		for (auto iter = _time_parmap.begin(); iter != _time_parmap.end(); iter++)
		{
			t_gtime epoch = iter->first;
			for (auto par : iter->second) {
				string obj = "";
				if (par->par.parType == t_gpar::CRD_X) {
					obj = par->par.site;
					dx[obj] = par->correct_value;
				}
				else if (par->par.parType == t_gpar::CRD_Y) {
					obj = par->par.site;
					dy[obj] = par->correct_value;
				}
				else if (par->par.parType == t_gpar::CRD_Z) {
					obj = par->par.site;
					dz[obj] = par->correct_value;
				}
				else {
					continue;
				}

			}
		}

		for (auto iter = dx.begin(); iter != dx.end(); iter++)
		{
			stadata.addsta(iter->first, iter->second, dy[iter->first], dz[iter->first]);
		}

	}

	void t_gallrecover::get_clk13data(t_gallprec& clkdata, t_gallprec::clk_type type)
	{
		set<t_gallprec::clk_type> type_list;
		if (type != t_gallprec::UNDEF) {
			type_list.insert(type);
		}
		else {
			type_list.insert(t_gallprec::AR);
			type_list.insert(t_gallprec::AS);
		}

		for (auto iter = _time_parmap.begin(); iter != _time_parmap.end(); iter++)
		{
			map<string, double> ifcb_data;
			t_gtime epoch = iter->first;
			int i = 0;
			sort(_time_parmap_pce[epoch].begin(), _time_parmap_pce[epoch].end());
			for (auto par : _time_parmap_pce[epoch]) {
				string obj = "";
				if (par.par.parType == t_gpar::CLK && type_list.find(t_gallprec::AR) != type_list.end()) {
					obj = par.par.site;
				}
				//else if ((par->par.parType == t_gpar::G15CLK_SAT || par->par.parType == t_gpar::E15CLK_SAT || 
				//	par->par.parType == t_gpar::C15CLK_SAT || par->par.parType == t_gpar::J15CLK_SAT) 
				//	&& type_list.find(t_gallprec::AS) != type_list.end()) {
				//	obj = par->par.prn;
				//}
				else if (par.par.parType == t_gpar::CLK_SAT && type_list.find(t_gallprec::AS) != type_list.end()) {
					obj = par.par.prn;
				}
				else if (par.par.parType == t_gpar::IFCB_F3 && type_list.find(t_gallprec::AS) != type_list.end()) {
					obj = par.par.prn;
				}
				else {
					continue;
				}

				if (par.correct_value == 0.0) {
					continue;
				}

				double clk[3] = { 0.0,0.0,0.0 };
				double clk_tri[4] = { 0.0,0.0,0.0,0.0 };//add by xiongyun
				double var[3] = { 0.0,0.0,0.0 };

				// xiongyun change for IFCB

				if (par.par.parType == t_gpar::IFCB_F3)
				{
					ifcb_data[par.par.prn] = (par.par.value() + par.correct_value) / CLIGHT;
					continue;
				}

				clk[0] = (par.par.value() + par.correct_value) / CLIGHT;
				clk_tri[0] = clk[0];
				if (ifcb_data.find(par.par.prn) != ifcb_data.end())
				{
					clk_tri[3] = ifcb_data[par.par.prn];
				}
				if (ifcb_data.size() == 0)
				{
					clkdata.addclk(obj, epoch, clk, var);
				}
				else
				{
					clkdata.addclk_tri(obj, epoch, clk_tri, var);
				}
			}
		}
	}

	void t_gallrecover::get_iondata(vector<t_tuple_ion>& ion_data)
	{
		for (auto iter = _time_parmap.begin(); iter != _time_parmap.end(); iter++)
		{
			t_gtime epoch = iter->first;
			for (auto par : iter->second) 
			{
				string type = "";
				if (par->par.parType == t_gpar::SION)
				{
					type = "SION";
				}
				else if (par->par.parType == t_gpar::VION)
				{
					type = "VION";
				}
				else
				{
					continue;
				}

				if (double_eq(par->correct_value, 0.0)) continue;

				double value = 0.0;
				double sigma = 0.5;
				
				value = (par->par.value() + par->correct_value);
				auto tmp = make_tuple(type, par->par.site, par->par.prn,
					value, sigma, par->par.beg, par->par.end);
				ion_data.push_back(tmp);
			}
		}
	}

	vector<t_grecover_equation> t_gallrecover::get_first_phase_recover_equation(string site, string sat, string freq)
	{

		map<t_gobscombtype, int> phase_order = 
		{
			{t_gobscombtype("LC"),10},
			{t_gobscombtype("L1"),9},
			{t_gobscombtype("L2"),8},
			{t_gobscombtype("L3"),7},
			{t_gobscombtype("L4"),6},
			{t_gobscombtype("L5"),5},
			{t_gobscombtype("L6"),4},
			{t_gobscombtype("L7"),3},
			{t_gobscombtype("L8"),2},
			{t_gobscombtype("L9"),1}
		};

		// jdhuang : 如果规定了freq。则把该频率的权限设为最优先
		if (phase_order.find(t_gobscombtype(freq)) != phase_order.end())
		{
			phase_order[t_gobscombtype(freq)] = 11;
		}

		vector<t_grecover_equation> first_phase_recover_equation;

		for(auto time_equ : _site_sat_time_equmap[site][sat]) 
		{		
			t_grecover_equation first_equ(time_equ.first,site,sat);
			int max_order = 0;
			for(auto equ : time_equ.second)
			{
				//add for slr
				if (equ->obstype.is_SLR())
				{
					first_equ = *equ;
					//first_phase_recover_equation.push_back(first_equ);
				}
				if (!equ->obstype.is_phase()) 
				{
					continue;
				}
				if (phase_order[equ->obstype] > max_order)
				{
					max_order = phase_order[equ->obstype];
					first_equ = *equ;
				}
			}
			first_phase_recover_equation.push_back(first_equ);

		}
		return first_phase_recover_equation;
	}

	vector<t_grecover_par> t_gallrecover::get_recover_par(t_gpar::t_type parType)
	{
		vector<t_grecover_par> tmp;
		for (auto type_par : _type_parmap[parType])
		{
			t_grecover_par par(type_par->par, type_par->correct_value);
			tmp.push_back(par);
		}
		return tmp;
	}


	t_gtime t_gallrecover::get_beg_time() const
	{
		return _recover_head.get_beg_time();
	}

	t_gtime t_gallrecover::get_end_time() const
	{
		return _recover_head.get_end_time();
	}

	t_gtime t_gallrecover::get_equ_beg_time() const
	{
		return _time_equmap.begin()->first;
	}

	double t_gallrecover::get_interval() const
	{
		return _recover_head.interval;
	}

	void t_gallrecover::set_interval(double intv)
	{
		_recover_head.interval = intv;
	}

	const vector<t_grecover_data*>& t_gallrecover::get_all_recover_data() const
	{
		return _recover_data;
	}

	double t_gallrecover::get_recover_data_value(t_gpar::t_type parType)
	{
		vector<t_grecover_data*>::iterator it;
		for (it = _recover_data.begin(); it != _recover_data.end(); it++)
		{
			double correct_value;
			double par_value;
			if (parType == dynamic_cast<t_grecover_par*>(*it)->par.parType)
			{
				par_value = dynamic_cast<t_grecover_par*>(*it)->par.value();
				correct_value = dynamic_cast<t_grecover_par*>(*it)->correct_value;
				if (double_eq(par_value, 0.0) || double_eq(correct_value, 0.0)) return 0;
				return par_value + correct_value;
			}
		}
		return 0;
	}

	const t_map_time_par& t_gallrecover::get_map_time_par() const
	{
		return _time_parmap;
	}

	const t_map_time_equ& t_gallrecover::get_map_time_equ() const
	{
		return _time_equmap;
	}
	const t_map_site_equ& t_gallrecover::get_map_site_equ() const
	{
		return _site_sat_time_equmap;
	}

	set<string> t_gallrecover::get_sat_list() const
	{
		return _recover_head.sat_list;
	}

	set<string> t_gallrecover::get_site_list() const
	{
		return _recover_head.site_list;
	}

	set<t_gtime> t_gallrecover::get_time_list() const
	{
		return _recover_head.time_list;
	}

	vector<t_gtime> t_gallrecover::get_all_obs_time_list() const
	{
		t_gtime beg = _time_equmap.begin()->first;
		t_gtime end = _time_equmap.rbegin()->first;

		int nepoch = floor((end - beg) / _recover_head.interval) + 1;
		vector<t_gtime> ans;
		for (int i = 0; i < nepoch; i++) {
			ans.push_back(beg + _recover_head.interval * i);
		}
		return ans;
	}

	t_gallpar t_gallrecover::get_all_pars()
	{
		t_gallpar tmp;

		t_gtime beg  = get_beg_time();
		t_gtime end  = get_end_time();
		double  intv = get_interval();
		while (beg <= end)
		{
			auto par_end = _time_parmap[beg].end();
			auto par_crt = _time_parmap[beg].begin();
			for (par_crt; par_crt != par_end; par_crt++)
			{
				// 加上改正数
				auto recover_par = (*par_crt)->par;
				recover_par.value(recover_par.value() + (*par_crt)->correct_value);
				tmp.addParam(recover_par);
			}
		}
		
		return  tmp;
	}

	map<t_gtime, int> t_gallrecover::get_sat_number_per_epo() const
	{
		map<t_gtime, int> nsat;
		for (auto it : _time_equmap)
		{
			nsat[it.first] = it.second.size() / 2;
		}
		return nsat;
	}

	void t_gallrecover::_add_common_data(t_grecover_data * data)
	{
		_recover_data.push_back(data);
		_recover_head.time_list.insert(data->get_recover_time());
	}

}