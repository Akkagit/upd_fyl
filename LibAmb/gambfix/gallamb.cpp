/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gallamb.cpp
* @brief	function for ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#include "gambfix/gambOw.h"
#include "gset/gsetamb.h"
#include "gset/gsetinp.h"
#include "gutils/gtime.h"
#include "gallamb.h"
#include "gall/gallbias.h"
#include "gambfix/gambcommon.h"
#include "gmodels/gprecisemodel.h"
#include "gutils/gstring.h"
#include "gutils/ginfolog.h"

#include "gio/gfile.h"
#include "gcoders/rinexo.h"

#include <thread>
namespace great
{
	t_gallambs::t_gallambs()
	{
	}

	t_gallambs::t_gallambs(t_gsetbase* gset, t_glog* log)
	{
		this->_glog = log;
		this->_gset = gset;
		_dd_mode = dynamic_cast<t_gsetamb*>(gset)->dd_mode();
		//_dcb_model = dynamic_cast<t_gsetproc*>(gset)->dcb_model();
		//_freq_model = dynamic_cast<t_gsetproc*>(gset)->freq_model();

		_frequency = dynamic_cast<t_gsetproc*>(gset)->frequency();
		_band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::GPS);
		_band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::GAL);
		_band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::GLO);
		_band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::BDS);
		_band_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::QZS);

		_freq_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::GPS);
		_freq_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::GAL);
		_freq_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::GLO);
		_freq_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::BDS);
		_freq_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::QZS);
#ifdef DEBUG_AMBFIX
		if (!_Ow_ambs.empty())
		{
			cout << endl;
			cout << "========================> The Ow amb in constuctor" << endl;
			int tmp1 = 0;
			for (auto iter : _Ow_ambs)
			{
				tmp1++;
				cout << setw(6) << "OW_AMB"
					<< setw(8) << tmp1
					<< setw(6) << iter.second.rec
					<< setw(4) << iter.second.sat
					<< setw(4) << iter.second.fixed
					<< setw(4) << iter.second.band()
					<< setw(16) << setprecision(8) << iter.second.LC_value
					<< setw(16) << setprecision(8) << iter.second.LC_sigma
					<< setw(16) << setprecision(8) << iter.second.freq()
					<< setw(16) << setprecision(8) << iter.second.wavelength()
					<< setw(16) << iter.second.beg.str_ymdhms("Beg", false)
					<< setw(16) << iter.second.end.str_ymdhms("End", false)
					<< endl;
			}
			cout << endl;
		}

		if (!_Cb_ambs.empty())
		{
			auto key_cb = make_tuple(t_gamb_Cb::COMB_TYPE::COMB_IF, GOBS_LC::LC_L1, GOBS_LC::LC_L2);
			cout << endl;
			cout << "========================> The Cb amb in constuctor" << endl;
			int tmp2 = 0;
			for (auto iter : _Cb_ambs)
			{
				tmp2++;
				cout << setw(6) << "Cb_AMB"
					<< setw(8) << tmp2
					<< setw(6) << iter.second.rec
					<< setw(4) << iter.second.sat
					<< setw(4) << iter.second.fixed(key_cb)
					<< setw(16) << setprecision(8) << iter.second.freq(GOBS_LC::LC_L1)
					<< setw(16) << setprecision(8) << iter.second.freq(GOBS_LC::LC_L2)
					<< setw(16) << setprecision(8) << iter.second.value(key_cb)
					<< setw(16) << setprecision(8) << iter.second.sigma(key_cb)
					<< setw(16) << setprecision(8) << iter.second.freq(key_cb)
					<< setw(16) << setprecision(8) << iter.second.length(key_cb)
					<< setw(30) << iter.second.beg.str_ymdhms("Beg : ", false)
					<< setw(30) << iter.second.end.str_ymdhms("End : ", false)
					<< endl;
			}
			cout << endl;
		}

		if (!_Sd_ambs.empty())
		{

		}

		if (!_Dd_ambs.empty())
		{

		}
#endif

	}

	t_gallambs::~t_gallambs()
	{
		if (_gprecise)
		{
			delete _gprecise;
			_gprecise = nullptr;

		}
	}

	void t_gallambs::add_set(t_gsetbase* gset)
	{
		_gset = gset;
		_dd_mode = dynamic_cast<t_gsetamb*>(gset)->dd_mode();
		//_dcb_model = dynamic_cast<t_gsetproc*>(gset)->dcb_model();
		//_freq_model = dynamic_cast<t_gsetproc*>(gset)->freq_model();

		_band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::GPS);
		_band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::GAL);
		_band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::GLO);
		_band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::BDS);
		_band_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::QZS);
														  
		_freq_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::GPS);
		_freq_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::GAL);
		_freq_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::GLO);
		_freq_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::BDS);
		_freq_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::QZS);
	}

	void t_gallambs::add_log(t_glog* log)
	{
		_glog = log;
	}

	void t_gallambs::add_data(t_gallpar* pars, t_gallobs* obs, t_gallbias* allbias,t_gallobj* obj)
	{
		string func_id = "t_gallambs::add_data";
		this->_gobs = obs;
		this->_gbias = allbias;
		this->_gpars = pars;
		this->_gobj = obj;
		// ADD DATA
		t_gallproc* data = new t_gallproc();
		data->Add_Data(this->_gbias);
		_gprecise = new t_gprecisemodel(data, _glog, _gset);


		// 从 par获取模糊度
		if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "_get_ows_table_from_pars");
		//if (_dd_mode != DD_MODEL::RAW_CB) t_gallambs::_get_ows_table_from_pars();
		t_gallambs::_get_ows_table_from_pars();

		// 从 obs获取MW。只有无电离层双差有用
		if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "_get_obs_table_from_gobs");
		if (_dd_mode != DD_MODEL::RAW_CB) t_gallambs::_get_obs_table_from_gobs();
	}

	void t_gallambs::add_data(t_gallpar* pars, t_gallobs* obs, t_gallbias* allbias, t_gupd* upd)
	{
		string func_id = "t_gallambs::add_data";
		this->_gobs = obs;
		this->_gbias = allbias;
		this->_gpars = pars;
		this->_gupd = upd;
		// ADD DATA
		t_gallproc* data = new t_gallproc();
		data->Add_Data(this->_gbias);
		_gprecise = new t_gprecisemodel(data, _glog, _gset);
		if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "_get_ows_table_from_pars");
		t_gallambs::_get_ows_table_from_pars();
		if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "_get_obs_table_from_gobs");
		t_gallambs::_get_obs_table_from_gobs();
	}

	void t_gallambs::add_id(const t_gamb::AMB_ID& ID)
	{
		_IDs.insert(ID);
	}

	void t_gallambs::add_WL_bdeci(t_gbdeci* WL)
	{
		_evaluator_MW = WL;
	}

	void t_gallambs::add_NL_bdeci(t_gbdeci* NL)
	{
		_evaluator_NL = NL;
	}

	t_gallambs::t_d t_gallambs::get_Sds()
	{
		return this->_nSd_ambs;
	}
	t_gallambs::t_dd t_gallambs::get_Dds()
	{
		return this->_nDd_ambs;
	}

	bool t_gallambs::cmb_amb(const AMB_COMB_OPTION& opt, const t_gamb::AMB_ID& ID, const vector<pair<string, string>>& baselines, double shortest_arc)
	{
		string func_id = "combine_Ow_ambs";

		// DD或者SD
		if (baselines.empty())
		{
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "the baselines is empty!");
			return false;
		}

		// 遍历所有基线
		for (const auto& iter : baselines)
		{
			if (_dd_mode == DD_MODEL::IF_CB_WN)
			{
				if (ID == t_gamb::AMB_ID::IF_L1_L2 || ID == t_gamb::AMB_ID::IF_L1_L3)
				{
					if (opt == AMB_COMB_OPTION::CMB_SAT_SD) this->_cmb_intersat_Sd_amb(iter, ID, shortest_arc);
					if (opt == AMB_COMB_OPTION::CMB_REC_DD) this->_cmb_interrec_Dd_amb(iter, ID, shortest_arc);
				}
			}
			else if (_dd_mode == DD_MODEL::RAW_CB_WN || _dd_mode == DD_MODEL::RAW_CB)
			{
				if (ID == t_gamb::AMB_ID::ALL_L1) continue;
				if (opt == AMB_COMB_OPTION::CMB_SAT_SD) this->_cmb_intersat_Sd_amb(iter, t_gamb::AMB_ID::ALL_L1, ID, shortest_arc);
				if (opt == AMB_COMB_OPTION::CMB_REC_DD) this->_cmb_interrec_Dd_amb(iter, t_gamb::AMB_ID::ALL_L1, ID, shortest_arc);
			}	
			else continue;
		}

		if (opt == AMB_COMB_OPTION::CMB_REC_DD)
		{
			// get Dd value
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "_get_Dds_value");

			if (_dd_mode == DD_MODEL::IF_CB_WN)
			{
				if (ID == t_gamb::AMB_ID::IF_L1_L2 || ID == t_gamb::AMB_ID::IF_L1_L3)
				{
					_get_Dds_value(ID);
					if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "_select_Dds_value");
					_select_Dds_value(ID);
					if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "_get_dependent_Dd");
					_get_dependent_IF_Dd(ID);
					_get_fix_rate(ID);
				}
			}
			else if (_dd_mode == DD_MODEL::RAW_CB_WN || _dd_mode == DD_MODEL::RAW_CB)
			{
				if (ID == t_gamb::AMB_ID::ALL_L2)
				{
					_get_Dds_value(t_gamb::AMB_ID::ALL_L1_L2);
					_select_Dds_value(t_gamb::AMB_ID::ALL_L1_L2);
					_get_dependent_1X_Dd(t_gamb::AMB_ID::ALL_L1_L2);
					_get_fix_rate(t_gamb::AMB_ID::ALL_L1_L2);

				}
				if (ID == t_gamb::AMB_ID::ALL_L3)
				{
					_get_Dds_value(t_gamb::AMB_ID::ALL_L1_L3);
					_select_Dds_value(t_gamb::AMB_ID::ALL_L1_L3);
					_get_dependent_1X_Dd(t_gamb::AMB_ID::ALL_L1_L3);
					_get_fix_rate(t_gamb::AMB_ID::ALL_L1_L3);
				}
			}
			else
			{
				return true;
			}
			
		}
		return true;
	}
	
	bool t_gallambs::cmb_amb(const AMB_COMB_OPTION& opt, const t_gamb::AMB_ID& ID, const vector<pair<string, string>>& baselines, 
		double shortest_arc, t_gsetamb::UPD_MODE upd_mode) {
		string func_id = "cmb_amb";
		if (baselines.empty())
		{
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "the baselines is empty!");
			return false;
		}
		for (const auto& iter : baselines)
		{
			if (opt == AMB_COMB_OPTION::CMB_SAT_SD)   this->_cmb_intersat_Sd_amb(iter, ID, shortest_arc);
			else continue;
		}
		if (opt == AMB_COMB_OPTION::CMB_SAT_SD)
		{
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "_get_Sds_value");
			_get_Sds_value(upd_mode);
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "_select_Sds_value");
			_select_Sds_value();
			_get_dependent_IF_Sd2(t_gamb::AMB_ID::IF_L1_L2);
			_get_SDfix_rate(ID);
		}
		return true;
	}

	bool t_gallambs::_cmb_intersat_Sd_amb(const pair<string, string>& sat_baseline, const t_gamb::AMB_ID& ID, double shortest_arc)
	{
		string func_id = "_cmb_intersat_Sd_amb";

#ifdef DEBUG_AMBFIX
		// 输出调试信息
		if (this->_glog)
		{
			this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "combine Sd amb for : begin combine Sd ambs with one rec and two sat");
		}
#endif // DEBUG_AMBFIX

		// 检查是一个系统
		if (sat_baseline.first.substr(0, 1) != sat_baseline.second.substr(0, 1))
		{
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_WARN, this->_class_id, func_id, "not the same system " +
				sat_baseline.first + " " + sat_baseline.second);
			return false;
		}

		// 检查有没有Ow数据
		// 判断系统有没有OW，没有OW
		t_dd_sys_index sys_key = make_pair(t_gsys::sat2gsys(sat_baseline.first), ID);
		if (this->_nOw_ambs.count(sys_key) < 1 || this->_nOw_ambs[sys_key].empty())
		{
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_WARN, this->_class_id, func_id, "no Ow " +
				sat_baseline.first + " " + sat_baseline.second);
			return false;
		}

		// 遍历测站
		const auto& sys_Ow = this->_nOw_ambs[sys_key];
		for (const auto& iter_rec : this->_rec_list)
		{
			// 检查测站有没有Ow
			if (sys_Ow.count(iter_rec) < 1)
			{
				if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_WARN, this->_class_id, func_id, "no Ow " +
					iter_rec);
				continue;
			}

			// 找到这个测站的所有模糊度
			const auto& Ow_Ambs = this->_nOw_ambs[sys_key][iter_rec];

			// 找到两个卫星的模糊度
			vector<int> satA;
			vector<int> satB;
			for (int amb_idx = 0; amb_idx < Ow_Ambs.size(); amb_idx++)
			{
				if (Ow_Ambs[amb_idx].sat == sat_baseline.first)  satA.push_back(amb_idx);
				if (Ow_Ambs[amb_idx].sat == sat_baseline.second) satB.push_back(amb_idx);
			}

			// 任意一个没有模糊度，直接剔除
			if (satA.empty() || satB.empty())
			{
				if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_WARN, this->_class_id, func_id, "not ambs for " +
					sat_baseline.first + " or " + sat_baseline.second);
				continue;
			}

			// 遍历两颗卫星的模糊度，组成基线
			// 原则上 
			// 1. 不同卫星系统不差分
			// 2. 共同弧度小于阈值不差分
			for (const auto& iter_satA : satA)
			{
				for (const auto& iter_satB : satB)
				{
					// 再检查一次测站
					if (Ow_Ambs[iter_satA].rec != Ow_Ambs[iter_satB].rec) continue;

					// 再判断有没有交集
					const auto& amb_1st = Ow_Ambs[iter_satA];
					const auto& amb_2nd = Ow_Ambs[iter_satB];

					// 没有任何交集
					auto MAX_BEG = MAX_TIME(amb_1st.beg, amb_2nd.beg);
					auto MIN_END = MIN_TIME(amb_1st.end, amb_2nd.end);
					if (MAX_BEG >= MIN_END) continue;

					// 有交集，再判断交集长度符合阈值
					double cross_arc = MIN_END - MAX_BEG;
					if (round(cross_arc) < shortest_arc) continue;

					// 符合所有条件，构建单差模糊度
					t_gamb_Sd tmp(iter_satA, iter_satB, ID);
					tmp.sgl_obj = Ow_Ambs[iter_satA].rec;
					tmp.dbl_obj_1st = Ow_Ambs[iter_satA].sat;
					tmp.dbl_obj_2nd = Ow_Ambs[iter_satB].sat;
					tmp.beg = (MAX_BEG);
					tmp.end = (MIN_END);

					// 以测站索引的形式存储
					string sat_rec = tmp.sgl_obj + "_" + tmp.dbl_obj_1st + "_" + tmp.dbl_obj_2nd;
					_glog->logInfo("t_gallambs", "_cmb_intersat_Sd_amb", "combine Sd for " + sat_rec);
					this->_nSd_ambs[sys_key][tmp.sgl_obj].push_back(tmp);
				}
			}
		}
		return true;
	}

	bool t_gallambs::_cmb_interrec_Dd_amb(const pair<string, string>& rec_baseline, const t_gamb::AMB_ID& ID, double shortest_arc)
	{
		string func_id = "_cmb_interrec_Dd_amb";

		// 取出这个key的所有模糊度
		for (const auto& iter_key : _nSd_ambs)
		{
			// find the bug, fix
			// 每次只需要处理一个频率的单差就好
			const auto& key = iter_key.first;
			const auto& Ambs = iter_key.second;
			if (key.second != ID) continue;

			// 检查测站是否有模糊度
			if ((Ambs.count(rec_baseline.first) < 1) ||
				(Ambs.count(rec_baseline.second) < 1))
				return false;

			// 找到第一，第二个测站的所有单差模糊度
			const auto& Sd_1st = this->_nSd_ambs[key][rec_baseline.first];
			const auto& Sd_2nd = this->_nSd_ambs[key][rec_baseline.second];

			// 任意一个测站没有任何模糊度。直接退出
			if (Sd_1st.empty() || Sd_2nd.empty()) return false;

			// 开始构建双差模糊度
			// 原则上 
			// 1. 不同卫星系统不差分
			// 2. 共同弧度小于阈值不差分
			// 遍历第一个测站的单差模糊度
			const auto& Ow_1st = this->_nOw_ambs[key][rec_baseline.first];
			const auto& Ow_2nd = this->_nOw_ambs[key][rec_baseline.second];

			for (const auto& iter_1st : Sd_1st)
			{
				// 第一个测站对应的基线
				const auto& idx1 = (iter_1st.idx2ow_1st.find(1) != iter_1st.idx2ow_1st.end()) ? iter_1st.idx2ow_1st.find(1)->second : -1;
				const auto& idx2 = (iter_1st.idx2ow_1st.find(2) != iter_1st.idx2ow_1st.end()) ? iter_1st.idx2ow_1st.find(2)->second : -1;
				if (idx1 == -1 || idx2 == -1) continue;
				const auto& line_1st = make_pair(Ow_1st[idx1].sat, Ow_1st[idx2].sat);

				for (const auto& iter_2nd : Sd_2nd)
				{
					const auto& idx3 = (iter_2nd.idx2ow_1st.find(1) != iter_2nd.idx2ow_1st.end()) ? iter_2nd.idx2ow_1st.find(1)->second : -1;
					const auto& idx4 = (iter_2nd.idx2ow_1st.find(2) != iter_2nd.idx2ow_1st.end()) ? iter_2nd.idx2ow_1st.find(2)->second : -1;
					if (idx3 == -1 || idx4 == -1) continue;
					const auto& line_2nd = make_pair(Ow_2nd[idx3].sat, Ow_2nd[idx4].sat);

					// 先判断是否同一条测站基线
					if (line_1st != line_2nd) continue;

					// 再判断有没有交集
					const auto& amb_1st = iter_1st;
					const auto& amb_2nd = iter_2nd;

					auto MAX_BEG = MAX_TIME(amb_1st.beg, amb_2nd.beg);
					auto MIN_END = MIN_TIME(amb_1st.end, amb_2nd.end);

					// 没有任何交集
					if (MAX_BEG >= MIN_END) continue;

					// 有交集，再判断交集长度符合阈值
					double cross_arc = MIN_END - MAX_BEG;
					if (round(cross_arc) < shortest_arc) continue;

					// 顺序为sat1，sat2，rec1，rec2
					pair<string, string> idx_key(Ow_1st[idx1].rec, Ow_2nd[idx3].rec);

					// 符合所有条件，构建单差模糊度
					t_gamb_Dd tmp(ID, idx1, idx2, idx3, idx4);
					tmp.beg = MAX_BEG;
					tmp.end = MIN_END;

					tmp.sat_1st = Ow_1st[idx1].sat;
					tmp.sat_2nd = Ow_1st[idx2].sat;
					tmp.rec_1st = Ow_1st[idx1].rec;
					tmp.rec_2nd = Ow_2nd[idx3].rec;

					// 以测站基线索引的形式存储
					string sat_rec = tmp.sat_1st + "_" + tmp.sat_2nd + "_" + tmp.rec_1st + "_" + tmp.rec_2nd;
					_glog->logInfo("t_gallambs", "_cmb_intersat_Sd_amb", "combine Sd for " + sat_rec);
					this->_nDd_ambs[key][idx_key].push_back(tmp);

				}
			}
		}
		return true;

		//
	}

	bool t_gallambs::_cmb_intersat_Sd_amb(const pair<string, string>& sat_baseline, const t_gamb::AMB_ID& ID1, const t_gamb::AMB_ID& ID2, double shortest_arc)
	{
		string func_id = "_cmb_intersat_Sd_amb";

#ifdef DEBUG_AMBFIX
		// 输出调试信息
		if (this->_glog)
		{
			this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "combine Sd amb for : begin combine Sd ambs with one rec and two sat");
		}
#endif // DEBUG_AMBFIX

		// 检查卫星基线是一个系统
		if (sat_baseline.first.substr(0, 1) != sat_baseline.second.substr(0, 1))
		{
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_WARN, this->_class_id, func_id, "not the same system " +
				sat_baseline.first + " " + sat_baseline.second);
			return false;
		}

		// ID1 必须是 L1
		if (ID1 != t_gamb::AMB_ID::ALL_L1)
		{
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_WARN, this->_class_id, func_id, "ID1 should be L1 " +
				sat_baseline.first + " " + sat_baseline.second);
			return false;
		}

		// 检查有没有Ow数据
		// 判断系统有没有OW，没有OW
		// 取出ID1和ID2的所有OW模糊度
		t_dd_sys_index sys_key_1st = make_pair(t_gsys::sat2gsys(sat_baseline.first), ID1);
		t_dd_sys_index sys_key_2nd = make_pair(t_gsys::sat2gsys(sat_baseline.first), ID2);
		// 检查OW模糊度是否都存在
		if (this->_nOw_ambs.count(sys_key_1st) < 1 || 
			this->_nOw_ambs.count(sys_key_2nd) < 1 ||
			this->_nOw_ambs[sys_key_1st].empty() ||
			this->_nOw_ambs[sys_key_2nd].empty())
		{
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_WARN, this->_class_id, func_id, "no Ow " +
				sat_baseline.first + " " + sat_baseline.second);
			return false;
		}

		// 遍历测站
		const auto& sys_Ow_1st = this->_nOw_ambs[sys_key_1st];
		const auto& sys_Ow_2nd = this->_nOw_ambs[sys_key_2nd];
		for (const auto& iter_rec : this->_rec_list)
		{
			// 检查测站有没有Ow
			if (sys_Ow_1st.count(iter_rec) < 1 || sys_Ow_2nd.count(iter_rec) < 1)
			{
				if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_WARN, this->_class_id, func_id, "no Ow " + iter_rec);
				continue;
			}

			// 找到这个测站的所有模糊度
			const auto& Ow_Ambs_1st = this->_nOw_ambs[sys_key_1st][iter_rec];
			const auto& Ow_Ambs_2nd = this->_nOw_ambs[sys_key_2nd][iter_rec];

			if (Ow_Ambs_1st.empty()) continue;
			if (Ow_Ambs_2nd.empty()) continue;

			// 找到测站模糊度中，两个卫星的所有模糊度
			vector<int> Amb_satA_1st;
			vector<int> Amb_satB_1st;
			for (int amb_idx = 0; amb_idx < Ow_Ambs_1st.size(); amb_idx++)
			{
				if (Ow_Ambs_1st[amb_idx].sat == sat_baseline.first)  Amb_satA_1st.push_back(amb_idx);
				if (Ow_Ambs_1st[amb_idx].sat == sat_baseline.second) Amb_satB_1st.push_back(amb_idx);
			}
			vector<int> Amb_satA_2nd;
			vector<int> Amb_satB_2nd;
			for (int amb_idx = 0; amb_idx < Ow_Ambs_2nd.size(); amb_idx++)
			{
				if (Ow_Ambs_2nd[amb_idx].sat == sat_baseline.first)  Amb_satA_2nd.push_back(amb_idx);
				if (Ow_Ambs_2nd[amb_idx].sat == sat_baseline.second) Amb_satB_2nd.push_back(amb_idx);
			}

			// 任意一个没有模糊度，直接剔除
			if (Amb_satA_1st.empty() || Amb_satB_1st.empty() || Amb_satA_2nd.empty() || Amb_satB_2nd.empty())
			{
				if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_WARN, this->_class_id, func_id, "not ambs for " +
					sat_baseline.first + " or " + sat_baseline.second);
				continue;
			}

			// 遍历两颗卫星的模糊度，组成基线
			// 原则上 
			// 1. 不同卫星系统不差分
			// 2. 共同弧度小于阈值不差分
			// 3. 构建基线的时候，以L1为准
			for (const auto& iter_satA : Amb_satA_1st)
			{
				for (const auto& iter_satB : Amb_satB_1st)
				{
					// 再检查一次测站
					if (Ow_Ambs_1st[iter_satA].rec != Ow_Ambs_1st[iter_satB].rec) continue;

					// 再判断有没有交集
					const auto& amb_1st = Ow_Ambs_1st[iter_satA];
					const auto& amb_2nd = Ow_Ambs_1st[iter_satB];

					// 没有任何交集
					auto MAX_BEG = MAX_TIME(amb_1st.beg, amb_2nd.beg);
					auto MIN_END = MIN_TIME(amb_1st.end, amb_2nd.end);
					if (MAX_BEG >= MIN_END) continue;

					// 有交集，再判断交集长度符合阈值
					double cross_arc = MIN_END - MAX_BEG;
					if (round(cross_arc) < shortest_arc) continue;

					// L1的模糊度符合条件
					// 找一下有没有ID2的模糊度，ID1和ID2同时有才会成立
					map<int, int> idx2ow_1st; 
					idx2ow_1st[1] = iter_satA; 
					idx2ow_1st[2] = iter_satB;

					map<int, int> idx2ow_2nd;
					for (const auto& find_A : Amb_satA_2nd)
					{
						if (Ow_Ambs_2nd[find_A].beg != amb_1st.beg) continue;
						if (Ow_Ambs_2nd[find_A].end != amb_1st.end) continue;
						idx2ow_2nd[1] = find_A;
						break;
					}

					for (const auto& find_B : Amb_satB_2nd)
					{
						if (Ow_Ambs_2nd[find_B].beg != amb_2nd.beg) continue;
						if (Ow_Ambs_2nd[find_B].end != amb_2nd.end) continue;
						idx2ow_2nd[2] = find_B;
						break;
					}

					if (idx2ow_1st.size() < 2) continue;
					if (idx2ow_2nd.size() < 2) continue;

					// 符合所有条件，构建单差模糊度
					t_gamb::AMB_ID tmp_ID = t_gamb::AMB_ID::UNDEF;
					     if (ID2 == t_gamb::AMB_ID::ALL_L2) tmp_ID = t_gamb::AMB_ID::ALL_L1_L2;
					else if (ID2 == t_gamb::AMB_ID::ALL_L3) tmp_ID = t_gamb::AMB_ID::ALL_L1_L3;
					else continue;

					t_dd_sys_index sys_key = make_pair(t_gsys::sat2gsys(sat_baseline.first), tmp_ID);
					t_gamb_Sd tmp(idx2ow_1st, idx2ow_2nd, tmp_ID);
					tmp.sgl_obj     = Ow_Ambs_1st[iter_satA].rec;
					tmp.dbl_obj_1st = Ow_Ambs_1st[iter_satA].sat;
					tmp.dbl_obj_2nd = Ow_Ambs_1st[iter_satB].sat;
					tmp.beg = (MAX_BEG);
					tmp.end = (MIN_END);

					// 以测站索引的形式存储
					string sat_rec = tmp.sgl_obj + "_" + tmp.dbl_obj_1st + "_" + tmp.dbl_obj_2nd;
					_glog->logInfo("t_gallambs", "_cmb_intersat_Sd_amb", "combine Sd for " + sat_rec);
					this->_nSd_ambs[sys_key][tmp.sgl_obj].push_back(tmp);
				}
			}
		}
		return true;
	}

	bool t_gallambs::_cmb_interrec_Dd_amb(const pair<string, string>& rec_baseline, const t_gamb::AMB_ID& ID1, const t_gamb::AMB_ID& ID2, double shortest_arc)
	{
		string func_id = "_cmb_interrec_Dd_amb";
		t_gamb::AMB_ID ID = t_gamb::AMB_ID::UNDEF;
		if (ID1 != t_gamb::AMB_ID::ALL_L1) return false;
		if (ID1 == t_gamb::AMB_ID::ALL_L1 && ID2 == t_gamb::AMB_ID::ALL_L2)
		{
			ID = t_gamb::AMB_ID::ALL_L1_L2;
		}
		else if (ID1 == t_gamb::AMB_ID::ALL_L1 && ID2 == t_gamb::AMB_ID::ALL_L3)
		{
			ID = t_gamb::AMB_ID::ALL_L1_L3;
		}
		else
		{
			throw(logic_error("wrong ID of ID1 and ID2"));
			//return false;
		}

		// 取出这个key的所有模糊度
		for (const auto& iter_key : _nSd_ambs)
		{
			// find the bug, fix
			// 每次只需要处理一个频率的单差就好
			const auto& sys_amb_key  = iter_key.first;
			const auto&     rec_ambs = iter_key.second;
			if (sys_amb_key.second != ID) continue;

			// 检查测站是否有模糊度
			if ((rec_ambs.count(rec_baseline.first) < 1) ||
				(rec_ambs.count(rec_baseline.second) < 1))
				return false;

			// 找到第一，第二个测站的所有单差模糊度
			const auto& Sd_1st = this->_nSd_ambs[sys_amb_key][rec_baseline.first];
			const auto& Sd_2nd = this->_nSd_ambs[sys_amb_key][rec_baseline.second];

			// 任意一个测站没有任何模糊度。直接退出
			if (Sd_1st.empty() || Sd_2nd.empty()) return false;

			// 开始构建双差模糊度
			// 原则上 
			// 1. 不同卫星系统不差分
			// 2. 共同弧度小于阈值不差分
			// 遍历第一个测站的单差模糊度
			for (const auto& iter_1st : Sd_1st)
			{
				// 第一个测站对应的基线
				const auto& line_1st = make_pair(iter_1st.dbl_obj_1st, iter_1st.dbl_obj_2nd);
				for (const auto& iter_2nd : Sd_2nd)
				{
					// 第二个测站对应的基线
					const auto& line_2nd = make_pair(iter_2nd.dbl_obj_1st, iter_2nd.dbl_obj_2nd);
					// 先判断是否同一条测站基线
					if (line_1st != line_2nd) continue;
					// 再判断有没有交集
					const auto& amb_1st = iter_1st;
					const auto& amb_2nd = iter_2nd;

					auto MAX_BEG = MAX_TIME(amb_1st.beg, amb_2nd.beg);
					auto MIN_END = MIN_TIME(amb_1st.end, amb_2nd.end);

					// 没有任何交集
					if (MAX_BEG >= MIN_END) continue;

					// 有交集，再判断交集长度符合阈值
					double cross_arc = MIN_END - MAX_BEG;
					if (round(cross_arc) < shortest_arc) continue;

					// 顺序为sat1，sat2，rec1，rec2
					pair<string, string> idx_key(iter_1st.sgl_obj, iter_2nd.sgl_obj);

					// 符合所有条件，构建单差模糊度
					if (iter_1st.idx2ow_1st.size() < 2)continue;
					if (iter_1st.idx2ow_2nd.size() < 2)continue;
					if (iter_2nd.idx2ow_1st.size() < 2)continue;
					if (iter_2nd.idx2ow_2nd.size() < 2)continue;

					map<int, int> idx2ow_1st;
					map<int, int> idx2ow_2nd;
					idx2ow_1st[1] = iter_1st.idx2ow_1st.find(1)->second;
					idx2ow_1st[2] = iter_1st.idx2ow_1st.find(2)->second;
					idx2ow_1st[3] = iter_2nd.idx2ow_1st.find(1)->second;
					idx2ow_1st[4] = iter_2nd.idx2ow_1st.find(2)->second;

					idx2ow_2nd[1] = iter_1st.idx2ow_2nd.find(1)->second;
					idx2ow_2nd[2] = iter_1st.idx2ow_2nd.find(2)->second;
					idx2ow_2nd[3] = iter_2nd.idx2ow_2nd.find(1)->second;
					idx2ow_2nd[4] = iter_2nd.idx2ow_2nd.find(2)->second;

					t_gamb_Dd tmp(ID, idx2ow_1st, idx2ow_2nd);
					tmp.beg = MAX_BEG;
					tmp.end = MIN_END;

					tmp.sat_1st = iter_1st.dbl_obj_1st;
					tmp.sat_2nd = iter_1st.dbl_obj_2nd;
					tmp.rec_1st = iter_1st.sgl_obj;
					tmp.rec_2nd = iter_2nd.sgl_obj;

					// 以测站基线索引的形式存储
					string sat_rec = tmp.sat_1st + "_" + tmp.sat_2nd + "_" + tmp.rec_1st + "_" + tmp.rec_2nd;
					_glog->logInfo("t_gallambs", "_cmb_intersat_Sd_amb", "combine Sd for " + sat_rec);
					this->_nDd_ambs[sys_amb_key][idx_key].push_back(tmp);
				}
			}
		}
		return true;
	}

	void t_gallambs::_get_ows_table_from_pars()
	{
		if (this->_gpars == nullptr) return;
		for (unsigned int i = 0; i < _gpars->parNumber(); ++i)
		{
			const auto& par = (*_gpars)[i];
			//bds geo do not paticipate in ambfix
			if (t_gsys::bds_geo(par.prn)) continue;
			if (_rec_list.find(par.site) == _rec_list.end()) _rec_list.insert(par.site);
			if (_sat_list.find(par.prn) == _sat_list.end()) _sat_list.insert(par.prn);

			// 无电离层组合
			t_gamb::AMB_ID ID = t_gamb::AMB_ID::UNDEF;
			switch (par.parType)
			{
			case t_gpar::AMB_IF:  ID = t_gamb::AMB_ID::IF_L1_L2; break;
			case t_gpar::AMB13_IF:ID = t_gamb::AMB_ID::IF_L1_L3; break;
			case t_gpar::AMB_L1:  ID = t_gamb::AMB_ID::ALL_L1; break;
			case t_gpar::AMB_L2:  ID = t_gamb::AMB_ID::ALL_L2; break;
			case t_gpar::AMB_L3:  ID = t_gamb::AMB_ID::ALL_L3; break;
			case t_gpar::AMB_L4:  ID = t_gamb::AMB_ID::ALL_L4; break;
			case t_gpar::AMB_L5:  ID = t_gamb::AMB_ID::ALL_L5; break;
			default:  break;
			}
			if (_IDs.find(ID) == _IDs.end()) _IDs.insert(ID);

			t_gamb_Ow amb(ID, par.prn, par.site, par.beg, par.end);
			amb.value = par.value();
			auto gsys = t_gsys::sat2gsys(par.prn);
			auto key = make_pair(gsys, ID);
			this->_nOw_ambs[key][par.site].push_back(amb);
		}


	}

	void t_gallambs::_get_obs_table_from_gobs()
	{
		if (_IDs.empty()) return;
		//if (_IDs.find(t_gamb::AMB_ID::IF_L1_L2) == _IDs.end()) return;
		if (!this->_gbias) return;
		if (!this->_gobs) return;
		if (!this->_gset) return;

		// 初始化time_table,存储所有可能的时间节点
		double  proc_len = 0.0;
		t_gtime proc_beg = FIRST_TIME;
		t_gtime proc_end = FIRST_TIME;
		set<string>  rec_list;
		READ_OFILE_MODE rom;
		try
		{
			rec_list = dynamic_cast<t_gsetgen*>(_gset)->recs();
			proc_beg = dynamic_cast<t_gsetgen*>(_gset)->beg();
			proc_end = dynamic_cast<t_gsetgen*>(_gset)->end();
			proc_len = dynamic_cast<t_gsetamb*>(_gset)->wl_interval();
			rom = dynamic_cast<t_gsetproc*>(_gset)->read_ofile_mode();
		}
		catch (...)
		{
			return;
		}
		if (double_eq(proc_len, 0.0)) return;

		// 获取测站列表，测站列表为空就直接返回
		if (rec_list.empty()) return;

		// 获取时间列表
		time_table.clear();
		auto tmp_now = proc_beg - 2 * proc_len;
		auto tmp_end = (proc_end + 2 * proc_len);
		do
		{
			time_table.insert(tmp_now);
			tmp_now = tmp_now + proc_len;
		} while (tmp_now <= tmp_end);

		if (time_table.empty()) return;

		// 获取MW观测值
		for (const auto& sys : this->_nOw_ambs)
		{
			for (const auto& ambs : sys.second)
			{
				const string& rec = ambs.first;
				//test
				//cout << rec << endl;
				if (rec_list.find(rec) == rec_list.end()) continue;

				//read one obs file
				t_gallobs  obs;  obs.glog(_glog); obs.gset(_gset);
				if (rom == READ_OFILE_MODE::REALTIME)
				{
					t_gdata*   gdata = nullptr;
					vector<string> inp = dynamic_cast<t_gsetinp*>(_gset)->inputs(IFMT::RINEXO_INP);
					auto itINP = inp.begin();
					for (size_t i = 0; i < inp.size() && itINP != inp.end(); ++i, ++itINP)
					{
						// Get the file format/path, which will be used in decoder
						string path(*itINP);
						string id("ID" + int2str(i));

						// check the obs file valid
						auto pos = path.find_last_of("o");
						string name = path.substr(pos - 11, 4);
						transform(name.begin(), name.end(), name.begin(), ::toupper);
						if (name != rec) continue;

						// Check the file path
						if (path.substr(0, 7) != "file://")
						{
							_glog->comment(0, "main", "path is not file (skipped)!");
							continue;
						}

						t_gio* gio = nullptr;
						t_gcoder* gcoder = nullptr;

						gdata = &obs; gcoder = new t_rinexo(_gset, "", 5000);

						// runepoch for the time costed each epoch (i guess)
						t_gtime runepoch(t_gtime::GPS);
						t_gtime lstepoch(t_gtime::GPS);
						// READ DATA FROM FILE
						if (gcoder) {
							// prepare gio for the file
							gio = new t_gfile;
							gio->glog(_glog);
							gio->path(path);

							// Put the file into gcoder
							gcoder->clear();
							gcoder->path(path);
							gcoder->glog(_glog);
							// Put the data container into gcoder
							gcoder->add_data(id, gdata);
							gcoder->add_data("OBJ", _gobj);
							gio->coder(gcoder);

							// Read the data from file here
							gio->run_read();
							// Write the information of reading process to log file
							_glog->comment(t_glog::LOG_LV::LOG_INFO, "main", "main", "READ: " + path + " time: "
								+ dbl2str(lstepoch.diff(runepoch)) + " sec");
							delete gio; gio = nullptr;
							delete gcoder; gcoder = nullptr;
						}
					}
					_gobs = &obs;
				}
			
				// 遍历找出测站观测值区间
				t_gtime min = LAST_TIME;
				t_gtime max = FIRST_TIME;
				for (const auto& amb : ambs.second)
				{
					if (min > amb.beg) min = amb.beg;
					if (max < amb.end) max = amb.end;
				}

				// 判断是否超限
				auto it_min = time_table.lower_bound(min);
				auto it_max = time_table.lower_bound(max);
				if (it_min == time_table.end()) it_min = time_table.begin();
				if (it_max == time_table.end()) it_max = time_table.end();

				// 循环获取观测值
				// 获取所有需要用到这个时间观测值的卫星, 如果没有卫星，直接跳过
				map<t_gtime, set<string>> sat_used;
				for (const auto& amb : ambs.second)
				{
					if (amb.beg > amb.end) continue;
					auto beg = time_table.lower_bound(amb.beg);
					auto end = time_table.lower_bound(amb.end);
					if (beg == time_table.end()) continue;
					if (end == time_table.end()) end = time_table.end();

					/*auto now = beg;
					for (now; *now <= *end && now != time_table.end(); now++)
					{
						sat_used[*now].insert(amb.sat);
					}*/
					auto now = beg;
					for (now; now != end; now++)
					{
						sat_used[*now].insert(amb.sat);
					}
				}

				// 开始循环
				auto loop_now = it_min;
				auto loop_end = (it_max == time_table.end()) ? it_max : (++it_max);
				
				for (loop_now; loop_now != loop_end; loop_now++)
				{
					// 获取这个测站在这个时间的所有观测值，没有观测值直接跳过
					vector<t_gobsgnss*> sat_datas = _gobs->obs(true, rec, *loop_now);
					if (sat_datas.empty())
					{
						continue;
					}

					if (sat_used.empty())
					{
						continue;
					}

					// 遍历之前获取的观测值，把有用的观测值加入存储变量中
					for (auto& sat_iter : sat_datas)
					{
						auto sat = sat_iter->sat();
						auto sys = t_gsys::sat2gsys(sat);

						if (sat_used[*loop_now].count(sat) <= 0) continue;

						pair<string, t_gtime> key(sat, *loop_now);
						if (this->_MW_used[sys][rec].count(key) < 1)
						{
							// 没有观测值，也组不出MW 
							if (sat_iter->obs_empty()) continue;

							GOBSBAND b1 = _band_index[sys][FREQ_1];
							GOBSBAND b2 = _band_index[sys][FREQ_2];
							GOBSBAND b3 = _band_index[sys][FREQ_3];

							t_gobs gobsP1;
							t_gobs gobsP2;
							t_gobs gobsP3;
							//cout << rec << "  " << sat << "  "<< setw(15)<<setprecision(8)<< loop_now->dmjd();
							t_gobs gobsL1;
							t_gobs gobsL2;
							t_gobs gobsL3;

							bool isRawAll = dynamic_cast<t_gsetproc*>(_gset)->obs_combin() == RAW_ALL;
							bool isDouble = dynamic_cast<t_gsetproc*>(_gset)->frequency()  == 2;
							// jdhuang:非组合的时候特殊处理
							if (!isDouble)
							{
								gobsP1 = t_gobs(sat_iter->select_range(b1, isRawAll));
								gobsP2 = t_gobs(sat_iter->select_range(b2, isRawAll));
								gobsP3 = t_gobs(sat_iter->select_range(b3, isRawAll));
								gobsL1 = t_gobs(sat_iter->select_phase(b1, isRawAll));
								gobsL2 = t_gobs(sat_iter->select_phase(b2, isRawAll));
								gobsL3 = t_gobs(sat_iter->select_phase(b3, isRawAll));
							}
							else
							{
								gobsP1 = t_gobs(sat_iter->select_range(b1, isRawAll));
								gobsP2 = t_gobs(sat_iter->select_range(b2, isRawAll));
								gobsL1 = t_gobs(sat_iter->select_phase(b1, isRawAll));
								gobsL2 = t_gobs(sat_iter->select_phase(b2, isRawAll));
							}
							
							double mw_obs_12 = 0.0; 
							double mw_obs_13 = 0.0; 
							//bool isIono13 = (isDouble && !isRawAll && b1 == t_gsys::band_priority(sys, FREQ_1) && b2 == t_gsys::band_priority(sys, FREQ_3));
							sat_iter->apply_bias(*_gbias);
							if(!isDouble)
							{
								mw_obs_12 = sat_iter->MW_cycle(gobsL1, gobsL2, gobsP1, gobsP2);
								mw_obs_13 = sat_iter->MW_cycle(gobsL1, gobsL3, gobsP1, gobsP3);
								if (double_eq(mw_obs_12, 0.0)) continue;
								if (double_eq(mw_obs_13, 0.0)) continue;

								//auto obs_type = OBSCOMBIN::MW_COMBIN;
								//mw_obs_12 = mw_obs_12 + _gprecise->cmbIntraFreqDcb("CAS", obs_type, this->_gbias, sat_iter, &gobsP1, &gobsP2);
								//mw_obs_13 = mw_obs_13 + _gprecise->cmbIntraFreqDcb("CAS", obs_type, this->_gbias, sat_iter, &gobsP1, &gobsP3);
							}
							else
							{
								mw_obs_12 = sat_iter->MW_cycle(gobsL1, gobsL2, gobsP1, gobsP2);
								if (double_eq(mw_obs_12, 0.0)) continue;

								//auto obs_type = OBSCOMBIN::MW_COMBIN;
								//mw_obs_12 = mw_obs_12 + _gprecise->cmbIntraFreqDcb("CAS", obs_type, this->_gbias, sat_iter, &gobsP1, &gobsP2);
							}

							this->_MW_used[sys][rec][key][2] = mw_obs_12;
							this->_MW_used[sys][rec][key][3] = mw_obs_13;
						}
					}// finsh time
				} // finsh time loop

			}// finsh the rec	
		}// finsh sys

	}
	void t_gallambs::_get_dependent_IF_Sd(const t_gamb::AMB_ID& ID)
	{
		string func_id = "_get_dependent_Sd";
		for (auto& sys : this->_nSd_ambs)
		{
			if (sys.first.second != ID) continue;
			if (sys.second.empty()) continue;
			for (auto& bl : sys.second)
			{
				Matrix vector_space;
				auto num_Ow_1st = this->_nOw_ambs[sys.first][bl.first].size();
				int ncol = num_Ow_1st;
				vector<t_gamb_Sd> result;
				for (auto& amb : bl.second)
				{
					int index1 = (amb.idx2ow_1st.find(1) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(1)->second + 1) : -1;
					int index2 = (amb.idx2ow_1st.find(2) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(2)->second + 1) : -1;
					if (index1 < 0 || index2 < 0 ) continue;
					if (vector_space.is_zero())
					{
						vector_space.resize(1, ncol);
						vector_space = 0.0;
						vector_space(1, index1) = 1;
						vector_space(1, index2) = -1;
						result.push_back(amb);
						if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "Base line is : " +
							int2str(result.size() + 1) + " " +
							amb.sgl_obj + " " +
							amb.dbl_obj_1st + " " +
							amb.dbl_obj_2nd);
						continue;
					}
					RowVector new_vector(0.0, ncol);
					new_vector(index1) = 1;
					new_vector(index2) = -1;
					Matrix A = vector_space & new_vector;
					DiagonalMatrix D;
					if (A.nrows() < A.ncols())
					{
						SVD(A.t(), D);
					}
					else
					{
						SVD(A, D);
					}
					bool isIndepend = (fabs(D(D.nrows(), D.ncols())) > 0.01);
					if (isIndepend)
					{
						vector_space = A;
						result.push_back(amb);
						if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "  Independ Base line is : " +
							int2str(result.size() + 1) + " " +
							amb.sgl_obj + " " +
							amb.dbl_obj_1st + " " +
							amb.dbl_obj_2nd);
					}
					else
					{
						if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "  Undepend Base line is : " +
							int2str(result.size() + 1) + " " +
							amb.sgl_obj + " " +
							amb.dbl_obj_1st + " " +
							amb.dbl_obj_2nd);
					}
				}
				bl.second = result;
			}
		}
	}

	bool t_gallambs::_checkAmbDepend(bool isFirst, int iNamb, int *iNdef, int iN_oneway, int *arriIpt2ow, int iMaxamb_ow, int iMaxamb_for_check)
	{

		const string strCprogName = "check_amb_depend";
		const double dEPS = 1e-12;
		static int iNdim_ow = 0;
		static int iNdim_for_check = 0;
		int i, j;
		//static double *pdE = 0;
		//static double *pdC = 0;
		const double dOper[4] = { 1.0, -1.0, -1.0, 1.0 };
		double dC_dot;
		bool isLdepend;

		try
		{
			if (isFirst)
			{
				if (_pdE) { delete _pdE; _pdE = nullptr; }
				if (_pdC) { delete _pdC; _pdC = nullptr; }
				//lvhb test in 202010
				_pdE = new double[iMaxamb_ow * (iMaxamb_for_check + 100)];
				_pdC = new double[iMaxamb_for_check + 100];
				if (iMaxamb_ow == 0 || iMaxamb_for_check == 0)
				{
					string err = "***ERROR: memory allocatation for e&c ";
					throw err;
				}
				iNdim_ow = iMaxamb_ow;
				iNdim_for_check = iMaxamb_for_check;
				return false;
			}

			if (*iNdef >= iNdim_for_check)
			{
				string err2 = "***ERROR: independent ones already reaches the allocated memory ";
				throw err2;
			}

			dC_dot = 0.0;
			for (i = 0; i < *iNdef; i++)
			{
				_pdC[i] = 0.0;
				for (j = 0; j < iN_oneway; j++)
				{
					if (arriIpt2ow[j] > iNdim_ow)
					{
						string err3 = "***ERROR: base element beyond the allocated memory ";
						throw err3;
					}
					_pdC[i] = _pdC[i] + _pdE[(arriIpt2ow[j] - 1)*(iNdim_for_check + 1) + i] * dOper[j];
				}
				dC_dot = dC_dot + _pdC[i] * _pdC[i];
			}

			isLdepend = false;
			if (abs(dC_dot - iN_oneway) <= dEPS)
			{
				isLdepend = true;
				return isLdepend;
			}

			for (j = 0; j < iNamb; j++)
			{
				_pdE[j * (iNdim_for_check + 1) + *iNdef] = 0.0;
			}
			for (j = 0; j < iN_oneway; j++)
			{
				_pdE[(arriIpt2ow[j] - 1) * (iNdim_for_check + 1) + *iNdef] = dOper[j];
			}

			if (dC_dot > dEPS)
			{
				for (i = 0; i < *iNdef; i++)
				{
					if (abs(_pdC[i]) < dEPS)
					{
						continue;
					}
					for (j = 0; j < iNamb; j++)
					{
						_pdE[j * (iNdim_for_check + 1) + *iNdef] = _pdE[j * (iNdim_for_check + 1) + *iNdef] - _pdC[i] * _pdE[j * (iNdim_for_check + 1) + i];
					}
				}

				dC_dot = 0.0;
				for (j = 0; j < iNamb; j++)
				{
					dC_dot = dC_dot + _pdE[j * (iNdim_for_check + 1) + *iNdef] * _pdE[j * (iNdim_for_check + 1) + *iNdef];
				}
				dC_dot = sqrt(dC_dot);
				for (j = 0; j < iNamb; j++)
				{
					_pdE[j * (iNdim_for_check + 1) + *iNdef] = _pdE[j * (iNdim_for_check + 1) + *iNdef] / dC_dot;
				}
			}
			else
			{
				dC_dot = sqrt(iN_oneway * 1.0);
				for (j = 0; j < iN_oneway; j++)
				{
					_pdE[(arriIpt2ow[j] - 1) * (iNdim_for_check + 1) + *iNdef] = _pdE[(arriIpt2ow[j] - 1) * (iNdim_for_check + 1) + *iNdef] / dC_dot;
				}
			}

			_pdC[*iNdef] = dC_dot;
			*iNdef = *iNdef + 1;
			return false;
		}
		catch (...)
		{
			write_log_info(_glog, 1, "ERROR[t_gambiguity::_checkAmbDepend] : throw exception");
			return false;
		}
	}

	void t_gallambs::_get_dependent_IF_Sd2(const t_gamb::AMB_ID& ID)
	{
		const double dEPS = 1e-12;
		static int iNdim_ow = 0;
		static int iNdim_for_check = 0;
		int i, j;
		//static double *pdE = 0;
		//static double *pdC = 0;
		const double dOper[4] = { 1.0, -1.0, -1.0, 1.0 };
		double dC_dot;
		bool isLdepend;
		int ndef=0;
		int iMaxamb_ow;
		int iMaxamb_for_check;
		bool isfirst = true;
		bool is_depend;
		string func_id = "_get_dependent_Sd";
		for (auto& sys : this->_nSd_ambs)
		{
			if (sys.first.second != ID) continue;
			if (sys.second.empty()) continue;
			for (auto& bl : sys.second)
			{
				auto& bl_ambs_fix_WL = this->_nSd_ambs_fix_WL[sys.first][bl.first];
				auto& bl_ambs_fix_NO = this->_nSd_ambs_fix_NO[sys.first][bl.first];
				Matrix vector_space;
				auto num_Ow_1st = this->_nOw_ambs[sys.first][bl.first].size();
				vector<t_gamb_Sd> result;
				for (auto& amb : bl.second)
				{
					int ipt2ow[2];
					ipt2ow[0] = amb.idx2ow_1st.find(1)->second + 1;
					ipt2ow[1] = amb.idx2ow_1st.find(2)->second + 1;
					
					is_depend = _checkAmbDepend(isfirst, num_Ow_1st, &ndef, 2, ipt2ow, num_Ow_1st, 4000);
					isfirst = false;

					if (!is_depend)
					{
						result.push_back(amb);
						if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "  Independ Base line is : " +
							int2str(result.size() + 1) + " " +
							amb.sgl_obj + " " +
							amb.dbl_obj_1st + " " +
							amb.dbl_obj_2nd);
					}
					else
					{
						if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "  Undepend Base line is : " +
							int2str(result.size() + 1) + " " +
							amb.sgl_obj + " " +
							amb.dbl_obj_1st + " " +
							amb.dbl_obj_2nd);
					}
				}
				bl.second = result;

				result.clear();
				for (auto& amb : bl_ambs_fix_WL)
				{
					int ipt2ow[2];
					ipt2ow[0] = amb.idx2ow_1st.find(1)->second + 1;
					ipt2ow[1] = amb.idx2ow_1st.find(2)->second + 1;

					is_depend = _checkAmbDepend(isfirst, num_Ow_1st, &ndef, 2, ipt2ow, num_Ow_1st, 4000);
					isfirst = false;

					if (!is_depend)
					{
						result.push_back(amb);
						if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "  Independ Base line is : " +
							int2str(result.size() + 1) + " " +
							amb.sgl_obj + " " +
							amb.dbl_obj_1st + " " +
							amb.dbl_obj_2nd);
					}
					else
					{
						if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "  Undepend Base line is : " +
							int2str(result.size() + 1) + " " +
							amb.sgl_obj + " " +
							amb.dbl_obj_1st + " " +
							amb.dbl_obj_2nd);
					}
				}

				bl_ambs_fix_WL = result;
				result.clear();

				for (auto& amb : bl_ambs_fix_NO)
				{
					int ipt2ow[2];
					ipt2ow[0] = amb.idx2ow_1st.find(1)->second + 1;
					ipt2ow[1] = amb.idx2ow_1st.find(2)->second + 1;

					is_depend = _checkAmbDepend(isfirst, num_Ow_1st, &ndef, 2, ipt2ow, num_Ow_1st, 4000);
					isfirst = false;

					if (!is_depend)
					{
						result.push_back(amb);
						if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "  Independ Base line is : " +
							int2str(result.size() + 1) + " " +
							amb.sgl_obj + " " +
							amb.dbl_obj_1st + " " +
							amb.dbl_obj_2nd);
					}
					else
					{
						if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "  Undepend Base line is : " +
							int2str(result.size() + 1) + " " +
							amb.sgl_obj + " " +
							amb.dbl_obj_1st + " " +
							amb.dbl_obj_2nd);
					}
					
				}
				bl_ambs_fix_NO = result;
			}
		}		
	}

	void t_gallambs::_get_dependent_IF_Dd(const t_gamb::AMB_ID& ID)
	{
		string func_id = "_get_dependent_Dd";
		// get the sys
		for (auto& sys : this->_nDd_ambs)
		{
			if (sys.first.second != ID) continue;
			if (sys.second.empty()) continue;
			// get the bl
			for (auto& bl : sys.second)
			{
				vector<t_gamb_Dd> result;
				vector<t_gamb_Dd> result_fix_WL;
				vector<t_gamb_Dd> result_fix_NO;
				Matrix vector_space;
				// 确定了双差模糊度的列表，判断独立性
				// check the depend
				auto num_Ow_1st = this->_nOw_ambs[sys.first][bl.first.first].size();
				auto num_Ow_2nd = this->_nOw_ambs[sys.first][bl.first.second].size();

				//auto num_Ow_1st_WL = this->_nDd_ambs_fix_WL[sys.first][bl.first].size();
				//auto num_Ow_2nd_WL = this->_nDd_ambs_fix_WL[sys.first][bl.first].size();

				//auto num_Ow_1st_NL = this->_nDd_ambs_fix_NO[sys.first][bl.first].size();
				//auto num_Ow_2nd_NL = this->_nDd_ambs_fix_NO[sys.first][bl.first].size();

				//int ncol = num_Ow_1st + num_Ow_2nd;

				auto& bl_ambs_fix_WL = this->_nDd_ambs_fix_WL[sys.first][bl.first];
				auto& bl_ambs_fix_NO = this->_nDd_ambs_fix_NO[sys.first][bl.first];

				math_dependent_IF_Dd(bl.second,      num_Ow_1st, num_Ow_2nd, result,        vector_space);
				math_dependent_IF_Dd(bl_ambs_fix_WL, num_Ow_1st, num_Ow_2nd, result_fix_WL, vector_space);
				math_dependent_IF_Dd(bl_ambs_fix_NO, num_Ow_1st, num_Ow_2nd, result_fix_NO, vector_space);
				//for (auto& amb : bl.second)
				//{
				//	int index1 = (amb.idx2ow_1st.find(1) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(1)->second + 1) : -1;
				//	int index2 = (amb.idx2ow_1st.find(2) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(2)->second + 1) : -1;
				//	int index3 = (amb.idx2ow_1st.find(3) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(3)->second + 1) : -1;
				//	int index4 = (amb.idx2ow_1st.find(4) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(4)->second + 1) : -1;

				//	if (index1 < 0 || index2 < 0 || index3 < 0 || index4 < 0) continue;

				//	// 构建向量空间
				//	if (vector_space.is_zero())
				//	{
				//		vector_space.resize(1, ncol);
				//		vector_space = 0.0;
				//		vector_space(1, index1) = 1;
				//		vector_space(1, index2) = -1;
				//		vector_space(1, index3 + num_Ow_1st) = -1;
				//		vector_space(1, index4 + num_Ow_1st) = 1;
				//		result.push_back(amb);

				//		if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "Base line is : " +
				//			int2str(result.size() + 1) + " " +
				//			amb.rec_1st + " " +
				//			amb.rec_2nd + " " +
				//			amb.sat_1st + " " +
				//			amb.sat_2nd);
				//		continue;
				//	}

				//	RowVector new_vector(0.0, ncol);
				//	new_vector(index1) = 1;
				//	new_vector(index2) = -1;
				//	new_vector(index3 + num_Ow_1st) = -1;
				//	new_vector(index4 + num_Ow_1st) = 1;

				//	// 数学原理：奇异矩阵经过SVD分解后，D矩阵对角线的非零元素个数等于A矩阵的矩阵秩
				//	Matrix A = vector_space & new_vector;
				//	DiagonalMatrix D;
				//	if (A.nrows() < A.ncols())
				//	{
				//		SVD(A.t(), D);
				//	}
				//	else
				//	{
				//		SVD(A, D);
				//	}
				//	bool isIndepend = (fabs(D(D.nrows(), D.ncols())) > 0.01);

				//	if (isIndepend)
				//	{
				//		vector_space = A;
				//		result.push_back(amb);
				//		if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "  Independ Base line is : " +
				//			int2str(result.size() + 1) + " " +
				//			amb.rec_1st + " " +
				//			amb.rec_2nd + " " +
				//			amb.sat_1st + " " +
				//			amb.sat_2nd);
				//	}
				//}

				bl.second = result;
				bl_ambs_fix_WL = result_fix_WL;
				bl_ambs_fix_NO = result_fix_NO;
			}
		}
	}

	void t_gallambs::_get_dependent_1X_Dd(const t_gamb::AMB_ID& ID)
	{
		string func_id = "_get_dependent_1X_Dd";

		// get the sys
		for (auto& sys_id_key : this->_nDd_ambs)
		{
			const auto& gsys = sys_id_key.first.first;
			const auto& gids = sys_id_key.first.second;
			auto& gambs = sys_id_key.second;
			if (gids != ID)    continue;
			if (gambs.empty()) continue;

			// get the bl
			for (auto& amb_for_one_rec_bl : gambs)
			{
				// 找到OW
				pair<GSYS, t_gamb::AMB_ID> index_1;
				pair<GSYS, t_gamb::AMB_ID> index_2;
				if (ID == t_gamb_Dd::AMB_ID::ALL_L1_L2)
				{
					index_1 = make_pair(gsys, t_gamb_Dd::AMB_ID::ALL_L1);
					index_2 = make_pair(gsys, t_gamb_Dd::AMB_ID::ALL_L2);
				}
				else if (ID == t_gamb_Dd::AMB_ID::ALL_L1_L3)
				{
					index_1 = make_pair(gsys, t_gamb_Dd::AMB_ID::ALL_L1);
					index_2 = make_pair(gsys, t_gamb_Dd::AMB_ID::ALL_L3);
				}
				else
				{
					throw(logic_error("wrong amb id"));
				}

				const auto& baseline = amb_for_one_rec_bl.first;
				auto& bl_ambs  = amb_for_one_rec_bl.second;
				auto& bl_ambs_fix_WL = this->_nDd_ambs_fix_WL[sys_id_key.first][amb_for_one_rec_bl.first];
				auto& bl_ambs_fix_NO = this->_nDd_ambs_fix_NO[sys_id_key.first][amb_for_one_rec_bl.first];
				const auto& rec_1st = baseline.first;
				const auto& rec_2nd = baseline.second;

				// 确定了双差模糊度的列表，判断独立性
				// check the depend
				// for L1
				auto num_Ow_1st = this->_nOw_ambs[index_1][rec_1st].size();
				auto num_Ow_2nd = this->_nOw_ambs[index_1][rec_2nd].size();

				// for LX
				auto num_Ow_3rd = this->_nOw_ambs[index_2][rec_1st].size();
				auto num_Ow_4th = this->_nOw_ambs[index_2][rec_2nd].size();

				int ncol     = num_Ow_1st + num_Ow_2nd + num_Ow_3rd + num_Ow_4th;
				int ncol_1st = num_Ow_1st + num_Ow_2nd;

				vector<t_gamb_Dd> result;	
				vector<t_gamb_Dd> result_fix_WL;	
				vector<t_gamb_Dd> result_fix_NO;	
				Matrix vector_space;

				math_dependent_1X_Dd(bl_ambs, num_Ow_1st, num_Ow_2nd, num_Ow_3rd, num_Ow_4th, result, vector_space);
				math_dependent_1X_Dd(bl_ambs_fix_WL, num_Ow_1st, num_Ow_2nd, num_Ow_3rd, num_Ow_4th, result, vector_space);
				math_dependent_1X_Dd(bl_ambs_fix_NO, num_Ow_1st, num_Ow_2nd, num_Ow_3rd, num_Ow_4th, result, vector_space);
				//for (const auto& amb : bl_ambs)
				//{
				//	int index_1st_1 = (amb.idx2ow_1st.find(1) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(1)->second + 1) : -1;
				//	int index_1st_2 = (amb.idx2ow_1st.find(2) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(2)->second + 1) : -1;
				//	int index_1st_3 = (amb.idx2ow_1st.find(3) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(3)->second + 1) : -1;
				//	int index_1st_4 = (amb.idx2ow_1st.find(4) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(4)->second + 1) : -1;

				//	int index_2nd_1 = (amb.idx2ow_2nd.find(1) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(1)->second + 1) : -1;
				//	int index_2nd_2 = (amb.idx2ow_2nd.find(2) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(2)->second + 1) : -1;
				//	int index_2nd_3 = (amb.idx2ow_2nd.find(3) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(3)->second + 1) : -1;
				//	int index_2nd_4 = (amb.idx2ow_2nd.find(4) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(4)->second + 1) : -1;

				//	if (index_1st_1 < 0 || index_1st_2 < 0 || index_1st_3 < 0 || index_1st_4 < 0) continue;
				//	if (index_2nd_1 < 0 || index_2nd_2 < 0 || index_2nd_3 < 0 || index_2nd_4 < 0) continue;

				//	// 构建向量空间
				//	if (vector_space.is_zero())
				//	{
				//		vector_space.resize(1, ncol);
				//		vector_space = 0.0;
				//		vector_space(1, index_1st_1) = 1;
				//		vector_space(1, index_1st_2) = -1;
				//		vector_space(1, index_1st_3 + num_Ow_1st) = -1;
				//		vector_space(1, index_1st_4 + num_Ow_1st) = 1;

				//		vector_space(1, ncol_1st + index_2nd_1) = 1;
				//		vector_space(1, ncol_1st + index_2nd_2) = -1;
				//		vector_space(1, ncol_1st + index_2nd_3 + num_Ow_3rd) = -1;
				//		vector_space(1, ncol_1st + index_2nd_4 + num_Ow_3rd) = 1;

				//		result.push_back(amb);

				//		if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "Base line is : " +
				//			int2str(result.size() + 1) + " " +
				//			amb.rec_1st + " " +
				//			amb.rec_2nd + " " +
				//			amb.sat_1st + " " +
				//			amb.sat_2nd);
				//		continue;
				//	}

				//	RowVector new_vector(0.0, ncol);
				//	new_vector(           index_1st_1             ) =  1;
				//	new_vector(           index_1st_2             ) = -1;
				//	new_vector(           index_1st_3 + num_Ow_1st) = -1;
				//	new_vector(           index_1st_4 + num_Ow_1st) =  1;
				//	new_vector(ncol_1st + index_2nd_1             ) =  1;
				//	new_vector(ncol_1st + index_2nd_2             ) = -1;
				//	new_vector(ncol_1st + index_2nd_3 + num_Ow_3rd) = -1;
				//	new_vector(ncol_1st + index_2nd_4 + num_Ow_3rd) =  1;

				//	// 数学原理：奇异矩阵经过SVD分解后，D矩阵对角线的非零元素个数等于A矩阵的矩阵秩
				//	Matrix A = vector_space & new_vector;
				//	DiagonalMatrix D;
				//	if (A.nrows() < A.ncols())
				//	{
				//		SVD(A.t(), D);
				//	}
				//	else
				//	{
				//		SVD(A, D);
				//	}
				//	bool isIndepend = (fabs(D(D.nrows(), D.ncols())) > 0.01);

				//	if (isIndepend)
				//	{
				//		vector_space = A;
				//		result.push_back(amb);
				//		if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "  Independ Base line is : " +
				//			int2str(result.size() + 1) + " " +
				//			amb.rec_1st + " " +
				//			amb.rec_2nd + " " +
				//			amb.sat_1st + " " +
				//			amb.sat_2nd);
				//	}
				//}

				amb_for_one_rec_bl.second = result;
				bl_ambs_fix_WL = result_fix_WL;
				bl_ambs_fix_NO = result_fix_NO;
			}
		}

	}

	void t_gallambs::_get_fix_rate(const t_gamb::AMB_ID& ID)
	{

		// fix ALL
		for (auto& sys : this->_nDd_ambs)
		{
			const auto& amb_id = sys.first.second;
			const auto& sys_key = sys.first;
			if (amb_id != ID) continue;
			for (auto& bl : sys.second)
			{
				const auto& rec_1st = bl.first.first;
				const auto& rec_2nd = bl.first.second;

				if (bl.second.empty()) continue;
				for (auto& amb : bl.second)
				{
					all_NL_num[sys.first] = all_NL_num[sys.first] + 1;
					all_WL_num[sys.first] = all_WL_num[sys.first] + 1;

					if (amb.MW_fixed) fix_WL_num[sys.first] = fix_WL_num[sys.first] + 1;
					if (amb.NL_fixed) fix_NL_num[sys.first] = fix_NL_num[sys.first] + 1;
				}
			}
		}

		// fix WL
		for (auto& sys : this->_nDd_ambs_fix_WL)
		{
			const auto& amb_id = sys.first.second;
			const auto& sys_key = sys.first;
			if (amb_id != ID) continue;
			for (auto& bl : sys.second)
			{
				const auto& rec_1st = bl.first.first;
				const auto& rec_2nd = bl.first.second;

				if (bl.second.empty()) continue;
				for (auto& amb : bl.second)
				{
					if (!double_eq(amb.NL_value, 0.0) && !double_eq(amb.MW_value, 0.0))
					{
						all_NL_num[sys.first] = all_NL_num[sys.first] + 1;
						all_WL_num[sys.first] = all_WL_num[sys.first] + 1;
					}

					if (amb.MW_fixed) fix_WL_num[sys.first] = fix_WL_num[sys.first] + 1;
					if (amb.NL_fixed) fix_NL_num[sys.first] = fix_NL_num[sys.first] + 1;
				}
			}
		}

		// fix NO
		for (auto& sys : this->_nDd_ambs_fix_NO)
		{
			const auto& amb_id = sys.first.second;
			const auto& sys_key = sys.first;
			if (amb_id != ID) continue;
			for (auto& bl : sys.second)
			{
				const auto& rec_1st = bl.first.first;
				const auto& rec_2nd = bl.first.second;

				if (bl.second.empty()) continue;
				for (auto& amb : bl.second)
				{
					if (!double_eq(amb.NL_value, 0.0) && !double_eq(amb.MW_value, 0.0))
					{
						all_NL_num[sys.first] = all_NL_num[sys.first] + 1;
						all_WL_num[sys.first] = all_WL_num[sys.first] + 1;
					}

					if (amb.MW_fixed) fix_WL_num[sys.first] = fix_WL_num[sys.first] + 1;
					if (amb.NL_fixed) fix_NL_num[sys.first] = fix_NL_num[sys.first] + 1;
				}
			}
		}
	}

	void t_gallambs::_get_SDfix_rate(const t_gamb::AMB_ID& ID)
	{

		// fix ALL
		for (auto& sys : this->_nSd_ambs)
		{
			if (sys.first.second != ID) continue;
			if (sys.second.empty()) continue;
			for (auto& bl : sys.second)
			{
				if (bl.second.empty()) continue;
				for (auto& amb : bl.second)
				{
					const auto& rec = amb.sgl_obj;
					all_NL_num[sys.first] = all_NL_num[sys.first] + 1;
					all_WL_num[sys.first] = all_WL_num[sys.first] + 1;

					if (amb.MW_fixed) fix_WL_num[sys.first] = fix_WL_num[sys.first] + 1;
					if (amb.NL_fixed) fix_NL_num[sys.first] = fix_NL_num[sys.first] + 1;
				}
			}
		}

		// fix WL
		for (auto& sys : this->_nSd_ambs_fix_WL)
		{
			if (sys.first.second != ID) continue;
			if (sys.second.empty()) continue;
			for (auto& bl : sys.second)
			{
				if (bl.second.empty()) continue;
				for (auto& amb : bl.second)
				{
					if (!double_eq(amb.NL_value, 0.0) && !double_eq(amb.MW_value, 0.0))
					{
						all_NL_num[sys.first] = all_NL_num[sys.first] + 1;
						all_WL_num[sys.first] = all_WL_num[sys.first] + 1;
					}

					if (amb.MW_fixed) fix_WL_num[sys.first] = fix_WL_num[sys.first] + 1;
					if (amb.NL_fixed) fix_NL_num[sys.first] = fix_NL_num[sys.first] + 1;
				}
			}
		}

		// fix NO
		for (auto& sys : this->_nSd_ambs_fix_NO)
		{
			if (sys.first.second != ID) continue;
			if (sys.second.empty()) continue;
			for (auto& bl : sys.second)
			{
				if (bl.second.empty()) continue;
				for (auto& amb : bl.second)
				{
					if (!double_eq(amb.NL_value, 0.0) && !double_eq(amb.MW_value, 0.0))
					{
						all_NL_num[sys.first] = all_NL_num[sys.first] + 1;
						all_WL_num[sys.first] = all_WL_num[sys.first] + 1;
					}

					if (amb.MW_fixed) fix_WL_num[sys.first] = fix_WL_num[sys.first] + 1;
					if (amb.NL_fixed) fix_NL_num[sys.first] = fix_NL_num[sys.first] + 1;
				}
			}
		}
	}

	void t_gallambs::math_dependent_IF_Dd(vector<t_gamb_Dd>& input, int num_Ow_1st, int num_Ow_2nd, vector<t_gamb_Dd>& result, Matrix& vector_space)
	{
		string func_id = "math_dependent_IF_Dd";
		// 确定了双差模糊度的列表，判断独立性
		// check the depend
		int ncol = num_Ow_1st + num_Ow_2nd;

		for (auto& amb : input)
		{
			int index1 = (amb.idx2ow_1st.find(1) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(1)->second + 1) : -1;
			int index2 = (amb.idx2ow_1st.find(2) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(2)->second + 1) : -1;
			int index3 = (amb.idx2ow_1st.find(3) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(3)->second + 1) : -1;
			int index4 = (amb.idx2ow_1st.find(4) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(4)->second + 1) : -1;

			if (index1 < 0 || index2 < 0 || index3 < 0 || index4 < 0) continue;

			// 构建向量空间
			if (vector_space.is_zero())
			{
				vector_space.resize(1, ncol);
				vector_space = 0.0;
				vector_space(1, index1) = 1;
				vector_space(1, index2) = -1;
				vector_space(1, index3 + num_Ow_1st) = -1;
				vector_space(1, index4 + num_Ow_1st) = 1;
				result.push_back(amb);

				if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "Base line is : " +
					int2str(result.size() + 1) + " " +
					amb.rec_1st + " " +
					amb.rec_2nd + " " +
					amb.sat_1st + " " +
					amb.sat_2nd);
				continue;
			}

			RowVector new_vector(0.0, ncol);
			new_vector(index1) = 1;
			new_vector(index2) = -1;
			new_vector(index3 + num_Ow_1st) = -1;
			new_vector(index4 + num_Ow_1st) = 1;

			// 数学原理：奇异矩阵经过SVD分解后，D矩阵对角线的非零元素个数等于A矩阵的矩阵秩
			Matrix A = vector_space & new_vector;
			DiagonalMatrix D;
			if (A.nrows() < A.ncols())
			{
				SVD(A.t(), D);
			}
			else
			{
				SVD(A, D);
			}
			bool isIndepend = (fabs(D(D.nrows(), D.ncols())) > 0.01);

			if (isIndepend)
			{
				vector_space = A;
				result.push_back(amb);
				if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "  Independ Base line is : " +
					int2str(result.size() + 1) + " " +
					amb.rec_1st + " " +
					amb.rec_2nd + " " +
					amb.sat_1st + " " +
					amb.sat_2nd);
			}
		}
	}

	void t_gallambs::math_dependent_1X_Dd(vector<t_gamb_Dd>& input,int num_Ow_1st, int num_Ow_2nd, int num_Ow_3rd, int num_Ow_4th, vector<t_gamb_Dd>& result, Matrix& vector_space)
	{
		string func_id = "math_dependent_1X_Dd";

		int ncol = num_Ow_1st + num_Ow_2nd + num_Ow_3rd + num_Ow_4th;
		int ncol_1st = num_Ow_1st + num_Ow_2nd;

		for (const auto& amb : input)
		{
			int index_1st_1 = (amb.idx2ow_1st.find(1) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(1)->second + 1) : -1;
			int index_1st_2 = (amb.idx2ow_1st.find(2) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(2)->second + 1) : -1;
			int index_1st_3 = (amb.idx2ow_1st.find(3) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(3)->second + 1) : -1;
			int index_1st_4 = (amb.idx2ow_1st.find(4) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(4)->second + 1) : -1;

			int index_2nd_1 = (amb.idx2ow_2nd.find(1) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(1)->second + 1) : -1;
			int index_2nd_2 = (amb.idx2ow_2nd.find(2) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(2)->second + 1) : -1;
			int index_2nd_3 = (amb.idx2ow_2nd.find(3) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(3)->second + 1) : -1;
			int index_2nd_4 = (amb.idx2ow_2nd.find(4) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(4)->second + 1) : -1;

			if (index_1st_1 < 0 || index_1st_2 < 0 || index_1st_3 < 0 || index_1st_4 < 0) continue;
			if (index_2nd_1 < 0 || index_2nd_2 < 0 || index_2nd_3 < 0 || index_2nd_4 < 0) continue;

			// 构建向量空间
			if (vector_space.is_zero())
			{
				vector_space.resize(1, ncol);
				vector_space = 0.0;
				vector_space(1, index_1st_1) = 1;
				vector_space(1, index_1st_2) = -1;
				vector_space(1, index_1st_3 + num_Ow_1st) = -1;
				vector_space(1, index_1st_4 + num_Ow_1st) = 1;

				vector_space(1, ncol_1st + index_2nd_1) = 1;
				vector_space(1, ncol_1st + index_2nd_2) = -1;
				vector_space(1, ncol_1st + index_2nd_3 + num_Ow_3rd) = -1;
				vector_space(1, ncol_1st + index_2nd_4 + num_Ow_3rd) = 1;

				result.push_back(amb);

				if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "Base line is : " +
					int2str(result.size() + 1) + " " +
					amb.rec_1st + " " +
					amb.rec_2nd + " " +
					amb.sat_1st + " " +
					amb.sat_2nd);
				continue;
			}

			RowVector new_vector(0.0, ncol);
			new_vector(index_1st_1) = 1;
			new_vector(index_1st_2) = -1;
			new_vector(index_1st_3 + num_Ow_1st) = -1;
			new_vector(index_1st_4 + num_Ow_1st) = 1;
			new_vector(ncol_1st + index_2nd_1) = 1;
			new_vector(ncol_1st + index_2nd_2) = -1;
			new_vector(ncol_1st + index_2nd_3 + num_Ow_3rd) = -1;
			new_vector(ncol_1st + index_2nd_4 + num_Ow_3rd) = 1;

			// 数学原理：奇异矩阵经过SVD分解后，D矩阵对角线的非零元素个数等于A矩阵的矩阵秩
			Matrix A = vector_space & new_vector;
			DiagonalMatrix D;
			if (A.nrows() < A.ncols())
			{
				SVD(A.t(), D);
			}
			else
			{
				SVD(A, D);
			}
			bool isIndepend = (fabs(D(D.nrows(), D.ncols())) > 0.01);

			if (isIndepend)
			{
				vector_space = A;
				result.push_back(amb);
				if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_INFO, this->_class_id, func_id, "  Independ Base line is : " +
					int2str(result.size() + 1) + " " +
					amb.rec_1st + " " +
					amb.rec_2nd + " " +
					amb.sat_1st + " " +
					amb.sat_2nd);
			}
		}
	}

	void t_gallambs::_select_Sds_value()
	{
		for (auto& sys : this->_nSd_ambs)
		{
			for (auto& bl : sys.second)
			{
				const auto& rec = bl.first;
				if (bl.second.empty()) continue;
				vector<t_gamb_Sd> result;
				vector<t_gamb_Sd> result_fix_WL;
				vector<t_gamb_Sd> result_fix_NO;
				for (auto& amb : bl.second)
				{
					if (amb.NL_fixed && amb.MW_fixed)
						result.push_back(amb);
					if (amb.MW_fixed && !amb.NL_fixed && !double_eq(amb.MW_value, 0.0) && !double_eq(amb.NL_value, 0.0))
					{
						result_fix_WL.push_back(amb);
					}
					if (!amb.MW_fixed && !amb.NL_fixed && !double_eq(amb.MW_value, 0.0) && !double_eq(amb.NL_value, 0.0))
					{
						result_fix_NO.push_back(amb);
					}
				}
				sort(result_fix_WL.begin(), result_fix_WL.end(), sort_Sd);
				sort(result_fix_NO.begin(), result_fix_NO.end(), sort_Sd);
				this->_nSd_ambs_fix_WL[sys.first][rec] = result_fix_WL;
				this->_nSd_ambs_fix_NO[sys.first][rec] = result_fix_NO;

				if (result.empty())
				{
					vector<t_gamb_Sd>().swap(bl.second);
					continue;
				}
				sort(result.begin(), result.end(), sort_Sd);
				bl.second = result;
			}
		}
	}
	
	void t_gallambs::_select_Dds_value(const t_gamb::AMB_ID& ID)
	{
		for (auto& sys : this->_nDd_ambs)
		{
			const auto& amb_id = sys.first.second;
			const auto& sys_key = sys.first;
			if (amb_id != ID) continue;
			for (auto& bl : sys.second)
			{
				const auto& rec_1st = bl.first.first;
				const auto& rec_2nd = bl.first.second;

				if (bl.second.empty()) continue;
				vector<t_gamb_Dd> result;
				vector<t_gamb_Dd> result_fix_WL;
				vector<t_gamb_Dd> result_fix_NO;

				for (auto& amb : bl.second)
				{		
					if (amb.NL_fixed && amb.MW_fixed) result.push_back(amb);

					if (amb.MW_fixed && !amb.NL_fixed && !double_eq(amb.MW_value,0.0) && !double_eq(amb.NL_value,0.0))
					{
						result_fix_WL.push_back(amb);
					}

					if (!amb.MW_fixed && !amb.NL_fixed && !double_eq(amb.MW_value, 0.0) && !double_eq(amb.NL_value, 0.0))
					{
						result_fix_NO.push_back(amb);
					}
				}

				sort(result_fix_WL.begin(), result_fix_WL.end(), sort_Dd);
				sort(result_fix_NO.begin(), result_fix_NO.end(), sort_Dd);
				this->_nDd_ambs_fix_WL[sys.first][bl.first] = result_fix_WL;
				this->_nDd_ambs_fix_NO[sys.first][bl.first] = result_fix_NO;
				if (result.empty())
				{
					vector<t_gamb_Dd>().swap(bl.second);
					continue;
				}
				sort(result.begin(), result.end(), sort_Dd);	
				bl.second = result;
			}
		}
	}

	//void t_gallambs::_check_UDUC_Dds(const GSYS& gsys, const pair<t_gamb::AMB_ID, t_gamb::AMB_ID>& IDs)
	//{
	//	auto key1 = make_pair(gsys, IDs.first);
	//	auto key2 = make_pair(gsys, IDs.second);

	//	double freq1 = t_gamb::LC2FREQ("G01", GOBS_LC::LC_L1);
	//	double freq2 = t_gamb::LC2FREQ("G01", GOBS_LC::LC_L2);
	//	double freq = freq1 - freq2;
	//	double length = CLIGHT / freq;
	//	double length1 = CLIGHT / freq1;
	//	double length2 = CLIGHT / freq2;

	//	if (this->_nDd_ambs.find(key1) == this->_nDd_ambs.end()) return;
	//	if (this->_nDd_ambs.find(key2) == this->_nDd_ambs.end()) return;

	//	auto& Dds_1st = this->_nDd_ambs[key1];
	//	auto& Dds_2nd = this->_nDd_ambs[key2];

	//	for (auto& BL_1st : Dds_1st)
	//	{

	//		auto baseline = BL_1st.first;
	//		if (Dds_2nd.find(baseline) == Dds_2nd.end()) continue;

	//		auto& Bl_Dd_1st = BL_1st.second;
	//		auto& Bl_Dd_2nd = Dds_2nd[baseline];

	//		for (auto& Dd_1st : Bl_Dd_1st)
	//		{
	//			if (!Dd_1st.fixed) continue;
	//			Dd_1st.fixed = false;
	//			//if (double_eq(Dd_1st.value, 0.0)) continue;

	//			for (auto& Dd_2nd : Bl_Dd_2nd)
	//			{
	//				if (!Dd_2nd.fixed) continue;
	//				Dd_2nd.fixed = false;
	//				//if (double_eq(Dd_2nd.value, 0.0)) continue;

	//				if (Dd_1st.rec_1st != Dd_2nd.rec_1st) continue;
	//				if (Dd_1st.rec_2nd != Dd_2nd.rec_2nd) continue;
	//				if (Dd_1st.sat_1st != Dd_2nd.sat_1st) continue;
	//				if (Dd_1st.sat_2nd != Dd_2nd.sat_2nd) continue;

	//				if (Dd_1st.beg != Dd_2nd.beg)continue;
	//				if (Dd_1st.end != Dd_2nd.end)continue;

	//				if (Dd_2nd.MW_fixed == false)
	//				{
	//					Dd_1st.fixed = false;
	//					Dd_2nd.fixed = false;
	//					cout << " Find one is not fix : WL get form obs is not fixed" << endl;
	//					continue;
	//				}

	//				double Dd_WL = (Dd_1st.Ow_value[1] / length1 - Dd_2nd.Ow_value[1] / length2) -
	//					(Dd_1st.Ow_value[2] / length1 - Dd_2nd.Ow_value[2] / length2) -
	//					(Dd_1st.Ow_value[3] / length1 - Dd_2nd.Ow_value[3] / length2) +
	//					(Dd_1st.Ow_value[4] / length1 - Dd_2nd.Ow_value[4] / length2);

	//				double  cycle_WL = Dd_WL;
	//				int    icycle_WL = round(cycle_WL);
	//				double dcycle_WL = cycle_WL - icycle_WL;

	//				if (fabs(dcycle_WL) > 0.14)
	//				{
	//					Dd_1st.fixed = false;
	//					Dd_2nd.fixed = false;
	//					cout << " Find one is not fix : fabs(dcycle) >= 0.14" << endl;
	//					continue;
	//				}

	//				if (icycle_WL != round(Dd_2nd.MW_value))
	//				{
	//					Dd_1st.fixed = false;
	//					Dd_2nd.fixed = false;
	//					cout << " Find one is not fix : fabs(MW) != fabs(L1-L2)" << endl;
	//					cout << " fabs(MW) - fabs(L1-L2) = " << cycle_WL - Dd_2nd.MW_value << endl;
	//					continue;
	//				}

	//				Dd_1st.fixed = true;
	//				Dd_2nd.fixed = true;
	//				break;
	//			}
	//		}
	//	}
	//}

	void t_gallambs::_get_Dds_value(const t_gamb::AMB_ID& ID)
	{
		for (auto& sys : this->_nDd_ambs)
		{
			const auto& sys_id_key = sys.first;
			const auto& sys_key    = sys.first.first;
			const auto& id_key     = sys.first.second;
			if (id_key == t_gamb::AMB_ID::UNDEF) continue;
			if (id_key != ID) continue;

			for (auto& bl_ambs : sys.second)
			{
				const auto& rec_1st = bl_ambs.first.first;
				const auto& rec_2nd = bl_ambs.first.second;
				const auto& bl = bl_ambs.first;

				for (auto& amb : bl_ambs.second)
				{
					if (ID == t_gamb::AMB_ID::IF_L1_L2 || ID == t_gamb::AMB_ID::IF_L1_L3)
					{
						if (!_get_IF_Dd_value(sys_id_key, bl, amb)) continue;
						if (!_get_MW_value(sys_id_key, bl, amb)) continue;
						if (!_get_NL_value(sys_id_key, bl, amb)) continue;
					}
					else if (ID == t_gamb::AMB_ID::ALL_L1_L2 || ID == t_gamb::AMB_ID::ALL_L1_L3)
					{
						if (!_get_1X_Dd_value(sys_id_key, bl, amb)) continue;
						if (_dd_mode != DD_MODEL::RAW_CB && !_get_MW_value(   sys_id_key, bl, amb)) continue;
						if (_dd_mode == DD_MODEL::RAW_CB && !_get_WL_value(   sys_id_key, bl, amb)) continue;
						if (!_get_NL_value(   sys_id_key, bl, amb)) continue;
					}
					else if (ID == t_gamb::AMB_ID::ALL_L1 || ID == t_gamb::AMB_ID::ALL_L2 || ID == t_gamb::AMB_ID::ALL_L3)
					{
						if (!_get_LX_Dd_value(sys_id_key, bl, amb)) continue;
					}
					else
					{
						throw logic_error("wrong type of ambs while _get_Dds_value");
					}
				}// loop baseline
			}// loop sys
		} // loop all
	}

	//void t_gallambs::_get_WL_value(const pair<GSYS, t_gamb::AMB_ID>& ID1, const pair<GSYS, t_gamb::AMB_ID>& ID2)
	//{
	//	// 这两个ID没有模糊度
	//	if (this->_nDd_ambs.find(ID1) == this->_nDd_ambs.end()) return;
	//	if (this->_nDd_ambs.find(ID2) == this->_nDd_ambs.end()) return;

	//	// 确保ID1 < ID2
	//	if (ID1.second >= ID2.second) return;

	//	for (auto& Dds_for_1st_bl : this->_nDd_ambs[ID1])
	//	{
	//		const auto& bl_1st = Dds_for_1st_bl.first;
	//		for (auto& Dds_for_2nd_bl : this->_nDd_ambs[ID2])
	//		{
	//			const auto& bl_2nd = Dds_for_2nd_bl.first;
	//			if (bl_1st != bl_2nd) continue;

	//			// get all the amb with the baseline;
	//			for (auto& Dd1 : Dds_for_1st_bl.second)
	//			{
	//				for (auto& Dd2 : Dds_for_2nd_bl.second)
	//				{
	//					_get_WL_value(ID1, ID2, bl_1st, Dd1, Dd2);
	//				}
	//			}
	//		}
	//	}
	//}

	bool t_gallambs::_get_NL_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& amb)
	{
		if (index_key.second == t_gamb::AMB_ID::IF_L1_L2 || index_key.second == t_gamb::AMB_ID::IF_L1_L3)
		{
			if (!_get_IF_NL_value(index_key, bl, amb)) return false;
		}
		else if (index_key.second == t_gamb::AMB_ID::ALL_L1_L2 ||
			index_key.second == t_gamb::AMB_ID::ALL_L1_L3)
		{
			if (!_get_1X_NL_value(index_key, bl, amb)) return false;
		}
		else if(index_key.second == t_gamb::AMB_ID::ALL_L1 ||
			index_key.second == t_gamb::AMB_ID::ALL_L2 ||
			index_key.second == t_gamb::AMB_ID::ALL_L3)
		{
			return true;
		}
		else
		{
			return false;
		}

		return true;
	}

	bool t_gallambs::_get_IF_NL_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& amb)
	{
		const auto& rec_1st = bl.first;
		const auto& rec_2nd = bl.second;

		if (double_eq(amb.value, 0.0)) return false;
		if (double_eq(amb.MW_value, 0.0)) return false;

		// jdhuang : fix for muti-GNSS IF combination
		t_gobsgnss gnss;
		gnss.sat(amb.sat_1st);
		auto gsys = t_gsys::sat2gsys(amb.sat_1st);

		double freq_1st = 0.0;
		double freq_2nd = 0.0;

		if (index_key.second == t_gamb::AMB_ID::IF_L1_L2)
		{
			freq_1st = gnss.frequency(_band_index[gsys][FREQ_1]);
			freq_2nd = gnss.frequency(_band_index[gsys][FREQ_2]);
		}
		else if (index_key.second == t_gamb::AMB_ID::IF_L1_L3)
		{
			freq_1st = gnss.frequency(_band_index[gsys][FREQ_1]);
			freq_2nd = gnss.frequency(_band_index[gsys][FREQ_3]);
		}
		else
		{
			// no support for others
			return false;
		}

	
		//double freq_1st = t_gamb::LC2FREQ(amb.sat_1st, GOBS_LC::LC_L1);
		//double freq_2nd = t_gamb::LC2FREQ(amb.sat_1st, GOBS_LC::LC_L2);

		double length_1st = CLIGHT / freq_1st;
		double length_2nd = CLIGHT / freq_2nd;

		double length_3rd = length_1st * length_2nd / (length_2nd - length_1st);
		double length_4th = length_1st * length_2nd / (length_2nd + length_1st);

		amb.NL_value = amb.value / (length_4th)-round(amb.MW_value) * length_3rd / length_2nd;
		amb.NL_sigma = amb.sigma / length_4th;

		if (double_eq(0.0, amb.NL_sigma)) amb.NL_sigma = 0.05;
		amb.NL_fixed = (_evaluator_NL) ? _evaluator_NL->bdeci(amb.NL_value, amb.NL_sigma, 1) : false;
		if (amb.NL_fixed)
		{
			amb.NL_fixed = (fabs(round(amb.NL_value) - amb.NL_value) < 0.15) ? true : false;
		}
		amb.NL_fixed = (amb.MW_fixed && amb.NL_fixed);
		return true;
	}

	bool t_gallambs::_get_1X_NL_value(const pair<GSYS, t_gamb::AMB_ID>& index_key_input, const pair<string, string>& bl, t_gamb_Dd& amb)
	{
		const auto& rec_1st = bl.first;
		const auto& rec_2nd = bl.second;
		pair<GSYS, t_gamb::AMB_ID> index_ID1;
		pair<GSYS, t_gamb::AMB_ID> index_ID2;

		if(index_key_input.second == t_gamb::AMB_ID::ALL_L1_L2) 
		{
			index_ID1 = make_pair(index_key_input.first, t_gamb::AMB_ID::ALL_L1);
			index_ID2 = make_pair(index_key_input.first, t_gamb::AMB_ID::ALL_L2);
		}
		else if (index_key_input.second == t_gamb::AMB_ID::ALL_L1_L3)
		{
			index_ID1 = make_pair(index_key_input.first, t_gamb::AMB_ID::ALL_L1);
			index_ID2 = make_pair(index_key_input.first, t_gamb::AMB_ID::ALL_L3);

		}
		else
		{
			return false;
		}

		if (amb.DD_LX_value.size() < 2) return false;
		if (double_eq(amb.MW_value, 0.0)) return false;
		if (double_eq(amb.DD_LX_value[index_ID1.second], 0.0)) return false;
		if (double_eq(amb.DD_LX_value[index_ID2.second], 0.0)) return false;

		//lijie 
		t_gobsgnss gnss;
		gnss.sat(amb.sat_1st);
		auto gsys = t_gsys::sat2gsys(amb.sat_1st);
		double freq_1st = gnss.frequency(_band_index[gsys][FREQ_1]);
		double freq_2nd;
			 if (index_ID2.second == t_gamb::AMB_ID::ALL_L2) freq_2nd = gnss.frequency(_band_index[gsys][FREQ_2]);
		else if (index_ID2.second == t_gamb::AMB_ID::ALL_L3) freq_2nd = gnss.frequency(_band_index[gsys][FREQ_3]);
		else return false; 

		double length_1st = CLIGHT / freq_1st;
		double length_2nd = CLIGHT / freq_2nd;

		double apha =  freq_1st * freq_1st / (freq_1st * freq_1st - freq_2nd * freq_2nd);
		double beta = -freq_2nd * freq_2nd / (freq_1st * freq_1st - freq_2nd * freq_2nd);
		
		amb.NL_value = (freq_1st + freq_2nd) / CLIGHT * (apha * amb.DD_LX_value[index_ID1.second] + beta * amb.DD_LX_value[index_ID2.second]) 
			- freq_1st / (freq_1st - freq_2nd) * round(amb.MW_value);

		amb.NL_sigma = 0.05;
		amb.NL_fixed = (_evaluator_NL) ? _evaluator_NL->bdeci(amb.NL_value, amb.NL_sigma, 1) : false;

		if (amb.NL_fixed)
		{
			amb.NL_fixed = (fabs(round(amb.NL_value) - amb.NL_value) < 0.15) ? true : false;
		}

		amb.NL_fixed = (amb.MW_fixed && amb.NL_fixed);

		return true;
	}

	//bool t_gallambs::_get_WL_value(const pair<GSYS, t_gamb::AMB_ID>& index_key1, const pair<GSYS, t_gamb::AMB_ID>& index_key2, 
	//	const pair<string, string>& bl, t_gamb_Dd& amb1, t_gamb_Dd& amb2)
	//{
	//	if (index_key1.first != index_key2.first) return false;
	//	if (!amb2.MW_fixed) return false;
	//	if (amb1.rec_1st != amb2.rec_1st)         return false;
	//	if (amb1.rec_2nd != amb2.rec_2nd)         return false;
	//	if (amb1.sat_1st != amb2.sat_1st)         return false;
	//	if (amb1.sat_2nd != amb2.sat_2nd)         return false;

	//	if (amb1.beg != amb2.beg) return false;
	//	if (amb1.end != amb2.end) return false;

	//	const auto& rec_1st = bl.first;
	//	const auto& rec_2nd = bl.second;

	//	int amb1_index1 = (amb1.idx2ow_1st.find(1) != amb1.idx2ow_1st.end()) ? (amb1.idx2ow_1st.find(1)->second) : -1;
	//	int amb1_index2 = (amb1.idx2ow_1st.find(2) != amb1.idx2ow_1st.end()) ? (amb1.idx2ow_1st.find(2)->second) : -1;
	//	int amb1_index3 = (amb1.idx2ow_1st.find(3) != amb1.idx2ow_1st.end()) ? (amb1.idx2ow_1st.find(3)->second) : -1;
	//	int amb1_index4 = (amb1.idx2ow_1st.find(4) != amb1.idx2ow_1st.end()) ? (amb1.idx2ow_1st.find(4)->second) : -1;

	//	int amb2_index1 = (amb2.idx2ow_1st.find(1) != amb2.idx2ow_1st.end()) ? (amb2.idx2ow_1st.find(1)->second) : -1;
	//	int amb2_index2 = (amb2.idx2ow_1st.find(2) != amb2.idx2ow_1st.end()) ? (amb2.idx2ow_1st.find(2)->second) : -1;
	//	int amb2_index3 = (amb2.idx2ow_1st.find(3) != amb2.idx2ow_1st.end()) ? (amb2.idx2ow_1st.find(3)->second) : -1;
	//	int amb2_index4 = (amb2.idx2ow_1st.find(4) != amb2.idx2ow_1st.end()) ? (amb2.idx2ow_1st.find(4)->second) : -1;

	//	if (amb1_index1 < 0 || amb1_index2 < 0 || amb1_index3 < 0 || amb1_index4 < 0) return false;
	//	if (amb2_index1 < 0 || amb2_index2 < 0 || amb2_index3 < 0 || amb2_index4 < 0) return false;

	//	if (double_eq(amb1.value, 0.0)) return false;
	//	if (double_eq(amb2.value, 0.0)) return false;

	//	//lijie
	//	t_gobsgnss gnss;
	//	gnss.sat(amb1.sat_1st);
	//	auto gsys = t_gsys::sat2gsys(amb1.sat_1st);
	//	double freq1 = gnss.frequency(_band_index[gsys][FREQ_1]);
	//	double freq2 = gnss.frequency(_band_index[gsys][FREQ_2]);

	//	/*double freq1 = t_gamb::LC2FREQ(amb1.sat_1st, GOBS_LC::LC_L1);
	//	double freq2 = t_gamb::LC2FREQ(amb1.sat_1st, GOBS_LC::LC_L2);*/
	//	double freq  = freq1 - freq2;
	//	double length1 = CLIGHT / freq1;
	//	double length2 = CLIGHT / freq2;

	//	if (double_eq(amb1.Ow_value[1], 0.0)) return false;
	//	if (double_eq(amb1.Ow_value[2], 0.0)) return false;
	//	if (double_eq(amb1.Ow_value[3], 0.0)) return false;
	//	if (double_eq(amb1.Ow_value[4], 0.0)) return false;

	//	if (double_eq(amb2.Ow_value[1], 0.0)) return false;
	//	if (double_eq(amb2.Ow_value[2], 0.0)) return false;
	//	if (double_eq(amb2.Ow_value[3], 0.0)) return false;
	//	if (double_eq(amb2.Ow_value[4], 0.0)) return false;

	//	double Dd_WL = (amb1.Ow_value[1] / length1 - amb2.Ow_value[1] / length2) -
	//		           (amb1.Ow_value[2] / length1 - amb2.Ow_value[2] / length2) -
	//		           (amb1.Ow_value[3] / length1 - amb2.Ow_value[3] / length2) +
	//		           (amb1.Ow_value[4] / length1 - amb2.Ow_value[4] / length2);

	//	bool WL_fix = fabs(Dd_WL - round(Dd_WL)) < 0.15 ? true : false;
	//	//if (WL_fix)
	//	//{
	//	//	WL_fix = ((fabs(amb2.MW_value - Dd_WL) < 0.15) && (round(amb2.MW_value) == round(Dd_WL)))? true : false;
	//	//}

	//	amb1.Mw_value = amb2.Mw_value;
	//	amb1.MW_value = amb2.MW_value;
	//	amb1.MW_sigma = amb2.MW_sigma;
	//	amb1.MW_fixed = amb2.MW_fixed;
	//	amb1.WL_value = Dd_WL;
	//	amb2.WL_value = Dd_WL;
	//	amb1.WL_sigma = 0.0;
	//	amb2.WL_sigma = 0.0;
	//	amb1.WL_fixed = WL_fix;
	//	amb2.WL_fixed = WL_fix;

	//	//cout << "WL value             : " << setw(20) << right << amb1.WL_value 
	//	//	<< " Dd_WL - round(Dd_WL) : " << setw(20) << right << Dd_WL - round(Dd_WL) 
	//	//	<< " Dd_WL - Dd_MW        : " << setw(20) << right << amb1.WL_value - amb1.MW_value
	//	//	<< " Dd_WL_amb1 - Dd_MW_1 : " << setw(20) << right << amb1.Ow_value[1] / length1 - amb2.Ow_value[1] / length2 - amb2.Mw_value[1]
	//	//	<< " Dd_WL_amb2 - Dd_MW_2 : " << setw(20) << right << amb1.Ow_value[2] / length1 - amb2.Ow_value[2] / length2 - amb2.Mw_value[2]
	//	//	<< " Dd_WL_amb3 - Dd_MW_3 : " << setw(20) << right << amb1.Ow_value[3] / length1 - amb2.Ow_value[3] / length2 - amb2.Mw_value[3]
	//	//	<< " Dd_WL_amb4 - Dd_MW_4 : " << setw(20) << right << amb1.Ow_value[4] / length1 - amb2.Ow_value[4] / length2 - amb2.Mw_value[4]
	//	//	<< endl;
	//	return true;
	//}

	void t_gallambs::_get_Sds_value(t_gsetamb::UPD_MODE upd_mode)
	{
		for (auto& sys : this->_nSd_ambs)
		{
			const auto& sys_key = sys.first;
			const auto& sys_ID = sys.first.first;
			const auto& ID_key = sys.first.second;
			if (ID_key == t_gamb::AMB_ID::UNDEF) continue;
			for (auto& bl : sys.second)
			{
				const auto& rec = bl.first;
				for (auto& amb : bl.second)
				{
					int index1 = (amb.idx2ow_1st.find(1) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(1)->second) : -1;
					int index2 = (amb.idx2ow_1st.find(2) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(2)->second) : -1;
					if (index1 < 0 || index2 < 0) continue;
					if (double_eq(this->_nOw_ambs[sys_key][rec][index1].value, 0.0)) continue;
					if (double_eq(this->_nOw_ambs[sys_key][rec][index2].value, 0.0)) continue;
					amb.Ow_value[1] = this->_nOw_ambs[sys_key][rec][index1].value;
					amb.Ow_value[2] = this->_nOw_ambs[sys_key][rec][index2].value;
					amb.Ow_time[1] = make_pair(this->_nOw_ambs[sys_key][rec][index1].beg, this->_nOw_ambs[sys_key][rec][index1].end);
					amb.Ow_time[2] = make_pair(this->_nOw_ambs[sys_key][rec][index2].beg, this->_nOw_ambs[sys_key][rec][index2].end);
					amb.value = this->_nOw_ambs[sys_key][rec][index1].value -
						this->_nOw_ambs[sys_key][rec][index2].value;
					amb.sigma = 0.0;
					if (double_eq(amb.value, 0.0)) continue;
					if (ID_key == t_gamb::AMB_ID::IF_L1_L2)
					{
						double mw_value_1st = 0.0; double mw_sigma_1st = 0.0;
						double mw_value_2nd = 0.0; double mw_sigma_2nd = 0.0;
						auto beg_1st = this->time_table.lower_bound(this->_nOw_ambs[sys_key][rec][index1].beg);
						auto beg_2nd = this->time_table.lower_bound(this->_nOw_ambs[sys_key][rec][index2].beg);
						auto end_1st = this->time_table.lower_bound(this->_nOw_ambs[sys_key][rec][index1].end);
						auto end_2nd = this->time_table.lower_bound(this->_nOw_ambs[sys_key][rec][index2].end);
						//cout << amb.dbl_obj_1st << "MW" << endl;
						_get_MW(amb.sgl_obj, amb.dbl_obj_1st, 2, this->time_table, beg_1st, end_1st, mw_value_1st, mw_sigma_1st);
						//cout << amb.sgl_obj << "  " << amb.dbl_obj_1st << "  " << "AMBC" << beg_1st->str_mjdsod() << end_1st->str_mjdsod() << "  " << setw(20) << setprecision(15) << mw_value_1st << "  " << mw_sigma_1st << endl;
						//cout << amb.dbl_obj_2nd << "MW" << endl;
						_get_MW(amb.sgl_obj, amb.dbl_obj_2nd, 2, this->time_table, beg_2nd, end_2nd, mw_value_2nd, mw_sigma_2nd);
						//cout << amb.sgl_obj << "  " << amb.dbl_obj_2nd << "  " << "AMBC" << beg_2nd->str_mjdsod() << end_2nd->str_mjdsod() << "  " << setw(20) << setprecision(15) << mw_value_2nd << "  " << mw_sigma_2nd << endl;
						if (double_eq(mw_value_1st, 0.0)) continue;
						if (double_eq(mw_value_2nd, 0.0)) continue;
						amb.MW_value = mw_value_1st -
							mw_value_2nd;
						amb.MW_sigma = sqrt(mw_sigma_1st * mw_sigma_1st +
							mw_sigma_2nd * mw_sigma_2nd);
						if (upd_mode != t_gsetamb::UPD_MODE::OSB) _add_upd("WL", amb);
						//cout << " WL " << amb.dbl_obj_1st << amb.dbl_obj_2nd << setw(20) << setprecision(15) << amb.WL_value <<" " << amb.WL_sigma << endl;
						
						//lijie
						t_gobsgnss gnss;
						gnss.sat(amb.dbl_obj_1st);
						auto gsys = t_gsys::sat2gsys(amb.dbl_obj_2nd);
						double freq_1st = gnss.frequency(_band_index[gsys][FREQ_1]);
						double freq_2nd = gnss.frequency(_band_index[gsys][FREQ_2]);
						/*double freq_1st = t_gamb::LC2FREQ(amb.dbl_obj_1st, GOBS_LC::LC_L1);
						double freq_2nd = t_gamb::LC2FREQ(amb.dbl_obj_2nd, GOBS_LC::LC_L2);*/
						double length_1st = CLIGHT / freq_1st;
						double length_2nd = CLIGHT / freq_2nd;
						double length_3rd = length_1st * length_2nd / (length_2nd - length_1st);
						double length_4th = length_1st * length_2nd / (length_2nd + length_1st);
						amb.NL_value = amb.value / (length_4th)-round(amb.MW_value) * length_3rd / length_2nd;
						amb.NL_sigma = amb.sigma / length_4th;
						if (upd_mode == t_gsetamb::UPD_MODE::UPD) _add_upd("NL", amb);
						/*cout << " NL fact" << setw(20) << setprecision(15) << length_4th << " " << length_3rd << " " << length_2nd << " value " << amb.value << " " << round(amb.WL_value) << endl;
						cout << " NL " << amb.dbl_obj_1st << amb.dbl_obj_2nd << setw(20) << setprecision(15) << amb.NL_value << " " << amb.NL_sigma << endl;
						cout << "WL ";*/
						if (double_eq(0.0, amb.MW_sigma)) amb.MW_sigma = 999.0;
						amb.MW_fixed = (_evaluator_MW) ? _evaluator_MW->bdeci(amb.MW_value, amb.MW_sigma, 1) : false;
						/*cout << "NL ";*/
						
						if (double_eq(0.0, amb.NL_sigma)) amb.NL_sigma = 0.05;
						amb.NL_fixed = (_evaluator_NL) ? _evaluator_NL->bdeci(amb.NL_value, amb.NL_sigma, 1) : false;
						/*if (amb.WL_fixed == true && amb.NL_fixed == true)
						{
							cout << "AMBC" << "true   "  << amb.dbl_obj_1st << " " << amb.dbl_obj_2nd << " " << amb.beg.str_mjdsod() << amb.end.str_mjdsod() << "  " << setw(20) << setprecision(15) << amb.NL_value << amb.NL_sigma << " " << setw(20) << setprecision(15) << amb.WL_value << " " << amb.WL_sigma << endl;
						}
						else
						{
							cout << "AMBC" << "false   "  << amb.dbl_obj_1st << " " << amb.dbl_obj_2nd << " " << amb.beg.str_mjdsod() << amb.end.str_mjdsod() << "  " << setw(20) << setprecision(15) << amb.NL_value << amb.NL_sigma << " " << setw(20) << setprecision(15) << amb.WL_value << " " << amb.WL_sigma << endl;
						}*/
					}
					else
					{
						const auto& sat = this->_nOw_ambs[sys_key][rec][index1].sat;
						double length = CLIGHT;
						//lijie
						t_gobsgnss gnss;
						gnss.sat(amb.dbl_obj_1st);
						auto gsys = t_gsys::sat2gsys(amb.dbl_obj_2nd);
						switch (ID_key)
						{
						case t_gamb::AMB_ID::ALL_L1: length = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_1]); break;
						case t_gamb::AMB_ID::ALL_L2: length = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_2]); break;
						case t_gamb::AMB_ID::ALL_L3: length = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_3]); break;
						case t_gamb::AMB_ID::ALL_L4: length = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_4]); break;
						case t_gamb::AMB_ID::ALL_L5: length = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_5]); break;
						default: break;
						}
						/*switch (ID_key)
						{
						case t_gamb::AMB_ID::ALL_L1: length = CLIGHT / t_gamb::LC2FREQ(sat, GOBS_LC::LC_L1); break;
						case t_gamb::AMB_ID::ALL_L2: length = CLIGHT / t_gamb::LC2FREQ(sat, GOBS_LC::LC_L2); break;
						case t_gamb::AMB_ID::ALL_L3: length = CLIGHT / t_gamb::LC2FREQ(sat, GOBS_LC::LC_L3); break;
						case t_gamb::AMB_ID::ALL_L4: length = CLIGHT / t_gamb::LC2FREQ(sat, GOBS_LC::LC_L4); break;
						case t_gamb::AMB_ID::ALL_L5: length = CLIGHT / t_gamb::LC2FREQ(sat, GOBS_LC::LC_L5); break;
						default: break;
						}*/
						if (double_eq(length, 0.0)) continue;
						double dcycle = amb.value / length;
						int    icycle = round(dcycle);
						amb.value = icycle * length;
						amb.fixed = fabs(dcycle - icycle) < 0.15 ? true : false;
						if (!amb.fixed) continue;
						string mark = " ALL";
						string sLC;
						if (amb.ID == t_gamb::AMB_ID::ALL_L1) sLC = "L1";
						else if (amb.ID == t_gamb::AMB_ID::ALL_L2) sLC = "L2";
						else if (amb.ID == t_gamb::AMB_ID::ALL_L3) sLC = "L3";
						else if (amb.ID == t_gamb::AMB_ID::ALL_L4) sLC = "L4";
						else if (amb.ID == t_gamb::AMB_ID::ALL_L5) sLC = "L5";
						else sLC = "UK";
						string fix = amb.fixed ? "true" : "false";
						cout << setw(6) << left << mark
							<< setw(6) << left << sLC
							<< setw(6) << amb.dbl_obj_1st
							<< setw(6) << amb.dbl_obj_2nd
							<< setw(6) << amb.sgl_obj;
						cout << setw(10) << right << amb.beg.mjd() << " "
							<< setw(16) << fixed << setprecision(8) << amb.beg.sod() + amb.beg.dsec();
						cout << setw(10) << right << amb.end.mjd() << " "
							<< setw(16) << fixed << setprecision(8) << amb.end.sod() + amb.end.dsec();
						cout << setw(16) << right << fixed << setprecision(8) << amb.Ow_value.find(1)->second
							<< setw(16) << right << fixed << setprecision(8) << amb.Ow_value.find(2)->second
							<< setw(16) << right << fixed << setprecision(8) << amb.value
							<< setw(16) << right << fixed << setprecision(8) << length
							<< setw(16) << right << fixed << setprecision(8) << dcycle
							<< setw(16) << right << fixed << setprecision(8) << icycle
							<< setw(16) << right << fixed << setprecision(8) << dcycle - icycle
							<< setw(16) << right << fixed << setprecision(8) << dcycle * length
							<< setw(16) << right << fixed << setprecision(8) << icycle * length
							<< setw(16) << right << fixed << setprecision(8) << amb.sigma << " " << fix << endl;
					}
				}// loop baseline
			}// loop sys
		} // loop all
	}

	bool t_gallambs::_get_MW_value(const pair<GSYS, t_gamb::AMB_ID>& index_key_input, const pair<string, string>& bl, t_gamb_Dd& amb)
	{
		const auto& rec_1st = bl.first;
		const auto& rec_2nd = bl.second;
		
		pair<GSYS, t_gamb::AMB_ID> index_key;
		if (index_key_input.second == t_gamb::AMB_ID::IF_L1_L2 || index_key_input.second == t_gamb::AMB_ID::IF_L1_L3)
		{
			index_key = index_key_input;
		}
		else if (index_key_input.second == t_gamb::AMB_ID::ALL_L1_L2) index_key = make_pair(index_key_input.first, t_gamb::AMB_ID::ALL_L2);
		else if (index_key_input.second == t_gamb::AMB_ID::ALL_L1_L3) index_key = make_pair(index_key_input.first, t_gamb::AMB_ID::ALL_L3);
		else throw logic_error("error amb ID");
			
		int index1 = -1;
		int index2 = -1;
		int index3 = -1;
		int index4 = -1;

		if (index_key_input.second == t_gamb::AMB_ID::ALL_L1_L2 || index_key_input.second == t_gamb::AMB_ID::ALL_L1_L3)
		{
			index1 = (amb.idx2ow_2nd.find(1) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(1)->second) : -1;
			index2 = (amb.idx2ow_2nd.find(2) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(2)->second) : -1;
			index3 = (amb.idx2ow_2nd.find(3) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(3)->second) : -1;
			index4 = (amb.idx2ow_2nd.find(4) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(4)->second) : -1;

		}
		else if (index_key_input.second == t_gamb::AMB_ID::IF_L1_L2 || index_key_input.second == t_gamb::AMB_ID::IF_L1_L3)
		{
			index1 = (amb.idx2ow_1st.find(1) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(1)->second) : -1;
			index2 = (amb.idx2ow_1st.find(2) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(2)->second) : -1;
			index3 = (amb.idx2ow_1st.find(3) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(3)->second) : -1;
			index4 = (amb.idx2ow_1st.find(4) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(4)->second) : -1;
		}
		else
		{
			throw logic_error("The amb ID is Wrong, it should be IF_L1_L2, ALL_L1_L2 and ALL_L1_L3.");
		}

		if (index1 < 0 || index2 < 0 || index3 < 0 || index4 < 0) return false;

		const auto& amb1 = this->_nOw_ambs[index_key][rec_1st][index1];
		const auto& amb2 = this->_nOw_ambs[index_key][rec_1st][index2];
		const auto& amb3 = this->_nOw_ambs[index_key][rec_2nd][index3];
		const auto& amb4 = this->_nOw_ambs[index_key][rec_2nd][index4];

		double mw_value_1st = 0.0; double mw_sigma_1st = 0.0;
		double mw_value_2nd = 0.0; double mw_sigma_2nd = 0.0;
		double mw_value_3rd = 0.0; double mw_sigma_3rd = 0.0;
		double mw_value_4th = 0.0; double mw_sigma_4th = 0.0;

		auto beg_1st = this->time_table.lower_bound(amb1.beg);
		auto beg_2nd = this->time_table.lower_bound(amb2.beg);
		auto beg_3rd = this->time_table.lower_bound(amb3.beg);
		auto beg_4th = this->time_table.lower_bound(amb4.beg);

		auto end_1st = this->time_table.lower_bound(amb1.end);
		auto end_2nd = this->time_table.lower_bound(amb2.end);
		auto end_3rd = this->time_table.lower_bound(amb3.end);
		auto end_4th = this->time_table.lower_bound(amb4.end);

		int freq = (index_key_input.second == t_gamb::AMB_ID::ALL_L1_L2 || index_key_input.second == t_gamb::AMB_ID::IF_L1_L2) ? 2 : 3;
		_get_MW(amb.rec_1st, amb.sat_1st, freq, this->time_table, beg_1st, end_1st, mw_value_1st, mw_sigma_1st);
		_get_MW(amb.rec_1st, amb.sat_2nd, freq, this->time_table, beg_2nd, end_2nd, mw_value_2nd, mw_sigma_2nd);
		_get_MW(amb.rec_2nd, amb.sat_1st, freq, this->time_table, beg_3rd, end_3rd, mw_value_3rd, mw_sigma_3rd);
		_get_MW(amb.rec_2nd, amb.sat_2nd, freq, this->time_table, beg_4th, end_4th, mw_value_4th, mw_sigma_4th);

		if (double_eq(mw_value_1st, 0.0)) return false;
		if (double_eq(mw_value_2nd, 0.0)) return false;
		if (double_eq(mw_value_3rd, 0.0)) return false;
		if (double_eq(mw_value_4th, 0.0)) return false;

		
		amb.MW_value = mw_value_1st -
			mw_value_2nd -
			mw_value_3rd +
			mw_value_4th;

		amb.MW_sigma = sqrt(mw_sigma_1st * mw_sigma_1st +
			mw_sigma_2nd * mw_sigma_2nd +
			mw_sigma_3rd * mw_sigma_3rd +
			mw_sigma_4th * mw_sigma_4th);

		amb.Mw_value[1] = mw_value_1st;
		amb.Mw_value[2] = mw_value_2nd;
		amb.Mw_value[3] = mw_value_3rd;
		amb.Mw_value[4] = mw_value_4th;
		amb.Mw_sigma[1] = mw_sigma_1st;
		amb.Mw_sigma[2] = mw_sigma_2nd;
		amb.Mw_sigma[3] = mw_sigma_3rd;
		amb.Mw_sigma[4] = mw_sigma_4th;


		if (double_eq(0.0, amb.MW_sigma)) amb.MW_sigma = 999.0;
		amb.MW_fixed = (_evaluator_MW) ? _evaluator_MW->bdeci(amb.MW_value, amb.MW_sigma, 1) : false;

		if (amb.MW_fixed)
		{
			amb.MW_fixed = (fabs(round(amb.MW_value) - amb.MW_value) < 0.15) ? true : false;
		}
	
		return true;
	}

	bool t_gallambs::_get_WL_value(const pair<GSYS, t_gamb::AMB_ID>& index_key_input, const pair<string, string>& bl, t_gamb_Dd & amb)
	{
		const auto& rec_1st = bl.first;
		const auto& rec_2nd = bl.second;

		pair<GSYS, t_gamb::AMB_ID> index_key1;
		pair<GSYS, t_gamb::AMB_ID> index_key2;

		if (index_key_input.second == t_gamb::AMB_ID::ALL_L1_L2)
		{
			index_key1 = make_pair(index_key_input.first, t_gamb::AMB_ID::ALL_L1);
			index_key2 = make_pair(index_key_input.first, t_gamb::AMB_ID::ALL_L2);
		}
		else if (index_key_input.second == t_gamb::AMB_ID::ALL_L1_L3)
		{
			index_key1 = make_pair(index_key_input.first, t_gamb::AMB_ID::ALL_L1);
			index_key2 = make_pair(index_key_input.first, t_gamb::AMB_ID::ALL_L3);
		}
		else throw logic_error("error amb ID");

		int amb1_index1 = -1;
		int amb1_index2 = -1;
		int amb1_index3 = -1;
		int amb1_index4 = -1;

		int amb2_index1 = -1;
		int amb2_index2 = -1;
		int amb2_index3 = -1;
		int amb2_index4 = -1;

		if (index_key_input.second == t_gamb::AMB_ID::ALL_L1_L2 || index_key_input.second == t_gamb::AMB_ID::ALL_L1_L3)
		{
			amb1_index1 = (amb.idx2ow_1st.find(1) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(1)->second) : -1;
			amb1_index2 = (amb.idx2ow_1st.find(2) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(2)->second) : -1;
			amb1_index3 = (amb.idx2ow_1st.find(3) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(3)->second) : -1;
			amb1_index4 = (amb.idx2ow_1st.find(4) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(4)->second) : -1;
					
			amb2_index1 = (amb.idx2ow_2nd.find(1) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(1)->second) : -1;
			amb2_index2 = (amb.idx2ow_2nd.find(2) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(2)->second) : -1;
			amb2_index3 = (amb.idx2ow_2nd.find(3) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(3)->second) : -1;
			amb2_index4 = (amb.idx2ow_2nd.find(4) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(4)->second) : -1;

		}
		else
		{
			throw logic_error("The amb ID is Wrong, it should be IF_L1_L2, ALL_L1_L2 and ALL_L1_L3.");
		}

		if (amb1_index1 < 0 || amb1_index2 < 0 || amb1_index3 < 0 || amb1_index4 < 0) return false;
		if (amb2_index1 < 0 || amb2_index2 < 0 || amb2_index3 < 0 || amb2_index4 < 0) return false;

		const auto& L1_amb1 = this->_nOw_ambs[index_key1][rec_1st][amb1_index1];
		const auto& L1_amb2 = this->_nOw_ambs[index_key1][rec_1st][amb1_index2];
		const auto& L1_amb3 = this->_nOw_ambs[index_key1][rec_2nd][amb1_index3];
		const auto& L1_amb4 = this->_nOw_ambs[index_key1][rec_2nd][amb1_index4];

		const auto& LX_amb1 = this->_nOw_ambs[index_key2][rec_1st][amb2_index1];
		const auto& LX_amb2 = this->_nOw_ambs[index_key2][rec_1st][amb2_index2];
		const auto& LX_amb3 = this->_nOw_ambs[index_key2][rec_2nd][amb2_index3];
		const auto& LX_amb4 = this->_nOw_ambs[index_key2][rec_2nd][amb2_index4];

		if (double_eq(L1_amb1.value, 0.0)) return false;
		if (double_eq(L1_amb2.value, 0.0)) return false;
		if (double_eq(L1_amb3.value, 0.0)) return false;
		if (double_eq(L1_amb4.value, 0.0)) return false;
		if (double_eq(LX_amb1.value, 0.0)) return false;
		if (double_eq(LX_amb2.value, 0.0)) return false;
		if (double_eq(LX_amb3.value, 0.0)) return false;
		if (double_eq(LX_amb4.value, 0.0)) return false;

		double mw_value_1st = 0.0; double mw_sigma_1st = 0.0;
		double mw_value_2nd = 0.0; double mw_sigma_2nd = 0.0;
		double mw_value_3rd = 0.0; double mw_sigma_3rd = 0.0;
		double mw_value_4th = 0.0; double mw_sigma_4th = 0.0;
	
		//使用xml给的band计算freq
		//GOBS_LC              LC1 = GOBS_LC::LC_L1;
		//GOBS_LC              LC2 = (index_key_input.second == t_gamb::AMB_ID::ALL_L1_L2) ? GOBS_LC::LC_L2 : GOBS_LC::LC_L3;
		//double length_L1 = CLIGHT / (t_gamb::LC2FREQ(amb.sat_1st, LC1));
		//double length_Lx = CLIGHT / (t_gamb::LC2FREQ(amb.sat_1st, LC2));
		//double length_WL = CLIGHT / (t_gamb::LC2FREQ(amb.sat_1st, LC1) - t_gamb::LC2FREQ(amb.sat_1st, LC2));
		//
		t_gobsgnss gnss;
		gnss.sat(amb.sat_1st);
		auto gsys = t_gsys::sat2gsys(amb.sat_1st);
		FREQ_SEQ freq1 = FREQ_1;
		FREQ_SEQ freq2 = (index_key_input.second == t_gamb::AMB_ID::ALL_L1_L2) ? FREQ_2 : FREQ_3;
		double length_L1 = CLIGHT / (gnss.frequency(_band_index[gsys][freq1]));
		double length_Lx = CLIGHT / (gnss.frequency(_band_index[gsys][freq2]));
		double length_WL = CLIGHT / (gnss.frequency(_band_index[gsys][freq1]) - gnss.frequency(_band_index[gsys][freq2]));


		mw_value_1st = L1_amb1.value/length_L1 - LX_amb1.value/length_Lx;
		mw_value_2nd = L1_amb2.value/length_L1 - LX_amb2.value/length_Lx;
		mw_value_3rd = L1_amb3.value/length_L1 - LX_amb3.value/length_Lx;
		mw_value_4th = L1_amb4.value/length_L1 - LX_amb4.value/length_Lx;

		if (double_eq(mw_value_1st, 0.0)) return false;
		if (double_eq(mw_value_2nd, 0.0)) return false;
		if (double_eq(mw_value_3rd, 0.0)) return false;
		if (double_eq(mw_value_4th, 0.0)) return false;

		amb.MW_value = mw_value_1st -
			mw_value_2nd -
			mw_value_3rd +
			mw_value_4th;

		amb.MW_sigma = sqrt(mw_sigma_1st * mw_sigma_1st +
			mw_sigma_2nd * mw_sigma_2nd +
			mw_sigma_3rd * mw_sigma_3rd +
			mw_sigma_4th * mw_sigma_4th);

		amb.Mw_value[1] = mw_value_1st;
		amb.Mw_value[2] = mw_value_2nd;
		amb.Mw_value[3] = mw_value_3rd;
		amb.Mw_value[4] = mw_value_4th;
		amb.Mw_sigma[1] = mw_sigma_1st;
		amb.Mw_sigma[2] = mw_sigma_2nd;
		amb.Mw_sigma[3] = mw_sigma_3rd;
		amb.Mw_sigma[4] = mw_sigma_4th;


		if (fabs(amb.MW_value - round(amb.MW_value) < 0.15))
		{
			amb.MW_fixed = (fabs(round(amb.MW_value) - amb.MW_value) < 0.15) ? true : false;
		}
		return true;
	}

	bool t_gallambs::_get_IF_Dd_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& amb)
	{
		const auto& rec_1st = bl.first;
		const auto& rec_2nd = bl.second;

		int index1 = -1;
		int index2 = -1;
		int index3 = -1;
		int index4 = -1;

		if (index_key.second == t_gamb::AMB_ID::ALL_L1 ||
			index_key.second == t_gamb::AMB_ID::IF_L1_L2 ||
			index_key.second == t_gamb::AMB_ID::IF_L1_L3)
		{
			index1 = (amb.idx2ow_1st.find(1) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(1)->second) : -1;
			index2 = (amb.idx2ow_1st.find(2) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(2)->second) : -1;
			index3 = (amb.idx2ow_1st.find(3) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(3)->second) : -1;
			index4 = (amb.idx2ow_1st.find(4) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(4)->second) : -1;
		}
		else
		{
			index1 = (amb.idx2ow_2nd.find(1) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(1)->second) : -1;
			index2 = (amb.idx2ow_2nd.find(2) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(2)->second) : -1;
			index3 = (amb.idx2ow_2nd.find(3) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(3)->second) : -1;
			index4 = (amb.idx2ow_2nd.find(4) != amb.idx2ow_2nd.end()) ? (amb.idx2ow_2nd.find(4)->second) : -1;
		}

		if (index1 < 0 || index2 < 0 || index3 < 0 || index4 < 0) return false;

		const auto& amb1 = this->_nOw_ambs[index_key][rec_1st][index1];
		const auto& amb2 = this->_nOw_ambs[index_key][rec_1st][index2];
		const auto& amb3 = this->_nOw_ambs[index_key][rec_2nd][index3];
		const auto& amb4 = this->_nOw_ambs[index_key][rec_2nd][index4];

		if (double_eq(amb1.value, 0.0)) return false;
		if (double_eq(amb2.value, 0.0)) return false;
		if (double_eq(amb3.value, 0.0)) return false;
		if (double_eq(amb4.value, 0.0)) return false;

		auto ID1 = this->_nOw_ambs[index_key][rec_1st][index1].ID;
		auto ID2 = this->_nOw_ambs[index_key][rec_1st][index2].ID;
		auto ID3 = this->_nOw_ambs[index_key][rec_2nd][index3].ID;
		auto ID4 = this->_nOw_ambs[index_key][rec_2nd][index4].ID;

		if (ID1 != index_key.second)  return false;
		if (ID1 != ID2) return false;
		if (ID1 != ID3) return false;
		if (ID1 != ID4) return false;

		amb.Ow_value[1] = amb1.value;
		amb.Ow_value[2] = amb2.value;
		amb.Ow_value[3] = amb3.value;
		amb.Ow_value[4] = amb4.value;

		amb.Ow_time[1] = make_pair(amb1.beg, amb1.end);
		amb.Ow_time[2] = make_pair(amb2.beg, amb2.end);
		amb.Ow_time[3] = make_pair(amb3.beg, amb3.end);
		amb.Ow_time[4] = make_pair(amb4.beg, amb4.end);

		if (index_key.second == t_gamb::AMB_ID::IF_L1_L2 || index_key.second == t_gamb::AMB_ID::IF_L1_L3)
		{
			amb.value = amb1.value - amb2.value - amb3.value + amb4.value;
			amb.sigma = 0.0;
		}
		else if (index_key.second == t_gamb::AMB_ID::ALL_L1 ||
			     index_key.second == t_gamb::AMB_ID::ALL_L2 ||
			     index_key.second == t_gamb::AMB_ID::ALL_L3)
		{
			amb.DD_LX_value[index_key.second] = amb1.value - amb2.value - amb3.value + amb4.value;
			amb.DD_LX_sigma[index_key.second] = 0.0;
		}
		else
		{
			amb.value = 0.0;
			amb.sigma = 0.0;
			amb.DD_LX_value[index_key.second] = 0.0;
			amb.DD_LX_sigma[index_key.second] = 0.0;
			amb.DD_LX_fixed[index_key.second] = false;
			return false;
		}


		return true;
	}

	bool t_gallambs::_get_1X_Dd_value(const pair<GSYS, t_gamb::AMB_ID>& index_key_input, const pair<string, string>& bl, t_gamb_Dd& amb)
	{
		pair<GSYS, t_gamb::AMB_ID> index_ID1;
		pair<GSYS, t_gamb::AMB_ID> index_ID2;

		if (index_key_input.second == t_gamb::AMB_ID::ALL_L1_L2)
		{
			index_ID1 = make_pair(index_key_input.first, t_gamb::AMB_ID::ALL_L1);
			index_ID2 = make_pair(index_key_input.first, t_gamb::AMB_ID::ALL_L2);
		}
		else if (index_key_input.second == t_gamb::AMB_ID::ALL_L1_L3)
		{
			index_ID1 = make_pair(index_key_input.first, t_gamb::AMB_ID::ALL_L1);
			index_ID2 = make_pair(index_key_input.first, t_gamb::AMB_ID::ALL_L3);

		}
		else
		{
			return false;
		}

		if (!_get_IF_Dd_value(index_ID1, bl, amb)) return false; 
		if (!_get_IF_Dd_value(index_ID2, bl, amb)) return false; 

		return true;
	}

	bool t_gallambs::_get_LX_Dd_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& amb)
	{
		const auto& rec_1st = bl.first;
		const auto& rec_2nd = bl.second;

		int index1 = -1;
		int index2 = -1;
		int index3 = -1;
		int index4 = -1;


		index1 = (amb.idx2ow_1st.find(1) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(1)->second) : -1;
		index2 = (amb.idx2ow_1st.find(2) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(2)->second) : -1;
		index3 = (amb.idx2ow_1st.find(3) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(3)->second) : -1;
		index4 = (amb.idx2ow_1st.find(4) != amb.idx2ow_1st.end()) ? (amb.idx2ow_1st.find(4)->second) : -1;


		if (index1 < 0 || index2 < 0 || index3 < 0 || index4 < 0) return false;

		const auto& amb1 = this->_nOw_ambs[index_key][rec_1st][index1];
		const auto& amb2 = this->_nOw_ambs[index_key][rec_1st][index2];
		const auto& amb3 = this->_nOw_ambs[index_key][rec_2nd][index3];
		const auto& amb4 = this->_nOw_ambs[index_key][rec_2nd][index4];

		if (double_eq(amb1.value, 0.0)) return false;
		if (double_eq(amb2.value, 0.0)) return false;
		if (double_eq(amb3.value, 0.0)) return false;
		if (double_eq(amb4.value, 0.0)) return false;

		auto ID1 = this->_nOw_ambs[index_key][rec_1st][index1].ID;
		auto ID2 = this->_nOw_ambs[index_key][rec_1st][index2].ID;
		auto ID3 = this->_nOw_ambs[index_key][rec_2nd][index3].ID;
		auto ID4 = this->_nOw_ambs[index_key][rec_2nd][index4].ID;

		if (ID1 != index_key.second)  return false;
		if (ID1 != ID2) return false;
		if (ID1 != ID3) return false;
		if (ID1 != ID4) return false;

		amb.Ow_value[1] = amb1.value;
		amb.Ow_value[2] = amb2.value;
		amb.Ow_value[3] = amb3.value;
		amb.Ow_value[4] = amb4.value;

		amb.Ow_time[1] = make_pair(amb1.beg, amb1.end);
		amb.Ow_time[2] = make_pair(amb2.beg, amb2.end);
		amb.Ow_time[3] = make_pair(amb3.beg, amb3.end);
		amb.Ow_time[4] = make_pair(amb4.beg, amb4.end);
		
		double wavelength = 0.0;
		//lijie
		t_gobsgnss gnss;
		gnss.sat(amb.sat_1st);
		auto gsys = t_gsys::sat2gsys(amb.sat_1st);
		switch (index_key.second)
		{
		case t_gamb::AMB_ID::ALL_L1:
			wavelength = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_1]);
			break;
		case t_gamb::AMB_ID::ALL_L2:
			wavelength = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_2]);
			break;
		case t_gamb::AMB_ID::ALL_L3:
			wavelength = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_3]);
			break;
		default:
			throw logic_error("error type of amb while get freq length.");
		}
		/*switch (index_key.second)
		{
		case t_gamb::AMB_ID::ALL_L1:
			wavelength = CLIGHT / t_gamb::LC2FREQ(amb.sat_1st, GOBS_LC::LC_L1);
			break;
		case t_gamb::AMB_ID::ALL_L2:
			wavelength = CLIGHT / t_gamb::LC2FREQ(amb.sat_1st, GOBS_LC::LC_L2);
			break;
		case t_gamb::AMB_ID::ALL_L3:
			wavelength = CLIGHT / t_gamb::LC2FREQ(amb.sat_1st, GOBS_LC::LC_L3);
			break;
		default:
			throw logic_error("error type of amb while get freq length.");
		}*/

		amb.value = amb1.value/ wavelength - amb2.value/ wavelength - amb3.value/ wavelength + amb4.value/ wavelength;
		amb.sigma = 0.05;


		string tmp = amb.sat_1st + "_" + amb.sat_2nd + "_" + amb.rec_1st + "_" + amb.rec_2nd;
		_glog->logInfo("t_gallambs", "_get_LX_Dd_value", "get Dd value for " + tmp);
		double dcycle = round(amb.value) - amb.value;
		if (fabs(dcycle) < 0.15)
		{
			amb.fixed = true;
			amb.MW_fixed = true;
			amb.WL_fixed = true;
			amb.NL_fixed = true;
			amb.MW_value = amb.value;
			amb.WL_value = 0.0;
			amb.NL_value = 0.0;
		}
		else
		{
			amb.fixed = false;
			amb.MW_value = false;
			amb.WL_value = false;
			amb.NL_fixed = false;
		}
	

		return true;
	}

	void t_gallambs::_get_MW(const string& rec, const string& sat, const int& freq, const set<t_gtime>& time_table, set<t_gtime>::iterator beg, set<t_gtime>::iterator end,
		double& value, double& sigma)
	{
		list<pair<double, double> > mwobs;

		for (auto now = beg; *now <= *end; now++)
		{
			if (now == time_table.end()) break;
			auto sys = t_gsys::sat2gsys(sat);
			auto key = make_pair(sat, *now);
			if (this->_MW_used[sys][rec].count(key) <= 0)
			{
				continue;
			}

			double mw_obs = this->_MW_used[sys][rec][key][freq];

			// 求出来的是0.0，白瞎
			if (double_eq(mw_obs, 0.0))
			{
				continue;
			}
			//cout << setw(25) << setprecision(14) << mw_obs << endl;

			mwobs.push_back(make_pair(mw_obs, 1.0));
		}

		// 所有采样点都没有观测值
		if (mwobs.empty())
		{
			value = 0.0;
			sigma = 0.0;
			return;
		};

		//平滑MW观测值
		double mean_mw = 0.0, _sigma = 0.0, mean_sigma = 0.0;
		getMeanWgt(true, mwobs, mean_mw, _sigma, mean_sigma);

		//平滑出0.0，白瞎
		if (double_eq(mean_sigma, 0.0)) mean_sigma = 999.0;

		value = mean_mw;
		sigma = mean_sigma;

	}

	set<pair<GSYS, t_gamb::AMB_ID>> t_gallambs::get_Dds_index_key(const t_gamb::AMB_ID& ID)
	{
		set<pair<GSYS, t_gamb::AMB_ID>> tmp;

		for (const auto& index_key : _nDd_ambs)
		{
			if (index_key.first.second == ID)
			{
				tmp.insert(index_key.first);
			}
		}
		return tmp;
	}

	set<t_gamb::AMB_ID> t_gallambs::IDs()
	{
		return _IDs;
	}


	bool sort_Dd(const t_gamb_Dd& A, const t_gamb_Dd& B)
	{
		int tmp_A = floor((A.end - A.beg) / 300.0);
		int tmp_B = floor((B.end - B.beg) / 300.0);

		if (tmp_A > tmp_B) return true;
		if (tmp_A < tmp_B) return false;

		if (A.sat_1st < B.sat_1st) return true;
		if (A.sat_1st > B.sat_1st) return false;

		if (A.sat_2nd < B.sat_2nd) return true;
		if (A.sat_2nd > B.sat_2nd) return false;

		if (A.rec_1st < B.rec_1st) return true;
		if (A.rec_1st > B.rec_1st) return false;

		if (A.rec_2nd < B.rec_2nd) return true;
		return false;
	}
	bool sort_Sd(const t_gamb_Sd& A, const t_gamb_Sd& B)
	{
		int tmp_A = floor((A.end - A.beg) / 10.0);
		int tmp_B = floor((B.end - B.beg) / 10.0);

		if (tmp_A > tmp_B) return true;
		if (tmp_A < tmp_B) return false;

		if (A.sgl_obj < B.sgl_obj) return true;

		if (A.dbl_obj_1st < B.dbl_obj_1st) return true;
		if (A.dbl_obj_1st > B.dbl_obj_1st) return false;

		if (A.dbl_obj_2nd < B.dbl_obj_2nd) return true;
		return false;
	}

	bool t_gallambs::_add_upd(string type, t_gamb_Sd& sd)
	{
		string func_id = "_add_upd";
		string sat1;
		string sat2;
		t_gtime beg;
		t_gtime end;
		double sd_rnl_cor;
		double sigcor;
		sat1 = sd.dbl_obj_1st;
		sat2 = sd.dbl_obj_2nd;
		beg = sd.beg;
		end = sd.end;
		double upd_wl1 = 0.0, upd_wl2 = 0.0, upd_nl1 = 0.0, upd_nl2 = 0.0, sig = 0.0, rnl1 = 0.0, rnl2 = 0.0;
		upd_nl1 = 0.0, upd_nl2 = 0.0, sig = 0.0;
		if (type == "WL")
		{
			if (!_getSingleUpd("WL", WL_IDENTIFY, sat1, upd_wl1, sig,rnl1)
				|| !_getSingleUpd("WL", WL_IDENTIFY, sat2, upd_wl2, sig,rnl2))
			{
				write_log_info(_glog, 1, "Warning[t_gallambs::_add_upd] : _getSingleUpd Wrong : WL, Sat: " + sat1 + " " + sat2);
				return false;
			}
			sd.MW_value += upd_wl2 - upd_wl1;
		}
		else
		{
			if (!_getSingleUpd("NL", end, sat1, upd_nl1, sig,rnl1)
				|| !_getSingleUpd("NL", end, sat2, upd_nl2, sig,rnl2))
			{
				write_log_info(_glog, 1, "Warning[t_gambiguity::_applyUpd] : _getSingleUpd Wrong : NL, Sat: " + sat1 + " " + sat2);
				return false;
			}
			
			sd.rnl += rnl2 - rnl1;
			sd.NL_value += sd.rnl;
		}
		//cout << "correct UPD" << type << sat1 << " " << sat2 << beg.str_mjdsod() << end.str_mjdsod() << "  " << upd_wl1 << upd_wl2 << setw(20) << setprecision(15) << sd.WL_value << endl;
		return true;
	}
		bool t_gallambs::_getSingleUpd(string mode, t_gtime t, string sat, double& value, double& sigma)
		{
			UPDTYPE upd_type = UPDTYPE::NONE;
			if (mode == "WL")     upd_type = UPDTYPE::WL;
			else if(mode == "NL") upd_type = UPDTYPE::NL;

			one_epoch_upd epoch_upd = _gupd->get_epo_upd(upd_type, t);
			value = 0.0;
			if (epoch_upd.find(sat) == epoch_upd.end()) return false;
			if (mode == "WL")
			{
				if (epoch_upd[sat]->sigma > 0.2 || epoch_upd[sat]->npoint <= 2) return false;
				value = epoch_upd[sat]->value;
			}
			else if (mode == "NL")
			{
				if (epoch_upd[sat]->sigma > 0.1 || epoch_upd[sat]->npoint <= 3) return false;
				value = epoch_upd[sat]->value;
				sigma += pow(epoch_upd[sat]->sigma, 2);
			}
			else
			{
				write_log_info(_glog, 1, "Warning[t_gambiguity::_applySingleUpd] : Undefined UPD Mode (WL/NL) : " + mode);
				return false;
			}
			return true;
		}


		bool t_gallambs::_getSingleUpd(string mode, t_gtime t, string sat, double& value, double& sigma, double& rnl)
		{
			UPDTYPE upd_type = UPDTYPE::NONE;
			if (mode == "WL")     upd_type = UPDTYPE::WL;
			else if (mode == "NL") upd_type = UPDTYPE::NL;

			rnl = 0;
			one_epoch_upd epoch_upd = _gupd->get_epo_upd(upd_type, t);
			value = 0.0;
			if (epoch_upd.find(sat) == epoch_upd.end()) return false;
			if (mode == "WL")
			{
				if (epoch_upd[sat]->sigma > 0.2 || epoch_upd[sat]->npoint <= 2) return false;
				value = epoch_upd[sat]->value;
			}
			else if (mode == "NL")
			{
				if (epoch_upd[sat]->sigma > 0.1 || epoch_upd[sat]->npoint <= 3) return false;
				value = epoch_upd[sat]->value;
				sigma += pow(epoch_upd[sat]->sigma, 2);
				rnl += epoch_upd[sat]->value;
			}
			else
			{
				write_log_info(_glog, 1, "Warning[t_gambiguity::_applySingleUpd] : Undefined UPD Mode (WL/NL) : " + mode);
				return false;
			}
			return true;
		}

		void t_gallambs::_setWaveLen(map<string, map<string, double>>& sys_wavelen)
		{
			t_gobsgnss* gnss = new t_gobsgnss();
			set<string> sys = dynamic_cast<t_gsetgen*>(_gset)->sys();
			for (auto sys_iter = sys.begin(); sys_iter != sys.end(); sys_iter++)
			{
				GOBSBAND b1, b2, b3;
				string sys_short = "NONE";
				switch (t_gsys::str2gsys(*sys_iter))
				{
				case GPS:
					b1 = BAND_1;
					b2 = BAND_2;
					b3 = BAND_5;
					sys_short = "G";
					break;
				case GAL:
					b1 = BAND_1;
					b2 = BAND_5;
					b3 = BAND_7;
					sys_short = "E";
					break;
				case BDS:
					b1 = BAND_2;
					b2 = BAND_7;
					b3 = BAND_6;
					sys_short = "C";
					break;
				case GLO:
					b1 = BAND_1;
					b2 = BAND_2;
					b3 = BAND;
					sys_short = "R";   // FDMA
					break;
				default:
					break;
				}
				if (sys_short == "NONE") return;
				if (*sys_iter == "GLO")
				{
				}
				else
				{
					gnss->sat(sys_short);
					sys_wavelen[sys_short]["L1"] = gnss->wavelength(b1);
					sys_wavelen[sys_short]["L2"] = gnss->wavelength(b2);
					sys_wavelen[sys_short]["L3"] = gnss->wavelength(b3);
					sys_wavelen[sys_short]["WL"] = gnss->wavelength_WL(b1, b2);
					sys_wavelen[sys_short]["NL"] = gnss->wavelength_NL(b1, b2);
					sys_wavelen[sys_short]["EWL"] = gnss->wavelength_WL(b2, b3);
				}
			}
			delete gnss; gnss = nullptr;
		}
}


