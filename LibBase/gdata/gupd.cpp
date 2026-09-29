/**
*
* @verbatim
History
-1.0 BoWong  2019-02-23  creat the file.
-1.1 BoWong  2019-04-08  Adding Doxygen Style Code Remarks
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gupd.cpp
* @brief	Storage the upd files' data include wide-lane, narrow-lane and ex-wide-lane
*
* @author   BoWong, Wuhan University
* @version	1.0.0
* @date		2019-04-08
*
*/

#include"gdata/gupd.h"

using namespace std;
namespace great {

	/** @brief default constructor. */
	t_updrec::t_updrec()
	{
		obj = " ";
		npoint = 0;
		value = 0.0;
		sigma = 1E4;
		isRef = false;
	}

	/** @brief default constructor. */
	t_gupd::t_gupd(){
		gtrace("t_gupdlane::constructor");
		id_type(t_gdata::UPD);

		_ewl_flag = EWL_IDENTIFY;
		_ewl24_flag = EWL24_IDENTIFY;
		_ewl25_flag = EWL25_IDENTIFY;
		_wl_flag = WL_IDENTIFY;

		_valid_beg[UPDTYPE::EWL_epoch] = FIRST_TIME;
		_valid_beg[UPDTYPE::EWL]       = FIRST_TIME;
		_valid_beg[UPDTYPE::EWL24]     = FIRST_TIME;
		_valid_beg[UPDTYPE::EWL25]     = FIRST_TIME;
		_valid_beg[UPDTYPE::WL]        = FIRST_TIME;
		_valid_beg[UPDTYPE::NL]        = FIRST_TIME;
		_valid_beg[UPDTYPE::IFCB]      = FIRST_TIME;

		_wait_stream = true;
	}

	/** @brief default destructor. */
	t_gupd::~t_gupd()
	{
		_upd.clear();
	}

	void t_gupd::add_sat_upd(UPDTYPE upd_type, t_gtime epoch, string prn, t_updrec one_sat_upd)
	{
		_upd[upd_type][epoch][prn] = make_shared<t_updrec>(one_sat_upd);
	}

	void t_gupd::add_epo_upd(UPDTYPE upd_type, t_gtime epoch, one_epoch_upd one_epo_upd){
		_upd[upd_type][epoch] = one_epo_upd;
	}

	void t_gupd::set_est_updtype(UPDTYPE mode)
	{
		_est_upd_type = mode;
	}

	UPDTYPE t_gupd::get_est_updtype()
	{
		return _est_upd_type;
	}

	bool t_gupd::upd_usable(const t_gtime& t, const string& str)
	{
		try
		{
			bool upd_usable;
			if (_est_upd_type == UPDTYPE::NL) {
				upd_usable = (_upd[_est_upd_type][t][str]->sigma < 0.100 && _upd[_est_upd_type][t][str]->npoint >= 3)
					|| _upd[_est_upd_type][t][str]->isRef || _upd[_est_upd_type][t][str]->sigma <= 1e-4;
			}
			else
			{
				upd_usable = (_upd[_est_upd_type][t][str]->sigma < 0.120 && _upd[_est_upd_type][t][str]->npoint >= 2)
					|| _upd[_est_upd_type][t][str]->isRef || _upd[_est_upd_type][t][str]->sigma <= 5e-4;
			}

			return upd_usable;
		}
		catch (...)
		{
			cout << "ERROR:t_gupd::upd_usable wrong!!!" << endl;
			throw(-1);
		}
	}

	void t_gupd::re_init_upd(const UPDTYPE& upd_type, const t_gtime& t, string str)
	{
		t_updrec one_upd;
		one_upd.npoint = 0;
		one_upd.value = 0.0;
		one_upd.sigma = 1E4;
		one_upd.isRef = false;
		t_gupd::add_sat_upd(upd_type, t, str,  one_upd);
	}

	one_epoch_upd& t_gupd::get_epo_upd(const UPDTYPE& upd_type, const t_gtime& t)
	{
		if (_upd[upd_type].find(t) != _upd[upd_type].end())
		{
			return _upd[upd_type][t];
		}
		else
		{
			auto it_epo = _upd[upd_type].lower_bound(t);
			if (it_epo != _upd[upd_type].end())
			{
				if (it_epo->first.diff(t) >= 30) {
					_wait_stream = false;
					return _null_epoch_upd;				
				}
				else {
					return it_epo->second;
				}
			}
			else
			{
				return _null_epoch_upd;
				
			}

			// Second
			//auto it_epo = _upd[upd_type].lower_bound(t);
			//if (it_epo == _upd[upd_type].end())
			//{
			//	if (it_epo != _upd[upd_type].begin())
			//	{
			//		it_epo--;
			//		if (t.diff(it_epo->first) < 30) return it_epo->second;
			//	}
			//	return _null_epoch_upd;
			//}
			//else
			//{
			//	auto it_epo2 = it_epo;
			//	if (it_epo2 != _upd[upd_type].begin()) it_epo2--;
			//	double diff1 = it_epo->first.diff(t);
			//	double diff2 = t.diff(it_epo2->first);
			//	if (diff1 >= 30 && diff2 >= 30) {
			//		return _null_epoch_upd;
			//	}
			//	else {
			//		it_epo = diff1 <= diff2 ? it_epo : it_epo2;
			//		return it_epo->second;
			//	}			
			//}

			// First
			//auto it_epo = _upd[upd_type].begin();
			//for (; it_epo != _upd[upd_type].end(); it_epo++)
			//{
			//	if (it_epo->first == _ewl_flag || it_epo->first == _wl_flag) continue;
			//	if (it_epo->first > t && it_epo->first.diff(t) < 30) break;  //不超出30秒才可以使用
			//	if (it_epo->first < t && t.diff(it_epo->first) < 30) break;  //不超出30秒才可以使用
			//}
			//if (it_epo == _upd[upd_type].end())
			//{
			//	return _null_epoch_upd;
			//}
			//else 
			//{
			//	return it_epo->second;			
			//}
		}
	}

	void t_gupd::reset_upd(const UPDTYPE& upd_type, const t_gtime& t, const string& str, const double& value,
		const double& sigma, const int& npoint)
	{
		_upd[upd_type][t][str]->npoint = npoint;
		_upd[upd_type][t][str]->value = value;
		_upd[upd_type][t][str]->sigma = sigma;
	}

	void t_gupd::reset_upd_value(const UPDTYPE& upd_type, const t_gtime& t, const string& str, const double& value)
	{
		_upd[upd_type][t][str]->value = value;
	}

	void t_gupd::copy_upd(const UPDTYPE& upd_type, const t_gtime& pre_t, const string& str, const t_gtime& current_t,
		const bool& is_first, const bool& is_site)
	{
		if (is_first == true)
		{
			t_gupd::re_init_upd(upd_type, current_t, str);
		}
		else
		{
			t_updrec one_upd;
			one_upd.npoint = _upd[upd_type][pre_t][str]->npoint;
			one_upd.value = _upd[upd_type][pre_t][str]->value;
			one_upd.sigma = _upd[upd_type][pre_t][str]->sigma;
			one_upd.isRef = _upd[upd_type][pre_t][str]->isRef;
			if (is_site) _upd[upd_type][pre_t].erase(str);
			t_gupd::add_sat_upd(upd_type, current_t, str, one_upd);

		}
		return;
	}

	void t_gupd::set_valid_beg(const UPDTYPE& upd_type, const t_gtime& t)
	{ 
		if(_valid_beg[upd_type] > t) _valid_beg[upd_type] = t;
	}

	bool t_gupd::upd_avail(const t_gtime& now)
	{
		one_epoch_upd epoch_upd_wl = this->get_epo_upd(UPDTYPE::WL, now);
		one_epoch_upd epoch_upd_nl = this->get_epo_upd(UPDTYPE::NL, now);
		if ((epoch_upd_wl.size() < 5 || epoch_upd_nl.size() < 5) && _wait_stream) {
			return false;
		}
		else
		{
			_wait_stream = true;
			return true;
		}
		
	}

	string updmode2str(UPDTYPE mode)
	{
		string str;
		switch (mode)
		{
		case UPDTYPE::EWL:
			str = "EWL";
			break;
		case UPDTYPE::EWL24:
			str = "EWL24";
			break;
		case UPDTYPE::EWL25:
			str = "EWL25";
			break;
		case UPDTYPE::EWL_epoch:
			str = "EWL_epoch";
			break;
		case UPDTYPE::NL:
			str = "NL";
			break;
		case UPDTYPE::WL:
			str = "WL";
			break;
		case UPDTYPE::IFCB:
			str = "IFCB";
			break;
		case UPDTYPE::RT_SAVE:
			str = "RT_SAVE";
			break;
		case UPDTYPE::NONE:
			str = "NONE";
			break;
		default:
			str = "NONE";
			break;
		}
		return str;
	}

	UPDTYPE str2updmode(string mode)
	{
		if (mode == "EWL_epoch")  return UPDTYPE::EWL_epoch;
		else if (mode == "EWL")   return UPDTYPE::EWL;
		else if (mode == "EWL24")   return UPDTYPE::EWL24;
		else if (mode == "EWL25")   return UPDTYPE::EWL25;
		else if (mode == "WL")    return UPDTYPE::WL;
		else if (mode == "NL")    return UPDTYPE::NL;
		else if (mode == "IFCB")  return UPDTYPE::IFCB;
		else if (mode == "RT_SAVE") return UPDTYPE::RT_SAVE;
		else                      return UPDTYPE::NONE;
	}

}//namespace