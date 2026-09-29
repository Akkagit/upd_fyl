/**
*
* @verbatim
History
	-1.0 Jie Li  2019-11-14  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambcon.cpp
* @brief	Storage the ambigulity constraint files' data
*
* @author   Jie Li, Wuhan University
* @version	1.0.0
* @date		2019-11-14
*
*/

#include "gdata/gambcon.h"
 
using namespace std;
namespace great {

	t_gambcon::t_gambcon()
	{
		gtrace("t_gambcon::constructor");
		id_type(t_gdata::AMBCON);
	}

	t_gambcon::~t_gambcon()
	{
	}

	void t_gambcon::set_header(t_gtime beg, double intv, int num_amb, set<string>& sats, set<string>& recs)
	{
		_beg = beg;
		_interval = intv;
		_num_fixamb = num_amb;
		_sats = sats;
		_sites = recs;
	}

	void t_gambcon::add_sdamb(t_gamb_Sd& sdamb)
	{
		_sdamb.push_back(sdamb);
		_nsdamb[sdamb.beg].push_back(sdamb);
	}
	void t_gambcon::add_ddamb(t_gamb_Dd& ddamb)
	{
		_ddamb.push_back(ddamb);
		_nddamb[ddamb.beg].push_back(ddamb);
	}

	vector<t_gamb_Dd> t_gambcon::get_Dds()
	{
		return this->_ddamb;
	}

	/*vector<t_gamb_Dd> t_gambcon::get_Dds(t_gtime time)
	{
		vector<t_gamb_Dd> ans;
		for (auto iter : _index_ddamb[time])
		{
			ans.push_back(*iter);
		}
		return ans;
	}
*/
	void t_gambcon::get_Dds(vector<t_gamb_Dd>& ddamb)
	{
		ddamb = _ddamb;
		std::sort(ddamb.begin(), ddamb.end(), ambcon_sort_Dd);
	}

	void t_gambcon::get_Sds(vector<t_gamb_Sd>& sdamb)
	{
		sdamb = _sdamb;
	}
	vector<t_gamb_Dd>& t_gambcon::get_Dds(const t_gtime& time)
	{
		// TODO: 在此处插入 return 语句
		return _nddamb[time];
	}

	vector<t_gamb_Sd>& t_gambcon::get_Sds(const t_gtime& time)
	{
		return _nsdamb[time];
	}

	void t_gambcon::get_header(t_gtime & beg, double & intv, int & num_amb, set<string>& sats, set<string>& recs)
	{
		beg = _beg;
		intv = _interval;
		num_amb = _num_fixamb;
		sats = _sats;
		recs = _sites;
	}

	bool ambcon_sort_Dd(const t_gamb_Dd& A, const t_gamb_Dd& B)
	{
		if (A.ID > B.ID) return true;
		if (A.ID < B.ID) return false;

		int tmp_A = floor((A.end - A.beg) / 30.0);
		int tmp_B = floor((B.end - B.beg) / 30.0);

		if (tmp_A > tmp_B) return true;
		if (tmp_A < tmp_B) return false;


		if (A.sat_1st < B.sat_1st) return true;
		if (A.sat_1st > B.sat_1st) return false;

		if (A.sat_2nd < B.sat_2nd) return true;
		if (A.sat_2nd > B.sat_2nd) return false;


		if (A.rec_1st < B.rec_1st) return true;
		if (A.rec_1st > B.rec_1st) return false;

					
		if (A.rec_2nd < B.rec_2nd) return true;
		if (A.rec_2nd > B.rec_2nd) return false;

		return false;
	}

	bool ambcon_sort_Sd(const t_gamb_Sd& A, const t_gamb_Sd& B)
	{
		if (A.ID > B.ID) return true;
		if (A.ID < B.ID) return false;

		int tmp_A = floor((A.end - A.beg) / 30.0);
		int tmp_B = floor((B.end - B.beg) / 30.0);
		if (tmp_A > tmp_B) return true;
		if (tmp_A < tmp_B) return false;
		if (A.dbl_obj_1st < B.dbl_obj_1st) return true;
		if (A.dbl_obj_1st > B.dbl_obj_1st) return false;
		if (A.dbl_obj_2nd < B.dbl_obj_2nd) return true;
		if (A.dbl_obj_2nd > B.dbl_obj_2nd) return false;
		if (A.sgl_obj < B.sgl_obj) return true;
		if (A.sgl_obj > B.sgl_obj) return false;
		return false;
	}
}//namespace

