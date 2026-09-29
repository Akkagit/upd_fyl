/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambBlRec.h
* @brief	site baseline of SD and DD ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#ifndef GAMBBLREC_H
#define GAMBBLREC_H 

#include "gutils/ExportAMB.h"
#include "gio/glog.h"
#include "gutils/gnss.h"
#include "gambfix/gambOw.h"
#include "gambfix/gambBlSat.h"
#include "gambfix/gambBl.h"
#include "gall/gallobj.h"
#include "gall/gallobs.h"

using namespace std;
using namespace gnut;

namespace great
{
	// 双差模糊度
	class AMB_LIBRARY_EXPORT t_gamb_baseline_recs : public t_gamb_baseline_sats
	{
	public:
		// 默认构造函数
		t_gamb_baseline_recs();
		t_gamb_baseline_recs(const t_gtime& time, double cut_length);

		// 直接给定所有测站构建出独立的基线
		virtual ~t_gamb_baseline_recs();

		void add_cut(const t_gtime& time, double cut_length);
		void add_gobj(t_gallobj* obj);
		void add_gset(t_gsetbase* set);
		void add_gobs(t_gallobs* obs);
		void add_glog(t_glog* log);
		void add_rec(const string& rec, const string& sat);

		vector<pair<string, string>> baselines();

	protected:
		t_gsetbase*             _gset = nullptr;
		t_gallobj*              _gobj = nullptr;
		t_glog*                 _glog = nullptr;

		t_gtime                 _cut_time;
		t_gtime                 _beg_time;
		t_gtime                 _end_time;
		double                  _cut_length = 3500; //km

		vector<string>                         _recs;
		map<string,  set<string>>              _sys;
		map<string,  t_gtriple>                _xyz;
	protected:
		bool  _get_rec_xyz();
		bool  _create_baseline();
		bool  _check_amb_depend();
	};

	bool  compare(const t_gamb_baseline& bl1, const t_gamb_baseline& bl2);
}

#endif