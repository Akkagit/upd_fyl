/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambBlSat.h
* @brief	sattelite baseline of SD and DD ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#ifndef GAMBBLSAT_H
#define GAMBBLSAT_H 

#include "gutils/ExportAMB.h"
#include "gio/glog.h"
#include "gutils/gnss.h"
#include "gambfix/gambBl.h"

using namespace std;
using namespace gnut;

namespace great
{
	// 双差模糊度
	class AMB_LIBRARY_EXPORT t_gamb_baseline_sats
	{
	public:
		// 默认构造函数
		t_gamb_baseline_sats();
		// 直接给定所有测站构建出独立的基线
		t_gamb_baseline_sats(set<string>& sats);
		virtual ~t_gamb_baseline_sats();

		// 独立基线的条数
		void add_sat(const string& sat);
		int  baseline_num();
		vector<pair<string, string>>  baselines();

		vector<tuple<string, string, string, string>>  get_depend_Dd_baseline(vector<tuple<string, string, string, string>>& new_baseline);
	protected:
		vector<t_gamb_baseline>  _baselines;
		vector<string>           _sats;
	};
}

#endif