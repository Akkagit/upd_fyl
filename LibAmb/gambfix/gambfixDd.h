/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambfixDd.h
* @brief	for the use of DD fix ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*
*/

#ifndef GAMBFIXDD_H
#define GAMBFIXDD_H 

#include "gutils/ExportAMB.h"
#include "gambfix/gambfix.h"
#include "gset/gsetamb.h"
#include "gdata/gambcon.h"

using namespace std;

namespace great
{

	/**
	*@ brief class for LAMADA search finding fixing ambiguity.
	*/
	class AMB_LIBRARY_EXPORT t_gambfixDd : public t_gambfix
	{
	public:
		t_gambfixDd();
		t_gambfixDd(t_glog* log, t_gsetbase* set);
		
		virtual ~t_gambfixDd();

		
		// 确定beg到end时间内的双差模糊度
		bool ProcessBatch();
		bool ProcessBatch_ALL();

		// 确定测站基线
		vector<pair<string, string>> baselines();

	protected:
		t_gamb_baseline_recs  _gRecBl;
		t_gamb_baseline_sats  _gSatBl;

		// 通过测站/卫星列表确定地面测站的基线
		bool _screen_rec_and_sat(t_gallpar* pars, t_gsetbase* set, t_gallobj* obj);
		void _write_ambcon();
	private:
		string _class_id = "t_gambfixDd";
	};
}

#endif