/**
* @verbatim
	History
	 -1.0 Guolong Feng	 2020-07-20  creat the file.
	 -1.1 Jie Li         2020-09-14  add new way for clock repairing 
  @endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file			gclockrepair.h
* @author       Guolong Feng, Wuhan University
* @version		1.0.0
* @date			2020-07-20
*
*/

#ifndef CLOCKREPAIR_H
#define CLOCKREPAIR_H

#include "gio/gfile.h"

#include "gall/gallobs.h"

#include "gcoders/rinexo.h"

#include "gset/gsetout.h"
#include "gset/gsetgen.h"
#include "gset/gsetgnss.h"

#include "gutils/gfileconv.h"
#include "gutils/ExportBase.h"
using namespace gnut;

namespace great
{
	// 钟差修复
	class BASE_LIBRARY_EXPORT t_gclock_repair
	{
	public:
		t_gclock_repair() {};
		t_gclock_repair(t_gsetbase* settings, t_gallobs* allobs);

		virtual ~t_gclock_repair() {};
		//isRealtime -- true->李盼博士博士论文提到的方法，边修改边探测（还存在一点问题）； false-> pannda,采用源数据探测，已调试通过
		void repairClockJump(string site, const t_gtime& begT, const t_gtime& endT, double sampling, bool isRealtime = false);

		int ProcessBatch(string site, const t_gtime& begT, const t_gtime& endT, double sampling);

		void setLog(t_glog* log) { _glog = log; };
		void setObj(t_gallobj* obj) { _allobj = obj; };
		void setOBS(t_gallobs* obs) { _allobs = obs; };

	private:

		void _repairClockJump(vector<t_spt_gobs> AllObsPre, vector<t_spt_gobs> AllObsCrt);
		void _repairClockJump(const vector<t_gsatdata>& AllOrigPre, vector<t_spt_gobs> AllObsCrt);
		int _clock_jump_detection(vector<t_spt_gobs> AllObsPre, vector<t_spt_gobs> AllObsCrt);
		int _clock_jump_repair(int numJump, vector<t_spt_gobs> AllObs);

		int _check_msec_jump(const vector<t_gsatdata>& epoOrigPre, vector<t_spt_gobs> AllObsCrt,bool isPhase);


		vector<shared_ptr<t_gobsgnss>>   _epoDataPre;
		vector<t_gsatdata>				 _epoOrigPre;//当前历元原始观测值
		map<GSYS, map<FREQ_SEQ, GOBSBAND>> _band_index;
		double _sumS;
		int    _jumpCor;

		t_gallobs*   _allobs;
		t_gallobj*   _allobj;
		t_gsetbase*  _set;
		t_glog*      _glog;
		set<string>  _site_list;

	};



} // namespace

#endif
