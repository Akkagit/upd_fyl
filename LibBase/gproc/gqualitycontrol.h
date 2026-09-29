
/**
* @verbatim
	History
	 -1.0 Guolong Feng	 2020-07-05  creat the file.
  @endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file			qualityControl.h
* @author       Guolong Feng, Wuhan University
* @version		1.0.0
* @date			2020-07-05
*
*/

#ifndef QUALITYCONTROL_H
#define QUALITYCONTROL_H

#include "gall/gallnav.h"
#include "gall/gallobs.h"

#include "gdata/gsatdata.h"

#include "gset/gsetgen.h"
#include "gset/gsetrec.h"
#include "gset/gsetgnss.h"
#include "gset/gsetproc.h"
#include "gset/gsetturboedit.h"

#include "gutils/ExportBase.h"

namespace great
{
	enum SMOOTH_MODEL
	{
		SMT_DOPPLER,
		SMT_PHASE,
		SMT_NONE
	};

	// 伪距平滑
	class BASE_LIBRARY_EXPORT t_gsmooth
	{
	public:
		t_gsmooth() {};
		t_gsmooth(t_gsetbase* settings);

		virtual ~t_gsmooth() {};

		void smooth_range_obs(vector<t_gsatdata>& obsdata, const t_gtime& now);

	private:

		void _doppler_smt_range(vector<t_gsatdata>& obsdata, const t_gtime& now);
		void _phase_smt_range(vector<t_gsatdata>& obsdata, const t_gtime& now);

		SMOOTH_MODEL                    _smoothModel;
		int                             _smoothWindow;
		double                          _smoothFactor;
		double                          _sampling;
		// site/sat/range_obs/pre_smt_value
		map<string, map<string, map<GOBS, double> > >                   _pre_smt_range;
		// 用于判断滑动窗口
		map<string, map<string, t_gtime > >                             _smt_beg_time;
		// site/sat/波段/(obs_type,phase/dolpper value)
		map<string, map<string, map<GOBSBAND, pair<GOBS, double> > > >  _pre_orig_val;

	};

	// 北斗星上多路径改正
	// BeiDou satellite-induced code pseudorange variations correct
	class BASE_LIBRARY_EXPORT t_gbds_codebias_cor
	{
	public:
		t_gbds_codebias_cor(t_gsetbase* settings);
		virtual ~t_gbds_codebias_cor() {};

		void apply_IGSO_MEO(const string& rec, t_gtriple& rec_crd, t_gallnav* gnav, vector<t_gsatdata>& obsdata);

	private:

		t_gsetbase*                          _set;
		map< GSYS, map<FREQ_SEQ, GOBSBAND> > _band_index;
		bool                                 _correct_bds_code_bias;
		// Wanninger & Beer : BeiDou satellite-induced code pseudorange variations: diagnosis and therapy [unit:m]
		map<GOBSBAND, map<string, map<int, double> > > _IGSO_MEO_Corr;

		bool _recAprCoordinate(const string& rec, t_gtriple& rec_crd, t_gallnav* gnav, vector<t_gsatdata>& obsdata);

	};


	// 粗差探测与剔除
	class BASE_LIBRARY_EXPORT t_goutliers_process
	{
	public:
		t_goutliers_process(t_gsetbase* settings, t_glog* glog);

		virtual ~t_goutliers_process();

		void setLog(t_glog* glog) { _log = glog; };

		void excludeBadObs(vector<t_gsatdata>& obsdata); // Low SNR < 10dB

		void flagRangeOutliers(shared_ptr<t_gobsgnss> ObsPre, shared_ptr<t_gobsgnss> Obs, double sampling);

	private:
		t_glog*            _log;
		t_gsetbase*        _set;
		t_giof*            _debug_outliers=nullptr;

	};

	//// 钟差修复
	//class BASE_LIBRARY_EXPORT t_gclock_repair
	//{
	//public:
	//	t_gclock_repair() {};
	//	t_gclock_repair(t_gsetbase* settings);

	//	virtual ~t_gclock_repair() {};

	//	void repairClockJump(vector<t_spt_gobs> AllObsPre, vector<t_spt_gobs> AllObs);

	//private:

	//	int _clock_jump_detection(vector<t_spt_gobs> AllObsPre, vector<t_spt_gobs> AllObs);
	//	int _clock_jump_repair(int numJump, vector<t_spt_gobs> AllObs);

	//	map<GSYS, map<FREQ_SEQ, GOBSBAND>> _band_index;
	//	double _sumS;
	//};

	// 组合
	class BASE_LIBRARY_EXPORT t_gqualitycontrol
	{
	public:
		t_gqualitycontrol(t_gsetbase* settings, t_gallnav* gnav);
		virtual ~t_gqualitycontrol();

		int processOneEpoch(const t_gtime& now, const string& rec, t_gtriple& rec_crd, vector<t_gsatdata>& obsdata);
		void setNav(t_gallnav* gnav) { _gnav = gnav; };

	protected:

		t_glog*            _log;
		t_gsetbase*        _set;
		t_gallnav*         _gnav;

		t_gbds_codebias_cor   _bds_codebias_cor;
		t_gsmooth             _smooth_range;
		t_goutliers_process   _outliers_proc;


	};


} // namespace

#endif
