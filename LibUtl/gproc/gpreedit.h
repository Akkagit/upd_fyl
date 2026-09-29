/**
*
* @verbatim
	History
	 -1.0 Jiande Huang  2019-09-20  created
	 -1.1 Jiande Huang  2020-09-12  Fixed memory leak
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gpreedit.h
* @brief		main control class for estimate clk	
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-09-20
*

* 功能设计：为了读取配置文件，准备lsq步骤前的所有文件
* 功能需求：
*	1. 生成ics文件
*	2. 生成tb 文件
*   3. 检查必要的数据文件是否存在
*
*/


#ifndef GPREEDIT_H
#define GPREEDIT_H

#include "gutils/ExportUtl.h"
#include "gset/gsetbase.h"
#include "gset/gsetinp.h"
#include "gio/glog.h"
#include "gall/gallnav.h"
#include "gall/gallics_nav.h"
#include "gutils/gtrs2crs.h"
#include "gproc/gsppflt.h"

#include <map>

namespace great
{

	/**
	* @brief main control class for caculate clk
	*/
	class UTL_LIBRARY_EXPORT t_gpreedit
	{
	public:
		
		t_gpreedit();
		t_gpreedit(t_gsetbase* set, t_glog* log);

		bool ProcessBatch(const t_gtime& beg, const t_gtime& end);
		
		virtual ~t_gpreedit();

	protected:
		t_gsetbase* _gset;
		t_glog*     _log;
		set<string> _rec_list;
		set<string> _sat_list;
		set<string> _sys_list;
		t_gtime		_beg_time;
		t_gtime		_end_time;
		double      _intv;

		t_gallnav* _gnav = nullptr;
		t_gallobj* _gobj  = nullptr;
		t_gallobs* _gobs  = nullptr;
		set<string> _inp_all;

		t_gallics_nav _gics;

		set<ICS_PAR> _def_ECOM_3;
		set<ICS_PAR> _def_ECOM_5;
		set<ICS_PAR> _def_ECOM_9;
		set<ICS_PAR> _def_ECOM2;
		FORCE_MODEL  _def_solar_model;
		t_gtrs2crs*  _trs2crs;

		map<string, t_gsppflt*>   _sppflt;		    ///< sppflt, provide init crd


	public:
		void set_nav(t_gallnav* nav);
		void set_erp(t_gpoleut1* erp);
		void set_obj(t_gallobj* obj);
		void set_obs(t_gallobs* obs);
	protected:
		bool _gnerate_ics_file(const t_gtime& beg, const t_gtime& end, double intv);
		bool _gnerate_clk_file(const t_gtime& beg, const t_gtime& end, double intv);
		bool _gnerate_tb_csh_file(const t_gtime& beg, const t_gtime& end, double intv);
		bool _gnerate_rec_crd_xml(const t_gtime& beg, const t_gtime& end);
		bool _check_inp_files(set<string>& _inp_all);

	};
	
}

#endif /*GCLKLSQ_H*/