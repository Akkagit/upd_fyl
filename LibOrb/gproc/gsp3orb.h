/**
*
* @verbatim
	History
	 -1.0 Yujie Qin		2019-10-10 creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gsp3orb.h
* @brief		convert sp3 file or kin file to orb file
*
* @author       Yujie Qin	, Wuhan University
* @version		1.0.0
* @date		    2019-10-10
*
*/


#ifndef GSP3ORB_H
#define GSP3ORB_H

#include "gutils/gtime.h"
#include "gutils/gmutex.h"
#include "gutils/gnss.h"
#include "gutils/gtrs2crs.h"
#include "gio/glog.h"
#include "gio/giof.h"
#include "gdata/gics_leo.h"
#include "gdata/gorb.h"
#include "gset/gsetinp.h"
#include "gset/gsetout.h"
#include "gall/gallprec.h"
#include "gall/gallics_leo.h"
#include "gall/gallics_nav.h"
#include "newmat/newmat.h"
#include "gset/gsetsp3orb.h"
#include "gset/gsetforce.h"
#include <vector>
#include "gutils/gleo.h"
#include "gio/gfile.h"
#include "gcoders/icsnav.h"
#include "gcoders/icsleo.h"
#include "gutils/ExportOI.h"
//#define NEXT_POINT_NUM 5

using namespace std;
using namespace gnut;

namespace great
{
	class OI_LIBRARY_EXPORT t_gsp3orb
	{

	public:
		t_gsp3orb();
		t_gsp3orb(set<string> mark, set<string> sys, t_gsetbase* set, t_glog * log);
		virtual ~t_gsp3orb();
	public:
		bool processBatch(const t_gtime& beg, const t_gtime& end);
		void setDataUsed(t_gallprec* gsp3);
		void setDataUsed(t_gpoleut1* poleut1);
		void GenerateProduct();
		t_gallics_leo*    _ics_leo;
		t_gallics_nav*    _ics_nav;
	protected:
		bool _initial_data();
		bool _get_settings();
		bool _extend_time_period(t_gtime beg, t_gtime end, int step_out, t_gtime& extend_beg, t_gtime& extend_end);
		void _generate_icsfile();
		void _generate_orbfile();
		map<string,map<string, double>>  _ref_state;
		set<string>       _sats;
		set<string>       _sys;
		t_glog*		      _glog;             ///< Genereal   log output. 
		t_gsetbase*		  _set;              ///< Base setting. 
		t_gorb*           _gorb;             ///< orb data (output)
		t_gtime           _epoch;          ///< current epoch
		t_gallprec*       _sp3;
		t_gpoleut1*       _poleut1;
		t_gtime           _beg;
		t_gtime           _end;
		t_gtime           _beg_extend;
		t_gtime           _end_extend;
		t_gtime           _ref;
		t_gtrs2crs*       _trs2crs;
		int               _interval;
		bool              _flag;
		string            _frame;
		map<string, int>               _integra_degree; ///< By zhangwei for multi-LEO processiing
		map<string, int>               _step_inp;
		map<string, int>               _step_out;
		map<string, int>               _neq_size;
		map<string, map<string, int>>  _drag_para;
		map<string, map<string, int>>  _srp_para;
		map<string, map<string, int>>  _cust_para;
		map<string, vector<string>>    _param;
		map<string, int>               _param_size;

	};
}
#endif // !GSP3ORB_H