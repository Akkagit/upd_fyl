/**
*
* @verbatim
	History
	 -1.0 Yujie Qin		2019-12-02 creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gorbsp3.h
* @brief		convert orb file to sp3 file
*
* @author       Yujie Qin	, Wuhan University
* @version		1.0.0
* @date		    2019-12-02
*
*/


#ifndef GORBSP3_H
#define GORBSP3_H

#include "gutils/gtime.h"
#include "gutils/gmutex.h"
#include "gutils/gnss.h"
#include "gutils/gtrs2crs.h"
#include "gio/glog.h"
#include "gio/giof.h"
#include "gdata/gorb.h"
#include "gdata/gleapsecond.h"
#include "gset/gsetinp.h"
#include "gset/gsetout.h"
#include "gall/gallprec.h"
#include "gall/gallics_leo.h"
#include "gall/gallics_nav.h"
#include "newmat/newmat.h"
#include "gset/gsetorbsp3.h"
#include <vector>
#include "gio/gfile.h"
#include "gutils/ExportOI.h"
#include "gcoders/sp3.h"
//#define NEXT_POINT_NUM 5

using namespace std;
using namespace gnut;

namespace great
{
	class OI_LIBRARY_EXPORT t_gorbsp3
	{

	public:
		t_gorbsp3();
		t_gorbsp3(vector<string> mark, t_gsetbase* set, t_glog * log);
		virtual ~t_gorbsp3();
	public:
		bool processBatch(const t_gtime& beg, const t_gtime& end);
		void setDataUsed(t_gorb* gorb);
		void setDataUsed(t_gpoleut1* poleut1);
		void GenerateProduct();
		t_gallics_leo*    _ics_leo;
		t_gallics_nav*    _ics_nav;
	protected:
		bool _initial_data();
		bool _get_settings();
		bool _init_sp3();
		bool _extend_time_period(t_gtime beg, t_gtime end, t_gtime& extend_beg, t_gtime& extend_end);
		void _generate_sp3file();
		double _intv;
		map<string,map<string, double>>  _ref_state;
		vector<string>    _sats;
		vector<string>    _sat3;
		string            _sattype;
		string            _frame;
		string            _PV;
		t_glog*		      _glog;             ///< Genereal   log output. 
		t_gsetbase*		  _set;              ///< Base setting. 
		t_gorb*           _gorb;             ///< orb data (output)
		t_gtime           _epoch;          ///< current epoch
		t_gallprec*       _sp3;
		t_gpoleut1*       _poleut1;
		t_gtime           _beg;
		t_gtime           _end;
		t_gtrs2crs*       _trs2crs;
		bool              _flag;
		int               _integra_degree;
		int               _step_inp;
		int               _step_out;
		int               _neq_size;
		int               _maxsat3;
		map<string, int>  _drag_para;
		map<string, int>  _cust_para;
		set<string>       _param;
		int               _param_size;

	};
}
#endif // !GSP3ORB_H