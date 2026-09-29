/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambfix.h
* @brief	for the use of fix ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*
*/

#ifndef GAMBFIX_H
#define GAMBFIX_H 

#include "gutils/ExportAMB.h"
#include "gset/gsetbase.h"
#include "gproc/glsq.h"
#include "gset/gsetamb.h"
#include "gall/gallpar.h"
#include "gall/gallobj.h"
#include "gdata/gupd.h"
#include "gambfix/gallamb.h"
#include "gambfix/gambOw.h"
#include "gambfix/gambBlRec.h"
#include "gambfix/gambBlSat.h"

using namespace std;
using namespace gnut;
namespace great
{
	/**
	*@ brief class for LAMADA search finding fixing ambiguity.
	*/
	class AMB_LIBRARY_EXPORT t_gambfix
	{
	
	public:

		/** @brief default constructor. */
		t_gambfix();
		t_gambfix(t_glog* log, t_gsetbase* set);

		void add_data(t_gallobs* obs, t_gallbias* bias, t_gallpar* pars, t_gallobj* objs);
		void add_data(t_gallobs* obs, t_gallbias* bias, t_gallpar* pars,  t_gupd* upd, t_gallobj* objs);
		void add_data(t_gallobs* obs);
		void add_data(t_gallbias* bias);
		void add_data(t_gallpar* pars);
		void add_data(t_gallobj* objs);

		/** @brief default destructor. */
		virtual ~t_gambfix();

		// 确定beg到end时间内的双差模糊度
		bool ProcessBatch(t_gtime& beg, t_gtime& end);
	protected:
		t_glog*      _glog = nullptr;
		t_gsetbase*  _gset = nullptr;
		t_gallobs*   _gobs = nullptr;
		t_gallbias*  _gbias = nullptr;
		t_gallpar*   _gpar  = nullptr;
		t_gallobj*   _gobj  = nullptr;
		t_gupd*      _gupd  = nullptr;
		t_gallambs   _gambs;

		std::set<string>     _rec_list;
		std::set<string>     _sat_list;

		t_gtime _beg;    ///< begin epoch
		t_gtime _end;    ///< end epoch
		t_gtime _crt;    ///< current epoch
		double  _intv;  ///< sampling interval /s
		
		t_gbdeci* _evaluator_EWL = nullptr;
		t_gbdeci* _evaluator_WL = nullptr;
		t_gbdeci* _evaluator_NL = nullptr;
		// bool      _apply_irc = false;
		t_gsetamb::FIX_MODE _fix_mode;                          ///< set ambiguity fixing mode
		t_gsetamb::UPD_MODE _upd_mode;							///< set ambiguity fixing mode
		double _ratio;                                          ///< threshold in LAMBDA method
		double _min_common_time;                                ///< the Minimum common time of two observation arc
		double _wl_interval;                                    ///< widelane ambiguity interval
		double _max_baseline_length;                            ///< the limit of baseline length
		map<string, double> _map_EWL_decision;                   ///< deriation, sigma in EWL-cycle
		map<string, double> _map_WL_decision;                   ///< deriation, sigma in WL/NL-cycle
		map<string, double> _map_NL_decision;

		double getWaveLen(string sat, string type);
		//bool _get_Dd_amb(Sd_amb);

		string _class_id = "t_gambfix";
	private:


	};
}

#endif