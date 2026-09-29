/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambfix.cpp
* @brief	function for ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#include "gambfix/gambfix.h"
#include <utility>

using namespace gnut;

namespace great
{
	t_gambfix::t_gambfix()
	{
	}

	t_gambfix::t_gambfix(t_glog* log, t_gsetbase* set)
	{
		cerr << "================> constructor of t_gambfixDd" << endl;

		this->_glog = log;
		this->_gset = set;
	
		if (_gset)
		{
			auto gsetamb = dynamic_cast<t_gsetamb*>(_gset);
			if (gsetamb)
			{
				this->_fix_mode            = gsetamb->fix_mode();
				this->_upd_mode            = gsetamb->upd_mode();
				this->_ratio               = gsetamb->lambda_ratio();
				this->_min_common_time     = gsetamb->min_common_time();
				this->_wl_interval         = gsetamb->wl_interval();
				this->_max_baseline_length = gsetamb->max_baseline_length();
				this->_map_EWL_decision    = gsetamb->get_amb_decision("EWL");
				this->_map_WL_decision     = gsetamb->get_amb_decision("WL");
				this->_map_NL_decision     = gsetamb->get_amb_decision("NL");
				//this->_apply_irc           = gsetamb->apply_irc();
			}

			auto gsetgen = dynamic_cast<t_gsetgen*>(_gset);
			if (gsetgen)
			{
				this->_beg = gsetgen->beg();
				this->_end = gsetgen->end();
				this->_intv = gsetgen->sampling();
			}

		}
		this->_evaluator_EWL = new t_gbdeci(this->_map_EWL_decision["maxdev"], this->_map_EWL_decision["maxsig"],
			0.0, this->_map_EWL_decision["alpha"]);
		this->_evaluator_WL = new t_gbdeci(this->_map_WL_decision["maxdev"], this->_map_WL_decision["maxsig"],
			0.0, this->_map_WL_decision["alpha"]);
		this->_evaluator_NL = new t_gbdeci(this->_map_NL_decision["maxdev"], this->_map_NL_decision["maxsig"],
			0.0, this->_map_NL_decision["alpha"]);
	}

	void t_gambfix::add_data(t_gallobs* obs, t_gallbias* bias, t_gallpar* pars, t_gallobj* objs)
	{
		cerr << " ============> t_gambfix::add_data <===============" << endl;

		this->_gobs  = obs;
		this->_gbias = bias;
		this->_gpar = pars;
		this->_gobj = objs;


	}
	void t_gambfix::add_data(t_gallobs* obs, t_gallbias* bias, t_gallpar* pars, t_gupd* upd, t_gallobj* objs)
	{
		cerr << " ============> t_gambfix::add_data <===============" << endl;
		this->_gobs = obs;
		this->_gbias = bias;
		this->_gpar = pars;
		this->_gupd = upd;
		this->_gobj = objs;
	}

	void t_gambfix::add_data(t_gallobs* obs)
	{
		this->_gobs = obs;
	}

	void t_gambfix::add_data(t_gallbias* bias)
	{
		this->_gbias = bias;
	}

	void t_gambfix::add_data(t_gallpar* pars)
	{
		this->_gpar = pars;
	}

	void t_gambfix::add_data(t_gallobj* objs)
	{
		this->_gobj = objs;
	}

	t_gambfix::~t_gambfix()
	{
		if (_evaluator_WL)
		{
			delete _evaluator_WL;
			_evaluator_WL = nullptr;
		}

		if (_evaluator_NL)
		{
			delete _evaluator_NL;
			_evaluator_NL = nullptr;
		}
	}


	bool t_gambfix::ProcessBatch(t_gtime& beg, t_gtime& end)
	{
		return false;
	}


	double t_gambfix::getWaveLen(string sat, string type)
	{
		t_gobsgnss* gnss = new t_gobsgnss();
		double l = 0;
		if (sat.substr(0, 1) == "G")
		{
			gnss->sat("G");
			if (type == "L1") l = gnss->wavelength(BAND_1);
			if (type == "L2") l = gnss->wavelength(BAND_2);
			if (type == "L3") l = gnss->wavelength(BAND_5);
			if (type == "WL") l = gnss->wavelength_WL(BAND_1, BAND_2);
			if (type == "NL") l = gnss->wavelength_NL(BAND_1, BAND_2);
			if (type == "EWL") l = gnss->wavelength_WL(BAND_2, BAND_5);
		}
		else if (sat.substr(0, 1) == "E")
		{
			gnss->sat("E");
			if (type == "L1") l = gnss->wavelength(BAND_1);
			if (type == "L2") l = gnss->wavelength(BAND_5);
			if (type == "L3") l = gnss->wavelength(BAND_6);
			if (type == "WL") l = gnss->wavelength_WL(BAND_1, BAND_5);
			if (type == "NL") l = gnss->wavelength_NL(BAND_1, BAND_5);
			if (type == "EWL") l = gnss->wavelength_WL(BAND_5, BAND_6);
		}
		else if (sat.substr(0, 1) == "C")
		{
			gnss->sat("C");
			if (type == "L1") l = gnss->wavelength(BAND_2);
			if (type == "L2") l = gnss->wavelength(BAND_7);
			if (type == "L3") l = gnss->wavelength(BAND_6);
			if (type == "WL") l = gnss->wavelength_WL(BAND_2, BAND_7);
			if (type == "NL") l = gnss->wavelength_NL(BAND_2, BAND_7);
			if (type == "EWL") l = gnss->wavelength_WL(BAND_7, BAND_6);
		}

		else
		{

		}
		delete gnss;
		return l;
	}

}

