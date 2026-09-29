/**
*
* @verbatim
	History
	 -1.0 Hongmin Zhang		2020-03-20  creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gkbr.cpp
* @brief		KBR
*
* @author       Hongmin Zhang, Wuhan University
* @version		1.0.0
* @date		    2020-03-20
*
*/


#include <algorithm>
#include <set>
#include "gkbr.h"

#include "gset/gsetout.h" 
#include "gset/gsetorbdif.h"
#include "gutils/gtrs2crs.h"
#include "gutils/gcrs2acr.h"
#include "gutils/gleo.h"
#include "gutils/gtrs2crs.h"

#include "gcoders/kbr.h"

using namespace std;

namespace great
{
	t_gkbr::t_gkbr()
	{
	}

	t_gkbr::t_gkbr(t_gsetbase* set, t_glog* log)
	{
		_set = set;
		_glog = log;
		_tide = new t_gtideIERS(log);
		//setOBJ(dynamic_cast<t_gallobj*>((*data)[t_gdata::ALLOBJ]));
		//_result = new t_gresult();
	}

	t_gkbr::~t_gkbr()
	{

	}

	void t_gkbr::setORB(t_gallorb* gorb)
	{
		_gorb = gorb;
	}

	void t_gkbr::setKBR(t_gallkbr* gkbr)
	{
		_gkbr = gkbr;
	}

	void t_gkbr::setLEO(set<string>& leo)
	{
		_leo = leo;
	}

	bool t_gkbr::processBatch(t_gtime& beg,t_gtime& end)
	{
		string function_name = "t_kbr::processingBatch...";

		t_gtime kbr_beg = _gkbr->getBeg();
		t_gtime kbr_end = _gkbr->getEnd();

		double interval = _gkbr->getIntv();
		vector<string> sats(2);
		sats = _gkbr->getSat();

		double _xsat1[3], _xsat2[3];
		t_gtriple xsat1, xsat2;
		double dist;
		double omc;

		map<int,map<t_gtime, double>> kbromc;
		
		int intv;

		if (beg > kbr_end || end < kbr_beg)
		{
			_glog->comment(0, function_name, "ERROR :The wrong input time!");
			return false;
		}

		_start = kbr_beg;
		_stop = kbr_end;

		if (beg > kbr_beg)
		{
			intv = (int)((beg - kbr_beg) / interval);
			_start = kbr_beg + interval * (intv + 1);
		}

		if (end < kbr_end)
		{
			intv = (int)((kbr_end - end) / interval);
			_stop = kbr_end - interval * (intv + 1);
		}

		int tab = 0;
		for (auto iter = _leo.begin(); iter != _leo.end(); iter++)
		{
			if (strs2strl(*iter) == sats[0] || strs2strl(*iter) == sats[1])
			{
				tab = 1;
			}
		}
		if (tab == 0)
		{
			_glog->comment(0, function_name, "ERROR :The wrong input leo names!");
			return false;
		}

		t_gtime TI;
		int i = 0;
		int tab_kbr = 0;
		int arcobs[5] = { 0 };
		stringstream str, str1;
		string mjd, sod;

		for (TI = _start; TI <= _stop;)
		{
			if (_gkbr->getkbr(TI)->flag == 1)
			{
				str << TI.mjd();
				str >> mjd;
				str1 << TI.sod();
				str1 >> sod;
				_glog->comment(0, function_name, "Warning : New Arc at  "+ mjd +"  "+ sod);
				tab_kbr += 1;
				if (tab_kbr >= 5)
				{
					_glog->comment(0, function_name, "Error: arc number is larger than maximum number 5");
					return false;
				}
			}
			if (_gkbr->getkbr(TI)->flag > 1)
			{
				str << TI.mjd();
				str >> mjd;
				str1 << TI.sod();
				str1 >> sod;
				_glog->comment(0, function_name, "Warning: Remove bad KBR observation at  " + mjd + "  " + sod+"  with quality flag> 1");
				TI = TI + interval;
				continue;
			}
			_gorb->get_crs_sat_crd(TI, sats[0], _xsat1[0], _xsat1[1], _xsat1[2]);
			_gorb->get_crs_sat_crd(TI, sats[1], _xsat2[0], _xsat2[1], _xsat2[2]);
		     xsat1.set(_xsat1);
			 xsat2.set(_xsat2);
			 dist = (xsat1.crd_cvect() - xsat2.crd_cvect()).NormFrobenius();
			 omc = _gkbr->getkbr(TI)->bias_range + _gkbr->getkbr(TI)->ltime_cor + _gkbr->getkbr(TI)->ant_cor - dist;
			 //_gkbromc->abc[i] = omc;
			 kbromc[tab_kbr][TI] = omc;
			 arcobs[tab_kbr] += 1;
			 TI = TI + interval;
		}

		double aver[5] = { 0.0 };
		map<t_gtime, double> kbr_temp;
		for (auto iter = kbromc.begin(); iter != kbromc.end(); iter++)
		{
			kbr_temp = iter->second;
			for (auto iterT = kbr_temp.begin(); iterT != kbr_temp.end(); iterT++)
			{
				aver[iter->first] += iterT->second / arcobs[iter->first];
			}
			
		}


		for (auto iter = kbromc.begin(); iter != kbromc.end(); iter++)
		{
			kbr_temp = iter->second;
			for (auto iterT = kbr_temp.begin(); iterT != kbr_temp.end(); iterT++)
			{
				iterT->second = iterT->second - aver[iter->first];
			}
			
		}

		for (auto iter = kbromc.begin(); iter != kbromc.end(); iter++)
		{
			for (auto iterT = kbr_temp.begin(); iterT != kbr_temp.end(); iterT++)
			{
				_gkbr->setOMC(iterT->first, iterT->second);
			}
			
		}

		return true;
	}

	void t_gkbr::GenerateProduct()
	{
		t_gfile gout;
		string kbr_path = "";

		kbr_path = dynamic_cast<t_gsetout*>(_set)->outputs("kbromc");

		if (kbr_path.empty()) {
			kbr_path = "kbromc";
		}

		gout.path(kbr_path);
		gout.glog(_glog);

		t_kbr kbr_coder(_set);
		kbr_coder.add_data("ID2", _gkbr);
		gout.coder(&kbr_coder);

		gout.run_write();
		return;
	}
}