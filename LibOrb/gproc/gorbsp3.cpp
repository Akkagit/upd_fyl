/**
*
* @verbatim
	History
	 -1.0 Yujie Qin		2019-12-02 creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gorbsp3.cpp
* @brief		convert orb file to sp3 file
*
* @author       Yujie Qin	, Wuhan University
* @version		1.0.0
* @date		    2019-12-02
*
*/

#include <iostream>
#include <iomanip>
#include <algorithm>
#include "newmat/newmatio.h"
#include "gproc/gorbsp3.h"
using namespace gnut;
namespace great
{
	t_gorbsp3::t_gorbsp3()
	{

	}
	t_gorbsp3::t_gorbsp3(vector<string> mark, t_gsetbase* set, t_glog * log)
	{
		_flag = false;
		_sats = mark;
		_set = set;
		_glog = log;
		_sp3 = new t_gallprec();
		_epoch.tsys(t_gtime::GPS);
		_maxsat3 = 119;
	}

	t_gorbsp3::~t_gorbsp3()
	{
		delete _sp3; _sp3 = NULL;
		
	}

	void t_gorbsp3::setDataUsed(t_gorb* orb)
	{
		_gorb = orb;
	}
	void t_gorbsp3::setDataUsed(t_gpoleut1* poleut1)
	{
		_poleut1 = poleut1;
	}
	
	bool t_gorbsp3::processBatch(const t_gtime& beg, const t_gtime& end)
	{
		_PV = "V";
		if (!_initial_data())
		{
			_glog->comment(0, "orbsp3", "ERROR : orbsp3 have no comeplete file !");
			return false;
		}
		if (!_get_settings())
		{
			_glog->comment(0, "orbsp3", "ERROR : orbsp3 have something wrong in xml-file !");
			return false;
		}
		_extend_time_period(beg, end, _beg, _end);
		if (!_init_sp3())
		{
			_glog->comment(0, "orbsp3", "ERROR : orbsp3 have something wrong in initial sp3-file !");
			return false;
		}

		_trs2crs = new t_gtrs2crs(false, _poleut1);
		
		_epoch = _beg;
		Matrix rot_trs2crs(0.0, 3, 3);
		Matrix rot_dgmat(0.0, 3, 3);
		while (_epoch <= _end)
		{
			t_gtime tdt_gtime = _epoch;
			tdt_gtime.tsys(t_gtime::TT);
			if (_frame == "trs")
			{
				_trs2crs->calcRotMat(tdt_gtime, true, true, true);
				// rotation matrix of TRS to CRS
				rot_trs2crs = _trs2crs->getRotMat();
				
				rot_dgmat = _trs2crs->getMatDu();
			}
			
			t_gtriple target_pos, target_vel;
			ColumnVector trs_pos(3), trs_vel(3), crs_pos(3), crs_vel(3);
			for (int i =0; i<_sats.size();i++)
			{
				string sat = _sats[i];
				if (!_gorb->get_xyz(_epoch, sat, crs_pos(1), crs_pos(2), crs_pos(3)))
				{
					continue;
				}
				_gorb->get_vel(_epoch, sat, crs_vel(1), crs_vel(2), crs_vel(3));
				crs_pos = crs_pos / 1000;
				crs_vel = crs_vel / 1000;

				if (_frame == "trs")
				{
					trs_pos = (crs_pos.t() * rot_trs2crs).t();
					trs_vel = (crs_vel.t() * rot_trs2crs).t();
					ColumnVector xhelp(3);
					xhelp = (trs_pos.t() * rot_dgmat / RAD2TSEC).t();
					trs_vel = trs_vel + xhelp;
					for (int i = 0; i < 3; i++)
					{
						target_pos[i] = trs_pos(i + 1);
						target_vel[i] = trs_vel(i + 1);
					}
				}
				else
				{
					for (int i = 0; i < 3; i++)
					{
						target_pos[i] = crs_pos(i + 1);
						target_vel[i] = crs_vel(i + 1);
					}
				}
				if (_PV == "V")
				{
					t_gtriple dxyz(0, 0, 0);
					double xyzt[4] = { target_vel[0],target_vel[1],target_vel[2],0};
					double var[4] = {0,0,0,0 };
					_sp3->add_pos_vel(_sat3[i], _epoch, target_pos, 0.0, dxyz, 0.0, xyzt, var);
				}
				else
				{
					t_gtriple dxyz(0, 0, 0);
					_sp3->addpos(_sat3[i], _epoch, target_pos, 0.0, dxyz, 0);
				}
				
			}
			
			
			_epoch = _epoch + _intv;
		}
		if (_trs2crs)
		{
			delete _trs2crs;
			_trs2crs = nullptr;
		}
		return true;
	}

	bool t_gorbsp3::_initial_data()
	{
		if (!_poleut1 || !_gorb)
			return false;
		return true;
	}

	bool t_gorbsp3::_get_settings()
	{
		_intv = dynamic_cast<t_gsetorbsp3*>(_set)->intv();
		_frame = dynamic_cast<t_gsetorbsp3*>(_set)->frame();
		if (_frame == "")
		{
			_frame = "trs";
		}
	    return true;
		
	}

	bool t_gorbsp3::_init_sp3()
	{
		t_orb_hd orb_hd;
		_gorb->get_head(orb_hd);
		t_gtime orb_beg = orb_hd.beg;
		t_gtime orb_end = orb_hd.end;
		/*if (_beg < orb_beg)_beg = orb_beg;
		if (_end > orb_beg)_end = orb_end;*/
		_beg = orb_beg;
		_end = orb_end;
		_sp3->add_interval(_intv);
		
		_sp3->set_beg(_beg);
		_sp3->set_end(_end);
		set<string> sats = orb_hd.sat_prn;
		_sattype = orb_hd.sat_group;
		vector<string> sats_used;
		vector<string>::const_iterator it_beg = _sats.begin();
		vector<string>::const_iterator it_end = _sats.end();
		vector<string>::const_iterator iter = it_beg;
		for (iter; iter != it_end; iter++)
		{
			
			if (sats.count(*iter) <= 0) continue;
			sats_used.push_back(*iter);

		}
		if (sats_used.size() == 0)
		{
			return false;
		}
		int nsize = sats_used.size();
		if (_sattype != "LEO")
		{
			_sat3 = sats_used;
			_sats = sats_used;
		}
		else
		{
			for (int i = 0; i < nsize; i++)
			{
				string sat3 = "L" + int2str(i + 1, 2);
				_sat3.push_back(sat3);
			}
		}
		for (int i = nsize; i < _maxsat3; i++)
		{
			_sat3.push_back("  0");
		}
		it_beg = sats_used.begin();
		it_end = sats_used.end();
		iter = it_beg;
		_sats = sats_used;
		_sp3->set_sat(_sats);
		_sp3->set_sat3(_sat3);
		_sp3->set_data_type("V");
		_sp3->set_sat_type(_sattype);
		_sp3->set_agency("GNUT");
		_sp3->add_interval(_sat3[0], _intv);
		return true;
		
	}
	bool t_gorbsp3::_extend_time_period(t_gtime beg, t_gtime end, t_gtime& extend_beg, t_gtime& extend_end)
	{
		try
		{
			int sec_beg = int((beg.dsec(false) + beg.sod(false)) / _intv) *  _intv;
			int sec_end = int((end.dsec(false) + end.sod(false)) / _intv) *  _intv;
			extend_beg.from_mjd(beg.mjd(false), sec_beg, sec_beg - floor(sec_beg), false);
			extend_end.from_mjd(end.mjd(false), sec_end, sec_end - floor(sec_end), false);
			return true;
		}
		catch (...)
		{
			return false;
		}

	}


	
	void t_gorbsp3::_generate_sp3file()
	{
		
	}

	void t_gorbsp3::GenerateProduct()
	{
		_generate_sp3file();
		string tmp = dynamic_cast<t_gsetout*>(_set)->outputs("sp3");
		if (tmp.empty())
		{
			tmp = "sp3_file";
		}

		t_gio* gout_sp3 = new t_gfile;  t_gcoder* gsp3_coder = new t_sp3(_set, "", 4096);
		string sp3_path = tmp;
		//string ics_path = (dynamic_cast<t_gsetinp*>(_set)->inputs("ics"))[0];
		gout_sp3->path(sp3_path); gout_sp3->glog(_glog);
		gsp3_coder->clear(); gsp3_coder->glog(_glog); gsp3_coder->path(sp3_path);
		gsp3_coder->add_data("ID1", _sp3);
		gout_sp3->coder(gsp3_coder);
		gout_sp3->run_write();
		delete gout_sp3; gout_sp3 = nullptr;
		delete gsp3_coder; gsp3_coder = nullptr;

	}
}
