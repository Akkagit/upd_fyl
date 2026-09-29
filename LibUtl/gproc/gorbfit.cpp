/**
*
* @verbatim
	History
	 -1.0 Jie Li		2019-09-06  creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gorbfit.cpp
* @brief		orbit fit processing for GNS satellite
* @author       Jie Li, Wuhan University
* @version		1.0.0
* @date		    2019-09-06
*
*/


#include <algorithm>
#include <set>
#include "gorbfit.h"

#include "gset/gsetinp.h"
#include "gset/gsetout.h" 
#include "gset/gsetorbdif.h"
#include "gutils/gtrs2crs.h"
#include "gutils/gcrs2acr.h"

#include "gio/gfile.h"
#include "gcoders/icsnav.h"
using namespace std;

namespace great {

	t_gorbfit::t_gorbfit():t_gorbdif()
	{
	}

	t_gorbfit::t_gorbfit(t_gsetbase * set, t_glog * log)
		:t_gorbdif(set,log)
	{
	
	}

	t_gorbfit::~t_gorbfit()
	{

	}

	bool t_gorbfit::processBatch()
	{
		string function_name = "t_gorbfit::processBatch";

		if (!_init_data())
		{
			_glog->comment(0, function_name, "ERROR : init data fail, check your xml file.");
			return false;
		}
		////check time
		if (!_check_time(""))
		{
			return false;
		}



		// loop by epoch
		ostringstream os;
		t_gtime epo_now = _beg;
		while (epo_now < _end)
		{
			os << " GPST" << setw(15) << right << epo_now.mjd()
				<< setw(20) << fixed << right << epo_now.sod() << endl;

			t_glsqEquationMatrix equ;
			Matrix B; DiagonalMatrix P; ColumnVector l;


			// loop by sat in xml
			for (auto iter = _sat.begin(); iter != _sat.end();)
			{
				
				// skip sat 
				set<string> sat_sp3 = _gsp3->satellites();
				if (sat_sp3.find(*iter) == sat_sp3.end())
				{
					_glog->comment(1, function_name, "WARNING : " + *iter + "can not be found in sp3 file!");
					_sat.erase(iter++);
					continue;
				}
				set<string> sat_ics = _gics->sat_list;
				if (sat_ics.find(*iter) == sat_ics.end())
				{
					_glog->comment(1, function_name, "WARNING : " + *iter + "can not be found in ics(orb) file!");
					_sat.erase(iter++);
					continue;
				}

				//get sat position in sp3
				double pos_trs_tmp[3];
				_gsp3->pos(*iter, epo_now, pos_trs_tmp);
		
				// jdhuang : sp3中没有位置信息，读rinex文件
				if (double_eq(pos_trs_tmp[0], 0.0) && double_eq(pos_trs_tmp[1], 0.0) && double_eq(pos_trs_tmp[2], 0.0))
				{
					_gsp3->t_gallnav::pos(*iter, epo_now, pos_trs_tmp);
				}

				ColumnVector pos_trs_sp3(0.0, 3);
				pos_trs_sp3(1) = pos_trs_tmp[0];
				pos_trs_sp3(2) = pos_trs_tmp[1];
				pos_trs_sp3(3) = pos_trs_tmp[2];
				//pos_trs_sp3 = pos_trs_sp3 / 1E3;
				//sp3值检查
				// 0 值或大值
				if ((abs(pos_trs_sp3(1)) + abs(pos_trs_sp3(2)) + abs(pos_trs_sp3(3))) == 0)
				{
					_sat.erase(iter++);
					continue;
				}
				if (pos_trs_sp3(1) >= 1E18)
				{
					_sat.erase(iter++);
					continue;
				}
				t_gtime epo_now_tdt = epo_now;
				epo_now_tdt.tsys(t_gtime::TT);
				//ChangeTimeSys(GPS2TDT, epo_now, epo_now_tdt);
				t_gtrs2crs*  trs2crs_cord = new t_gtrs2crs(false, _gpoleut1);
				trs2crs_cord->calcRotMat(epo_now_tdt, true, true, true);
				Matrix rot_trs2crs = trs2crs_cord->getRotMat();

				double diff = (epo_now - _ref_orb) / 86400.0;
				Matrix dxdpole(0.0, 3, 4);
				dxdpole.Column(1) = (trs2crs_cord->getMatDx()*pos_trs_sp3) / RAD2SEC;
				dxdpole.Column(2) = dxdpole.Column(1) * diff;
				dxdpole.Column(3) = (trs2crs_cord->getMatDy()*pos_trs_sp3) / RAD2SEC;
				dxdpole.Column(4) = dxdpole.Column(3) * diff;
				Matrix dxdut1(0.0, 3, 2);
				dxdut1.Column(1) = (trs2crs_cord->getMatDu()*pos_trs_sp3) / RAD2TSEC;
				dxdut1.Column(2) = dxdut1.Column(1)*diff;

				ColumnVector pos_crs_sp3(0.0, 3);
				pos_crs_sp3 = rot_trs2crs * pos_trs_sp3;

				// get sat position in orbfile
				ColumnVector pos_crs_orb(0.0, 3);
				_gorb->get_crs_sat_crd(epo_now, *iter, pos_crs_orb(1), pos_crs_orb(2), pos_crs_orb(3));
				//pos_crs_orb = pos_crs_orb / 1E3;
				Matrix partial;
				_gorb->getGnsPartial(epo_now, *iter, partial);
				//cout << "partial  " << partial << endl;

				int nPar = _lsqEstimator->_x_solve.parNumber();
				// resize the Matrix B l P
				B.ReSize(3, nPar); B = 0.0;
				l.ReSize(3);	l = 0.0;
				P.ReSize(3);	P = 1.0;

				l = pos_crs_sp3 - pos_crs_orb;

				// 广播星历容错
				if ((abs(l(1)) + abs(l(2)) + abs(l(3))) >= 1e4)
				{
					cout << *iter<<"   "<<"abs 1 : " << abs(l(1)) << "abs 2 : " << abs(l(2)) << "abs 3 : " << abs(l(3)) << endl;
					iter++;
					continue;
				}

				//cout << "l" << l << endl;
				// B for ICS and ERP
				for (int i = 1; i <= nPar; i++)
				{
					t_gpar par = _lsqEstimator->_x_solve[i - 1];
					if (par.prn != "")
					{
						if (*iter == par.prn)
						{
							auto ics = _gics->get_nav_ics(*iter);
							auto ics_list =ics->get_par_list();
							int ipar = 0;
							for (auto itpar = ics_list.begin(); itpar != ics_list.end(); itpar++)
							{
								ipar++;
								if (par.parType ==ics->_icspar2gpar((*itpar).first))
								{
									//cout << "partial  " << partial.row(ipar).t() * 1E3 << endl;
									B.column(i) = partial.row(ipar).t()*1E3;
									//cout << "B.column(i)  " << B.column(i) * 1E3 << endl;
									break;
								}
							}
						}
						continue;
					}
					
					switch (par.parType)
					{
					case t_gpar::SCALE:
						B.column(i) = pos_crs_orb / 1E6;
						break;
					case t_gpar::GEOCX:
						B(1, i) = 1000.0;
						break;
					case t_gpar::GEOCY:
						B(2, i) = 1000.0;
						break;
					case t_gpar::GEOCZ:
						B(3, i) = 1000.0;
						break;
					case t_gpar::XPOLE:
						B.column(i) = -1.0 * (dxdpole.Column(1));
						break;
					case t_gpar::YPOLE:
						B.column(i) = -1.0 * (dxdpole.Column(3));
						break;
					case t_gpar::UT1:
						B.column(i) = -1.0 * (dxdut1.Column(1));
						break;
					case t_gpar::DXPOLE:
						B.column(i) = -1.0 * (dxdpole.Column(2));
						break;
					case t_gpar::DYPOLE:
						B.column(i) = -1.0 * (dxdpole.Column(4));
						break;
					case t_gpar::DUT1:
						B.column(i) = -1.0 * (dxdut1.Column(2));
						break;
					default:
						break;
					}
				}

				ColumnVector cor(0.0, 3);
				if (_lsqEstimator->dx().size() != 0)
				{
					cor = B * _lsqEstimator->dx();
				}
				l = l - cor;


				vector<string> site_info = { "X","Y","Z" };
				equ.add_equ(B, P, l, site_info,vector<string>(3,*iter),vector<t_gobscombtype>(3,t_gobscombtype()));

				// output APRI to orbdif file 
				os << "APRI" << setw(5) << right << *iter
					<< setw(20) << right << epo_now.mjd()
					<< setw(10) << fixed << setprecision(3) << right << epo_now.sod();
				for (int i = 1; i <= 3; i++)
				{
					os << setw(15) << fixed << setprecision(3) << right << l(i);
				}
				for (int i = 1; i <= 3; i++)
				{
					os << setw(15) << fixed << setprecision(3) << right << cor(i);
				}
				os << endl;

				if (trs2crs_cord)
				{
					delete trs2crs_cord;
					trs2crs_cord = NULL;
				}

				iter++;
			}
			//_lsqEstimator->write_equation(equ, epo_now);
			_lsqEstimator->add_equation(equ, epo_now);
			epo_now = epo_now + _sampling;
		}

		_lsqEstimator->solve_NEQ();


		//write the parameters value and correction to differ file
		os << setw(20) << left << "Paramters Estimated"
			<< setw(20) << right << "Initial Value"
			<< setw(20) << right << "Correction"
			<< setw(20) << right << "Result" << endl;
		int npar = _lsqEstimator->_x_solve.parNumber();
		ColumnVector dx = _lsqEstimator->dx();
		for (int i = 0; i < npar; i++)
		{
			t_gpar par = _lsqEstimator->_x_solve[i];
			os << setw(20) << left << par.str_type()
				<< setw(20) << right << scientific << setprecision(10) << par.value()
				<< setw(20) << right << scientific << setprecision(10) << dx(i + 1)
				<< setw(20) << right << scientific << setprecision(10) << par.value() + dx(i + 1) << endl;
		}
		_orbdif->open(_orbdif->mask(), ios::out);
		_orbdif->write(os.str().c_str(), os.str().size());
		_orbdif->flush();

		// 残差统计并写入文件
		_recover_resid();

		//更新ICS文件
		_update_icsfile();
		if (_lsqEstimator)
		{
			delete _lsqEstimator;
			_lsqEstimator = NULL;
		}

		return true;
	}

	bool t_gorbfit::_init_data()
	{
		if (!_gics|| !_gorb || !_gsp3 || !_gpoleut1)
		{
			_glog->comment(0, "orbfit", "ERROR : orbfit ==> t_gorbfit::no comeplete file !");
			return false;
		}
		// get begin/end time,sampling, sat list and trans mode from xml
		//if (!_get_setting())
		if (!_get_setting_orbfit())			 
		{
			_glog->comment(0, "orbfit", "ERROR : orbfit ==> t_gorbfit::there is something wrong witn xml !");
			return false;
		}

		// init lsq
		if (_lsqEstimator) {
			delete _lsqEstimator;
		}
		_lsqEstimator = new t_glsq(NULL);


		_init_SAT_par();
		_init_ERP_par();

		_set_out();

		return true;
	}

	void t_gorbfit::_init_ERP_par()
	{
		int ipar_num = _lsqEstimator->_x_solve.parNumber();

		for (int i = 0; i < _trans.length(); i++)
		{
			map<string, double> data1, data2;
			_gpoleut1->getEopData(_beg, data1);
			_gpoleut1->getEopData(_end, data2);
			if (_trans.substr(i, 1) == "S")
			{
				t_gpar par_scale("", t_gpar::SCALE, ++ipar_num, "");
				par_scale.setTime(_beg, _end);
				par_scale.value(0.0);
				par_scale.apriori(1E20);

				_lsqEstimator->add_parameter(par_scale);
			}
			else if (_trans.substr(i, 1) == "T")
			{
				t_gpar par_geocx("", t_gpar::GEOCX, ++ipar_num, "");
				t_gpar par_geocy("", t_gpar::GEOCY, ++ipar_num, "");
				t_gpar par_geocz("", t_gpar::GEOCZ, ++ipar_num, "");

				par_geocx.setTime(_beg, _end);
				par_geocx.value(0.0);
				par_geocx.apriori(1E20);
				par_geocy.setTime(_beg, _end);
				par_geocy.value(0.0);
				par_geocy.apriori(1E20);
				par_geocz.setTime(_beg, _end);
				par_geocz.value(0.0);
				par_geocz.apriori(1E20);

				_lsqEstimator->add_parameter(par_geocx);
				_lsqEstimator->add_parameter(par_geocy);
				_lsqEstimator->add_parameter(par_geocz);
			}
			else if (_trans.substr(i, 1) == "R")
			{
				t_gpar par_xpole("", t_gpar::XPOLE, ++ipar_num, "");
				t_gpar par_ypole("", t_gpar::YPOLE, ++ipar_num, "");
				t_gpar par_ut1("", t_gpar::UT1, ++ipar_num, "");

				par_xpole.setTime(_beg, _end);
				par_xpole.value(data1["XPOLE"]);
				par_xpole.apriori(1E20);
				par_ypole.setTime(_beg, _end);
				par_ypole.value(data1["YPOLE"]);
				par_ypole.apriori(1E20);
				par_ut1.setTime(_beg, _end);
				par_ut1.value(data1["UT1-TAI"]);
				par_ut1.apriori(1E20);

				_lsqEstimator->add_parameter(par_xpole);
				_lsqEstimator->add_parameter(par_ypole);
				_lsqEstimator->add_parameter(par_ut1);
			}
			else if (_trans.substr(i, 1) == "D")
			{
				t_gpar par_dxpole("", t_gpar::DXPOLE, ++ipar_num, "");
				t_gpar par_dypole("", t_gpar::DYPOLE, ++ipar_num, "");
				t_gpar par_dut1("", t_gpar::DUT1, ++ipar_num, "");

				par_dxpole.setTime(_beg, _end);
				par_dxpole.value(data2["XPOLE"] - data1["XPOLE"]);
				par_dxpole.apriori(1E20);
				par_dypole.setTime(_beg, _end);
				par_dypole.value(data2["YPOLE"] - data1["YPOLE"]);
				par_dypole.apriori(1E20);
				par_dut1.setTime(_beg, _end);
				par_dut1.value(data2["UT1-TAI"] - data1["UT1-TAI"]);
				par_dut1.apriori(1E20);

				_lsqEstimator->add_parameter(par_dxpole);
				_lsqEstimator->add_parameter(par_dypole);
				_lsqEstimator->add_parameter(par_dut1);
			}
		}
	}

	void t_gorbfit::_init_SAT_par()
	{
		set<string> ics_sat = _gics->sat_list;
		int ipar_num = _lsqEstimator->_x_solve.parNumber();
		
		for (auto itsat = ics_sat.begin(); itsat != ics_sat.end(); ++itsat)
		{
			if (_sat.find(*itsat) == _sat.end())
			{
				_glog->comment(1, "orbfit:init_SAT_par", "WARNING : " + *itsat + "is not in the sat list of xml");
				continue;
			}


			auto ics = _gics->get_nav_ics(*itsat);
			auto par_list = ics->get_par_list();

			// cycle par list
			auto it_beg = par_list.begin();
			auto it_end = par_list.end();
			for (auto par_iter = it_beg; par_iter != it_end; par_iter++)
			{
				t_gpar::t_type par_type = t_gics_nav::_icspar2gpar(par_iter->first);

				// init 
				t_gpar par_temp("", par_type, ++ipar_num, *itsat);
				par_temp.setTime(par_iter->second.beg_time, par_iter->second.end_time);
				par_temp.value(par_iter->second.value);
				par_temp.apriori(1E20);

				_lsqEstimator->add_parameter(par_temp);
			}
		}
	}

	void t_gorbfit::_set_out()
	{
		string tmp = dynamic_cast<t_gsetout*>(_set)->outputs("orbfit");
		if (tmp.empty())
		{
			tmp = "orbfit" + int2str(_beg.year()) + int2str(_beg.doy(), 3);
		}
		else
		{
			substitute(tmp, "$(date)", int2str(_beg.year()) + int2str(_beg.doy(), 3), false);
		}
		_orbdif = new t_giof;
		_orbdif->mask(tmp);
		_orbdif->append(false);

	}

	void t_gorbfit::_update_icsfile()
	{
		set<string> sats_ics = _gics->sat_list;
		set<string> sat_list = sats_ics;											   					  
		// update ics data
		for (auto sat_iter = sats_ics.begin(); sat_iter != sats_ics.end(); sat_iter++)
		{
			if (_sat.find(*sat_iter) == _sat.end())
			{
				_gics->del_nav_ics(*sat_iter);
				sat_list.erase(*sat_iter);						 
				continue;
			}
			auto ics = _gics->get_nav_ics(*sat_iter);
			auto& par_list = ics->get_par_list();

			// cycle par list
			for (auto par_iter = par_list.begin(); par_iter != par_list.end(); par_iter++)
			{
				t_gpar::t_type par_type = t_gics_nav::_icspar2gpar(par_iter->first);

				int par_idx = _lsqEstimator->_x_solve.getParam("", par_type, *sat_iter);

				ColumnVector dx = _lsqEstimator->dx();
				ColumnVector stdx = _lsqEstimator->stdx();

				par_iter->second.value += dx(par_idx + 1);
				par_iter->second.sigma = stdx(par_idx + 1);

			}
		}
		//updata ics sat_list
		_gics->sat_list = _sat;
		_gics->sat_num = _sat.size();
		//output new icsfile
		t_gtime beg  = dynamic_cast<t_gsetgen*>(_set)->beg();
		string  path = dynamic_cast<t_gsetout*>(_set)->outputs("ics");
		if (path.empty()) {
			path = "grt$(date).ics";
		}
		if (path.find("$(gwkd)") != string::npos)
		{
			substitute(path, "$(gwkd)", beg.str_gwkd(), false);
		}
		else if (path.find("$(date)") != string::npos)
		{
			substitute(path, "$(date)", beg.str_yyyydoy(), false);
		}

		t_gfile gout; t_ics_nav gics_coder(_set, "", 4096);
		gout.path(path); gout.glog(_glog);
		gics_coder.clear(); gics_coder.glog(_glog); gics_coder.path(path);
		gics_coder.add_data("ID1", _gics);
		gout.coder(&gics_coder);
		gout.run_write();
	}

}

