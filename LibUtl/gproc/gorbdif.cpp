/**
*
* @verbatim
	History
	 -1.0 Jie Li		2019-09-06  creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gorbdif.cpp
* @brief		orbit differ processing for GNS satellite
* @author       Jie Li, Wuhan University
* @version		1.0.0
* @date		    2019-09-06
*
*/


#include <algorithm>
#include <set>
#include "gorbdif.h"

#include "gset/gsetout.h" 
#include "gset/gsetorbdif.h"
#include "gutils/gtrs2crs.h"
#include "gutils/gcrs2acr.h"
using namespace std;

namespace great
{

	t_gorbdif::t_gorbdif()
	{
	}

	t_gorbdif::t_gorbdif(t_gsetbase* set, t_glog* log) :
		_lsqEstimator(nullptr)
	{
		_set = set;
		_glog = log;
	}

	t_gorbdif::~t_gorbdif()
	{

	}

	bool t_gorbdif::processBatch()
	{
		string function_name = "t_gorbdif::processBatch";

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
		t_gtime epo_now = _beg;
		ostringstream os;
		// output ERP parameters that needed tio be estimated
		/*os << "Parameters to be estimated:" << endl;
		for (int i = 0; i < _lsqEstimator->_x_solve.parNumber(); i++)
		{
			os << " " << setw(3) << i << setw(15) << right << _lsqEstimator->_x_solve[i].str_type()
				<< setw(15) << right << _beg.mjd(false) << setw(20) << setprecision(10)<<_beg.sod()
				<< setw(15) << right << _end.mjd(false) << setw(20) << setprecision(10) << right << _end.sod()
				<< setw(15) << right << setprecision(10) << 0.0
				<< setw(15) << right << setprecision(10) << _lsqEstimator->_x_solve[i].apriori() << endl;
		}*/

		while (epo_now < _end)
		{
			os << " GPST" << setw(15) << right << epo_now.mjd()
				<< setw(20) << fixed << right << epo_now.sod() << endl;

			t_glsqEquationMatrix equ;
			Matrix B; DiagonalMatrix P; ColumnVector l;
			vector<tuple<string, string, string> > equ_info;

			// loop by sat
			set<string> sat_sp3 = _gsp3->satellites();
			set<string> sat_orb = _gorb->getSatlist();
			for (auto iter = _sat.begin(); iter != _sat.end();)
			{
				// skip sat 
				if (sat_sp3.find(*iter) == sat_sp3.end())
				{
					_glog->comment(1, function_name, "WARNING : " + *iter + "can not be found in sp3 file!");
					_sat.erase(iter++);
					continue;
				}
				if (sat_orb.find(*iter) == sat_orb.end())
				{
					_glog->comment(1, function_name, "WARNING : " + *iter + "can not be found in orb(ics) file!");
					_sat.erase(iter++);
					continue;
				}
				
				//get sat position in sp3
				double pos_trs_tmp[3];
				_gsp3->pos(*iter, epo_now, pos_trs_tmp);
				if (double_eq(pos_trs_tmp[0], 0.0) && double_eq(pos_trs_tmp[1], 0.0) && double_eq(pos_trs_tmp[2], 0.0))
				{
					_gsp3->t_gallnav::pos(*iter, epo_now, pos_trs_tmp);
				}

				ColumnVector pos_trs_sp3(0.0, 3);
				pos_trs_sp3(1) = pos_trs_tmp[0];
				pos_trs_sp3(2) = pos_trs_tmp[1];
				pos_trs_sp3(3) = pos_trs_tmp[2];

				//sp3值检查
				// 0 值或大值
				if ((abs(pos_trs_sp3(1)) + abs(pos_trs_sp3(2)) + abs(pos_trs_sp3(3))) == 0 || pos_trs_sp3(1) >= 1E18)
				{
					_sat.erase(iter++);
					continue;
				}
				t_gtime epo_now_tdt= epo_now;
				epo_now_tdt.tsys(t_gtime::TT);
				//ChangeTimeSys(GPS2TDT, epo_now, epo_now_tdt);
				t_gtrs2crs* trs2crs_cord = new t_gtrs2crs(false, _gpoleut1);
				trs2crs_cord->calcRotMat(epo_now_tdt, true, true, true);
				Matrix rot_trs2crs = trs2crs_cord->getRotMat();

				double diff = (epo_now - _ref_orb) / 86400.0;
				Matrix dxdpole(0.0, 3, 4);
				dxdpole.Column(1) = (trs2crs_cord->getMatDx() * pos_trs_sp3) / RAD2SEC;
				dxdpole.Column(2) = dxdpole.Column(1) * diff;
				dxdpole.Column(3) = (trs2crs_cord->getMatDy() * pos_trs_sp3) / RAD2SEC;
				dxdpole.Column(4) = dxdpole.Column(3) * diff;
				Matrix dxdut1(0.0, 3, 2);
				dxdut1.Column(1) = (trs2crs_cord->getMatDu() * pos_trs_sp3) / RAD2TSEC;
				dxdut1.Column(2) = dxdut1.Column(1) * diff;

				ColumnVector pos_crs_sp3(0.0, 3);
				pos_crs_sp3 = rot_trs2crs * pos_trs_sp3;

				// get sat position in orbfile
				ColumnVector pos_crs_orb(0.0, 3);
				if (!_gorb->get_crs_sat_crd(epo_now, *iter, pos_crs_orb(1), pos_crs_orb(2), pos_crs_orb(3)))
				{
					_glog->comment(1, function_name, "WARNING : " + *iter + "can not be get xyz from orb file!");
					iter++;
					continue;
				}

				int nPar = _lsqEstimator->_x_solve.parNumber();
				// resize the Matrix B l P
				B.ReSize(3, nPar); B = 0.0;
				l.ReSize(3);	l = 0.0;
				P.ReSize(3);	P = 1.0;

				l = pos_crs_sp3 - pos_crs_orb;

				for (int i = 1; i <= nPar; i++)
				{
					t_gpar par = _lsqEstimator->_x_solve[i - 1];
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

				if (_excsat.find(*iter) != _excsat.end() ||
					_excsys.find(t_gsys::gsys2str(t_gsys::sat2gsys(*iter))) != _excsys.end())
				{
					t_glsqEquationMatrix equ_tmp;
					equ_tmp.add_equ(B, P, l, vector<string>({ "X","Y","Z" }),
						vector<string>(3, *iter),
						vector<t_gobscombtype>(3, t_gobscombtype()));
					_lsqEstimator->write_equation(equ_tmp, epo_now);
					iter++;
					continue;
				}

				equ.add_equ(B, P, l, vector<string>({ "X","Y","Z" }),
					vector<string>(3, *iter),
					vector<t_gobscombtype>(3, t_gobscombtype()));

				if (trs2crs_cord)
				{
					delete trs2crs_cord;
					trs2crs_cord = NULL;
				}

				iter++;
			}
			/*if (B.is_zero())
			{
				_glog->comment(1, function_name, "WARNING : " + epo_now.str_mjdsod("no useful equations!"));
				epo_now = epo_now + _sampling;
				continue;
			}*/
			
			_lsqEstimator->add_equation(equ, epo_now);
			epo_now = epo_now + _sampling;
		}

		_lsqEstimator->solve_NEQ();

		/*os << "Parameters estimated result:" << endl;
		for (int i = 0; i < _lsqEstimator->_x_solve.parNumber(); i++)
		{
			os << " " << setw(3) << i << setw(15) << right << _lsqEstimator->_x_solve[i].str_type()
				<< setw(15) << right << _beg.mjd(false) << setw(20) << right << setprecision(10) << _beg.sod()
				<< setw(15) << right << _end.mjd(false) << setw(20) << right << setprecision(10) << _end.sod()
				<< setw(20) << right << setw(15) << setprecision(10) << _lsqEstimator->dx()(i + 1) << endl;
		}*/
		_orbdif->open(_orbdif->mask(), ios::out);
		_orbdif->write(os.str().c_str(), os.str().size());
		_orbdif->flush();

		// 残差统计并写入文件
		_recover_resid();

		if (_lsqEstimator)
		{
			delete _lsqEstimator;
			_lsqEstimator = NULL;
		}

		return true;
	}

	void t_gorbdif::setICS(t_gallics_nav* gics)
	{
		_gics = gics;
	}

	void t_gorbdif::setORB(t_gallorb* gorb)
	{
		_gorb = gorb;
	}

	void t_gorbdif::setSP3(t_gallnav* gsp3)
	{
		_gsp3 = gsp3;
	}

	void t_gorbdif::setPOLEUT1(t_gpoleut1* gpoleut1)
	{
		_gpoleut1 = gpoleut1;
	}

	bool t_gorbdif::_init_data()
	{
		if (!_gorb || !_gsp3 || !_gpoleut1)
		{
			_glog->comment(0, "orbdif", "ERROR : orbdif ==> t_gorbdif::no comeplete file !");
			return false;
		}
		// get begin/end time,sampling, sat list and trans mode from xml
		if (!_get_setting())
		{
			_glog->comment(0, "orbdif", "ERROR : orbdif ==> t_gorbdif::there is something wrong witn xml !");
			return false;
		}

		// init lsq
		if (_lsqEstimator) {
			delete _lsqEstimator;
		}
		_lsqEstimator = new t_glsq(NULL);

		_init_par();

		_set_out();

		return true;
	}

	bool t_gorbdif::_get_setting()
	{
		_beg = dynamic_cast<t_gsetgen*>(_set)->beg();
		_end = dynamic_cast<t_gsetgen*>(_set)->end();
		_sampling = dynamic_cast<t_gsetgen*>(_set)->sampling();

		_end = _beg + ((int)((_end - _beg) / _sampling)) * _sampling;

		// get satlist from xml
		set<string> sys = dynamic_cast<t_gsetgen*>(_set)->sys();
		_excsys = dynamic_cast<t_gsetorbdif*>(_set)->excsys();
		/*set<string> tmp;
		set_difference(sys.begin(), sys.end(), excsys.begin(), excsys.end(), inserter(tmp, tmp.begin()));
		sys = tmp;
		tmp.clear();
		if (sys.empty())
		{
			_glog->comment(0, "orbdif", "ERROR : orbdif ==> No satellite systerm, check your xml file(sys and excsys) !");
			return false;
		}*/
		// sats list from gset
		_sat = dynamic_cast<t_gsetgnss*>(_set)->sat();
		if (_sat.empty())
		{
			t_map_sats gnss_sats = GNSS_SATS();
			for (auto itGNS = sys.begin(); itGNS != sys.end(); ++itGNS) {
				GSYS gsys = t_gsys::str2gsys(*itGNS);
				_sat.insert(gnss_sats[gsys].begin(), gnss_sats[gsys].end());
			}
		}
		_excsat = dynamic_cast<t_gsetorbdif*>(_set)->excsat();
		/*set_difference(_sat.begin(), _sat.end(), excsat.begin(), excsat.end(), inserter(tmp, tmp.begin()));
		_sat = tmp;
		tmp.clear();*/

		// get trans from xml
		_trans = dynamic_cast<t_gsetorbdif*>(_set)->trans();
		return true;
	}

	// added by yqyuan for orbfit
	bool t_gorbdif::_get_setting_orbfit()
	{
		_beg = dynamic_cast<t_gsetgen*>(_set)->beg();
		_end = dynamic_cast<t_gsetgen*>(_set)->end();
		_sampling = dynamic_cast<t_gsetgen*>(_set)->sampling();

		_end = _beg + ((int)((_end - _beg) / _sampling)) * _sampling;

		// get satlist from xml
		set<string> sys = dynamic_cast<t_gsetgen*>(_set)->sys();
		set<string> excsys = dynamic_cast<t_gsetorbdif*>(_set)->excsys_orbfit();
		set<string> tmp;
		set_difference(sys.begin(), sys.end(), excsys.begin(), excsys.end(), inserter(tmp, tmp.begin()));
		sys = tmp;
		tmp.clear();
		if (sys.empty())
		{
			_glog->comment(0, "orbfit", "ERROR : orbfit ==> No satellite systerm, check your xml file !");
			return false;
		}
		// sats list from gset
		_sat = dynamic_cast<t_gsetgnss*>(_set)->sat();
		if (_sat.empty())
		{
			t_map_sats gnss_sats = GNSS_SATS();
			for (auto itGNS = sys.begin(); itGNS != sys.end(); ++itGNS) {
				GSYS gsys = t_gsys::str2gsys(*itGNS);
				_sat.insert(gnss_sats[gsys].begin(), gnss_sats[gsys].end());
			}
		}
		set<string> excsat = dynamic_cast<t_gsetorbdif*>(_set)->excsat_orbfit();
		set_difference(_sat.begin(), _sat.end(), excsat.begin(), excsat.end(), inserter(tmp, tmp.begin()));
		_sat = tmp;
		tmp.clear();

		// get trans from xml
		_trans = dynamic_cast<t_gsetorbdif*>(_set)->trans_orbfit();
		return true;
	}
	bool t_gorbdif::_check_time(string sat)
	{
		t_gtime beg_orb, end_orb;
		_gorb->getTime(sat, beg_orb, end_orb, _ref_orb);

		t_gtime beg_sp3 = _gsp3->beg_time(sat);
		t_gtime end_sp3 = _gsp3->end_time(sat);
		//rinexn 文件无数据，读sp3文件
		if (beg_sp3 == LAST_TIME || end_sp3 == FIRST_TIME)
		{
			beg_sp3 = dynamic_cast<t_gallprec*>(_gsp3)->beg_data(sat);
			end_sp3 = dynamic_cast<t_gallprec*>(_gsp3)->end_data(sat);
		}

		// jdhuang：修改条件判断阈值
		/*if (double_eq(fabs(_beg - beg_orb),3600) || 
			double_eq(fabs(_beg - beg_sp3),3600) || 
			double_eq(fabs(_end - end_orb),3600) || 
			double_eq(fabs(_end - end_sp3),3600))
		{
			_glog->comment(0, "orbdif", "ERROR : orbdif ==> time is otrange, please chelk time in xml file or sp3 file or orb file  !");
			return false;
		}*/
		
		if (_beg < beg_orb) _beg = beg_orb;
		if (_beg < beg_sp3) _beg = beg_sp3;
		if (_end > end_orb) _end = end_orb;
		if (_end > end_sp3) _end = end_sp3;
		if (_beg >= _end)
		{
			_glog->comment(0, "orbdif", "ERROR : orbdif ==> time is otrange, please chelk time in xml file or sp3 file or orb file  !");

			return false;
		}

		return true;
	}

	void t_gorbdif::_init_par()
	{
		//init erp par
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
		if (ipar_num == 0)
		{
			t_gpar par_scale("", t_gpar::SCALE, ++ipar_num, "");
			par_scale.setTime(_beg, _end);
			par_scale.value(0.0);
			par_scale.apriori(1E-20);

			_lsqEstimator->add_parameter(par_scale);
		}
	}

	void t_gorbdif::_set_out()
	{
		string tmp = dynamic_cast<t_gsetout*>(_set)->outputs("orbdif");

		if (tmp.empty())
		{
			string tmp = dynamic_cast<t_gsetout*>(_set)->outputs("dif");
		}

		if (tmp.empty())
		{
			tmp = "orbdif" + int2str(_beg.year()) + int2str(_beg.doy(), 3);
		}
		else
		{
			substitute(tmp, "$(date)", int2str(_beg.year()) + int2str(_beg.doy(), 3), false);
		}
		_orbdif = new t_giof;
		_orbdif->mask(tmp);
		_orbdif->append(false);

	}

	void t_gorbdif::_recover_resid()
	{
		_read_resfile();
		if (_res.size() == 0)
		{
			_glog->comment(0, "orbdif", "ERROR : orbdif ==> res record is empty!");
			return;
		}
		ostringstream os;
		os << setw(19) << right << "SAT";
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat)
		{
			os << "        " << *itsat << "       ";
		}
		os << endl;
		os << "      MJD       SOD";
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat)
		{
			os << " -A/X- -C/Y- -R/Z-";
		}
		os << endl;

		map<string, t_gtriple> max_acr;
		map<string, t_gtriple> min_acr;
		map<string, t_gtriple> fit_acr;
		map<string, double>	rms;
		double rmsmean = 0.0;
		_rms_statistics(max_acr, min_acr, fit_acr, rms, rmsmean);

		map<t_gtime, t_gtriple> tmp_for_epo = (*_res.begin()).second;
		auto it_epo = tmp_for_epo.begin();
		for (it_epo; it_epo != tmp_for_epo.end(); ++it_epo)
		{
			t_gtime epo_now = (*it_epo).first;
			os << "ACR" << setw(6) << right << epo_now.mjd()
				<< setw(9) << right << epo_now.sod();
			for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat)
			{
				os << setw(6) << right << (int)_res[(*itsat)][epo_now][0]
					<< setw(6) << (int)_res[(*itsat)][epo_now][1]
					<< setw(6) << (int)_res[(*itsat)][epo_now][2];
			}
			os << endl;
		}
		os << setw(18) << left << "FITRMS(ACR)";
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat)
		{
			os << setw(6) << right << (int)fit_acr[*(itsat)][0]
				<< setw(6) << (int)fit_acr[*(itsat)][1]
				<< setw(6) << (int)fit_acr[*(itsat)][2];
		}
		os << endl;
		os << setw(18) << left << "MAXRMS(ACR)";
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat)
		{
			os << setw(6) << right << (int)max_acr[*(itsat)][0]
				<< setw(6) << (int)max_acr[*(itsat)][1]
				<< setw(6) << (int)max_acr[*(itsat)][2];
		}
		os << endl;
		os << setw(18) << left << "MINRMS(ACR)";
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat)
		{
			os << setw(6) << right << (int)min_acr[*(itsat)][0]
				<< setw(6) << (int)min_acr[*(itsat)][1]
				<< setw(6) << (int)min_acr[*(itsat)][2];
		}
		os << endl;
		os << setw(18) << right << "SAT";
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat)
		{
			os << "        " << *itsat << "       ";
		}
		os << endl;
		double sigma0 = _lsqEstimator->sigma0() * 1E3;
		os << "Unit weight sigma0 = " << fixed << setprecision(1) << setw(10) << sigma0
			<< "   Unit = " << setprecision(3) << 1E3 / 1E6 << "m" << endl;

		os << "PRN mean";
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat)
		{
			os << setw(5) << right << *itsat;
		}
		os << endl;
		os << "RMS " << setw(4) << left << (int)rmsmean;
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat)
		{
			os << setw(5) << right << (int)rms[*itsat];
		}
		os << endl;

		_orbdif->write(os.str().c_str(), os.str().size());
		_orbdif->flush();
		if (_orbdif)
		{
			delete _orbdif;
			_orbdif = NULL;
		}
	}

	void t_gorbdif::_read_resfile()
	{

		t_gallrecover recover_data;
		_lsqEstimator->recover_parameters(recover_data);

		string tmp;
		int i = 0;
		auto map_time_equ = recover_data.get_map_time_equ();
		
		for (auto iter : map_time_equ) {
			t_gtime epoch = iter.first;
			for (auto equ: iter.second)
			{
				int idx = 0;
				if (equ->site_name == "X") {
					idx = 0;
				}
				else if (equ->site_name == "Y") {
					idx = 1;
				}
				else if (equ->site_name == "Z") {
					idx = 2;
				}
				else {
					throw("Recover equation is Wrong!");
				}
				_res[equ->sat_name][epoch][idx] = equ->resuidal;

			}
		}
	}

	void t_gorbdif::_rms_statistics(map<string, t_gtriple>& max_acr, map<string, t_gtriple>& min_acr,
		map<string, t_gtriple>& fit_acr, map<string, double>& rms, double& rmsmean)
	{

		rmsmean = 0.0;
		int nsat = 0;
		for (auto itres = _res.begin(); itres != _res.end(); ++itres)
		{
			string sat = (*itres).first;
			for (int i = 0; i < 3; i++)
			{
				max_acr[sat][i] = -1E10;
				min_acr[sat][i] = 1E10;
				fit_acr[sat][i] = 0.0;
			}
			rms[sat] = 0.0;
			map<t_gtime, t_gtriple > tmp_res = (*itres).second;
			double nepo = 0;
			for (auto it_t = tmp_res.begin(); it_t != tmp_res.end(); ++it_t)
			{
				nepo++;
				t_gcrs2acr* crs2acr_cord = new t_gcrs2acr();
				ColumnVector sat_pos(0.0, 3);
				ColumnVector sat_vel(0.0, 3);
				_gorb->get_crs_sat_crd((*it_t).first, sat, sat_pos(1), sat_pos(2), sat_pos(3));
				_gorb->get_crs_sat_vel((*it_t).first, sat, sat_vel(1), sat_vel(2), sat_vel(3));
				Matrix rotmat(0.0, 3, 3);
				crs2acr_cord->calcRotMat(sat_pos, sat_vel, rotmat);
				ColumnVector xyz = _res[sat][(*it_t).first].crd_cvect();
				ColumnVector acr = rotmat * xyz * 1E3;

				_res[sat][(*it_t).first][0] = acr(1);
				_res[sat][(*it_t).first][1] = acr(2);
				_res[sat][(*it_t).first][2] = acr(3);
				for (int i = 0; i < 3; i++)
				{
					max_acr[sat][i] = max(acr(i + 1), max_acr[sat][i]);
					min_acr[sat][i] = min(acr(i + 1), min_acr[sat][i]);
					fit_acr[sat][i] = fit_acr[sat][i] + acr(i + 1) * acr(i + 1);
				}
				if (crs2acr_cord)
				{
					delete crs2acr_cord;
					crs2acr_cord = NULL;
				}
			}
			if (nepo != 0)
			{
				nsat++;
				for (int i = 0; i < 3; i++)
				{
					fit_acr[sat][i] = sqrt(fit_acr[sat][i] / nepo);

				}
				rms[sat] = sqrt((pow(fit_acr[sat][0], 2) + pow(fit_acr[sat][1], 2) + pow(fit_acr[sat][2], 2)) / 3.0);
				rmsmean += rms[sat];
			}
		}
		rmsmean = rmsmean / nsat;
	}
}

