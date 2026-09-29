/**
*
* @verbatim
	History
	 -1.0 Yujie Qin		2019-12-12 creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gorbfitleo.cpp
* @brief		orbit fit processing for leos
*
* @author       Yujie Qin, Wuhan University
* @version		1.0.0
* @date		    2019-12-12
*
*/


#include <algorithm>
#include <set>
#include "gorbfitleo.h"

#include "gset/gsetinp.h"
#include "gset/gsetout.h" 
#include "gset/gsetorbfitleo.h"
#include "gutils/gtrs2crs.h"
#include "gutils/gcrs2acr.h"
#include "gutils/gleo.h"
#include "gio/gfile.h"
#include "gcoders/icsleo.h"
using namespace std;

namespace great {

	t_gorbfitleo::t_gorbfitleo()
	{
	}

	t_gorbfitleo::t_gorbfitleo(t_gsetbase * set, t_glog * log) :
		_lsqEstimator(nullptr)
	{
		_set = set;
		_glog = log;
	}

	t_gorbfitleo::~t_gorbfitleo()
	{

	}

	void t_gorbfitleo::setICS(t_gallics_leo* gics)
	{
		_gicsleo = gics;
	}
	void t_gorbfitleo::setORB(t_gallorb* gorb)
	{
		_gorb = gorb;
	}

	void t_gorbfitleo::setSP3(t_gallprec* gsp3)
	{
		_gsp3 = gsp3;
	}

	void t_gorbfitleo::setPOLEUT1(t_gpoleut1* gpoleut1)
	{
		_gpoleut1 = gpoleut1;
	}

	bool t_gorbfitleo::processBatch() {
		string function_name = "t_gorbfitleo::processBatch";
		if (!_init_data()) {
			_glog->comment(0, function_name, "ERROR : init data fail, check your xml file.");
			return false;
		}
		// get satellite list
		set<string> sat_sp3 = _gsp3->satellites();
		set<string> sat_orb = _gorb->getSatlist();
		int intv = _gsp3->intv();
		string producer = _gsp3->get_agency();
		for (auto iter = _sat.begin(); iter != _sat.end();) {
			if (sat_sp3.find(*iter) == sat_sp3.end()) {
				_glog->comment(1, function_name, "WARNING : " + *iter + "can not be found in sp3 file!");
				_sat.erase(iter++);
				continue;
			}
			if (sat_orb.find(*iter) == sat_orb.end()) {
				_glog->comment(1, function_name, "WARNING : " + *iter + "can not be found in orb(ics) file!");
				_sat.erase(iter++);
				continue;
			}
			if (_isfit) {
				set<string> sat_ics = _gicsleo->get_leo_sat_list();
				if (sat_ics.find(*iter) == sat_ics.end()) {
					_glog->comment(1, function_name, "WARNING : " + *iter + "can not be found in ics(orb) file!");
					_sat.erase(iter++);
					continue;
				}
			}
			++iter;
		}
		// loop by epoch
		ostringstream os;
		for (t_gtime epo_now = _beg; epo_now < _end; epo_now = epo_now + _sampling) {
			t_gtime epo_UTC = epo_now + 19.0;
			// loop by sat in xml
			for (auto iter = _sat.begin(); iter != _sat.end(); ++iter) {
				if (!_check_time(*iter)) return false;
				//get sat position in sp3
				double pos_trs_tmp[3];
				double vel_trs_tmp[3];
				int obsnum;
				if (_isinter)
				{
					if (!_isfit && producer != "GNUT" && (*iter == "JASON-2" || *iter == "JASON-3" || *iter == "HAIYANG-2")) {
						if (!_gsp3->posnew(*iter, epo_UTC, pos_trs_tmp)) continue;
					}
					else {
						if (!_gsp3->posnew(*iter, epo_now, pos_trs_tmp)) continue;
					}
				}
				else
				{
					if (!_isfit && producer != "GNUT" && (*iter == "JASON-2" || *iter == "JASON-3" || *iter == "HAIYANG-2")) {
						if (!_gsp3->get_pos_vel(*iter, epo_UTC, pos_trs_tmp, vel_trs_tmp, obsnum)) continue;
					}
					else {
						if (!_gsp3->get_pos_vel(*iter, epo_now, pos_trs_tmp, vel_trs_tmp, obsnum)) continue;
					}
				}
				
				// Check 0 or extreme numbers in sp3 records
				if (double_eq(pos_trs_tmp[0] * pos_trs_tmp[1] * pos_trs_tmp[2], 0.0) || pos_trs_tmp[0] >= 1E18) {
					// _sat.erase(iter++); // is it necessary?
					continue;
				}
				// get sat position in orbfile
				ColumnVector pos_crs_orb(0.0, 3);
				if (! _gorb->get_crs_sat_crd(epo_now, *iter, pos_crs_orb(1), pos_crs_orb(2), pos_crs_orb(3)))
					continue;
				//pos_crs_orb = pos_crs_orb / 1E3;
				os << " GPST" << setw(7) << right << epo_now.mjd()
					<< setw(7) << fixed << right << epo_now.sod() << "   ";
				ColumnVector pos_trs_sp3(0.0, 3);
				pos_trs_sp3(1) = pos_trs_tmp[0];
				pos_trs_sp3(2) = pos_trs_tmp[1];
				pos_trs_sp3(3) = pos_trs_tmp[2];
				// pos_trs_sp3 = pos_trs_sp3 / 1E3;
				
				// ChangeTimeSys(GPS2TDT, epo_now, epo_now_tdt);
				t_gtime epo_now_tdt = epo_now;
				epo_now_tdt.tsys(t_gtime::TT);
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


				//cout << "partial  " << setw(22) << setprecision(10) << partial << endl;

				t_glsqEquationMatrix equ;
				Matrix B; DiagonalMatrix P; ColumnVector l;
				int nPar = _lsqEstimator->_x_solve.parNumber();
				// resize the Matrix B l P
				B.ReSize(3, nPar); B = 0.0;
				l.ReSize(3);	   l = 0.0;
				P.ReSize(3);	   P = 1.0;

				l = pos_crs_sp3 - pos_crs_orb;
				Matrix partial;
				int site_index = _lsqEstimator->_x_solve.getParam("", t_gpar::PXSAT, *iter);
				if (site_index >= 0) {
					vector<t_gpar> par = _lsqEstimator->_x_solve.getOrbPar(*iter);
					partial.ReSize(par.size(), 3); partial = 0.0;
					_gorb->getLeoPartial(epo_now, *iter, partial, par);
				}
				for (int i = 1; i <= nPar; i++) {
					t_gpar par = _lsqEstimator->_x_solve[i-1];
					if (par.prn != "") {
						if (*iter == par.prn) {
							auto ics = _gicsleo->get_leo_ics(*iter);
							auto ics_list = ics->leoics();
							B.column(i) = partial.row(i - site_index).t() * 1E3;
							//cout << setw(25) << setprecision(12) << B.column(i);
						}
						continue;
					}
					switch (par.parType) {
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
					cor = B * _lsqEstimator->dx();
				l = l - cor;
				//cout << setw(15) << setprecision(8) << l << endl;
				equ.add_equ(B, P, l, vector<string>(3, *iter),
					vector<string>({ "X","Y","Z" }),
					vector<t_gobscombtype>(3, t_gobscombtype()));

				// output APRI to log 
				os << "APRI" << setw(15) << right << *iter;
				for (int i = 1; i <= 3; i++)
					os << setw(15) << fixed << setprecision(3) << right << l(i) * _unit;
				for (int i = 1; i <= 3; i++)
					os << setw(15) << fixed << setprecision(3) << right << cor(i) * _unit;
				os << endl;

				if (trs2crs_cord) { delete trs2crs_cord; trs2crs_cord = NULL; }
				_lsqEstimator->add_equation(equ, epo_now);
			} // loop: satellite 
		}     // loop: epoch
		//_lsqEstimator->write_equation(equ, epo_now);
		_lsqEstimator->solve_NEQ();
		//write the parameters value and correction to differ file
		os << setw(20) << left  << "Paramters Estimated"
		   << setw(20) << right << "Initial Value"
		   << setw(20) << right << "Correction"
		   << setw(20) << right << "Result" << endl;
		int npar = _lsqEstimator->_x_solve.parNumber();
		ColumnVector dx = _lsqEstimator->dx();
		for (int i = 0; i < npar; i++) {
			t_gpar par = _lsqEstimator->_x_solve[i];
			os  << setw(20) << left  << par.str_type()
				<< setw(20) << right << scientific << setprecision(10) << par.value() * _unit
				<< setw(20) << right << scientific << setprecision(10) << dx(i + 1) * _unit
				<< setw(20) << right << scientific << setprecision(10) << (par.value() + dx(i + 1)) * _unit << endl;
		}
		_glog->comment(1, function_name, "\n" + os.str());
		_recover_resid_leo();	
		// update ICS files
		if (_isfit) _update_icsfile();
		if (_lsqEstimator) { delete _lsqEstimator; _lsqEstimator = NULL;}
		return true;
	}

	void t_gorbfitleo::_read_resfile()
	{
		t_gallrecover recover_data;
		_lsqEstimator->recover_parameters(recover_data);

		string tmp;
		int i = 0;
		auto map_time_equ = recover_data.get_map_time_equ();

		for (auto iter : map_time_equ) {
			t_gtime epoch = iter.first;
			for (auto equ : iter.second)
			{
				int idx = 0;
				if (equ->sat_name == "X") {
					idx = 0;
				}
				else if (equ->sat_name == "Y") {
					idx = 1;
				}
				else if (equ->sat_name == "Z") {
					idx = 2;
				}
				else {
					throw("Recover equation is Wrong!");
				}
				_res[equ->site_name][epoch][idx] = equ->resuidal*_unit;

			}
		}
	}


	bool t_gorbfitleo::_check_time(string sat)
	{
		t_gtime beg_orb, end_orb;
		_gorb->getTime(sat, beg_orb, end_orb, _ref_orb);

		t_gtime beg_sp3 = _gsp3->beg_data(sat);
		t_gtime end_sp3 = _gsp3->end_data(sat);

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


	bool t_gorbfitleo::_init_data()
	{
		if (!_gicsleo || !_gorb || !_gsp3 || !_gpoleut1)
		{
			_glog->comment(0, "orbfitleo", "ERROR : orbfit ==> t_gorbfit::no comeplete file !");
			return false;
		}
		// get begin/end time,sampling, sat list and trans mode from xml
		if (!_get_setting())
		{
			_glog->comment(0, "orbfitleo", "ERROR : orbfit ==> t_gorbfit::there is something wrong witn xml !");
			return false;
		}
		if (_lsqEstimator) {
			delete _lsqEstimator;
		}
		_lsqEstimator = new t_glsq(NULL);
		_init_leo();
		
		_init_ERP_par();

		_set_out();
		return true;
	}

	bool t_gorbfitleo::_get_setting()
	{
		_beg = dynamic_cast<t_gsetgen*>(_set)->beg();
		_end = dynamic_cast<t_gsetgen*>(_set)->end();
		// get satlist from xml
		_sat = dynamic_cast<t_gsetorbfitleo*>(_set)->leo();

		if (_sat.empty())
		{
			_glog->comment(0, "orbdifleo", "ERROR : orbdif ==> No leos, check your xml file !");
			return false;
		}

		_isfit = dynamic_cast<t_gsetorbfitleo*>(_set)->isfit();
		_isinter= dynamic_cast<t_gsetorbfitleo*>(_set)->isinter();
		//_sampling = dynamic_cast<t_gsetgen*>(_set)->sampling();
		_sampling = _gsp3->intv(*_sat.begin());
		t_gtime beg_sp3 = _gsp3->beg_data(*_sat.begin()) + _sampling * 5;//to do yjqin
		t_gtime end_sp3 = _gsp3->end_data(*_sat.begin());
		/*if (_beg < beg_sp3)
			_beg = beg_sp3;
		if (_end > end_sp3)
			_end = end_sp3;*/
		if (abs(_sampling) < 5e-10)
		{
			_sampling = 30;
			//to do,yjqin
		}
		_end = _beg + ((int)((_end - _beg) / _sampling))*_sampling;
		// get trans from xml
		_trans = dynamic_cast<t_gsetorbfitleo*>(_set)->trans();
		_unit = dynamic_cast<t_gsetorbfitleo*>(_set)->unit();
		return true;
	}

	bool t_gorbfitleo::_init_leo()
	{
		// init lsq
		if (_isfit)
		{
			set<string> leo = _gicsleo->get_leo_sat_list();
			for (auto itsat = leo.begin(); itsat != leo.end(); ++itsat)
			{
				if (_sat.find(*itsat) == _sat.end())
				{
					_glog->comment(1, "orbfitleo:init_leo", "WARNING : " + *itsat + "is not in the sat list of xml");
					continue;
				}
				_icspar = 0;
				int ipar_num = _lsqEstimator->_x_solve.parNumber();
				vector<t_gpar> pars;
				

				auto ics = _gicsleo->get_leo_ics(*itsat);
				auto par_list = ics->leoics();

				for (auto iter_satpar = par_list->begin(); iter_satpar != par_list->end(); iter_satpar++)
				{
					t_gpar::t_type par_type = t_gics_leo::icspar2gpar(iter_satpar->first);
					vector<t_icsleo_record> onepar = iter_satpar->second;
					for (int k = 0; k < onepar.size(); k++)
					{
						//int indexPWC = 0;
						int npwc = onepar.size();
						t_gpar par_temp("", par_type, ++ipar_num, *itsat, false);
						par_temp.setTime(onepar[k].beg, onepar[k].end);
						par_temp.value(onepar[k].value);
						par_temp.apriori(1E20);

						if (onepar.size() > 1)  // PWC parameters
						{
							par_temp.nPWC = npwc;
							par_temp.lPWC = true;
							par_temp.piecewise_time = onepar[0].piecewise_time;
						}
						_lsqEstimator->add_parameter(par_temp);
						_icspar++;
						pars.push_back(par_temp);
					}
				}
				_lsqEstimator->_x_solve.setOrbPar(*itsat, pars);
			}
		}

		return true;
	}

	void t_gorbfitleo::_init_ERP_par()
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
		if (_trans.size() == 0)
		{
			t_gpar par_scale("", t_gpar::SCALE, ++ipar_num, "");
			par_scale.setTime(_beg, _end);
			par_scale.value(0);
			par_scale.apriori(1E-20);

			_lsqEstimator->add_parameter(par_scale);
		}
	}

	// void t_gorbfitleo::_init_SAT_par(string  satname) {}

	void t_gorbfitleo::_num2sstream(ostringstream& os, int num, int width) {
		if (abs(num) > 0.99*pow(10, width-1))
			os << setw(width) << setfill('*') << right << '*';
		else
			os << setw(width) << setfill(' ') << right << num;
	}

	void t_gorbfitleo::_recover_resid_leo() {
		_res.clear();
		_read_resfile();
		if (_res.size() == 0) {
			_glog->comment(0, "orbfitleo", "ERROR : orbdif ==> res record is empty!");
			return;
		}
		ostringstream os;
		os << setw(19) << right << "SAT";
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat)
			os << "        " << setw(10) << left << strl2strs(*itsat);
		os << endl;
		os << "      MJD       SOD";
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat)
			os << " -A/X- -C/Y- -R/Z-";
		os << endl;

		map<string, t_gtriple> max_acr, min_acr, fit_acr, mean_acr;
		map<string, double>	rms;
		double rmsmean = 0.0;
		_rms_statistics(max_acr, min_acr, fit_acr, mean_acr, rms, rmsmean);

		map<t_gtime, t_gtriple> tmp_for_epo = (*_res.begin()).second;
		auto it_epo = tmp_for_epo.begin();
		for (it_epo; it_epo != tmp_for_epo.end(); ++it_epo) {
			t_gtime epo_now = (*it_epo).first;
			os << "ACR" << setw(6) << setfill(' ') <<  right << epo_now.mjd()
			   << setw(10) << setfill(' ') << right << epo_now.sod();
			for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat) {
				_num2sstream(os, _res[*itsat][epo_now][0], 6);
				_num2sstream(os, _res[*itsat][epo_now][1], 6);
				_num2sstream(os, _res[*itsat][epo_now][2], 6);
			}
			os << endl;
		}
		os << setw(19) << left << "FITRMS(ACR)";
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat) {
			_num2sstream(os, fit_acr[*itsat][0], 6);
			_num2sstream(os, fit_acr[*itsat][1], 6);
			_num2sstream(os, fit_acr[*itsat][2], 6);
		}
		os << endl;
		os << setw(19) << left << "MAXRMS(ACR)";
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat) {
			_num2sstream(os, max_acr[*itsat][0], 6);
			_num2sstream(os, max_acr[*itsat][1], 6);
			_num2sstream(os, max_acr[*itsat][2], 6);
		}	
		os << endl;
		os << setw(19) << left << "MINRMS(ACR)";
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat) {
			_num2sstream(os, min_acr[*itsat][0], 6);
			_num2sstream(os, min_acr[*itsat][1], 6);
			_num2sstream(os, min_acr[*itsat][2], 6);
		}
		os << endl;
		os << setw(19) << left << "MEANRMS(ACR)";
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat) {
			_num2sstream(os, mean_acr[*itsat][0], 6);
			_num2sstream(os, mean_acr[*itsat][1], 6);
			_num2sstream(os, mean_acr[*itsat][2], 6);
		}
		os << endl;
		os << setw(19) << right << "SAT";
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat)
			os << "        " << setw(10) << left << strl2strs(*itsat);
		os << endl;
		double sigma0 = _lsqEstimator->sigma0() * 1E3*_unit;
		os << "Unit weight sigma0 = " << fixed << setprecision(1) << setw(10) << sigma0
		   << "   Unit = " << setprecision(3) << 1E3 / 1E6 / _unit << "m" << endl;
		os << "PRN mean";
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat)
			os << setw(5) << right << strl2strs(*itsat);
		os << endl;
		os << "RMS ";
		_num2sstream(os, rmsmean, 5);
		for (auto itsat = _sat.begin(); itsat != _sat.end(); ++itsat)
			_num2sstream(os, rms[*itsat], 5);
		os << endl;

		_orbdif->write(os.str().c_str(), os.str().size());
		_orbdif->flush();
		if (_orbdif) { delete _orbdif; _orbdif = NULL; }
	}

	void t_gorbfitleo::_set_out()
	{
		string tmp = dynamic_cast<t_gsetout*>(_set)->outputs("orbdif");
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

	void t_gorbfitleo::_update_icsfile()
	{
		// update
		for (auto sat_iter = _sat.begin(); sat_iter != _sat.end(); sat_iter++)
		{
			auto ics = _gicsleo->get_leo_ics(*sat_iter);
			map<ICS_PAR, vector<t_icsleo_record>>* par_list = ics->leoics();

			for (auto par_iter = par_list->begin(); par_iter != par_list->end(); par_iter++)
			{
				t_gpar::t_type par_type = t_gics_leo::icspar2gpar(par_iter->first);
				vector<t_icsleo_record> &pars = par_iter->second;
				for (int k = 0; k < pars.size(); k++)
				{
					t_gtime begtime = pars[k].beg;
					t_gtime endtime = pars[k].end;
					int par_idx = _lsqEstimator->_x_solve.getParam("", par_type, *sat_iter, begtime, endtime);
					if (par_idx >= 0)
					{
						ColumnVector dx = _lsqEstimator->dx();
						ColumnVector stdx = _lsqEstimator->stdx();

						// update value && sig
						pars[k].value += dx(par_idx + 1);
						pars[k].sig = stdx(par_idx + 1);

					}

				}

			}

		}

		//output new icsfile
		t_gtime beg  = dynamic_cast<t_gsetgen*>(_set)->beg();
		string  path = dynamic_cast<t_gsetout*>(_set)->outputs("icsleo");
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

		t_gfile gout; t_ics_leo gics_coder(_set, "", 4096);
		gout.path(path); gout.glog(_glog);
		gics_coder.clear(); gics_coder.glog(_glog); gics_coder.path(path);
		gics_coder.add_data("ID1", _gicsleo);
		gout.coder(&gics_coder);
		gout.run_write();
	}

	void t_gorbfitleo::_rms_statistics(map<string, t_gtriple>& max_acr, map<string, t_gtriple>& min_acr,
		map<string, t_gtriple>& fit_acr, map<string, t_gtriple>& mean_acr, map<string, double>& rms, double& rmsmean)
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
				mean_acr[sat][i] = 0.0;
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
					mean_acr[sat][i] = mean_acr[sat][i] + acr(i + 1);
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
					mean_acr[sat][i] = mean_acr[sat][i] / nepo;

				}
				rms[sat] = sqrt((pow(fit_acr[sat][0], 2) + pow(fit_acr[sat][1], 2) + pow(fit_acr[sat][2], 2)) / 3.0);
				rmsmean += rms[sat];
			}
		}
		rmsmean = rmsmean / nsat;
	}

}

