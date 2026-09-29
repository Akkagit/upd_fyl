/**
*
* @verbatim
	History
	 -1.0 Xiongyun  2019-12-10  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gpceprecisemodel.cpp
* @brief		precise model for computing  correction
*
* @author       Xiongyun, Wuhan University
* @version		1.0.0
* @date		    2019-12-10
*
*/

#include "gmodels/gpceprecisemodel.h"
#include "gmodels/gprecisemodel.h"
#include "gutils/ginfolog.h"

namespace great
{
	t_gpceprecisemodel::t_gpceprecisemodel(t_gallproc* data, t_glog* log, t_gsetbase* setting) :
		t_gprecisemodel(data, log, setting)
	{
		_carrier_range = dynamic_cast<t_gsetproc*>(setting)->carrier_range();
		_IRC = dynamic_cast<t_gsetproc*>(setting)->IRC();
	}

	t_gpceprecisemodel::~t_gpceprecisemodel()
	{

	}

	bool t_gpceprecisemodel::_cmb_equ(const t_gtime& epoch, t_gsatdata& crt_obs, t_gallpar& pars, t_glsqEquationMatrix& equ)
	{
		_isObs13Valid = false;
		auto gsys = crt_obs.gsys();
		bool obs_valid = _check_obs_valid(crt_obs, _observ, _band_index[gsys][FREQ_1], _band_index[gsys][FREQ_2]);
		if (!obs_valid)
		{
			write_log_info(_log, 1, "NOTE", "obs is not valid" + epoch.str_ymdhms("", false, false));
			return obs_valid;
		}

		bool update_valid = _update_obs_info(epoch, crt_obs, pars);
		if (!update_valid)
		{
			write_log_info(_log, 1, "NOTE", "update obs information failed" + epoch.str_ymdhms("", false, false));
			return update_valid;
		}

		// at fisrt prepare the obs
		bool pre_valid = _prepare_obs(epoch, pars);
		if (!pre_valid)
		{
			write_log_info(_log, 1, "NOTE", "prepare obs information failed" + epoch.str_ymdhms("", false, false));
			return pre_valid;
		}

		bool cmb_valid = false;
		switch (_observ)
		{
		case IONO_FREE:
			cmb_valid = _cmb_equ_IF(epoch, pars, equ);
			break;
		case RAW_ALL:
			cmb_valid = _cmb_equ_ALL(epoch, pars, equ);
			break;
		default:
			return false;
		}
		if (!cmb_valid)
		{
			write_log_info(_log, 1, "NOTE", "cmb failed");
			return false;
		}

		t_gprecisemodel::_update_obs_info(crt_obs);
		return true;
	}

	bool t_gpceprecisemodel::_update_obs_info(const t_gtime& epoch, t_gsatdata& obsdata, t_gallpar& pars)
	{
		if (!_gall_nav || !_gallobj)
		{
			write_log_info(_log, 1, "NOTE", "no pars");
			return false;
		}

		// Some frequently used variables are also defined here.
		_crt_epo = epoch;
		_crt_obs = obsdata;
		_crt_obs.tb12(obsdata.tb12()); // 因为拷贝构造函数无法拷贝这两个值，暂时如此处理
		_crt_obs.tb13(obsdata.tb13()); // 因为拷贝构造函数无法拷贝这两个值，暂时如此处理
		_crt_sat = obsdata.sat();
		_crt_rec = obsdata.site();
		_crt_sys = obsdata.gsys();
		_crt_obj = _gallobj->obj(_crt_rec);
		_crt_obj_type = _crt_obj->id_type();
		_site = _crt_rec;
		_grec = _crt_obj;

		bool rec_clk_valid = _update_obj_clk("rec" + _crt_rec, epoch, pars, _crt_rec_clk, _obj_clk, _crt_sys);
		bool sat_clk_valid = _update_obj_clk("sat" + _crt_sat, epoch, pars, _crt_sat_clk, _obj_clk, _crt_sys);
		bool sat_ifcb_valid = _update_obj_clk("sat" + _crt_sat + "_ifcb", epoch, pars, _crt_sat_ifcb, _obj_clk, _crt_sys);
		if (!rec_clk_valid || !sat_clk_valid)
		{
			write_log_info(_log, 1, "NOTE", "no rec or sat clk for epoch " + epoch.str_ymdhms());
			return false;
		}
		return true;
	}

	bool t_gpceprecisemodel::_update_obj_clk(const string& obj, const t_gtime& crt_epoch,
		t_gallpar& pars, double& clk, map<string, pair<t_gtime, double> >& obj_clk, GSYS& ctr_sys)
	{
		if (!_gall_nav || obj.substr(3).empty())
		{
			write_log_info(_log, 1, "NOTE", "no navigation files or obj is NONE");
			return false;
		}

		// clk[s] from clk files
		double clk_rms = 0.0;
		double dclk = 0.0;
		string type = obj.substr(0, 3);
		string name = obj.substr(3);
		// update clk
		if (obj_clk[name].first != crt_epoch)
		{
			int clk_valid = _gall_nav->clk(name, crt_epoch, &clk, &clk_rms, &dclk);
			if (clk_valid < 0) clk = 0.0;
			if (name.size() <= 4)
			{
				obj_clk[name].first = crt_epoch;
				obj_clk[name].second = clk;
			}
			else
			{
				int idx = pars.getParam("", t_gpar::IFCB_F3, name.substr(0, 3));
				obj_clk[name].first = crt_epoch;
				obj_clk[name].second += pars[idx].value();
			}
		}
		else
		{
			clk = obj_clk[name].second;
		}
		// get clk from clk files
		//--can not get any clk informtion from file or gallpar
		int idx_clk = -1;
		if (type == "rec")
		{
			if (_IRC)
			{
				t_gpar::t_type recclk_type;
				switch (ctr_sys)
				{
				case GSYS::GPS:
					recclk_type = t_gpar::CLK_G;
					break;
				case GSYS::GLO:
					recclk_type = t_gpar::CLK_R;
					break;
				case GSYS::BDS:
					recclk_type = t_gpar::CLK_C;
					break;
				case GSYS::GAL:
					recclk_type = t_gpar::CLK_E;
					break;
				case GSYS::QZS:
					recclk_type = t_gpar::CLK_J;
					break;
				default:
					write_log_info(_log, 1, "NOTE", "unknown satellite system(GRECJ) ");
					return false;
				}
				idx_clk = pars.getParam(name, recclk_type, "");
			}
			else
			{
				idx_clk = pars.getParam(name, t_gpar::CLK, "");
			}
		}
		else if (type == "sat" && name.size() == 3) idx_clk = pars.getParam("", t_gpar::CLK_SAT, name);
		else if (type == "sat" && name.size() > 3) idx_clk = pars.getParam("", t_gpar::IFCB_F3, name.substr(0, 3));
		else idx_clk = -1;
		if (double_eq(clk, 0.0) && idx_clk < 0)
		{
			write_log_info(_log, 1, "NOTE", "no rec clk for " + crt_epoch.str_ymdhms());
			return false;
		}
		//-- no clk file but clk par exist
		if (double_eq(clk, 0.0) && idx_clk >= 0)
		{
			//clk = pars[idx_clk].value() / CLIGHT;
			if (_IRC && type == "rec" && _rec_clk[name])
			{
				clk = _rec_clk[name];
			}
			else
			{
				clk = pars[idx_clk].value() / CLIGHT;
			}
		}
		//-- update clk par
		if (idx_clk >= 0)
		{
			pars[idx_clk].value(clk * CLIGHT);
			if (_frequency == 3 && type == "rec")
			{
				idx_clk = pars.getParam(name, t_gpar::G15CLK, "");
				if (idx_clk >= 0) pars[idx_clk].value(clk * CLIGHT);
				idx_clk = pars.getParam(name, t_gpar::E15CLK, "");
				if (idx_clk >= 0) pars[idx_clk].value(clk * CLIGHT);
				idx_clk = pars.getParam(name, t_gpar::C15CLK, "");
				if (idx_clk >= 0) pars[idx_clk].value(clk * CLIGHT);
				idx_clk = pars.getParam(name, t_gpar::J15CLK, "");
				if (idx_clk >= 0) pars[idx_clk].value(clk * CLIGHT);
			}
		}

		return true;
	}

	bool t_gpceprecisemodel::_cmb_equ_IF(const t_gtime& epoch, t_gallpar& pars, t_glsqEquationMatrix& equ)
	{
		vector<pair<int, double> > coefP, coefP13;
		vector<pair<int, double> > coefL, coefL13;
		GOBSBAND b1 = t_gsys::band_priority(_crt_sys, FREQ_1);
		GOBSBAND b2 = t_gsys::band_priority(_crt_sys, FREQ_2);
		GOBSBAND b3 = t_gsys::band_priority(_crt_sys, FREQ_3);
		t_gobs obsP = _crt_obs.select_range(b1);
		t_gobs obsL = _crt_obs.select_phase(b1);
		t_gobs obsP2 = _crt_obs.select_range(b2);
		t_gobs obsL2 = _crt_obs.select_phase(b2);
		t_gobs obsP3 = _crt_obs.select_range(b3);
		t_gobs obsL3 = _crt_obs.select_phase(b3);

		double omcP = 0.0;
		double omcL = 0.0;
		double wgtP = 0.0;
		double wgtL = 0.0;

		_isTriple = false;
		if (!_omc_obs_IF(epoch, obsP, pars, omcP) ||
			!_wgt_obs(obsP, _factorP, wgtP) ||
			//!_prt_obs(epoch, obsP, pars, coefP) merge
			!_prt_obs(epoch, pars, obsP2, coefP))
		{
			return false;
		}
		//if (_carrier_range)
		//{
		//	equ.add_equ(coefP, wgtP, omcP, _crt_rec, _crt_sat, t_gobscombtype(obsP, _observ));
		//}
		//else
		//{
		//	equ.add_equ(coefP, wgtP, omcP, _crt_rec, _crt_sat, t_gobscombtype(obsP, _observ));
		//}
		t_gobscombtype typeP(obsP, _observ);
		equ.add_equ(coefP, wgtP, omcP, _crt_rec, _crt_sat, typeP, false);
		if (_frequency == 3 && _isObs13Valid)
		{
			_isTriple = true;
			if (!_omc_obs_IF13(epoch, obsP3, pars, omcP) ||
				!_wgt_obs(obsP3, _factorP, wgtP, true) ||
				//!_prt_obs(epoch, obsP3, pars, coefP13)
				!_prt_obs(epoch, pars, obsP3, coefP13))
			{
				return false;
			}
			t_gobscombtype typeP3(obsP3, _observ);
			equ.add_equ(coefP13, wgtP, omcP, _crt_rec, _crt_sat, typeP3,false);
		}

		_isTriple = false;
		if (!_omc_obs_IF(epoch, obsL, pars, omcL) ||
			!_wgt_obs(obsL, _factorP, wgtL) ||
			//!_prt_obs(epoch, obsL, pars, coefL)
			!_prt_obs(epoch, pars, obsL2, coefL))
		{
			equ.remove_last_equation();
			return false;
		}

		// set new amb flag 
		if (!_carrier_range)
		{
			bool is_newamb = false;
			int idx = pars.getParam(_crt_rec, t_gpar::AMB_IF, _crt_sat);
			if (pars[idx].beg == epoch) {
				is_newamb = true;
			}

			t_gobscombtype type(obsL, _observ);
			equ.add_equ(coefL, wgtL, omcL, _crt_rec, _crt_sat, type, is_newamb);
		}
		else
		{
			t_gobscombtype type(obsL, _observ);
			equ.add_equ(coefL, wgtL, omcL, _crt_rec, _crt_sat, type, false);
		}

		if (_frequency == 3 && _isObs13Valid)
		{
			_isTriple = true;
			if (!_omc_obs_IF13(epoch, obsL3, pars, omcL) ||
				!_wgt_obs(obsL3, _factorP, wgtL, true) ||
				//!_prt_obs(epoch, obsL3, pars, coefL13)
				!_prt_obs(epoch, pars, obsL3, coefL))
			{
				equ.remove_last_equation();
				return false;
			}
			bool is_newamb = false;
			int idx = pars.getParam(_crt_rec, t_gpar::AMB13_IF, _crt_sat);
			if (pars[idx].beg == epoch) {
				is_newamb = true;
			}
			equ.add_equ(coefL13, wgtL, omcL, _crt_rec, _crt_sat, t_gobscombtype(obsL3, _observ), is_newamb);
		}

		//P(obs_num) = wgtL; l(obs_num) = omcL;
#ifdef DEBUG_LOG
		string tmp2;
		string time2(_crt_epo.str_ymdhms("", false, false));
		char tmp_num2[8];   sprintf(tmp_num2, "%8d", obs_num);   string num2(tmp_num2);
		char tmp_wgt2[18];  sprintf(tmp_wgt2, "%16.8f", wgtL);   string wgt2(tmp_wgt2);
		char tmp_rsi2[18];  sprintf(tmp_rsi2, "%16.8f", omcL);   string rsi2(tmp_rsi2);
		tmp2 = time2 +
			" recL : " + _crt_rec +
			" satL : " + _crt_sat +
			" numL : " + num2 +
			" wgtL : " + wgt2 +
			" rsiL : " + rsi2;
#endif 

#ifdef DEBUG_LOG
		write_log_info(_debug, 0, "DEBUG", tmp2);
#endif // DEBUG_LOG
		return true;
	}

	bool t_gpceprecisemodel::_cmb_equ_ALL(const t_gtime& epoch, t_gallpar& pars, t_glsqEquationMatrix& equ)
	{
		// first get the sys and the bands
		//std::vector<GOBSBAND> crt_bands;
		//if (_SYS_bands.count(_crt_sys) > 0) crt_bands = _SYS_bands[_crt_sys];
		std::map<FREQ_SEQ, GOBSBAND> crt_bands;
		if (_band_index.count(_crt_sys) > 0) crt_bands = _band_index[_crt_sys];
		if (crt_bands.empty())
		{
			write_log_info(_log, 1, "NOTE", "no band for " + t_gsys::gsys2str(_crt_sys));
			return false;
		}
		_isTriple = false;
		if (_frequency == 3 && _isObs13Valid)
		{
			_isTriple = true;
		}

		for (const auto& iter : crt_bands)
		{
			vector<pair<int, double> > coefP;
			vector<pair<int, double> > coefL;

			double omcP = 0.0;
			double omcL = 0.0;
			double wgtP = 0.0;
			double wgtL = 0.0;
			//GOBSBAND band = iter; merge
			GOBSBAND band = iter.second;
			t_gobs  obsP = _crt_obs.select_range(band, true);
			t_gobs  obsL = _crt_obs.select_phase(band, true);

			// jdhuang add for freq3 check
			auto freq = t_gsys::band2freq(_crt_sys, band);
			if (_frequency == 2 && freq >= FREQ_3) continue; // 双频模式，只处理freq 1 和 freq 2
			if (_frequency == 3 && freq == FREQ_2 && !_crt_obs.tb12()) continue;
			if (_frequency == 3 && freq == FREQ_3 && !_crt_obs.tb13()) continue;

			if (!_omc_obs_ALL(epoch, obsP, pars, omcP))
			{
				write_log_info(_log, 1, "NOTE", "omc range obs failed");
				return false;
			};

			if (!_wgt_obs(obsP, _factorP, wgtP))
			{
				write_log_info(_log, 1, "NOTE", "weight range obs failed");
				return false;
			};

			//if (!_prt_obs(epoch, obsP, pars, coefP)) merge
			if (!_prt_obs(epoch, pars, obsP, coefP))
			{
				write_log_info(_log, 1, "NOTE", "partial range obs failed");
				return false;
			};
			t_gobscombtype typeP(obsP, _observ);
			equ.add_equ(coefP, wgtP, omcP, _crt_rec, _crt_sat, typeP, false);

#ifdef DEBUG_LOG

			string tmp;
			string time(_crt_epo.str_ymdhms("", false));
			char tmp_num[8];   sprintf(tmp_num, "%8d", obs_num); string num(tmp_num);
			char tmp_wgt[18];  sprintf(tmp_wgt, "%16.8f", wgtP); string wgt(tmp_wgt);
			char tmp_rsi[18];  sprintf(tmp_rsi, "%16.8f", omcP); string rsi(tmp_rsi);
			tmp = time +
				"Range rec : " + _crt_rec +
				" satP : " + _crt_sat +
				" numP : " + num +
				" wgtP : " + wgt +
				" rsiP : " + rsi;
#endif // DEBUG_LSQPROCIF

			if (!_omc_obs_ALL(epoch, obsL, pars, omcL))
			{
				write_log_info(_log, 1, "NOTE", "omc phase obs failed");
				return false;
			};

			if (!_wgt_obs(obsL, _factorP, wgtL))
			{
				write_log_info(_log, 1, "NOTE", "weight phase obs failed");
				return false;
			};

			//if (!_prt_obs(epoch, obsL, pars, coefL)) merge
			if (!_prt_obs(epoch, pars, obsL, coefL))
			{
				write_log_info(_log, 1, "NOTE", "partial phase obs failed");
				return false;
			};

			int idx = -1;
			switch (t_gsys::band2freq(_crt_sys, band))
			{
			case FREQ_1:
				idx = pars.getParam(_crt_rec, t_gpar::AMB_L1, _crt_sat);
				break;
			case FREQ_2:
				idx = pars.getParam(_crt_rec, t_gpar::AMB_L2, _crt_sat);
				break;
			case FREQ_3:
				idx = pars.getParam(_crt_rec, t_gpar::AMB_L3, _crt_sat);
				break;
			case FREQ_4:
				idx = pars.getParam(_crt_rec, t_gpar::AMB_L4, _crt_sat);
				break;
			case FREQ_5:
				idx = pars.getParam(_crt_rec, t_gpar::AMB_L5, _crt_sat);
				break;
			}

			bool is_newamb = (pars[idx].beg == epoch) ? true : false;
			equ.add_equ(coefL, wgtL, omcL, _crt_rec, _crt_sat, t_gobscombtype(obsL, _observ), is_newamb);

#ifdef DEBUG_LOG
			string tmp2;
			string time2(_crt_epo.str_ymdhms("", false));
			char tmp_num2[8];   sprintf(tmp_num2, "%8d", obs_num);   string num2(tmp_num2);
			char tmp_wgt2[18];  sprintf(tmp_wgt2, "%16.8f", wgtL);   string wgt2(tmp_wgt2);
			char tmp_rsi2[18];  sprintf(tmp_rsi2, "%16.8f", omcL);   string rsi2(tmp_rsi2);
			tmp2 = time2 +
				"Phase rec : " + _crt_rec +
				" satL : " + _crt_sat +
				" numL : " + num2 +
				" wgtL : " + wgt2 +
				" rsiL : " + rsi2;
#endif // DEBUG_LSQPROCIF

		}
		return true;
	}

	bool t_gpceprecisemodel::_omc_obs_IF(const t_gtime& crt_epo, t_gobs& gobs, t_gallpar& pars, double& omc)
	{
		double obs = 0.0;
		t_gobs gobs1, gobs2, gobs3;
		GOBSBAND b1, b2;
		b1 = t_gsys::band_priority(_crt_sys, FREQ_1);
		b2 = t_gsys::band_priority(_crt_sys, FREQ_2);
		//sort(_SYS_freq[_crt_sys].begin(), _SYS_freq[_crt_sys].end());
		//b1 = t_gsys::band_priority(_crt_sys, _SYS_freq[_crt_sys][0]);
		//b2 = t_gsys::band_priority(_crt_sys, _SYS_freq[_crt_sys][1]);

		if (gobs.is_code())
		{
			if (_freq_index[_crt_sys][b2] != FREQ_3)
			{
				gobs1 = t_gobs(_crt_obs.select_range(b1));
				gobs2 = t_gobs(_crt_obs.select_range(b2));
				obs = _crt_obs.P3(gobs1, gobs2);
				//apply_code_DCB(_gallbias, _crt_obs, obs, _observ, &gobs1, &gobs2);
				correct_DCB("COD_R", _gallbias, _crt_obs, obs, _observ, &gobs1, &gobs2, nullptr);

			}
			else
			{
				gobs1 = t_gobs(_crt_obs.select_range(t_gsys::band_priority(_crt_sys, FREQ_1)));
				gobs2 = t_gobs(_crt_obs.select_range(t_gsys::band_priority(_crt_sys, FREQ_2)));
				gobs3 = t_gobs(_crt_obs.select_range(t_gsys::band_priority(_crt_sys, FREQ_3)));
				obs = _crt_obs.P3(gobs1, gobs3);
				//apply_cas_dcb(_gallbias, _crt_obs, obs, _observ, &gobs, nullptr);
				correct_DCB("CAS_R", _gallbias, _crt_obs, obs, _observ, &gobs1, nullptr, nullptr);

			}
		}
		else if (gobs.is_phase())
		{
			gobs1 = t_gobs(_crt_obs.select_phase(b1));
			gobs2 = t_gobs(_crt_obs.select_phase(b2));
			obs = _crt_obs.L3(gobs1, gobs2);
		}
		else
		{
			write_log_info(_log, 1, "NOTE", "error obs type");
			return false;
		}

		t_gtime epo = crt_epo;
		double ModelObs;
		ModelObs = cmpObs(epo, pars, _crt_obs, gobs, true);
		if (ModelObs < 0)
		{
			write_log_info(_log, 1, "NOTE", "Compute ModelObs error!");
			return false; //glfeng
		}

		double Amb = 0.0;
		if (gobs.is_phase())
		{
			if (!_carrier_range)
			{
				int idx = pars.getParam(_crt_rec, t_gpar::AMB_IF, _crt_sat);
				if (idx < 0)
				{
					write_log_info(_log, 1, "NOTE", "can't get the amb param!");
					return false;
				}

				// if amb is new  then init value with Obs - Modelobs
				if (double_eq(pars[idx].value(), 0.0) || pars[idx].beg == crt_epo)
				{
					pars[idx].value(obs - ModelObs);
				}
				Amb = pars[idx].value();
			}
		}

		omc = (obs - ModelObs - Amb);
#ifdef DEBUG_GLO		
		cout << fixed << setprecision(5) << setw(20) << Amb << setw(20) << omc;
		//if (gobs.is_phase()) { cout << "Sat OMC LC " << _crt_obs.sat() << "  " << omc << endl; }
		//else { cout << "Sat OMC PC " << _crt_obs.sat() << "  " << omc << endl; }
#endif
		return true;
	}

	bool t_gpceprecisemodel::_omc_obs_IF13(const t_gtime& crt_epo, t_gobs& gobs, t_gallpar& pars, double& omc)
	{
		double obs = 0.0;
		t_gobs gobs1, gobs2, gobs3;
		GOBSBAND b1 = t_gsys::band_priority(_crt_sys, FREQ_1);
		GOBSBAND b2 = t_gsys::band_priority(_crt_sys, FREQ_2);
		GOBSBAND b3 = t_gsys::band_priority(_crt_sys, FREQ_3);
		if (gobs.is_code())
		{
			gobs1 = t_gobs(_crt_obs.select_range(b1));
			gobs3 = t_gobs(_crt_obs.select_range(b3));
			obs = _crt_obs.P3(gobs1, gobs3);
		}
		else if (gobs.is_phase())
		{
			gobs1 = t_gobs(_crt_obs.select_phase(b1));
			gobs3 = t_gobs(_crt_obs.select_phase(b3));
			obs = _crt_obs.L3(gobs1, gobs3);
		}
		else
		{
			write_log_info(_log, 1, "NOTE", "error obs type");
			return false;
		}

		t_gtime epo = crt_epo;
		double ModelObs = cmpObs(epo, pars, _crt_obs, gobs, true);
		if (ModelObs < 0)
		{
			write_log_info(_log, 1, "NOTE", "Compute ModelObs error!");
			return false;
		}

		double Amb13 = 0.0;
		if (gobs.is_phase())
		{
			int idx = pars.getParam(_crt_rec, t_gpar::AMB13_IF, _crt_sat);
			if (idx < 0)
			{
				write_log_info(_log, 1, "NOTE", "can't get the amb param!");
				return false;
			}

			// if amb is new  then init value with Obs - Modelobs
			if (double_eq(pars[idx].value(), 0.0) || pars[idx].beg == crt_epo)
			{
				pars[idx].value(obs - ModelObs);
			}
			Amb13 = pars[idx].value();
		}

		omc = (obs - ModelObs - Amb13);

		return true;
	}

	bool t_gpceprecisemodel::_omc_obs_ALL(const t_gtime& crt_epo, t_gobs& gobs, t_gallpar& pars, double& omc)
	{
		double obs = 0.0;
		if (gobs.is_code())
		{
			obs = _crt_obs.obs_C(gobs);
			//apply_code_DCB(_gallbias, _crt_obs, obs, _observ, &gobs);

		}
		else if (gobs.is_phase())
		{
			obs = _crt_obs.obs_L(gobs);
		}
		else
		{
			write_log_info(_log, 1, "NOTE", "error obs type");
			return false;
		}

		t_gtime epo = crt_epo;
		double ModelObs = cmpObs(epo, pars, _crt_obs, gobs, true);
		if (ModelObs < 0)
		{
			write_log_info(_log, 1, "NOTE", "Compute ModelObs error!");
			return false;  //glfeng
		}

		double Amb = 0.0;
		if (gobs.is_phase())
		{
			int idx = -1;
			switch (t_gsys::band2freq(_crt_sys, gobs.band()))
			{
			case FREQ_1:
				idx = pars.getParam(_crt_rec, t_gpar::AMB_L1, _crt_sat);
				break;
			case FREQ_2:
				idx = pars.getParam(_crt_rec, t_gpar::AMB_L2, _crt_sat);
				break;
			case FREQ_3:
				idx = pars.getParam(_crt_rec, t_gpar::AMB_L3, _crt_sat);
				break;
			case FREQ_4:
				idx = pars.getParam(_crt_rec, t_gpar::AMB_L4, _crt_sat);
				break;
			case FREQ_5:
				idx = pars.getParam(_crt_rec, t_gpar::AMB_L5, _crt_sat);
				break;
			}

			if (idx < 0)
			{
				write_log_info(_log, 1, "NOTE", "can't get the amb param!");
				return false;
			}

			// if amb is new  then init value with Obs - Modelobs
			if (double_eq(pars[idx].value(), 0.0) || pars[idx].beg == crt_epo)
			{
				// jdhuang change for ambfix
				//pars[idx].value(obs - ModelObs);
				t_gobs gobs_P = t_gobs(_crt_obs.select_range(gobs.band(), true));
				double obs_P = _crt_obs.obs_C(gobs_P);
				//apply_DCB(_gallbias, _crt_obs, obs_P, _observ, &gobs_P);
				if (_frequency == 2)
				{
					//apply_code_DCB(_gallbias, _crt_obs, obs_P, _observ, &gobs_P);
					correct_DCB("COD_R", _gallbias, _crt_obs, obs, _observ, &gobs_P, nullptr, nullptr);
				}
				else if (_frequency == 3)
				{
					//apply_cas_dcb(_gallbias, _crt_obs, obs_P, _observ, &gobs_P, nullptr);
					correct_DCB("CAS_R", _gallbias, _crt_obs, obs, _observ, &gobs_P, nullptr, nullptr);
				}

				// jdhuang : without sion pars, a test
				pars[idx].value(obs - obs_P);
			}
			Amb = pars[idx].value();
		}

		omc = (obs - ModelObs - Amb);
		return true;
	}

	bool t_gpceprecisemodel::_wgt_obs(t_gobs& gobs, double factorP, double& wgt, bool isTriple)
	{
		GOBSTYPE type = gobs.type();
		if (type != TYPE_C &&
			type != TYPE_P &&
			type != TYPE_L)
		{
			return false;
		}

		double factor = 0.0;
		double sigRange = 0.0;
		double sigPhase = 0.0;

		if (_crt_sys == GSYS::GLO || _crt_sys == GSYS::BDS) factor = 1.5; // same with panda xiongyun
		//else if (_crt_sys == GSYS::BDS && _crt_sat <= "C05") factor = 5.0;          
		else factor = factorP;

		// get the sys sigRange
		if (type == TYPE_C ||
			type == TYPE_P)
		{
			if (_crt_obj_type == t_gdata::REC)  //add if by zhangwei
			{
				switch (_crt_sys)
				{

				case GPS: sigRange = _sigCodeGPS; break;
				case GLO: sigRange = _sigCodeGLO; break;
				case GAL: sigRange = _sigCodeGAL; break;
				case BDS: sigRange = _sigCodeBDS; break;
				case QZS: sigRange = _sigCodeQZS; break;
				default: sigRange = 0.0; return false;
				}
			}
			else if (_crt_obj_type == t_gdata::REC_LEO)
			{
				switch (_crt_sys)
				{

				case GPS: sigRange = _sigCodeGPSLEO; break;
				case GLO: sigRange = _sigCodeGLOLEO; break;
				case GAL: sigRange = _sigCodeGALLEO; break;
				case BDS: sigRange = _sigCodeBDSLEO; break;
				case QZS: sigRange = _sigCodeQZSLEO; break;
				default: sigRange = 0.0; return false;
				}
			}
			else
			{
				sigRange = 0.0;
				return false;
			}
		}
		else if (type == TYPE_L)
		{
			if (_crt_obj_type == t_gdata::REC)  //add if by zhangwei
			{
				switch (_crt_sys)
				{
				case GPS: sigPhase = _sigPhaseGPS; break;
				case GLO: sigPhase = _sigPhaseGLO; break;
				case GAL: sigPhase = _sigPhaseGAL; break;
				case BDS: sigPhase = _sigPhaseBDS; break;
				case QZS: sigPhase = _sigPhaseQZS; break;
				default: sigPhase = 0.0; return false;
				}
			}
			else if (_crt_obj_type == t_gdata::REC_LEO)
				switch (_crt_sys)
				{
				case GPS: sigPhase = _sigPhaseGPSLEO; break;
				case GLO: sigPhase = _sigPhaseGLOLEO; break;
				case GAL: sigPhase = _sigPhaseGALLEO; break;
				case BDS: sigPhase = _sigPhaseBDSLEO; break;
				case QZS: sigPhase = _sigPhaseQZSLEO; break;
				default: sigPhase = 0.0; return false;
				}
			else
			{
				sigPhase = 0.0;
				return false;
			}
		}
		else
		{
			write_log_info(_log, 1, "NOTE", "no obs type");
			return	false;
		}

		if (_crt_obj_type == t_gdata::REC)
		{
			double leo_ele = _crt_obs.ele_leo();
			double leo_ele_deg = _crt_obs.ele_leo_deg();
			double sin_leo_ele = sin(leo_ele);
			// get the _weight factor
			switch (_weight)
			{
			case DEF_OBSWEIGHT:
				cerr << "gspplsq: WeightObs (default) should not happened!\n";
				break;
			case EQUAL:
				factor *= 1;
				break;
			case SINEL:
				factor *= 1.0 / 2.0 / sin(leo_ele);
				break;
			case SINEL2:
				factor *= 1.0 / 2.0 / pow(sin_leo_ele, 2);
				break;
			case SINEL4:
				factor *= 1.0 / 2.0 / pow(sin_leo_ele, 4);
				break;
			case PARTELE:
				factor = (leo_ele_deg <= 30.0) ? factor * (1.0 / 2.0 / sin_leo_ele) : factor;
				break;
			default:
				cerr << "gspplsq: we can't deal with this WeightObs method!";
				return false;
			}
		}

		// get the combine obs factor
	// TODO finish other weight method 
		double coef1 = 0.0, coef2 = 0.0;
		GOBSBAND b1, b2;
		switch (_observ)
		{
		case gnut::IONO_FREE:
			b1 = t_gsys::band_priority(_crt_sys, FREQ_1);
			b2 = t_gsys::band_priority(_crt_sys, FREQ_2);
			//b1 = t_gsys::band_priority(_crt_sys, _SYS_freq[_crt_sys][0]);
			//b2 = t_gsys::band_priority(_crt_sys, _SYS_freq[_crt_sys][1]);
			if (isTriple)
			{
				b1 = t_gsys::band_priority(_crt_sys, FREQ_1);
				b2 = t_gsys::band_priority(_crt_sys, FREQ_3);
			}
			if (type == TYPE_C || type == TYPE_P)  //by zhangwei
			{
				_crt_obs.coef_ionofree(b1, coef1, b2, coef2);
				sigRange *= sqrt(coef1 * coef1 + coef2 * coef2);
			}
			else if (type == TYPE_L)
			{
				_crt_obs.coef_ionofree_phi(b1, coef1, b2, coef2);
				//sigPhase = sqrt(pow(sigPhase, 2)*coef1*coef1 + pow(sigPhase, 4)*coef2*coef2);
				sigPhase = sqrt(pow(sigPhase, 2) * coef1 * coef1 + pow(sigPhase, 2) * coef2 * coef2);  //by zhangwei
				sigPhase = sigPhase * _crt_obs.wavelength(b1);
			}
			else
			{
				sigRange = 0.0;
				sigPhase = 0.0;
			}
			break;
		case gnut::RAW_SINGLE:
		case gnut::RAW_DOUBLE:
		case gnut::RAW_ALL:
			b1 = gobs.band();
			//if (type == TYPE_C)
			if (type == TYPE_C || type == TYPE_P)  //by zhangwei
			{
				sigRange = sigRange;
			}
			else if (type == TYPE_L)
			{
				sigPhase = sigPhase * _crt_obs.wavelength(b1);
			}
			else
			{
				sigRange = 0.0;
				sigPhase = 0.0;
			}
			break;
		case gnut::WL_COMBIN://add by BoWong 2019/11/13
			b1 = t_gsys::band_priority(_crt_sys, FREQ_1);
			if (type == TYPE_C)
			{
				sigRange = sigRange;
			}
			else if (type == TYPE_L)
			{
				sigPhase = sigPhase * _crt_obs.wavelength(b1);
			}
			else
			{
				sigRange = 0.0;
				sigPhase = 0.0;
			}
			break;
		default:
			cerr << "gspplsq: No other weighting method have been implmented!";
			return false;
		}

		if (sigPhase == 0.0 && sigRange == 0.0) return false;

		switch (type)
		{
		case TYPE_L:
			if (sigPhase == 0.0) return false;
			wgt = 1.0 / pow(factor * sigPhase, 2);
#ifdef DEBUG_GLO
			cout << fixed << setprecision(5) << setw(20) << wgt << endl;
#endif
			break;
		case TYPE_C:
		case TYPE_P:
			if (sigRange == 0.0) return false;
			wgt = 1.0 / pow(factor * sigRange, 2);
#ifdef DEBUG_GLO
			cout << fixed << setprecision(5) << setw(20) << wgt << endl;
#endif
			break;
		default:
			return false;
		}
		//if (_crt_sys == GSYS::BDS && type == TYPE_L) wgt *= 100;
		return true;
	}

	double t_gpceprecisemodel::cmpObs(t_gtime& epoch, t_gallpar& param, t_gsatdata& obsdata, t_gobs& gobs, bool com)
	{

		//return (t_gpppmodel::cmpObs(epoch, param, obsdata, gobs, com));

		gtrace("t_gpppmodel::cmpObs");

		double spp_model;
		if (_IRC && _carrier_range)
		{
			spp_model = t_gsppmodel::cmpObsIRC(epoch, param, obsdata, gobs);
		}
		else
		{
			spp_model = t_gsppmodel::cmpObs(epoch, param, obsdata, gobs);
		}
		if (spp_model < 0) return -1;

		// Cartesian coordinates to ellipsodial coordinates
		t_gtriple ell;
		xyz2ell(_trs_rec_crd, ell, false);

		// Wind up correction 
		double wind = 0.0;
		if (gobs.is_phase())
		{
			double wavelength = 0.0;
			if (_observ == IONO_FREE)
			{
				GSYS gs = obsdata.gsys();
				if (!_isTriple)
				{
					GOBSBAND b1 = t_gsys::band_priority(gs, FREQ_1);
					GOBSBAND b2 = t_gsys::band_priority(gs, FREQ_2);
					wavelength = obsdata.wavelength_L3(b1, b2);
				}
				if (_frequency == 3 && _isTriple)
				{
					GOBSBAND b1 = t_gsys::band_priority(gs, FREQ_1);
					GOBSBAND b3 = t_gsys::band_priority(gs, FREQ_3);
					wavelength = obsdata.wavelength_L3(b1, b3);
				}
			}
			else
			{
				wavelength = obsdata.wavelength(gobs.band());
			}

			if (fabs(obsdata.wind()) > 0) wind = obsdata.wind() * wavelength;
			else wind = windUp(obsdata, _trs_rec_crd.crd_cvect()) * wavelength;
		}

		// Phase center variation correction
		double pcv_R = 0.0;
		double pcv_S = 0.0;
		double pco_R = 0.0;
		double pco_S = 0.0;
		if (_gallobj)
		{
			shared_ptr<t_gobj>  sat_obj = _gallobj->obj(_crt_sat);
			shared_ptr<t_gobj>  rec_obj = _gallobj->obj(_crt_rec);

			shared_ptr<t_gpcv>  sat_pcv;
			shared_ptr<t_gpcv>  rec_pcv;

			if (sat_obj != 0)  sat_pcv = sat_obj->pcv(epoch);
			if (rec_obj != 0)  rec_pcv = rec_obj->pcv(epoch);

			GOBS_LC lc = LC_IF;
			if (_observ != IONO_FREE)
			{
				if (t_gsys::band2freq(obsdata.gsys(), gobs.band()) == FREQ_1) lc = LC_L1;
				if (t_gsys::band2freq(obsdata.gsys(), gobs.band()) == FREQ_2) lc = LC_L2;
				if (t_gsys::band2freq(obsdata.gsys(), gobs.band()) == FREQ_3) lc = LC_L3;
				if (t_gsys::band2freq(obsdata.gsys(), gobs.band()) == FREQ_4) lc = LC_L4;
				if (t_gsys::band2freq(obsdata.gsys(), gobs.band()) == FREQ_5) lc = LC_L5;
			}


			GOBSBAND b1 = t_gsys::band_priority(obsdata.gsys(), FREQ_1);
			GOBSBAND b2 = t_gsys::band_priority(obsdata.gsys(), FREQ_2);
			if (sat_pcv && com && !_realtime)
			{
				// Satellite phase center variation
				// -- Satellite phase center offset
				t_gtriple pco(0, 0, 0);
				if (sat_pcv->pcoS(obsdata, pco, lc, b1, b2) > 0)
				{
					string antenna = sat_pcv->anten();
					t_gtriple dx(0, 0, 0);
					ColumnVector i(3);
					ColumnVector j(3);
					ColumnVector k(3);
					if (_attitudes == YAW_NOMI)	   attitude(obsdata, "", i, j, k);				//nominal modeling
					else if (_attitudes == YAW_RTCM) attitude(obsdata, obsdata.yaw(), i, j, k);	//value from RTCM used
					else							   attitude(obsdata, antenna, i, j, k);		//default
					dx[0] = pco[0] * i(1) + pco[1] * j(1) + pco[2] * k(1);
					dx[1] = pco[0] * i(2) + pco[1] * j(2) + pco[2] * k(2);
					dx[2] = pco[0] * i(3) + pco[1] * j(3) + pco[2] * k(3);
					//sat_pcv->pco_proj(pco_S, obsdata, _trs_rec_crd, dx, lc);
					sat_pcv->pco_proj(pco_S, obsdata, _trs_rec_crd, dx);//merge
					obsdata.addpco(dx);
				}

			}

			if (rec_pcv != 0)
			{
				// Receiver phase center offset
				t_gtriple pco(0.0, 0.0, 0.0);
				if (rec_pcv->pcoR(obsdata, pco, lc, b1, b2) > 0)
				{
					Matrix _rot_matrix = _RotMatrix_Ant(obsdata, _rec_epo, rec_obj, false);
					t_gtriple dx(_rot_matrix * (pco.crd_cvect()));
					//rec_pcv->pco_proj(pco_R, obsdata, _trs_rec_crd, dx, lc);
					rec_pcv->pco_proj(pco_R, obsdata, _trs_rec_crd, dx);//merge
				}
				pco_R *= -1;
			}

			if (gobs.is_phase())
			{
				if (sat_pcv != 0 && com && !_realtime)
				{
					// Satellite phase center variation
					sat_pcv->pcvS(pcv_S, obsdata, _trs_rec_crd, lc, b1, b2);
				}
				if (rec_pcv != 0) {
					// Receiver phase center variation
					//rec_pcv->is_noazi(true);
					rec_pcv->pcvR(pcv_R, obsdata, lc, b1, b2);
				}
			}
		}

#ifdef DEBUG_GLO
		cout << fixed << setprecision(5) << setw(10) << pco_R << setw(10) << pco_S << setw(10) << pcv_R << setw(10) << pcv_S
			<< setw(10) << wind;
#endif
		return spp_model +
			wind +
			pco_R +
			pco_S +
			pcv_R +
			pcv_S;
	}

}