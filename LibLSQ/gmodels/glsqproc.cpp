/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2019-04-04  created
	 -1.1 glfeng   2019-09-05  add function _batch_turboedit and variable _turboedit_software
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    glsqproc.cpp
* @brief		virtual class for process
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-18
*
*/

#include "gmodels/glsqproc.h"
#include "gset/gsetpar.h"
#include "gset/gsetamb.h"
#include "gset/gsetout.h"
#include "gset/gsetturboedit.h"
#include "gio/gfile.h"
#include "gutils/ginfolog.h"
#include "gutils/gcycleslip.h"
#include "gutils/gturboedit.h"
#include "gutils/gtime.h"
#include "gcoders/recover.h"
#include "gcoders/rinexc.h"
#include "gcoders/rb.h"
#include "gcoders/poleut1.h"
#include "gcoders/sta.h"
#include "gcoders/eop.h"
#include "gcoders/ambupd.h"
#include "gall/galleop.h"
#include "gutils/gfileconv.h"
#include "gutils/gtime.h"
#include "gmodels/glsqprocIF.h"
#include "gproc/gupdateparIRC.h"
#include "gproc/gupdateparIF.h"
#include "gutils/ginfolog.h"
#include "gproc/gupdateparALL.h"
#include "gcoders/ion.h"
#include "gutils/gstring.h"
#include "gutils//gleo.h"
#include <algorithm>
#include <thread>


//#define DEBUG_LOG
//#define DEBUG_BPL_LOG
// Judge whether file exsit
#if  defined _WIN32 || defined _WIN64
#include <io.h>
#include <direct.h> 
#define ACCESS(fileName,accessMode) _access(fileName,accessMode)
#define MKDIR(path) _mkdir(path)
#else
#include <unistd.h>
#include <sys/stat.h>
#define ACCESS(fileName,accessMode) access(fileName,accessMode)
#define MKDIR(path) mkdir(path,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif


namespace great
{
	t_glsqproc::t_glsqproc()
	{}

	t_glsqproc::t_glsqproc(t_gsetbase* set, t_gallproc* data, t_glog* log) :
		_glog(log),
		_gset(set),
		_gall_proc(data),
		_recclk_threshold(1e-8)
	{
		// check log
		if (_glog == nullptr)
		{
			_myLog = _create_log(_glog);
		}

		
		// get settings from xml
		_init_xml_settings(set);
		_init_lsq_settings(set, data, log);

		// used for updating parameters in kin mode (for PPP or LEO kinematic mode)
		//yjqin
		if (_leo_list.size() == 0 && _crd_est == KIN)
		{
			_read_kinfile();
		}

		
	}

	t_glsqproc::~t_glsqproc()
	{
		// Init Paramter in Constructor???
		if (_lsq)
		{
			delete _lsq;
			_lsq = nullptr;
		}

		if (_myLog)
		{
			_delete_log(_glog);
		}
	}

	//added by liyuhao
	void t_glsqproc::add_coder(const vector<t_gcoder*>& coder)
	{
		_gcoder = coder;
	}

	bool t_glsqproc::ProcessBatch(t_gallproc* data, const t_gtime& beg, const t_gtime& end)
	{
		return true;
	}

	bool t_glsqproc::GenerateProduct()
	{
		return false;
	}

	bool t_glsqproc::_init_crs_erp_pars(t_glsq* lsq)
	{
		int ipar_num = _lsq->_x_solve.parNumber();

		// init EOP parmeter
		// init the Eop dx dy vdx vdy
		t_gpar par_xpole("", t_gpar::XPOLE, ++ipar_num, "");
		t_gpar par_ypole("", t_gpar::YPOLE, ++ipar_num, "");
		t_gpar par_dxpole("", t_gpar::DXPOLE, ++ipar_num, "");
		t_gpar par_dypole("", t_gpar::DYPOLE, ++ipar_num, "");
		t_gpar par_ut1("", t_gpar::UT1, ++ipar_num, "");
		t_gpar par_dut1("", t_gpar::DUT1, ++ipar_num, "");

		// TODO:init value time sigma		
		map<string, double> data1, data2;

		try
		{
			if (!_gdata_erp) _gdata_erp = dynamic_cast<t_gpoleut1*>((*_gall_proc)[t_gdata::ALLPOLEUT1]);
			if (!_gdata_erp) return false;
			_gdata_erp->getEopData(_beg_time, data1);
			_gdata_erp->getEopData(_end_time, data2);
		}
		catch (exception e)
		{
			write_log_info(_glog, 0, "ERROR", e.what());
			return false;
		}

		double sig_xpole = dynamic_cast<t_gsetpar*>(_gset)->sigXpole();
		double sig_ypole = dynamic_cast<t_gsetpar*>(_gset)->sigYpole();
		double sig_dxpole = dynamic_cast<t_gsetpar*>(_gset)->sigDxpole();
		double sig_dypole = dynamic_cast<t_gsetpar*>(_gset)->sigDypole();
		double sig_ut1 = dynamic_cast<t_gsetpar*>(_gset)->sigUt1();
		double sig_dut1 = dynamic_cast<t_gsetpar*>(_gset)->sigDut1();

		par_xpole.setTime(_beg_time, _end_time);
		par_xpole.value(data1["XPOLE"]);
		par_xpole.apriori(sig_xpole);
		par_ypole.setTime(_beg_time, _end_time);
		par_ypole.value(data1["YPOLE"]);
		par_ypole.apriori(sig_ypole);
		par_dxpole.setTime(_beg_time, _end_time);
		par_dxpole.value(data2["XPOLE"] - data1["XPOLE"]);
		par_dxpole.apriori(sig_dxpole);
		par_dypole.setTime(_beg_time, _end_time);
		par_dypole.value(data2["YPOLE"] - data1["YPOLE"]);
		par_dypole.apriori(sig_dypole);
		par_ut1.setTime(_beg_time, _end_time);
		par_ut1.value(data1["UT1-TAI"]);
		par_ut1.apriori(sig_ut1);
		par_dut1.setTime(_beg_time, _end_time);
		par_dut1.value(data2["UT1-TAI"] - data1["UT1-TAI"]);
		par_dut1.apriori(sig_dut1);

		// add the Estimator
		lsq->add_parameter(par_xpole);
		lsq->add_parameter(par_ypole);
		lsq->add_parameter(par_dxpole);
		lsq->add_parameter(par_dypole);
		lsq->add_parameter(par_ut1);
		lsq->add_parameter(par_dut1);

		return true;
	}

	bool t_glsqproc::_init_rec_trop_pars(t_glsq* lsq)
	{
		if (_rec_list.empty() || !_lsq || !_gset) return false;

		int ipar_num = _lsq->_x_solve.parNumber();
		auto it_beg = _rec_list.begin();
		auto it_end = _rec_list.end();
		for (auto it = it_beg; it != it_end; it++)
		{
			// init the ztd
			string rec = *it;
			t_gpar par_ztd(rec, t_gpar::TRP, ++ipar_num, "");
			shared_ptr<t_gobj> _grec = _gall_obj->obj(rec);
			if (_grec->id_type() == t_gdata::REC_LEO)
			{
				//it++;  //by zhangwei
				continue;
			}
			double sig_ztd = dynamic_cast<t_gsetpar*>(_gset)->sigZtd(rec);
			double sig_TropPd = dynamic_cast<t_gsetpar*>(_gset)->sigTropPd(rec);
			par_ztd.apriori(sig_ztd);
			par_ztd.value(0.0);
			par_ztd.setTime(_beg_time, _beg_time + _pwc_intv * 60.0);
			lsq->add_parameter(par_ztd);
			if (_ztd_model == PWC)
			{
				lsq->add_par_state_equ(par_ztd.parType, 1, _pwc_intv / 60.0, sig_TropPd);
			}
			if (_ztd_model == STO)
			{
				lsq->add_par_state_equ(par_ztd.parType, 1, _obs_intv / 3600.0, sig_TropPd);
			}

			// init tropospheric gradients
			if (_tropo_grad)
			{
				double sig_grd = dynamic_cast<t_gsetpar*>(_gset)->sigGRD(rec);
				double sig_GrdPd = dynamic_cast<t_gsetpar*>(_gset)->sigGrdPd(rec);
				t_gpar par_grdN(rec, t_gpar::GRD_N, ++ipar_num, "");
				par_grdN.setTime(_beg_time, _beg_time + _grd_pwc_intv * 60.0);
				par_grdN.apriori(sig_grd);
				par_grdN.value(0.0);
				lsq->add_parameter(par_grdN);
				
				t_gpar par_grdE(rec, t_gpar::GRD_E, ++ipar_num, "");
				par_grdE.setTime(_beg_time, _beg_time + _grd_pwc_intv * 60.0);
				par_grdE.value(0.0);
				par_grdE.apriori(sig_grd);
				lsq->add_parameter(par_grdE);

				if (_grd_model == GRD_PWC)
				{
					lsq->add_par_state_equ(par_grdN.parType, 1, _grd_pwc_intv / 60.0, sig_GrdPd);
					lsq->add_par_state_equ(par_grdE.parType, 1, _grd_pwc_intv / 60.0, sig_GrdPd);
				}
				if (_grd_model == GRD_STO)
				{
					lsq->add_par_state_equ(par_grdN.parType, 1, _obs_intv / 3600.0, sig_GrdPd);
					lsq->add_par_state_equ(par_grdE.parType, 1, _obs_intv / 3600.0, sig_GrdPd);
				}
			}
		}
		return true;
	}

	bool t_glsqproc::_init_gns_vion_pars(t_glsq* lsq)
	{
		if (_rec_list.empty() || _sat_list.empty() || !_lsq || !_gset) return false;

		int ipar_num = _lsq->_x_solve.parNumber();

		auto it_beg = _rec_list.begin();
		auto it_end = _rec_list.end();
		for (auto it = it_beg; it != it_end; it++)
		{
			string rec = *it;
			// init the ztd
			t_gpar par_ion(rec, t_gpar::VION, ++ipar_num, "");
			double sig_ion = dynamic_cast<t_gsetpar*>(_gset)->sigVion(rec);
			par_ion.apriori(sig_ion);

			// init in precisemodel
			double dt = dynamic_cast<t_gsetgen*>(_gset)->sampling(); // minute
			if (dt == 0.0) return false;
			par_ion.value(0.0);
			par_ion.setTime(_beg_time, _beg_time + dt);
			lsq->add_parameter(par_ion);
			double sigIonoPd = dynamic_cast<t_gsetpar*>(_gset)->sigIonoPd(rec);
			lsq->add_par_state_equ(par_ion.parType, 1, dt / 60.0, sigIonoPd);
		}
		return true;
	}

	bool t_glsqproc::_init_gns_sion_pars(t_glsq* lsq)
	{
		if (_rec_list.empty() || _sat_list.empty() || !_lsq || !_gset) return false;

		for (const auto& it_rec : _rec_list)
		{
			for (const auto& it_sat : _sat_list)
			{
				string rec = it_rec;
				string sat = it_sat;
				//double sig_ion = dynamic_cast<t_gsetpar*>(_gset)->sigVion(rec);  glfeng
				double sig_ion = dynamic_cast<t_gsetpar*>(_gset)->sigSion(rec);

				int ipar_num = _lsq->_x_solve.parNumber();
				t_gpar par_ion(rec, t_gpar::SION, ++ipar_num, sat);
				par_ion.apriori(sig_ion);
				par_ion.value(0.0);
				par_ion.setTime(_beg_time, _beg_time);
				lsq->add_parameter(par_ion);
			}
		}
		return true;
	}

	bool t_glsqproc::_init_rec_clk_pars(t_glsq* lsq)
	{

		if (_rec_list.empty() || !lsq || !_gset)
		{
			return false;
		}

		string site_ref = dynamic_cast<t_gsetproc*>(_gset)->ref_clk();
		auto it_beg = _rec_list.begin();
		auto it_end = _rec_list.end();
		for (auto it = it_beg; it != it_end; it++)
		{

			int ipar_num = lsq->_x_solve.parNumber();

			string site = *it;
			double sigclk = dynamic_cast<t_gsetpar*>(_gset)->sigRecCLK(site);
			double sig_refclk = dynamic_cast<t_gsetproc*>(_gset)->sig_ref_clk();

			if (site == site_ref)
			{
				if (sig_refclk == 0.0)
				{
					//it++;
					continue;
				}
				else
				{
					sigclk = sig_refclk;
				}
			}

			//intit the Clock
			t_gpar par_clk(site, t_gpar::CLK, ++ipar_num, "");
			par_clk.value(0.0);
			par_clk.setTime(_beg_time, _beg_time);
			par_clk.apriori(sigclk);
			lsq->add_parameter(par_clk);
		}
		return true;
	}

	bool t_glsqproc::_init_rec_clk13_pars(t_glsq* lsq)
	{

		if (_rec_list.empty() || !lsq || !_gset)
		{
			return false;
		}

		string site_ref = dynamic_cast<t_gsetproc*>(_gset)->ref_clk();
		auto it_beg = _rec_list.begin();
		auto it_end = _rec_list.end();
		for (auto it = it_beg; it != it_end; it++)
		{

			int ipar_num = lsq->_x_solve.parNumber();

			string site = *it;
			double sigclk = dynamic_cast<t_gsetpar*>(_gset)->sigRecCLK(site);
			double sig_refclk = dynamic_cast<t_gsetproc*>(_gset)->sig_ref_clk();

			if (site == site_ref)
			{
				if (sig_refclk == 0.0)
				{
					//it++;
					continue;
				}
				else
				{
					sigclk = sig_refclk;
				}
			}

			//intit the Clock
			set<string>::const_iterator iter_sys = _sys_list.begin();
			t_gpar par_clk;
			for (iter_sys; iter_sys != _sys_list.end(); iter_sys++)
			{
				GSYS   crt_sys = t_gsys::str2gsys(*iter_sys);
				switch (crt_sys)
				{
				case gnut::GPS:
					par_clk=t_gpar(site, t_gpar::G15CLK, ++ipar_num, "");
					par_clk.value(0.0);
					par_clk.setTime(_beg_time, _beg_time);
					par_clk.apriori(sigclk);
					lsq->add_parameter(par_clk);
					break;
				case gnut::GAL:
					par_clk = t_gpar(site, t_gpar::E15CLK, ++ipar_num, "");
					par_clk.value(0.0);
					par_clk.setTime(_beg_time, _beg_time);
					par_clk.apriori(sigclk);
					lsq->add_parameter(par_clk);
					break;
				case gnut::BDS:
					par_clk = t_gpar(site, t_gpar::C15CLK, ++ipar_num, "");
					par_clk.value(0.0);
					par_clk.setTime(_beg_time, _beg_time);
					par_clk.apriori(sigclk);
					lsq->add_parameter(par_clk);
					break;
				case gnut::QZS:
					par_clk = t_gpar(site, t_gpar::J15CLK, ++ipar_num, "");
					par_clk.value(0.0);
					par_clk.setTime(_beg_time, _beg_time);
					par_clk.apriori(sigclk);
					lsq->add_parameter(par_clk);
					break;
				default:
					break;
				}
			}
		}
		return true;
	}

	bool t_glsqproc::_init_rec_clk_pars_irc(t_glsq* lsq)
	{

		if (_rec_list.empty() || !lsq || !_gset)
		{
			return false;
		}

		string site_ref = dynamic_cast<t_gsetproc*>(_gset)->ref_clk();
		auto it_beg = _rec_list.begin();
		auto it_end = _rec_list.end();
		set<string>::const_iterator itsys_beg = _sys_list.begin();
		set<string>::const_iterator itsys_end = _sys_list.end();
		for (auto it = it_beg; it != it_end; it++)
		{
			for (set<string>::const_iterator itsys = itsys_beg; itsys != itsys_end; itsys++)
			{
				int ipar_num = lsq->_x_solve.parNumber();

				string site = *it;
				double sigclk = dynamic_cast<t_gsetpar*>(_gset)->sigRecCLK(site);
				double sig_refclk = dynamic_cast<t_gsetproc*>(_gset)->sig_ref_clk();

				if (site == site_ref)
				{
					if (sig_refclk == 0.0)
					{
						continue;
					}
					else
					{
						sigclk = sig_refclk;
					}
				}

				//intit the Clock
				t_gpar::t_type clktype = t_gpar::t_type::NO_DEF;
				GSYS crt_sys = t_gsys::str2gsys(*itsys);
				switch (crt_sys)
				{
					case GSYS::GPS:
						clktype = t_gpar::CLK_G;
						break;
					case GSYS::BDS:
						clktype = t_gpar::CLK_C;
						break;
					case GSYS::GAL:
						clktype = t_gpar::CLK_E;
						break;
					case GSYS::GLO:
						clktype = t_gpar::CLK_R;
						break;
					case GSYS::QZS:
						clktype = t_gpar::CLK_J;
						break;
					default:
						break;
				}
				t_gpar par_clk(site, clktype, ++ipar_num, "");
				par_clk.value(0.0);
				par_clk.setTime(_beg_time, _beg_time);
				par_clk.apriori(sigclk);
				lsq->add_parameter(par_clk);
			}
		}
		return true;
	}

	bool t_glsqproc::_init_sat_clk_pars(t_glsq* lsq)
	{
		if (!lsq || !_gset || _sat_list.empty()) return false;

		int ipar_num = _lsq->_x_solve.parNumber();
		auto it_beg = _sat_list.begin();
		auto it_end = _sat_list.end();
		string sat_ref = dynamic_cast<t_gsetproc*>(_gset)->ref_clk();
		double sig_refclk = dynamic_cast<t_gsetproc*>(_gset)->sig_ref_clk();
		
		for (auto iter = it_beg; iter != it_end; iter++)
		{
			string satID = *iter;
			double sigclk = dynamic_cast<t_gsetpar*>(_gset)->sigSatCLK(satID);
			if (satID == sat_ref && sig_refclk == 0.0)
			{
				continue;
			}
			else if (satID == sat_ref && sig_refclk != 0.0)
			{
				sigclk = sig_refclk;
			}
			else
			{
				sigclk = dynamic_cast<t_gsetpar*>(_gset)->sigSatCLK(satID);
			}

			t_gpar par_satclk("", t_gpar::CLK_SAT, ++ipar_num, satID);
			par_satclk.setTime(_beg_time, _beg_time);
			par_satclk.apriori(sigclk);
			par_satclk.value(0.0);
			lsq->add_parameter(par_satclk);
		}
		return true;
	}

	bool t_glsqproc::_init_sat_ifcb_pars(t_glsq* lsq)
	{
		if (!lsq || !_gset || _sat_list.empty()) return false;

		int ipar_num = _lsq->_x_solve.parNumber();
		auto it_beg = _sat_list.begin();
		auto it_end = _sat_list.end();
		//string sat_ref = dynamic_cast<t_gsetproc*>(_gset)->ref_clk();
		//double sig_refclk = dynamic_cast<t_gsetproc*>(_gset)->sig_ref_clk();

		for (auto iter = it_beg; iter != it_end; iter++)
		{
			string satID = *iter;
			double sigclk = dynamic_cast<t_gsetpar*>(_gset)->sigSatCLK(satID);
			//if (satID == sat_ref && sig_refclk == 0.0)
			//{
			//	continue;
			//}
			//else if (satID == sat_ref && sig_refclk != 0.0)
			//{
			//	sigclk = sig_refclk;
			//}
			//else
			//{
			//	sigclk = dynamic_cast<t_gsetpar*>(_gset)->sigSatCLK(satID);
			//}

			//t_gpar::t_type clktype;
			//GSYS crt_sys = t_gsys::str2gsys(satID.substr(0, 1));
			//switch (crt_sys)
			//{
			//	case GSYS::GPS:
			//		clktype = t_gpar::G15CLK_SAT;
			//		break;
			//	case GSYS::BDS:
			//		clktype = t_gpar::C15CLK_SAT;
			//		break;
			//	case GSYS::GAL:
			//		clktype = t_gpar::E15CLK_SAT;
			//		break;
			//	case GSYS::QZS:
			//		clktype = t_gpar::J15CLK_SAT;
			//		break;
			//	default:
			//		break;
			//}

			//t_gpar par_satclk("", clktype, ++ipar_num, satID);
			t_gpar par_satclk("", t_gpar::IFCB_F3, ++ipar_num, satID);
			par_satclk.setTime(_beg_time, _beg_time);
			par_satclk.apriori(3000);
			par_satclk.value(0.0);
			lsq->add_parameter(par_satclk);
		}
		return true;
	}

	bool t_glsqproc::_init_gns_ics_pars(t_glsq* lsq)
	{
		if (!_gall_gns) return false;

		int ipar_num = _lsq->_x_solve.parNumber();
		// cycle sat list

		set<string>::const_iterator sat_end = _sat_list.end();
		set<string>::const_iterator sat_beg = _sat_list.begin();
		for (auto sat_iter = sat_beg; sat_iter != sat_end; sat_iter++)
		{
			auto ics = _gall_gns->get_nav_ics(*sat_iter);
			auto par_list = ics->get_par_list();
			if (par_list.empty())
			{
				return false;
			}

			// cycle par list
			auto it_beg = par_list.begin();
			auto it_end = par_list.end();
			for (auto par_iter = it_beg; par_iter != it_end; par_iter++)
			{
				t_gpar::t_type par_type = t_gics_nav::_icspar2gpar(par_iter->first);

				// init 
				t_gpar par_temp("", par_type, ++ipar_num, *sat_iter);
				par_temp.setTime(par_iter->second.beg_time, par_iter->second.end_time);
				par_temp.value(par_iter->second.value);

				auto sig_pos = dynamic_cast<t_gsetpar*>(_gset)->sigSatPos(*sat_iter);
				auto sig_vel = dynamic_cast<t_gsetpar*>(_gset)->sigSatVel(*sat_iter);
				auto sig_force = dynamic_cast<t_gsetpar*>(_gset)->sigSatEcom(*sat_iter);
				string id = "sig" + par_temp.str_type().substr(0, 2);
				if (par_type >= t_gpar::PXSAT && par_type <= t_gpar::PZSAT)
				{
					par_temp.apriori(sig_pos[id] / 1000.0); // change unit to km
				}
				else if (par_type >= t_gpar::VXSAT && par_type <= t_gpar::VZSAT) 
				{
					par_temp.apriori(sig_vel[id] / 1000.0); // change unit to km
				}
				else if (par_type >= t_gpar::ECOM_D0 && par_type <= t_gpar::ECOM2_X1s)
				{
					double tmp = sig_force[id];
					if (double_eq(tmp, 0.0)) tmp = 0.1;
					par_temp.apriori(tmp);
				}
				lsq->add_parameter(par_temp);
			}
		}

		return true;
	}

	bool t_glsqproc::_init_leo_ics_pars(t_glsq* lsq)
	{

		int ipar_num = _lsq->_x_solve.parNumber();
		// cycle sat list
		for (auto sat_iter = _leo_list.begin(); sat_iter != _leo_list.end(); sat_iter++)
		{
			string leosat = strs2strl(*sat_iter);
			auto ics = _gall_leo->get_leo_ics(leosat);
			map<ICS_PAR, vector<t_icsleo_record>>* leopar;
			leopar = ics->leoics();
			vector<t_gpar> pars;
			int npar = 0;
			// cycle par list
			for (auto par_iter = leopar->begin(); par_iter != leopar->end(); par_iter++)
			{
				auto sig_pos = dynamic_cast<t_gsetpar*>(_gset)->sigLeoPos(*sat_iter);
				auto sig_vel = dynamic_cast<t_gsetpar*>(_gset)->sigLeoVel(*sat_iter);
				double sig_dyn = dynamic_cast<t_gsetpar*>(_gset)->sigLeoDyn(*sat_iter);

				t_gpar::t_type par_type = t_gics_leo::icspar2gpar(par_iter->first);
				vector<t_icsleo_record> onepar = par_iter->second;
				int npwc = round((_end_time - _beg_time) / (onepar[0].end - onepar[0].beg));
				int nparics = onepar.size();
				if (npwc > nparics)
				{
					write_log_info(_glog, 0, "ERROR", "The end time: " + _end_time.str_hms() + " is beyond the ics time table: from " + onepar[0].beg.str_hms() +  " to " + onepar[(int64_t)npar - 1].end.str_hms());
					return false;
				}
				for (int k = 0; k < npwc; k++) //nparics to npwc
				{
					//int indexPWC = 0;
					if (onepar[k].beg >= _beg_time && onepar[k].beg < _end_time)
					{
						npar++;
						t_gpar par_temp(*sat_iter, par_type, ++ipar_num, "", false);
						if (k == npwc - 1)
						{
							par_temp.setTime(onepar[k].beg, _end_time);
						}
						else
						{
							par_temp.setTime(onepar[k].beg, onepar[k].end);
						}
						par_temp.value(onepar[k].value);
						//par_temp.apriori(onepar[k].sig);
						par_temp.value(onepar[k].value);
						
						string id = "sig" + par_temp.str_type().substr(0, 2);
						if (par_type >= t_gpar::PXSAT && par_type <= t_gpar::PZSAT) {
							par_temp.apriori(sig_pos[id] / 1000.0); // change unit to km
						}
						else if (par_type >= t_gpar::VXSAT && par_type <= t_gpar::VZSAT) {
							par_temp.apriori(sig_vel[id] / 1000.0); // change unit to km
						}
						else {
							par_temp.apriori(sig_dyn);
						}

						if (onepar.size() > 1)  // PWC parameters
						{
							par_temp.nPWC = npwc;
							par_temp.lPWC = true;
							par_temp.piecewise_time = onepar[0].piecewise_time;
							//par_temp.indPWC = indexPWC;
							//indexPWC++;
						}
						_lsq->add_parameter(par_temp);
						pars.push_back(par_temp);
						
					}
				}
				
			}
			_lsq->_x_solve.setOrbPar(leosat, pars);
		}
		return true;
	}


	// set sys bia pars
	bool t_glsqproc::_init_sys_bia_pars(t_glsq* lsq)
	{
		// 单系统，无电离层组合
		if (_sys_list.size() < 2 && _obs_mode == IONO_FREE && _frequency == 2) return true;
		if (_sys_list.empty())  return false;

		// get the ref sys
		//string ref_sys = _sys_list.count("GPS") > 0 ? "GPS" : *_sys_list.begin();
		string ref_sys;
		if      (_sys_list.count("GPS") > 0)  ref_sys = "GPS";
		else if (_sys_list.count("LEO") > 0)  ref_sys = "LEO";
		else                                  ref_sys = *_sys_list.begin();
		for (const auto& sys : _sys_list)
		{
			if (sys == ref_sys) continue;
			for (const auto& crt_rec : _rec_list)
			{
				t_gpar par_ISB;
				par_ISB.value(0.0);
				par_ISB.apriori(3000.0);

				GSYS   crt_sys = t_gsys::str2gsys(sys);
				switch (crt_sys)
				{
				case gnut::GPS:
					break;
				case gnut::LEO:  // add for PPPLEO settings, regard LEO/GPS as the same system, glfeng
					break;
				case gnut::GAL:
					par_ISB = t_gpar(crt_rec, t_gpar::GAL_ISB, lsq->_x_solve.parNumber() + 1, "");
					break;
				case gnut::BDS:
					par_ISB = t_gpar(crt_rec, t_gpar::BDS_ISB, lsq->_x_solve.parNumber() + 1, "");
					break;
				case gnut::QZS:
					par_ISB = t_gpar(crt_rec, t_gpar::QZS_ISB, lsq->_x_solve.parNumber() + 1, "");
					break;
				case gnut::GLO:
					//par_ISB = t_gpar(crt_rec, t_gpar::GLO_ISB, lsq->_x_solve.parNumber() + 1, "");
					break;
				default:
					break;
				}
				
				// ISB par type
				if (crt_sys != GLO && crt_sys != LEO)
				{
					_set_isb_par_setting(lsq, par_ISB, _sysbias_model, 0.0, 3000.0);
				}	

				// IFB, special for GLONASS
				if (crt_sys != gnut::GLO) continue;
				if (_sysbias_model != AUTO_CON && _sysbias_model != AUTO_RWK && _sysbias_model != AUTO_WHIT) continue;
				for (const auto& iter_sat : _sat_list)
				{
					if ((iter_sat).substr(0, 1) != "R") continue;
					t_gpar par_IFB;
					par_IFB = t_gpar(crt_rec, t_gpar::GLO_IFB, lsq->_x_solve.parNumber() + 1, iter_sat);
					_set_isb_par_setting(lsq, par_IFB, _sysbias_model, 0.0, 3000.0);
				} // next glonass sat

			} // next rec
		}// next system
		
		return true;
	}

	bool t_glsqproc::_init_sys_rec_ifb_pars(t_glsq * lsq)
	{
		if (_ifb_model != IFB_MODEL::EST_REC_IFB) return false;

		// if raw all mode, init IFB par for each satellite
		// one sat one ifb
		map<GSYS,set<t_gpar::t_type>> par_list;

		// 最多支持处理五频数据
		switch (_frequency)
		{
		case 5:
			par_list[GPS].insert(t_gpar::t_type::GPS_REC_IFB_C5);
			par_list[GAL].insert(t_gpar::t_type::GAL_REC_IFB_C5);
			par_list[BDS].insert(t_gpar::t_type::BDS_REC_IFB_C5);
		case 4:
			par_list[GPS].insert(t_gpar::t_type::GPS_REC_IFB_C4);
			par_list[GAL].insert(t_gpar::t_type::GAL_REC_IFB_C4);
			par_list[BDS].insert(t_gpar::t_type::BDS_REC_IFB_C4);
		case 3:
			par_list[GPS].insert(t_gpar::t_type::GPS_REC_IFB_C3);
			par_list[GAL].insert(t_gpar::t_type::GAL_REC_IFB_C3);
			par_list[BDS].insert(t_gpar::t_type::BDS_REC_IFB_C3);
		case 2:
			break;
		default:
			return false;
		}

		// one rec one ifb
		for (const auto& rec : _rec_list)
		{
			t_gpar par_IFB;
			par_IFB.value(0.0);
			par_IFB.apriori(3000.0);

			for (const auto& sys : _sys_list)
			{
				auto gsys = t_gsys::str2gsys(sys);
				switch (gsys)
				{
				case GPS:
				{
					for (const auto& par_type : par_list[GPS])
					{
						par_IFB = t_gpar(rec, par_type, lsq->_x_solve.parNumber() + 1, "");
						_set_isb_par_setting(lsq, par_IFB, _sysbias_model, 0.0, 9000.0);
						//par_IFB.setTime(_beg_time, _beg_time);
						//par_IFB.apriori(9000.0);
						//lsq->add_parameter(par_IFB);
					}
					break;
				}
				case GAL:
				{
					for (const auto& par_type : par_list[GAL])
					{
						par_IFB = t_gpar(rec, par_type, lsq->_x_solve.parNumber() + 1, "");

						_set_isb_par_setting(lsq, par_IFB, _sysbias_model, 0.0, 9000.0);

						//par_IFB.setTime(_beg_time, _beg_time);
						//par_IFB.apriori(9000.0);
						//lsq->add_parameter(par_IFB);

						//// IFB par type
						//if (_sysbias_model == AUTO_CON || _sysbias_model == ISB_CON)
						//{
						//	par_IFB.setTime(_beg_time, _end_time);
						//	lsq->add_parameter(par_IFB);
						//}
						//else if (_sysbias_model == AUTO_WHIT)
						//{
						//	par_IFB.setTime(_beg_time, _beg_time);
						//	_lsq->add_parameter(par_IFB);
						//}
						//else
						//{
						//	par_IFB.setTime(_beg_time, _beg_time);
						//	lsq->add_parameter(par_IFB);
						//	lsq->add_par_state_equ(par_IFB.parType, 1, 24, 3000.0);
						//}
					}
					break;
				}
				case BDS:
					for (const auto& par_type : par_list[BDS])
					{
						par_IFB = t_gpar(rec, par_type, lsq->_x_solve.parNumber() + 1, "");
						_set_isb_par_setting(lsq, par_IFB, _sysbias_model, 0.0, 9000.0);
						//par_IFB.setTime(_beg_time, _beg_time);
						//par_IFB.apriori(9000.0);
						//lsq->add_parameter(par_IFB);

						//// IFB par type
						//if (_sysbias_model == AUTO_CON || _sysbias_model == ISB_CON)
						//{
						//	par_IFB.setTime(_beg_time, _end_time);
						//	lsq->add_parameter(par_IFB);
						//}
						//else if (_sysbias_model == AUTO_WHIT)
						//{
						//	par_IFB.setTime(_beg_time, _beg_time);
						//	_lsq->add_parameter(par_IFB);
						//}
						//else
						//{
						//	par_IFB.setTime(_beg_time, _beg_time);
						//	lsq->add_parameter(par_IFB);
						//	lsq->add_par_state_equ(par_IFB.parType, 1, 24, 3000.0);
						//}
					}
					break;
				default:
					throw logic_error("can not support sys : " + sys);
				}
			}
		}

		return true;
	}

	bool t_glsqproc::_init_sys_sat_ifb_pars(t_glsq * lsq)
	{
		if (_ifb_model != IFB_MODEL::EST_SAT_IFB) return false;


		// if raw all mode, init IFB par for each satellite
		// one sat one ifb
		set<t_gpar::t_type> par_list;

		// 最多支持处理五频数据
		switch (_frequency)
		{
		case 5:
			par_list.insert(t_gpar::t_type::SAT_IFB_C5);
		case 4:
			par_list.insert(t_gpar::t_type::SAT_IFB_C4);
		case 3:
			par_list.insert(t_gpar::t_type::SAT_IFB_C3);
		case 2:
			break;
		default:
			return false;
		}

		for (const auto& par_type : par_list)
		{
			for (const auto& sat : _sat_list)
			{
				if (_ifb_model != IFB_MODEL::EST_SAT_IFB) continue;
				// pass GLONASS
				if (sat.find("R") != string::npos) continue;
				t_gpar par_IFB;
				//par_IFB.value(0.0);
				//par_IFB.apriori(3000.0);

				par_IFB = t_gpar("", par_type, lsq->_x_solve.parNumber() + 1, sat);
				//lsq->add_parameter(par_IFB);

				_set_isb_par_setting(lsq, par_IFB, _sysbias_model, 0.0, 3000.0);
			}
		}

		return true;
	}

	bool t_glsqproc::_init_sys_ifb_pars(t_glsq * lsq)
	{
		if (_ifb_model != IFB_MODEL::EST_IFB) return false;

		// if raw all mode, init IFB par for each satellite
		// one sat one ifb
		set<t_gpar::t_type> par_list;

		// 最多支持处理五频数据
		switch (_frequency)
		{
		case 5:
			par_list.insert(t_gpar::t_type::IFB_C5);
		case 4:
			par_list.insert(t_gpar::t_type::IFB_C4);
		case 3:
			par_list.insert(t_gpar::t_type::IFB_C3);
		case 2:
			break;
		default:
			return false;
		}

		// one sat-rec one ifb
		for (const auto& sat : _sat_list)
		{
			// pass GLONASS
			if (sat.find("R") != string::npos) continue;

			for (const auto& rec : _rec_list)
			{
				t_gpar par_IFB;
				par_IFB.value(0.0);
				par_IFB.apriori(3000.0);

				for (const auto& sys : _sys_list)
				{
					for (const auto& par_type : par_list)
					{
						par_IFB = t_gpar(rec, par_type, lsq->_x_solve.parNumber() + 1, sat);

						_set_isb_par_setting(lsq, par_IFB, _sysbias_model, 0.0, 3000.0);

						//// ISB par type
						//if (_sysbias_model == AUTO_CON || _sysbias_model == ISB_CON)
						//{
						//	par_IFB.setTime(_beg_time, _end_time);
						//	lsq->add_parameter(par_IFB);
						//}
						//else
						//{
						//	par_IFB.setTime(_beg_time, _beg_time);
						//	lsq->add_parameter(par_IFB);
						//	lsq->add_par_state_equ(par_IFB.parType, 1, 24, 3000.0);
						//}
					}
				}
			}
		}
		return true;
	}
	
	// set sys bia pars
	bool t_glsqproc::_init_sys_bia_pars_pce(t_glsq* lsq)
	{
		if (_sys_list.size() < 2) return true;
		if (_sys_list.empty()) return false;

		// get the ref sys
		string ref_sys;
		set<string>::const_iterator iter_sys = _sys_list.begin();
		if (_sys_list.count("GPS") > 0) ref_sys = "GPS";
		else ref_sys = *iter_sys;


		for (iter_sys; iter_sys != _sys_list.end(); iter_sys++)
		{
			if (*iter_sys == ref_sys) continue;

			set<string>::const_iterator rec_beg = _rec_list.begin();
			set<string>::const_iterator rec_end = _rec_list.end();
			set<string>::const_iterator rec_iter = rec_beg;

			for (rec_iter; rec_iter != rec_end; ++rec_iter)
			{
				t_gpar par_ISB, par_IFB;
				int ipar_num = _lsq->_x_solve.parNumber();
				GSYS   crt_sys = t_gsys::str2gsys(*iter_sys);
				string crt_rec = *rec_iter;

				switch (crt_sys)
				{
				case gnut::GPS:
					break;
				case gnut::GAL:
					par_ISB = t_gpar(crt_rec, t_gpar::GAL_ISB, ++ipar_num, "");
					par_ISB.value(0.0);
					par_ISB.setTime(_beg_time, _end_time);
					par_ISB.apriori(3000.0);
					_lsq->add_parameter(par_ISB);
					break;
				case gnut::BDS:
					par_ISB = t_gpar(crt_rec, t_gpar::BDS_ISB, ++ipar_num, "");
					par_ISB.value(0.0);
					par_ISB.setTime(_beg_time, _end_time);
					par_ISB.apriori(3000.0);
					_lsq->add_parameter(par_ISB);
					break;
				case gnut::QZS:
					par_ISB = t_gpar(crt_rec, t_gpar::QZS_ISB, ++ipar_num, "");
					par_ISB.value(0.0);
					par_ISB.setTime(_beg_time, _end_time);
					par_ISB.apriori(3000.0);
					_lsq->add_parameter(par_ISB);
					break;
				case gnut::GLO:
					par_ISB = t_gpar(crt_rec, t_gpar::GLO_ISB, ++ipar_num, "");
					par_ISB.value(0.0);
					par_ISB.setTime(_beg_time, _end_time);
					par_ISB.apriori(3000.0);
					_lsq->add_parameter(par_ISB);
					break;
				default:
					break;
				}
			}
		}
		return true;
	}
	
	bool t_glsqproc::_init_IF_amb_pars(t_glsq* lsq)
	{
		if (_rec_list.empty() ||
			_sat_list.empty())
		{
			string rec_empty("false");
			string sat_empty("false");
			if (_rec_list.empty()) rec_empty = "true";
			if (_sat_list.empty()) sat_empty = "true";
			write_log_info(_glog, 0, "ERROR", "can not init lsq-IF amb pars for  rec list empty is " + rec_empty + " and sat list empty is " + sat_empty);
			return false;
		}


		double sigAmb = dynamic_cast<t_gsetpar*>(_gset)->sigAmb();
		set<string>::const_iterator rec_beg = _rec_list.begin();
		set<string>::const_iterator rec_end = _rec_list.end();
		set<string>::const_iterator sat_beg = _sat_list.begin();
		set<string>::const_iterator sat_end = _sat_list.end();
		for (set<string>::const_iterator it_rec = rec_beg; it_rec != rec_end; ++it_rec)
		{
			for (set<string>::const_iterator it_sat = sat_beg; it_sat != sat_end; ++it_sat)
			{
				int ipar_num = _lsq->_x_solve.parNumber();
				t_gpar par_amb = t_gpar(*it_rec, t_gpar::AMB_IF, ipar_num + 1, *it_sat);
				par_amb.value(0.0);
				par_amb.setTime(_beg_time, _end_time);
				par_amb.apriori(sigAmb);
				lsq->add_parameter(par_amb);
			}
		}
		return true;
	}

	bool t_glsqproc::_init_IF_amb13_pars(t_glsq* lsq)
	{
		if (_rec_list.empty() ||
			_sat_list.empty())
		{
			string rec_empty("false");
			string sat_empty("false");
			if (_rec_list.empty()) rec_empty = "true";
			if (_sat_list.empty()) sat_empty = "true";
			write_log_info(_glog, 0, "ERROR", "can not init lsq-IF amb13 pars for  rec list empty is " + rec_empty + " and sat list empty is " + sat_empty);
			return false;
		}


		double sigAmb = dynamic_cast<t_gsetpar*>(_gset)->sigAmb();
		set<string>::const_iterator rec_beg = _rec_list.begin();
		set<string>::const_iterator rec_end = _rec_list.end();
		set<string>::const_iterator sat_beg = _sat_list.begin();
		set<string>::const_iterator sat_end = _sat_list.end();
		for (set<string>::const_iterator it_rec = rec_beg; it_rec != rec_end; ++it_rec)
		{
			for (set<string>::const_iterator it_sat = sat_beg; it_sat != sat_end; ++it_sat)
			{
				int ipar_num = _lsq->_x_solve.parNumber();
				t_gpar par_amb = t_gpar(*it_rec, t_gpar::AMB13_IF, ipar_num + 1, *it_sat);
				par_amb.value(0.0);
				par_amb.setTime(_beg_time, _end_time);
				par_amb.apriori(sigAmb);
				lsq->add_parameter(par_amb);
			}
		}
		return true;
	}

	bool t_glsqproc::_InitSLRFData(t_gslrf* _gslrf)
	{
		if (!_gslrf)
		{
			write_log_info(_glog, 1, "ERROR", "Init SLR data failed");
			return false;
		}
		int size_sitslr = _gslrf->getSize();
		for (int j = 0; j < size_sitslr; j++)
		{
			for (int jj = 0; jj < 3; jj++)
			{
				_gslrf->getOnegslrf(j)._xsta[jj] = _gslrf->getOnegslrf(j)._xsta[jj] + _gslrf->getOnegslrf(j)._xsta[jj + 3] * ((_beg_time.dmjd(false) + _end_time.dmjd(false)) / 2.0 - _gslrf->getOnegslrf(j).now.dmjd(false)) / 365.25;
			}
			xyz2ell(_gslrf->getOnegslrf(j)._xsta, _gslrf->getOnegslrf(j)._blh, false);
			_gslrf->getOnegslrf(j)._blh[2] = _gslrf->getOnegslrf(j)._blh[2] / 1000;
			for (int jj = 0; jj < 3; jj++)
			{
				_gslrf->getOnegslrf(j)._xsta[jj] = _gslrf->getOnegslrf(j)._xsta[jj] / 1000;    //m->km
			}

		}
		return true;
	}

	bool t_glsqproc::_initSLRProcPars(t_glsq* lsq)
	{
		if (!_gall_preproc)
		{
			write_log_info(_glog, 1, "ERROR", "ALLPREPROC data is empty.");
			return false;
		}
		_slr_list = _gall_preproc->getStaCode();
		_sat_sta = _gall_preproc->getSatSta();
		if (_slr_list.empty())
		{
			write_log_info(_glog, 1, "ERROR", "no slr station in POD_LSQ.");
			return false;
		}

		bool slr_valid = _init_slr_pars(lsq);
		bool rb_valid = _init_slr_rb(lsq);
		//bool rb_valid = true;
		bool erp_valid = _init_crs_erp_pars_slr(lsq);

		if (slr_valid || rb_valid || erp_valid)
		{
			write_log_info(_glog, 1, "INIT", "init the pars of slr estimated in LSQ.");
			return true;
		}
		else
		{
			write_log_info(_glog, 1, "ERROR", "can not init the pars of slr estimated in LSQ.");
			return false;
		}
	}

	bool t_glsqproc::_init_crs_erp_pars_slr(t_glsq* lsq)
	{
		int ipar_num = _lsq->_x_solve.parNumber();

		t_gtime erp_beg, erp_end;

		// init EOP parmeter
		// init the Eop dx dy vdx vdy

		// TODO:init value time sigma		
		map<string, double> data1, data2;

		double sig_xpole = dynamic_cast<t_gsetpar*>(_gset)->sigXpole();
		double sig_ypole = dynamic_cast<t_gsetpar*>(_gset)->sigYpole();
		double sig_dxpole = dynamic_cast<t_gsetpar*>(_gset)->sigDxpole();
		double sig_dypole = dynamic_cast<t_gsetpar*>(_gset)->sigDypole();
		double sig_ut1 = dynamic_cast<t_gsetpar*>(_gset)->sigUt1();
		double sig_dut1 = dynamic_cast<t_gsetpar*>(_gset)->sigDut1();

		erp_beg.from_mjd(_beg_time.mjd(), 0, 0.0);

		do
		{
			t_gpar par_xpole("", t_gpar::XPOLE, ++ipar_num, "", false);
			t_gpar par_ypole("", t_gpar::YPOLE, ++ipar_num, "", false);
			t_gpar par_dxpole("", t_gpar::DXPOLE, ++ipar_num, "", false);
			t_gpar par_dypole("", t_gpar::DYPOLE, ++ipar_num, "", false);
			t_gpar par_ut1("", t_gpar::UT1, ++ipar_num, "", false);
			t_gpar par_dut1("", t_gpar::DUT1, ++ipar_num, "", false);

			if (erp_beg + 86400 > _end_time)
			{
				erp_end = _end_time;
			}
			else
			{
				erp_end.from_mjd(erp_beg.mjd() + 1, 0, 0.0);
			}

			try
			{
				if (!_allerp) _allerp = dynamic_cast<t_gpoleut1*>((*_gall_proc)[t_gdata::ALLPOLEUT1]);
				if (!_allerp) return false;
				_allerp->getEopData(erp_beg, data1);
				_allerp->getEopData(erp_end, data2);
			}
			catch (exception e)
			{
				write_log_info(_glog, 0, "ERROR", e.what());
				return false;
			}

			par_xpole.setTime(erp_beg, erp_end);
			par_xpole.value(data1["XPOLE"]);
			par_xpole.apriori(sig_xpole);
			par_ypole.setTime(erp_beg, erp_end);
			par_ypole.value(data1["YPOLE"]);
			par_ypole.apriori(sig_ypole);
			par_dxpole.setTime(erp_beg, erp_end);
			par_dxpole.value(data2["XPOLE"] - data1["XPOLE"]);
			par_dxpole.apriori(sig_dxpole);
			par_dypole.setTime(erp_beg, erp_end);
			par_dypole.value(data2["YPOLE"] - data1["YPOLE"]);
			par_dypole.apriori(sig_dypole);
			par_ut1.setTime(erp_beg, erp_end);
			par_ut1.value(data1["UT1-TAI"]);
			par_ut1.apriori(sig_ut1);
			par_dut1.setTime(erp_beg, erp_end);
			par_dut1.value(data2["UT1-TAI"] - data1["UT1-TAI"]);
			//par_dut1.value(data1["LOD"]);
			par_dut1.apriori(sig_dut1);

			// add the Estimator
			lsq->add_parameter(par_xpole);
			lsq->add_parameter(par_ypole);
			lsq->add_parameter(par_dxpole);
			lsq->add_parameter(par_dypole);
			lsq->add_parameter(par_ut1);
			lsq->add_parameter(par_dut1);

			erp_beg = erp_beg + 86400;

		} while (erp_beg < _end_time);

		return true;
	}

	bool t_glsqproc::_init_slr_rb(t_glsq* lsq)
	{
		int ipar_num = _lsq->_x_solve.parNumber();
		int ista;

		if (_sat_sta.empty())
		{
			return false;
		}

		map<string, map<string,double>> RB;
		map<string, map<string,double>>* std;
		RB = _gall_rb->getrb();
		std = _gall_rb->getstd();
		double range_bias;

		auto it_beg = _sat_sta.begin();
		auto it_end = _sat_sta.end();
		double sigrb;

		for (auto it = it_beg; it != it_end; it++)
		{
			for (auto itsta = it->second.begin(); itsta != it->second.end(); itsta++)
			{
				t_gpar par_rb(*itsta, t_gpar::RB, ++ipar_num, it->first);
				if (RB.empty())
				{
					range_bias = 0.0;
				}
				else
				{
					range_bias = (RB)[it->first][*itsta];
				}

				if (std->empty())
				{
					sigrb = dynamic_cast<t_gsetpar*>(_gset)->sigRB();
				}
				else
				{
					sigrb = (*std)[it->first][*itsta];
				}

				par_rb.value(range_bias);
				// set the inital value
				t_gtime endT = _end_time;
				par_rb.setTime(_beg_time, endT);
				par_rb.apriori(sigrb);
				lsq->add_parameter(par_rb);
			}
			
		}

		return true;
	}

	bool t_glsqproc::_init_slr_pars(t_glsq* lsq)
	{
		int ipar_num = _lsq->_x_solve.parNumber();
		int ista;
		//vector<int> core_station = { 7821,7839,7941,7825,8834,7840,7119,7810,7090,7841,7501,1886,7105,7124,7080,7407,1889,7827 };
		ColumnVector x_f(0.0, 3), v_f(0.0, 3);
		auto it_beg = _slr_list.begin();
		auto it_end = _slr_list.end();
		t_gtriple rec_std = dynamic_cast<t_gsetpar*>(_gset)->sigRecPosTriple("SLR");
		for (auto it = it_beg; it != it_end; it++)
		{
			for (ista = 0; ista < _gdata_slrf->getSize(); ista++)
			{
				if (it->first == _gdata_slrf->getCode(ista))
				{
					break;
				}
			}
			if (ista >= _gdata_slrf->getSize())
			{
				write_log_info(_glog, 1, "NOTE", "The code is not found in SLRF2008");
				return false;
			}

			// init xyz 
			t_gpar par_x(it->second, t_gpar::CRD_X, ++ipar_num, "");
			t_gpar par_y(it->second, t_gpar::CRD_Y, ++ipar_num, "");
			t_gpar par_z(it->second, t_gpar::CRD_Z, ++ipar_num, "");

			for (int m = 0; m < 3; m++)
			{
				double ttmp = _gdata_slrf->getOnegslrf(ista)._xsta[m] + _gdata_slrf->getOnegslrf(ista)._offset[m] * 0.001;
				x_f(m + 1) = _gdata_slrf->getOnegslrf(ista)._xsta[m] + _gdata_slrf->getOnegslrf(ista)._offset[m] * 0.001;
				x_f(m + 1) = x_f(m + 1)*1000.0;
				v_f(m + 1) = 0.0;
			}

			// set the parameter tim
			par_x.value(x_f(1));
			par_y.value(x_f(2));
			par_z.value(x_f(3));

			// set the inital value
			t_gtime endT = _end_time;

			par_x.setTime(_beg_time, endT);
			par_y.setTime(_beg_time, endT);
			par_z.setTime(_beg_time, endT);

			//auto ret = find(core_station.begin(), core_station.end(), it->first);        
			//if (ret == core_station.end())
			//{
			//	// set the apriori  [crd read from rinex]
			//	par_x.apriori(rec_std[0]);
			//	par_y.apriori(rec_std[1]);
			//	par_z.apriori(rec_std[2]);
			//}
			//else
			//{
			//	par_x.apriori(1e-5);
			//	par_y.apriori(1e-5);
			//	par_z.apriori(1e-5);
			//}

			par_x.apriori(rec_std[0]);
			par_y.apriori(rec_std[1]);
			par_z.apriori(rec_std[2]);

			//  add parameters
			lsq->add_parameter(par_x);
			lsq->add_parameter(par_y);
			lsq->add_parameter(par_z);

		}

		return true;
	}
	
	bool t_glsqproc::_initLsqProdData(t_gallprod* data)
	{
		return false;
	}

	bool t_glsqproc::_initLsqProcPars(t_glsq* lsq)
	{
		return false;
	}

	bool t_glsqproc::_initLsqProcData(t_gallproc* data)
	{
		return false;
	}
	
	bool t_glsqproc::_select_rec_obs(const string& rec, const t_gtime& epoch, vector<t_gsatdata>& all_obs, vector<t_gsatdata>& rec_obs)
	{
		int  index_cout = -1;                        // add for PPP-RTK
		map<string, map <string, int> > valid_index; // add for PPP-RTK
		
		if (!_slip12)
		{
			_glog->logError("t_glsqproc", "_select_rec_obs", "_slip 12 is nullptr");
			return false;
		}

		if ((_frequency == 3) && !_slip13)
		{
			_glog->logError("t_glsqproc", "_select_rec_obs", "_slip 13 is nullptr");
			return false;
		}

		if ((_frequency == 5) && (!_slip13 || !_slip14 || !_slip15))
		{
			_glog->logError("t_glsqproc", "_select_rec_obs", "_slip 13/14/15 is nullptr");
			return false;
		}

		// get all the obs for the rec
		for (auto& iter : all_obs)
		{
			index_cout++;                            // add for PPP-RTK
			const string& obs_rec = iter.site();
			const string& obs_sat = iter.sat();

			if (obs_rec != rec) continue;

			// xiongyun add for realtime
			bool obs_valid = true;
			if (/*_realtime*/_lite_turboedit)
			{
			    // fix by glfeng
				//map<GOBS, int> slip = iter.slip();
				//for (auto it = slip.begin(); it != slip.end(); it++)
				//{
				//	if (it->second != 0)
				//	{
				//		obs_valid = false;
				//		break;
				//	}
				//}
				// 对于使用lli模式的周跳探测，需要手动设置, glfeng
				set<GOBSBAND> band_avail = iter.band_avail();
				bool frq1_ok = false;
				for (auto itband : band_avail)
				{
					FREQ_SEQ frq2 = _freq_index[iter.gsys()][itband];
					switch (frq2)
					{
					case FREQ_SEQ::FREQ_1:
						frq1_ok = true;
						break;
					case FREQ_SEQ::FREQ_2:
						iter.tb12(true);
						break;
					case FREQ_SEQ::FREQ_3:
						iter.tb13(true);
						break;
					case FREQ_SEQ::FREQ_4:
						iter.tb14(true);
						break;
					case FREQ_SEQ::FREQ_5:
						iter.tb15(true);
						break;
					default:
						break;
					}
				}
				if (!frq1_ok)
				{
					iter.tb12(false);
					iter.tb13(false);
					iter.tb14(false);
					iter.tb15(false);
				}
				
			}
			else
			{
				obs_valid = _slip12->use_of_obs(obs_rec, obs_sat, epoch);
			}
			if (obs_valid)
			{
				rec_obs.push_back(iter);
				valid_index[obs_rec][obs_sat] = index_cout;      // add for PPP-RTK
			}
		}

		if (rec_obs.empty())
		{
			_glog->logDebug("t_glsqproc", "_select_rec_obs", "This epoch have no useful data : " + rec);
			return false;
		}

		if (/*_realtime*/_lite_turboedit)
		{
			return true;
		}

		// set tb log in satdata
		for (auto& iter : rec_obs)
		{
			const string& obs_rec = iter.site();
			const string& obs_sat = iter.sat();

			if(_use_carrier_range) _apply_carrier_range(obs_rec, obs_sat, epoch, iter);
			// flag to marking useing of obs in tb
			bool obs_12 = _slip12->use_of_obs(obs_rec, obs_sat, epoch);
			iter.tb12(obs_12);
			all_obs[valid_index[obs_rec][obs_sat]].tb12(obs_12);              // add for PPP-RTK

			if (_frequency < 3) continue;
			bool obs_13 = _slip13->use_of_obs(obs_rec, obs_sat, epoch);
			iter.tb13(obs_13);
			all_obs[valid_index[obs_rec][obs_sat]].tb13(obs_13);              // add for PPP-RTK

			if (_frequency < 5) continue;
			bool obs_14 = _slip14->use_of_obs(obs_rec, obs_sat, epoch);
			iter.tb14(obs_14);
			all_obs[valid_index[obs_rec][obs_sat]].tb14(obs_14);              // add for PPP-RTK

			bool obs_15 = _slip15->use_of_obs(obs_rec, obs_sat, epoch);
			iter.tb15(obs_15);
			all_obs[valid_index[obs_rec][obs_sat]].tb15(obs_15);              // add for PPP-RTK
		}
		return true;
	}

	bool t_glsqproc::_apply_carrier_range(const string& rec, const string& sat, const t_gtime& epoch, t_gsatdata& obs)
	{
		double c1 = 0;
		double c2 = 0;
		double c3 = 0;
		GSYS gsys = obs.gsys();
		if (_slip12->is_carrier_range(rec, sat, epoch, c1, c2)) {
			GOBS gobsL1 = obs.select_phase(_band_index[gsys][FREQ_1]);
			GOBS gobsL2 = obs.select_phase(_band_index[gsys][FREQ_2]);
			double value1 = obs.getobs(gobsL1) - c1;
			double value2 = obs.getobs(gobsL2) - c2;
			obs.resetobs(gobsL1, value1);
			obs.resetobs(gobsL2, value2);
			obs.is_carrier_range(_band_index[gsys][FREQ_1], true);
			obs.is_carrier_range(_band_index[gsys][FREQ_2], true);
		}
		if (_frequency > 2 && _slip13) {
			if (_slip13->is_carrier_range(rec, sat, epoch, c1, c3)) {
				GOBS gobsL3 = obs.select_phase(_band_index[gsys][FREQ_3]);
				double value3 = obs.getobs(gobsL3) - c3;
				obs.resetobs(gobsL3, value3);
				obs.is_carrier_range(_band_index[gsys][FREQ_3], true);
			}
		}
		// 目前只考虑到3频
		return true;
	}

	bool t_glsqproc::_processOneEpoch(const t_gtime& crt_epoch, std::vector<t_gsatdata>& crt_obs)
	{
		// check obs size()
		
		if (crt_obs.empty())
		{
			_glog->logError("t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod("crt obs is empty"));
			return false;
		}
		//check glonass id--lijie
		for (auto &obs_tmp : crt_obs)
		{
			if (t_gsys::sat2gsys(obs_tmp.sat()) == GSYS::GLO)
			{
				if (obs_tmp.channel() == 255)
				{
					if(!_gdata_satinfo)
					{
						_glog->logInfo("t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod("no useful frequency id for: " + obs_tmp.sat()));
						continue;
					}
					else
					{
						obs_tmp.channel(_gdata_satinfo->getGloFid(obs_tmp.sat(), crt_epoch));
					}
				}
			}
		}
		//map<GSYS, vector<GOBSBAND> > tmp;
		bool updata_valid = _lsq->update_parameter(crt_epoch, crt_obs, false);
		if (!updata_valid)
		{
			_glog->logError("t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod("update_parameter failed"));
			return false;
		}
		//if (_lsq->_x_solve.amb_prns().size() <= 3)return false;

		// init the value of SION or VION
		// only used in pod all
		if (_gdata_ion && _useIon_cmd)
		{
			_add_sion_pars_value(crt_epoch, crt_epoch, _lsq, _gdata_ion);
		}

		// loop while recs
		for (const auto& crt_rec : _rec_list)
		{
			// claer equ before
			_crt_equ.clear_allequ();
			//_crt_equ.swap_allequ();

			// reset init value of coord parameters for kin mode
			if (_crd_est == KIN && _kin_xyz[_crt_rec].find(crt_epoch) != _kin_xyz[_crt_rec].end())
			{
				int idx_x = _lsq->_x_solve.getParam(crt_rec, t_gpar::CRD_X, "");
				int idx_y = _lsq->_x_solve.getParam(crt_rec, t_gpar::CRD_Y, "");
				int idx_z = _lsq->_x_solve.getParam(crt_rec, t_gpar::CRD_Z, "");

				_lsq->_x_solve.setParValue(idx_x, _kin_xyz[crt_rec][crt_epoch].crd(0));
				_lsq->_x_solve.setParValue(idx_y, _kin_xyz[crt_rec][crt_epoch].crd(1));
				_lsq->_x_solve.setParValue(idx_z, _kin_xyz[crt_rec][crt_epoch].crd(2));
			}

			// get all the obs for the rec
			vector<t_gsatdata> crt_rec_obs;
			bool select_valid = _select_rec_obs(crt_rec, crt_epoch, crt_obs, crt_rec_obs);
			if (!select_valid)
			{
				_glog->logDebug("t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod("no useful data : " + crt_rec));
				continue;
			}

			// Process this site and get the equations
			bool proc_site = _processOneRec(crt_epoch, crt_rec, crt_rec_obs);
			if (!proc_site)
			{
				_glog->logInfo("t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod("no useful equations : " + crt_rec));
				continue;
			}

			// add the new equations
			try
			{
				_crt_equ.set_newamb(dynamic_cast<t_gturboedit*>(_slip12.get()));
				_lsq->add_equation(_crt_equ, crt_epoch);
			}
			catch (...)
			{
				_glog->comment(t_glog::LOG_LV::LOG_ERROR, "t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod("add_equation throw error"));
				return false;
			}

			// optional function
			if (_gdata_ambcon && (_ambfix_cmd || _ambfix_xml))
			{
				_glog->logInfo(" t_glsqproc", "_processOneEpoch", "add ambcon for epoch " + crt_epoch.str_mjdsod());
				bool ambcon_valid = false;
				if (_obs_mode == RAW_ALL)
				{
					ambcon_valid = _pod_add_RAW_CB_Dd_ambcon(crt_epoch, _lsq, _gdata_ambcon, true);
					ambcon_valid = ambcon_valid && _pod_add_RAW_CB_Dd_ambcon(crt_epoch, _lsq, _gdata_ambcon, false);
				}
				else if(_obs_mode == IONO_FREE)
				{
					ambcon_valid = _pod_add_Dd_ambcon(this->_lsq, _gdata_ambcon);
				}

				if (!ambcon_valid)
				{
					_glog->comment(t_glog::LOG_LV::LOG_ERROR, "t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod(": no ambcon added"));
				}
			}
		}


		return true;
	}

	bool t_glsqproc::_processOneRec(const t_gtime& crt_epoch, const std::string& crt_rec, std::vector<t_gsatdata>& crt_obs)
	{
		if (crt_obs.empty() || crt_rec.empty())
		{
			_glog->logInfo("t_glsqproc", "_processOneRec", crt_epoch.str_mjdsod("have no useful data : " + crt_rec));
			return false;
		}

		// save as mark
		_crt_rec = crt_rec;

		if (crt_epoch.mjd() == 58610 && crt_epoch.sod() == 27540)
		{
			cout << "anb";
		}

		// loop for Combine equation
		vector<double> codeOmc;
		bool while_valid = false;
		do
		{
			// process one rec one sat and get equ
			// ========================================================================================================
			_crt_equ.clear_allequ();

			std::vector<t_gsatdata>::iterator it_beg = crt_obs.begin();
			std::vector<t_gsatdata>::iterator it_end = crt_obs.end();
			string crt_sat("NONE");
			for (auto it = it_beg; it != it_end; it++)
			{
				crt_sat = it->sat();
				if (_sat_list.empty() || _sat_list.find(crt_sat) == _sat_list.end())
				{
					continue;
				}
				bool proc_valid = _processOneRecOneSat(crt_epoch, crt_rec, crt_sat, *it);
				if (!proc_valid)
				{
					_glog->logInfo("t_glsqproc", "_processOneRec", crt_epoch.str_mjdsod("can not process the rec : " + crt_rec + " and the sat : " + crt_sat + " "));
					continue;
				}
			}

			// glfeng
			if (_is_prefix_res)
			{
				if (this->_prefix_residuals(crt_obs))
				{
					while_valid = true;
					continue;
				}
			}

			codeOmc = _crt_equ.get_codeomc();
			// =======================================================================================================
			/*if (_frequency >= 3 && _obs_mode == IONO_FREE)
			{
					vector<double> tmp;
					for (int i = 0; i < codeOmc.size();)
					{
						tmp.push_back(codeOmc[i]);
						i = i + 2;
					}
					codeOmc.clear();
					codeOmc = tmp;
			}*/

			//reset recclk omc, check range
			// =======================================================================================================
			int idx = _lsq->_x_solve.getParam(crt_rec, t_gpar::CLK, "");
			double recclk = 0.0;
			double sigma = 0.0;
			int k = Check_Range(codeOmc, recclk, sigma);
			
			// fix according LX
			if (idx >= 0 && !while_valid && (fabs(recclk / CLIGHT) > _recclk_threshold/*1e-8*/ ||
				(k >= 1 && fabs(_lsq->_x_solve[idx].value() / CLIGHT) < 1e-15)))
			{
				// reset recclk
				recclk += _lsq->_x_solve[idx].value();
				_lsq->_x_solve[idx].value(recclk);
				_bias_model->update_obj_clk(crt_rec, crt_epoch, recclk / CLIGHT);

				codeOmc.clear();
				while_valid = true;
			}
			else
			{
				while_valid = false;
			}

			// =======================================================================================================
		} while (while_valid);

		// give the satclk initial value
		if (!_init_sat_clk())
		{
			_glog->logInfo("t_glsqproc", "_processOneRec", crt_epoch.str_mjdsod("_init_sat_clk failed" + _crt_sat));
			return false;
		}

		return true;
	}

	bool t_glsqproc::_processOneRecOneSat(const t_gtime& crt_epoch, const std::string& rec, const std::string& sat, t_gsatdata& crt_obs)
	{
		t_gtime epoch = crt_epoch;
		bool code_valid = _base_model->cmb_equ(epoch, _lsq->_x_solve, crt_obs, _crt_equ);
		if (!code_valid)
		{
			_glog->logInfo("t_glsqproc", "_processOneRecOneSat", crt_epoch.str_mjdsod("can not combine equ for code in "));
			return false;
		}
		return true;
	}

	bool t_glsqproc::_init_sat_clk_IF(vector<t_gsatdata>& crt_obs, ColumnVector& l, t_glsq* lsq)
	{
		if (!_lsq || crt_obs.empty()) return false;
		for (int i = 0; i < crt_obs.size(); i++)
		{
			string sat = crt_obs[i].sat();
			int idx_satclk = _lsq->_x_solve.getParam("", t_gpar::CLK_SAT, sat);
			if (idx_satclk >= 0 && _lsq->_x_solve[idx_satclk].value() == 0.0)
			{
				try
				{
					// satclk delay = - code_omc
					double satclk = -(l(2 * i + 1) + 1.0);
					// code omc
					l(2 * i + 1) = l(2 * i + 1) + satclk;
					// phase omc
					l(2 * i + 2) = l(2 * i + 2) + satclk;
					// clk par init
					_lsq->_x_solve[idx_satclk].value(satclk);

				}
				catch (exception e)
				{
					write_log_info(_glog, 1, "ERROR", e.what());
					return false;
				}
			}
		}

		return true;
	}

	bool t_glsqproc::_init_sat_clk_ALL(vector<t_gsatdata>& crt_obs, ColumnVector& l, t_glsq* lsq)
	{
		return true;
		//return false;
	}

	bool t_glsqproc::_init_sat_clk()
	{
		set<string> sat_list = _crt_equ.get_sitelist(_crt_rec);
		for (auto sat_iter = sat_list.begin(); sat_iter != sat_list.end(); sat_iter++)
		{
			// get satclk par idx
			int idx_satclk = _lsq->_x_solve.getParam("", t_gpar::CLK_SAT, *sat_iter);
			if (idx_satclk >= 0 && double_eq(_lsq->_x_solve[idx_satclk].value(),0.0))
			{
				try
				{
					// get all phase code observ equations;
					vector<int> equ_list = _crt_equ.find_equ(_crt_rec, *sat_iter);
					double satclk = 0.0;
					// get satclk delay initial value 
					// TODO: use average function(now use the first code omc as the initial value)
					for (int idx : equ_list)
					{
						if (_crt_equ.get_obscombtype(idx).is_code())
						{
							satclk = -(_crt_equ.l[idx] + 1.0);
							break;
						}
					}
					// clk par init
					_lsq->_x_solve[idx_satclk].value(satclk);
					_bias_model->update_obj_clk(*sat_iter, _crt_time, satclk/CLIGHT);

					// omc = omc + satclk;
					for (int idx : equ_list)
					{
						_crt_equ.l[idx] += satclk;
					}
				}
				catch (exception e)
				{
					write_log_info(_glog, 1, "ERROR", e.what());
					return false;
				}
			}
		}
		return true;
	}

	bool t_glsqproc::_init_xml_settings(t_gsetbase* set)
	{
		// for all values with time
		_beg_time = dynamic_cast<t_gsetgen*>(set)->beg();
		_end_time = dynamic_cast<t_gsetgen*>(set)->end();
		_obs_intv = dynamic_cast<t_gsetgen*>(set)->sampling();

		// for all values with mode
		_obs_mode = dynamic_cast<t_gsetproc*>(set)->obs_combin();
		_lsq_mode = dynamic_cast<t_gsetproc*>(set)->lsq_mode();
		_fix_mode = dynamic_cast<t_gsetamb*>(set)->fix_mode();
		_upd_mode = dynamic_cast<t_gsetamb*>(set)->upd_mode();

		_slip_model = dynamic_cast<t_gsetproc*>(_gset)->slip_model();
		_iono_order = dynamic_cast<t_gsetproc*>(_gset)->iono_order();
		_ztd_model = dynamic_cast<t_gsetproc*>(_gset)->ztd_model(&_pwc_intv);
		_crd_est = dynamic_cast<t_gsetproc*>(_gset)->crd_est();
		_sysbias_model = dynamic_cast<t_gsetproc*>(_gset)->sysbias_model();
		_ifb_model = dynamic_cast<t_gsetproc*>(_gset)->ifb_model();
		_frequency = dynamic_cast<t_gsetproc*>(set)->frequency();
		_ambfix_xml = dynamic_cast<t_gsetproc*>(_gset)->ambfix();
		_tropo_grad = dynamic_cast<t_gsetproc*>(_gset)->tropo_grad();    // this two line add by glfeng
		if (_tropo_grad) _grd_model = dynamic_cast<t_gsetproc*>(_gset)->grd_model(&_grd_pwc_intv);
		

		// for all settings with sat and rec
		_rec_list = dynamic_cast<t_gsetgen*>(set)->rec_all();
		_leo_list = dynamic_cast<t_gsetgen*>(set)->rec_leo();
		_sat_list = dynamic_cast<t_gsetgnss*>(set)->sat();
		_sys_list = dynamic_cast<t_gsetgen*>(set)->sys();
		_rover_list = dynamic_cast<t_gsetgen*>(set)->list_rover();

		// 这两种模式都是直接对lli进行判断 _lite_turboedit是为了那些非实时，又只想采用GF/MW/Gap 进行判断的场景
		_realtime = dynamic_cast<t_gsetproc*>(set)->realtime();
		_lite_turboedit = dynamic_cast<t_gsetturboedit*>(set)->liteMode();
		if (_realtime || _slip_model == SLIPMODEL::DEF_DETECT_MODEL) _lite_turboedit = true;

		//POD mode
		_mode = dynamic_cast<t_gsetgen*>(set)->mode();
		_slr = bool(dynamic_cast<t_gsetgen*>(set)->slr());

		_num_threads = dynamic_cast<t_gsetproc*>(set)->num_threads();
#ifdef USE_OPENMP
		omp_set_num_threads(_num_threads);
#endif

		// add by glfeng
		_is_prefix_res = dynamic_cast<t_gsetproc*>(set)->pre_fix_residuals();
		// add by xiongyun
		_IRC = dynamic_cast<t_gsetproc*>(set)->IRC();
		_carrier_range = dynamic_cast<t_gsetproc*>(set)->carrier_range();
		_use_carrier_range = dynamic_cast<t_gsetproc*>(set)->apply_carrier_range();
		_ultrasp3 = dynamic_cast<t_gsetproc*>(set)->ultrasp3(); // only used in realtime
		_maxres_norm = dynamic_cast<t_gsetproc*>(set)->max_res_norm();
		sys_support = dynamic_cast<t_gsetgnss*>(set)->SYS_SUPPORTED();
		for (auto it = sys_support.begin(); it != sys_support.end(); it++) {
			_maxres_C[*it] = dynamic_cast<t_gsetgnss*>(set)->maxres_C(*it);
			_maxres_L[*it] = dynamic_cast<t_gsetgnss*>(set)->maxres_L(*it);
		}
		_band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::GPS);
		_band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::GAL);
		_band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::GLO);
		_band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::BDS);
		_band_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::QZS);

		_freq_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(set)->freq_index(gnut::GPS);
		_freq_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(set)->freq_index(gnut::GAL);
		_freq_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(set)->freq_index(gnut::GLO);
		_freq_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(set)->freq_index(gnut::BDS);
		_freq_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(set)->freq_index(gnut::QZS);

		// 比较gen的sys，剔除sat list
		std::set<string> tmp_sat_list;
		for (const auto& sat : _sat_list)
		{
			auto gsys = t_gsys::sat2gsys(sat);
			auto ssys = t_gsys::gsys2str(gsys);
			if (_sys_list.count(ssys) == 0) continue;
			tmp_sat_list.insert(sat);
		}
		_sat_list = tmp_sat_list;

		// 如果为空，按照系统补充卫星
		if (_sat_list.empty())
		{
			for (auto sys_iter = _sys_list.begin(); sys_iter != _sys_list.end(); sys_iter++)
			{
				auto sats_temp = dynamic_cast<t_gsetgnss*>(set)->sat(t_gsys::str2gsys(*sys_iter));
				_sat_list.insert(sats_temp.begin(), sats_temp.end());
			}
		}

		// 根据sat_rm删除
		std::set<string> sat_rm = dynamic_cast<t_gsetgen*>(set)->sat_rm();
		for (auto iter : sat_rm) { _sat_list.erase(iter); }


		// Get site coordinates from XML/rinex O files
		std::set<string> rec_rm;
		for (auto it : _rec_list)
		{
			if (_leo_list.find(it) != _leo_list.end())  // ignore LEO satellites
				continue;
			_rec_crds[it] = dynamic_cast<t_gsetrec*>(_gset)->get_crd_xyz(it);
			_rec_stds[it] = dynamic_cast<t_gsetrec*>(_gset)->get_std_xyz(it);
			if (_rec_crds[it].zero())
			{
				rec_rm.insert(it);
			}
		}
		if (_crd_est != FIX) {
			for (auto iter : rec_rm) { _rec_list.erase(iter); }
		}



		return true;
	}

	bool t_glsqproc::_init_lsq_settings(t_gsetbase* set, t_gallproc* data, t_glog* log)
	{
		// cycle slip
		if (_frequency == 2)
		{
			_slip12 = shared_ptr<t_gcycleslip>(new t_gturboedit(set, log));
		}
		else if (_frequency == 3)
		{
			_slip12 = shared_ptr<t_gcycleslip>(new t_gturboedit(set, log));
			_slip13 = shared_ptr<t_gcycleslip>(new t_gturboedit(set, log, 3));
		}
		else if (_frequency == 4)
		{
			_slip12 = shared_ptr<t_gcycleslip>(new t_gturboedit(set, log));
			_slip13 = shared_ptr<t_gcycleslip>(new t_gturboedit(set, log, 3));
			_slip14 = shared_ptr<t_gcycleslip>(new t_gturboedit(set, log, 4));
		}
		else if (_frequency == 5)
		{
			_slip12 = shared_ptr<t_gcycleslip>(new t_gturboedit(set, log));
			_slip13 = shared_ptr<t_gcycleslip>(new t_gturboedit(set, log, 3));
			_slip14 = shared_ptr<t_gcycleslip>(new t_gturboedit(set, log, 4));
			_slip15 = shared_ptr<t_gcycleslip>(new t_gturboedit(set, log, 5));
		}
		else
		{
			throw logic_error("you should set the frequency in xml as 1/2/3/4/5.");
		}

		// Nothing but cout a string here
		// Init Paramter in Constructor, make sure lsq is new class
		if (_lsq)
		{
			delete _lsq;
			_lsq = nullptr;
		}
		_lsq = new t_glsq(_gset);

		// Nothing but cout a string here
		if (!_IRC)
		{
			_bias_model = shared_ptr<t_gbiasmodel>(new t_gprecisebias(data, log, set));
		}
		else
		{
			_bias_model = shared_ptr<t_gbiasmodel>(new t_gIRCbias(data, log, set));
		}
		
		if (_frequency == 2 && _obs_mode == IONO_FREE)
		{
			if (!_IRC)
			{
				//shared_ptr<t_gupdatepar> updateIF(new t_gupdateparIF(_slip12, set));
				shared_ptr<t_gupdatepar> updateIF(new t_gupdateparIF(_slip12, _band_index));
				updateIF->set_amb_update_way(/*_realtime*/_lite_turboedit);
				updateIF->set_interval(_obs_intv);
				_lsq->set_update_par(updateIF);

				_base_model = shared_ptr<t_gbasemodel>(new t_gcombIF(set, log, _bias_model, data));
			}
			else
			{
				shared_ptr<t_gupdatepar> updateIRC(new t_gupdateparIRC(_slip12));
				updateIRC->set_amb_update_way(/*_realtime*/_lite_turboedit);
				updateIRC->set_interval(_obs_intv);
				_lsq->set_update_par(updateIRC);
			
				_base_model = shared_ptr<t_gbasemodel>(new t_gcombIF(set, log, _bias_model, data));
			}
		}
		else if (_frequency >= 3 && _obs_mode == IONO_FREE)
		{
			shared_ptr<t_gupdatepar> update_if13(new t_gupdateparIF13(_slip12, _slip13, _band_index));
			update_if13->set_amb_update_way(/*_realtime*/_lite_turboedit);
			update_if13->set_interval(_obs_intv);
			_lsq->set_update_par(update_if13);

			_base_model = shared_ptr<t_gbasemodel>(new t_gcombIF13(set, log, _bias_model, data));
		}
		else if (_obs_mode == RAW_ALL || _obs_mode == RAW_MIX)
		{
			shared_ptr<t_gupdatepar> update_all(new t_gupdateparALL(_slip12, _slip13, _slip14, _slip15, _band_index, _iono_order));
			update_all->set_amb_update_way(/*_realtime*/_lite_turboedit);
			update_all->set_interval(_obs_intv);
			if (!_rover_list.empty()) update_all->set_sig_ion(dynamic_cast<t_gsetpar*>(_gset)->sigSion(_rover_list[0]));
			_lsq->set_update_par(update_all);

			_base_model = shared_ptr<t_gbasemodel>(new t_gcombALL(set, log, _bias_model, data));
		}
		else if (_frequency == 2 && _obs_mode == WL_COMBIN)
		{
			shared_ptr<t_gupdatepar> update_wl(new t_updateparWL(_slip12, _band_index));
			update_wl->set_amb_update_way(/*_realtime*/_lite_turboedit);
			update_wl->set_interval(_obs_intv);
			_lsq->set_update_par(update_wl);

			_base_model = shared_ptr<t_gbasemodel>(new t_gcombWL(set, log, _bias_model, data));
		}
		else
		{
			throw logic_error("Support IF[DOUBLE/TRIPLE]/UDUC[5 FREQ]");
		}


		return true;
	}

	/*void t_glsqproc::logBPL(const Matrix & B, const DiagonalMatrix & P, const ColumnVector & l)
	{
		Matrix tmp = B;
		for (int i = 1; i <= tmp.Nrows(); i++)
		{
			char row[4]; sprintf(row, "%4d", i); string rows(row);
			string tmp_str;
			for (int j = 1; j <= tmp.Ncols(); j++)
			{
				char col[4];        sprintf(col, "%4d", j); string cols(col);
				char tmp_array[20]; sprintf(tmp_array, "%18.8f", tmp(i, j));

				tmp_str = tmp_str + rows + cols + tmp_array;
			}
#ifdef DEBUG_LOG
			write_log_info(_B_log, 0, "DEBUG", tmp_str);
#endif
		}

		DiagonalMatrix tmp2 = P;
		for (int i = 1; i <= tmp2.Nrows(); i++)
		{
			char row[4]; sprintf(row, "%4d", i); string rows(row);
			string tmp_str;

			char col[4];        sprintf(col, "%4d", 1); string cols(col);
			char tmp_array[20]; sprintf(tmp_array, "%18.8f", tmp2(i));
			tmp_str = tmp_str + rows + cols + tmp_array;
#ifdef DEBUG_LOG
			write_log_info(_P_log, 0, "DEBUG", tmp_str);
#endif
		}

		ColumnVector tmp3 = l;
		for (int i = 1; i <= tmp3.Nrows(); i++)
		{
			char row[4]; sprintf(row, "%4d", i); string rows(row);
			string tmp_str;
			char col[4];        sprintf(col, "%4d", 1); string cols(col);
			char tmp_array[20]; sprintf(tmp_array, "%18.8f", tmp3(i));
			tmp_str = tmp_str + rows + cols + tmp_array;
#ifdef DEBUG_LOG
			write_log_info(_L_log, 0, "DEBUG", tmp_str);
#endif
		}
	}*/

	bool t_glsqproc::_extract_resfile(t_gallrecover& recover_data,string resfile)
	{

		// encode recover file
		t_gfile gout; 
		t_resfile grecover_coder(_gset);
		string recover_path = dynamic_cast<t_gsetout*>(_gset)->outputs("recover");
		if (recover_path.empty())
			recover_path = dynamic_cast<t_gsetout*>(_gset)->outputs("rcv");
		if (recover_path.empty())
			recover_path = resfile;
		if (recover_path.empty())
			recover_path = "resfile_temp_$(date)";

		substitute(recover_path, "$(date)", _beg_time.str_yyyydoy(), false);
		substitute(recover_path, "$(rec)", _crt_rec, false);
		//set I/O
		gout.path(recover_path); 
		gout.glog(_glog);
		//set coder
		grecover_coder.glog(_glog);
		grecover_coder.path(recover_path);
		grecover_coder.add_data("ID1", &recover_data);

		gout.coder(&grecover_coder);
		// write
		gout.run_write();

		return true;
	}

	bool t_glsqproc::_extract_clkfile(t_gallrecover& recover_data,t_gallprec::clk_type type)
	{
		// get allprec data
		t_gallprec clk_data;

		recover_data.get_clkdata(clk_data,type);

		// encoder data
		t_gfile gout;
		string clk_path = "";
		if (type == t_gallprec::AS) 
		{
			if(clk_path.empty()) clk_path = dynamic_cast<t_gsetout*>(_gset)->outputs("satclk");
			if(clk_path.empty()) clk_path = dynamic_cast<t_gsetout*>(_gset)->outputs("cas");
		}
		else if (type == t_gallprec::AR) 
		{
		    if (clk_path.empty()) clk_path = dynamic_cast<t_gsetout*>(_gset)->outputs("recclk");
			if (clk_path.empty()) clk_path = dynamic_cast<t_gsetout*>(_gset)->outputs("car");
		}
		
		if (clk_path.empty()) 
		{
			clk_path = "clkfiletemp_$(date)";
		}
		substitute(clk_path, "$(date)",_beg_time.str_yyyydoy(), false);
		substitute(clk_path, "$(rec)", _crt_rec, false);
		
		gout.path(clk_path);
		gout.glog(_glog);

		t_rinexc clk_coder(_gset);
		clk_coder.add_data("ID1", &clk_data);
		gout.coder(&clk_coder);

		gout.run_write();
		return true;
	}
	
	bool t_glsqproc::_extract_clkfile(t_gallrecover& recover_data, t_gallprec::clk_type type, bool isTriple)
	{
		if (!isTriple)
		{
			return _extract_clkfile(recover_data, type);
		}
		// get allprec data
		t_gallprec clk_data;

		recover_data.get_clk13data(clk_data, type);

		// encoder data
		t_gfile gout;
		string clk_path = "";
		if (type == t_gallprec::AS) {
			clk_path = dynamic_cast<t_gsetout*>(_gset)->outputs("satclk");
		}
		else if (type == t_gallprec::AR) {
			clk_path = dynamic_cast<t_gsetout*>(_gset)->outputs("recclk13");
		}
		if (clk_path.empty()) {
			clk_path = "clkfiletemp_$(date)";
		}
		substitute(clk_path, "$(date)", _beg_time.str_yyyydoy(), false);

		gout.path(clk_path);
		gout.glog(_glog);

		t_rinexc clk_coder(_gset);
		clk_coder.add_data("ID1", &clk_data);
		gout.coder(&clk_coder);

		gout.run_write();
		return true;
	}

	bool t_glsqproc::_extract_rbfile(t_gallrecover& recover_data)
	{
		// get allprec data
		t_gallrb rb_data;

		recover_data.get_rbdata(rb_data);

		_get_rb_std(rb_data);

		// encoder data
		t_gfile gout;
		string rb_path = "";
		
        rb_path = dynamic_cast<t_gsetout*>(_gset)->outputs("rb");
		
		if (rb_path.empty()) {
			rb_path = "RangeBias";
		}

		gout.path(rb_path);
		gout.glog(_glog);

		t_rb rb_coder(_gset);
		rb_coder.add_data("ID2", &rb_data);
		gout.coder(&rb_coder);

		gout.run_write();
		return true;
	}

	bool t_glsqproc::_extract_eopfile(t_gallrecover& recover_data)
	{

		t_galleop eop_data;

		recover_data.get_eopdata(eop_data);

		// encoder data
		t_gfile gout;
		string eop_path = "";

		eop_path = dynamic_cast<t_gsetout*>(_gset)->outputs("eop");

		if (eop_path.empty()) {
			eop_path = "EOP";
		}

		gout.path(eop_path);
		gout.glog(_glog);

		t_eop eop_coder(_gset);
		eop_coder.add_data("ID4", &eop_data);
		gout.coder(&eop_coder);

		gout.run_write();
		return true;
	}

	bool t_glsqproc::_extract_polefile(t_gallrecover& recover_data)
	{
		t_gpoleut1 pole;

		recover_data.get_poledata(pole);

		// encoder data
		t_gfile gout;
		string pole_path = "";

		pole_path = dynamic_cast<t_gsetout*>(_gset)->outputs("pole");

		if (pole_path.empty()) {
			pole_path = "poleut1_out";
		}

		gout.path(pole_path);
		gout.glog(_glog);

		t_poleut1 pole_coder(_gset);
		pole_coder.add_data("ID5", &pole);
		gout.coder(&pole_coder);

		gout.run_write();
		return true;

	}

	bool t_glsqproc::_extract_xmlcrd(t_gallrecover& recover_data)
	{
		// get allprec data
		t_gallsta sta_data;

		recover_data.get_stadata(sta_data);

		// encoder data
		t_gfile gout;
		string sta_path = "";

		sta_path = dynamic_cast<t_gsetout*>(_gset)->outputs("sta");

		if (sta_path.empty()) {
			sta_path = "StaCor";
		}

		gout.path(sta_path);
		gout.glog(_glog);

		t_sta sta_coder(_gset);
		sta_coder.add_data("ID3", &sta_data);
		gout.coder(&sta_coder);

		gout.run_write();

		xml_document xdoc;

		auto root = xdoc.load_file(_gset->xml_name().c_str());
		xdoc.save_file(string(_gset->xml_name() + t_gtime::current_time().str_hms()).c_str());

		auto config = xdoc.child("config");

		string obj = "SNX";
		// remove for case
		if (config.child("receiver")) config.remove_child("receiver");

		config.append_child("receiver");
		xml_node pars = config.child("receiver");

		pars.append_child(pugi::node_comment).set_value("The part for receiver settings. ");
		pars.append_child(pugi::node_comment).set_value("id    : for index, Rec name.     ");
		pars.append_child(pugi::node_comment).set_value("X     :  X pos.                   ");
		pars.append_child(pugi::node_comment).set_value("Y     :  Y pos.                   ");
		pars.append_child(pugi::node_comment).set_value("Z     :  Z pos.                   ");
		pars.append_child(pugi::node_comment).set_value("dX    : dX pos.                   ");
		pars.append_child(pugi::node_comment).set_value("dY    : dY pos.                   ");
		pars.append_child(pugi::node_comment).set_value("dZ    : dZ pos.                   ");

		map<string, sta_cor> rec_crds = sta_data.getsta();

		for (auto it : rec_crds)
		{

			string rec = it.first;

			t_gtriple rec_crd(it.second.dx, it.second.dy, it.second.dz);
			t_gtriple rec_std(0.001, 0.001, 0.001);

			xml_node par_sta = pars.append_child("rec");
			pugi::xml_attribute tmp_attri;
			tmp_attri = par_sta.append_attribute("id");
			tmp_attri.set_value(rec.c_str());

			string crd_x = format("%16.6lf", rec_crd[0]);
			tmp_attri = par_sta.append_attribute("X");
			tmp_attri.set_value(crd_x.c_str());

			string crd_y = format("%16.6lf", rec_crd[1]);
			tmp_attri = par_sta.append_attribute("Y");
			tmp_attri.set_value(crd_y.c_str());

			string crd_z = format("%16.6lf", rec_crd[2]);
			tmp_attri = par_sta.append_attribute("Z");
			tmp_attri.set_value(crd_z.c_str());

			crd_x = format("%16.6lf", rec_std[0]);
			tmp_attri = par_sta.append_attribute("dX");
			tmp_attri.set_value(crd_x.c_str());

			crd_y = format("%16.6lf", rec_std[1]);
			tmp_attri = par_sta.append_attribute("dY");
			tmp_attri.set_value(crd_y.c_str());

			crd_z = format("%16.6lf", rec_std[2]);
			tmp_attri = par_sta.append_attribute("dZ");
			tmp_attri.set_value(crd_z.c_str());

			tmp_attri = par_sta.append_attribute("obj");
			tmp_attri.set_value(obj.c_str());
		}

		xdoc.save_file(string(_gset->xml_name()).c_str());
		return true;

	}

	bool t_glsqproc::_extract_stafile(t_gallrecover& recover_data)
	{
		// get allprec data
		t_gallsta sta_data;

		recover_data.get_stadata(sta_data);

		// encoder data
		t_gfile gout;
		string sta_path = "";

		sta_path = dynamic_cast<t_gsetout*>(_gset)->outputs("sta");

		if (sta_path.empty()) {
			sta_path = "StaCor";
		}

		gout.path(sta_path);
		gout.glog(_glog);

		t_sta sta_coder(_gset);
		sta_coder.add_data("ID3", &sta_data);
		gout.coder(&sta_coder);

		gout.run_write();
		return true;
	}

	bool t_glsqproc::_extract_ambupd(t_gallambupd* allambupd)
	{
		// Encoder ambupd files
		string path = dynamic_cast<t_gsetout*>(_gset)->outputs("ambupd");
		set<string> ambupd_sites = allambupd->get_sites();
		for (const auto site : ambupd_sites) {
			if (path.empty()) path = "file://" + site + "_ambupd_" + _beg_time.str_yyyydoy();
			else substitute(path, "$(rec)", site, false);

			t_gallambupd* outambupd = new t_gallambupd();
			outambupd->addAmbUpd(site, allambupd->getOneSiteAmbUpd(site));
			t_gfile* gout = new t_gfile;
			gout->glog(_glog);
			gout->path(path);

			t_gcoder* gcoder = new t_ambupd(_gset, "", 4096);
			gcoder->clear();
			gcoder->glog(_glog);
			gcoder->path(path);
			gcoder->add_data("IDX", outambupd);

			gout->coder(gcoder);

			_glog->comment(1, "t_glsqproc::_extract_ambupd", " Writing AMBUPD file : " + path);
			gout->run_write();
			delete gout; gout = NULL;
			delete gcoder; gcoder = NULL;
			delete outambupd; outambupd = NULL;
		}
		delete allambupd; allambupd = NULL;

		return true;
	}

	bool t_glsqproc::_get_rb_std(t_gallrb& rb_data)
	{
		map<string, map<string, double>> RB = rb_data.getrb();
		for (auto iter = RB.begin(); iter != RB.end(); iter++)
		{
			for (auto itersite = iter->second.begin(); itersite != iter->second.end(); itersite++)
			{
				int irb = _lsq->_x_solve.getParam(itersite->first, t_gpar::RB, iter->first);
				double stdx = _lsq->stdx(irb);
				rb_data.addstd(itersite->first, iter->first, stdx);
			}
			
		}
		
		return true;
	}

	bool t_glsqproc::_init_rec_crd_pars(t_glsq* lsq, bool isPPP)
	{
		int ipar_num = _lsq->_x_solve.parNumber();

		//
		set<string> rec_tmp = _gall_obj->objects_set(t_gdata::REC);
		set<string> rec_final;
		set_intersection(rec_tmp.begin(), rec_tmp.end(), _rec_list.begin(), _rec_list.end(), insert_iterator<set<string>>(rec_final, rec_final.begin()));

		_rec_list = rec_final;
		if (_rec_list.empty())
		{
			write_log_info(_glog, 0, "ERROR : The rec list is empty !");
			return false;
		}

		auto it_beg = _rec_list.begin();
		auto it_end = _rec_list.end();
		for (auto it = it_beg; it != it_end; )
		{
			string site = *it;
			shared_ptr<t_gobj> _grec = _gall_obj->obj(site);

			// get the init value
			// @TODO : consider the crd and std is zeros

			t_gtriple rec_crd;
			t_gtriple rec_std;

			if (isPPP && _realtime != true)
			{
				_grec->get_recent_crd(_beg_time, 0.1, rec_crd, rec_std);
			}
			else
			{
				rec_crd = _rec_crds[site];
				rec_std = _rec_stds[site];
			}

			if (rec_std.zero() || isPPP)
			{
				rec_std = dynamic_cast<t_gsetpar*>(_gset)->sigRecPosTriple(site);
			}

			// jdhuang : est crd
			if (_crd_est == FIX) 
			{
				// only need sinex file 
				t_gtriple snx_crd, snx_std;
				if (!_grec->get_recent_crd(_beg_time,0.1,snx_crd,snx_std)) 
				{					
					// add for station no snx crd
					if (dynamic_cast<t_gsetproc*>(_gset)->realtime())//add by xiongyun for realtime
					{
						// add by xiongyun: to get site crd from resfile
						if (_siteres.size() != 0)
						{
							if (_siteres.find(str2upper(site)) != _siteres.end())
							{
								snx_crd[0] = _siteres[str2upper(site)]->get_recover_data_value(t_gpar::t_type::CRD_X);
								snx_crd[1] = _siteres[str2upper(site)]->get_recover_data_value(t_gpar::t_type::CRD_Y);
								snx_crd[2] = _siteres[str2upper(site)]->get_recover_data_value(t_gpar::t_type::CRD_Z);
								if (double_eq(snx_crd[0], 0.0) || double_eq(snx_crd[1], 0.0) || double_eq(snx_crd[2], 0.0))
								{
									_glog->comment(t_glog::LOG_LV::LOG_INFO, "glsqproc", "+init_rec_crd_pars", site + "'s crd is 0 in resfile, Erase the site!");
									it = _rec_list.erase(it);
								}
								else
								{
									snx_std = t_gtriple(0.0001, 0.0001, 0.0001);
									_grec->crd(snx_crd, snx_std, FIRST_TIME, LAST_TIME, true);
									it++;
								}
							}
							else
							{
								_glog->comment(t_glog::LOG_LV::LOG_INFO, "glsqproc", "+init_rec_crd_pars", site + " has no crd in resfile, Erase the site!");
								it = _rec_list.erase(it);
							}
						}
						else
						{
							_glog->comment(t_glog::LOG_LV::LOG_INFO, "glsqproc", "+init_rec_crd_pars", site + " has no crd in resfile, Erase the site!");
							it = _rec_list.erase(it);
						}
					}
					else
					{
						//DEBUG by ZHJ(Don't use the station with no snx crd）
					//	_glog->comment(t_glog::LOG_LV::LOG_INFO, "glsqproc", "+init_rec_crd_pars", site + "'s std precision is low, Erase the site!");
					//	it = _rec_list.erase(it);
					//	continue;

						snx_crd = dynamic_cast<t_gsetrec*>(_gset)->get_crd_xyz(site);
						snx_std = dynamic_cast<t_gsetrec*>(_gset)->get_std_xyz(site);
						if (snx_std[0] > 0.1 || snx_std[1] > 0.1 || snx_std[2] > 0.1)
						{
							_glog->comment(t_glog::LOG_LV::LOG_INFO, "glsqproc", "+init_rec_crd_pars", site + "'s std precision is low, Erase the site!");
							it = _rec_list.erase(it);
							continue;
						}
						if (double_eq(snx_crd[0], 0.0) || double_eq(snx_crd[1], 0.0) || double_eq(snx_crd[2], 0.0))
						{
							_glog->comment(t_glog::LOG_LV::LOG_INFO, "glsqproc", "+init_rec_crd_pars", site + "'s crd precision is absent, Erase the site!");
							it = _rec_list.erase(it);
						}
						else
						{
							snx_std = t_gtriple(0.0001, 0.0001, 0.0001);
							// change by glfeng make sure all the time the crd is the same;
							_grec->crd(snx_crd, snx_std, FIRST_TIME, LAST_TIME, true);
							it++;
						}
					}
				}
				else 
				{
					// change by ZHJ make sure all the time the crd is the same;
					_grec->crd(snx_crd, snx_std,FIRST_TIME,LAST_TIME,true);
					it++;
				}
				continue;
			}


			// init xyz 
			t_gpar par_x(site, t_gpar::CRD_X, ++ipar_num, "");
			t_gpar par_y(site, t_gpar::CRD_Y, ++ipar_num, "");
			t_gpar par_z(site, t_gpar::CRD_Z, ++ipar_num, "");

			// set the parameter tim
			par_x.value(rec_crd[0]);
			par_y.value(rec_crd[1]);
			par_z.value(rec_crd[2]);

			// set the inital value
			t_gtime endT = _end_time;
			if (_crd_est == CONSTRPAR::KIN || _crd_est == CONSTRPAR::SIMU_KIN) 
			{
				endT = _beg_time;
			}
			par_x.setTime(_beg_time, endT);
			par_y.setTime(_beg_time, endT);
			par_z.setTime(_beg_time, endT);

			// set the apriori  [crd read from rinex]
			par_x.apriori(rec_std[0]);
			par_y.apriori(rec_std[1]);
			par_z.apriori(rec_std[2]);

			//  add parameters
			lsq->add_parameter(par_x);
			lsq->add_parameter(par_y);
			lsq->add_parameter(par_z);

			it++;
		}
		if (_leo_list.size() != 0)
		{
			_rec_list.insert(_leo_list.begin(), _leo_list.end());
		}
		return true;
	}

	bool t_glsqproc::_get_site_crd(string site, t_gtriple& site_crd)
	{
		shared_ptr<t_gobj> _grec = _gall_obj->obj(site);
		t_gtriple snx_std;
		if (!_grec->get_recent_crd(_beg_time, 0.1, site_crd, snx_std))
		{
			// add by xiongyun: to get site crd from resfile
			if (_siteres.size() != 0)
			{
				if (_siteres.find(str2upper(site)) != _siteres.end())
				{
					site_crd[0] = _siteres[str2upper(site)]->get_recover_data_value(t_gpar::t_type::CRD_X);
					site_crd[1] = _siteres[str2upper(site)]->get_recover_data_value(t_gpar::t_type::CRD_Y);
					site_crd[2] = _siteres[str2upper(site)]->get_recover_data_value(t_gpar::t_type::CRD_Z);
					if (double_eq(site_crd[0], 0.0) || double_eq(site_crd[1], 0.0) || double_eq(site_crd[2], 0.0))
					{
						return false;
					}

					return true;
				}

				return false;
			}

			return false;
		}

		return true;
	}


	void t_glsqproc::_read_kinfile()
	{
		vector<string> kin_list = dynamic_cast<t_gsetinp*>(_gset)->inputs("kin");

		if (kin_list.empty())
		{
			for (auto crt_rec : _rec_list)
			{
				string kin_site;
				transform(crt_rec.begin(), crt_rec.end(), back_inserter(kin_site), ::tolower);
				kin_site =  "kin_" + _beg_time.str_yyyydoy() + "_" + kin_site;
				kin_list.push_back(kin_site);
			}
		}

		fstream kinfile;
		for (auto crt_site_file : kin_list)
		{
			string site_upper;
			string site_low = crt_site_file.substr(crt_site_file.find_last_of("_") + 1);
			if (IS_SUPPORTED_LEO(site_low))
			{
				site_low = strl2strs(site_low);
			}
			transform(site_low.begin(), site_low.end(), back_inserter(site_upper), ::toupper);

			if (crt_site_file.find(GFILE_PREFIX) != string::npos) crt_site_file = crt_site_file.substr(7);
			kinfile.open(crt_site_file, ios::in);

			string temp;
			t_gtime time_t;
			while (getline(kinfile, temp)) {
				if (temp.substr(0, 3) == "*  ") {
					string ymd = temp.substr(3, 28);
					stringstream time_inp(ymd);
					int yr, m, d, h, min;
					double s;
					time_inp >> yr >> m >> d >> h >> min >> s;
					time_t.from_ymdhms(yr, m, d, h, min, s);
					continue;
				}
				if (temp.substr(0, 4) == "PK01") {
					//string xyz = temp.substr(4, 42);
					string xyz = temp.substr(4);
					stringstream xxx(xyz);
					t_gtriple tmp;
					xxx >> tmp[0] >> tmp[1] >> tmp[2];

					_kin_xyz[site_upper][time_t] = tmp * 1000;
					continue;
				}
			}
			kinfile.close();
		}
		
	}

	bool t_glsqproc::_create_log(t_glog* log)
	{
		if (!log)
		{
			this->_glog = new t_glog();
			this->_glog->mask("podlsqALL.log");
			this->_glog->cache_size(99);
			this->_glog->tsys(t_gtime::GPS);
			this->_glog->time_stamp(true);
			this->_glog->verb(1);
			return true;
		}
		return false;
	}

	bool t_glsqproc::_delete_log(t_glog* log)
	{
		if (log)
		{
			delete log;
			log = nullptr;
			return true;
		}
		return false;
	}

	bool t_glsqproc::_add_sion_pars_value(const t_gtime& beg, const t_gtime& end, t_glsq* lsq, t_gion* gion)
	{

		for (const auto& it_rec : _rec_list)
		{
			for (const auto& it_sat : _sat_list)
			{
				int idx = _lsq->_x_solve.getParam(it_rec, t_gpar::SION, it_sat);
				if (idx < 0) continue;
				else
				{
					double value = gion ? gion->get_value(it_rec, it_sat, beg, end) : 0.0;
					_lsq->_x_solve[idx].value(value);
				}
			}
		}
		return false;
	}

	void t_glsqproc::_extract_ion_file(t_gallrecover& recover_data)
	{
		t_gion ion_data;
		vector<t_tuple_ion> ions;
		recover_data.get_iondata(ions);

		for (const auto& iter : ions)
		{
			t_gion_record tmp(std::get<0>(iter),
				std::get<1>(iter), std::get<2>(iter),
				std::get<3>(iter), std::get<4>(iter),
				std::get<5>(iter), std::get<6>(iter));

			ion_data.add_record(tmp);
		}

		// encoder data
		t_gfile gout;
		string ion_path = dynamic_cast<t_gsetout*>(_gset)->outputs("ion");

		if (ion_path.empty())
		{
			ion_path = "grt$(date).ion";
		}
		substitute(ion_path, "$(date)", _beg_time.str_yyyydoy(), false);

		gout.path(ion_path);
		gout.glog(_glog);

		t_ion coder(_gset);
		coder.add_data("ID1", &ion_data);
		gout.coder(&coder);

		gout.run_write();
		return;
	}

	bool t_glsqproc::_pod_add_RAW_CB_Dd_ambcon(const t_gtime& crt_epoch, t_glsq* lsq, t_gambcon* ambcon, const bool& isL12)
	{
		if (lsq == nullptr || ambcon == nullptr)
		{
			_glog->logError("t_gpodlsqALL", "_pod_add_RAW_CB_Dd_ambcon", "lsq or ambcon is empty");
			return false;
		}

		const auto& ddamb = ambcon->get_Dds(crt_epoch);
		if (ddamb.empty())
		{
			_glog->logError("t_gpodlsqALL", "_pod_add_RAW_CB_Dd_ambcon", "dd amb for " + crt_epoch.str_mjdsod() + " is empty");
			return false;
		}

		// 双差系数
		double op[5] = { 0.0, 1.0, -1.0, -1.0, 1.0 };
		for (const auto& amb : ddamb)
		{
			const auto& Dd_amb = amb;
			if (Dd_amb.ID == t_gamb::AMB_ID::ALL_L1_L2 && !isL12) continue;
			if (Dd_amb.ID == t_gamb::AMB_ID::ALL_L1_L3 && isL12) continue;

			// 将双差模糊度的测站和卫星提取出来
			string rec_1st = Dd_amb.rec_1st;
			string rec_2nd = Dd_amb.rec_2nd;
			string sat_1st = Dd_amb.sat_1st;
			string sat_2nd = Dd_amb.sat_2nd;

			const t_gtime& beg = Dd_amb.beg;
			const t_gtime& end = Dd_amb.end;

			// 确保四个都有
			string strType = "";
			t_gpar::t_type par_type1 = t_gpar::t_type::AMB_L1;
			t_gpar::t_type par_type2 = t_gpar::t_type::NO_DEF;
			FREQ_SEQ FQ1 = FREQ_1;
			FREQ_SEQ FQ2 = FREQ_1;
			if (Dd_amb.ID == t_gamb::AMB_ID::ALL_L1_L2)
			{
				par_type2 = t_gpar::t_type::AMB_L2;
				FQ2 = FREQ_2;
				strType = "ALL_L1_L2";
			}
			else if (Dd_amb.ID == t_gamb::AMB_ID::ALL_L1_L3)
			{
				par_type2 = t_gpar::t_type::AMB_L3;
				FQ2 = FREQ_3;
				strType = "ALL_L1_L3";
			}
			else
			{
				par_type2 = t_gpar::t_type::NO_DEF;
				continue;
			}
			//lijie
			t_gobsgnss gnss;
			gnss.sat(sat_1st);
			auto gsys = t_gsys::sat2gsys(amb.sat_1st);
			double length1 = CLIGHT / gnss.frequency(_band_index[gsys][FQ1]);
			double length2 = CLIGHT / gnss.frequency(_band_index[gsys][FQ2]);
			
			int amb1_1st = lsq->_x_solve.getAmbParam(rec_1st, sat_1st, par_type1, beg, end);
			int amb1_2nd = lsq->_x_solve.getAmbParam(rec_1st, sat_2nd, par_type1, beg, end);
			int amb1_3rd = lsq->_x_solve.getAmbParam(rec_2nd, sat_1st, par_type1, beg, end);
			int amb1_4th = lsq->_x_solve.getAmbParam(rec_2nd, sat_2nd, par_type1, beg, end);

			int amb2_1st = lsq->_x_solve.getAmbParam(rec_1st, sat_1st, par_type2, beg, end);
			int amb2_2nd = lsq->_x_solve.getAmbParam(rec_1st, sat_2nd, par_type2, beg, end);
			int amb2_3rd = lsq->_x_solve.getAmbParam(rec_2nd, sat_1st, par_type2, beg, end);
			int amb2_4th = lsq->_x_solve.getAmbParam(rec_2nd, sat_2nd, par_type2, beg, end);

			if (amb1_1st < 0 || amb1_2nd < 0 || amb1_3rd < 0 || amb1_4th < 0) continue;
			if (amb2_1st < 0 || amb2_2nd < 0 || amb2_3rd < 0 || amb2_4th < 0) continue;


			// 有可能确实取整后为0
			if (double_eq(Dd_amb.MW_value, 0.0)) continue;
			if (double_eq(Dd_amb.NL_value, 0.0)) continue;

			// 构建约束方程
			// constraint X(ix1)-X(ix2)-X(ix3)+X(ix4) = nc_fix - nc
			const auto& par_L1_1 = lsq->_x_solve.getAmbParam(amb1_1st);
			const auto& par_L1_2 = lsq->_x_solve.getAmbParam(amb1_2nd);
			const auto& par_L1_3 = lsq->_x_solve.getAmbParam(amb1_3rd);
			const auto& par_L1_4 = lsq->_x_solve.getAmbParam(amb1_4th);

			const auto& par_LX_1 = lsq->_x_solve.getAmbParam(amb2_1st);
			const auto& par_LX_2 = lsq->_x_solve.getAmbParam(amb2_2nd);
			const auto& par_LX_3 = lsq->_x_solve.getAmbParam(amb2_3rd);
			const auto& par_LX_4 = lsq->_x_solve.getAmbParam(amb2_4th);

			if (double_eq(par_L1_1.value(), 0.0)) continue;
			if (double_eq(par_L1_2.value(), 0.0)) continue;
			if (double_eq(par_L1_3.value(), 0.0)) continue;
			if (double_eq(par_L1_4.value(), 0.0)) continue;
			if (double_eq(par_LX_1.value(), 0.0)) continue;
			if (double_eq(par_LX_2.value(), 0.0)) continue;
			if (double_eq(par_LX_3.value(), 0.0)) continue;
			if (double_eq(par_LX_4.value(), 0.0)) continue;

			double nf1 = (round(Dd_amb.MW_value) + round(Dd_amb.NL_value)) * length1;
			double nf2 = (round(Dd_amb.NL_value)) * length2;

			double nc1 = par_L1_1.value() * op[1] +
				         par_L1_2.value() * op[2] +
				         par_L1_3.value() * op[3] +
				         par_L1_4.value() * op[4];

			double nc2 = par_LX_1.value() * op[1] +
				         par_LX_2.value() * op[2] +
				         par_LX_3.value() * op[3] +
				         par_LX_4.value() * op[4];


			double nl_L1 = nf1 - nc1;
			double nl_LX = nf2 - nc2;

			string str_nc = format(" %s %16.8lf %16.8lf", " NC = ", nc1, nc2);
			string str_nf = format(" %s %16.8lf %16.8lf", " NF = ", nf1, nf2);
			string str_nl = format(" %s %16.8lf %16.8lf", " NL = ", nl_L1, nl_LX);
			string str_beg = Dd_amb.beg.str_mjdsod(" beg ");
			string str_end = Dd_amb.end.str_mjdsod(" end ");
			string str_idx1 = " " + int2str(amb1_1st, 6) + " " + int2str(amb1_2nd, 6) + " " + int2str(amb1_3rd, 6) + " " + int2str(amb1_4th, 6) + " ";
			string str_idxX = " " + int2str(amb2_1st, 6) + " " + int2str(amb2_2nd, 6) + " " + int2str(amb2_3rd, 6) + " " + int2str(amb2_4th, 6) + " ";

			string str_line = strType + "_" + sat_1st + "_" + sat_2nd + "_" + rec_1st + "_" + rec_2nd;
			string strOut = str_line + str_idx1 + str_idxX + str_nc + str_nf + str_nl + str_beg + str_end;

			// 添加约束方程
			// 按照升序排好位置
			map<int, double> par_index_L1;
			par_index_L1[amb1_1st] = op[1];
			par_index_L1[amb1_2nd] = op[2];
			par_index_L1[amb1_3rd] = op[3];
			par_index_L1[amb1_4th] = op[4];

			map<int, double> par_index_LX;
			par_index_LX[amb2_1st] = op[1];
			par_index_LX[amb2_2nd] = op[2];
			par_index_LX[amb2_3rd] = op[3];
			par_index_LX[amb2_4th] = op[4];

			vector< pair<int, double> > B_L1;
			vector< pair<int, double> > B_LX;
			for (const auto& B_iter : par_index_L1)
			{
				B_L1.push_back(make_pair(B_iter.first + 1, B_iter.second));
			}
			for (const auto& B_iter : par_index_LX)
			{
				B_LX.push_back(make_pair(B_iter.first + 1, B_iter.second));
			}

			t_glsqEquationMatrix virtual_equ;
			t_gobscombtype combtype;

			if (Dd_amb.ID == t_gamb::AMB_ID::ALL_L1_L2 && isL12)
			{
				// 输出添加约束情况
				//cout << strOut << endl;

				t_amb_index key1 = make_tuple(amb1_1st, par_L1_1.site, par_L1_1.prn, par_L1_1.beg, par_L1_1.end);
				t_amb_index key2 = make_tuple(amb1_2nd, par_L1_2.site, par_L1_2.prn, par_L1_2.beg, par_L1_2.end);
				t_amb_index key3 = make_tuple(amb1_3rd, par_L1_3.site, par_L1_3.prn, par_L1_3.beg, par_L1_3.end);
				t_amb_index key4 = make_tuple(amb1_4th, par_L1_4.site, par_L1_4.prn, par_L1_4.beg, par_L1_4.end);
				auto key = make_tuple(key1, key2, key3, key4);

				if (_Dd_amb_L13_used.count(key) > 0)
				{
					if (!double_eq(_Dd_amb_L13_used[key], nl_L1))
					{
						cout << "remove for without L1 " + strOut << endl;
						continue;
					}

					virtual_equ.add_equ(B_LX, 1e7, nl_LX, rec_1st + "_" + rec_2nd, sat_1st + "_" + sat_2nd, combtype, false);
				}
				else
				{
					virtual_equ.add_equ(B_L1, 1e7, nl_L1, rec_1st + "_" + rec_2nd, sat_1st + "_" + sat_2nd, combtype, false);
					virtual_equ.add_equ(B_LX, 1e7, nl_LX, rec_1st + "_" + rec_2nd, sat_1st + "_" + sat_2nd, combtype,false);
				}

				_Dd_amb_L12_used[key] = nl_L1;

			}
			else if (Dd_amb.ID == t_gamb::AMB_ID::ALL_L1_L3 && !isL12)
			{
				t_amb_index key1 = make_tuple(amb1_1st, par_L1_1.site, par_L1_1.prn, par_L1_1.beg, par_L1_1.end);
				t_amb_index key2 = make_tuple(amb1_2nd, par_L1_2.site, par_L1_2.prn, par_L1_2.beg, par_L1_2.end);
				t_amb_index key3 = make_tuple(amb1_3rd, par_L1_3.site, par_L1_3.prn, par_L1_3.beg, par_L1_3.end);
				t_amb_index key4 = make_tuple(amb1_4th, par_L1_4.site, par_L1_4.prn, par_L1_4.beg, par_L1_4.end);
				auto key = make_tuple(key1, key2, key3, key4);

				if (_Dd_amb_L12_used.count(key) > 0)
				{
					if (!double_eq(_Dd_amb_L12_used[key], nl_L1))
					{
						cout << "remove for without L12 " + strOut << endl;
						continue;
					}
					virtual_equ.add_equ(B_LX, 1e7, nl_LX, rec_1st + "_" + rec_2nd, sat_1st + "_" + sat_2nd, combtype, false);
				}
				else
				{
					virtual_equ.add_equ(B_L1, 1e7, nl_L1, rec_1st + "_" + rec_2nd, sat_1st + "_" + sat_2nd, combtype, false);
					virtual_equ.add_equ(B_LX, 1e7, nl_LX, rec_1st + "_" + rec_2nd, sat_1st + "_" + sat_2nd, combtype, false);
				}

				_Dd_amb_L13_used[key] = nl_L1;

				// 输出添加约束情况
				//cout << strOut << endl;
			}
			else
			{
				throw logic_error("wrong type of dd amb while add ambcon for pod all in lsq");
			}
			//_glog->comment(t_glog::LOG_LV::LOG_INFO, "t_gpodlsqALL", "_pod_add_Dd_ambcon", "add ambcon in " + crt_epoch.str_ymdhms(strOut));
			lsq->add_equation(virtual_equ,_crt_time,false);
		}
		return true;
	}

	bool t_glsqproc::_pod_add_Dd_ambcon(t_glsq* lsq, t_gambcon* ambcon)
	{
		if (lsq == nullptr || ambcon == nullptr) return false;

		auto ddambs = _gdata_ambcon->get_Dds(_crt_time);
		if (ddambs.empty()) return false;
		// 双差系数
		double op[5] = { 0.0, 1.0, -1.0, -1.0, 1.0 };
		for (const auto& iter : ddambs)
		{
			const auto& Dd_amb = iter;

			// 如果两个没有同时固定
			if (!(Dd_amb.MW_fixed && Dd_amb.NL_fixed))
			{
				continue;
			}
			//ddamb=0
			if (double_eq(Dd_amb.MW_value, 0.0) || double_eq(Dd_amb.NL_value, 0.0)) continue;

			// 将双差模糊度的测站和卫星提取出来
			string rec_1st = Dd_amb.rec_1st;
			string rec_2nd = Dd_amb.rec_2nd;
			string sat_1st = Dd_amb.sat_1st;
			string sat_2nd = Dd_amb.sat_2nd;

			// add for single difference
			if (rec_2nd == "----")
			{
				_pod_add_one_Sd_ambcon(lsq, &Dd_amb);
				continue;
			}

			const t_gtime& beg = Dd_amb.beg;
			const t_gtime& end = Dd_amb.end;

			// 确保四个都有
			int par_1st = -1;
			int par_2nd	= -1;
			int par_3rd	= -1;
			int par_4th	= -1;

			if (Dd_amb.ID == t_gamb::AMB_ID::IF_L1_L2)
			{
				par_1st = lsq->_x_solve.getAmbParam(rec_1st, sat_1st, t_gpar::AMB_IF, beg, end);
				par_2nd = lsq->_x_solve.getAmbParam(rec_1st, sat_2nd, t_gpar::AMB_IF, beg, end);
				par_3rd = lsq->_x_solve.getAmbParam(rec_2nd, sat_1st, t_gpar::AMB_IF, beg, end);
				par_4th = lsq->_x_solve.getAmbParam(rec_2nd, sat_2nd, t_gpar::AMB_IF, beg, end);
			}
			else if (Dd_amb.ID == t_gamb::AMB_ID::IF_L1_L3)
			{
				par_1st = lsq->_x_solve.getAmbParam(rec_1st, sat_1st, t_gpar::AMB13_IF, beg, end);
				par_2nd = lsq->_x_solve.getAmbParam(rec_1st, sat_2nd, t_gpar::AMB13_IF, beg, end);
				par_3rd = lsq->_x_solve.getAmbParam(rec_2nd, sat_1st, t_gpar::AMB13_IF, beg, end);
				par_4th = lsq->_x_solve.getAmbParam(rec_2nd, sat_2nd, t_gpar::AMB13_IF, beg, end);
			}
			else
			{
				return false;
			}


			if (par_1st < 0 || par_2nd < 0 || par_3rd < 0 || par_4th < 0)
			{
				continue;
			}


			// 构建约束方程
			// constraint X(ix1)-X(ix2)-X(ix3)+X(ix4)=nc_fix-nc

			if (double_eq(lsq->_x_solve.getAmbParam(par_1st).value(), 0.0)) continue;
			if (double_eq(lsq->_x_solve.getAmbParam(par_2nd).value(), 0.0)) continue;
			if (double_eq(lsq->_x_solve.getAmbParam(par_3rd).value(), 0.0)) continue;
			if (double_eq(lsq->_x_solve.getAmbParam(par_4th).value(), 0.0)) continue;

			double nc = lsq->_x_solve.getAmbParam(par_1st).value() * op[1] +
				lsq->_x_solve.getAmbParam(par_2nd).value() * op[2] +
				lsq->_x_solve.getAmbParam(par_3rd).value() * op[3] +
				lsq->_x_solve.getAmbParam(par_4th).value() * op[4];

			// jdhuang : fix for muti-GNSS
			t_gobsgnss gnss;
			gnss.sat(sat_1st);
			auto gsys = t_gsys::sat2gsys(sat_1st);
			double wlength_1st = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_1]);
			double wlength_2nd = 0.0;
			if (Dd_amb.ID == t_gamb::AMB_ID::IF_L1_L2)
			{
				wlength_2nd = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_2]);
			}
			else if (Dd_amb.ID == t_gamb::AMB_ID::IF_L1_L3)
			{
				wlength_2nd = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_3]);
			}
			else
			{
				return false;
			}

			//double wlength_1st = CLIGHT / t_gamb::LC2FREQ(sat_1st, GOBS_LC::LC_L1);
			//double wlength_2nd = CLIGHT / t_gamb::LC2FREQ(sat_1st, GOBS_LC::LC_L2);
			double lambdaw = wlength_1st / (1 - wlength_1st / wlength_2nd);
			double lambdan = wlength_1st / (1 + wlength_1st / wlength_2nd);
			double nc_fix = round(Dd_amb.NL_value) + round(Dd_amb.MW_value) * lambdaw / wlength_2nd;
			nc_fix = nc_fix * lambdan;
			double nl = nc_fix - nc;
			//cout << fixed << nc_fix << "     " << nc << endl;
			//cout << sat_1st << "  " << sat_2nd << "  " << rec_1st << "  " << rec_2nd << "   " << beg.sod() << "  " << end.sod() << endl;
			/*string amb1 = format("%s,%16.8lf", "AMB1", lsq->_x_solve.getAmbParam(par_1st).value());
			string amb2 = format("%s,%16.8lf", "AMB2", lsq->_x_solve.getAmbParam(par_2nd).value());
			string amb3 = format("%s,%16.8lf", "AMB3", lsq->_x_solve.getAmbParam(par_3rd).value());
			string amb4 = format("%s,%16.8lf", "AMB4", lsq->_x_solve.getAmbParam(par_4th).value());
			string str_nc = format("%s,%16.8lf", "NC  ", nc);
			string str_nf = format("%s,%16.8lf", "NF  ", nc_fix);
			string str_nl = format("%16.8lf", nl);
			string strOut = amb1 + amb2 + amb3 + amb4 + str_nc + str_nf + str_nl;*/

			// 添加约束方程
			// 按照升序排好位置
			map<int, double> par_index;
			par_index[par_1st] = op[1];
			par_index[par_2nd] = op[2];
			par_index[par_3rd] = op[3];
			par_index[par_4th] = op[4];

			vector< pair<int, double> > B;
			for (const auto& B_iter : par_index)
			{
				B.push_back(make_pair(B_iter.first + 1, B_iter.second));
			}

			t_glsqEquationMatrix virtual_equ;
			t_gobscombtype type;

			// jdhuang : todo 
			// For IGSO and GEO, 1e6 or 1e5 maybe better
			virtual_equ.add_equ(B, 1e7, nl, rec_1st + "_" + rec_2nd, sat_1st + "_" + sat_2nd, type, false);

			_glog->comment(t_glog::LOG_LV::LOG_WARN, "t_gpodlsqIF", "_pod_add_Dd_ambcon", _crt_time.str_mjdsod("one ambcon add success"));

			lsq->add_equation(virtual_equ, _crt_time, false);

		}
		return true;
	}
	bool t_glsqproc::_pod_add_one_Sd_ambcon(t_glsq* lsq, const t_gamb_Dd* Dd_amb)
	{
		string rec = Dd_amb->rec_1st;
		string sat_1st = Dd_amb->sat_1st;
		string sat_2nd = Dd_amb->sat_2nd;
		double op[3] = { 0.0, 1.0, -1.0 };
		if (double_eq(Dd_amb->MW_value, 0.0) || double_eq(Dd_amb->NL_value, 0.0)) return false;
		const t_gtime& beg = Dd_amb->beg;
		const t_gtime& end = Dd_amb->end;
		int par_1st = lsq->_x_solve.getAmbParam(rec, sat_1st, t_gpar::AMB_IF, beg, end);
		int par_2nd = lsq->_x_solve.getAmbParam(rec, sat_2nd, t_gpar::AMB_IF, beg, end);
		if (par_1st < 0 || par_2nd < 0)
		{
			return false;
		}
		if (double_eq(lsq->_x_solve.getAmbParam(par_1st).value(), 0.0)) return false;
		if (double_eq(lsq->_x_solve.getAmbParam(par_2nd).value(), 0.0)) return false;
		double nc = lsq->_x_solve.getAmbParam(par_1st).value() * op[1] +
			lsq->_x_solve.getAmbParam(par_2nd).value() * op[2];

		t_gobsgnss gnss;
		gnss.sat(sat_1st);
		auto gsys = t_gsys::sat2gsys(sat_1st);
		double wlength_1st = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_1]);
		double wlength_2nd = CLIGHT / gnss.frequency(_band_index[gsys][FREQ_2]);
		//double wlength_1st = CLIGHT / t_gamb::LC2FREQ(sat_1st, GOBS_LC::LC_L1);
		//double wlength_2nd = CLIGHT / t_gamb::LC2FREQ(sat_1st, GOBS_LC::LC_L2);
		double lambdaw = wlength_1st / (1 - wlength_1st / wlength_2nd);
		double lambdan = wlength_1st / (1 + wlength_1st / wlength_2nd);
		double nc_fix = round(Dd_amb->NL_value) + round(Dd_amb->MW_value) * lambdaw / wlength_2nd;
		nc_fix = nc_fix * lambdan;
		nc_fix = nc_fix - Dd_amb->rnl*lambdan;
		double nl = nc_fix - nc;
		map<int, double> par_index;
		par_index[par_1st] = op[1];
		par_index[par_2nd] = op[2];
		vector< pair<int, double> > B;
		for (const auto& B_iter : par_index)
		{
			B.push_back(make_pair(B_iter.first + 1, B_iter.second));
		}
		t_glsqEquationMatrix virtual_equ;
		t_gobscombtype type;
		virtual_equ.add_equ(B, 1e7, nl, rec, sat_1st + "_" + sat_2nd, type, false);
		lsq->add_equation(virtual_equ);
		return true;
	}

	void t_glsqproc::_set_isb_par_setting(t_glsq* lsq, t_gpar& par, SYSBIASMODEL model, double value, double apriori)
	{
		par.value(value);
		par.apriori(apriori);

		if (model == AUTO_CON || model == ISB_CON)
		{
			par.setTime(_beg_time, _end_time);
			lsq->add_parameter(par);
		}
		else if (model == AUTO_RWK || model == ISB_RWK)
		{
			par.setTime(_beg_time, _beg_time);
			lsq->add_parameter(par);
			lsq->add_par_state_equ(par.parType, 1, 24, 3000.0);
		}
		else if (model == AUTO_WHIT)
		{
			par.setTime(_beg_time, _beg_time);
			lsq->add_parameter(par);
		}
	}


	bool t_glsqproc::_prefix_residuals(vector<t_gsatdata>& vec_obs)
	{
		set<string> bad_sats;
		string sat, type;
		GSYS gs;
		for (int iequ = 0; iequ < _crt_equ.num_equ(); iequ++)
		{
			sat = _crt_equ.get_satname(iequ);
			type = _crt_equ.get_obscombtype2str(iequ);
			gs = t_gsys::str2gsys(sat.substr(0, 1));

			if (!type.empty() && type[0] == 'L')
			{
				if (fabs(_crt_equ.l[iequ]) > _maxres_L[gs] * 600)          // 验前残差 与 验后残差 600倍
				{
					bad_sats.insert(sat);
					_glog->logInfo("t_glsqproc", "_prefix_residuals:", " outlier (v_orig: " + type + 
						" " + _crt_equ.get_sitename(iequ) + " " + sat + ":" + dbl2str(fabs(_crt_equ.l[iequ])) + ")");
				}
			}
			else if (!type.empty() && (type[0] == 'P' || type[0] == 'C'))  // 验前残差 与 验后残差 300倍
			{
				if (fabs(_crt_equ.l[iequ]) > _maxres_C[gs] * 300)
				{
					bad_sats.insert(sat);
					_glog->logInfo("t_glsqproc", "_prefix_residuals:", " outlier (v_orig: " + type + 
						" " + _crt_equ.get_sitename(iequ) + " " + sat + ":" + dbl2str(fabs(_crt_equ.l[iequ])) + ")");
				}
			}
			else
			{
				continue;
			}
		}

		bool eraseSat = false;
		auto data = vec_obs.begin();
		while (data != vec_obs.end())
		{
			string sat = data->sat();
			if (bad_sats.find(sat) != bad_sats.end())
			{
				data = vec_obs.erase(data);
				eraseSat = true;
			}
			else
			{
				data++;
			}
		}

		return eraseSat;
	}

}