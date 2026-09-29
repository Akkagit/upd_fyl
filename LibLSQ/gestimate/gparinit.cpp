//#include "gparinit.h"
//#include "gset/gsetpar.h"
//#include "gdata/gpoleut1.h"
//
//namespace gnut
//{
//
//	t_glsqparInit::t_glsqparInit(t_glsq* lsq,t_gsetbase* set):
//		_lsq(lsq),_set(set)
//	{
//		_beg_time = dynamic_cast<t_gsetgen*>(_set)->beg();
//		_end_time = dynamic_cast<t_gsetgen*>(_set)->end();
//	}
//
//	t_glsqparInit::~t_glsqparInit()
//	{
//
//	}
//
//	bool t_glsqparInit::initRecCrdPars(set<string>& sites, t_gallobj* allobj, double std_threshold)
//	{
//
//		int ipar_num = _lsq->_x_solve.parNumber();
//		for (auto it = sites.begin(); it != sites.end(); )
//		{
//			string site = *it;
//			shared_ptr<t_gobj> _grec = allobj->obj(site);
//			if (!_grec) 
//			{
//				it = sites.erase(it);
//				continue;
//			}
//			CONSTRPAR crd_est = dynamic_cast<t_gsetproc*>(_set)->crd_est();
//
//			// get the init value
//			t_gtriple rec_crd = _grec->crd(_beg_time);
//			t_gtriple rec_std = _grec->std(_beg_time);
//
//			// judge threshold and crd intial
//			if ((std_threshold > 0 && rec_std.norm() > std_threshold) || rec_crd.zero()) 
//			{
//				cerr << site + "'s crd precision is low, erase!!!" << endl;
//				it = sites.erase(it);
//				continue;
//			}
//			// init xyz 
//			if (crd_est != FIX) 
//			{
//				t_gpar par_x(site, t_gpar::CRD_X, ++ipar_num, "");
//				t_gpar par_y(site, t_gpar::CRD_Y, ++ipar_num, "");
//				t_gpar par_z(site, t_gpar::CRD_Z, ++ipar_num, "");
//
//				// set the parameter tim
//				par_x.value(rec_crd[0]);
//				par_y.value(rec_crd[1]);
//				par_z.value(rec_crd[2]);
//				// set the inital value
//				par_x.setTime(_beg_time, _end_time);
//				par_y.setTime(_beg_time, _end_time);
//				par_z.setTime(_beg_time, _end_time);
//				// set the apriori  [crd read from rinex]
//				par_x.apriori(rec_std[0]);
//				par_y.apriori(rec_std[1]);
//				par_z.apriori(rec_std[2]);
//				//  add parameters
//				_lsq->add_parameter(par_x);
//				_lsq->add_parameter(par_y);
//				_lsq->add_parameter(par_z);
//			}
//			it++;
//		}
//
//		return true;
//	}
//
//	bool t_glsqparInit::initRecTroPars(set<string>& sites)
//	{
//		int ipar_num = _lsq->_x_solve.parNumber();
//		for (auto it = sites.begin(); it != sites.end(); it++)
//		{
//			string site = *it;
//
//			// init the ztd
//			t_gpar par_ztd(site, t_gpar::TRP, ++ipar_num, "");
//
//			//double sig_ztd = 0.201;
//			double sig_ztd = dynamic_cast<t_gsetpar*>(_set)->getSigZtd(site);
//			par_ztd.apriori(sig_ztd);
//		
//			// init in precisemodel
//			double dt = 0.0; // minute
//			ZTDMODEL ztd_model = dynamic_cast<t_gsetproc*>(_set)->ztd_model(&dt);
//			par_ztd.value(0.0);
//			par_ztd.setTime(_beg_time, _beg_time + dt * 60.0);
//
//			_lsq->add_parameter(par_ztd);
//
//			// PWC:120min  dt_hour=120/60=2  pd=0.015
//			double sig_TropPd = dynamic_cast<t_gsetpar*>(_set)->getSigTropPd(site);
//			if (ztd_model == PWC) 
//			{
//				_lsq->add_par_state_equ(par_ztd.parType, 1, dt / 60.0, sig_TropPd);
//			}
//		}
//		return true;
//	}
//
//	bool t_glsqparInit::initRecClkPars(set<string>& sites)
//	{
//		int ipar_num = _lsq->_x_solve.parNumber();
//		for (auto it = sites.begin(); it != sites.end(); it++) {
//
//			string site = *it;
//			//double sigclk = 9000;
//			double sigclk = dynamic_cast<t_gsetpar*>(_set)->getRecSigClk(site);
//			string site_ref = dynamic_cast<t_gsetproc*>(_set)->ref_site();
//			double sig_refclk = dynamic_cast<t_gsetproc*>(_set)->sig_ref_clk();
//			if (site == site_ref) {
//				if (sig_refclk == 0.0) {
//					it++;
//					continue;
//				}
//				else {
//					sigclk = sig_refclk;
//				}
//			}
//			//intit the Clock
//			t_gpar par_clk(site, t_gpar::CLK, ++ipar_num, "");
//			par_clk.value(0.0);
//			par_clk.setTime(_beg_time, _beg_time);
//			par_clk.apriori(sigclk);
//			_lsq->add_parameter(par_clk);
//		}
//		return true;
//	}
//
//	bool t_glsqparInit::initSysBiaPars(set<string>& sites, set<string>& sys)
//	{
//		int ipar_num = _lsq->_x_solve.parNumber();
//		// init the ISB/IFB 
//		// TODO  set in estimate equations
//		if (sys.size() > 1)
//		{
//			for (auto it = sites.begin(); it != sites.end(); it++)
//			{
//				string site = *it;
//				set<string>::iterator iter_sys = sys.begin();
//				// for scene that haven't GPS
//				if (sys.find("GPS") == sys.end())
//					iter_sys++;
//				// set ISB
//				for (; iter_sys != sys.end(); iter_sys++)
//				{
//					t_gpar* par_ISB;
//					GSYS sys = t_gsys::str2gsys(*iter_sys);
//					switch (sys)
//					{
//					case gnut::GPS:
//						break;
//					case gnut::GAL:
//						par_ISB = new t_gpar(site, t_gpar::GAL_ISB, ++ipar_num, "");
//						par_ISB->value(0.0);
//						par_ISB->setTime(_beg_time, _end_time);
//						par_ISB->apriori(3000.0);
//						_lsq->add_parameter(*par_ISB);
//						break;
//					case gnut::GLO:
//						par_ISB = new t_gpar(site, t_gpar::GLO_ISB, ++ipar_num, "");
//						par_ISB->value(0.0);
//						par_ISB->setTime(_beg_time, _end_time);
//						par_ISB->apriori(3000.0);
//						_lsq->add_parameter(*par_ISB);
//						break;
//					case gnut::BDS:
//						par_ISB = new t_gpar(site, t_gpar::BDS_ISB, ++ipar_num, "");
//						par_ISB->value(0.0);
//						par_ISB->setTime(_beg_time, _end_time);
//						par_ISB->apriori(3000.0);
//						_lsq->add_parameter(*par_ISB);
//						break;
//					case gnut::QZS:
//						par_ISB = new t_gpar(site, t_gpar::QZS_ISB, ++ipar_num, "");
//						par_ISB->value(0.0);
//						par_ISB->setTime(_beg_time, _end_time);
//						par_ISB->apriori(3000.0);
//						_lsq->add_parameter(*par_ISB);
//						break;
//					default:
//						break;
//					}
//					delete par_ISB;
//				}
//			}
//		}
//
//		// set GLO IFB
//		if (sys.find("GLO") != sys.end())
//		{
//			for (auto it = sites.begin(); it != sites.end(); it++)
//			{
//				string site = *it;
//				set<string> glo_sats = GNSS_SATS()[GLO];
//				for (auto sat_iter = glo_sats.begin(); sat_iter != glo_sats.end(); sat_iter++)
//				{
//					t_gpar par_GLOIFB(site, t_gpar::GLO_IFPB, ++ipar_num, *sat_iter);
//					par_GLOIFB.value(0.0);
//					par_GLOIFB.setTime(_beg_time, _end_time);
//					par_GLOIFB.apriori(3000.0);
//					_lsq->add_parameter(par_GLOIFB);
//				}
//			}
//		}
//		return true;
//	}
//
//	bool t_glsqparInit::initSatClkPars(set<string>& sats)
//	{
//		int ipar_num = _lsq->_x_solve.parNumber();
//		for (auto iter = sats.begin(); iter != sats.end(); iter++) {
//			string satID = *iter;
//			//double sigclk = 5000;
//			double sigclk = dynamic_cast<t_gsetpar*>(_set)->getSatSigClk(satID);
//			string sat_ref = dynamic_cast<t_gsetproc*>(_set)->ref_sat();
//			double sig_refclk = dynamic_cast<t_gsetproc*>(_set)->sig_ref_clk();
//			if (satID == sat_ref && sig_refclk == 0.0) continue;
//			else if (satID == sat_ref && sig_refclk != 0.0) sigclk = sig_refclk;
//			t_gpar par_satclk("", t_gpar::CLK_SAT, ++ipar_num, satID);
//			par_satclk.setTime(_beg_time, _beg_time);
//			par_satclk.apriori(sigclk);
//			par_satclk.value(0.0);
//			_lsq->add_parameter(par_satclk);
//		}
//		return true;
//	}
//
//	bool t_glsqparInit::initGnsIcsPars(set<string>& sats, t_gdata * icsdata)
//	{
//		int ipar_num = _lsq->_x_solve.parNumber();
//		// cycle sat list
//		for (auto sat_iter = sats.begin(); sat_iter != sats.end(); sat_iter++)
//		{
//			auto ics = dynamic_cast<t_gallics_nav*>(icsdata)->getIcs(*sat_iter);
//			auto par_list = ics->get_par_list();
//			// cycle par list
//			for (auto par_iter = par_list.begin(); par_iter != par_list.end(); par_iter++)
//			{
//				t_gpar::t_type par_type = t_gics_nav::_icspar2gpar(par_iter->first);
//				// init 
//				t_gpar par_temp("", par_type, ++ipar_num, *sat_iter);
//				par_temp.setTime(par_iter->second.beg, par_iter->second.end);
//				par_temp.value(par_iter->second.value);
//				// par_temp.apriori(par_iter->second.sig);
//				// according to the panda apriori value
//				auto sig_pos = dynamic_cast<t_gsetpar*>(_set)->getSatSigPos(*sat_iter);
//				auto sig_vel = dynamic_cast<t_gsetpar*>(_set)->getSatSigVel(*sat_iter);
//				auto sig_force = dynamic_cast<t_gsetpar*>(_set)->getSatSigEcom(*sat_iter);
//				string id = "sig" + par_temp.str_type().substr(0, 2);
//				// pos
//				if (par_type >= t_gpar::PXSAT && par_type <= t_gpar::PZSAT) {
//					//par_temp.apriori(1E-2);
//					par_temp.apriori(sig_pos[id]/1000.0); // change unit to km
//				}
//				// vel
//				else if (par_type >= t_gpar::VXSAT && par_type <= t_gpar::VZSAT) {
//					//par_temp.apriori(1E-4);
//					par_temp.apriori(sig_vel[id] / 1000.0); // change unit to km
//				}
//				// force
//				else if (par_type >= t_gpar::D0 && par_type <= t_gpar::Bs) {
//					//par_temp.apriori(0.1);
//					par_temp.apriori(sig_force[id]);
//				}
//				//add
//				_lsq->add_parameter(par_temp);
//			}
//		}
//		return true;
//	}
//
//	bool t_glsqparInit::initCrsEopPars(t_gdata * poledata)
//	{
//		int ipar_num = _lsq->_x_solve.parNumber();
//		// init EOP parmeter
//		// init the Eop dx dy vdx vdy
//		t_gpar par_xpole("", t_gpar::XPOLE, ++ipar_num, "");
//		t_gpar par_ypole("", t_gpar::YPOLE, ++ipar_num, "");
//		t_gpar par_dxpole("", t_gpar::DXPOLE, ++ipar_num, "");
//		t_gpar par_dypole("", t_gpar::DYPOLE, ++ipar_num, "");
//		t_gpar par_ut1("", t_gpar::UT1, ++ipar_num, "");
//		t_gpar par_dut1("", t_gpar::DUT1, ++ipar_num, "");
//
//
//		// TODO:init value time sigma		
//		map<string, double> data1, data2;
//
//		try
//		{
//			dynamic_cast<t_gpoleut1*>(poledata)->getEopData(_beg_time, data1);
//			dynamic_cast<t_gpoleut1*>(poledata)->getEopData(_end_time, data2);
//		}
//		catch (exception e)
//		{
//			cerr << e.what();
//			return false;
//		}
//
//		double sig_xpole = dynamic_cast<t_gsetpar*>(_set)->getSigXpole();
//		double sig_ypole = dynamic_cast<t_gsetpar*>(_set)->getSigYpole();
//		double sig_dxpole = dynamic_cast<t_gsetpar*>(_set)->getSigDxpole();
//		double sig_dypole = dynamic_cast<t_gsetpar*>(_set)->getSigDypole();
//		double sig_ut1 = dynamic_cast<t_gsetpar*>(_set)->getSigUt1();
//		double sig_dut1 = dynamic_cast<t_gsetpar*>(_set)->getSigDut1();
//
//		par_xpole.setTime(_beg_time, _end_time);
//		par_xpole.value(data1["XPOLE"]);
//		par_xpole.apriori(sig_xpole);
//		par_ypole.setTime(_beg_time, _end_time);
//		par_ypole.value(data1["YPOLE"]);
//		par_ypole.apriori(sig_ypole);
//		par_dxpole.setTime(_beg_time, _end_time);
//		par_dxpole.value(data2["XPOLE"] - data1["XPOLE"]);
//		par_dxpole.apriori(sig_dxpole);
//		par_dypole.setTime(_beg_time, _end_time);
//		par_dypole.value(data2["YPOLE"] - data1["YPOLE"]);
//		par_dypole.apriori(sig_dypole);
//		par_ut1.setTime(_beg_time, _end_time);
//		par_ut1.value(data1["UT1-TAI"]);
//		par_ut1.apriori(sig_ut1);
//		par_dut1.setTime(_beg_time, _end_time);
//		par_dut1.value(data2["UT1-TAI"] - data1["UT1-TAI"]);
//		par_dut1.apriori(sig_dut1);
//
//		// add the Estimator
//		_lsq->add_parameter(par_xpole);
//		_lsq->add_parameter(par_ypole);
//		_lsq->add_parameter(par_dxpole);
//		_lsq->add_parameter(par_dypole);
//		_lsq->add_parameter(par_ut1);
//		_lsq->add_parameter(par_dut1);
//
//		return true;
//	}
//}