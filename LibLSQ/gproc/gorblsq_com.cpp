//#include "gproc/gorblsq_com.h"
//#include "gutils/gleo.h"
//#include "gset/gsetout.h"
//#include "gio/gfile.h"
//#include "gset/gsetinp.h"
//#include "gcoders/icsnav.h"
//#include "gcoders/icsleo.h"
//
//
//namespace gnut {
//
//	t_gorblsq_com::t_gorblsq_com(t_gsetbase* set, t_glog* log):
//		t_glsqprocIF(set,log),_recclk_file("recclk"),_satclk_file("satclk")
//	{
//		try
//		{
//			//_sites = dynamic_cast<t_gsetgen*>(set)->recskd();
//			/*for (auto iter = _allsites.begin(); iter != _allsites.end(); iter++)
//			{
//				if (IS_SUPPORTED_LEO(*iter))
//				{
//					string satleo = strs2strl(*iter);
//					_satleo.insert(satleo);
//				}
//				
//			}*/
//			_site_ref = dynamic_cast<t_gsetproc*>(set)->ref_site();
//			_sat_ref = dynamic_cast<t_gsetproc*>(set)->ref_sat();
//			_sig_refclk = dynamic_cast<t_gsetproc*>(set)->sig_ref_clk();
//			_crd_est = dynamic_cast<t_gsetproc*>(_set)->crd_est();
//			_setOut();
//		}
//		catch (...)
//		{
//			if (_log)
//			{
//				_log->comment(1, "ERROR:  t_gorblsq_leo::t_gorblsq_leo() constructor  unknown mistake");
//			}
//			else
//			{
//				cout << "ERROR:  t_gorblsq_leo::t_gorblsq_leo() constructor  unknown mistake" << endl;
//			}
//		}
//	}
//
//	t_gorblsq_com::~t_gorblsq_com()
//	{
//		
//	}
//
//	bool t_gorblsq_com::GenerateProduct()
//	{
//		// update clk file
//		_lsq->extract_satclk(_satclk_file);
//		_lsq->extract_recclk(_recclk_file);
//
//		_update_icsfile();
//		_update_icsleo();
//
//		// Update ICS file
//		t_gio* gout = new t_gfile; t_gcoder* gics_coder = new t_icsnav(_set, "", 4096);
//		string ics_path = (dynamic_cast<t_gsetinp*>(_set)->inputs("ics"))[0] + "TEST";
//		gout->path(ics_path); gout->glog(_log);
//		gics_coder->clear(); gics_coder->glog(_log); gics_coder->path(ics_path);
//		gics_coder->add_data("ID1", _icsdata);
//		gics_coder->add_data("OBJ", _allobj);
//		gout->coder(gics_coder);
//		gout->run_write();
//		
//		
//		// Update ICSleo file
//		gout = new t_gfile; gics_coder = new t_icsleo(_set, "", 4096);
//		ics_path = (dynamic_cast<t_gsetinp*>(_set)->inputs("icsleo"))[0] + "TEST";
//		gout->path(ics_path); gout->glog(_log);
//		gics_coder->clear(); gics_coder->glog(_log); gics_coder->path(ics_path);
//		gics_coder->add_data("ID1", _icsleo);
//		gics_coder->add_data("OBJ", _allobj);
//		gout->coder(gics_coder);
//		gout->run_write();
//		if (gout)
//		{
//			delete gout; gout = NULL;
//		}
//		if (gics_coder)
//		{
//			delete gics_coder; gics_coder = NULL;
//		}
//		return true;
//	}
//
//	bool t_gorblsq_com::_InitLsq(t_gallproc* data)
//	{
//		_poledata = dynamic_cast<t_gpoleut1*>((*data)[t_gdata::ALLPOLEUT1]);
//		_icsdata = dynamic_cast<t_gallics_nav*>((*data)[t_gdata::ALLICS]);
//		_icsleo = dynamic_cast<t_gallics_leo*>((*data)[t_gdata::ALLICSLEO]);
//		if (!_poledata || !_icsdata||!_icsleo) 
//		{
//			if (_log) {
//				_log->comment(1, "poleut1 or icsnav or icsleo data not found!!!");
//			}
//			return false;
//		}
//
//		//_allsats = _icsdata->sat_list;
//		//_sys = _icsdata->sys_list;
//		set<string> newsats;
//		set_intersection(_sat_list.begin(), _sat_list.end()
//			, _icsdata->sat_list.begin(), _icsdata->sat_list.end(),
//			inserter(newsats, newsats.begin()));
//		_sat_list = newsats;
//		_satleo= _icsleo->getLeoList();
//		/*set<string> longsite = _icsleo->leo_list();
//		set<string> shortsite;
//		for (auto site = longsite.begin(); site != longsite.end(); site++)
//		{
//			shortsite.insert(strl2strs(*site));
//		}
//		_allsites = shortsite;*/
//
//		if (_sys_list.size() == 0 || _sat_list.size() == 0 || _rec_list.size() == 0) {
//			
//			if (_log) {
//				_log->comment(1, "ics file have no data!");
//			}
//			return false;
//		}
//
//		
//		// init CLK_SAT
//			//GSYS sys = t_gsys::str2gsys(*iter_sys);
//			//set<string> sats = GNSS_SATS()[sys];
//		int ipar_num = _lsq->_x_solve.parNumber();
//		for (auto iter = _sat_list.begin(); iter != _sat_list.end(); iter++) {
//			string satID = *iter;
//			double sigclk = 5000;
//			if (satID == _sat_ref && _sig_refclk == 0.0) continue;
//			else if (satID == _sat_ref && _sig_refclk != 0.0) sigclk = _sig_refclk;
//			t_gpar par_satclk("", t_gpar::CLK_SAT, ++ipar_num, satID);
//			par_satclk.setTime(_beg_time, _beg_time);
//			par_satclk.apriori(sigclk);
//			par_satclk.value(0.0);
//			_lsq->add_parameter(par_satclk);
//		}
//
//		_init_ORB_parameter();
//
//
//		// count par number
//		ipar_num = _lsq->_x_solve.parNumber();
//		// init xyz recclk ztd isb/ifb
//		// loop loop over sites from settings
//		for (set<string>::iterator it = _rec_list.begin(); it != _rec_list.end();)
//		{
//			string site = *it;
//			shared_ptr<t_gobj> _grec = _allobj->obj(site);
//			if (_grec->id_type() == t_gdata::REC_LEO)
//			{
//				it++;
//				continue;
//			}
//			// init xyz
//			if (_crd_est != FIX) {
//				
//
//				if (_grec == nullptr) {
//					it = _rec_list.erase(it);
//					continue;
//				}
//				// get the init value
//				t_gtriple rec_crd = _grec->crd(_beg_time);
//				t_gtriple rec_std = _grec->std(_beg_time);
//				// only need sinex file 
//				//if ( _grec->id_type() == t_gdata::REC && (rec_std.norm() > 0.1 || rec_crd.zero())) {
//				if (rec_std.norm() > 0.1 || rec_crd.zero()) {
//					if (_log)
//					{
//						_log->comment(1, site + "'s crd precision is low, erase!!!");
//					}
//					it = _rec_list.erase(it);
//					continue;
//				}
//				t_gpar par_x(site, t_gpar::CRD_X, ++ipar_num, "");
//				t_gpar par_y(site, t_gpar::CRD_Y, ++ipar_num, "");
//				t_gpar par_z(site, t_gpar::CRD_Z, ++ipar_num, "");
//				// set the parameter tim
//				par_x.value(rec_crd[0]);
//				par_y.value(rec_crd[1]);
//				par_z.value(rec_crd[2]);
//				// set the inital value
//				par_x.setTime(_beg_time, _end_time);
//				par_y.setTime(_beg_time, _end_time);
//				par_z.setTime(_beg_time, _end_time);
//
//				// set the apriori  [crd read from rinex]
//				par_x.apriori(rec_std[0]);
//				par_y.apriori(rec_std[1]);
//				par_z.apriori(rec_std[2]);
//
//				//  add parameters
//				_lsq->add_parameter(par_x);
//				_lsq->add_parameter(par_y);
//				_lsq->add_parameter(par_z);
//			}
//			else
//			{
//				if (_grec == nullptr) {
//					it = _rec_list.erase(it);
//					continue;
//				}
//				// get the init value
//				t_gtriple rec_crd = _grec->crd(_beg_time);
//				t_gtriple rec_std = _grec->std(_beg_time);
//				// only need sinex file 
//				if (rec_std.norm() > 0.1 || rec_crd.zero()) {
//					if (_log)
//					{
//						_log->comment(1, site + "'s crd precision is low, erase!!!");
//					}
//					it = _rec_list.erase(it);
//					continue;
//				}
//			}
//			if (site == _site_ref && _sig_refclk == 0.0) {
//				it++;
//				continue;
//			}
//			
//			
//			double sigclk = 9000;
//			if (site == _site_ref && _sig_refclk != 0.0)
//				sigclk = _sig_refclk;
//			//intit the Clock
//			t_gpar par_clk(site, t_gpar::CLK, ++ipar_num, "");
//			par_clk.value(0.0);
//			par_clk.setTime(_beg_time, _beg_time);
//			par_clk.apriori(sigclk);
//			_lsq->add_parameter(par_clk);
//
//			// init the ztd
//			t_gpar par_ztd(site, t_gpar::TRP, ++ipar_num, "");
//			t_gtriple xyz, ell;
//			if (_lsq->_x_solve.getCrdParam(site, xyz, _beg_time, _end_time) < 0) {
//				xyz = _grec->crd(_beg_time);
//			}
//			//_lsqEstimator->_x_solve.getCrdParam(site, xyz, _beg, _end);
//			xyz2ell(xyz, ell, false);
//			par_ztd.zhd = _gModel->tropoModel()->getZHD(ell, _beg_time);
//			par_ztd.apriori(0.201);
//			par_ztd.value(_gModel->tropoModel()->getZWD(ell, _beg_time));
//			par_ztd.setTime(_beg_time, _beg_time + 3600);
//			_lsq->add_parameter(par_ztd);
//			// PWC:120min  dt_hour=120/60=2  pd=0.015
//			_lsq->add_par_state_equ(par_ztd.parType, 1, 2, 0.015);
//			
//
//
//			// init the ISB/IFB 
//			// TODO  set in estimate equations
//			if (_sys_list.size() > 1)
//			{
//				set<string>::iterator iter_sys = _sys_list.begin();
//				// for scene that haven't GPS
//				if (_sys_list.find("GPS") == _sys_list.end())
//					iter_sys++;
//				// set ISB
//				for (; iter_sys != _sys_list.end(); iter_sys++)
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
//
//			// set GLO IFB
//			if (_sys_list.find("GLO") != _sys_list.end())
//			{
//				set<string> glo_sats = GNSS_SATS()[GLO];
//				for (auto sat_iter = glo_sats.begin(); sat_iter != glo_sats.end(); sat_iter++)
//				{
//					t_gpar par_GLOIFB(*it, t_gpar::GLO_IFPB, ++ipar_num, *sat_iter);
//					par_GLOIFB.value(0.0);
//					par_GLOIFB.setTime(_beg_time, _end_time);
//					par_GLOIFB.apriori(3000.0);
//					_lsq->add_parameter(par_GLOIFB);
//				}
//			}
//
//			if (site == _site_ref && _sig_refclk == 0.0) {
//				it++;
//				continue;
//			}
//			// cycle
//			it++;
//		}
//
//
//		//_init_parameter();
//		
//		_init_EOP_parameter();
//		_init_SAT_parameter();
//		
//		
//		
//
//
//		_lsq->reset_resfile("orb_resfile");
//
//		return true;
//	}
//
//
//	void t_gorblsq_com::_setOut()
//	{
//		try
//		{
//
//			string tmp = dynamic_cast<t_gsetout*>(_set)->outputs("recclk");
//			if (tmp.empty())
//			{
//				tmp = "rec_" + int2str(_beg_time.year()) + int2str(_beg_time.doy(false), 3);
//			}
//			else
//			{
//				substitute(tmp, "$(date)", int2str(_beg_time.year()) + int2str(_beg_time.doy(false), 3), false);
//			}
//			_recclk_file = tmp;
//			tmp = dynamic_cast<t_gsetout*>(_set)->outputs("satclk");
//
//			if (tmp.empty())
//			{
//				tmp = "sat_" + int2str(_beg_time.year()) + int2str(_beg_time.doy(false), 3);
//			}
//			else
//			{
//				substitute(tmp, "$(date)", int2str(_beg_time.year()) + int2str(_beg_time.doy(false), 3), false);
//			}
//			_satclk_file = tmp;
//
//			tmp = dynamic_cast<t_gsetout*>(_set)->outputs("leoclk");
//
//			/*if (tmp.empty())
//			{
//				tmp = "leo_" + int2str(_beg.year()) + int2str(_beg.doy(false), 3);
//			}
//			else
//			{
//				substitute(tmp, "$(date)", int2str(_beg.year()) + int2str(_beg.doy(false), 3), false);
//			}
//			_leoclk_file = tmp;*/
//		}
//		catch (...)
//		{
//			if (_log)
//			{
//				_log->comment(1, "ERROR:  t_gorblsq_com::_setOut  unknown mistake");
//			}
//			else
//			{
//				cout << "ERROR:  t_gorblsq_com::_setOut  unknown mistake" << endl;
//			}
//		}
//	}
//
//	bool t_gorblsq_com::_init_EOP_parameter()
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
//			_poledata->getEopData(_beg_time, data1);
//			_poledata->getEopData(_end_time, data2);
//		}
//		catch (exception e)
//		{
//			if (_log)
//			{
//				_log->comment(1, e.what());
//			}
//			return false;
//		}
//
//		par_xpole.setTime(_beg_time, _end_time);
//		par_xpole.value(data1["XPOLE"]);
//		par_xpole.apriori(0.3);
//		par_ypole.setTime(_beg_time, _end_time);
//		par_ypole.value(data1["YPOLE"]);
//		par_ypole.apriori(0.3);
//		par_dxpole.setTime(_beg_time, _end_time);
//		par_dxpole.value(data2["XPOLE"] - data1["XPOLE"]);
//		par_dxpole.apriori(0.03);
//		par_dypole.setTime(_beg_time, _end_time);
//		par_dypole.value(data2["YPOLE"] - data1["YPOLE"]);
//		par_dypole.apriori(0.03);
//		par_ut1.setTime(_beg_time, _end_time);
//		par_ut1.value(data1["UT1-TAI"]);
//		par_ut1.apriori(1E-4);
//		par_dut1.setTime(_beg_time, _end_time);
//		par_dut1.value(data2["UT1-TAI"] - data1["UT1-TAI"]);
//		par_dut1.apriori(2E-3);
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
//
//	}
//
//	bool t_gorblsq_com::_init_SAT_parameter()
//	{
//		int ipar_num = _lsq->_x_solve.parNumber();
//		// cycle sat list
//		for (auto sat_iter = _sat_list.begin(); sat_iter != _sat_list.end(); sat_iter++)
//		{
//			auto ics = _icsdata->getIcs(*sat_iter);
//			auto par_list = ics->get_par_list();
//			// cycle par list
//			for (auto par_iter = par_list.begin(); par_iter != par_list.end();par_iter++)
//			{
//				t_gpar::t_type par_type = t_gics_nav::_icspar2gpar(par_iter->first);
//				// init 
//				t_gpar par_temp("", par_type, ++ipar_num, *sat_iter);
//				par_temp.setTime(par_iter->second.beg,par_iter->second.end);
//				par_temp.value(par_iter->second.value);
//				// par_temp.apriori(par_iter->second.sig);
//				// according to the panda apriori value
//				if (par_type >= t_gpar::PXSAT && par_type <= t_gpar::PZSAT) {
//					par_temp.apriori(1E-2);
//				}
//				else if (par_type >= t_gpar::VXSAT && par_type <= t_gpar::VZSAT) {
//					par_temp.apriori(1E-4);
//				}
//				else if (par_type >= t_gpar::D0 && par_type <= t_gpar::Bs) {
//					par_temp.apriori(0.1);
//				}
//				//add
//				_lsq->add_parameter(par_temp);
//			}
//		}
//		return true;
//	}
//
//	/** @brief init LEO orbit parameters. */
//	void t_gorblsq_com::_init_ORB_parameter()
//	{
//		try
//		{
//			int ipar_num = _lsq->_x_solve.parNumber();
//			// cycle sat list
//			for (auto sat_iter = _satleo.begin(); sat_iter != _satleo.end(); sat_iter++)
//			{
//				auto ics = _icsleo->getIcsLeo(*sat_iter);
//				map<ICS_PAR, vector<t_icsleo_record>>* leopar;
//				leopar = ics->leoics();
//				int npar = 0;
//				// cycle par list
//				for (auto par_iter = leopar->begin(); par_iter != leopar->end(); par_iter++)
//				{
//					t_gpar::t_type par_type = t_gicsleo::icspar2gpar(par_iter->first);
//					vector<t_icsleo_record> onepar = par_iter->second;
//					for (int k = 0; k < onepar.size(); k++)
//					{
//						//int indexPWC = 0;
//						int npwc = ceil((_end_time - _beg_time) / (onepar[0].end - onepar[0].beg));
//						if (onepar[k].beg >= _beg_time && onepar[k].beg < _end_time || onepar[k].end > _beg_time && onepar[k].end <= _end_time)
//						{
//							npar++;
//							t_gpar par_temp(strl2strs(*sat_iter), par_type, ++ipar_num, "", false);
//							par_temp.setTime(onepar[k].beg, onepar[k].end);
//							par_temp.value(onepar[k].value);
//							//par_temp.apriori(onepar[k].sig);
//							par_temp.value(onepar[k].value);
//							/*if (double_eq(onepar[k].sig, 0.0))
//							{*/
//							if (par_type >= t_gpar::PXSAT && par_type <= t_gpar::PZSAT) {
//								par_temp.apriori(100.0);
//							}
//							else if (par_type >= t_gpar::VXSAT && par_type <= t_gpar::VZSAT) {
//								par_temp.apriori(1.0);
//							}
//							else {
//								par_temp.apriori(0.1);
//							}
//							/*switch (par_type)
//							{
//							case PXSAT:
//							case PYSAT:
//							case PZSAT:
//								par_temp.apriori(100.0); break;
//							case VXSAT:
//							case VYSAT:
//							case VZSAT:
//								par_temp.apriori(1.0); break;
//							default:
//								par_temp.apriori(0.1);
//								break;
//							}*/
//							/*}
//							else
//							{
//								par_temp.apriori(onepar[k].sig);
//							}*/
//							if (onepar.size() > 1)  // PWC parameters
//							{
//								par_temp.nPWC = npwc;
//								par_temp.lPWC = true;
//								//par_temp.indPWC = indexPWC;
//								//indexPWC++;
//							}
//							_lsq->add_parameter(par_temp);
//							_lsq->_x_solve.setOrbPar(par_temp);
//						}
//
//					}
//
//					
//				}
//				double sigclk = 9000;
//				string site = strl2strs(*sat_iter);
//				if (site == _site_ref && _sig_refclk != 0.0)
//					sigclk = _sig_refclk;
//				t_gpar par_clk(site, t_gpar::CLK, ++ipar_num, "");
//				par_clk.value(0.0);
//				par_clk.setTime(_beg_time, _beg_time);
//				par_clk.apriori(sigclk);
//				_lsq->add_parameter(par_clk);
//			}
//		}
//		catch (...)
//		{
//			if (_log)
//			{
//				_log->comment(1, "ERROR:  t_gorblsq_com::_init_ORB_parameter  unknown mistake");
//			}
//			else
//			{
//				cout << "ERROR:  t_gorblsq_com::_init_ORB_parameter  unknown mistake" << endl;
//			}
//		}
//	}
//
//	
//
//	void t_gorblsq_com::_update_icsfile()
//	{
//		// update
//		for (auto sat_iter = _sat_list.begin(); sat_iter != _sat_list.end(); sat_iter++)
//		{
//
//			auto ics = _icsdata->getIcs(*sat_iter);
//			auto& par_list = ics->get_par_list();
//
//			// cycle par list
//			for (auto par_iter = par_list.begin(); par_iter != par_list.end(); par_iter++)
//			{
//				t_gpar::t_type par_type = t_gics_nav::_icspar2gpar(par_iter->first);
//
//				int par_idx = _lsq->_x_solve.getParam("", par_type, *sat_iter);
//
//				ColumnVector dx = _lsq->dx();
//				ColumnVector stdx = _lsq->stdx();
//
//
//				// update value && sig
//				cout << setprecision(15) << scientific
//					<< setw(20) << par_iter->second.value << setw(20) << dx(par_idx + 1) << endl;
//				par_iter->second.value += dx(par_idx + 1);
//				par_iter->second.sig = stdx(par_idx + 1);
//
//			}
//
//		}
//	}
//	/** @brief upsate LEO ics information. */
//	void t_gorblsq_com::_update_icsleo()
//	{
//		try
//		{
//			// update
//			for (auto sat_iter = _satleo.begin(); sat_iter != _satleo.end(); sat_iter++)
//			{
//
//				auto ics = _icsleo->getIcsLeo(*sat_iter);
//				map<ICS_PAR, vector<t_icsleo_record>>* par_list = ics->leoics();
//
//				// cycle par list
//				for (auto par_iter = par_list->begin(); par_iter != par_list->end(); par_iter++)
//				{
//					t_gpar::t_type par_type = t_gicsleo::icspar2gpar(par_iter->first);
//					vector<t_icsleo_record> &pars = par_iter->second;
//					for (int k = 0; k < pars.size(); k++)
//					{
//						t_gtime begtime = pars[k].beg;
//						t_gtime endtime = pars[k].end;
//						int par_idx = _lsq->_x_solve.getParam(strl2strs(*sat_iter), par_type, "", begtime, endtime);
//						if (par_idx >= 0)
//						{
//							ColumnVector dx = _lsq->dx();
//							ColumnVector stdx = _lsq->stdx();
//
//							// update value && sig
//							pars[k].value += dx(par_idx + 1);
//							pars[k].sig = stdx(par_idx + 1);
//
//						}
//
//					}
//
//				}
//
//			}
//		}
//		catch (...)
//		{
//			if (_log)
//			{
//				_log->comment(1, "ERROR:  t_gorblsq_leo::_update_icsleo  unknown mistake");
//			}
//			else
//			{
//				cout << "ERROR:  t_gorblsq_leo::_update_icsleo  unknown mistake" << endl;
//			}
//		}
//	}
//}