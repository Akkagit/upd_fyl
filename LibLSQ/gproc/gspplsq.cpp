/*

*/

#include "gproc/gspplsq.h"
#include "gset/gsetpar.h"


namespace great
{

	t_gspplsq::t_gspplsq(string mark, t_gsetbase* set, t_gallobj* obj, bool dynamic)
		:t_gspp(mark, set),
		_minsat(static_cast<size_t>(SPP_MINSAT)),
		_lsq(nullptr),
		_result(0),
		_dynamics(dynamic)
	{
		gtrace("t_spplsq::tgspplsq");

		// set obj
		setOBJ(obj);

		// set site && inono free
		dynamic_cast<t_gsppmodel*>(_gModel)->setrec(obj->obj(_site));
		if (_observ != RAW_MIX && _observ!= RAW_SINGLE)
		{
			dynamic_cast<t_gsppmodel*>(_gModel)->reset_observ(IONO_FREE);
			_observ = IONO_FREE;
		}
		_sig_ion = dynamic_cast<t_gsetpar*>(set)->sigSion(_site);

		// get set from setfile
		_get_own_settings();

		// set output
		_setOut();

	}

	t_gspplsq::~t_gspplsq()
	{
		gtrace("t_gspplsq::~g_gspplsq");

		if (_lsq)
			delete _lsq;

		if (_result) { if (_result->is_open()) { _result->close(); }; delete _result; }
	}

	int t_gspplsq::processBatch(const t_gtime& beg, const t_gtime& end)
	{
		_obs_num.clear();
		gtrace("t_gspplsq::processBatch");
		int suc_nepoch = 0;
		if (beg > end) { cerr << _site << " - processing not started [beg > end]\n"; return -1; }

		// set the _dop
		_dop.set_data(_gnav, _gobs, _site);

		// define the time begin/end/intv/now
		t_gtime	begT(beg);
		t_gtime endT(end);
		t_gtime now(beg);

		// time loop
		while (now <= end)
		{
			int invalid_equation = ProcessOneEpoch(now, NULL);
			if (invalid_equation <= 0) {
				now = now + _sampling;
				continue;
			}
			// update time
			now = now + _sampling;
			suc_nepoch++;
		}

		return suc_nepoch;
	}

	int t_gspplsq::ProcessOneEpoch(const t_gtime& now, vector<t_gsatdata>* vec_obs)
	{
		// get the obs !! This is the Raw obs !!!
		_data.erase(_data.begin(), _data.end());

		if (vec_obs) _data = *vec_obs;
		else         _data = _gobs->obs(_site, now);

		// no obs
		if (_data.size() == 0) {
			if (_log) { _log->comment(1, "gspplsq", _site + now.str_ymdhms(" no observation found at epoch: ")); }
			return -1;
		}

		// apply dcb by lvhb 20200803
		if (_gallbias)
		{
			for (auto& itdata : _data)
			{
				itdata.apply_bias(*_gallbias);
			}
		}

		int invalid_equation = _processEpoch(now);
		if (invalid_equation <= 0) {
			if (_log) {
				_log->comment(1, "gspplsq", now.str_ymdhms("havo no equations at the epoch "));
			}
			return -1;
		}

		return 1;
	}

	double t_gspplsq::getRecClk(const t_gtime& time)
	{
		auto param = _lsq->_x_solve;
		int i = param.getParam(_site, t_gpar::CLK, "");
		return i<0 ? 0.0: param[i].value();
	}

	t_gtriple t_gspplsq::getvel(const t_gtime& time)
	{
		auto param = _lsq->_x_solve;
		int i = param.getParam(_site, t_gpar::VEL_X, "");
		int j = param.getParam(_site, t_gpar::VEL_Y, "");
		int k = param.getParam(_site, t_gpar::VEL_Z, "");
		t_gtriple vel(param[i].value(), param[j].value(), param[k].value());
		return vel;
	}

	t_gtriple t_gspplsq::getCrd(const t_gtime& time)
	{
		if (_map_crd.find(time) != _map_crd.end()) {
			return _map_crd[time];
		}
		auto left_crd = _map_crd.lower_bound(time);
		auto right_crd = left_crd;
		if (_map_crd.size() > 0)
		{
			left_crd--;
		}
		if (left_crd != _map_crd.end() && right_crd != _map_crd.end()) {
			return (left_crd->second + right_crd->second) / 2;
		}
		else if (left_crd != _map_crd.end()) {
			return left_crd->second;
		}
		else if (right_crd != _map_crd.end()) {
			return right_crd->second;
		}
		else {
			return t_gtriple(0.0, 0.0, 0.0);
		}
	}

	t_gtriple t_gspplsq::getAvergeCrd()
	{
		t_gtriple crd_aver(0.0, 0.0, 0.0);
		for (auto iter_crd = _map_crd.begin(); iter_crd != _map_crd.end(); iter_crd++)
		{
			crd_aver += iter_crd->second;
		}
		crd_aver /= _map_crd.size();
		return crd_aver;
	}

	set<string> t_gspplsq::outcomsat(map<string, t_gsatdata>& satdata)
	{
		set<string> com_sat;

		for (vector<t_gsatdata>::iterator it = _data.begin(); it != _data.end(); it++) {
			string comsat = it->sat();
			com_sat.insert(comsat);
			satdata[comsat] = *it;
		}

		return com_sat;
	}

	map<t_gtime, int> t_gspplsq::getobsnum()
	{
		return _obs_num;
	}

	int t_gspplsq::_processEpoch(const t_gtime& runEpoch)
	{


		// construct lsq estimator
		t_gtriple crd_init(0.0, 0.0, 0.0);
		double  clk_init = 0.0;
		//crd_init = _grec->crd(runEpoch);
		int invalid_equation = 0;
		int count_iter = 0;
		t_gtriple vel_init(0.0, 0.0, 0.0);
		double  vclk_init = 0.0;
		// iterator
		while (count_iter <= 10)
		{
			// construct lsq
			if (_lsq)
				delete _lsq;
			_lsq = new t_glsq(nullptr);
			t_gpar parx(_site, t_gpar::CRD_X, 0, "");
			t_gpar pary(_site, t_gpar::CRD_Y, 1, "");
			t_gpar parz(_site, t_gpar::CRD_Z, 2, "");
			t_gpar parclk(_site, t_gpar::CLK, 3, "");
			parx.value(crd_init[0]);
			pary.value(crd_init[1]);
			parz.value(crd_init[2]);
			parclk.value(clk_init);

			_lsq->add_parameter(parx);
			_lsq->add_parameter(pary);
			_lsq->add_parameter(parz);
			_lsq->add_parameter(parclk);
			if (_dynamics) {
				t_gpar parv_x(_site, t_gpar::VEL_X, 4, "");
				t_gpar parv_y(_site, t_gpar::VEL_Y, 5, "");
				t_gpar parv_z(_site, t_gpar::VEL_Z, 6, "");
				t_gpar parv_clk(_site, t_gpar::CLK_RAT, 7, "");
				parv_x.value(vel_init[0]);
				parv_y.value(vel_init[1]);
				parv_z.value(vel_init[2]);
				parv_clk.value(vclk_init);
				_lsq->add_parameter(parv_x);
				_lsq->add_parameter(parv_y);
				_lsq->add_parameter(parv_z);
				_lsq->add_parameter(parv_clk);
			}
			if (_sys.size() >= 2)
			{
				string ref_sys;
				if      (_sys.count("GPS") > 0) ref_sys = "GPS";
				else if (_sys.count("GAL") > 0) ref_sys = "GAL";
				else if (_sys.count("BDS") > 0) ref_sys = "BDS";
				else if (_sys.count("QZS") > 0) ref_sys = "QZS";
				else                            ref_sys = "GLO";

				t_gpar par_isbifb;
				for (auto it_sys : _sys)
				{	
					GSYS   crt_sys = t_gsys::str2gsys(it_sys);
					switch (crt_sys)
					{
					case gnut::GPS:
						break;
					case gnut::GAL:
						par_isbifb = t_gpar(_site, t_gpar::GAL_ISB, _lsq->_x_solve.parNumber() + 1, "");
						par_isbifb.value(0.0);
						_lsq->add_parameter(par_isbifb);
						break;
					case gnut::BDS:
						par_isbifb = t_gpar(_site, t_gpar::BDS_ISB, _lsq->_x_solve.parNumber() + 1, "");
						par_isbifb.value(0.0);
						_lsq->add_parameter(par_isbifb);
						break;
					case gnut::QZS:
						par_isbifb = t_gpar(_site, t_gpar::QZS_ISB, _lsq->_x_solve.parNumber() + 1, "");
						par_isbifb.value(0.0);
						_lsq->add_parameter(par_isbifb);
						break;
					case gnut::GLO:
						for (auto itdata : _data)
						{
							string sat = itdata.sat();
							if (sat.substr(0, 1) != "R") continue;
							par_isbifb = t_gpar(_site, t_gpar::GLO_IFB, _lsq->_x_solve.parNumber() + 1, sat);
							par_isbifb.value(0.0);
							_lsq->add_parameter(par_isbifb);
						}
						break;
					default:
						break;
					}
				}
			}

			// Process this Epoch and get the equations
			Matrix B; DiagonalMatrix P; ColumnVector l;
			// Check whether the receiver is valid or not
			if (_grec == nullptr) {
				if (_log) _log->comment(2, "gspplsq", "No receiver settings available!!!");
				return -1;
			}
			// Get the current epoch
			_epoch = runEpoch;

			//  the fist estimator of _data (erase some unused data)
			if (_prepareData() == 0) {
				if (_log) _log->comment(2, "gspplsq", "This epoch have no useful data");
				return -1;
			}

			// for obs set log && set IONO
			for (auto iter = _data.begin(); iter != _data.end(); iter++) {
				iter->glog(_log);
				if (_observ == RAW_MIX)
				{
					t_gpar parSION(_site, t_gpar::SION, _lsq->_x_solve.parNumber() + 1, iter->sat());
					parSION.value(0.0);
					_lsq->add_parameter(parSION);
				}
			}

			// Get the parameter 
			t_gallpar& param = _lsq->_x_solve;
			t_gtriple groundEll, groundXYZ;
			param.getCrdParam(_site, groundXYZ);
			xyz2ell(groundXYZ.crd_array(), groundEll.crd_array(), false);

			// define a number of measurements 
			unsigned int nObs = _data.size();
			unsigned int mult = 1;
			//if (_observ == RAW_DOUBLE) { mult = 2; nObs *= 2; }
			//if (_observ == RAW_ALL) { mult = 2; nObs *= 5; } // reservation for 5 freq - not used raws will be removed
			if (_observ == RAW_MIX) { mult = 2; nObs *= 5; } // reservation for 5 freq - not used raws will be removed
			if (_phase) { mult *= 2; nObs *= 2; }  // code + phase

			unsigned int nPar = param.parNumber();

			// resize the Matrix B l P
			if (_dynamics) nObs *= 5;
			B.ReSize(nObs, nPar); B = 0.0;
			l.ReSize(nObs); l = 0.0;
			P.ReSize(nObs);	P = 0.0;

			// count for useful obs
			unsigned int iobs = 1;

			// loop over the obs data generate Matrix 
			// TODO Maybe the obstype should be optional?
			for (auto it = _data.begin(); it != _data.end();) {
				if (_addObsP(*it, iobs, groundEll, B, l, P) <= 0) {
					it = _data.erase(it);
					continue;
				}
				invalid_equation++;
				if (_addObsL(*it, iobs, groundEll, B, l, P) <= 0) {
					it = _data.erase(it);
					continue;
				}
				if (_dynamics) {
					if (_addObsD(*it, iobs, groundEll, B, l, P) <= 0) {
						it = _data.erase(it);
						continue;
					}
				}
				it++;
			}

			// clean the zero row or column
			iobs--;
			//invalid_equation = iobs;
			B = B.Rows(1, iobs);
			P = P.SymSubMatrix(1, iobs);
			l = l.Rows(1, iobs);

			if (invalid_equation < 4) {
				if (_log) {
					_log->comment(1, "gspplsq ", _site + runEpoch.str_ymdhms("no enough equations:"));
				}
				break;
			}

			// add the new equations
			_lsq->add_equation(B, P, l, runEpoch);
			

			if (_log) {
				_log->comment(1, "gspplsq", _site + runEpoch.str_ymdhms("processing epoch:"));
			}

			try
			{
				_lsq->solve_NEQ();
			}
			catch (NPDException)
			{
				if (_log) {
					_log->comment(1, _site + _epoch.str_ymdhms() + "_lsq->solve_NEQ failed");
				}
				invalid_equation = -1;
				break;
			}
			t_gtriple dx((_lsq->dx()).Rows(1, 3));
			t_gtriple crd_temp(0.0, 0.0, 0.0);
			_lsq->_x_solve.getCrdParam(_site, crd_temp);
			crd_temp += dx;
			if ((_lsq->dx()).Rows(1, 4).NormFrobenius() < 1E-4) {
				if (_dynamics && count_iter != 10)
				{
					if ((_lsq->dx()).Rows(5, 8).NormFrobenius() < 1E-6) {
						break;
					}
				}
				else {
					break;
				}

			}
			crd_init = crd_temp;
			clk_init += _lsq->dx()(4);

			if (_dynamics) {
				t_gtriple dv((_lsq->dx()).Rows(5, 7));
				int i = param.getParam(_site, t_gpar::VEL_X, "");
				int j = param.getParam(_site, t_gpar::VEL_Y, "");
				int k = param.getParam(_site, t_gpar::VEL_Z, "");
				vel_init[0] = param[i].value();
				vel_init[1] = param[j].value();
				vel_init[2] = param[k].value();
				vel_init += dv;
				vclk_init += _lsq->dx()(8);
			}
			count_iter++;
			_obs_num[runEpoch] = P.nrows();
		}

		if (invalid_equation >= 4 && count_iter <= 10) {
			_map_crd[runEpoch] = crd_init;
			_prepareData();
			return invalid_equation;
		}
		else {
			if (_log) {
				_log->comment(1, runEpoch.str_ymdhms() + ": solve fail!!");
			}
			return -1;
		}

		return -1;
	}

	int	t_gspplsq::_get_own_settings()
	{
		_sys = dynamic_cast<t_gsetgen*>(_set)->sys();
		_beg = dynamic_cast<t_gsetgen*>(_set)->beg();
		_end = dynamic_cast<t_gsetgen*>(_set)->end();
		_pos_kin = dynamic_cast<t_gsetproc*>(_set)->pos_kin();
		_crd_est = dynamic_cast<t_gsetproc*>(_set)->crd_est();

		return 1;
	}

	int t_gspplsq::_prepareData()
	{
		gtrace("t_gsppflt::_prepareData");
		vector<t_gsatdata>::iterator it = _data.begin();

		// erase the unused data of gnss system unused
		it = _data.begin();
		while (it != _data.end()) {
			switch (_gnss) {
			case GNS: break;
			case GPS: if (it->gsys() != GPS) { it = _data.erase(it); continue; } break;
			case GLO: if (it->gsys() != GLO) { it = _data.erase(it); continue; } break;
			case GAL: if (it->gsys() != GAL) { it = _data.erase(it); continue; } break;
			case BDS: if (it->gsys() != BDS) { it = _data.erase(it); continue; } break;
			case QZS: if (it->gsys() != QZS) { it = _data.erase(it); continue; } break;
			case IRN: if (it->gsys() != IRN) { it = _data.erase(it); continue; } break;
			case SBS: if (it->gsys() != SBS) { it = _data.erase(it); continue; } break;
			}
			it++;
		}

		// ============= add for spp parameters initialation glfeng
		Matrix BB;

		BB.ReSize(_data.size(), 4);
		BB = 0.0;
		int iobs = 0;
		// =============
		// filter the obs according to the obs && the sat

		// if obs is used ?
		// TODO: the obs type is LC PC;it should be optional
		it = _data.begin();
		while (it != _data.end()) {
			GSYS gs = it->gsys();
			GOBSBAND b1, b2;

			b1 = t_gsys::band_priority(gs, FREQ_1);
			b2 = t_gsys::band_priority(gs, FREQ_2);

			double obsP = it->P3(b1, b2);

			if (double_eq(obsP, 0.0) && (_observ == RAW_MIX || _observ == RAW_SINGLE))
			{
				if(!double_eq(it->obs_L(t_gband(b1, GOBSATTR::ATTR)), 0.0))
					obsP = it->obs_C(t_gband(b1, GOBSATTR::ATTR));
			}

			if (double_eq(obsP, 0.0)) {
				if (_log) {
					_log->comment(2, "gspplsq", "have no PIF erase the sat:" + it->sat());
				}
				it = _data.erase(it);
				continue;
			}

			if (_satPos(_epoch, *it) < 0) {
				if (_log) {
					ostringstream str;
					str << "prepareData: erasing data since _satPos failed, "
						<< "epo: " << _epoch.str_hms() << ", " <<
						"prn: " << it->sat();
					_log->comment(2, "gsppflt", str.str());
				}
				it = _data.erase(it);
				continue;
			}

			iobs++;
			BB(iobs, 1) = it->satcrd().crd(0);
			BB(iobs, 2) = it->satcrd().crd(1);
			BB(iobs, 3) = it->satcrd().crd(2);
			//BB(iobs, 4) = it->P3(b1, b2) + it->clk();
			BB(iobs, 4) = obsP + it->clk();

			it++;
		}
		//lvhb
		//cout << "bb" << endl;
		//cout << fixed << setprecision(6) << setw(15) << BB << endl;

		int i = _lsq->_x_solve.getParam(_site, t_gpar::CRD_X, "");
		int j = _lsq->_x_solve.getParam(_site, t_gpar::CRD_Y, "");
		int k = _lsq->_x_solve.getParam(_site, t_gpar::CRD_Z, "");
		int m = _lsq->_x_solve.getParam(_site, t_gpar::CLK, "");

		if (double_eq(_lsq->_x_solve[i].value(), 0.0) && double_eq(_lsq->_x_solve[j].value(), 0.0) &&
            double_eq(_lsq->_x_solve[k].value(), 0.0))
		{
			if (_data.size() < _minsat) {
				if (_log)
					_log->comment(2, "gspplsq", _epoch.str_ymdhms(_site + " epoch ")
						+ " skipped (Bancroft not calculated: " + int2str(_data.size()) + " < _minsat)");
				return -1;
			}

			BB = BB.Rows(1, iobs);    // delete zero rows

			if (BB.Nrows() < static_cast<int>(_minsat)) 
			{
				if (_log && _log->verb() >= 2)
					_log->comment(2, "gsppflt", _epoch.str_ymdhms(_site + " epoch ")
						+ " skipped (Bancroft not calculated: BB.Nrows < _minsat)");
				return -1;
			}

			ColumnVector vBanc;
			gbancroft(BB, vBanc);

			_lsq->_x_solve[i].value(vBanc(1));
			_lsq->_x_solve[j].value(vBanc(2));
			_lsq->_x_solve[k].value(vBanc(3));
			_lsq->_x_solve[m].value(vBanc(4));
		}

		// if sattelite is used?
		//	relative have been thought in clk
		it = _data.begin();
		while (it != _data.end()) {
			if (it->addprd(_gnav) < 0) {
				if (_log) {
					_log->comment(2, "gspplsq", "have no eph erase the sat:" + it->sat());
				}
				it = _data.erase(it);
				continue;
			}
			it++;
		}

		// elvation
		// compute sat elevation and rho
		// TODO  real time should change 
		it = _data.begin();
		while (it != _data.end()) {
			t_gtriple xyz_r, xyz_s, neu_s;
			xyz_s = it->satcrd();
			// ARP Correction
			//xyz_r = _grec->crd_arp(_epoch);
			_lsq->_x_solve.getCrdParam(_site, xyz_r);

			t_gtriple xyz_rho = xyz_s - xyz_r;
			t_gtriple ell_r;
			xyz2ell(xyz_r, ell_r, false);

			xyz2neu(ell_r, xyz_rho, neu_s);

			// Earth rotation correction
			t_gtriple xRec;
			double rho0 = sqrt(pow(xyz_r[0] - xyz_s[0], 2) + pow(xyz_r[1] - xyz_s[1], 2) + pow(xyz_r[2] - xyz_s[2], 2));
			double dPhi = OMEGA * rho0 / CLIGHT;
			xRec[0] = xyz_r[0] * cos(dPhi) - xyz_r[1] * sin(dPhi);
			xRec[1] = xyz_r[1] * cos(dPhi) + xyz_r[0] * sin(dPhi);
			xRec[2] = xyz_r[2];

			// Apply tides 
			_apply_tides(_epoch, xRec);

			double tmp = (it->satcrd().crd_cvect() - xRec.crd_cvect()).norm_Frobenius();

			it->addrho(tmp);

			double NE2 = neu_s[0] * neu_s[0] + neu_s[1] * neu_s[1];
			double ele = acos(sqrt(NE2) / it->rho());
			if (neu_s[2] < 0.0) {
				ele *= -1.0;
			}

			if (sqrt(NE2) / it->rho() > 1.0)
				it->addele(0.0);
			else it->addele(ele);

			double azi = atan2(neu_s[1], neu_s[0]);
			if (azi < 0) azi += 2 * G_PI;
			it->addazi(azi);

			// check elevation cut-off 
			//if (it->ele_deg() < _minElev) {
			//	if (_log) {
			//		_log->comment(2, "gspplsq", "elevation too small erase the sat:" + it->sat());
			//	}
			//	it = _data.erase(it);
			//	continue;
			//}

			it++;
		}

		// TODO add the res

		return _data.size();

	}


	int t_gspplsq::_addObsP(t_gsatdata& satdata, unsigned int& iobs, t_gtriple& ell, Matrix& B, ColumnVector& OMC, DiagonalMatrix& P)
	{
		gtrace("t_gspplsq::_addObsP");

		// TODO obstype should be optional !!
		// TODO Maybeneed our own obstype
		// just use the PC here
		//  This part  for getting the obs according to the OBS
		t_gsys gsys(satdata.gsys());
		GSYS gs = gsys.gsys();

		GOBSBAND b1, b2;
		b1 = t_gsys::band_priority(gs, FREQ_1);
		b2 = t_gsys::band_priority(gs, FREQ_2);

		t_gobs gobs1(TYPE_C, b1, ATTR);
		t_gobs gobs2(TYPE_C, b2, ATTR);

		int addObs = 0;
		double Pobs = satdata.P3(gobs1, gobs2);
		double P1 = satdata.obs_C(gobs1);
		double P2 = satdata.obs_C(gobs2);

		if (_observ == IONO_FREE && double_eq(Pobs, 0.0)) return -1;
		if (_observ == RAW_SINGLE && double_eq(P1, 0.0)) return -1;
		if (_observ == RAW_MIX && double_eq(P1, 0.0) && double_eq(P2, 0.0)) return -1;

		bool com = true;
		if (_gnav) com = _gnav->com();

		// this part for getting the resuids l (gpsmod)
		// TODO:: applyDCB
		//_applyDCB(satdata, Pobs, _observ);
		if (_observ == RAW_SINGLE || _observ == RAW_MIX)
		{
			double modelPobs = _gModel->cmpObs(_epoch, _lsq->_x_solve, satdata, gobs1, com);
			if (modelPobs < 0) return -1;
			OMC(iobs) = P1 - modelPobs;

			if (_observ == RAW_MIX && !double_eq(P2, 0.0))
			{
				addObs = 1;
				modelPobs = _gModel->cmpObs(_epoch, _lsq->_x_solve, satdata, gobs2, com);
				if (modelPobs < 0) addObs = 0;
				else               OMC(iobs + addObs) = P2 - modelPobs;
			}
		}
		else
		{
			double modelPobs = _gModel->cmpObs(_epoch, _lsq->_x_solve, satdata, gobs1, com);
			if (modelPobs < 0) return -1;
			OMC(iobs) = Pobs - modelPobs;
		}


		// This part for getting the  Weight P
		// TODO Maybe edit the _weight
		if (_observ == RAW_SINGLE || _observ == RAW_MIX)
		{
			double weightObs = _weightObs(satdata, _observ, &gobs1);
			if (weightObs < 0) return -1;
			P(iobs, iobs) = weightObs;
			P(iobs + addObs, iobs + addObs) = weightObs;
		}
		else
		{
			double weightObs = _weightObs(satdata, _observ, &gobs1, &gobs2);
			if (weightObs < 0) return -1;
			P(iobs, iobs) = weightObs;
		}


		// This part for getting the  Design Matrix B
		t_gallpar& param = _lsq->_x_solve;
		for (unsigned int ipar = 1; ipar <= param.parNumber(); ipar++) {
			B(iobs, ipar) = param[ipar - 1].partial(satdata, _epoch, ell, gobs1);
			if(_observ == RAW_MIX && addObs == 1) 
				B(iobs + 1, ipar) = param[ipar - 1].partial(satdata, _epoch, ell, gobs2);

		}

		// Only have PC
		iobs++;
		iobs += addObs; // greater than 0 only RAW_MIX.

		return 1;
	}

	// SPP have no use of phase
	int	 t_gspplsq::_addObsL(t_gsatdata& satdata, unsigned int& iobs, t_gtriple& ell, Matrix& B, ColumnVector& OMC, DiagonalMatrix& P)
	{
		gtrace("t_gspplsq::_addObsL");

		return 1;
	}

	int t_gspplsq::_addObsD(t_gsatdata& satdata, unsigned int& iobs, t_gtriple& XYZ, Matrix& A, ColumnVector& l, DiagonalMatrix& P)
	{
		gtrace("t_gpvtflt::_addObsD");

		t_gsys gsys(satdata.gsys());
		GSYS gs = gsys.gsys();

		GOBSBAND b1, b2, b3, b4, b5;
		b1 = b2 = b3 = b4 = b5 = BAND;

		if (_auto_band) {    // automatic dual band selection -> for Anubis purpose
			set<GOBSBAND> bands = satdata.band_avail();
			auto itBAND = bands.begin();
			if (bands.size() < 2) return -1;
			b1 = *itBAND;
			itBAND++;
			b2 = *itBAND;
		}
		else {           // select fix defined band according the table
			b1 = t_gsys::band_priority(gs, FREQ_1);
			b2 = t_gsys::band_priority(gs, FREQ_2);
			b3 = t_gsys::band_priority(gs, FREQ_3);
			b4 = t_gsys::band_priority(gs, FREQ_4);
			b5 = t_gsys::band_priority(gs, FREQ_5);
		}

		//t_gobs gobs1, gobs2, gobs3, gobs4, gobs5;
		t_gobs gobs1(TYPE_D, b1, ATTR);
		t_gobs gobs2(TYPE_D, b2, ATTR);
		t_gobs gobs3(TYPE_D, b3, ATTR);
		t_gobs gobs4(TYPE_D, b4, ATTR);
		t_gobs gobs5(TYPE_D, b5, ATTR);

		double D1, D2, D3, D4, D5;
		D1 = D2 = D3 = D4 = D5 = 0.0;
		vector<int> index_freq(5, 0);
		D1 = satdata.obs_D(gobs1);
		D2 = satdata.obs_D(gobs2);
		if (b3 != BAND) D3 = satdata.obs_D(gobs3);
		if (b4 != BAND) D4 = satdata.obs_D(gobs4);
		if (b5 != BAND) D5 = satdata.obs_D(gobs5);

		double weightObs = _weightObs(satdata, _observ, &gobs1, &gobs2);
		//double weight_coef = _weightObs(satdata, gobs1);
		int obsnum = 0;
		if (!double_eq(D1, 0)) { obsnum++; index_freq[0] = 1; }
		if (!double_eq(D2, 0)) { obsnum++; index_freq[1] = 1; }
		if (!double_eq(D3, 0)) { obsnum++; index_freq[2] = 1; }
		if (!double_eq(D4, 0)) { obsnum++; index_freq[3] = 1; }
		if (!double_eq(D5, 0)) { obsnum++; index_freq[4] = 1; }

		if (obsnum == 0)return -1;

		//Create reduced measurements (prefit residuals)
		double modObsD = 0.0;
		int n = 0;
		auto& param = _lsq->_x_solve;
		for (int i = 0; i < index_freq.size(); i++)
		{
			double Di = 0.0;
			if (!index_freq[i]) continue;
			if (i == 0) { modObsD = _gModel->cmpObsD(_epoch, param, satdata, gobs1); Di = D1; }
			if (i == 1) { modObsD = _gModel->cmpObsD(_epoch, param, satdata, gobs2); Di = D2; }
			if (i == 2) { modObsD = _gModel->cmpObsD(_epoch, param, satdata, gobs3); Di = D3; }
			if (i == 3) { modObsD = _gModel->cmpObsD(_epoch, param, satdata, gobs4); Di = D4; }
			if (i == 4) { modObsD = _gModel->cmpObsD(_epoch, param, satdata, gobs5); Di = D5; }
			/*if (fabs(-Di - modObsD) > 0.1)
			{
				obsnum--;
				continue;
			}*/
			if (double_eq(modObsD, 0.0) || double_eq(Di, 0.0)) continue;
			l(iobs + n) = -Di - modObsD; n++;
		}

		// Create weight matrix
		for (int i = 0; i < obsnum; i++)
			P(iobs + i, iobs + i) = weightObs;

		// Create first design matrix
		double A_conf = 0.0;
		int i = param.getParam(_site, t_gpar::VEL_X, "");
		int j = param.getParam(_site, t_gpar::VEL_Y, "");
		int k = param.getParam(_site, t_gpar::VEL_Z, "");
		t_gtriple groundVel(param[i].value(), param[j].value(), param[k].value());

		n = 0;
		for (int i = 0; i < index_freq.size(); i++)
		{
			if (!index_freq[i]) continue;
			for (unsigned int ipar = 1; ipar <= param.parNumber(); ipar++)
			{
				if (i == 0) A_conf = param[ipar - 1].partial_doppler(satdata, XYZ, groundVel);
				if (i == 1) A_conf = param[ipar - 1].partial_doppler(satdata, XYZ, groundVel);
				if (i == 2) A_conf = param[ipar - 1].partial_doppler(satdata, XYZ, groundVel);
				if (i == 3) A_conf = param[ipar - 1].partial_doppler(satdata, XYZ, groundVel);
				if (i == 4) A_conf = param[ipar - 1].partial_doppler(satdata, XYZ, groundVel);
				if (double_eq(A_conf, 0.0))continue;
				if (param[ipar - 1].parType == t_gpar::CRD_X || param[ipar - 1].parType == t_gpar::CRD_Y
					|| param[ipar - 1].parType == t_gpar::CRD_Z)A_conf = 0.0;
				A(iobs + n, ipar) = A_conf;
			}
			n++;
		}

		iobs = iobs + obsnum;

		return 1;
	}

	// Satelite position
    // ----------
	int t_gspplsq::_satPos(t_gtime& epo, t_gsatdata& gsatdata)
	{
		gtrace("t_gsppflt::_satPos");

		string satname = gsatdata.sat();

		int i;

		int base_size = dynamic_cast<t_gsetgen*>(_set)->list_base().size();
		auto nppmodel = dynamic_cast<t_gsetproc*>(_set)->npp_model();//lvhb added
		// modified by zhShen
		if (dynamic_cast<t_gsetproc*>(_set)->realtime() && !base_size && nppmodel != NPP_MODEL::URTK)
			i = gsatdata.addprd_realtime(_gnav); //add sat crd and clk
		else
			i = gsatdata.addprd(_gnav); //add sat crd and clk

		if (i < 0) return i;

		return 1;
	}


	int t_gspplsq::_apply_tides(t_gtime& _epoch, t_gtriple& xRec)
	{
		// spp don't apply tides
		return 1;
	}

	double t_gspplsq::_weightObs(t_gsatdata& satdata, OBSCOMBIN obstype, t_gobs* go1, t_gobs* go2)
	{
		// TODO according to obs type
		// This weight is computing according to the Panda

		GSYS gs = satdata.gsys();
		GOBSTYPE type = go1->type();

		double factor = 1.0;
		double sigRange = 0.0;
		// get the sys sigRange
		if (type == TYPE_C) {
			switch (gs)
			{
			case GPS: sigRange = _sigCodeGPS; break;
			case GLO: sigRange = _sigCodeGLO; break;
			case GAL: sigRange = _sigCodeGAL; break;
			case BDS: sigRange = _sigCodeBDS; break;
			case QZS: sigRange = _sigCodeQZS; break;
			default: sigRange = 0.0; return -1;
			}
		}
		else if (type == TYPE_L) {
			switch (gs)
			{
			case GPS: sigRange = _sigPhaseGPS; break;
			case GLO: sigRange = _sigPhaseGLO; break;
			case GAL: sigRange = _sigPhaseGAL; break;
			case BDS: sigRange = _sigPhaseBDS; break;
			case QZS: sigRange = _sigPhaseQZS; break;
			default: sigRange = 0.0; return -1;
			}
		}
		else if (type == TYPE_D) {
			switch (gs) {
			case GPS: sigRange = _sigDopplerGPS; break;
			case GLO: sigRange = _sigDopplerGLO; break;
			case GAL: sigRange = _sigDopplerGAL; break;
			case BDS: sigRange = _sigDopplerBDS; break;
			case QZS: sigRange = _sigDopplerQZS; break;
			default: sigRange = 0.0;
			}
		}
		else {
			return -1;
		}

		// get the _weight factor
		switch (_weight) {
		case DEF_OBSWEIGHT: cerr << "gspplsq: WeightObs (default) should not happened!\n"; break;
		case EQUAL:	factor *= 1; break;
		case SINEL: factor *= 1.0 / 2.0 / sin(satdata.ele()); break;
		case SINEL2: factor *= 1.0 / 2.0 / pow(sin(satdata.ele()), 2); break;
		case SINEL4: factor *= 1.0 / 2.0 / pow(sin(satdata.ele()), 4); break;
		case PARTELE:
			if (satdata.ele_deg() <= 30.0) {
				factor *= (1.0 / 2.0 / sin(satdata.ele()));
			}
			else factor *= 1.0;
			break;

		default:
			cerr << "gspplsq: we can't deal with this WeightObs method!\n";
			return -1;
		}

		// get the combine obs factor
		// TODO finish other weight method 
		double coef1 = 0.0, coef2 = 0.0;
		switch (obstype)
		{
		case gnut::IONO_FREE:
			satdata.coef_ionofree(go1->band(), coef1, go2->band(), coef2);
			sigRange *= sqrt(coef1 * coef1 + coef2 * coef2);
			break;
		default:
			//cerr << "gspplsq: No other weighting method have been implmented!";
			return 1.0 / pow(factor * sigRange, 2);
		}


		return 1.0 / pow(factor * sigRange, 2);
	}




	void t_gspplsq::_setOut()
	{
		// set res file
		t_gspp::_setOut();

		// set output file
		string tmp;
		tmp = dynamic_cast<t_gsetout*>(_set)->outputs("flt");
		if (!tmp.empty()) {
			substitute(tmp, "$(rec)", _site, false);
			_result = new t_giof;
			_result->tsys(t_gtime::GPS);
			_result->mask(tmp);
			_result->append(dynamic_cast<t_gsetout*>(_set)->append());
		}

	}

	void t_gspplsq::_printout(const t_glsq& lsqestimator)
	{
		ostringstream os;

		t_gtime epoch = lsqestimator.epo();
		t_gallpar par = lsqestimator._x_solve;
		SymmetricMatrix Qx = lsqestimator.Qx();
		ColumnVector dx = lsqestimator.dx();
		string str_dsec = dbl2str(epoch.dsec());

		int idx_x, idx_y, idx_z;
		idx_x = par.getParam(_site, t_gpar::CRD_X, "");
		idx_y = par.getParam(_site, t_gpar::CRD_Y, "");
		idx_z = par.getParam(_site, t_gpar::CRD_Z, "");

		t_gtriple xyz, ell, xyz_standard;
		xyz[0] = par[idx_x].value() + dx(idx_x + 1);
		xyz[1] = par[idx_y].value() + dx(idx_y + 1);
		xyz[2] = par[idx_z].value() + dx(idx_z + 1);

		// just for test by ZHJ
		xyz_standard[0] = -3.60766537057412e+06;
		xyz_standard[1] = 4.14786796276078e+06;
		xyz_standard[2] = 3.22371704084420e+06;

		xyz2ell(xyz_standard, ell, false);

		t_gtriple dxyz, neu;
		dxyz = xyz - xyz_standard;

		xyz2neu(ell, dxyz, neu);

		os << "LSQ:= " << fixed << setprecision(4)
			<< epoch.str_ymdhms() << str_dsec.substr(2) << setprecision(4);
		//os << setw(20) << neu[0]
		//	<< setw(20) << neu[1]
		//	<< setw(20) << neu[2];
		//os << endl;
		for (unsigned int ipar = 0; ipar < par.parNumber(); ipar++)
		{
			os << setw(15) << par[ipar].value() + dx(ipar + 1);
		}
		os << endl;

		// Print flt results
		if (_result) {
			_result->write(os.str().c_str(), os.str().size());
			_result->flush();
		}
		else {
			_log->comment(1, "t_gspplsq::_printout", "have no output file!");
		}
	}

}