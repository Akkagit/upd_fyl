#include "gpvtflt_old.h"
#include "gutils/gtimesync.h"
#include "gproc/gpreproc.h"
#include "gproc/gexeturboedit.h"
#include "gmodels/gprecisemodel.h"
#include "gcoders/bnccorr.h"
#include "gcoders/bncobs.h"

gsins::t_gpvtflt_old::t_gpvtflt_old(string mark, string mark_base, t_gsetbase* gset)
	: t_gpppflt(mark, gset),
	t_gspp(mark, gset),
	_fix_mode(t_gsetamb::FIX_MODE::NO),
	_isFirstFix(true),
	_amb_state(false),
	_site_base(mark_base),
	_gModel_base(nullptr),
	_pos_constrain(false)
{
	gtrace("t_gpvtflt::t_gpvtflt");
	_vel = t_gtriple(0, 0, 0);  _Qx_vel.ReSize(4); _Qx_vel = 0.0;
	_dclkStoModel = new t_whitenoise(dynamic_cast<t_gsetflt*>(_set)->noise_dclk());
	_velStoModel = new t_whitenoise(dynamic_cast<t_gsetflt*>(_set)->noise_vel());
	_fix_mode = dynamic_cast<t_gsetamb*>(_set)->fix_mode();
	_max_res_norm = dynamic_cast<t_gsetproc*>(_set)->max_res_norm();
	_sd_sat = dynamic_cast<t_gsetproc*>(_set)->sd_sat();
	_isBase = false;
	if (!_site_base.empty())_isBase = true;
	if (_isBase) _rtkinit();

	// amb fix
	if (_ambfix)
	{
		delete _ambfix;
		_ambfix = nullptr;
	}
	if (_fix_mode != t_gsetamb::FIX_MODE::NO)
	{
		_ambfix = new t_gambiguity(_site, _set);
	}

}

gsins::t_gpvtflt_old::~t_gpvtflt_old()
{
	gtrace("t_gpvtflt::~t_gpvtflt");
	if (_ambfix)
	{
		delete _ambfix;
		_ambfix = nullptr;
	}
	if (_gModel_base)
	{
		delete _gModel_base;
		_gModel_base = nullptr;
	}

}

int gsins::t_gpvtflt_old::_addObsD(t_gsatdata& satdata, unsigned int& iobs, t_gallpar& param, t_gtriple& XYZ, Matrix& A, ColumnVector& l, DiagonalMatrix& P)
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

	t_gobs gobs1, gobs2, gobs3, gobs4, gobs5;

	_getgobs(satdata.sat(), TYPE_D, b1, gobs1);
	_getgobs(satdata.sat(), TYPE_D, b2, gobs2);
	_getgobs(satdata.sat(), TYPE_D, b3, gobs3);
	_getgobs(satdata.sat(), TYPE_D, b4, gobs4);
	_getgobs(satdata.sat(), TYPE_D, b5, gobs5);

	double D1, D2, D3, D4, D5;
	D1 = D2 = D3 = D4 = D5 = 0.0;
	vector<int> index_freq(5, 0);
	D1 = satdata.obs_D(gobs1);
	D2 = satdata.obs_D(gobs2);
	if (b3 != BAND) D3 = satdata.obs_D(gobs3);
	if (b4 != BAND) D4 = satdata.obs_D(gobs4);
	if (b5 != BAND) D5 = satdata.obs_D(gobs5);

	double weight_coef = _weightObs(satdata, gobs1);
	int obsnum = 0;
	if (!double_eq(D1, 0)) { obsnum++; index_freq[0] = 1; }
	if (!double_eq(D2, 0)) { obsnum++; index_freq[1] = 1; }
	if (!double_eq(D3, 0)) { obsnum++; index_freq[2] = 1; }
	if (!double_eq(D4, 0)) { obsnum++; index_freq[3] = 1; }
	if (!double_eq(D5, 0)) { obsnum++; index_freq[4] = 1; }

	if (obsnum == 0)return -1;

	double sigDoppler = 0.0;
	switch (gs) {
	case GPS: sigDoppler = _sigDopplerGPS; break;
	case GLO: sigDoppler = _sigDopplerGLO; break;
	case GAL: sigDoppler = _sigDopplerGAL; break;
	case BDS: sigDoppler = _sigDopplerBDS; break;
	case QZS: sigDoppler = _sigDopplerQZS; break;
	default: sigDoppler = 0.0;
	}

	//Create reduced measurements (prefit residuals)
	double modObsD = 0.0;
	int n = obsnum;
	for (int i = 0; i < obsnum; i++)
	{
		double Di = 0.0;
		if (i == 0) { modObsD = _gModel->cmpObsD(_epoch, param, satdata, gobs1); Di = D1; }
		if (i == 1) { modObsD = _gModel->cmpObsD(_epoch, param, satdata, gobs2); Di = D2; }
		if (i == 2) { modObsD = _gModel->cmpObsD(_epoch, param, satdata, gobs3); Di = D3; }
		if (i == 3) { modObsD = _gModel->cmpObsD(_epoch, param, satdata, gobs4); Di = D4; }
		if (i == 4) { modObsD = _gModel->cmpObsD(_epoch, param, satdata, gobs5); Di = D5; }
		if (double_eq(modObsD, 0.0) || double_eq(Di, 0.0)) continue;
		l(iobs + i) = -Di - modObsD;
	}

	// Create weight matrix
	for (int i = 0; i < obsnum; i++)
		P(iobs + i, iobs + i) = weight_coef * (1 / (sigDoppler * sigDoppler));

	// Create first design matrix
	double A_conf;
	int i = param.getParam(_site, t_gpar::VEL_X, "");
	int j = param.getParam(_site, t_gpar::VEL_Y, "");
	int k = param.getParam(_site, t_gpar::VEL_Z, "");
	t_gtriple groundVel(param[i].value(), param[j].value(), param[k].value());

	for (int i = 0; i < obsnum; i++)
	{
		for (unsigned int ipar = 1; ipar <= param.parNumber(); ipar++)
		{
			if (i == 0) A_conf = param[ipar - 1].partial_doppler(satdata, XYZ, groundVel);
			if (i == 1) A_conf = param[ipar - 1].partial_doppler(satdata, XYZ, groundVel);
			if (i == 2) A_conf = param[ipar - 1].partial_doppler(satdata, XYZ, groundVel);
			if (i == 3) A_conf = param[ipar - 1].partial_doppler(satdata, XYZ, groundVel);
			if (i == 4) A_conf = param[ipar - 1].partial_doppler(satdata, XYZ, groundVel);
			if (double_eq(A_conf, 0.0))continue;
			A(iobs + i, ipar) = A_conf;
		}
	}

	iobs = iobs + obsnum;

	return 1;
}

int gsins::t_gpvtflt_old::_rtkinit()
{

	/* setting initiating */
	_param.delAllParam();
	_Qx.Release();
	int ipar = 0;

	// Add coordinates parameters
	if (_crd_est != FIX) {
		_param.addParam(t_gpar(_site, t_gpar::CRD_X, ++ipar, ""));
		_param.addParam(t_gpar(_site, t_gpar::CRD_Y, ++ipar, ""));
		_param.addParam(t_gpar(_site, t_gpar::CRD_Z, ++ipar, ""));
	}

	// Add tropospheric wet delay parameter
	if (_tropo_est) {
		//rover 
		t_gpar trp_rover(_site, t_gpar::TRP, ++ipar, "");
		trp_rover.setMF(_ztd_mf);
		_param.addParam(trp_rover);
		//base
		t_gpar trp_base(_site_base, t_gpar::TRP, ++ipar, "");
		trp_base.setMF(_ztd_mf);
		_param.addParam(trp_base);
	}

	// Filling init parameter covariance matrix
	_Qx.ReSize(_param.parNumber()); _Qx = 0.0;
	double crdInit = _sig_init_crd;
	double ztdInit = _sig_init_ztd;
	for (unsigned int i = 1; i <= _param.parNumber(); i++)
	{
		if (_param[i - 1].parType == t_gpar::CRD_X)
			_Qx(i, i) = crdInit * crdInit;
		else if (_param[i - 1].parType == t_gpar::CRD_Y)
			_Qx(i, i) = crdInit * crdInit;
		else if (_param[i - 1].parType == t_gpar::CRD_Z)
			_Qx(i, i) = crdInit * crdInit;
		else if (_param[i - 1].parType == t_gpar::TRP)
			_Qx(i, i) = ztdInit * ztdInit;
	}

	return 1;
}

int gsins::t_gpvtflt_old::_setCrd()
{
	//save base coordinate to grec
	shared_ptr<gnut::t_gobj> grec_base = _gallobj->obj(_site_base);
	/* rover/reference fixed position */
	if (grec_base == nullptr) return -1;

	_gModel_base = new t_gpppmodel(_site_base, _log, _set);
	_gModel_base->reset_observ(_observ);
	_gModel_base->setOBJ(_gallobj);

	//get base coordinate
	BASEPOS basepos = dynamic_cast<t_gsetproc*>(_set)->basepos();
	t_gtime beg = dynamic_cast<t_gsetgen*>(_set)->beg();
	t_gtime end = dynamic_cast<t_gsetgen*>(_set)->end();

	if (basepos == SPP)
	{// from spp
		/*
		t_gspplsq* _spplsq = new t_gspplsq(_site_base, _set, _gallobj);
		_spplsq->setDAT(_gobs, _gnav);

		//if (_realtime == true)beg = gobs->beg_obs(_site);
		if (_spplsq->processBatch(beg, end) <= 0) {
			if (_log) {
				_log->comment(0, _site_base + "spplsq failed");
			}
			delete _spplsq; _spplsq = nullptr;
			return -1;
		}
		t_gtriple crd = _spplsq->getAvergeCrd();
		if (_spplsq) {
			delete _spplsq;
			_spplsq = nullptr;
		}
		grec_base->crd(crd, t_gtriple(1, 1, 1), beg, end);*/
	}
	else if (basepos == CFILE)
	{
		// set the value [read from xml/snx/rnxo]
		//vector<double> crd = dynamic_cast<t_gsetgen*>(_set)->crd(_site_base);
		auto crd = dynamic_cast<t_gsetrec*>(_set)->get_crd_xyz(_site_base);
		grec_base->crd(t_gtriple(crd[0], crd[1], crd[2]), t_gtriple(0.01, 0.01, 0.01), beg, end, true);
	}

	if (_crd_est == FIX && !_pos_kin) {
		// set the value [read from xml/snx/rnxo]
		//vector<double> crd = dynamic_cast<t_gsetgen*>(_set)->crd(_site);
		auto crd = dynamic_cast<t_gsetrec*>(_set)->get_crd_xyz(_site);
		_grec->crd(t_gtriple(crd[0], crd[1], crd[2]), t_gtriple(0.01, 0.01, 0.01), beg, end, true);
	}

	return 1;
}

/* detect valid zero-difference residuals*/
bool gsins::t_gpvtflt_old::_valid_residual(bool phase_process, string sat_name, enum FREQ_SEQ &f, map<string, map<FREQ_SEQ, pair<int, int>>>& index_l) {
	try
	{
		/* if no phase observable, psudorange is also unusable */
		if (phase_process)
		{
			auto tmp_phase = index_l.at(sat_name).at(_observ == IONO_FREE ? FREQ_1 : f).second;
		}
		else {
			auto tmp_phase = index_l.at(sat_name).at(_observ == IONO_FREE ? FREQ_1 : f).second;
			auto tmp_code = index_l.at(sat_name).at(_observ == IONO_FREE ? FREQ_1 : f).first;
		}
	}
	catch (exception e)
	{
		return false;
	}
	return true;
}

int gsins::t_gpvtflt_old::_combineDD(Matrix & A, SymmetricMatrix & P, ColumnVector & l)
{
	gtrace("t_gpvtflt_newmodel::_combineDD");
	try
	{
		// ����index_l todo 
		map<string, map<FREQ_SEQ, pair<int, int>>> index_l;
		for (int i = 0; i < _obs_index.size(); i++)
		{
			string sat = _obs_index[i].first;
			FREQ_SEQ f = _obs_index[i].second.first;
			GOBSTYPE obstype = _obs_index[i].second.second;
			if (index_l.find(sat) == index_l.end())
			{
				index_l.insert(make_pair(sat, map<FREQ_SEQ, pair<int, int>>()));
			}
			if (index_l[sat].find(f) == index_l[sat].end())
			{
				index_l[sat].insert(make_pair(f, make_pair(-1, -1)));
			}
			if (obstype == TYPE_C) index_l[sat][f].first = i + 1;
			else if (obstype == TYPE_L) index_l[sat][f].second = i + 1;
		}

		_obs_index.clear();
		int iobs = 1;
		set<string> sysall = dynamic_cast<t_gsetgen*>(_set)->sys();

		int nrows = A.Nrows();
		Matrix DD; DD.ReSize(nrows, nrows); DD = 0.0;

		bool isPhaseProcess = true;
		for (auto sys_iter = sysall.begin(); sys_iter != sysall.end(); sys_iter++)
		{
			enum GSYS sys = t_gsys::str2gsys(*sys_iter);
			vector < GOBSBAND > band = dynamic_cast<t_gsetgnss*>(_set)->band(sys);
			int nf = 5;
			if (band.size()) nf = band.size();
			if (_observ == IONO_FREE) nf = 1;
			FREQ_SEQ f;
			for (FREQ_SEQ freq = FREQ_1; freq <= 2 * nf; freq = (FREQ_SEQ)(freq + 1)) {
				if (freq <= nf) {//phase equations
					isPhaseProcess = true;
					f = freq;
				}
				else {//code equations
					isPhaseProcess = false;
					f = (FREQ_SEQ)(freq - nf);
				}
				string sat, sat_ref;
				t_gsatdata obs_sat_ref; enum GSYS gs;
				for (auto it = _data.begin(); it != _data.end(); it++)
				{
					string sat = it->sat();

					gs = it->gsys();
					if (gs == QZS) gs = GPS;
					if (gs != sys) continue;
					//if (gs == BDS && sat <= "C05") continue;

					if (_valid_residual(isPhaseProcess, sat, f, index_l) == false)continue;

					if (sat_ref.empty()) {
						sat_ref = sat;
						obs_sat_ref = *it;
						continue;
					}

					if (it->ele_deg() >= obs_sat_ref.ele_deg())
					{
						sat_ref = sat;
						obs_sat_ref = *it;
					}
				}//end select sat_ref

				if (sat_ref.empty()) continue;

				int index_ref;
				if (isPhaseProcess)	index_ref = index_l[sat_ref][f].second;
				else index_ref = index_l[sat_ref][f].first;
				/* make double difference */
				for (auto it = _data.begin(); it != _data.end(); it++)
				{
					sat = it->sat();
					if (sat == sat_ref)continue;
					gs = it->gsys();
					if (gs == QZS) gs = GPS;
					if (gs != sys) continue;

					if (_valid_residual(isPhaseProcess, sat, f, index_l) == false) continue;

					int index_sat;
					if (isPhaseProcess)	index_sat = index_l[sat][f].second;
					else index_sat = index_l[sat][f].first;
					DD(iobs, index_ref) = -1;
					DD(iobs, index_sat) = 1;
					GOBSTYPE obstype = TYPE_C;
					if (isPhaseProcess)obstype = TYPE_L;
					_obs_index.push_back(make_pair(it->sat(), make_pair(f, obstype)));
					iobs++;
				}//end sat
			}//end f
		}//end sys

		// delete zero rows/cols
		iobs--;
		DD = DD.Rows(1, iobs);

		A = DD * A;
		l = DD * l;
		P << (DD * P.i() * DD.t()).i();
		// cout << _epoch.str_ymdhms()<< endl;
		// t_out("A", A);
		// t_out("L", l);
		// t_out("p", P);

		return 1;
	}
	catch (...)
	{
		if (_log) _log->comment(1, "t_gpvtflt_newmodel", "generate double-differented equation failed!");
		return -1;
	}
}

int gsins::t_gpvtflt_old::_prepareData()
{
	gtrace("t_gpvtflt::_prepareData");
	int flag = 1;
	if (_isBase) {
		if (_prepareData(_site_base, _data_base) < 0) {
			if (_log) _log->comment(2, _site + "gpvtflt", "Preparing data failed!");
			flag = -1;
		}
	}
	if (_prepareData(_site, _data) < 0) {
		if (_log) _log->comment(2, _site + "gpvtflt", "Preparing data failed!");
		flag = -1;
	}

	if (_isBase) {
		if ((_nSat = _selcomsat(_data_base, _data)) < _minsat) {
			if (_log) _log->comment(2, "gpvtflt", "no common satellite");
			flag = -1;
		}
	}

	return flag;
}

int gsins::t_gpvtflt_old::_prepareData(const string & ssite, vector<t_gsatdata>& sdata)
{
	gtrace("t_gpvtflt::_prepareData");

	vector<t_gsatdata>::iterator it = sdata.begin();
	while (it != sdata.end()) {
		switch (_gnss) {
		case GNS: break;
		case GPS: if (it->gsys() != GPS) { it = sdata.erase(it); continue; } break;
		case GLO: if (it->gsys() != GLO) { it = sdata.erase(it); continue; } break;
		case GAL: if (it->gsys() != GAL) { it = sdata.erase(it); continue; } break;
		case BDS: if (it->gsys() != BDS) { it = sdata.erase(it); continue; } break;
		case QZS: if (it->gsys() != QZS) { it = sdata.erase(it); continue; } break;
		case IRN: if (it->gsys() != IRN) { it = sdata.erase(it); continue; } break;
		case SBS: if (it->gsys() != SBS) { it = sdata.erase(it); continue; } break;
		}
		it++;
	}

	Matrix BB;

	BB.ReSize(sdata.size(), 4);
	BB = 0.0;
	int iobs = 0;
	vector<t_gsatdata>::iterator iter = sdata.begin();

	_nSat = sdata.size();  // configured GNSS is considered (see erase in above swich)

	while (iter != sdata.end()) {

		GSYS gs = iter->gsys();

		GOBSBAND b1, b2;

		if (_auto_band) {
			set<GOBSBAND> bands;
			bands = iter->band_avail(_phase);
			auto itBAND = bands.begin();
			if (bands.size() < 2) { iter++; continue; }
			b1 = *itBAND;
			itBAND++;
			b2 = *itBAND;
			//std::cout << iter->sat() << " " << b1 << " " << b2 << std::endl;
		}
		else {
			vector < GOBSBAND > band = dynamic_cast<t_gsetgnss*>(_set)->band(gs);
			if (band.size())
			{
				b1 = band[0];
				b2 = band[1];
			}
			else {
				b1 = t_gsys::band_priority(gs, FREQ_1);
				b2 = t_gsys::band_priority(gs, FREQ_2);
			}

		}

		iter->glog(_log);

		// check data availability lvhb added in 20200525
		auto l1 = iter->select_phase(b1);
		auto l2 = iter->select_phase(b2);
		auto p1 = iter->select_range(b1);
		auto p2 = iter->select_range(b2);
		if (p1 == X || p2 == X || l1 == X || l2 == X) {
			iter = sdata.erase(iter);
			continue;
		}

		// check data availability
		double P3 = iter->P3(b1, b2);
		double L3 = iter->L3(b1, b2);
		if (_observ == RAW_SINGLE) {
			P3 = iter->obs_C(t_gband(b1, GOBSATTR::ATTR));
			L3 = iter->obs_L(t_gband(b1, GOBSATTR::ATTR));
		}

#ifdef DEBUG
		cerr << "iter->epoch/_epoch: " << iter->epoch().str_hms() << _epoch.str_hms(" / ")
			<< " " << t_gsys::gsys2str(gs) << " " << iter->sat() << " " << iter->channel()
			<< " " << fixed << setprecision(3) << P3 << " " << L3 << endl;
#endif

		if (double_eq(L3, 0.0) && _phase) {
			if (_log) {
				ostringstream str;
				str << "prepareData: erasing data due to no phase double bands observation, " <<
					"epo: " << _epoch.str_hms() << ", " <<
					"prn: " << iter->sat() << " (" << iter->channel() << ")";
				_log->comment(2, "gsppflt", str.str());
			}
			iter = sdata.erase(iter); // !!!! OK u vektoru funguje bez ++ (erase zvysuje pointer !)      
			continue;
		}

		if (double_eq(P3, 0.0)) {
			if (_log) {
				ostringstream str;
				str << "prepareData: erasing data due to no code double bands observation, " <<
					"epo: " << _epoch.str_hms() << ", " <<
					"prn: " << iter->sat();
				_log->comment(2, "gsppflt", str.str());
			}
			iter = sdata.erase(iter); // !!!! OK u vektoru funguje bez ++ (erase zvysuje pointer !)      
			continue;
		}

		if (_satPos(_epoch, *iter) < 0) {
			if (_log) {
				ostringstream str;
				str << "prepareData: erasing data since _satPos failed, "
					<< "epo: " << _epoch.str_hms() << ", " <<
					"prn: " << iter->sat();
				_log->comment(2, "gsppflt", str.str());
			}
#ifdef DEBUG
			cerr << "prepareData: sat pos not calculated "
				<< _epoch.str_ymdhms(ssite + " " + iter->sat() + " " + int2str(sdata.size()) + " ")
				<< endl; cerr.flush();
#endif       
			iter = sdata.erase(iter); // !!!! OK u vektoru funguje bez ++ (erase zvysuje pointer !)

		}
		else {
			iobs++;
			BB(iobs, 1) = iter->satcrd().crd(0);
			BB(iobs, 2) = iter->satcrd().crd(1);
			BB(iobs, 3) = iter->satcrd().crd(2);
			BB(iobs, 4) = P3 + iter->clk();

			// compute the MW combine for ambfix (added by zhshen)
			//if(_observ == IONO_FREE)	_combineMW(*iter);

			iter++;
			continue;
		}
	}//end sdata

	if (sdata.size() < _minsat) {
		if (_log)
			_log->comment(2, "gsppflt", _epoch.str_ymdhms(ssite + " epoch ")
				+ " skipped (Bancroft not calculated: " + int2str(sdata.size()) + " < _minsat)");
		//      cout << t_gsys::gsys2str(_gnss) << " " << _epoch.str_ymdhms( _site + " epoch ") << " " << int2str(_data.size()) <<  " skipped (Bancroft not calculated: _data.size < _minsat)" << endl;
		return -1;
	}

	BB = BB.Rows(1, iobs);    // delete zero rows

	if (BB.Nrows() < static_cast<int>(_minsat)) {
		if (_log && _log->verb() >= 2)
			_log->comment(2, "gsppflt", _epoch.str_ymdhms(_site + " epoch ")
				+ " skipped (Bancroft not calculated: BB.Nrows < _minsat)");
		//      cout << _epoch.str_ymdhms( _site + " epoch ") << int2str(_data.size()) <<  " skipped (Bancroft not calculated: _data.size < _minsat)" << endl;
		return -1;
	}

	_vBanc = 0.0;
	gbancroft(BB, _vBanc);

	shared_ptr<t_gobj> grec = _gallobj->obj(ssite);
	if (grec == nullptr) return -1;

	if (!_initialized && (_vBanc.Rows(1, 3) - grec->crd_arp(_epoch).crd_cvect()).NormFrobenius() > 1000) _valid_crd_xml = false;

	if (_valid_crd_xml && !_initialized && (_crd_est == FIX && ssite == _site)) {
		_vBanc(1) = grec->crd_arp(_epoch)[0];
		_vBanc(2) = grec->crd_arp(_epoch)[1];
		_vBanc(3) = grec->crd_arp(_epoch)[2];
	}
	if (ssite == _site && _crd_est != FIX) {
		if (_pos_constrain)
		{
			_vBanc(1) = _extn_pos[0];
			_vBanc(2) = _extn_pos[1];
			_vBanc(3) = _extn_pos[2];
			_pos_constrain = false;
		}
		_grec->crd(t_gtriple(_vBanc), t_gtriple(1, 1, 1), _epoch, LAST_TIME);//modified by zhshen
	}
	//if (_isBase && ssite == _site_base)	_vBanck_base = _vBanc;

	t_gtriple test_xyz(_vBanc);
	t_gtriple test_ell;
	xyz2ell(test_xyz, test_ell, true);
	double radius = test_xyz.norm();
	if (radius < B_WGS - 500) {
		string warning = "WARNING: Unexpected site (" + ssite + ") coordinates from Bancroft. Orbits/clocks or code observations should be checked.";
		if (_log)  _log->comment(1, "gsppflt", warning + _epoch.str_ymdhms(" Epoch "));
		if (_glog) _glog->comment(2, "gsppflt", warning + _epoch.str_ymdhms(" Epoch "));
		if (!_phase) { return -1; }
	}

#ifdef DEBUG
	cout << _epoch.str_ymdhms() << endl; cout.flush();
	//   cout << "\nBB\n" << fixed << setprecision(10) << BB << endl;
	cout << "vBanc\n" << fixed << setprecision(10) << _vBanc << endl; cout.flush();
	//   cout << "Clk [mili sec]" << fixed << setprecision(10) << _vBanc(4)*1000/CLIGHT << endl;
	//   cout << "_Qx\n"  << fixed << setprecision(3)<< _Qx << endl;    cout.flush();
	//   int ooo; cin >> ooo;
#endif

	//Compute sat elevation and rho
	iter = sdata.begin();
	while (iter != sdata.end()) {
		t_gtriple xyz_r, xyz_s;
		xyz_s = iter->satcrd();

		if (_isBase && ssite == _site_base) {//lvhb added in 20200523
			xyz_r = _gallobj->obj(ssite)->crd_arp(_epoch);
		}
		else {
			if (_crd_est == FIX && _valid_crd_xml) {
				xyz_r = _gallobj->obj(ssite)->crd_arp(_epoch);
			}
			else if (_crd_est == EST) {
				if (!_initialized) {
					xyz_r.set(_vBanc.Rows(1, 3));
				}
				else {
					if (_pos_kin) {
						xyz_r.set(_vBanc.Rows(1, 3));
					}
					else _param.getCrdParam(ssite, xyz_r);
				}

			}
			else {
				cerr << "Undefined CRD constraining" << endl;
				return -1;
			}
		}
		Add_rho_azel(/*runEpoch,*/ ssite, xyz_s, xyz_r, *iter); //lvhb added in 20200523

		// check elevation cut-off 
		if (iter->ele_deg() < _minElev) {
			ostringstream os;
			os << "Erasing " << iter->sat() << " data due to low elevation angle (ele = " << fixed << setprecision(1) << iter->ele_deg()
				<< ") " << iter->epoch().str_ymdhms();
			if (_log) _log->comment(2, "gpvtflt", os.str());
			iter = sdata.erase(iter); // !!!! zveda iterator !!!!
			if (sdata.size() < _minsat) {
				if (_log) _log->comment(2, "", ssite + _epoch.str_ymdhms(" epoch ") + int2str(sdata.size()) + " skipped (data.size < _rtk_set->minsat)");
				return -1;
			}
			continue;
		}

		// Printing beta and orbit angles for deep verbosity
		ostringstream os;
		os << iter->sat() << " " << iter->epoch().str_ymdhms() << " " << fixed << setprecision(1) << iter->beta() << " " << iter->orb_angle() << endl;
		if (_log) _log->comment(4, "grtkflt", os.str());

		// check satellite eclipsing
		if (iter->ecl()) {
			ostringstream os;
			os << "Erasing " << iter->sat() << " data due to satellite eclipsing (beta = " << fixed << setprecision(1) << iter->beta()
				<< " ,orbit angle = " << iter->orb_angle() << ") " << iter->epoch().str_ymdhms();
			if (_log) _log->comment(2, "grtkflt", os.str());
			iter = sdata.erase(iter); // !!!! zveda iterator !!!!
			if (sdata.size() < _minsat) {
				if (_log) _log->comment(2, "", ssite + _epoch.str_ymdhms(" epoch ") + int2str(sdata.size()) + " skipped (data.size < _rtk_set->minsat)");
				return -1;
			}
		}
		else iter++;

	}//end sdata

	return 0;

}

void gsins::t_gpvtflt_old::_predict()
{
	gtrace("t_gpvtflt::_predict");

	// Predict coordinates, clock and troposphere
	// state vector is the same, covariance matrix is predicted with noise   
	int i;

	_cntrep++;

	// add/remove inter system bias
	if (!_isBase) _syncSys();

	// add/remove ionosphere delay
	_syncIono();

	// add/remove inter-frequency biases
	if (!_isBase) _syncIFB();

	_Noise.ReSize(_Qx.Nrows()); _Noise = 0;

	// Predict coordinates, clock and troposphere
	// state vector is the same, covariance matrix is predicted with noise   

	double crdInit = _sig_init_crd;
	double ztdInit = _sig_init_ztd;

	i = _param.getParam(_site, t_gpar::CRD_X, "");
	if (i >= 0)
	{
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0)
		{
			_param[i].value(_vBanc(1));
			_Qx(i + 1, i + 1) = crdInit * crdInit;
		}
		else
		{
			if (_epoch > _crd_begStat && _epoch < _crd_endStat)
			{
				if (_smooth) _Noise(i + 1, i + 1) = 0.0;
			}
			else
			{
				if (_pos_kin) _param[i].value(_vBanc(1));
				if (_cntrep == 1 && _success) _Qx(i + 1, i + 1) += _crdStoModel->getQ() * _crdStoModel->getQ();
				if (_smooth) _Noise(i + 1, i + 1) = _crdStoModel->getQ() * _crdStoModel->getQ();
			}
		}
	}

	i = _param.getParam(_site, t_gpar::CRD_Y, "");
	if (i >= 0)
	{
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0)
		{
			_param[i].value(_vBanc(2));
			_Qx(i + 1, i + 1) = crdInit * crdInit;
		}
		else
		{
			if (_epoch > _crd_begStat && _epoch < _crd_endStat)
			{
				if (_smooth) _Noise(i + 1, i + 1) = 0.0;
			}
			else
			{
				if (_pos_kin) _param[i].value(_vBanc(2));
				if (_cntrep == 1 && _success) _Qx(i + 1, i + 1) += _crdStoModel->getQ() * _crdStoModel->getQ();
				if (_smooth) _Noise(i + 1, i + 1) = _crdStoModel->getQ() * _crdStoModel->getQ();
			}
		}
	}

	i = _param.getParam(_site, t_gpar::CRD_Z, "");
	if (i >= 0)
	{
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0)
		{
			_param[i].value(_vBanc(3));
			_Qx(i + 1, i + 1) = crdInit * crdInit;
		}
		else
		{
			if (_epoch > _crd_begStat && _epoch < _crd_endStat)
			{
				if (_smooth) _Noise(i + 1, i + 1) = 0.0;
			}
			else
			{
				if (_pos_kin) _param[i].value(_vBanc(3));
				if (_cntrep == 1 && _success) _Qx(i + 1, i + 1) += _crdStoModel->getQ() * _crdStoModel->getQ();
				if (_smooth) _Noise(i + 1, i + 1) = _crdStoModel->getQ() * _crdStoModel->getQ();
			}
		}
	}

	i = _param.getParam(_site, t_gpar::CLK, "");
	if (i >= 0)
	{
		_param[i].value(_vBanc(4));
		for (unsigned int jj = 1; jj <= _param.parNumber(); jj++) _Qx(i + 1, jj) = 0.0;
		_Qx(i + 1, i + 1) = _clkStoModel->getQ() * _clkStoModel->getQ();
		if (_smooth) _Noise(i + 1, i + 1) = _clkStoModel->getQ() * _clkStoModel->getQ();
	}

	i = _param.getParam(_site, t_gpar::GLO_ISB, "");
	if (i >= 0)
	{
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0)
		{
			_Qx(i + 1, i + 1) = _sig_init_glo * _sig_init_glo;
		}
		else
		{
			if (_cntrep == 1) _Qx(i + 1, i + 1) += _gloStoModel->getQ();
			if (_smooth) _Noise(i + 1, i + 1) = _gloStoModel->getQ();
		}
	}

	i = _param.getParam(_site, t_gpar::GLO_IFCB, "");
	if (i >= 0) {
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
			_Qx(i + 1, i + 1) = _sig_init_glo * _sig_init_glo;
		}
		else {
			if (_cntrep == 1) _Qx(i + 1, i + 1) += _gloStoModel->getQ();
			if (_smooth) _Noise(i + 1, i + 1) = _gloStoModel->getQ();
		}
	}

	i = _param.getParam(_site, t_gpar::GAL_ISB, "");
	if (i >= 0) {
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
			_Qx(i + 1, i + 1) = _sig_init_gal * _sig_init_gal;
		}
		else {
			if (_cntrep == 1) _Qx(i + 1, i + 1) += _galStoModel->getQ();
			if (_smooth) _Noise(i + 1, i + 1) = _galStoModel->getQ();
		}
	}

	i = _param.getParam(_site, t_gpar::BDS_ISB, "");
	if (i >= 0) {
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
			_Qx(i + 1, i + 1) = _sig_init_bds * _sig_init_bds;
		}
		else {
			if (_cntrep == 1) _Qx(i + 1, i + 1) += _bdsStoModel->getQ();
			if (_smooth) _Noise(i + 1, i + 1) = _bdsStoModel->getQ();
		}
	}

	i = _param.getParam(_site, t_gpar::QZS_ISB, "");
	if (i >= 0) {
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
			_Qx(i + 1, i + 1) = _sig_init_qzs * _sig_init_qzs;
		}
		else {
			if (_cntrep == 1) _Qx(i + 1, i + 1) += _qzsStoModel->getQ();
			if (_smooth) _Noise(i + 1, i + 1) = _qzsStoModel->getQ();
		}
	}

	i = _param.getParam(_site, t_gpar::TRP, "");
	if (i >= 0) {
		t_gtriple Ell, XYZ;

		if (_param.getCrdParam(_site, XYZ) > 0) {
		}
		else if (_valid_crd_xml) { XYZ = _grec->crd_arp(_epoch); }
		xyz2ell(XYZ, Ell, false);
		// ZTD apriori (=ZHD) is updated every epoch
		if (_gModel->tropoModel() != 0) {
			_param[i].apriori(_gModel->tropoModel()->getZHD(Ell, _epoch));
		
		}

		if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
			if (_valid_ztd_xml) {
				_param[i].value(_aprox_ztd_xml - _gModel->tropoModel()->getZHD(Ell, _epoch));
			}
			else {
				if (_gModel->tropoModel() != 0) {
					_param[i].value(_gModel->tropoModel()->getZWD(Ell, _epoch));
					
				}
			}

			_Qx(i + 1, i + 1) = ztdInit * ztdInit;
		}
		else {
			if (_epoch > _ztd_begStat && _epoch < _ztd_endStat) {
				if (_smooth) _Noise(i + 1, i + 1) = 0.0;
			}
			else {
				if (_cntrep == 1)  _Qx(i + 1, i + 1) += _trpStoModel->getQ();
				if (_smooth) _Noise(i + 1, i + 1) = _trpStoModel->getQ();
			}
		}
	}

	// predict ionosphere delay  
	vector<t_gsatdata>::iterator it;
	for (it = _data.begin(); it != _data.end(); it++)
	{
		i = _param.getParam(_site, t_gpar::VION, it->sat());
		if (i >= 0)
		{
			if (_gion)
			{   // Introduce new apriory value from IONEX if available
				t_gtriple site_xyz(0.0, 0.0, 0.0);
				t_gtriple site_ell(0.0, 0.0, 0.0);
				t_gtriple ipp_ell(0.0, 0.0, 0.0);

				_param.getCrdParam(_site, site_xyz);
				if (site_xyz.zero()) site_xyz.set(_vBanc);
				xyz2ell(site_xyz, site_ell, false);
				ell2ipp(*it, site_ell, ipp_ell);

				// use iono-free combination instead ! //
		//      double ionomodel = _gion->iono(ipp_ell[0] * R2D, ipp_ell[1] * R2D, it->epoch());
				double ionomodel = 1.0; // ( ionomodel * 40.28 * 1e16 ) / (G01_F*G01_F);

				_param[i].apriori(ionomodel);
			}

			if (_cntrep == 1 &&
				!double_eq(_Qx(i + 1, i + 1), _sig_init_vion*_sig_init_vion)) {
				_Qx(i + 1, i + 1) += _ionStoModel->getQ() * _ionStoModel->getQ();
			}

			// if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
			//   cout << "DEBUG " << it->sat() << " " << it->epoch().str_ymdhms() << endl;
			//   _Qx(i + 1, i + 1) =  _sig_init_vion * _sig_init_vion;
			// }
			// else {                
			//   if (_cntrep == 1) _Qx(i + 1, i + 1) += _ionStoModel->getQ() * _ionStoModel->getQ(); 
			//   if (_smooth) _Noise(i + 1, i + 1) = _ionStoModel->getQ() * _ionStoModel->getQ();
			// }
		}

		i = _param.getParam(_site, t_gpar::SION, it->sat());
		if (i >= 0)
		{
			if (_cntrep == 1 && !double_eq(_Qx(i + 1, i + 1), _sig_init_vion*_sig_init_vion)) {
				_Qx(i + 1, i + 1) += _ionStoModel->getQ() * _ionStoModel->getQ();
			}
		}
	}

	// predict inter-frequency clock bias
	for (it = _data.begin(); it != _data.end(); it++)
	{
		i = _param.getParam(_site, t_gpar::IFCB_F3, it->sat());
		if (i >= 0)
		{
			if (_cntrep == 1) _Qx(i + 1, i + 1) += 0.001;
			if (_smooth) _Noise(i + 1, i + 1) = 0.001;
		}
		i = _param.getParam(_site, t_gpar::IFCB_F4, it->sat());
		if (i >= 0)
		{
			if (_cntrep == 1) _Qx(i + 1, i + 1) += 0.001;
			if (_smooth)      _Noise(i + 1, i + 1) = 0.001;
		}
		i = _param.getParam(_site, t_gpar::IFCB_F5, it->sat());
		if (i >= 0)
		{
			if (_cntrep == 1) _Qx(i + 1, i + 1) += 0.001;
			if (_smooth)      _Noise(i + 1, i + 1) = 0.001;
		}
	}


	double velInit = _sig_init_vel;

	i = _param.getParam(_site, t_gpar::VEL_X, "");
	if (i >= 0) {
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
			_param[i].value(0);
			_Qx(i + 1, i + 1) = velInit * velInit;
		}
		else {
			if (_epoch > _crd_begStat&& _epoch < _crd_endStat) {
				if (_smooth) _Noise(i + 1, i + 1) = 0.0;
			}
			else {
				if (_cntrep == 1 && _success) _Qx(i + 1, i + 1) += _velStoModel->getQ() * _velStoModel->getQ();
				if (_smooth) _Noise(i + 1, i + 1) = _velStoModel->getQ() * _velStoModel->getQ();
			}
		}
	}

	i = _param.getParam(_site, t_gpar::VEL_Y, "");
	if (i >= 0) {
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
			_param[i].value(0);
			_Qx(i + 1, i + 1) = velInit * velInit;
		}
		else {
			if (_epoch > _crd_begStat&& _epoch < _crd_endStat) {
				if (_smooth) _Noise(i + 1, i + 1) = 0.0;
			}
			else {
				if (_cntrep == 1 && _success) _Qx(i + 1, i + 1) += _velStoModel->getQ() * _velStoModel->getQ();
				if (_smooth) _Noise(i + 1, i + 1) = _velStoModel->getQ() * _velStoModel->getQ();
			}
		}
	}

	i = _param.getParam(_site, t_gpar::VEL_Z, "");
	if (i >= 0) {
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
			_param[i].value(0);
			_Qx(i + 1, i + 1) = velInit * velInit;
		}
		else {
			if (_epoch > _crd_begStat&& _epoch < _crd_endStat) {
				if (_smooth) _Noise(i + 1, i + 1) = 0.0;
			}
			else {
				if (_cntrep == 1 && _success) _Qx(i + 1, i + 1) += _velStoModel->getQ() * _velStoModel->getQ();
				if (_smooth) _Noise(i + 1, i + 1) = _velStoModel->getQ() * _velStoModel->getQ();
			}
		}
	}

	i = _param.getParam(_site, t_gpar::CLK_RAT, "");
	if (i >= 0) {
		//_param[i].value(0);
		for (unsigned int jj = 1; jj <= _param.parNumber(); jj++) _Qx(i + 1, jj) = 0.0;
		_Qx(i + 1, i + 1) = _dclkStoModel->getQ() * _dclkStoModel->getQ();
		if (_smooth) _Noise(i + 1, i + 1) = _dclkStoModel->getQ() * _dclkStoModel->getQ();
	}


	i = _param.getParam(_site, t_gpar::GRD_N, "");
	if (i >= 0) {
		if (_cntrep == 1 && _initialized) _Qx(i + 1, i + 1) += _grdStoModel->getQ();
		if (_smooth) _Noise(i + 1, i + 1) = _grdStoModel->getQ();
	}

	i = _param.getParam(_site, t_gpar::GRD_E, "");
	if (i >= 0) {
		if (_cntrep == 1 && _initialized) _Qx(i + 1, i + 1) += _grdStoModel->getQ();
		if (_smooth) _Noise(i + 1, i + 1) = _grdStoModel->getQ();
	}

	// Predict ambiguity 
	// Add or remove ambiguity parameter and appropriate rows/columns covar. matrix
	if (_phase) {
		if (_isBase)_udsdAmb();
		else _syncAmb();
	}

	if (_isBase) {
		//tropspheric parameter
		if (_tropo_est) {
			//rover: sppflt::predict()
			//base
			//_site,_gmodel,_grec,_epoch;
			t_gpppmodel* gModel = new t_gpppmodel(_site_base, _log, _set);
			gModel->reset_observ(_observ);
			gModel->setOBJ(_gallobj);
			shared_ptr<t_gobj> grec = _gallobj->obj(_site_base);
			if (grec == nullptr) return;

			double ztdInit = _sig_init_ztd;
			int i = _param.getParam(_site_base, t_gpar::TRP, "");
			if (i >= 0) {
				t_gtriple Ell, XYZ;

				if (_param.getCrdParam(_site_base, XYZ) > 0) {
				}
				else { XYZ = grec->crd_arp(_epoch); }
				xyz2ell(XYZ, Ell, false);
				// ZTD apriori (=ZHD) is updated every epoch
				if (gModel->tropoModel() != 0) {
					_param[i].apriori(gModel->tropoModel()->getZHD(Ell, _epoch));
				}

				if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
					if (_valid_ztd_xml) {
						_param[i].value(_aprox_ztd_xml - gModel->tropoModel()->getZHD(Ell, _epoch));
					}
					else {
						if (gModel->tropoModel() != 0) {
							_param[i].value(gModel->tropoModel()->getZWD(Ell, _epoch));
						}
					}

					_Qx(i + 1, i + 1) = ztdInit * ztdInit;
				}
				else {
					if (_epoch > _ztd_begStat && _epoch < _ztd_endStat) {
						if (_smooth) _Noise(i + 1, i + 1) = 0.0;
					}
					else {
						if (_cntrep == 1)  _Qx(i + 1, i + 1) += _trpStoModel->getQ();
						if (_smooth) _Noise(i + 1, i + 1) = _trpStoModel->getQ();
					}
				}
			}

			delete gModel;
		}
	}

	// Create _Noise matrix (DiagonalMatrix)
	if (_smooth) {
		Matrix tmp = _Noise;
		_Noise.ReSize(_Qx.Nrows()); _Noise = 0;
		for (int i = 1; i <= _Noise.Ncols(); i++) {
			if (i > tmp.Ncols()) break;
			_Noise(i, i) = tmp(i, i);
		}
	}
}

int gsins::t_gpvtflt_old::_addObsP(t_gsatdata & satdata, unsigned int & iobs, Matrix & A, ColumnVector & l, SymmetricMatrix & P)
{
	gtrace("t_gpvtflt_old::_addObsP");

	t_gsatdata satdata_other;

	if (_isBase)
	{
		for (int i = 0; i < _data_base.size(); i++)
		{
			if (_data_base[i].sat() == satdata.sat())
				satdata_other = _data_base[i];
		}
		if (satdata_other.site() == "") return -1;
	}

	// 获取流动站和基站的坐�?
	t_gtriple xyz, ell;
	if (_param.getCrdParam(_site, xyz) < 0)
	{
		xyz = _gallobj->obj(_site)->crd_arp(_epoch);
	}
	xyz2ell(xyz.crd_array(), ell.crd_array(), false);

	t_gtriple xyz_base, ell_base;
	if (_isBase)
	{
		if (_param.getCrdParam(_site_base, xyz_base) < 0)
		{
			xyz_base = _gallobj->obj(_site_base)->crd_arp(_epoch);
		}
		xyz2ell(xyz_base.crd_array(), ell_base.crd_array(), false);
	}

	// 根据流动站数据选择频率和观测值类�?
	t_gsys gsys(satdata.gsys());
	GSYS gs = gsys.gsys();

	GOBSBAND b1, b2, b3, b4, b5;
	b1 = b2 = b3 = b4 = b5 = BAND;

	if (_auto_band) {    // automatic dual band selection -> for Anubis purpose
		set<GOBSBAND> bands = satdata.band_avail(_phase);
		auto itBAND = bands.begin();
		if (bands.size() < 2) return -1;
		b1 = *itBAND;
		itBAND++;
		b2 = *itBAND;
	}
	else
	{
		// select fix defined band according the table
		vector < GOBSBAND > band = dynamic_cast<t_gsetgnss*>(_set)->band(gs);
		if (band.size())
		{
			b1 = band[0];
			b2 = band[1];
			if (band.size() > 2) b3 = band[2];
			if (band.size() > 3) b4 = band[3];
			if (band.size() > 4) b5 = band[4];
		}
		else {
			b1 = t_gsys::band_priority(gs, FREQ_1);
			b2 = t_gsys::band_priority(gs, FREQ_2);
			b3 = t_gsys::band_priority(gs, FREQ_3);
			b4 = t_gsys::band_priority(gs, FREQ_4);
			b5 = t_gsys::band_priority(gs, FREQ_5);
		}
	}

	t_gobs gobs1, gobs2, gobs3, gobs4, gobs5;
	_getgobs(satdata.sat(), TYPE_C, b1, gobs1);
	_getgobs(satdata.sat(), TYPE_C, b2, gobs2);
	_getgobs(satdata.sat(), TYPE_C, b3, gobs3);
	_getgobs(satdata.sat(), TYPE_C, b4, gobs4);
	_getgobs(satdata.sat(), TYPE_C, b5, gobs5);


	// 收集流动站的观测�?
	double PIF, P1, P2, P3, P4, P5;
	PIF = P1 = P2 = P3 = P4 = P5 = 0.0;

	PIF = satdata.P3(gobs1, gobs2);
	P1 = satdata.obs_C(gobs1);
	P2 = satdata.obs_C(gobs2);
	if (b3 != BAND) P3 = satdata.obs_C(gobs3);
	if (b4 != BAND) P4 = satdata.obs_C(gobs4);
	if (b5 != BAND) P5 = satdata.obs_C(gobs5);

	if (_observ == IONO_FREE && double_eq(PIF, 0.0)) return -1;
	if (_observ == RAW_SINGLE && double_eq(P1, 0.0)) return -1;
	if ((_observ == RAW_DOUBLE || _observ == RAW_ALL) && (double_eq(P1, 0.0) || double_eq(P2, 0.0))) return -1;

	// 收集基站的观测�?
	double PIF_base, P1_base, P2_base, P3_base, P4_base, P5_base;
	PIF_base = P1_base = P2_base = P3_base = P4_base = P5_base = 0.0;

	if (_isBase)
	{
		PIF_base = satdata_other.P3(gobs1, gobs2);
		P1_base = satdata_other.obs_C(gobs1);
		P2_base = satdata_other.obs_C(gobs2);
		if (b3 != BAND) P3_base = satdata_other.obs_C(gobs3);
		if (b4 != BAND) P4_base = satdata_other.obs_C(gobs4);
		if (b5 != BAND) P5_base = satdata_other.obs_C(gobs5);

		if (_observ == IONO_FREE && double_eq(PIF_base, 0.0)) return -1;
		if (_observ == RAW_SINGLE && double_eq(P1_base, 0.0)) return -1;
		if ((_observ == RAW_DOUBLE || _observ == RAW_ALL) && (double_eq(P1_base, 0.0) || double_eq(P2_base, 0.0))) return -1;
	}

	// 根据流动站的数据计算权重
	double weight_coef_b1 = _weightObs(satdata, gobs1);
	double weight_coef_b2 = _weightObs(satdata, gobs2);
	double weight_coef_b3 = _weightObs(satdata, gobs3);
	double weight_coef_b4 = _weightObs(satdata, gobs4);
	double weight_coef_b5 = _weightObs(satdata, gobs5);
	double weight_coef = 0;

	if (_observ == IONO_FREE)
	{
		double koef1 = 0.0;
		double koef2 = 0.0;
		satdata.coef_ionofree(b1, koef1, b2, koef2);
		weight_coef = (weight_coef_b1 + weight_coef_b2) / 2;
	}
	else weight_coef = weight_coef_b1;

	if (_isBase)
	{
		double weight_coef_b1_base = _weightObs(satdata_other, gobs1);
		double weight_coef_b2_base = _weightObs(satdata_other, gobs2);
		double weight_coef_b3_base = _weightObs(satdata_other, gobs3);
		double weight_coef_b4_base = _weightObs(satdata_other, gobs4);
		double weight_coef_b5_base = _weightObs(satdata_other, gobs5);
		double weight_coef_base = 0;

		if (_observ == IONO_FREE)
		{
			double koef1 = 0.0;
			double koef2 = 0.0;
			satdata_other.coef_ionofree(b1, koef1, b2, koef2);
			weight_coef_base = (weight_coef_b1_base + weight_coef_b2_base) / 2;
		}
		else weight_coef_base = weight_coef_b1_base;

		weight_coef_b1 = (double_eq(weight_coef_b1*weight_coef_b1_base, 0.0) ? 0.0 : (1.0 / (1.0 / weight_coef_b1 + 1.0 / weight_coef_b1_base)));
		weight_coef_b2 = (double_eq(weight_coef_b2*weight_coef_b2_base, 0.0) ? 0.0 : (1.0 / (1.0 / weight_coef_b2 + 1.0 / weight_coef_b2_base)));
		weight_coef_b3 = (double_eq(weight_coef_b3*weight_coef_b3_base, 0.0) ? 0.0 : (1.0 / (1.0 / weight_coef_b3 + 1.0 / weight_coef_b3_base)));
		weight_coef_b4 = (double_eq(weight_coef_b4*weight_coef_b4_base, 0.0) ? 0.0 : (1.0 / (1.0 / weight_coef_b4 + 1.0 / weight_coef_b4_base)));
		weight_coef_b5 = (double_eq(weight_coef_b5*weight_coef_b5_base, 0.0) ? 0.0 : (1.0 / (1.0 / weight_coef_b5 + 1.0 / weight_coef_b5_base)));
		weight_coef = (double_eq(weight_coef   *weight_coef_base, 0.0) ? 0.0 : (1.0 / (1.0 / weight_coef + 1.0 / weight_coef_base)));

	}

	int addIobs = 1;
	if ((_isBase && P3 > 0 && P3_base > 0) || ((!_isBase) && P3 > 0)) addIobs++;
	if ((_isBase && P4 > 0 && P4_base > 0) || ((!_isBase) && P4 > 0)) addIobs++;
	if ((_isBase && P5 > 0 && P5_base > 0) || ((!_isBase) && P5 > 0)) addIobs++;

	int usedFrq = addIobs + 1;
	if (_observ == IONO_FREE || _observ == RAW_SINGLE) _frqNum[satdata.sat()] = 1;
	if (_observ == RAW_DOUBLE) _frqNum[satdata.sat()] = 2;
	if (_observ == RAW_ALL)    _frqNum[satdata.sat()] = usedFrq;

	double sigCode = 0.0;
	switch (gs) {
	case GPS: sigCode = _sigCodeGPS; break;
	case GLO: sigCode = _sigCodeGLO; break;
	case GAL: sigCode = _sigCodeGAL; break;
	case BDS: sigCode = _sigCodeBDS; break;
	case QZS: sigCode = _sigCodeQZS; break;
	default: sigCode = 0.0;
	}

	// Create weight matrix
	P(iobs, iobs) = (weight_coef * (1 / (sigCode * sigCode)));

	if (_observ == RAW_DOUBLE ||
		_observ == RAW_ALL)
	{
		P(iobs + 1, iobs + 1) = (weight_coef_b2 * (1 / (sigCode * sigCode)));
	}
	if (_observ == RAW_ALL)
	{
		for (int i = 2; i <= addIobs; i++)
		{
			if (i == 2) P(iobs + i, iobs + i) = (weight_coef_b3 * (1 / (sigCode * sigCode)));
			if (i == 3) P(iobs + i, iobs + i) = (weight_coef_b4 * (1 / (sigCode * sigCode)));
			if (i == 4) P(iobs + i, iobs + i) = (weight_coef_b5 * (1 / (sigCode * sigCode)));
		}
	}

	bool com = true;
	if (_gnav) com = _gnav->com();


	// 组成流动站和基站的参数——用于计算模型�?
	t_gallpar param_rover;
	t_gallpar param_base;

	for (int iPar = 0; iPar < _param.parNumber(); iPar++)
	{
		if (_param[iPar].parType == t_gpar::CRD_X)
		{
			t_gpar par_rover = _param[iPar];
			param_rover.addParam(par_rover);
			if (_isBase)
			{
				t_gpar par_base = _param[iPar];
				par_base.site = _site_base;
				par_base.value(xyz_base[0]);
				param_base.addParam(par_base);
			}
		}
		else if (_param[iPar].parType == t_gpar::CRD_Y)
		{
			t_gpar par_rover = _param[iPar];
			param_rover.addParam(par_rover);
			if (_isBase)
			{
				t_gpar par_base = _param[iPar];
				par_base.site = _site_base;
				par_base.value(xyz_base[1]);
				param_base.addParam(par_base);
			}
		}
		else if (_param[iPar].parType == t_gpar::CRD_Z)
		{
			t_gpar par_rover = _param[iPar];
			param_rover.addParam(par_rover);
			if (_isBase)
			{
				t_gpar par_base = _param[iPar];
				par_base.site = _site_base;
				par_base.value(xyz_base[2]);
				param_base.addParam(par_base);
			}
		}
		else if (_param[iPar].parType == t_gpar::TRP)
		{
			if (_param[iPar].site == _site)
				param_rover.addParam(_param[iPar]);
			else if (_param[iPar].site == _site_base)
				param_base.addParam(_param[iPar]);
		}
		else if (_param[iPar].parType == t_gpar::VION
			|| _param[iPar].parType == t_gpar::SION)
		{
			// do nothing
		}
		else
		{
			// 考虑PPP的情�?
			if (_param[iPar].site == _site)
				param_rover.addParam(_param[iPar]);
		}
	}

	//Create reduced measurements (prefit residuals)
	double modObsP = 0.0;
	double modObsP_base = 0.0;

	if (_observ == IONO_FREE)
	{
		modObsP = _gModel->cmpObs(_epoch, param_rover, satdata, gobs1, com); if (modObsP < 0) return -1;
		_applyDCB(satdata, PIF, &gobs1, &gobs2);

		if (_isBase) {
			modObsP_base = _gModel_base->cmpObs(_epoch, param_base, satdata_other, gobs1, com); if (modObsP_base < 0) return -1;
			_applyDCB(satdata_other, PIF_base, &gobs1, &gobs2);
		}

		l(iobs) = PIF - modObsP;
		if (_isBase) l(iobs) = l(iobs) - (PIF_base - modObsP_base);

		_obs_index.push_back(make_pair(satdata.sat(), make_pair(FREQ_1, TYPE_C)));
	}
	else if (_observ == RAW_SINGLE)
	{
		modObsP = _gModel->cmpObs(_epoch, param_rover, satdata, gobs1, com); if (modObsP < 0) return -1;
		_applyDCB(satdata, P1, &gobs1);

		if (_isBase) {
			modObsP_base = _gModel_base->cmpObs(_epoch, param_base, satdata_other, gobs1, com); if (modObsP_base < 0) return -1;
			_applyDCB(satdata_other, P1_base, &gobs1);
		}

		l(iobs) = P1 - modObsP;
		if (_isBase) l(iobs) = l(iobs) - (P1_base - modObsP_base);

		_obs_index.push_back(make_pair(satdata.sat(), make_pair(FREQ_1, TYPE_C)));
	}
	else if (_observ == RAW_DOUBLE)
	{
		modObsP = _gModel->cmpObs(_epoch, param_rover, satdata, gobs1, com); if (modObsP < 0) return -1;
		_applyDCB(satdata, P1, &gobs1);

		if (_isBase) {
			modObsP_base = _gModel_base->cmpObs(_epoch, param_base, satdata_other, gobs1, com); if (modObsP_base < 0) return -1;
			_applyDCB(satdata_other, P1_base, &gobs1);
		}

		l(iobs) = P1 - modObsP;
		if (_isBase) l(iobs) = l(iobs) - (P1_base - modObsP_base);
		_obs_index.push_back(make_pair(satdata.sat(), make_pair(FREQ_1, TYPE_C)));


		modObsP = _gModel->cmpObs(_epoch, param_rover, satdata, gobs2, com); if (modObsP < 0) return -1;
		_applyDCB(satdata, P2, &gobs2);

		if (_isBase) {
			modObsP_base = _gModel_base->cmpObs(_epoch, param_base, satdata_other, gobs2, com); if (modObsP_base < 0) return -1;
			_applyDCB(satdata_other, P2_base, &gobs2);
		}

		l(iobs + 1) = P2 - modObsP;
		if (_isBase) l(iobs + 1) = l(iobs + 1) - (P2_base - modObsP_base);
		_obs_index.push_back(make_pair(satdata.sat(), make_pair(FREQ_2, TYPE_C)));

	}
	else if (_observ == RAW_ALL) {
		for (int i = 0; i <= addIobs; i++) {
			double Pi = 0.0;
			double Pi_base = 0.0;
			double modObsP_base = 0.0;
			FREQ_SEQ fi;
			if (i == 0) { modObsP = _gModel->cmpObs(_epoch, param_rover, satdata, gobs1, com); _applyDCB(satdata, P1, &gobs1); Pi = P1; fi = FREQ_1; }
			if (i == 1) { modObsP = _gModel->cmpObs(_epoch, param_rover, satdata, gobs2, com); _applyDCB(satdata, P2, &gobs2); Pi = P2; fi = FREQ_2; }
			if (i == 2) { modObsP = _gModel->cmpObs(_epoch, param_rover, satdata, gobs3, com);                                 Pi = P3; fi = FREQ_3; }
			if (i == 3) { modObsP = _gModel->cmpObs(_epoch, param_rover, satdata, gobs4, com);                                 Pi = P4; fi = FREQ_4; }
			if (i == 4) { modObsP = _gModel->cmpObs(_epoch, param_rover, satdata, gobs5, com);                                 Pi = P5; fi = FREQ_5; }

			if (_isBase)
			{
				if (i == 0) { modObsP_base = _gModel_base->cmpObs(_epoch, param_base, satdata_other, gobs1, com); _applyDCB(satdata_other, P1_base, &gobs1); Pi_base = P1_base; }
				if (i == 1) { modObsP_base = _gModel_base->cmpObs(_epoch, param_base, satdata_other, gobs2, com); _applyDCB(satdata_other, P2_base, &gobs2); Pi_base = P2_base; }
				if (i == 2) { modObsP_base = _gModel_base->cmpObs(_epoch, param_base, satdata_other, gobs3, com);                                            Pi_base = P3_base; }
				if (i == 3) { modObsP_base = _gModel_base->cmpObs(_epoch, param_base, satdata_other, gobs4, com);                                            Pi_base = P4_base; }
				if (i == 4) { modObsP_base = _gModel_base->cmpObs(_epoch, param_base, satdata_other, gobs5, com);                                            Pi_base = P5_base; }
			}
			if ((!_isBase && (double_eq(modObsP, -1.0) || double_eq(Pi, 0.0)))
				|| (_isBase && (double_eq(modObsP, -1.0) || double_eq(Pi, 0.0) || double_eq(modObsP_base, -1.0) || double_eq(Pi_base, 0.0))))
				continue;  // double_eq(modObsP, -1.0) modified by zhshen
			l(iobs + i) = Pi - modObsP;
			if (_isBase) l(iobs + i) = l(iobs + i) - (Pi_base - modObsP_base);
			_obs_index.push_back(make_pair(satdata.sat(), make_pair(fi, TYPE_C)));
		}
	}

	// 处理电离�?
	int iParIono = _param.getParam(_site, t_gpar::SION, satdata.sat());
	if (iParIono >= 0)
	{
		for (int i = iobs - 1; i < _obs_index.size(); i++)
		{
			double f1 = G01_F;
			double fk = satdata.frequency(t_gsys::band_priority(gs, _obs_index[i].second.first));
			double alfa = 0.0;
			if (_obs_index[i].second.second == TYPE_L) alfa = -(f1*f1) / (fk*fk);
			if (_obs_index[i].second.second == TYPE_C)  alfa = (f1*f1) / (fk*fk);
			double iono_delay = alfa * _param[iParIono].value();
			l(i + 1) = l(i + 1) - iono_delay;
		}
	}

	// Create first design matrix
	if (gs == GLO ||
		gs == GPS ||
		gs == GAL ||
		gs == BDS ||
		gs == QZS) {
		for (unsigned int ipar = 1; ipar <= _param.parNumber(); ipar++) {
			if (_param[ipar - 1].site == _site)
			{
				A(iobs, ipar) = _param[ipar - 1].partial(satdata, _epoch, ell, gobs1);
				if (_observ == RAW_DOUBLE ||
					_observ == RAW_ALL)       A(iobs + 1, ipar) = _param[ipar - 1].partial(satdata, _epoch, ell, gobs2);
				if (_observ == RAW_ALL) {
					for (int i = 2; i <= addIobs; i++) {
						if (i == 2) A(iobs + i, ipar) = _param[ipar - 1].partial(satdata, _epoch, ell, gobs3);
						if (i == 3) A(iobs + i, ipar) = _param[ipar - 1].partial(satdata, _epoch, ell, gobs4);
						if (i == 4) A(iobs + i, ipar) = _param[ipar - 1].partial(satdata, _epoch, ell, gobs5);
					}
				}
			}
			else if (_param[ipar - 1].site == _site_base)
			{
				A(iobs, ipar) = -_param[ipar - 1].partial(satdata_other, _epoch, ell_base, gobs1);
				if (_observ == RAW_DOUBLE ||
					_observ == RAW_ALL)       A(iobs + 1, ipar) = -_param[ipar - 1].partial(satdata_other, _epoch, ell_base, gobs2);
				if (_observ == RAW_ALL) {
					for (int i = 2; i <= addIobs; i++) {
						if (i == 2) A(iobs + i, ipar) = -_param[ipar - 1].partial(satdata_other, _epoch, ell_base, gobs3);
						if (i == 3) A(iobs + i, ipar) = -_param[ipar - 1].partial(satdata_other, _epoch, ell_base, gobs4);
						if (i == 4) A(iobs + i, ipar) = -_param[ipar - 1].partial(satdata_other, _epoch, ell_base, gobs5);
					}
				}
			}
		}
	}
	else return -1;

#ifdef DEBUG       
	cout << satdata.sat() << fixed << setprecision(10)
		<< "  PIF: (" << BAND_1 << "," << BAND_2 << ") :" << setw(12) << satdata.P3(b1, b2)
		<< "  cmpObs: " << setw(12) << modObsP
		<< endl;
#endif   

	if (_observ == RAW_DOUBLE) iobs += 2;
	else if (_observ == RAW_ALL) iobs = iobs + addIobs + 1;
	else iobs++;

	return 1;
}

int gsins::t_gpvtflt_old::_addObsL(t_gsatdata & satdata, unsigned int & iobs, Matrix & A, ColumnVector & l, SymmetricMatrix & P)
{
	gtrace("t_gpvtflt_old::_addObsP");

	t_gsatdata satdata_other;

	if (_isBase)
	{
		for (int i = 0; i < _data_base.size(); i++)
		{
			if (_data_base[i].sat() == satdata.sat())
				satdata_other = _data_base[i];
		}
		if (satdata_other.site() == "") return -1;
	}

	// 获取流动站和基站的坐�?
	t_gtriple xyz, ell;
	if (_param.getCrdParam(_site, xyz) < 0)
	{
		xyz = _gallobj->obj(_site)->crd_arp(_epoch);
	}
	xyz2ell(xyz.crd_array(), ell.crd_array(), false);

	t_gtriple xyz_base, ell_base;
	if (_isBase)
	{
		if (_param.getCrdParam(_site_base, xyz_base) < 0)
		{
			xyz_base = _gallobj->obj(_site_base)->crd_arp(_epoch);
		}
		xyz2ell(xyz_base.crd_array(), ell_base.crd_array(), false);
	}

	// 根据流动站数据选择频率和观测值类�?
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
	else {
		// select fix defined band according the table
		vector < GOBSBAND > band = dynamic_cast<t_gsetgnss*>(_set)->band(gs);
		if (band.size())
		{
			b1 = band[0];
			b2 = band[1];
			if (band.size() > 2) b3 = band[2];
			if (band.size() > 3) b4 = band[3];
			if (band.size() > 4) b5 = band[4];
		}
		else {
			b1 = t_gsys::band_priority(gs, FREQ_1);
			b2 = t_gsys::band_priority(gs, FREQ_2);
			b3 = t_gsys::band_priority(gs, FREQ_3);
			b4 = t_gsys::band_priority(gs, FREQ_4);
			b5 = t_gsys::band_priority(gs, FREQ_5);
		}
	}

	t_gobs gobs1, gobs2, gobs3, gobs4, gobs5;

	_getgobs(satdata.sat(), TYPE_L, b1, gobs1);
	_getgobs(satdata.sat(), TYPE_L, b2, gobs2);
	_getgobs(satdata.sat(), TYPE_L, b3, gobs3);
	_getgobs(satdata.sat(), TYPE_L, b4, gobs4);
	_getgobs(satdata.sat(), TYPE_L, b5, gobs5);

	if (_phase) {

		bool com = true;
		if (_gnav) com = _gnav->com();
		// 收集流动站的数据
		double LIF, L1, L2, L3, L4, L5;
		LIF = L1 = L2 = L3 = L4 = L5 = 0.0;

		LIF = satdata.L3(gobs1, gobs2);
		L1 = satdata.obs_L(gobs1);
		L2 = satdata.obs_L(gobs2);
		if (b3 != BAND) L3 = satdata.obs_L(gobs3);
		if (b4 != BAND) L4 = satdata.obs_L(gobs4);
		if (b5 != BAND) L5 = satdata.obs_L(gobs5);

		if (_observ == IONO_FREE && double_eq(LIF, 0.0)) return -1;
		if (_observ == RAW_SINGLE && double_eq(L1, 0.0)) return -1;
		if ((_observ == RAW_DOUBLE || _observ == RAW_ALL) && (double_eq(L1, 0.0) || double_eq(L2, 0.0))) return -1;

		// 收集基站的数�?
		double LIF_base, L1_base, L2_base, L3_base, L4_base, L5_base;
		LIF_base = L1_base = L2_base = L3_base = L4_base = L5_base = 0.0;

		if (_isBase)
		{
			LIF_base = satdata_other.L3(gobs1, gobs2);
			L1_base = satdata_other.obs_L(gobs1);
			L2_base = satdata_other.obs_L(gobs2);
			if (b3 != BAND) L3_base = satdata_other.obs_L(gobs3);
			if (b4 != BAND) L4_base = satdata_other.obs_L(gobs4);
			if (b5 != BAND) L5_base = satdata_other.obs_L(gobs5);

			if (_observ == IONO_FREE && double_eq(LIF_base, 0.0)) return -1;
			if (_observ == RAW_SINGLE && double_eq(L1_base, 0.0)) return -1;
			if ((_observ == RAW_DOUBLE || _observ == RAW_ALL) && (double_eq(L1_base, 0.0) || double_eq(L2_base, 0.0))) return -1;
		}

		// 根据流动站数据计算权�?
		double weight_coef = _weightObs(satdata, gobs1);

		if (_isBase)
		{
			double weight_coef_base = _weightObs(satdata_other, gobs1);
			weight_coef = (double_eq(weight_coef*weight_coef_base, 0.0) ? 0.0 : (1.0 / (1.0 / weight_coef + 1.0 / weight_coef_base)));
		}

		int addIobs = 1;
		if (L3 > 0) addIobs++;
		if (L4 > 0) addIobs++;
		if (L5 > 0) addIobs++;

		double sigPhase = 0.0;
		switch (gs) {
		case GPS: sigPhase = _sigPhaseGPS; break;
		case GLO: sigPhase = _sigPhaseGLO; break;
		case GAL: sigPhase = _sigPhaseGAL; break;
		case BDS: sigPhase = _sigPhaseBDS; break;
		case QZS: sigPhase = _sigPhaseQZS; break;
		default: sigPhase = 0.0;
		}

		// Create weight matrix
		P(iobs, iobs) = weight_coef * (1 / (sigPhase * sigPhase));
		if (_observ == RAW_DOUBLE ||
			_observ == RAW_ALL)      P(iobs + 1, iobs + 1) = weight_coef * (1 / (sigPhase * sigPhase));
		if (_observ == RAW_ALL)
			for (int i = 2; i <= addIobs; i++) P(iobs + i, iobs + i) = weight_coef * (1 / (sigPhase * sigPhase));


		GOBS type1, type2, type3, type4, type5;

		if (gobs1.gattr().attr() == ATTR) type1 = satdata.id_phase(gobs1.gband().band());
		else type1 = gobs1.gobs();
		if (gobs2.gattr().attr() == ATTR) type2 = satdata.id_phase(gobs2.gband().band());
		else type2 = gobs2.gobs();
		if (gobs3.gattr().attr() == ATTR) type3 = satdata.id_phase(gobs3.gband().band());
		else type3 = gobs3.gobs();
		if (gobs4.gattr().attr() == ATTR) type4 = satdata.id_phase(gobs4.gband().band());
		else type4 = gobs4.gobs();
		if (gobs5.gattr().attr() == ATTR) type5 = satdata.id_phase(gobs5.gband().band());
		else type5 = gobs5.gobs();

		// 组成流动站和基站的参数——用于计算模型�?
		t_gallpar param_rover;
		t_gallpar param_base;

		for (int iPar = 0; iPar < _param.parNumber(); iPar++)
		{
			if (_param[iPar].parType == t_gpar::CRD_X)
			{
				t_gpar par_rover = _param[iPar];
				param_rover.addParam(par_rover);
				if (_isBase)
				{
					t_gpar par_base = _param[iPar];
					par_base.site = _site_base;
					par_base.value(xyz_base[0]);
					param_base.addParam(par_base);
				}
			}
			else if (_param[iPar].parType == t_gpar::CRD_Y)
			{
				t_gpar par_rover = _param[iPar];
				param_rover.addParam(par_rover);
				if (_isBase)
				{
					t_gpar par_base = _param[iPar];
					par_base.site = _site_base;
					par_base.value(xyz_base[1]);
					param_base.addParam(par_base);
				}
			}
			else if (_param[iPar].parType == t_gpar::CRD_Z)
			{
				t_gpar par_rover = _param[iPar];
				param_rover.addParam(par_rover);
				if (_isBase)
				{
					t_gpar par_base = _param[iPar];
					par_base.site = _site_base;
					par_base.value(xyz_base[2]);
					param_base.addParam(par_base);
				}
			}
			else if (_param[iPar].parType == t_gpar::TRP)
			{
				if (_param[iPar].site == _site)
					param_rover.addParam(_param[iPar]);
				else if (_param[iPar].site == _site_base)
					param_base.addParam(_param[iPar]);
			}
			else if (_param[iPar].parType == t_gpar::VION
				|| _param[iPar].parType == t_gpar::SION)
			{
				// do nothing
			}
			else
			{
				// 考虑PPP的情�?
				if (_param[iPar].site == _site)
					param_rover.addParam(_param[iPar]);
			}
		}

		for (size_t iPar = 0; iPar < _param.parNumber(); iPar++) {
			if (_param[iPar].prn.compare(satdata.sat()) == 0) {

				// 根据参数中有的模糊度来选择对应的数�?
				// 计算流动站的模型�?
				double modObsL = 0.0;
				double obs_L = 0.0;
				if (_param[iPar].parType == t_gpar::AMB_IF) { obs_L = LIF; modObsL = _gModel->cmpObs(_epoch, param_rover, satdata, gobs1, com); }
				else if (_param[iPar].parType == t_gpar::AMB_L1) { obs_L = L1;  modObsL = _gModel->cmpObs(_epoch, param_rover, satdata, gobs1, com); }
				else if (_param[iPar].parType == t_gpar::AMB_L2) { obs_L = L2;  modObsL = _gModel->cmpObs(_epoch, param_rover, satdata, gobs2, com); }
				else if (_param[iPar].parType == t_gpar::AMB_L3) { obs_L = L3;  modObsL = _gModel->cmpObs(_epoch, param_rover, satdata, gobs3, com); }
				else if (_param[iPar].parType == t_gpar::AMB_L4) { obs_L = L4;  modObsL = _gModel->cmpObs(_epoch, param_rover, satdata, gobs4, com); }
				else if (_param[iPar].parType == t_gpar::AMB_L5) { obs_L = L5;  modObsL = _gModel->cmpObs(_epoch, param_rover, satdata, gobs5, com); }

				// 计算基站的模型�?
				double modObsL_base = 0.0;
				double obs_L_base = 0.0;
				if (_isBase)
				{
					if (_param[iPar].parType == t_gpar::AMB_IF) { obs_L_base = LIF_base; modObsL_base = _gModel_base->cmpObs(_epoch, param_base, satdata_other, gobs1, com); }
					else if (_param[iPar].parType == t_gpar::AMB_L1) { obs_L_base = L1_base;  modObsL_base = _gModel_base->cmpObs(_epoch, param_base, satdata_other, gobs1, com); }
					else if (_param[iPar].parType == t_gpar::AMB_L2) { obs_L_base = L2_base;  modObsL_base = _gModel_base->cmpObs(_epoch, param_base, satdata_other, gobs2, com); }
					else if (_param[iPar].parType == t_gpar::AMB_L3) { obs_L_base = L3_base;  modObsL_base = _gModel_base->cmpObs(_epoch, param_base, satdata_other, gobs3, com); }
					else if (_param[iPar].parType == t_gpar::AMB_L4) { obs_L_base = L4_base;  modObsL_base = _gModel_base->cmpObs(_epoch, param_base, satdata_other, gobs4, com); }
					else if (_param[iPar].parType == t_gpar::AMB_L5) { obs_L_base = L5_base;  modObsL_base = _gModel_base->cmpObs(_epoch, param_base, satdata_other, gobs5, com); }

				}

				if ((!_isBase &&modObsL < 0) || (_isBase && (modObsL < 0 || modObsL_base < 0))) return -1;

				// Reset ambiguity in case of cycle slip
				if (!_isBase)
				{
					satdata.addslip(false);
					if ((_param[iPar].parType == t_gpar::AMB_L1 && satdata.getlli(type1) >= 1) ||
						(_param[iPar].parType == t_gpar::AMB_L2 && satdata.getlli(type2) >= 1) ||
						(_param[iPar].parType == t_gpar::AMB_L3 && satdata.getlli(type3) >= 1) ||
						(_param[iPar].parType == t_gpar::AMB_L4 && satdata.getlli(type4) >= 1) ||
						(_param[iPar].parType == t_gpar::AMB_L5 && satdata.getlli(type5) >= 1) ||
						(_param[iPar].parType == t_gpar::AMB_IF && (satdata.getlli(type1) >= 1 || satdata.getlli(type2) >= 1))) {
						satdata.addslip(true);
						string prn = satdata.sat();
						if (_param[iPar].parType == t_gpar::AMB_IF) { _param[iPar].value(obs_L - satdata.L3(b1, b2)); }
						else if (_param[iPar].parType == t_gpar::AMB_L1) { _param[iPar].value(obs_L - satdata.obs_C(t_gband(b1, ATTR))); }
						else if (_param[iPar].parType == t_gpar::AMB_L2) { _param[iPar].value(obs_L - satdata.obs_C(t_gband(b2, ATTR))); }
						else if (_param[iPar].parType == t_gpar::AMB_L3) { _param[iPar].value(obs_L - satdata.obs_C(t_gband(b3, ATTR))); }
						else if (_param[iPar].parType == t_gpar::AMB_L4) { _param[iPar].value(obs_L - satdata.obs_C(t_gband(b4, ATTR))); }
						else if (_param[iPar].parType == t_gpar::AMB_L5) { _param[iPar].value(obs_L - satdata.obs_C(t_gband(b5, ATTR))); }
						//cout << "cycle slip !!!!" << endl;
						_Qx(iPar + 1, iPar + 1) = _sigAmbig * _sigAmbig;
						_newAMB[prn] = 1;
						if (_smooth) _slips.insert(prn);
						//              cout << "Cycle slip: " << prn << " " << gobs2str(type1) << " " << gobs2str(type2) << " " << _epoch.str_hms() << fixed << setprecision(3) << " " << _param[iPar].value() << " " << LIF << " " << modObsL << endl;
						if (_smooth) {
							SymmetricMatrix Qp = _fltdata.Qp();
							Qp(iPar + 1, iPar + 1) = _sigAmbig * _sigAmbig;
							_fltdata.Qp(Qp);
						}
					}
				}
				else
				{
					satdata.addslip(false);
					satdata_other.addslip(false);
					if ((_param[iPar].parType == t_gpar::AMB_L1 && (satdata.getlli(type1) >= 1 || satdata_other.getlli(type1) >= 1)) ||
						(_param[iPar].parType == t_gpar::AMB_L2 && (satdata.getlli(type2) >= 1 || satdata_other.getlli(type2) >= 1)) ||
						(_param[iPar].parType == t_gpar::AMB_L3 && (satdata.getlli(type3) >= 1 || satdata_other.getlli(type3) >= 1)) ||
						(_param[iPar].parType == t_gpar::AMB_L4 && (satdata.getlli(type4) >= 1 || satdata_other.getlli(type4) >= 1)) ||
						(_param[iPar].parType == t_gpar::AMB_L5 && (satdata.getlli(type5) >= 1 || satdata_other.getlli(type5) >= 1)) ||
						(_param[iPar].parType == t_gpar::AMB_IF && (satdata.getlli(type1) >= 1 || satdata.getlli(type2) >= 1 || satdata_other.getlli(type1) >= 1 || satdata_other.getlli(type2) >= 1)))
					{
						satdata.addslip(true);
						satdata_other.addslip(true);
						string prn = satdata.sat();
						if (_param[iPar].parType == t_gpar::AMB_IF) { _param[iPar].value(obs_L - satdata.L3(b1, b2) - (obs_L_base - satdata_other.L3(b1, b2))); }
						else if (_param[iPar].parType == t_gpar::AMB_L1) { _param[iPar].value(obs_L - satdata.obs_C(t_gband(b1, ATTR)) - (obs_L_base - satdata_other.obs_C(t_gband(b1, ATTR)))); }
						else if (_param[iPar].parType == t_gpar::AMB_L2) { _param[iPar].value(obs_L - satdata.obs_C(t_gband(b2, ATTR)) - (obs_L_base - satdata_other.obs_C(t_gband(b2, ATTR)))); }
						else if (_param[iPar].parType == t_gpar::AMB_L3) { _param[iPar].value(obs_L - satdata.obs_C(t_gband(b3, ATTR)) - (obs_L_base - satdata_other.obs_C(t_gband(b3, ATTR)))); }
						else if (_param[iPar].parType == t_gpar::AMB_L4) { _param[iPar].value(obs_L - satdata.obs_C(t_gband(b4, ATTR)) - (obs_L_base - satdata_other.obs_C(t_gband(b4, ATTR)))); }
						else if (_param[iPar].parType == t_gpar::AMB_L5) { _param[iPar].value(obs_L - satdata.obs_C(t_gband(b5, ATTR)) - (obs_L_base - satdata_other.obs_C(t_gband(b5, ATTR)))); }						//cout << "cycle slip !!!!" << endl;
						_Qx(iPar + 1, iPar + 1) = _sigAmbig * _sigAmbig;
						_newAMB[prn] = 1;
						if (_smooth) _slips.insert(prn);
						//              cout << "Cycle slip: " << prn << " " << gobs2str(type1) << " " << gobs2str(type2) << " " << _epoch.str_hms() << fixed << setprecision(3) << " " << _param[iPar].value() << " " << LIF << " " << modObsL << endl;
						if (_smooth) {
							SymmetricMatrix Qp = _fltdata.Qp();
							Qp(iPar + 1, iPar + 1) = _sigAmbig * _sigAmbig;
							_fltdata.Qp(Qp);
						}
					}
				}

				// Create reduced measurement (prefit residuals)                        
				if (_param[iPar].parType == t_gpar::AMB_L1)
				{
					l(iobs) = obs_L - _param[iPar].value() - modObsL;
					_obs_index.push_back(make_pair(satdata.sat(), make_pair(FREQ_1, TYPE_L)));
				}
				else if (_param[iPar].parType == t_gpar::AMB_IF)
				{
					l(iobs) = obs_L - _param[iPar].value() - modObsL;
					_obs_index.push_back(make_pair(satdata.sat(), make_pair(FREQ_1, TYPE_L)));
				}
				if ((_observ == RAW_DOUBLE || _observ == RAW_ALL) && _param[iPar].parType == t_gpar::AMB_L2) {
					l(iobs + 1) = obs_L - _param[iPar].value() - modObsL;
					_obs_index.push_back(make_pair(satdata.sat(), make_pair(FREQ_2, TYPE_L)));
				}
				if (_observ == RAW_ALL) {
					if (_param[iPar].parType == t_gpar::AMB_L3) { l(iobs + 2) = obs_L - _param[iPar].value() - modObsL, _obs_index.push_back(make_pair(satdata.sat(), make_pair(FREQ_3, TYPE_L))); }
					if (_param[iPar].parType == t_gpar::AMB_L4) { l(iobs + 2) = obs_L - _param[iPar].value() - modObsL, _obs_index.push_back(make_pair(satdata.sat(), make_pair(FREQ_4, TYPE_L))); }
					if (_param[iPar].parType == t_gpar::AMB_L5) { l(iobs + 2) = obs_L - _param[iPar].value() - modObsL, _obs_index.push_back(make_pair(satdata.sat(), make_pair(FREQ_5, TYPE_L))); }
				}

				if (_isBase)
				{
					if (_param[iPar].parType == t_gpar::AMB_L1 || _param[iPar].parType == t_gpar::AMB_IF) {
						l(iobs) = l(iobs) - (obs_L_base - modObsL_base);
					}
					if ((_observ == RAW_DOUBLE || _observ == RAW_ALL) && _param[iPar].parType == t_gpar::AMB_L2) {
						l(iobs + 1) = l(iobs + 1) - (obs_L_base - modObsL_base);
					}
					if (_observ == RAW_ALL) {
						if (_param[iPar].parType == t_gpar::AMB_L3) l(iobs + 2) = l(iobs + 2) - (obs_L_base - modObsL_base);
						if (_param[iPar].parType == t_gpar::AMB_L4) l(iobs + 3) = l(iobs + 3) - (obs_L_base - modObsL_base);
						if (_param[iPar].parType == t_gpar::AMB_L5) l(iobs + 4) = l(iobs + 4) - (obs_L_base - modObsL_base);
					}
				}

				// check satellite post shadow period
				if (satdata.ecl()) return -1;
			}
		}

		// 处理电离�?
		int iParIono = _param.getParam(_site, t_gpar::SION, satdata.sat());
		if (iParIono >= 0)
		{
			for (int i = iobs - 1; i < _obs_index.size(); i++)
			{
				double f1 = G01_F;
				double fk = satdata.frequency(t_gsys::band_priority(gs, _obs_index[i].second.first));
				double alfa = 0.0;
				if (_obs_index[i].second.second == TYPE_L) alfa = -(f1*f1) / (fk*fk);
				if (_obs_index[i].second.second == TYPE_C)  alfa = (f1*f1) / (fk*fk);
				double iono_delay = alfa * _param[iParIono].value();
				l(i + 1) = l(i + 1) - iono_delay;
			}
		}

		// Create first design matrix 
		if (gs == GLO ||
			gs == GPS ||
			gs == GAL ||
			gs == BDS ||
			gs == QZS) {
			for (unsigned int ipar = 1; ipar <= _param.parNumber(); ipar++) {
				if (_param[ipar - 1].site == _site)
				{
					A(iobs, ipar) = _param[ipar - 1].partial(satdata, _epoch, ell, gobs1);
					if (_observ == RAW_DOUBLE ||
						_observ == RAW_ALL)       A(iobs + 1, ipar) = _param[ipar - 1].partial(satdata, _epoch, ell, gobs2);
					if (_observ == RAW_ALL) {
						for (int i = 2; i <= addIobs; i++) {
							if (i == 2) A(iobs + i, ipar) = _param[ipar - 1].partial(satdata, _epoch, ell, gobs3);
							if (i == 3) A(iobs + i, ipar) = _param[ipar - 1].partial(satdata, _epoch, ell, gobs4);
							if (i == 4) A(iobs + i, ipar) = _param[ipar - 1].partial(satdata, _epoch, ell, gobs5);
						}
					}
					if (_param[ipar - 1].parType == t_gpar::AMB_IF && satdata.sat() == _param[ipar - 1].prn && gobs1.is_phase())
						A(iobs, ipar) = 1.0;
					else if (_param[ipar - 1].parType == t_gpar::AMB_L1&& satdata.sat() == _param[ipar - 1].prn && gobs1.is_phase())
						A(iobs, ipar) = 1.0;
					else if (_param[ipar - 1].parType == t_gpar::AMB_L2&& satdata.sat() == _param[ipar - 1].prn && gobs2.is_phase())
						A(iobs + 1, ipar) = 1.0;
					else if (_param[ipar - 1].parType == t_gpar::AMB_L3&& satdata.sat() == _param[ipar - 1].prn && gobs3.is_phase())
						A(iobs + 2, ipar) = 1.0;
					else if (_param[ipar - 1].parType == t_gpar::AMB_L4&& satdata.sat() == _param[ipar - 1].prn && gobs4.is_phase())
						A(iobs + 3, ipar) = 1.0;
					else if (_param[ipar - 1].parType == t_gpar::AMB_L5&& satdata.sat() == _param[ipar - 1].prn && gobs5.is_phase())
						A(iobs + 4, ipar) = 1.0;
				}
				else if (_param[ipar - 1].site == _site_base)
				{
					A(iobs, ipar) = -_param[ipar - 1].partial(satdata_other, _epoch, ell_base, gobs1);
					if (_observ == RAW_DOUBLE ||
						_observ == RAW_ALL)       A(iobs + 1, ipar) = -_param[ipar - 1].partial(satdata_other, _epoch, ell_base, gobs2);
					if (_observ == RAW_ALL) {
						for (int i = 2; i <= addIobs; i++) {
							if (i == 2) A(iobs + i, ipar) = -_param[ipar - 1].partial(satdata_other, _epoch, ell_base, gobs3);
							if (i == 3) A(iobs + i, ipar) = -_param[ipar - 1].partial(satdata_other, _epoch, ell_base, gobs4);
							if (i == 4) A(iobs + i, ipar) = -_param[ipar - 1].partial(satdata_other, _epoch, ell_base, gobs5);
						}
					}
				}
			}
		}
		else return -1;
		if (_observ == RAW_DOUBLE) iobs += 2;
		else if (_observ == RAW_ALL) iobs = iobs + addIobs + 1;
		else iobs++;

	}

	return 1;
}

int gsins::t_gpvtflt_old::_combineMW(t_gsatdata& satdata)
{
	gtrace("t_gpvtflt::_combineMW");

	string sat = satdata.sat();
	double obs_intv = _sampling;
	t_gtime crt_time = _epoch;
	if (crt_time == _beg_time)
	{
		_MW[crt_time - obs_intv][sat][1] = 0.0;
		_MW[crt_time - obs_intv][sat][2] = 0.0;
		_MW[crt_time - obs_intv][sat][3] = 0.0;
		_MW[crt_time - obs_intv][sat][4] = 0.0;
		_MW[crt_time - obs_intv][sat][5] = 0.0; // ele
	}

	GOBSBAND b1, b2;
	if (_auto_band) {    // automatic dual band selection -> for Anubis purpose
		set<GOBSBAND> bands = satdata.band_avail();
		auto itBAND = bands.begin();
		b1 = *itBAND;
		itBAND++;
		b2 = *itBAND;
	}
	else {            // select fix defined band according the table 
		vector < GOBSBAND > band = dynamic_cast<t_gsetgnss*>(_set)->band(satdata.gsys());
		if (band.size())
		{
			b1 = band[0];
			b2 = band[1];
		}
		else {
			b1 = t_gsys::band_priority(satdata.gsys(), FREQ_1);
			b2 = t_gsys::band_priority(satdata.gsys(), FREQ_2);
		}
	}
	t_gobs gobsP1 = t_gobs(satdata.select_range(b1));
	t_gobs gobsP2 = t_gobs(satdata.select_range(b2));
	t_gobs gobsL1 = t_gobs(satdata.select_phase(b1));
	t_gobs gobsL2 = t_gobs(satdata.select_phase(b2));
	double mw_obs = satdata.MW_cycle(gobsL1, gobsL2, gobsP1, gobsP2);

	//dcb correct  �Ѿ��ڹ۲�֮���޸���dcb
	//OBSCOMBIN type = OBSCOMBIN::MW_COMBIN;
	//apply_DCB(_gallbias, satdata, mw_obs, type, &gobsP1, &gobsP2);
	if (_isBase)
	{
		t_gsatdata satdata_base;
		for (int i = 0; i < _data_base.size(); i++)
		{
			if (_data_base[i].sat() == satdata.sat())
				satdata_base = _data_base[i];
		}
		if (satdata_base.site() == "") mw_obs = 0.0;
		else {
			gobsP1 = t_gobs(satdata_base.select_range(b1));
			gobsP2 = t_gobs(satdata_base.select_range(b2));
			gobsL1 = t_gobs(satdata_base.select_phase(b1));
			gobsL2 = t_gobs(satdata_base.select_phase(b2));
			double mw_obs_base = satdata_base.MW_cycle(gobsL1, gobsL2, gobsP1, gobsP2);

			//apply_DCB(_gallbias, satdata_base, mw_obs_base, type, &gobsP1, &gobsP2);  �Ѿ��ڹ۲�֮���޸���dcb

			if (double_eq(mw_obs, 0.0) || double_eq(mw_obs_base, 0.0)) mw_obs = 0.0;
			else mw_obs -= mw_obs_base;
		}

	}

	if (!double_eq(mw_obs, 0.0))
	{
		// calculate smooth MW
		if (_MW[crt_time - obs_intv].find(sat) != _MW[crt_time - obs_intv].end())
		{
			_MW[crt_time][sat][1] = _MW[crt_time - obs_intv][sat][1] + 1;
			_MW[crt_time][sat][4] = mw_obs;
			if (!double_eq(_MW[crt_time][sat][1], 1.0) && _MW[crt_time][sat][1] > 1.0)
			{
				_MW[crt_time][sat][3] = _MW[crt_time - obs_intv][sat][3] * (_MW[crt_time][sat][1] - 2) /
					(_MW[crt_time][sat][1] - 1) + pow(mw_obs - _MW[crt_time - obs_intv][sat][2], 2)
					/ _MW[crt_time][sat][1];
			}
			else
			{
				_MW[crt_time][sat][3] = _MW[crt_time - obs_intv][sat][3];
			}
			_MW[crt_time][sat][2] = _MW[crt_time - obs_intv][sat][2] + (mw_obs - _MW[crt_time - obs_intv][sat][2])
				/ _MW[crt_time][sat][1];
			_MW[crt_time][sat][5] = satdata.ele_deg();
		}
		// exclude mw begining ,it is different from last all
		else
		{
			_MW[crt_time][sat][1] = 1;
			_MW[crt_time][sat][4] = mw_obs;
			_MW[crt_time][sat][3] = 0.0;
			_MW[crt_time][sat][2] = mw_obs / _MW[crt_time][sat][1];
			_MW[crt_time][sat][5] = satdata.ele_deg();
		}
	}
	// exclude zero wl
	else
	{
		_MW[crt_time][sat][1] = 0.0;
		_MW[crt_time][sat][2] = 0.0;
		_MW[crt_time][sat][3] = 0.0;
		_MW[crt_time][sat][4] = 0.0;
		_MW[crt_time][sat][5] = 0.0;
		return -1;
	}


	return 1;
}

bool gsins::t_gpvtflt_old::_external_pos(const t_gtriple& xyz_r, const t_gtriple& rms)
{
	_extn_pos = xyz_r;
	_extn_rms = rms;
	_pos_constrain = true;
	return true;
}

int gsins::t_gpvtflt_old::_processEpochVel()
{
	gtrace("t_gpvtflt::_processEpochVel");

	double threshold = 0.0;

	Matrix A;
	DiagonalMatrix P;
	ColumnVector l, dx;
	vector<t_gsatdata>::iterator it;
	t_gallpar param_vel;

	int ipar = 0;
	param_vel.addParam(t_gpar(_site, t_gpar::VEL_X, ++ipar, ""));
	param_vel.addParam(t_gpar(_site, t_gpar::VEL_Y, ++ipar, ""));
	param_vel.addParam(t_gpar(_site, t_gpar::VEL_Z, ++ipar, ""));
	param_vel.addParam(t_gpar(_site, t_gpar::CLK_RAT, ++ipar, ""));

	if (!_initialized || _Qx_vel(1, 1) == 0) {
		_Qx_vel = _sig_init_vel * _sig_init_vel;
	}
	else {
		for (int i = 1; i <= param_vel.parNumber(); i++) {
			if (param_vel[i - 1].parType == t_gpar::VEL_X)
			{
				param_vel[i - 1].value(_vel[0]);
				_Qx_vel(i, i) += _velStoModel->getQ() * _velStoModel->getQ();
			}
			else if (param_vel[i - 1].parType == t_gpar::VEL_Y)
			{
				param_vel[i - 1].value(_vel[1]);
				_Qx_vel(i, i) += _velStoModel->getQ() * _velStoModel->getQ();
			}
			else if (param_vel[i - 1].parType == t_gpar::VEL_Z)
			{
				param_vel[i - 1].value(_vel[2]);
				_Qx_vel(i, i) += _velStoModel->getQ() * _velStoModel->getQ();
			}
			else if (param_vel[i - 1].parType == t_gpar::CLK_RAT)
			{
				param_vel[i - 1].value(0);
				_Qx_vel(i, i) = _dclkStoModel->getQ() * _dclkStoModel->getQ();
			}
		}
	}


	do {
		// define a number of measurements 
		unsigned int nObs = _data.size();
		if (_doppler) { nObs *= 5; }

		unsigned int nPar = param_vel.parNumber();
		A.ReSize(nObs, nPar); A = 0.0;
		l.ReSize(nObs);       l = 0.0;
		P.ReSize(nObs);       P = 0.0;
		dx.ReSize(nPar);     dx = 0.0;
		unsigned int iobs = 1;

		t_gtriple groundXYZ;
		_param.getCrdParam(_site, groundXYZ);

		// Create matrices and vectors for estimation
		// ----------------------------------------------------------------------------------
		// loop over all measurement in current epoch

		for (it = _data.begin(); it != _data.end(); ) {
			if (_addObsD(*it, iobs, param_vel, groundXYZ, A, l, P) > 0) {}
			else { it = _data.erase(it); continue; } // just return next iterator

			++it; // increase iterator if not erased
		}
		// delete zero rows/cols
		iobs--;
		A = A.Rows(1, iobs);
		P = P.SymSubMatrix(1, iobs);
		l = l.Rows(1, iobs);

		_filter->update(A, P, l, dx, _Qx_vel);

		int freedom = A.Nrows() - A.Ncols();
		if (freedom < 1) {
			if (_log) _log->comment(1, "t_gpvtflt", "No redundant observations!");
			//    _restore(QsavBP, XsavBP);     
			//    return -1;
			freedom = 1;
		}

		int index = param_vel.getParam(_site, t_gpar::VEL_X, "");

		for (unsigned int iPar = 0; iPar < param_vel.parNumber(); iPar++) {
			param_vel[iPar].value(param_vel[iPar].value() + dx(param_vel[iPar].index));
		}
		threshold = sqrt(dx(index + 1) * dx(index + 1) + dx(index + 2) * dx(index + 2) + dx(index + 3) * dx(index + 3));
	} while (threshold > 1e-3);

	if (param_vel.getVelParam(_site, _vel) < 0)return -1;

	return 1;
}


int gsins::t_gpvtflt_old::_processEpoch_new(const t_gtime & runEpoch)
{
	gtrace("t_gpvtflt::_processEpoch");

	//int irc = t_gsppflt::_processEpoch(runEpoch);

#ifdef DEBUG
	cout << "gpvtflt: Epoch - full data size: " << _data.size() << " "
		<< runEpoch.str_hms() << endl;
#endif
	//cout << "gpvtflt: Epoch - full data size: " << _data.size() << " "
	//	<< runEpoch.str_hms() << endl;

	vector<t_gsatdata>::iterator ItException = _data.end();
	ColumnVector Threshold;
	double threshold;
	if (_grec == nullptr) {
		if (_log) _log->comment(1, "gpvtflt", "No receiver settings available!!!");
		return -1;
	}

	_epoch = runEpoch;
	t_gtriple crdapr = _grec->crd_arp(_epoch);

	if (double_eq(crdapr[0], 0.0) && double_eq(crdapr[1], 0.0) &&
		double_eq(crdapr[2], 0.0)) {
		_valid_crd_xml = false;
	}
	else {
		_valid_crd_xml = true;
	}

	if (!_valid_crd_xml) _sig_init_crd = 100.0;

	Matrix A;
	SymmetricMatrix P;
	ColumnVector l, dx;
	ColumnVector v_orig, v_norm;
	SymmetricMatrix Qsav, QsavBP;
	vector<t_gsatdata>::iterator it;
	t_gallpar XsavBP;
	double vtpv;
	int nobs_total, npar_number;
	string outlier = "";

	_cntrep = 0;    // number of iterations caused by outliers

	do {

		// erase sat satellite because of outlier
		it = _data.begin();
		while (it != _data.end()){
			if (it->sat() == outlier){
				it = _data.erase(it); continue;
			}
			it++;
		}

		if (_prepareData() < 0){
			_predict();//todo:add time gap
			return -1;
		}

		QsavBP = _Qx; XsavBP = _param;
		_predict();

		_initialized = true;

		if (_data.size() < _minsat) {
			if (_log) _log->comment(2, "gpvtflt", "Not enough visible satellites!");
			_restore(QsavBP, XsavBP);
			return -1;
		}

		// define a number of measurements 
		unsigned int nObs = _data.size();
		unsigned int mult = 1;
		if (_observ == RAW_DOUBLE) { mult = 2; nObs *= 2; }
		if (_observ == RAW_ALL) { mult = 2; nObs *= 5; } // reservation for 5 freq - not used raws will be removed
		if (_phase) { mult *= 2; nObs *= 2; }  // code + phase
		// if (_doppler) { mult; nObs *= 2; }
		// nObs *= 2;
#ifdef DEBUG
		cout << "gpvtflt: nObs: " << _data.size() << " " << nObs << " " << runEpoch.str_hms() << endl; cout.flush();
#endif       

		/*G01 P1 P2 P3... L1 L2L3*/
		unsigned int nPar = _param.parNumber();

		A.ReSize(nObs, nPar); A = 0.0;
		l.ReSize(nObs);       l = 0.0;
		P.ReSize(nObs);       P = 0.0;
		dx.ReSize(nPar);     dx = 0.0;
		unsigned int iobs = 1;
		_frqNum.clear(); _obs_index.clear();
		it = _data.begin();
		while (it != _data.end())
		{
			if (_addObsP((*it), iobs, A, l, P) > 0) {}
			else { it = _data.erase(it); continue; }
			if (_addObsL((*it), iobs, A, l, P) > 0) {}
			else { it = _data.erase(it); continue; }
			if (_observ == IONO_FREE) {
				if (_combineMW((*it)) < 0) { it = _data.erase(it); continue; }
			}
			++it;
		}

		if (iobs <= _minsat * mult) {
			if (_log) _log->comment(2, "gpvtflt", "Not enough processable observations!");
			_restore(QsavBP, XsavBP);
			return -1;
		}

		// delete zero rows/cols
		iobs--;
		A = A.Rows(1, iobs);
		P = P.SymSubMatrix(1, iobs);
		l = l.Rows(1, iobs);

		// rtk mode or single diff(sat) mode
		if (_isBase || _sd_sat)
		{
			if (_combineDD(A, P, l) < 0) return -1;
		}

		Qsav = _Qx;

		_filter->update(A, P, l, dx, _Qx);

		// increasing variance after update in case of introducing new ambiguity
		for (size_t iPar = 0; iPar < _param.parNumber(); iPar++) {
			if (_param[iPar].parType == t_gpar::AMB_IF) {
				string sat = _param[iPar].prn;
				if (_newAMB.find(sat) != _newAMB.end() && _cntrep == 1) {
					if (_newAMB[sat] == 1) _Qx(iPar + 1, iPar + 1) += 10;
					if (_newAMB[sat] == 2 && _Qx(iPar + 1, iPar + 1) > 0.01) _Qx(iPar + 1, iPar + 1) += 1;
					_newAMB[sat]++;
				}
			}
		}

		// post-fit residuals
		v_orig = l - A * dx;

		// normalized post-fit residuals
		Matrix Qv = A * _Qx * A.t() + P.i();
		v_norm.ReSize(v_orig.Nrows());
		for (int i = 1; i <= v_norm.Nrows(); i++) {
			v_norm(i) = sqrt(1 / Qv(i, i)) * v_orig(i);
		}

		int freedom = A.Nrows() - A.Ncols();
		if (freedom < 1) {
			if (_log) _log->comment(1, "t_gpvtflt", "No redundant observations!");
			//    _restore(QsavBP, XsavBP);     
			//    return -1;
			freedom = 1;
		}

		ColumnVector vtPv = v_orig.t() * P * v_orig;
		_sig_unit = vtPv(1) / freedom;
		vtpv = vtPv(1);
		nobs_total = A.Nrows();
		npar_number = A.ncols();

#ifdef DEBUG
		//    SymmetricMatrix Corr = cov2corr(_Qx);
		cout << "Used satellites: ";
		for (unsigned int i = 0; i < _data.size(); i++) cout << _data[i].sat() << " ";
		cout << endl;
		cout << "Used parameters: ";
		for (unsigned int i = 0; i < _param.parNumber(); i++) cout << _param[i].str_type() << "_" << _param[i].prn << " ";
		cout << endl;
		cout << fixed << setprecision(3)
			<< _epoch.str_ymdhms() << endl
			//         << " m0: " << setw(7) << _sig_unit << endl
			<< " dx: " << setw(7) << dx
			//         << " l.Nrows(): " << l.Nrows() << endl
			//         << setw(7) << l 
			//         << " v_orig.Nrows(): " << v_orig.Nrows() 
			//         << " v_orig: "  << setw(7) << v_orig
			//         << " v_norm: "  << setw(7) << v_norm
			//         << " A(" << A.Nrows() << "," << A.Ncols() << "):" << setw(7) << endl << A
			//         << " P: "  << setw(7) << P
			//         << " Q: "  << setw(7) << _Qx
			//         << " C: "  << setw(7) << Corr    
			<< " Param " << setw(7) << _param
			<< endl; cout.flush();

		// for(int i = 1; i <= _param.parNumber(); i++) cout << _param[_param.getParam(i)].str_type() << " "; cout << endl;
		// for(int i = 1; i <= _param.parNumber(); i++) cout << _param[_param.getParam(i)].str_type() << " " << Corr.row(i) << " "; cout << endl;

		//int ooo; cin >> ooo;
#endif
		// save post-fit residuals
		//if (!_isBase) 
		//{
		//	_save_residuals(v_norm, _data, RES_NORM);
		//	_save_residuals(v_orig, _data, RES_ORIG);
		//}

		// std::cout << _epoch.sod() << " " << _sig_unit << std::endl;

	} while (_outlierDetect(v_norm, Qsav, outlier) != 0);

	_nSat_excl = _nSat - _data.size();   // number of excluded satellites due to different reasons

	try {
		_n_ALL_flt++;
		Cholesky(_Qx);
	}
	catch (NPDException) {
		//     cout << _site << " " << _epoch.str_ymdhms() << " NPDException " << endl;
		_n_NPD_flt++;
	}
	// end of covariance matrix testing     

	if (_data.size() < _minsat) {
		if (_log) _log->comment(2, "gpvtflt", _site + _epoch.str_ymdhms(" epoch ")
			+ " skipped: " + int2str(_data.size()) + " < _minsat)");
		//    cout << "Epoch " << _epoch.str_hms() << " not calculated due to _minsat." << endl;     
		_restore(QsavBP, XsavBP);
		return -1;
	}

	// the velocity is not estimated by the real-time model (suggested by zhshen !!!)
	if (_doppler)
	{
		int irc = _processEpochVel();
		if (irc < 0)
		{
			if (_log) _log->comment(1, "t_gpvtflt", "estimate velocity by doppler failed!");
			_vel = 0; _Qx_vel = 0;
		}
	}

	ColumnVector dx_tmp = dx;
	SymmetricMatrix   Qx_tmp = _Qx;
	t_gallpar _param_tmp = _param;
	// ambfix added by zhshen
	if (_fix_mode != t_gsetamb::FIX_MODE::NO)
	{
		if (_isFirstFix)
		{
			// set glofrq number 
			_glofrq_num = _gobs->glo_freq_num();
			if (_glofrq_num.empty() && _gnav)
			{
				_glofrq_num = _gnav->glo_freq_num();
			}
			_ambfix->setWaveLength(_glofrq_num);
			_ambfix->setUPD(_gupd);
			_ambfix->setLOG(_log);
			//_ambfix->setFreq(_gnav);
			_isFirstFix = false;
		}

		_filter->add_data(_param, dx, _Qx, _sig_unit, Qsav);
		_filter->add_data(vtpv, nobs_total, npar_number);
		_ambfix->setObsType(_observ);
		if (_observ == IONO_FREE)	_ambfix->setMW(_MW[_epoch]);
		t_gamb_cmn amb_cmn(_epoch, _filter);
		map<string, int> active_amb; active_amb[_site] = npar_number;
		//amb_cmn.active_amb(active_amb);
#ifdef DEBUG
		ColumnVector tmp = _filter->dx();
		string out = "Before amb: \n";
		for (int i = 1; i <= tmp.Nrows(); i++)
		{
			out += dbl2str(tmp(i), 8);
		}
		if (_log) _log->comment(2, out);

#endif

		if (_ambfix->processBatch(_epoch, _filter, nullptr, "NL") < 0)
		{
			amb_cmn.amb_fixed(false);
			if (_log) _log->comment(2, "gpvtflt", _site + _epoch.str_ymdhms(" epoch ")
				+ ": fix ambiguity failed !");
		}
		else amb_cmn.amb_fixed(true);
#ifdef DEBUG
		ColumnVector tmp1 = _filter->dx();
		string out1 = "After amb: \n";
		for (int i = 1; i <= tmp1.Nrows(); i++)
		{
			out1 += dbl2str(tmp1(i), 8);
		}
		if (_log) _log->comment(2, out1);
#endif
		dx_tmp = amb_cmn.dx();
		Qx_tmp = amb_cmn.Qx();

		//if (amb_cmn.amb_fixed())_amb_state = true;
		//else _amb_state = false;
		_amb_state = amb_cmn.amb_fixed();
		//_fixed_pct = amb_cmn.fixed_percent();
	}

	// edit and save ZHD - for post-fit residual
	// (Must be before updating crd due to maintain consistency with a priory)
	for (unsigned int iPar = 0; iPar < _param.parNumber(); iPar++) {
		if (_param[iPar].parType == t_gpar::TRP)
		{
			t_gtriple Ell, XYZ;
			string site = _param[iPar].site;
			if (_param.getCrdParam(site, XYZ) < 0)
			{
				XYZ = _gallobj->obj(site)->crd_arp(_epoch);
			}
			xyz2ell(XYZ, Ell, false);
			if (site == _site && _gModel->tropoModel() != 0)
				_param[iPar].apriori(_gModel->tropoModel()->getZHD(Ell, _epoch));
			else if (site == _site_base && _gModel_base->tropoModel() != 0)
				_param[iPar].apriori(_gModel_base->tropoModel()->getZHD(Ell, _epoch));
		}
	}

	// Edit and save parematres (transmit the float info!!!)
	for (unsigned int iPar = 0; iPar < _param.parNumber(); iPar++) {
		_param[iPar].value(_param[iPar].value() + dx(_param[iPar].index));
	}

	// output the fixed result
	ostringstream os;
	for (unsigned int iPar = 0; iPar < _param.parNumber(); iPar++) {
		_param_tmp[iPar].value(_param_tmp[iPar].value() + dx_tmp(_param_tmp[iPar].index));
	}
	_prtOut(_epoch, _param_tmp, Qx_tmp, _data, os, line, true);
	// Print flt results
	if (_flt) {
		_flt->write(os.str().c_str(), os.str().size());
		_flt->flush();
	}

	return 1;
}

int gsins::t_gpvtflt_old::_outlierDetect(const ColumnVector & v, const SymmetricMatrix& Qsav, string & sat)
{
	gtrace("t_gpvtflt::outlierDetect");

	int nobs = v.Nrows(); double max = 0.0; int idx = 0;
	for (int i = 1; i <= nobs; i++)
	{
		if (fabs(v(i)) > max && fabs(v(i)) > _max_res_norm) {
			max = fabs(v(i));
			idx = i;
		}
	}
	if (idx > 0) {
		_Qx = Qsav;
		sat = _obs_index[idx - 1].first;
		string obsType = gobstype2str(_obs_index[idx - 1].second.second);
		ostringstream os;
		os << _site << " outlier (" << obsType << ") " << sat
			<< " v: " << fixed << setw(16) << right << setprecision(3) << max;
		if (_log)
			_log->comment(1, "gpvtflt", os.str());
	}
	else sat = "";

	return idx;
}

int gsins::t_gpvtflt_old::_outlierDetect(const Eigen::VectorXd & v, const SymmetricMatrix & Qsav, string & sat)
{
	gtrace("t_gpvtflt::outlierDetect");

	int nobs = v.size(); double max = 0.0; int idx = 0;
	for (int i = 0; i < nobs; i++)
	{
		if (fabs(v(i)) > max && fabs(v(i)) > _max_res_norm) {
			max = fabs(v(i));
			idx = i;
		}
	}
	if (idx > 0) {
		_Qx = Qsav;
		sat = _obs_index[idx].first;
		string obsType = gobstype2str(_obs_index[idx - 1].second.second);
		ostringstream os;
		os << _site << " outlier (" << obsType << ") " << sat
			<< " v: " << fixed << setw(16) << right << setprecision(3) << max;
		if (_log)
			_log->comment(1, "gpvtflt", os.str());
	}
	else sat = "";

	return idx;
}

string gsins::t_gpvtflt_old::_getsat_from_obs(int idx)
{
	gtrace("t_gpvtflt::_getsat_from_obs");



	return string();
}

void gsins::t_gpvtflt_old::_print_ele()
{
	gtrace("t_gpvtflt::_print_ele");

	ofstream file("ele.txt", ios::app);

	streambuf *stream_buffer_cout = std::cout.rdbuf();
	streambuf *stream_buffer_file = file.rdbuf();
	std::cout.rdbuf(stream_buffer_file);

	std::cout << _epoch.sow() + _epoch.dsec() << "  " << _data.size() << "  ";
	for (vector<t_gsatdata>::iterator it = _data.begin(); it != _data.end(); it++)
	{
		std::cout << it->sat() << ":" << it->ele_deg() << "  ";
	}
	std::cout << endl;
	std::cout.rdbuf(stream_buffer_cout);

}

int gsins::t_gpvtflt_old::processBatch(const t_gtime& beg_r, const t_gtime& end_r, bool prtOut)
{
	gtrace("t_gpvtflt::processBatch");

#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif 
	_gmutex.lock();

	if (_grec == nullptr) {
		ostringstream os;
		os << "ERROR: No object found (" << _site << "). Processing terminated!!! " << beg_r.str_ymdhms() << " -> " << end_r.str_ymdhms() << endl;
		if (_log) _log->comment(0, "t_gpppflt", os.str());
		return -1;
	}

	int sign = 1;
	t_gtime begT;
	t_gtime endT;

	if (!_beg_end) {
		begT = end_r;
		endT = beg_r;
		sign = -1;
		if (_log) _log->comment(2, "t_gpppflt", "Filtering in end -> begin direction!");
	}
	else {
		begT = beg_r;
		endT = end_r;
		if (_log) _log->comment(2, "t_gpppflt", "Filtering in begin -> end direction!");
	}

	// get the time
	_beg_time = (_beg_time < begT) ? begT : _beg_time;
	_end_time = (_end_time > endT) ? endT : _end_time;

	if (_isBase) {
		if (_setCrd() < 0) return -1;
	}
	//grec->crd(t_gtriple(0), t_gtriple(1, 1, 1), _beg_time, _end_time);

	_dop.set_data(_gnav, _gobs, _site);
	_running = true;
	for (int i = 0; i < _param.parNumber(); i++)
		if (_param[i].site == "") _param.setSite(_site);

	double subint = 0.1;
	if (_scale > 0) subint = 1.0 / _scale;
	if (_sampling > 1) subint = pow(10, floor(log10(_sampling)));

	string smtModStr(dynamic_cast<t_gsetflt*>(_set)->method_smt());
	string fltModStr(dynamic_cast<t_gsetflt*>(_set)->method_flt());
	int smt_delay = dynamic_cast<t_gsetflt*>(_set)->smt_delay();

	t_gtime now(begT);

#ifdef DEBUG
	if (_log) _log->comment(1, "t_gpppflt", begT.str_ymdhms("Filtering: ") + endT.str_ymdhms(" "));
#endif

	std::cout << _site << ": Start filtering: " << _site << " " << now.str_ymdhms() << " " << endT.str_ymdhms() << endl;
	_prtOutHeader();
	bool time_loop = true;
	//   high_resolution_clock::time_point chrono_t1 = high_resolution_clock::now();
	while (time_loop)
	{
		if (_beg_end && (now < endT || now == endT)) time_loop = true;
		else if (_beg_end && now > endT) { time_loop = false; break; }

		if (!_beg_end && (now > endT || now == endT)) time_loop = true;
		else if (!_beg_end && now < endT) { time_loop = false; break; }

		// synchronization
		if (now != endT) {
			if (!time_sync(now, _sampling, _scale, _log)) { // now.sod()%(int)_sampling != 0 ){      
				now.add_dsec(sign * subint / 100); // add_dsec used for synchronization!

				continue;
			}
			if (_sampling > 1) now.reset_dsec();

		}

		// added by zhShen
		if (dynamic_cast<t_gsetproc*>(_set)->realtime())
		{
			while (!((t_bncobs*)_gcoder[0])->available(now) || !((t_bnccorr*)_gcoder[1])->available(now))
				t_gtime::gmsleep(1);

			t_gtime epo_prev = now - 10.0;
			if (dynamic_cast<t_gsetproc*>(_set)->slip_model() == SLIPMODEL::TURBO_EDIT)
			{
				t_gexeturboedit gturbo(_set, _log, _gobs, _gnav);
				gturbo.ProcessBatch(_site, epo_prev, now, _sampling);
				if (_isBase)	gturbo.ProcessBatch(_site_base, epo_prev, now, _sampling);
			}
			else
			{
				t_gpreproc prep(_gobs, _set);
				prep.glog(_log);
				prep.preprocess(_site, epo_prev, now, _sampling, false);
				if (_isBase) prep.preprocess(_site_base, epo_prev, now, _sampling, false);
			}


			_gobs->erase(_site, now);
		}

		// clean/collect/filter epoch data
		_data.erase(_data.begin(), _data.end());
		_data = _gobs->obs(_site, now);

		if (_data.size() == 0) {
			if (_log && _log->verb() >= 1) {
				_log->comment(1, "gpppflt", _site + now.str_ymdhms(" no observation found at epoch: "));
			}
			if (_sampling > 1) now.add_secs(int(sign * _sampling));  // =<1Hz data
			else               now.add_dsec(sign * _sampling);       //  >1Hz data

			continue;
		}
		if (_isBase)
		{
			_data_base.erase(_data_base.begin(), _data_base.end());
			_data_base = _gobs->obs(_site_base, now);

			if (_data_base.size() == 0) {
				if (_log && _log->verb() >= 1) {
					_log->comment(1, "grtkflt", _site_base + now.str_ymdhms(" no observation found at epoch: "));
				}
				if (_sampling > 1) now.add_secs(int(sign * _sampling));  // =<1Hz data
				else               now.add_dsec(sign * _sampling);       //  >1Hz data

				continue;
			}

			//pick common satellites

		}

		t_gtime obsEpo = _data.begin()->epoch();

		_timeUpdate(obsEpo);

		// reseting params according to conf
		if (_reset_amb > 0) {
			if (now.sod() % _reset_amb == 0) {
				_reset_ambig();
			}
		}
		if (_reset_par > 0) {
			if (now.sod() % _reset_par == 0) {
				_reset_param();
			}
		}

		// save apriory coordinates     
		if(_crd_est != FIX)
			_saveApr(obsEpo, _param, _Qx);

		//#ifdef _WIN32
		//		system("cls");
		//#else
		//		system("clear");
		//#endif
		//		cerr << now.str() << endl;

				//int irc_epo = this->_processEpoch(obsEpo);
		int irc_epo = this->_processEpoch_new(obsEpo);

		if (irc_epo < 0) {
			_success = false;
			_removeApr(obsEpo);
			if (_log && _log->verb() >= 1)
				_log->comment(1, "gpppflt", _site + now.str_ymdhms(" epoch ") + " was not calculated");

			if (_sampling > 1) now.add_secs(int(sign * _sampling));  // =<1Hz data
			else               now.add_dsec(sign * _sampling);       //  >1Hz data


			if (dynamic_cast<t_gsetproc*>(_set)->slip_model() == SLIPMODEL::TURBO_EDIT)
			{
				t_gexeturboedit gturbo(_set, _log, _gobs, _gnav);
				t_gtime epo_prev = now;
				if (_sampling > 1) epo_prev.add_secs(-int(sign * 2 * _sampling));  // =<1Hz data
				else               epo_prev.add_dsec(-sign * 2 * _sampling);       //  >1Hz data     

				gturbo.ProcessBatch(_site, epo_prev, now, _sampling * 2);
			}
			else
			{
				// repeat preprocessing in case of skipped epoch -> sampling changed
				t_gpreproc prep(_gobs, _set);
				prep.glog(_log);
				t_gtime epo_prev = now;
				if (_sampling > 1) epo_prev.add_secs(-int(sign * 2 * _sampling));  // =<1Hz data
				else               epo_prev.add_dsec(-sign * 2 * _sampling);       //  >1Hz data     
				prep.preprocess(_site, epo_prev, now, _sampling * 2, false);
				//       cout << "New preprocessing due to skipping epoch: " << epo_prev.str_hms() << " " << now.str_hms() << " " << _sampling*2 << endl;
			}

			continue;
		}
		else _success = true;

		if (_log && _log->verb() >= 1)
			_log->comment(1, "gpppflt", _site + now.str_ymdhms(" processing epoch: "));

		//ostringstream os;
		//bool savePrd = true;
		//if (_smooth) savePrd = false;
		//_prtOut(obsEpo, _param, _Qx, _data, os, line, savePrd);

		//// Print flt results
		//if (_flt && prtOut) {
		//	_flt->write(os.str().c_str(), os.str().size());
		//	_flt->flush();
		//}

		//_initialized = true;

		double percent = now.diff(begT) / endT.diff(begT) * 100.0;
		cerr << "\r" << now.str_ymdhms() << setw(5) << " Q = " << (_amb_state ? 1 : 2) << fixed << setprecision(1) << setw(6) << percent << "%";

		if (_sampling > 1) now.add_secs(int(sign * _sampling));  // =<1Hz data
		else               now.add_dsec(sign * _sampling);       //  >1Hz data

															   // RTS Smooth
		if (_smooth) {
			if (smt_delay != 0 && now.sod() % smt_delay == 0) {
				_backflt(prtOut);
				_vfltdata.clear();
			}
		}
	}

	_running = false;
	//   high_resolution_clock::time_point chrono_t2 = high_resolution_clock::now();
	//   auto duration = duration_cast<microseconds>( chrono_t2 - chrono_t1 ).count();
	//   cout << "Duration FLT: " << duration << endl;

	//   chrono_t1 = high_resolution_clock::now();
	if (_smooth && smt_delay == 0) {
		_backflt(prtOut);
	}
	//   chrono_t2 = high_resolution_clock::now();
	//   duration = duration_cast<microseconds>( chrono_t2 - chrono_t1 ).count();
	//   cout << "Duration SMT: " << duration << endl;

	if (beg_r != end_r) {   // do not print in real-time
		double npd_perc = 0;
		npd_perc = (double(_n_NPD_flt) / double(_n_ALL_flt)) * 100;
		if (_log) _log->comment(1, "gpppflt", _site + ": Forward filter -     Not positive definite covariance matrices: " + dbl2str(npd_perc, 0) + " %   ("
			+ int2str(_n_NPD_flt) + ", " + int2str(_n_ALL_flt) + ")");

		if (_smooth) {
			npd_perc = (double(_n_NPD_smt) / double(_n_ALL_smt)) * 100;
			if (_log) _log->comment(1, "gpppflt", _site + ": Backward smoothing - Not positive definite covariance matrices: " + dbl2str(npd_perc, 0) + " %   ("
				+ int2str(_n_NPD_smt) + ", " + int2str(_n_ALL_smt) + ")");
		}
	}

	_gmutex.unlock();

	return 1;
}

int gsins::t_gpvtflt_old::processBatchFB(const t_gtime& beg_r, const t_gtime& end_r)
{
	gtrace("t_gpvtflt::processBatchFB");

	this->processBatch(beg_r, end_r, false);

	_beg_end = false;

	this->processBatch(beg_r, end_r, true);

	_beg_end = true;

	return 1;
}

void gsins::t_gpvtflt_old::Add_UPD(t_gupd * gupd)
{
	_gupd = gupd;
}

int gsins::t_gpvtflt_old::_selcomsat(vector<t_gsatdata>& data_base, vector<t_gsatdata>& data_rover)
{
	int nSat = 0;
	//com_sat.clear(); //mapobsdata.clear();
	vector<t_gsatdata> bdata, rdata;
	//first to recycle on base-site observations
	for (auto iter_base = data_base.begin(); iter_base != data_base.end(); iter_base++)
	{
		string comsat = iter_base->sat();
		//second to recycle on rover-site obserations
		for (auto iter_rover = data_rover.begin(); iter_rover != data_rover.end(); iter_rover++)
		{
			if (iter_rover->sat() == comsat /*&& _set->exsats.find(comsat) == _set->exsats.end()*/)
			{
				//mapobsdata[_site_base][comsat] = *iter_base;
				//mapobsdata[_site][comsat] = *iter_rover;
				bdata.push_back(*iter_base);
				rdata.push_back(*iter_rover);
				//com_sat.insert(comsat);
				nSat++;
				break;
			}
		}
	}
	if (nSat)
	{
		data_base = bdata;
		data_rover = rdata;
	}
	return nSat;
}

void gsins::t_gpvtflt_old::_udsdAmb()
{
	gtrace("t_gpppflt::_syncAmb");

	for (map<string, int>::iterator it = _newAMB.begin(); it != _newAMB.end(); ) {
		if (it->second > 2) _newAMB.erase(it++);
		else ++it;
	}

	// Add ambiguity parameter and appropriate rows/columns covar. matrix   
	set<string> mapPRN;
	for (int i = 0; i < _data.size(); i++)
	{
		auto rsatdata = _data[i];
		auto bsatdata = _data_base[i];

		GSYS gs = rsatdata.gsys();
		string sat = rsatdata.sat();
		mapPRN.insert(rsatdata.sat());

		bool com = true;
		if (_gnav) com = _gnav->com();
		vector < GOBSBAND > band = dynamic_cast<t_gsetgnss*>(_set)->band(gs);

		if (_observ == IONO_FREE) {
			GOBSBAND b1, b2;
			if (_auto_band) {    // automatic dual band selection -> for Anubis purpose
				set<GOBSBAND> bands = rsatdata.band_avail();
				auto itBAND = bands.begin();
				if (bands.size() < 2) continue;
				b1 = *itBAND;
				itBAND++;
				b2 = *itBAND;
			}
			else {            // select fix defined band according the table 
				if (band.size())
				{
					b1 = band[0];
					b2 = band[1];
				}
				else {
					b1 = t_gsys::band_priority(gs, FREQ_1);
					b2 = t_gsys::band_priority(gs, FREQ_2);
				}
			}
			bool index = true; double amb[2] = { 0 };
			t_gobs gobs1, gobs2;
			for (int j = 0; j < 2; j++)
			{
				t_gsatdata tmpsatdata = rsatdata;
				if (j == 1) tmpsatdata = bsatdata;
				gobs1 = tmpsatdata.select_phase(b1);
				gobs2 = tmpsatdata.select_phase(b2);
				if (gobs1.gobs() == X || gobs2.gobs() == X) {
					index = false; break;
				}
				amb[j] = tmpsatdata.L3(gobs1, gobs2);
				if (double_eq(amb[j], 0.0)) {
					index = false; break;
				}
				if (_param.getParam(_site, t_gpar::AMB_IF, sat) < 0) {
					t_gpppmodel* gModel;
					gModel = new t_gpppmodel(tmpsatdata.site(), _log, _set);
					gModel->reset_observ(_observ);
					gModel->setOBJ(_gallobj);
					t_gallpar param;
					double cmpObs = gModel->cmpObs(_epoch, param, tmpsatdata, gobs1, com);
					delete gModel;
					if (cmpObs < 0) {
						index = false; break;
					}
					amb[j] -= tmpsatdata.P3(gobs1, gobs2);
				}
				else { index = false; break; }
			}
			if (index == false) continue;
			double sdamb = amb[0] - amb[1];
			t_gpar newPar(_site, t_gpar::AMB_IF, _param.parNumber() + 1, sat);
			newPar.value(sdamb);      // first ambiguity value

			_param.addParam(newPar);
			_newAMB[sat] = 1;
			Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
			_Qx(_param.parNumber(), _param.parNumber()) = _sigAmbig * _sigAmbig;
			//cout << "ADD AMB IF " << it->sat() << " " << it->epoch().str_ymdhms() << endl;
			//if (_log) _log->comment(3, "gppp", "AMB_IF was added! For Sat PRN " + sat + " Epoch: " + _epoch.str_ymdhms());
		}
		else if (_observ == RAW_SINGLE || _observ == RAW_DOUBLE || _observ == RAW_ALL) {
			int newAmb = 0;
			int nf = 5;
			if (_auto_band) {    // automatic dual band selection -> for Anubis purpose
				GOBSBAND b1, b2;
				set<GOBSBAND> bands = rsatdata.band_avail();
				auto itBAND = bands.begin();
				if (bands.size() < 2) continue;
				b1 = *itBAND;
				itBAND++;
				b2 = *itBAND;
				band.clear(); band.push_back(b1); band.push_back(b2);
			}
			if (band.size()) nf = band.size();
			for (FREQ_SEQ f = FREQ_1; f <= nf; f = (FREQ_SEQ)(f + 1))
			{
				GOBSBAND b;
				t_gpar::t_type amb_type;
				switch (f)
				{
					/*case gnut::FREQ_X:	amb = t_gpar::AMB_IF; break;*/
				case gnut::FREQ_1:	amb_type = t_gpar::AMB_L1; break;
				case gnut::FREQ_2:	amb_type = t_gpar::AMB_L2; break;
				case gnut::FREQ_3:	amb_type = t_gpar::AMB_L3; break;
				case gnut::FREQ_4:	amb_type = t_gpar::AMB_L4; break;
				case gnut::FREQ_5:	amb_type = t_gpar::AMB_L5; break;
				default:	break;
				}
				// '>' -> '>=' . Modified by YuxuanZhou. 2020/05/31
				if (band.size() >= f) {    // automatic dual band selection -> for Anubis purpose
					b = band[f - 1];
				}
				else b = t_gsys::band_priority(gs, f);


				bool index = true; double amb[2] = { 0 };
				amb[0] = amb[1] = 0.0;
				t_gobs gobs;
				for (int j = 0; j < 2; j++)
				{
					t_gsatdata tmpsatdata = rsatdata;
					if (j == 1) tmpsatdata = bsatdata;
					//_getgobs(sat, TYPE_L, b, gobs);
					gobs = tmpsatdata.select_phase(b);
					if (gobs.gobs() == X) { index = false; break; }
					amb[j] = tmpsatdata.obs_L(gobs);
					if (double_eq(amb[j], 0.0)) { index = false; break; }
					if (_observ == RAW_SINGLE && f > FREQ_1) { index = false; break; }
					if (_observ == RAW_DOUBLE && f > FREQ_2) { index = false; break; }

					// add L1 amb - everytime when RAW observations
					if (_param.getParam(_site, amb_type, sat) < 0) {
						t_gpppmodel* gModel;
						gModel = new t_gpppmodel(tmpsatdata.site(), _log, _set);
						gModel->reset_observ(_observ);
						gModel->setOBJ(_gallobj);
						t_gallpar param;
						double cmpObs = gModel->cmpObs(_epoch, param, tmpsatdata, gobs, com);
						delete gModel;
						if (cmpObs < 0) {
							index = false; break;
						}
						amb[j] -= tmpsatdata.obs_C(gobs);
					}
					else { index = false; break; }
				}//end base+rove
				if (index == false) continue;
				double sdamb = amb[0] - amb[1];

				t_gpar newPar(_site, amb_type, _param.parNumber() + 1, sat);
				newPar.value(sdamb);
				_param.addParam(newPar);
				Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
				_Qx(_param.parNumber(), _param.parNumber()) = _sigAmbig * _sigAmbig;
				//         cout << "ADD AMB L1 " << it->sat() << " " << it->epoch().str_ymdhms()  << " " << L1 - cmpObs << endl;
				//if (_log) _log->comment(3, "gppp", "RAW AMB_L1 was added! For Sat PRN " + it->sat() + " Epoch: " + _epoch.str_ymdhms());
				newAmb = 1;
			}//end f
			if (newAmb) _newAMB[sat] = 1;
		}//end uc

	}//end data

  // Remove ambiguity parameter and appropriate rows/columns covar. matrix

	for (unsigned int i = 0; i < _param.parNumber(); i++) {
		if (_param[i].parType == t_gpar::AMB_IF ||
			_param[i].parType == t_gpar::AMB_L1 ||
			_param[i].parType == t_gpar::AMB_L2 ||
			_param[i].parType == t_gpar::AMB_L3 ||
			_param[i].parType == t_gpar::AMB_L4 ||
			_param[i].parType == t_gpar::AMB_L5) {

			set<string>::iterator prnITER = mapPRN.find(_param[i].prn);
			if (prnITER == mapPRN.end()) {

				if (_log)
					_log->comment(3, "gppp", "AMB will be removed! For Sat PRN " + _param[i].prn + " Epoch: " + _epoch.str_ymdhms());

#ifdef DEBUG
				cout << "AMB will be removed! For Sat PRN " << _param[i].prn
					<< " Epoch: " << _epoch.str_ymdhms() << endl;
#endif        

				_newAMB.erase(_param[i].prn);

				Matrix_remRC(_Qx, _param[i].index, _param[i].index);
				_param.delParam(i);
				_param.reIndex();
				i--;
			}
		}
	}

	return;
}

void gsins::t_gpvtflt_old::Add_rho_azel(/*const t_gtime & runEpoch,*/ const string & site_name, const t_gtriple & xyz_s, const t_gtriple & xyz_r, t_gsatdata & obs_sat)
{
	t_gtriple xyz_rho = xyz_s - xyz_r;
	t_gtriple ell_r, neu_s;
	gnut::xyz2ell(xyz_r, ell_r, false);

	xyz2neu(ell_r, xyz_rho, neu_s);

	// Earth rotation correction
	t_gtriple xRec;
	double rho0 = sqrt(pow(xyz_r[0] - xyz_s[0], 2) + pow(xyz_r[1] - xyz_s[1], 2) + pow(xyz_r[2] - xyz_s[2], 2));
	double dPhi = OMEGA * rho0 / CLIGHT;
	xRec[0] = xyz_r[0] * cos(dPhi) - xyz_r[1] * sin(dPhi);
	xRec[1] = xyz_r[1] * cos(dPhi) + xyz_r[0] * sin(dPhi);
	xRec[2] = xyz_r[2];

	// Apply tides 
	if (_isBase && site_name == _site_base)
	{
		string tmp = _site;
		_site = site_name;
		_apply_tides(_epoch, xRec);
		_site = tmp;
	}
	else _apply_tides(_epoch, xRec);

	double tmp = (obs_sat.satcrd().crd_cvect() - xRec.crd_cvect()).norm_Frobenius();

	obs_sat.addrho(tmp);
	//cout << iter->sat() << " " << iter->epoch().str_ymdhms() << fixed << setprecision(3) << " " << neu_s[0] << " " << neu_s[1] << " " << neu_s[2] << endl;
	double NE2 = neu_s[0] * neu_s[0] + neu_s[1] * neu_s[1];
	double ele = acos(sqrt(NE2) / tmp);
	if (neu_s[2] < 0.0) {
		ele *= -1.0;
	}

	if (sqrt(NE2) / tmp > 1.0)
		obs_sat.addele(0.0);
	else obs_sat.addele(ele);

	double azi = atan2(neu_s[1], neu_s[0]);
	if (azi < 0) azi += 2 * G_PI;
	obs_sat.addazi(azi);
	//if (!_use_ecl) obs_sat.addecl(_lastEcl);

	return;
}

void gsins::t_gpvtflt_old::_prtOut(t_gtime& epoch, t_gallpar& X, const SymmetricMatrix& Q, vector<t_gsatdata>& data, ostringstream& os, xml_node& node, bool saveProd)
{
	gtrace("t_gpvtflt::_prtOut");

	// get CRD params
	t_gtriple xyz, ell;
	X.getCrdParam(_site, xyz);
	xyz2ell(xyz, ell, false);

	// CRD using eccentricities
	t_gtriple xyz_ecc = xyz - _grec->eccxyz(epoch); // MARKER + ECC = ARP

													// get CRD rms  (XYZ)
	double Xrms = 0.0, Yrms = 0.0, Zrms = 0.0,
		Vxrms = 0.0, Vyrms = 0.0, Vzrms = 0.0;
	double cov_xy = 0.0, cov_xz = 0.0, cov_yz = 0.0;
	if (Q(1, 1) < 0) Xrms = -1; else Xrms = sqrt(Q(1, 1));
	if (Q(2, 2) < 0) Yrms = -1; else Yrms = sqrt(Q(2, 2));
	if (Q(3, 3) < 0) Zrms = -1; else Zrms = sqrt(Q(3, 3));
	t_gtriple crd_rms(Xrms, Yrms, Zrms);
	t_gtriple vRec(0, 0, 0);
	cov_xy = Q(2, 1);
	cov_xz = Q(3, 1);
	cov_yz = Q(3, 2);

	if (_doppler)
	{
		vRec = _vel;
		if (_Qx_vel(1, 1) < 0) Vxrms = -1; else Vxrms = sqrt(_Qx_vel(1, 1));
		if (_Qx_vel(2, 2) < 0) Vyrms = -1; else Vyrms = sqrt(_Qx_vel(2, 2));
		if (_Qx_vel(3, 3) < 0) Vzrms = -1; else Vzrms = sqrt(_Qx_vel(3, 3));
	}
	else
	{
		vRec = t_gtriple(0.0, 0.0, 0.0);
	}
	// get CLK param
	double clk = 0.0;
	int iclk = X.getParam(_site, t_gpar::CLK, "");
	if (iclk >= 0) clk = X[iclk].value();

	// get ZTD params
	double ZHD = 0.0, ZWD = 0.0, ZTD = 0.0, ZWDrms = 0.0;
	int iztd = X.getParam(_site, t_gpar::TRP, "");
	if (iztd >= 0) {
		ZHD = X[iztd].apriori();//_gModel->tropoModel()->getZHD(ell, epoch);
		ZWD = X[iztd].value();
		ZTD = ZHD + ZWD;
		if (Q(iztd + 1, iztd + 1) < 0) ZWDrms = -1;
		else ZWDrms = sqrt(Q(iztd + 1, iztd + 1));          // _param is indexed from 0, but Q from 1 !
	}
	else {
		// test purpose only (if ZTD not estimated)
		if (_gModel->tropoModel() != 0) {
			ZWD = _gModel->tropoModel()->getZWD(ell, epoch);
			ZHD = _gModel->tropoModel()->getZHD(ell, epoch);
			ZTD = ZHD + ZWD;
		}
	}

	// get Tropo gradient
	int ign = X.getParam(_site, t_gpar::GRD_N, "");
	int ige = X.getParam(_site, t_gpar::GRD_E, "");
	double grdN = 0.0, grdE = 0.0;
	//  double rmsN = 0.0, rmsE = 0.0;
	if (ign >= 0 && ige >= 0) {
		grdN = X[ign].value();
		grdE = X[ige].value();
		//    rmsN = sqrt(Q(ign+1,ign+1));                  // _param is indexed from 0, but Q from 1 !
		//    rmsE = sqrt(Q(ige+1,ige+1));                  // _param is indexed from 0, but Q from 1 !
	}

	set<string> sat_list;
	for (auto it = _data.begin(); it != _data.end(); it++) sat_list.insert(it->sat());
	_dop.set_sats(sat_list);
	double gdop = -1, pdop = -1;
	if (_dop.calculate(epoch, xyz) >= 0) gdop = _dop.gdop();
	if (_dop.calculate(epoch, xyz) >= 0) pdop = _dop.pdop();

	set<string> ambs = X.amb_prns();
	int nsat = ambs.size();

	// get amb status
	string amb = "Float";
	if (_amb_state)amb = "Fixed";

	if (_allprod != 0 && saveProd) {
		shared_ptr<t_gprod> prdcrd = _allprod->get(_site, t_gdata::POS, epoch);
		if (prdcrd) {
			dynamic_pointer_cast<t_gprodcrd>(prdcrd)->xyz(xyz_ecc);
			dynamic_pointer_cast<t_gprodcrd>(prdcrd)->xyz_rms(crd_rms);
			dynamic_pointer_cast<t_gprodcrd>(prdcrd)->cov(COV_XY, cov_xy);
			dynamic_pointer_cast<t_gprodcrd>(prdcrd)->cov(COV_XZ, cov_xz);
			dynamic_pointer_cast<t_gprodcrd>(prdcrd)->cov(COV_YZ, cov_yz);
		}
		else {
			prdcrd = make_shared<t_gprodcrd>(epoch);
			dynamic_pointer_cast<t_gprodcrd>(prdcrd)->xyz(xyz_ecc);
			dynamic_pointer_cast<t_gprodcrd>(prdcrd)->xyz_rms(crd_rms);
			dynamic_pointer_cast<t_gprodcrd>(prdcrd)->cov(COV_XY, cov_xy);
			dynamic_pointer_cast<t_gprodcrd>(prdcrd)->cov(COV_XZ, cov_xz);
			dynamic_pointer_cast<t_gprodcrd>(prdcrd)->cov(COV_YZ, cov_yz);
			_allprod->add(prdcrd, _site);
		}

		t_gtriple Ell, XYZ;
		if (_param.getCrdParam(_site, XYZ) > 0) {
		}
		else if (_valid_crd_xml) { XYZ = _grec->crd_arp(_epoch); }
		xyz2ell(XYZ, Ell, false);
		int itrp = _param.getParam(_site, t_gpar::TRP, "");
		if (itrp >= 0) _param[itrp].apriori(_gModel->tropoModel()->getZHD(Ell, _epoch));
	}
	t_gtriple blh;
	xyz2ell(xyz_ecc, blh, true);
	// write params
	if (_kml) {
		string val = dbl2str(ell[1] * R2D, 11) + "," + dbl2str(ell[0] * R2D, 11);
		this->_default_node(node, "coordinates", val.c_str(), false);
	}

	string str_dsec = dbl2str(epoch.dsec());

	os << fixed << setprecision(4) << " "
		// << epoch.str_ymdhms() << str_dsec.substr(2) << setprecision(4)
		<< " " << epoch.sow() + epoch.dsec()
		<< fixed << setprecision(4)
		<< " " << setw(15) << xyz_ecc[0]          // [m]
		<< " " << setw(15) << xyz_ecc[1]          // [m]
		<< " " << setw(15) << xyz_ecc[2]          // [m]
		<< " " << setw(10) << vRec[0]          // [m/s]
		<< " " << setw(10) << vRec[1]          // [m/s]
		<< " " << setw(10) << vRec[2]          // [m/s]
		//<< " " << setw(12) << dclk_Rec
		//<< " " << setw(9) << ZTD             // [m]
		//<< " " << setw(9) << grdN * 1000       // [mm]
		//<< " " << setw(9) << grdE * 1000       // [mm]
		//<< " " << setw(17) << clk //<< (clk/CLIGHT)*1e6             // [m]
		<< fixed << setprecision(4)
		<< " " << setw(9) << Xrms            // [m]
		<< " " << setw(9) << Yrms            // [m]
		<< " " << setw(9) << Zrms            // [m]
		<< " " << setw(9) << Vxrms           // [m/s]
		<< " " << setw(9) << Vyrms           // [m/s]
		<< " " << setw(9) << Vzrms           // [m/s]
		//<< " " << setw(7) << ZWDrms          // [m]
		<< fixed << setprecision(0)
		<< " " << setw(5) << nsat            // nsat
		<< fixed << setprecision(1)
		//<< " " << setw(3) << gdop            // gdop
		<< fixed << setprecision(2)
		<< " " << setw(5) << pdop            // pdop
		<< fixed << setprecision(2)
		<< " " << setw(8) << _sig_unit       // m0
		<< " " << setw(10) << amb
		//<< " " << setw(10) << _fixed_pct
		<< endl;

	return;
}

void gsins::t_gpvtflt_old::_get_result(t_gtime& epoch, t_gposdata::data_pos & pos)
{
	gtrace("t_gpvtflt::_get_result");

	double crt = epoch.sow() + epoch.dsec();
	t_gallpar X = _param;
	// get CRD params
	t_gtriple xyz, xyz_ecc;
	// get CRD rms  (XYZ)
	double Xrms = 0.0, Yrms = 0.0, Zrms = 0.0,
		Vxrms = 0.0, Vyrms = 0.0, Vzrms = 0.0;

	if (X.getCrdParam(_site, xyz) > 0)
	{
		// CRD using eccentricities
		xyz_ecc = xyz - _grec->eccxyz(epoch); // MARKER + ECC = ARP

		double cov_xy = 0.0, cov_xz = 0.0, cov_yz = 0.0;
		if (_Qx(1, 1) < 0) Xrms = -1; else Xrms = sqrt(_Qx(1, 1));
		if (_Qx(2, 2) < 0) Yrms = -1; else Yrms = sqrt(_Qx(2, 2));
		if (_Qx(3, 3) < 0) Zrms = -1; else Zrms = sqrt(_Qx(3, 3));
	}
	else {
		xyz_ecc = _grec->crd_arp(epoch);
		Xrms = Yrms = Zrms = 30;
	}
	// get VEL rms
	t_gtriple vRec(0, 0, 0);
	if (_doppler) {
		vRec = _vel;
		if (_Qx_vel(1, 1) < 0) Vxrms = -1; else Vxrms = sqrt(_Qx_vel(1, 1));
		if (_Qx_vel(2, 2) < 0) Vyrms = -1; else Vyrms = sqrt(_Qx_vel(2, 2));
		if (_Qx_vel(3, 3) < 0) Vzrms = -1; else Vzrms = sqrt(_Qx_vel(3, 3));
	}
	else {
		vRec = t_gtriple(0.0, 0.0, 0.0);
		Vxrms = Vyrms = Vzrms = 0;
	}

	set<string> sat_list;
	for (auto it = _data.begin(); it != _data.end(); it++) sat_list.insert(it->sat());
	_dop.set_sats(sat_list);
	double pdop = -1;
	if (_dop.calculate(epoch, xyz) >= 0) pdop = _dop.pdop();

	set<string> ambs = X.amb_prns();
	int nsat = ambs.size();

	Eigen::Vector3d Qpos(Xrms, Yrms, Zrms), Qvel(Vxrms, Vyrms, Vzrms);
	Eigen::Vector3d position(xyz_ecc[0], xyz_ecc[1], xyz_ecc[2]), velocity(vRec[0], vRec[1], vRec[2]);

	pos = t_gposdata::data_pos{ crt,position,velocity,Qpos,Qvel,pdop,nsat,_amb_state };
}

// added by zhshen
void gsins::t_gpvtflt_old::_prtOutHeader()
{
	gtrace("t_gpvtflt::_prtOutHeader");
	ostringstream os;

	// first line
	os << "#" <<
		setw(15) << "Seconds of Week" <<
		" " << setw(12) << "X-ECEF " <<          // [m]
		" " << setw(15) << "Y-ECEF" <<            // [m]
		" " << setw(15) << "Z-ECEF" <<           // [m]
		" " << setw(10) << "Vx-ECEF" <<         // [m/s]
		" " << setw(10) << "Vy-ECEF" <<         // [m/s]
		" " << setw(10) << "Vz-ECEF" <<         // [m/s]
											   //<< " " << setw(12) << dclk_Rec
											   //<< " " << setw(9) << ZTD             // [m]
											   //<< " " << setw(9) << grdN * 1000       // [mm]
											   //<< " " << setw(9) << grdE * 1000       // [mm]
											   //<< " " << setw(17) << clk //<< (clk/CLIGHT)*1e6             // [m]
		" " << setw(9) << "X-RMS" <<
		" " << setw(9) << "Y-RMS" <<
		" " << setw(9) << "Z-RMS" <<
		" " << setw(9) << "Vx-RMS" <<
		" " << setw(9) << "Vy-RMS" <<
		" " << setw(9) << "Vz-RMS" <<
		" " << setw(5) << "NSat" <<
		" " << setw(5) << "PDOP" <<
		" " << setw(8) << "sigma0" <<
		" " << setw(10) << "ambStatus" <<
		//" " << setw(10) << "fixed-num" <<
		endl;

	// second line
	os << "#" <<
		setw(15) << "(s)" <<
		" " << setw(12) << "(m)" <<         // [m]
		" " << setw(15) << "(m)" <<         // [m]
		" " << setw(15) << "(m)" <<         // [m]
		" " << setw(10) << "(m/s)" <<         // [m/s]
		" " << setw(10) << "(m/s)" <<          // [m/s]
		" " << setw(10) << "(m/s)" <<         // [m/s]
										 //<< " " << setw(12) << dclk_Rec
										 //<< " " << setw(9) << ZTD             // [m]
										 //<< " " << setw(9) << grdN * 1000       // [mm]
										 //<< " " << setw(9) << grdE * 1000       // [mm]
										 //<< " " << setw(17) << clk //<< (clk/CLIGHT)*1e6             // [m]
		" " << setw(9) << "(m)" <<          // [m]
		" " << setw(9) << "(m)" <<          // [m]
		" " << setw(9) << "(m)" <<       // [m]
		" " << setw(9) << "(m/s)" <<          // [m/s]
		" " << setw(9) << "(m/s)" <<          // [m/s]
		" " << setw(9) << "(m/s)" <<          // [m/s]
		" " << setw(5) << "(#)" <<
		" " << setw(5) << "(#)" <<
		" " << setw(8) << "(m)" <<
		" " << setw(10) << " " <<
		//" " << setw(10) << " " <<
		endl;

	// Print flt results
	if (_flt) {
		_flt->write(os.str().c_str(), os.str().size());
		_flt->flush();
	}
}

void gsins::t_gpvtflt_old::_combineEWL(t_gsatdata & satdata)
{
}

