#include "gpvtflt.h"
#include "gutils/gtimesync.h"
#include "gmodels/gprecisemodel.h"
#include "gcoders/bnccorr.h"
#include "gcoders/bncobs.h"
#include "gproc/gspplsq.h"
#include <algorithm>
//#define DEBUG 1

great::t_gpvtflt::t_gpvtflt(string mark, string mark_base, t_gsetbase* gset, t_gallproc *allproc)
	: t_gspp(mark, gset),
	t_gpppflt(mark, gset),
	_fix_mode(t_gsetamb::FIX_MODE::NO),
	_isFirstFix(true),
	_amb_state(false),
	_site_base(mark_base),
	_gModel_base(nullptr),
	_pos_constrain(false),
	_allproc(allproc),
	_gquality_control(gset, nullptr)
{
	gtrace("t_gpvtflt::t_gpvtflt");
	_vel = t_gtriple(0, 0, 0);  _Qx_vel.ReSize(4); _Qx_vel = 0.0;
	_dclkStoModel = new t_whitenoise(dynamic_cast<t_gsetflt*>(_set)->noise_dclk());
	_velStoModel = new t_whitenoise(dynamic_cast<t_gsetflt*>(_set)->noise_vel());
	_fix_mode = dynamic_cast<t_gsetamb*>(_set)->fix_mode();
	_upd_mode = dynamic_cast<t_gsetamb*>(_set)->upd_mode();
	_max_res_norm = dynamic_cast<t_gsetproc*>(_set)->max_res_norm();
	_sd_sat = dynamic_cast<t_gsetproc*>(_set)->sd_sat();
	_minsat = dynamic_cast<t_gsetproc*>(_set)->minsat();
	_realtime = dynamic_cast<t_gsetproc*>(_set)->realtime();
	_isBase = false;
	if (!_site_base.empty())_isBase = true;
	if (_isBase) _rtkinit();

	if (!_pos_kin)
	{
		xyz_standard = dynamic_cast<t_gsetrec*>(gset)->get_crd_xyz(_site);
		neu_sum_nepo = make_pair(0, t_gtriple(0.0, 0.0, 0.0));
	}
	//freq_index and band_index
	_band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::GPS);
	_band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::GAL);
	_band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::GLO);
	_band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::BDS);
	_band_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(gset)->band_index(gnut::QZS);

	_freq_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::GPS);
	_freq_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::GAL);
	_freq_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::GLO);
	_freq_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::BDS);
	_freq_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(gset)->freq_index(gnut::QZS);

	_gobs = dynamic_cast<t_gallobs*>((*_allproc)[t_gdata::ALLOBS]);
	_gnav = dynamic_cast<t_gallprec*>((*_allproc)[t_gdata::GRP_EPHEM]);
	_gionex_GIM = dynamic_cast<t_gionex*>((*_allproc)[t_gdata::IONEX]);
	_gupd = dynamic_cast<t_gupd*>((*_allproc)[t_gdata::UPD]);
	//_allaug = dynamic_cast<t_gaug*>((*_allproc)[t_gdata::AUG]);

	this->setOTL(dynamic_cast<t_gallotl*>((*_allproc)[t_gdata::ALLOTL]));
	this->setOBJ(dynamic_cast<t_gallobj*>((*_allproc)[t_gdata::ALLOBJ]));
	this->setDCB(dynamic_cast<t_gallbias*>((*_allproc)[t_gdata::ALLBIAS]));
	this->setFCB(dynamic_cast<t_gallbias*>((*_allproc)[t_gdata::ALLBIAS]));

	// init AMB update mode
	_frequency = dynamic_cast<t_gsetproc*>(gset)->frequency();
	_slip_model = dynamic_cast<t_gsetproc*>(gset)->slip_model();
	_turbo_liteMode = dynamic_cast<t_gsetturboedit*>(gset)->liteMode();
	if (_slip_model == SLIPMODEL::TURBO_EDIT)
	{
		_gturboedit = make_shared<t_gexeturboedit>(gset, _log, _allproc);
		if (!_turbo_liteMode)
		{
			_gturboedit->ProcessBatch(_site, _gobs->beg_obs(_site), _gobs->end_obs(_site), _sampling);
			switch (_frequency)
			{
			case 1:
			case 2:
				_slip12 = shared_ptr<t_gcycleslip>(new t_gturboedit(gset, nullptr));
				break;
			case 3:
				_slip12 = shared_ptr<t_gcycleslip>(new t_gturboedit(gset, nullptr));
				_slip13 = shared_ptr<t_gcycleslip>(new t_gturboedit(gset, nullptr, 3));
				break;
			case 4:
			case 5:
				_slip12 = shared_ptr<t_gcycleslip>(new t_gturboedit(gset, nullptr));
				_slip13 = shared_ptr<t_gcycleslip>(new t_gturboedit(gset, nullptr, 3));
				_slip14 = shared_ptr<t_gcycleslip>(new t_gturboedit(gset, nullptr, 4));
				_slip15 = shared_ptr<t_gcycleslip>(new t_gturboedit(gset, nullptr, 5));
				break;
			default:
				throw logic_error("you should set the frequency in xml as 1/2/3/4/5.");
				break;		
			}
		}
	}
	else  // DEFAULT
	{
		_gpre = make_shared<t_gpreproc>(_gobs, gset);
		_gpre->glog(_log);
		_gpre->setNav(_gnav);
	}

	if (_observ == RAW_SINGLE)   // extern ionosphere constraint for single-frequency PPP
	{
		if (_gionex_GIM) _iono_constraint = true;
		else             _iono_constraint = false;

		_epoch_cout = 0;
	}

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

	shared_ptr<t_gbiasmodel> precise_bias(new t_gprecisebias(_allproc, _log, gset));

	if (!_isBase)
	{
		switch (_observ)
		{
		case RAW_ALL:
		case RAW_DOUBLE:
		case RAW_MIX:
			_base_model = shared_ptr<t_gbasemodel>(new t_gcombALL(gset, _log, precise_bias, _allproc));	
			if (!_turbo_liteMode)
			{
				_update_AMB = shared_ptr<t_gupdatepar>(new t_gupdateparALL(_slip12, _slip13, _slip14, _slip15, _band_index, IONO_ORDER::NONE));
			}
			break;
		case IF_P1:	
			_base_model = shared_ptr<t_gbasemodel>(new t_gcombMIX(gset, _log, precise_bias, _allproc));	
			break;
		case IONO_FREE:		
			_base_model = shared_ptr<t_gbasemodel>(new  t_gcombIF(gset, _log, precise_bias, _allproc));	
			if (!_turbo_liteMode)
			{
				if (_frequency == 2) _update_AMB = shared_ptr<t_gupdatepar>(new t_gupdateparIF(_slip12, _band_index));
				if (_frequency == 3) _update_AMB = shared_ptr<t_gupdatepar>(new t_gupdateparIF13(_slip12, _slip13, _band_index));
			}
			break;
		default:			
			_base_model = shared_ptr<t_gbasemodel>(new t_gcombALL(gset, _log, precise_bias, _allproc));	
			if (!_turbo_liteMode)
			{
				_update_AMB = shared_ptr<t_gupdatepar>(new t_gupdateparALL(_slip12, _slip13, _slip14, _slip15, _band_index, IONO_ORDER::NONE));
			}
			break;
		}
	}
	else
	{
		_base_model = shared_ptr<t_gbasemodel>(new t_gcombDD(gset, _log, precise_bias, _allproc));
		dynamic_cast<t_gcombDD*>(&(*_base_model))->set_observ(_observ);
		dynamic_cast<t_gcombDD*>(&(*_base_model))->set_site(_site, _site_base);
	}

	// use to estimate bias parameter(ISB/IFB/..) as white noise
	_isClient = dynamic_cast<t_gsetgen*>(gset)->isClient();
	if(_isClient)	_gInterpol = make_shared<t_gInterpol>(gset, allproc, _log);	//lvhb added in 20201105

	///< real time
	_wl_Upd_time = WL_IDENTIFY;
	_ewl_Upd_time = EWL_IDENTIFY;

}

great::t_gpvtflt::~t_gpvtflt()
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

int great::t_gpvtflt::_addObsD(t_gsatdata& satdata, unsigned int& iobs, t_gallpar& param, t_gtriple& XYZ, Matrix& A, ColumnVector& l, DiagonalMatrix& P)
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

int great::t_gpvtflt::_rtkinit()
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

int great::t_gpvtflt::_setCrd()
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

	if (basepos == BASEPOS::SPP)
	{// from spp
		shared_ptr<t_gspplsq> spplsq = make_shared<t_gspplsq>(_site_base, _set, _gallobj);

		spplsq->setDAT(_gobs, _gnav);
		if (spplsq->processBatch(beg, end) <= 0) {
			if (_log) {
				_log->comment(0, _site_base + "spplsq failed");
			}
			return -1;
		}
		t_gtriple crd = spplsq->getAvergeCrd();
		grec_base->crd(crd, t_gtriple(1, 1, 1), beg, end,true);
	}
	else if (basepos == BASEPOS::CFILE)
	{
		// set the value [read from xml/snx/rnxo]
		auto crd = dynamic_cast<t_gsetrec*>(_set)->get_crd_xyz(_site_base);
		grec_base->crd(crd, t_gtriple(0.01, 0.01, 0.01), beg, end, true);
	}

	//rover
	if (_crd_est == FIX && !_pos_kin) {
		// set the value [read from xml/snx/rnxo]
		auto crd = dynamic_cast<t_gsetrec*>(_set)->get_crd_xyz(_site);
		_grec->crd(crd, t_gtriple(0.01, 0.01, 0.01), beg, end, true);
	}

	return 1;
}

/* detect valid zero-difference residuals*/
bool great::t_gpvtflt::_valid_residual(bool phase_process, string sat_name, enum FREQ_SEQ &f, map<string, map<FREQ_SEQ, pair<int, int>>>& index_l) {
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

int great::t_gpvtflt::_combineDD(Matrix & A, SymmetricMatrix & P, ColumnVector & l)
{
	gtrace("t_gpvtflt_newmodelmodel::_combineDD");
	try
	{
		// 生成index_l todo : sat,freq,pair<P,L>
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
			else if (obstype == GOBSTYPE::TYPE_L) index_l[sat][f].second = i + 1;
		}

		_obs_index.clear();
		int iobs = 1;

		set<string> sysall = dynamic_cast<t_gsetgen*>(_set)->sys();

		int nrows = A.Nrows();
		Matrix DD; DD.ReSize(nrows, nrows); DD = 0.0;
		_sat_ref.clear();
		bool isSetRefSat = dynamic_cast<t_gsetamb*>(_set)->isSetRefSat();

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
				if(f > _frequency) continue;// modified by lvhb in 20201211
				string sat, sat_ref;
				t_gsatdata obs_sat_ref; enum GSYS gs;
				for (auto it = _data.begin(); it != _data.end(); it++)
				{
					string sat = it->sat();

					gs = it->gsys();
					if (gs == QZS) gs = GPS;
					if (gs != sys) continue;
					if (gs == BDS && t_gsys::bds_geo(sat)) continue;
					if (it->islip()) continue;

					/* select reference satellite,lvhb modified in 20200826 */
					if (isSetRefSat && _observ == RAW_ALL) {
						FREQ_SEQ f1 = FREQ_1;
						FREQ_SEQ f2 = FREQ_2;
						if(_valid_residual(isPhaseProcess,sat,f1,index_l) == false
							|| _valid_residual(isPhaseProcess, sat, f2, index_l) == false
							) continue;
					}
					else if (_valid_residual(isPhaseProcess, sat, f, index_l) == false)continue;

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
				if ( _ipSatRep[sys] != "" && sat_ref != _ipSatRep[sys])
				{
					cerr << "refsat bug" << endl;
					continue;
				}
				if (f == FREQ_1) _sat_ref.insert(sat_ref);
			
				if (_log) _log->comment(2, "gpvtflt", _site + _epoch.str_ymdhms(" epoch ")
					+ sat_ref + ": is reference satellite!");

				int index_ref;
				if (isPhaseProcess)	index_ref = index_l[sat_ref][f].second;
				else index_ref = index_l[sat_ref][f].first;
				/* make double difference */
				for (auto it = _data.begin(); it != _data.end(); it++)
				{
					sat = it->sat();
					if (sat == sat_ref)	continue;

					gs = it->gsys();
					if (gs == QZS) gs = GPS;
					if (gs != sys) continue;

					if (_valid_residual(isPhaseProcess, sat, f, index_l) == false) continue;
					int index_sat;
					if (isPhaseProcess)	index_sat = index_l[sat][f].second;
					else index_sat = index_l[sat][f].first;
					
					/* lvhb added in 202020815. Actually it should divided in two threshold : phase and code,
					 maybe it can be optimized in future.*/
					//if (fabs(l(index_sat) - l(index_ref) > 30)) continue;

					DD(iobs, index_ref) = -1;
					DD(iobs, index_sat) = 1;
					GOBSTYPE obstype = TYPE_C;
					if (isPhaseProcess)obstype = GOBSTYPE:: TYPE_L;
					_obs_index.push_back(make_pair(it->sat(), make_pair(f, obstype)));
					iobs++;
				}//end sat
			}//end f
		}//end sys

		// delete zero rows/cols
		iobs--;
		if (iobs <= 0) return -1;

		DD = DD.Rows(1, iobs);

		A = DD * A;
		l = DD * l;
		P << (DD * P.i() * DD.t()).i();

#if DEBUG
		 cout << _epoch.str_ymdhms()<< endl;
		 t_out("ddA", A);
		 t_out("ddL", l);
		 t_out("ddpr", P);
		 t_out("sdQ", _Qx);
#endif

		return 1;
	}
	catch (...)
	{
		if (_log) _log->comment(1, "t_gpvtflt_newmodelmodel", "generate double-differented equation failed!");
		return -1;
	}
}

int great::t_gpvtflt::_postRes(const Matrix& A,const SymmetricMatrix& P,const ColumnVector& l, Matrix& pA, SymmetricMatrix& pP, ColumnVector& pl)
{
	if (_sat_ref.empty())	return -1;
	if (!_amb_state) { _post_A = A; _post_P = P; _post_l = l; return 0; }

	// index_l todo : sat,freq,pair<P,L>
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
		else if (obstype == GOBSTYPE::TYPE_L) index_l[sat][f].second = i + 1;
	}

	int nobs = P.Nrows(); int ncols = A.Ncols();
	pA.ReSize(nobs, A.Ncols()); pA = A;
	pP = P;
	pl.ReSize(nobs); pl = l;
	DiagonalMatrix delP; delP.ReSize(nobs); delP = 0;

	t_DD_amb DD = _ambfix->getDD();
	int iobs = 0;
	vector<int> ind;
	for (auto it = 1; it <= nobs; it++)ind.push_back(it);

	for (auto itref = _sat_ref.begin(); itref != _sat_ref.end(); itref++)
	{
		string satref = *itref;
		for (auto itdd = DD.begin(); itdd != DD.end(); itdd++)
		{
			if (itdd->isWlFixed && itdd->isNlFixed)
			{
				string sat1 = get<0>(itdd->ddSats[0]);
				string sat2 = get<0>(itdd->ddSats[1]);

				string crt_sat = sat1;
				if (crt_sat.substr(0, 1) != satref.substr(0, 1))continue;
				if (sat1 != satref && sat2 != satref) continue;
				if (sat1 == satref && sat2 != satref)
				{
					crt_sat = sat2;
					itdd->rlc *= -1.0;
					itdd->rwl *= -1.0;
					itdd->inl *= -1.0;
				}
				FREQ_SEQ f = FREQ_X; double dif = 0.0;
				if (_observ == OBSCOMBIN::IONO_FREE)
				{
					t_gobsgnss gnss(satref);
					auto gs = gnss.gsys();
					double lambda_1 = gnss.wavelength(_band_index[gs][FREQ_1]);
					double lambda_2 = gnss.wavelength(_band_index[gs][FREQ_2]);

					double wl = round(itdd->rwl);
					double mw_coef = lambda_2 * SQR(lambda_1)/(SQR(lambda_2)-SQR(lambda_1));
					dif = itdd->rlc - (itdd->inl * itdd->factor + wl * mw_coef);
					f = FREQ_1;
				}
				else if (_observ == OBSCOMBIN::RAW_ALL)
				{
					if(itdd->ambtype == "AMB_L1") f = FREQ_1;
					else if (itdd->ambtype == "AMB_L2") f = FREQ_2;
					else if (itdd->ambtype == "AMB_L3") f = FREQ_3;
					else if (itdd->ambtype == "AMB_L4") f = FREQ_4;
					else if (itdd->ambtype == "AMB_L5") f = FREQ_5;
					else continue;

					dif = itdd->rlc - itdd->inl * itdd->factor;
				}
				if (f != FREQ_X)
				{
					auto index = index_l[crt_sat][f].second;
					pl(index) -= dif;
					vector<int>::const_iterator it = find(ind.begin(), ind.end(), index);
					//if (it != ind.end())ind.erase(it);
				}
			}
		}//end dd
	}//end satref

	vector<int> ind1 = ind, ind2 = ind;
	int pobs = nobs - ind.size();
	Matrix_rem(pP, ind);
	Matrix_remR(pA, ind1);
	Matrix_remR(pl, ind2);
	pA.SubMatrix(1, pobs, 4, ncols) = 0.0;
	for (int i = 0; i < pobs; i++)pP(i + 1, i + 1) *= 10;
	//t_out("A:", A);
	//t_out("pA:", pA);
	//t_out("P:", P);
	//t_out("pP:", pP);
	//t_out("l:", l.t());
	//t_out("pl:", pl.t());

	return iobs;
}

int great::t_gpvtflt::_prepareData()
{
	gtrace("t_gpvtflt::_prepareData");
	int flag = 1;
	if (_isBase) {
		if (_prepareData(_site_base, _data_base) < 0) {
			if (_log) _log->comment(2, _site + "gpvtflt", "Preparing data failed!");
			flag = -1;
		}
		_vBanc_base = _vBanc;
	}
	if (_prepareData(_site, _data) < 0) {
		if (_log) _log->comment(2, _site + "gpvtflt", "Preparing data failed!");
		flag = -1;
	}
	if (flag == -1) {
		_nSat = 0;
		return -1;
	}

	if (_cntrep == 0){
		t_gtriple xyz = t_gtriple(_vBanc.Rows(1, 3)); // modified by zhshen (linux)
		_gquality_control.processOneEpoch(_epoch, _site, xyz, _data);
		if (_isBase) {
			xyz = t_gtriple(_vBanc_base.Rows(1, 3));
			_gquality_control.processOneEpoch(_epoch, _site_base, xyz, _data_base);
		}
	}

	if (_isBase) {
		if ((_nSat = _selcomsat(_data_base, _data)) < _minsat) {
			if (_log) _log->comment(2, "gpvtflt", "no common satellite");
			flag = -1;
		}
	}
	if (_data.size() < _minsat) {
		if (_log) _log->comment(2, "", _site + _epoch.str_ymdhms(" epoch ") + int2str(_data.size()) + " skipped (data.size < minsat)");
		return -1;
	}
	return flag;
}

int great::t_gpvtflt::_prepareData(const string & ssite, vector<t_gsatdata>& sdata)
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
	set<string> sat_rm = dynamic_cast<t_gsetgen*>(_set)->sat_rm();
	while (iter != sdata.end()) {

		// except sat from config file
		string satname = iter->sat();
		if (sat_rm.find(satname) != sat_rm.end()) {
			iter = sdata.erase(iter);
			continue;
		}

		GSYS gs = iter->gsys();

		//GOBSBAND b1, b2;
		GOBSBAND b1 = _band_index[gs][FREQ_1];
		GOBSBAND b2 = _band_index[gs][FREQ_2];

		iter->glog(_log);

		// check data availability lvhb added in 20200525
		auto l1 = iter->select_phase(b1);
		auto l2 = iter->select_phase(b2);

		auto p1 = iter->select_range(b1);
		auto p2 = iter->select_range(b2);

		// check data availability
		if (p1 == X || l1 == X) {
			iter = sdata.erase(iter);
			continue;
		}

		double P3, L3;
		if (_observ == RAW_SINGLE){
			P3 = iter->obs_C(p1);
			L3 = iter->obs_L(l1);
		}
		else if (_observ == IF_P1 || _observ == RAW_MIX){
			if (p2 == X || l2 == X) {
				P3 = iter->obs_C(p1);
				L3 = iter->obs_L(l1);
			}
			else{
				P3 = iter->P3(p1, p2);
				L3 = iter->L3(l1, l2);
			}
		}
		else{
			if (p2 == X || l2 == X) {
				iter = sdata.erase(iter);
				continue;
			}
			P3 = iter->P3(p1, p2);
			L3 = iter->L3(l1, l2);
		}

		if (double_eq(L3, 0.0) && _phase) {
			if (_log) {
				ostringstream str;
				str << "prepareData: erasing data due to no phase double bands observation, " <<
					"epo: " << _epoch.str_hms() << ", " <<
					"prn: " << iter->sat() << " (" << iter->channel() << ")";
				_log->comment(2, "gpvtflt", str.str());
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
				_log->comment(2, "gpvtflt", str.str());
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
				_log->comment(2, "gpvtflt", str.str());
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

			iter++;
			continue;
		}
	}//end sdata

	if (sdata.size() < _minsat) {
		if (_log)
			_log->comment(2, "gpvtflt", _epoch.str_ymdhms(ssite + " epoch ")
				+ " skipped (Bancroft not calculated: " + int2str(sdata.size()) + " < _minsat)");
		//      cout << t_gsys::gsys2str(_gnss) << " " << _epoch.str_ymdhms( _site + " epoch ") << " " << int2str(_data.size()) <<  " skipped (Bancroft not calculated: _data.size < _minsat)" << endl;
		return -1;
	}

	BB = BB.Rows(1, iobs);    // delete zero rows

	if (BB.Nrows() < static_cast<int>(_minsat)) {
		if (_log && _log->verb() >= 2)
			_log->comment(2, "gpvtflt", _epoch.str_ymdhms(_site + " epoch ")
				+ " skipped (Bancroft not calculated: BB.Nrows < _minsat)");
		//      cout << _epoch.str_ymdhms( _site + " epoch ") << int2str(_data.size()) <<  " skipped (Bancroft not calculated: _data.size < _minsat)" << endl;
		return -1;
	}

	_vBanc = 0.0;

	if (BB.Nrows() >= 4)
		gbancroft(BB, _vBanc);
	//if(_isClient)	cout << "sppflt 2  :" << fixed << fixed << _vBanc(1) << "  " << _vBanc(2) << "  " << _vBanc(3) << "  " <<ssite << endl;
	shared_ptr<t_gobj> grec = _gallobj->obj(ssite);
	if (grec == nullptr) return -1;

	if (!_initialized && (_vBanc.Rows(1, 3) - grec->crd_arp(_epoch).crd_cvect()).NormFrobenius() > 1000) {
		if (ssite == _site) _valid_crd_xml = false;
		else {//base setting, no set _valid_crd_xml value,lvhb modified in 20200925
			_vBanc(1) = grec->crd_arp(_epoch)[0];
			_vBanc(2) = grec->crd_arp(_epoch)[1];
			_vBanc(3) = grec->crd_arp(_epoch)[2];
		}
	}

	if (_valid_crd_xml && !_initialized) {
		_vBanc(1) = grec->crd_arp(_epoch)[0];
		_vBanc(2) = grec->crd_arp(_epoch)[1];
		_vBanc(3) = grec->crd_arp(_epoch)[2];
		//if (_isClient)	cout << "sppflt 2.1  :" << fixed << fixed << _vBanc(1) << "  " << _vBanc(2) << "  " << _vBanc(3) << "  " << ssite << endl;
	}

	//lvhb: just to judge spp/bancroft coordinates
	t_gtriple test_xyz(_vBanc);
	t_gtriple test_ell;
	xyz2ell(test_xyz, test_ell, true);
	double radius = test_xyz.norm();
	if (radius < B_WGS - 500) {
		string warning = "WARNING: Unexpected site (" + ssite + ") coordinates from Bancroft. Orbits/clocks or code observations should be checked.";
		if (_log)  _log->comment(1, "gpvtflt", warning + _epoch.str_ymdhms(" Epoch "));
		if (_glog) _glog->comment(2, "gpvtflt", warning + _epoch.str_ymdhms(" Epoch "));
		if (!_phase) { return -1; }
		return -1;
	}

	if (ssite == _site && _crd_est != FIX ) {
		/********* spplsq init **********/
		//t_gtriple crd(0, 0, 0); double clk = 0;
		//if (_get_spplsq_result(_site, _epoch, crd, clk) > 0)
		//{
		//	_vBanc(1) = crd[0];
		//	_vBanc(2) = crd[1];
		//	_vBanc(3) = crd[2];
		//	_vBanc(4) = clk;
		//}
		if (_pos_constrain)
		{
			_vBanc(1) = _extn_pos[0];
			_vBanc(2) = _extn_pos[1];
			_vBanc(3) = _extn_pos[2];
			int iclk = _param.getParam(_site, t_gpar::CLK, "");
			if (BB.Nrows() < 4)_vBanc(4) = _param[iclk].value();
		}
	}

	//Compute sat elevation and rho
	iter = sdata.begin();
	while (iter != sdata.end()) {
		t_gtriple xyz_r, xyz_s;
		xyz_s = iter->satcrd();

		if (_isBase && ssite == _site_base) {//lvhb added in 20200523
			xyz_r = grec->crd_arp(_epoch);
			if (xyz_r.norm() == 0) {
				cerr << endl;
				cerr << xyz_r[0] << " " << xyz_r[1] << " " << xyz_r[2] << endl;
			}
			_vBanc(1) = xyz_r[0];
			_vBanc(2) = xyz_r[1];
			_vBanc(3) = xyz_r[2];
		}
		else {
			if (_crd_est == FIX && _valid_crd_xml) {
				xyz_r = grec->crd_arp(_epoch);
				//lvhb added in 202008
				_vBanc(1) = xyz_r[0];
				_vBanc(2) = xyz_r[1];
				_vBanc(3) = xyz_r[2];
			}
			else if (_crd_est == EST || _crd_est == KIN || _crd_est == SIMU_KIN) {  //用于适配
				if (!_initialized) {
					xyz_r.set(_vBanc.Rows(1, 3));
				}
				else {
					if (_pos_kin || _crd_est == KIN || _crd_est == SIMU_KIN) {
						xyz_r.set(_vBanc.Rows(1, 3));
					}
					else _param.getCrdParam(ssite, xyz_r);
				}
				//modified by lvhb,20200925
				if (_pos_kin){
					t_gtriple crd_r = xyz_r - grec->eccxyz(_epoch);
					grec->crd(crd_r, t_gtriple(1, 1, 1), _epoch, LAST_TIME, true);
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
		
		// Changed by LX: add for PPP-RTK
		_crt_ele[iter->sat()] = iter->ele_deg();

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

void great::t_gpvtflt::_predict()
{
	gtrace("t_gpvtflt::_predict");

	// Predict coordinates, clock and troposphere
	// state vector is the same, covariance matrix is predicted with noise   
	int i;
	double bl = 0.0;

	_cntrep++;

	// add/remove inter system bias
	if (!_isBase) _syncSys();
	
	// add/remove ionosphere delay
	//Actually,we only estimate sion/sdsion by using white noise here.	lvhb 202007
	_syncIono();
	if (_isBase && 1)//lvhb 20200822
	{
		if ( _iono_est) {
			auto crd = _gallobj->obj(_site_base)->crd_arp(_epoch);
			t_gtriple tmpcrd(_vBanc);
			t_gtriple tmpell, tmpdxyz, tmpneu;
			xyz2ell(crd, tmpell, false);
			tmpdxyz = tmpcrd - crd;
			xyz2neu(tmpell, tmpdxyz, tmpneu);
			bl = tmpneu.norm();

			if (bl < 10000 && _pos_kin) {
				_delPar(t_gpar::VION);
				_delPar(t_gpar::SION);
			}
		}
	}

	// add/remove inter-frequency biases
	// LX changed IFB estimation
	//if (!_isBase) _syncIFB();
	if (!_isBase&&_frequency >=3) _syncIFB();

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
	//cout << "Pridavam CRDX " <<  " " << _param[i].value() << " " << _Qx(i + 1, i + 1) << endl;

	i = _param.getParam(_site, t_gpar::CRD_Y, "");
	if (i >= 0)
	{
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0 )
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
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0 )
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
	//cout << "Pridavam CLK " << " " << _param[i].value() << " " << _Qx(i + 1, i + 1) << endl;

	i = _param.getParam(_site, t_gpar::GLO_ISB, "");
	if (i >= 0)
	{
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0||_isClient)
		{
			_param[i].value(0.0);
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
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0||_isClient) {
			_param[i].value(0.0);
			_Qx(i + 1, i + 1) = _sig_init_glo * _sig_init_glo;
		}
		else {
			if (_cntrep == 1) _Qx(i + 1, i + 1) += _gloStoModel->getQ();
			if (_smooth) _Noise(i + 1, i + 1) = _gloStoModel->getQ();
		}
	}

	i = _param.getParam(_site, t_gpar::GAL_ISB, "");
	if (i >= 0) {
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0||_isClient) {
			_param[i].value(0.0);
			_Qx(i + 1, i + 1) = _sig_init_gal * _sig_init_gal;
		}
		else {
			if (_cntrep == 1) _Qx(i + 1, i + 1) += _galStoModel->getQ();
			if (_smooth) _Noise(i + 1, i + 1) = _galStoModel->getQ();
		}
	}

	i = _param.getParam(_site, t_gpar::BDS_ISB, "");
	if (i >= 0) {
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0||_isClient) {
			_param[i].value(0.0);
			_Qx(i + 1, i + 1) = _sig_init_bds * _sig_init_bds;
		}
		else {
			if (_cntrep == 1) _Qx(i + 1, i + 1) += _bdsStoModel->getQ();
			if (_smooth) _Noise(i + 1, i + 1) = _bdsStoModel->getQ();
		}
	}

	i = _param.getParam(_site, t_gpar::QZS_ISB, "");
	if (i >= 0) {
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0||_isClient) {
			_param[i].value(0.0);
			_Qx(i + 1, i + 1) = _sig_init_qzs * _sig_init_qzs;
		}
		else {
			if (_cntrep == 1) _Qx(i + 1, i + 1) += _qzsStoModel->getQ();
			if (_smooth) _Noise(i + 1, i + 1) = _qzsStoModel->getQ();
		}
	}

	// predict ionosphere delay  
	if (_iono_est)
	{
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
					_Qx(i + 1, i + 1) += _ionStoModel->getQ();// *_ionStoModel->getQ();
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
				// for (unsigned int jj = 1; jj <= _param.parNumber(); jj++) _Qx(i + 1, jj) = 0.0;
				double var = _sig_init_vion * _sig_init_vion;
				if (_isBase && double_eq(_Qx(i + 1, i + 1), var)) {
					var *= SQR(bl / (1e4));
					_Qx(i + 1, i + 1) = var;
					_param[i].value(1e-6);
				}
				else if (_isClient && !_isBase){
					_param[i].value(0.0);
					_Qx(i + 1, i + 1) = var;
				}

				if (_cntrep == 1 && !double_eq(_Qx(i + 1, i + 1), var)) {
					
					if (_isBase) {
						
						_Qx(i + 1, i + 1) += SQR(bl / (1e4) * cos(it->ele()))* _ionStoModel->getQ();
					}
					else {
						_Qx(i + 1, i + 1) += _ionStoModel->getQ();
					}
				}
			}
		}

	}

	i = _param.getParam(_site, t_gpar::IFB_GPS, "");
	if (i >= 0) {
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0||_isClient) {
			_param[i].value(0.0);
			_Qx(i + 1, i + 1) = 3000 * 3000;
		}
		else {
			if (_cntrep == 1) _Qx(i + 1, i + 1) += _gpsStoModel->getQ();
			if (_smooth) _Noise(i + 1, i + 1) = _gpsStoModel->getQ();
		}
	}
	i = _param.getParam(_site, t_gpar::IFB_GAL, "");
	if (i >= 0) {
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0||_isClient) {
			_param[i].value(0.0);
			_Qx(i + 1, i + 1) = 3000 * 3000;
		}
		else {
			if (_cntrep == 1) _Qx(i + 1, i + 1) += _galStoModel->getQ();
			if (_smooth) _Noise(i + 1, i + 1) = _galStoModel->getQ();
		}
	}
	i = _param.getParam(_site, t_gpar::IFB_BDS, "");
	if (i >= 0) {
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0||_isClient) {
			_param[i].value(0.0);
			_Qx(i + 1, i + 1) = 3000 * 3000;
		}
		else {
			if (_cntrep == 1) _Qx(i + 1, i + 1) += _bdsStoModel->getQ();
			if (_smooth) _Noise(i + 1, i + 1) = _bdsStoModel->getQ();
		}
	}

	i = _param.getParam(_site, t_gpar::IFB_QZS, "");
	if (i >= 0) {
		if (!_initialized || _Qx(i + 1, i + 1) == 0.0 || _isClient) {
			_Qx(i + 1, i + 1) = 3000 * 3000;
		}
		else {
			if (_cntrep == 1) _Qx(i + 1, i + 1) += _qzsStoModel->getQ();
			if (_smooth) _Noise(i + 1, i + 1) = _qzsStoModel->getQ();
		}
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
		else
		{
			if(!_turbo_liteMode && _slip_model == SLIPMODEL::TURBO_EDIT) _post_turbo_syncAmb();
			else                                                         _syncAmb();
		}
	}

	// ambiguity randomwalk
	if (_slip_model != SLIPMODEL::TURBO_EDIT || (_slip_model == SLIPMODEL::TURBO_EDIT && _turbo_liteMode))
	{
		for (unsigned int i = 0; i < _param.parNumber(); i++) {
			if (_param[i].parType == t_gpar::AMB_IF ||
				_param[i].parType == t_gpar::AMB_L1 ||
				_param[i].parType == t_gpar::AMB_L2 ||
				_param[i].parType == t_gpar::AMB_L3 ||
				_param[i].parType == t_gpar::AMB_L4 ||
				_param[i].parType == t_gpar::AMB_L5) {

				if (_cntrep == 1) _Qx(i + 1, i + 1) += _ambStoModel->getQ();
				if (_smooth) _Noise(i + 1, i + 1) = _ambStoModel->getQ();
				_param[i].stime = _param[i].end = _epoch;
			}
		}
	}
	//tropspheric parameter,modified by lvhb 20200917
	if (_tropo_est) {
		//rover: sppflt::predict()
		//base:
		double ztdInit = _sig_init_ztd;
		for (int j = 0;j<1+(_isBase?1:0);j++)
		{
			auto tmpsite = _site;
			auto tmpgrec = _grec;
			auto tmpgmodel = _gModel;

			if (j)
			{
				tmpsite = _site_base;
				tmpgrec = _gallobj->obj(_site_base);
				tmpgmodel = _gModel_base;
			}
			i = _param.getParam(tmpsite, t_gpar::TRP, "");
			if (i >= 0) {
				t_gtriple Ell, XYZ;

				if (_param.getCrdParam(tmpsite, XYZ) > 0) {
				}
				else { XYZ = tmpgrec->crd_arp(_epoch); }
				xyz2ell(XYZ, Ell, false);
				if (tmpgmodel->tropoModel() != 0) {
					_param[i].apriori(_sig_init_ztd);
				}

				if (!_initialized || _Qx(i + 1, i + 1) == 0.0||_isClient) {
					if (_valid_ztd_xml) {
						_param[i].value(_aprox_ztd_xml - tmpgmodel->tropoModel()->getZHD(Ell, _epoch));
					}
					else {
						if (tmpgmodel->tropoModel() != 0) {
							_param[i].value(tmpgmodel->tropoModel()->getZWD(Ell, _epoch));
						}
					}
					if (_isClient) _param[i].value(0.0);
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

void great::t_gpvtflt::_delPar(const t_gpar::t_type par)
{
	// Remove params and appropriate rows/columns covar. matrix
	for (unsigned int i = 0; i <= _param.parNumber() - 1; i++)
	{
		if (_param[i].parType == par)
		{
			Matrix_remRC(_Qx, _param[i].index, _param[i].index);
			_param.delParam(i);
			_param.reIndex();
			i--;
		}
	}

	return;
}

int great::t_gpvtflt::_combineMW(t_gsatdata& satdata)
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
	bool islip = (satdata.getlli(gobsL1.gobs()) >= 1 || satdata.getlli(gobsL2.gobs()) >= 1);

	//dcb correct  已经在观测之上修改了dcb
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

			if (double_eq(mw_obs, 0.0) || double_eq(mw_obs_base, 0.0)) mw_obs = 0.0;
			else mw_obs -= mw_obs_base;

			if (!islip) islip = (satdata_base.getlli(gobsL1.gobs()) >= 1 || satdata_base.getlli(gobsL2.gobs()) >= 1);
		}

	}

	if (!double_eq(mw_obs, 0.0))
	{
		// calculate smooth MW
		t_gtime pre_time = crt_time - obs_intv;
		if (_MW.find(pre_time) != _MW.end() &&_MW[pre_time].find(sat) != _MW[pre_time].end()
			&& !double_eq(_MW[pre_time][sat][4], 0.0)&& !islip)
		{
			_MW[crt_time][sat][1] = _MW[pre_time][sat][1] + 1;
			_MW[crt_time][sat][4] = mw_obs;
			if (!double_eq(_MW[crt_time][sat][1], 1.0) && _MW[crt_time][sat][1] > 1.0)
			{
				_MW[crt_time][sat][3] = _MW[pre_time][sat][3] * (_MW[crt_time][sat][1] - 2) /
					(_MW[crt_time][sat][1] - 1) + pow(mw_obs - _MW[pre_time][sat][2], 2)
					/ _MW[crt_time][sat][1];
			}
			else
			{
				_MW[crt_time][sat][3] = _MW[pre_time][sat][3];
			}
			_MW[crt_time][sat][2] = _MW[pre_time][sat][2] + (mw_obs - _MW[pre_time][sat][2])
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

bool great::t_gpvtflt::_external_pos(const t_gtriple& xyz_r, const t_gtriple& rms)
{
	_extn_pos = xyz_r;
	_extn_rms = rms;
	_pos_constrain = true;
	return true;
}

bool great::t_gpvtflt::_pos_virtual_obs(const t_gtriple & xyz_r, const t_gtriple & rms, Matrix & A, ColumnVector & l, SymmetricMatrix & P)
{
	try
	{
		int obs_num = A.Nrows(), par_num = _param.parNumber();
		int virtual_obs_num = 3;
		A.resize_keep(obs_num + virtual_obs_num, par_num);
		P.resize_keep(obs_num + virtual_obs_num);
		l.resize_keep(obs_num + virtual_obs_num);
		int icrdx = _param.getParam(_site, t_gpar::CRD_X, "");
		int icrdy = _param.getParam(_site, t_gpar::CRD_Y, "");
		int icrdz = _param.getParam(_site, t_gpar::CRD_Z, "");
		A(obs_num + 1, icrdx + 1) = 1; P(obs_num + 1, obs_num + 1) = SQR(rms[0]); l(obs_num + 1) = 0;
		A(obs_num + 2, icrdy + 1) = 1; P(obs_num + 2, obs_num + 2) = SQR(rms[1]); l(obs_num + 2) = 0;
		A(obs_num + 3, icrdz + 1) = 1; P(obs_num + 3, obs_num + 3) = SQR(rms[2]); l(obs_num + 3) = 0;
		return true;
	}
	catch (...)
	{
		if (_log) _log->comment(1, "t_gpvtflt", "Add position virtual obs failed!");
		return false;
	}
}

int great::t_gpvtflt::_processEpochVel()
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
	} while (threshold > 1e-2);

	if (param_vel.getVelParam(_site, _vel) < 0)return -1;

	return 1;
}

int great::t_gpvtflt::_processEpoch(const t_gtime & runEpoch)
{
	gtrace("t_gpvtflt::_processEpoch");

#ifdef DEBUG
	cout << "gpvtflt: Epoch - full data size: " << _data.size() << " "
		<< runEpoch.str_hms() << endl;
#endif

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
	ColumnVector v_orig, v_norm, v_test;
	SymmetricMatrix Qsav, QsavBP;
	vector<t_gsatdata>::iterator it;
	t_gallpar XsavBP;
	double vtpv;
	int nobs_total, npar_number;
	string outlier = "";
	t_gtime begT;
	_cntrep = 0;    // number of iterations caused by outliers

	// select obs for tb log (rover)
	if (_slip_model == SLIPMODEL::TURBO_EDIT && !_turbo_liteMode){
		if (!this->_post_turbo_select_obs(runEpoch, _data)) return -1;
	}

	do {
		// erase sat satellite because of outliers
		it = _data.begin();
		if (outlier != ""){		//lvhb modified
			while (it != _data.end()) {
				if (it->sat() == outlier) {
					it = _data.erase(it); continue;
				}
				it++;
			}
		}

		if (_prepareData() < 0) {
			if (_initialized){
				_predict();//todo:add time gap
	//			_reset_param();
	//			_reset_ambig();b
			}
			return -1;
		}
		
		QsavBP = _Qx; XsavBP = _param;
		_predict();
		//_prt_info(_epoch);

		_initialized = true;

		if (_data.size() < _minsat) {
			if (_log) _log->comment(2, "gpvtflt", "Not enough visible satellites!");
			_restore(QsavBP, XsavBP);
			return -1;
		}

		if (_pos_constrain) {
			int icrdx = _param.getParam(_site, t_gpar::CRD_X, "");
			for (int i = 0; i < 3; i++)
			{
				_Qx(_param[icrdx].index + i, _param[icrdx].index + i) = SQR(_extn_rms[i]);
			}
		}

		if (_observ == RAW_SINGLE && _iono_constraint && _cntrep == 1){
			_extern_ion_constraint(runEpoch);
		}
		
		// define a number of measurements 
		unsigned int nObs = _data.size();
		unsigned int mult = 1;
		if (_observ == RAW_DOUBLE) { mult = 2; nObs *= 2; }
		if (_observ == RAW_ALL || _observ == RAW_MIX) { mult = 2; nObs *= 5; } // reservation for 5 freq - not used raws will be removed
		if (_phase) { mult *= 2; nObs *= 2; }  // code + phase
		if (_observ == IF_P1) { mult = 2; nObs *= 3; }  // IF + P1(不增加参数，多一倍观测值）
		// if (_doppler) { mult; nObs *= 2; }
		// nObs *= 2;    

		/*G01 P1 P2 P3... L1 L2L3*/
		unsigned int nPar = _param.parNumber();
		unsigned int iobs = 1;

		//_print_ele();

		Matrix A_new;
		SymmetricMatrix P_new;
		ColumnVector l_new, dx_new;
		dx_new.ReSize(nPar);    dx_new = 0.0;
		unsigned int iobs_new = 1;
		_frqNum.clear(); _obs_index.clear();
		if (_isBase){
			dynamic_cast<t_gcombDD*>(&(*_base_model))->set_base_data(&_data_base);
			dynamic_cast<t_gcombDD*>(&(*_base_model))->set_rec_info(_gallobj->obj(_site_base)->crd_arp(_epoch), _vBanc(4), _vBanc_base(4));
			//dynamic_cast<t_gcombDD*>(&(*_base_model))->set_rec_info(_gallobj->obj(_site_base)->crd_arp(_epoch), 0, 0);
		}
		// use combmodel
		t_glsqEquationMatrix equ;
		it = _data.begin();
		for (it = _data.begin(); it != _data.end(); ) {
			//cout << fixed << setprecision(5) << "before  " << it->sat() << "  crd: " << it->satcrd().crd_cvect().t() << "   vel: " << it->satvel().crd_cvect().t() << "   Send time:  " << it->satTime().str_hms() << endl;
			int irc = _base_model->cmb_equ(_epoch, _param, *it, equ);
			if (irc < 0) { it = _data.erase(it); continue; }
			if(_observ==IONO_FREE) _combineMW(*it);
			//cout << fixed << setprecision(5) << "after  " << it->sat() << "  crd: " << it->satcrd().crd_cvect().t() << "   vel: " << it->satvel().crd_cvect().t() << "   Send time:  " << it->satTime().str_hms() << endl;
			++it;
		}
		iobs_new = equ.num_equ();
		// tranform equ into A, P, l
		_equationToMatrix(equ, A_new, P_new, l_new, nPar);

		// generate obs_index
		_obs_index.clear();
		_generateObsIndex(equ);

		if (iobs_new < _minsat * mult) {
			if (_log) _log->comment(2, "gpvtflt", "Not enough processable observations!");
			_restore(QsavBP, XsavBP);
			return -1;
		}

		A = A_new;
		P = P_new;
		l = l_new;
		dx = dx_new;
		//t_out("l", l); t_out("qx", _Qx);
		if (_isBase || _sd_sat){
			if (_combineDD(A, P, l) < 0) return -1;
		}
		
		Qsav = _Qx;
		try
		{
			_filter->update(A, P, l, dx, _Qx);
		}
		catch (...)
		{
			if (_log) _log->comment(1, "t_gpvtflt", " filter update failed!");
			_Qx = Qsav;
			return -1;
		}

		// increasing variance after update in case of introducing new ambiguity
		if (_cntrep == 1 && !_reset_amb && !_reset_par) {
			for (size_t iPar = 0; iPar < _param.parNumber(); iPar++) {
				//if (_param[iPar].parType == t_gpar::AMB_IF) {
				if (_param[iPar].str_type().find("AMB") != string::npos) {
					string sat = _param[iPar].prn;
					if (_newAMB.find(sat) != _newAMB.end() /*&& _cntrep == 1*/) {
						if (_newAMB[sat] == 1) _Qx(iPar + 1, iPar + 1) += 10;
						if (_newAMB[sat] == 2 && _Qx(iPar + 1, iPar + 1) > 0.01) _Qx(iPar + 1, iPar + 1) += 1;
						//_newAMB[sat]++;
					}
				}
			}
			auto amb_prn = _param.amb_prns();
			for (auto& sat : amb_prn) {
				if (_newAMB.find(sat) != _newAMB.end()) _newAMB[sat]++;
			}
		}
		//t_out("increasing qx:", _Qx);
		ColumnVector dxxx = dx;
		for (int i = 0; i < 3; i++)dxxx(i + 1) = 0;
		// post-fit residuals
		v_orig = l - A * dx;
		v_test = l - A*dxxx;

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

	_filter->add_data(_param, dx, _Qx, _sig_unit, Qsav);
	_filter->add_data(A, P, l);
	_filter->add_data(vtpv, nobs_total, npar_number);
	
	t_gallpar param_after = _param;

	//t_out("B:", A);
	//t_out("P:", P);
	//t_out("l:", l.t());

	_amb_resolution();
	for (unsigned int iPar = 0; iPar < _param.parNumber(); iPar++) {
			_param[iPar].value(_param[iPar].value() + dx(_param[iPar].index));
			if(_isClient && !_isBase) _param[iPar].amb_ini = !_amb_state;
	}

	return _amb_state?1:0;
}

int great::t_gpvtflt::_amb_resolution()
{
	gtrace("t_gpvtflt::_amb_resolution");

	ColumnVector dx_tmp = _filter->dx();
	SymmetricMatrix   Qx_tmp = _filter->Qx();
	_param_fixed = _filter->param();
	_amb_state = false;
#if DEBUG
	cout << _epoch.str_ymdhms() << endl;
	for (int i = 0; i < _param_fixed.parNumber(); i++) 
	{
		cout << fixed << setw(20) << " Float EPO  " << setw(20) << _filter->param()[i].str_type() + "  " <<
			setw(20) << setprecision(5) << _filter->param()[i].value() <<
			setw(15) << setprecision(5) << _filter->dx()(i + 1) <<
			setw(20) << setprecision(5) << _filter->param()[i].value() + _filter->dx()(i + 1) << endl;
	}
# endif
	if (_fix_mode != t_gsetamb::FIX_MODE::NO/* && _epoch - _beg_time>1700*/)
	{
		if (dynamic_cast<t_gsetproc*>(_set)->realtime() == true) {
			_ewl_Upd_time = _epoch;
			_wl_Upd_time = _epoch;
		}

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
		_ambfix->setObsType(_observ);
		_ambfix->setActiveAmb(_filter->npar_number());
		if (_observ == IONO_FREE)	_ambfix->setMW(_MW[_epoch]);

		// getting the reference satellite
		bool isSetRefSat = dynamic_cast<t_gsetamb*>(_set)->isSetRefSat();
		if (isSetRefSat && !_isBase) {//ppp mode
			bool ref_valid = _getSatRef();
			if (!ref_valid) return 0;
		}
		if (!isSetRefSat) _sat_ref.clear();
		_ambfix->setSatRef(_sat_ref);

		if (_observ == RAW_ALL /*&& !_isBase*/&&_frequency >= 3)
		{
			_ambfix->processBatch(_epoch, _filter, nullptr, "EWL");
		}
		if (_observ == RAW_ALL /*&& !_isBase*/)
		{
			_ambfix->processBatch(_epoch, _filter, nullptr, "WL");
		}
		if (_ambfix->processBatch(_epoch, _filter, nullptr, "NL") < 0)
		{
			_amb_state = false;
			if (_log) _log->comment(2, "gpvtflt", _site + _epoch.str_ymdhms(" epoch ")
				+ ": fix ambiguity failed !");
		}
		else _amb_state = _ambfix->amb_fixed();

	//	dx_tmp = _filter->dx();
	//	Qx_tmp = _filter->Qx();
	}
#ifdef DEBUG
	cout << _epoch.str_ymdhms() << endl;
	for (int i = 0; i < _param_fixed.parNumber(); i++) {
		cout << fixed << setw(20) << "EPO  " << setw(20) << _filter->param()[i].str_type() + "  " <<
			setw(20) << setprecision(5) << _filter->param()[i].value() <<
			setw(15) << setprecision(5) << _filter->dx()(i + 1) <<
			setw(20) << setprecision(5) << _filter->param()[i].value() + _filter->dx()(i + 1) << endl;
	}
#endif
	// output the fixed result
	ostringstream os;
	if (_amb_state)//fixed 
	{
		_param_fixed = _ambfix->getFinalParams();
		_prtOut(_epoch, _param_fixed, _filter->Qx(), _data, os, line, true);
		_prt_port(_epoch, _param_fixed, _filter->Qx(), _data);
		//	for (unsigned int iPar = 0; iPar < _param_fixed.parNumber(); iPar++)
		//{
		//	cout << "EPO_after        " << _param_fixed[iPar].str_type() << "         " << fixed << setw(12) << setprecision(3) << _param_fixed[iPar].value() << endl;
		//}
	}
	else
	{
		for (unsigned int iPar = 0; iPar < _param_fixed.parNumber(); iPar++) {
			_param_fixed[iPar].value(_param_fixed[iPar].value() + dx_tmp(_param_fixed[iPar].index));
		}
		_prtOut(_epoch, _param_fixed, Qx_tmp, _data, os, line, true);
		_prt_port(_epoch, _param_fixed, Qx_tmp, _data);
	}

	// Print flt results
	if (_flt) {
		_flt->write(os.str().c_str(), os.str().size());
		_flt->flush();
	}

	return 1;
}

bool t_gpvtflt::_getSatRef()
{
	_sat_ref.clear();
	//one_epoch_upd ewl_upd, wl_upd, nl_upd;
	//if (_frequency >= 3) ewl_upd = _gupd->get_epo_upd(EWL_IDENTIFY);
	//wl_upd = _gupd->get_epo_upd(WL_IDENTIFY);
	//nl_upd = _gupd->get_epo_upd(_epoch);
	one_epoch_upd wl_upd, nl_upd, ewl_upd;
	if (_gupd)
	{
		wl_upd = _gupd->get_epo_upd(UPDTYPE::WL, _wl_Upd_time);
		ewl_upd = _gupd->get_epo_upd(UPDTYPE::EWL, _ewl_Upd_time);
		if (_upd_mode == t_gsetamb::UPD_MODE::UPD) nl_upd = _gupd->get_epo_upd(UPDTYPE::NL, _epoch);
	}

	t_gallpar  params_ALL = _filter->param();
	set<string> sysall = dynamic_cast<t_gsetgen*>(_set)->sys();
	for (auto  sys_iter = sysall.begin(); sys_iter != sysall.end(); sys_iter++)
	{
		GSYS sys = t_gsys::str2gsys(*sys_iter);
		double max_mw = 0.0;
		string ref;

		for (int i = 0; i < params_ALL.parNumber(); i++)
		{
			if (params_ALL[i].str_type().find("AMB") == string::npos)continue;

			if (t_gsys::sat2gsys(params_ALL[i].prn) != sys)continue;

			if (sys == GSYS::BDS && stoi(params_ALL[i].prn.substr(1, 3)) <= 5)continue; // exclude BDS C01~C05

			// reference satellite must have WL/NL UPD
			//LX 注释窄巷因为整数钟时没有窄巷upd
			// elev >= 15 angle && continuous observation epoch > 600s
			//if (wl_upd.find(params_ALL[i].prn) == wl_upd.end() || _frequency >= 3 && ewl_upd.find(params_ALL[i].prn) == ewl_upd.end()) continue;
			//if (_frequency >= 3 && ewl_upd[params_ALL[i].prn]->npoint <= 2) continue;
			//if (!((wl_upd[params_ALL[i].prn]->npoint > 2 && _epoch - params_ALL[i].beg) >= 0 && _gupd->get_epo_upd(_epoch).find(params_ALL[i].prn) != _gupd->get_epo_upd(_epoch).end()\
			//	&& _crt_ele[params_ALL[i].prn] > 15))continue;

			//if (_upd_mode == t_gsetamb::UPD_MODE::UPD)
			//{
			//	if (_gupd && !(wl_upd.find(params_ALL[i].prn) != wl_upd.end() && nl_upd.find(params_ALL[i].prn) != nl_upd.end())) continue;
			//}
			//else if (_upd_mode == t_gsetamb::UPD_MODE::IRC)
			//{
			//	if (_gupd && !(wl_upd.find(params_ALL[i].prn) != wl_upd.end()))continue;
			//}
			//if (_crt_ele[params_ALL[i].prn] < 15 || !_isClient && (_epoch - params_ALL[i].beg < 300)) continue;

				// modified 0->300 for Server !!!!!
			if (_upd_mode == t_gsetamb::UPD_MODE::UPD)
			{
				if (_gupd && _frequency >= 3 && (ewl_upd.find(params_ALL[i].prn) == ewl_upd.end() || ewl_upd[params_ALL[i].prn]->npoint <= 2)) continue;
				if (_gupd && !(wl_upd.find(params_ALL[i].prn) != wl_upd.end() && nl_upd.find(params_ALL[i].prn) != nl_upd.end())) continue;
				if (_gupd && !(wl_upd[params_ALL[i].prn]->npoint > 2)) continue;
			}
			else if (_upd_mode == t_gsetamb::UPD_MODE::IRC)
			{
				if (_gupd && !(wl_upd.find(params_ALL[i].prn) != wl_upd.end()))continue;
			}

			if (_crt_ele[params_ALL[i].prn] < 15 || !_isClient && (_epoch - params_ALL[i].beg <= 0)) continue;


			double temp = (_epoch - params_ALL[i].beg + 1) * _crt_ele[params_ALL[i].prn];
			// max  elev * continuous observation epoch 
			if (temp > max_mw)
			{
				max_mw = temp;
				ref = params_ALL[i].prn;
			}
		}
		//	cerr << "ref_sat    " << (_gupd->get_epo_upd(EWL_IDENTIFY).find(ref) != _gupd->get_epo_upd(EWL_IDENTIFY).end()) << endl;
		//	cerr << "ref_sat    " << ref << endl;
		if (!ref.empty()) {
			_sat_ref.insert(ref);
		}
	}
	if (_sat_ref.empty())return false;
	else return true;
}

int great::t_gpvtflt::_outlierDetect(const ColumnVector & v, const SymmetricMatrix& Qsav, string & sat)
{
	gtrace("t_gpvtflt::_outlierDetect");

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
			_log->comment(1, "gpvtflt", _epoch.str_ymdhms(" epoch ") + os.str());
	}
	else sat = "";

	return idx;
}
int great::t_gpvtflt::_outlierDetect(const Eigen::VectorXd & v, const SymmetricMatrix & Qsav, string & sat)
{
	gtrace("t_gpvtflt::outlierDetect");

	int nobs = v.size(); double max = 0.0; int idx = -1;
	for (int i = 0; i < nobs; i++)
	{
		if (fabs(v(i)) > max && fabs(v(i)) > _max_res_norm) {
			max = fabs(v(i));
			idx = i;
		}
	}
	if (idx >= 0) {
		_Qx = Qsav;
		sat = _obs_index[idx].first;
		string obsType = gobstype2str(_obs_index[idx - 1].second.second);
		ostringstream os;
		os << _site << " outlier (" << obsType << ") " << sat
			<< " v: " << fixed << setw(16) << right << setprecision(3) << max;
		if (_log)
			_log->comment(1, "gpvtflt", _epoch.str_ymdhms(" epoch ") + os.str());
	}
	else sat = "";

	return idx;
}


string great::t_gpvtflt::_getsat_from_obs(int idx)
{
	return string();
}

void great::t_gpvtflt::_print_ele()
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

int great::t_gpvtflt::_prt_info(const t_gtime& runEpoch)
{
	cout << runEpoch.sow() + runEpoch.dsec() << "   " << _param.parNumber() << "  ";
	// Edit and save parematres
	for (unsigned int iPar = 0; iPar < _param.parNumber(); iPar++) {
		cout << "  " + _param[iPar].str_type() << " ";
		cout << fixed << setw(18) << setprecision(6) << _param[iPar].value() << "  " << _Qx(iPar + 1, iPar + 1) << "  " ;
	}
	cout << endl;
	return 0;
}

int great::t_gpvtflt::processBatch(const t_gtime& beg_r, const t_gtime& end_r, bool prtOut)
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
		_gmutex.unlock();
		return -1;
	}

	int sign = 1;
	double subint = 0.1;

	if (!_beg_end) sign = -1;

	InitProc(beg_r, end_r, &subint);

	int smt_delay = dynamic_cast<t_gsetflt*>(_set)->smt_delay();
	t_gtime now(_beg_time);

	std::cerr << _site << ": Start GNSS Processing filtering: " << now.str_ymdhms() << " " << _end_time.str_ymdhms() << endl;
	bool time_loop = true;
	//   high_resolution_clock::time_point chrono_t1 = high_resolution_clock::now();
	while (time_loop)
	{
		if (_beg_end && (now < _end_time || now == _end_time)) time_loop = true;
		else if (_beg_end && now > _end_time) { time_loop = false; break; }

		if (!_beg_end && (now > _end_time || now == _end_time)) time_loop = true;
		else if (!_beg_end && now < _end_time) { time_loop = false; break; }

		// synchronization
		if (now != _end_time) {
			if (!time_sync(now, _sampling, _scale, _log)) { // now.sod()%(int)_sampling != 0 ){      
				now.add_dsec(sign * subint / 100); // add_dsec used for synchronization!

				continue;
			}
			if (_sampling > 1) now.reset_dsec();

		}

		// For real time procesing, added by zhShen
		if (dynamic_cast<t_gsetproc*>(_set)->realtime()) {
			while (now > _gobs->end_obs(_site))t_gtime::gmsleep(10);
			if (_isBase) {
				while (now > _gobs->end_obs(_site_base))t_gtime::gmsleep(10);
			}
			else {
				while(!dynamic_cast<t_gallprec*>(_gnav)->corr_avali(now))t_gtime::gmsleep(10);

				if (_upd_mode == t_gsetamb::UPD_MODE::UPD && _fix_mode != t_gsetamb::FIX_MODE::NO)
				{
					while (!_gupd->upd_avail(now)) t_gtime::gmsleep(10);
				}
			}

			//while (!((t_bncobs*)_gcoder[0])->available(now) || (!_isBase && !((t_bnccorr*)_gcoder[1])->available(now)))
			//	t_gtime::gmsleep(1);
		}
		/* Preprocessing for one epoch, modified by lvhb */
		if (_isClient) {
			if (!_gInterpol->interpolAug(now, _data, _data_base)) {
				if (_sampling > 1) now.add_secs(int(sign * _sampling));  // =<1Hz data
				else               now.add_dsec(sign * _sampling);       //  >1Hz data
				continue;
			}
		}
		if (dynamic_cast<t_gsetproc*>(_set)->slip_model() == SLIPMODEL::TURBO_EDIT) {
			if (_isClient) _gturboedit->setSingleEpoData(&_data);
			_gturboedit->ProcessBatch(_site, now, now, _sampling);
			if (_isBase && !_isClient) {
				_gturboedit->ProcessBatch(_site_base, now, now, _sampling);
			}
		}
		else {
			if (_isClient) _gpre->setSingleEpoData(&_data);
			_gpre->preprocess(_site, now, now, _sampling, false);
			if (_isBase && !_isClient) {
				_gpre->preprocess(_site_base, now, now, _sampling, false);
			}
		}

		int irc_epo = ProcessOneEpoch(now);
		if (irc_epo < 0) {
			if (_sampling > 1) now.add_secs(int(sign * _sampling));  // =<1Hz data
			else               now.add_dsec(sign * _sampling);       //  >1Hz data

			continue;
		}
		else _success = true;
		
		if (_log && _log->verb() >= 1)
			_log->comment(1, "gpvtflt", _site + now.str_ymdhms(" processing epoch: "));

		double percent = now.diff(_beg_time) / _end_time.diff(_beg_time) * 100.0;
		std::cerr << "\r" << now.str_ymdhms() << setw(5) << " Q = " << (_amb_state ? 1 : 2) << fixed << setprecision(1) << setw(6) << percent << "%";

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
	if (_smooth && smt_delay == 0) {
		_backflt(prtOut);
	}

	if (beg_r != end_r) {   // do not print in real-time
		double npd_perc = 0;
		npd_perc = (double(_n_NPD_flt) / double(_n_ALL_flt)) * 100;
		if (_log) _log->comment(1, "gpvtflt", _site + ": Forward filter -     Not positive definite covariance matrices: " + dbl2str(npd_perc, 0) + " %   ("
			+ int2str(_n_NPD_flt) + ", " + int2str(_n_ALL_flt) + ")");

		if (_smooth) {
			npd_perc = (double(_n_NPD_smt) / double(_n_ALL_smt)) * 100;
			if (_log) _log->comment(1, "gpvtflt", _site + ": Backward smoothing - Not positive definite covariance matrices: " + dbl2str(npd_perc, 0) + " %   ("
				+ int2str(_n_NPD_smt) + ", " + int2str(_n_ALL_smt) + ")");
		}
	}
	
	//output 3drms for all epochs,fenggl 202010
	if (_crd_est != KIN && !_pos_kin && _enufile){
		neu_sum_nepo.second /= neu_sum_nepo.first;
		ostringstream enu_info;
		enu_info << " RMS = " << fixed << setw(12) << setprecision(4) << sqrt(neu_sum_nepo.second.crd(1))
			<< setw(12) << setprecision(4) << sqrt(neu_sum_nepo.second.crd(0))
			<< setw(12) << setprecision(4) << sqrt(neu_sum_nepo.second.crd(2)) << endl;
		_enufile->write(enu_info.str().c_str(), enu_info.str().size());
		_enufile->flush();
	}

	_gmutex.unlock();

	return 1;
}

bool great::t_gpvtflt::InitProc(const t_gtime& begT,const t_gtime& endT,double* subint)
{
	if (_beg_end) {
		_beg_time = begT;
		_end_time = endT;
		if (_log) _log->comment(2, "t_gpvtflt", "Filtering in begin -> end direction!");
	}
	else {
		_beg_time = endT;
		_end_time = begT;
		if (_log) _log->comment(2, "t_gpvtflt", "Filtering in end -> begin direction!");
	}

	// get the time
	//_beg_time = (_beg_time < begT) ? begT : _beg_time;
	//_end_time = (_end_time > endT) ? endT : _end_time;

	if (_isBase) {
		if (_setCrd() < 0) return -1;
	}
	//grec->crd(t_gtriple(0), t_gtriple(1, 1, 1), _beg_time, _end_time);

	_gquality_control.setNav(_gnav);
	_dop.set_data(_gnav, _gobs, _site);
	_running = true;

	if (!_isBase){
		for (int i = 0; i < _param.parNumber(); i++)
			if (_param[i].site == "") _param.setSite(_site);
	}

	//double subint = 0.1;
	if (subint){
		if (_scale > 0) *subint = 1.0 / _scale;
		if (_sampling > 1) *subint = pow(10, floor(log10(_sampling)));
	}


	_prtOutHeader();

	return true;
}

int great::t_gpvtflt::ProcessOneEpoch(const t_gtime& now, vector<t_gsatdata>* data_rover, vector<t_gsatdata>* data_base)
{
	// clean/collect/filter epoch data
	if (_isClient || data_rover != NULL){
		if (data_rover != NULL) {
			_data.erase(_data.begin(), _data.end());
			_data = *data_rover;
		}
	}
	else {
		_data.erase(_data.begin(), _data.end());
		_data = _gobs->obs(_site, now);
	}

	if (_data.size() == 0) {
		if (_log && _log->verb() >= 1) {
			_log->comment(1, "gpvtflt", _site + now.str_ymdhms(" no observation found at epoch: "));
		}
		return -1;
	}
	
	// apply dcb
	if (_gallbias){
		for (auto& itdata : _data){
			itdata.apply_bias(*_gallbias);
		}
	}


	if (_isBase)
	{
		if (_isClient || data_base != NULL){
			if (data_base != NULL){
				_data_base.erase(_data_base.begin(), _data_base.end());
				_data_base = *data_base;
			}
		}
		else {
			_data_base.erase(_data_base.begin(), _data_base.end());
			_data_base = _gobs->obs(_site_base, now);
		}

		if (_data_base.size() == 0) {
			if (_log && _log->verb() >= 1) {
				_log->comment(1, "grtkflt", _site_base + now.str_ymdhms(" no observation found at epoch: "));
			}
			return -1;
		}
		// apply dcb
		if (_gallbias) {
			for (auto& itdata_base : _data_base){
				itdata_base.apply_bias(*_gallbias);
			}
		}
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
	
	if (_realtime)
	{
		_gobs->erase(_site, obsEpo);
		if (_isBase) _gobs->erase(_site_base, obsEpo);
	}
	
	// save apriory coordinates     
	if (_crd_est != FIX)
		_saveApr(obsEpo, _param, _Qx);

	int irc_epo = t_gpvtflt::_processEpoch(obsEpo);
	
	// Constraint finish at this epoch
	_pos_constrain = false;

	if (irc_epo < 0) {
		_success = false;
		_removeApr(obsEpo);
		if (_log && _log->verb() >= 1)
			_log->comment(1, "gpvtflt", _site + now.str_ymdhms(" epoch ") + " was not calculated");
	}
	else _success = true;

	return irc_epo;
}

t_gallpar& great::t_gpvtflt::getUCparam(set<string>* satref)
{
	return _ambfix->updateUCAmb(satref);
}


int great::t_gpvtflt::processBatchFB(const t_gtime& beg_r, const t_gtime& end_r)
{
	gtrace("t_gpvtflt::processBatchFB");

	//this->processBatch(beg_r, end_r, false);

	_beg_end = false;

	//_reset_param();

	this->processBatch(beg_r, end_r, true);

	_beg_end = true;

	return 1;
}

void great::t_gpvtflt::Add_UPD(t_gupd * gupd)
{
	_gupd = gupd;
}

int great::t_gpvtflt::_selcomsat(vector<t_gsatdata>& data_base, vector<t_gsatdata>& data_rover)
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

void great::t_gpvtflt::_udsdAmb()
{
	gtrace("t_gpppflt::_udsdAmb");

	for (map<string, int>::iterator it = _newAMB.begin(); it != _newAMB.end(); ) {
		if (it->second > 2) _newAMB.erase(it++);
		else ++it;
	}
	if (_nSat == 0) return;
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

		tuple<GOBS, GOBS, GOBS, GOBS> amb_obs_identifier = make_tuple(X, X, X, X);

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
			bool update_amb = false;
			double amb[2] = { 0 };
			t_gobs gobs1, gobs2;
			t_gobs gobs1_P, gobs2_P;
			for (int j = 0; j < 2; j++)
			{
				t_gsatdata tmpsatdata = rsatdata;
				if (j == 1) tmpsatdata = bsatdata;
				gobs1 = tmpsatdata.select_phase(b1);
				gobs1_P = tmpsatdata.select_range(b1);
				gobs2 = tmpsatdata.select_phase(b2);
				gobs2_P = tmpsatdata.select_range(b2);
				if (j == 0)
				{
					get<0>(amb_obs_identifier) = gobs1.gobs();
					get<1>(amb_obs_identifier) = gobs2.gobs();
				}
				else if (j == 1)
				{
					get<2>(amb_obs_identifier) = gobs1.gobs();
					get<3>(amb_obs_identifier) = gobs2.gobs();
				}

				if (gobs1.gobs() == X || gobs2.gobs() == X||gobs1_P.gobs() == X || gobs2_P.gobs() == X) {
					update_amb = false; break;
				}
				amb[j] = tmpsatdata.L3(gobs1, gobs2);
				if (double_eq(amb[j], 0.0)) {
					update_amb = false; break;
				}
				amb[j] -= tmpsatdata.P3(gobs1_P, gobs2_P);
				if (_param.getParam(_site, t_gpar::AMB_IF, sat) < 0 )
				{
					update_amb = true;
				}
				if (tmpsatdata.getlli(gobs1.gobs()) >= 1
					|| tmpsatdata.getlli(gobs2.gobs()) >= 1)
				{
					update_amb = true; 
					//tmpsatdata.addslip(true);
					rsatdata.addslip(true);//lvhb save slip information in rover's data
				}
			}
			if (update_amb == false)
			{
				if (_amb_obs.find(make_pair(sat, t_gpar::AMB_IF)) == _amb_obs.end())
				{
					//std::cerr << "Error: amb_obs not correct!" << sat << " " << _epoch.str_hms() << std::endl;
					if (_log) _log->comment(1, "gpvtflt", "Error: amb_obs not correct!" + sat + " " + _epoch.str_hms());
				}
				else if (_amb_obs[make_pair(sat, t_gpar::AMB_IF)] != amb_obs_identifier)
				{
					//std::cerr << "Warning: amb_obs switched silently!" << sat << " " << _epoch.str_hms() << std::endl;
					if (_log) _log->comment(1, "gpvtflt", "Warning: amb_obs switched silently!" + sat + " " + _epoch.str_hms());
					_amb_obs[make_pair(sat, t_gpar::AMB_IF)] = amb_obs_identifier;
				}
				continue;
			}
			double sdamb = amb[0] - amb[1];
			int idx = _param.getParam(_site, t_gpar::AMB_IF, sat);
			if (idx<0)
			{
				t_gpar newPar(_site, t_gpar::AMB_IF, _param.parNumber() + 1, sat);
				newPar.value(sdamb);      // first ambiguity value
				_param.addParam(newPar);
				_newAMB[sat] = 1;
				Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
				_Qx(_param.parNumber(), _param.parNumber()) = _sigAmbig * _sigAmbig;
			}
			else
			{
				_param[idx].value(sdamb);
				_Qx(idx + 1, idx + 1) = _sigAmbig * _sigAmbig;
				_newAMB[sat] = 1;
			}

			_amb_obs[make_pair(sat, t_gpar::AMB_IF)] = amb_obs_identifier;

			//cout << "ADD AMB IF " << it->sat() << " " << it->epoch().str_ymdhms() << endl;
			//if (_log) _log->comment(3, "gppp", "AMB_IF was added! For Sat PRN " + sat + " Epoch: " + _epoch.str_ymdhms());
		}
		else if (_observ == RAW_SINGLE || _observ == RAW_DOUBLE || _observ == RAW_ALL || _observ == RAW_MIX) {
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
				if( f > _frequency) continue;// modified by lvhb in 20201211
				// '>' -> '>=' . Modified by YuxuanZhou. 2020/05/31
				if (band.size() >= f) {    // automatic dual band selection -> for Anubis purpose
					b = band[f - 1];
				}
				else b = t_gsys::band_priority(gs, f);

				bool update_amb = false;
				bool skip = false;
				double amb[2] = { 0 };
				amb[0] = amb[1] = 0.0;
				t_gobs gobs;
				for (int j = 0; j < 2; j++)
				{
					t_gsatdata tmpsatdata = rsatdata;
					if (j == 1) tmpsatdata = bsatdata;
					//_getgobs(sat, TYPE_L, b, gobs);
					gobs = tmpsatdata.select_phase(b);
					if (gobs.gobs() == X) { update_amb = false; skip = true; break; }
					amb[j] = tmpsatdata.obs_L(gobs);
					if (double_eq(amb[j], 0.0)) { update_amb = false; skip = true; break; }
					if (_observ == RAW_SINGLE && f > FREQ_1) { update_amb = false; skip = true; break; }
					if (_observ == RAW_DOUBLE && f > FREQ_2) { update_amb = false; skip = true; break; }
					t_gobs gobsP = tmpsatdata.select_range(b);
					if (gobsP.gobs() == X) { update_amb = false; skip = true; break; }
					amb[j] -= tmpsatdata.obs_C(gobsP);

					if (j == 0)
					{
						get<0>(amb_obs_identifier) = gobs.gobs();
					}
					else if (j == 1)
					{
						get<2>(amb_obs_identifier) = gobs.gobs();
					}
					// add L1 amb - everytime when RAW observations
					if (_param.getParam(_site, amb_type, sat) < 0)
					{
						update_amb = true;
					}
					if (tmpsatdata.getlli(gobs.gobs()) >= 1)
					{
						update_amb = true;
						//tmpsatdata.addslip(true);
						rsatdata.addslip(true);// save slip information in rover's data
					}

				}//end base+rove
				if (skip) continue;
				if (update_amb == false)
				{
					if (_amb_obs.find(make_pair(sat, amb_type)) == _amb_obs.end())
					{
						//std::cerr << "Error: amb_obs not correct!" << sat << " " << _epoch.str_hms() << std::endl;
						if (_log) _log->comment(1, "gpvtflt", "Error: amb_obs not correct!" + sat + " " + _epoch.str_hms());
					}
					else if (_amb_obs[make_pair(sat, amb_type)] != amb_obs_identifier)
					{
						//std::cerr << "Warning: amb_obs switched silently!" << sat << " " << _epoch.str_hms() << std::endl;
						if (_log) _log->comment(1, "gpvtflt", "Warning: amb_obs switched silently!" + sat + " " + _epoch.str_hms());
						_amb_obs[make_pair(sat, amb_type)] = amb_obs_identifier;
					}
					continue;
				}
				double sdamb = amb[0] - amb[1];
				int idx = _param.getParam(_site, amb_type, sat);
				if (idx < 0)
				{
					t_gpar newPar(_site, amb_type, _param.parNumber() + 1, sat);
					newPar.value(sdamb);
					_param.addParam(newPar);
					Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
					_Qx(_param.parNumber(), _param.parNumber()) = _sigAmbig * _sigAmbig;
					//         cout << "ADD AMB L1 " << it->sat() << " " << it->epoch().str_ymdhms()  << " " << L1 - cmpObs << endl;
					//if (_log) _log->comment(3, "gppp", "RAW AMB_L1 was added! For Sat PRN " + it->sat() + " Epoch: " + _epoch.str_ymdhms());
					newAmb = 1;
				}
				else
				{
					_param[idx].value(sdamb);
					_Qx(idx + 1, idx + 1) = _sigAmbig * _sigAmbig;
					newAmb = 1;
				}

				_amb_obs[make_pair(sat, amb_type)] = amb_obs_identifier;

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

				_amb_obs.erase(make_pair(_param[i].prn, _param[i].parType));

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
void great::t_gpvtflt::_syncAmb()
{
	gtrace("t_gpppflt::_syncAmb");

	for (map<string, int>::iterator it = _newAMB.begin(); it != _newAMB.end(); ) {
		if (it->second > 2) _newAMB.erase(it++);
		else ++it;
	}

	// Add ambiguity parameter and appropriate rows/columns covar. matrix   
	set<string> mapPRN;

	vector<t_gsatdata>::iterator it;
	for (it = _data.begin(); it != _data.end(); it++) {       // loop over all satellites
		mapPRN.insert(it->sat());

		tuple<GOBS, GOBS, GOBS, GOBS> amb_obs_identifier = make_tuple(X, X, X, X);

		GSYS gs = it->gsys();

		GOBSBAND b1 = _band_index[gs][FREQ_1];
		GOBSBAND b2 = _band_index[gs][FREQ_2];
		GOBSBAND b3 = _band_index[gs][FREQ_3];
		GOBSBAND b4 = _band_index[gs][FREQ_4];
		GOBSBAND b5 = _band_index[gs][FREQ_5];


		t_gobs gobs1, gobs2, gobs3, gobs4, gobs5;
		t_gobs gobs1_P, gobs2_P, gobs3_P, gobs4_P, gobs5_P;

		gobs1 = it->select_phase(b1);
		gobs2 = it->select_phase(b2);
		gobs3 = it->select_phase(b3);
		gobs4 = it->select_phase(b4);
		gobs5 = it->select_phase(b5);
		gobs1_P = it->select_range(b1);
		gobs2_P = it->select_range(b2);
		gobs3_P = it->select_range(b3);
		gobs4_P = it->select_range(b4);
		gobs5_P = it->select_range(b5);
		//std::cout << gobs1.attr() << " " << gobs1_P.attr() << " " << "gobs:" << (gobs1.attr() == gobs1_P.attr()) << std::endl;
		//std::cout << gobs2.attr() << " " << gobs2_P.attr() << " " << "gobs:" << (gobs2.attr() == gobs2_P.attr()) << std::endl;

		double LIF, L1, L2, L3, L4, L5;
		double PIF, P1, P2, P3, P4, P5;
		LIF = L1 = L2 = L3 = L4 = L5 = 0.0;
		PIF = P1 = P2 = P3 = P4 = P5 = 0.0;

		LIF = it->L3(gobs1, gobs2); PIF = it->P3(gobs1_P, gobs2_P);
		L1 = it->obs_L(gobs1); P1 = it->obs_C(gobs1_P);
		if ((_observ == RAW_SINGLE || _observ == IF_P1 || _observ == RAW_MIX) && !double_eq(L1, 0.0)) it->tb12(true);
		L2 = it->obs_L(gobs2); P2 = it->obs_C(gobs2_P);
		if (!double_eq(L2,0.0))it->tb12(true);
		L3 = it->obs_L(gobs3); P3 = it->obs_C(gobs3);
		if (!double_eq(L3, 0.0))it->tb13(true);
		L4 = it->obs_L(gobs4); P4 = it->obs_C(gobs4);
		if (!double_eq(L4, 0.0))it->tb14(true);
		L5 = it->obs_L(gobs5); P5 = it->obs_C(gobs5);
		if (!double_eq(L5, 0.0))it->tb15(true);


		bool com = true;
		if (_gnav) com = _gnav->com();

		if (_observ == IONO_FREE || _observ == IF_P1)
		{
			int idx = _param.getParam(_site, t_gpar::AMB_IF, it->sat());
			if (idx < 0)
			{
				if (double_eq(LIF, 0.0) || double_eq(PIF, 0.0)) continue;

				t_gpar newPar(it->site(), t_gpar::AMB_IF, _param.parNumber() + 1, it->sat());

				newPar.value(LIF - PIF);      // first ambiguity value
				newPar.setTime(_epoch, LAST_TIME); // beg -> end
				_param.addParam(newPar);
				_newAMB[it->sat()] = 1;

				Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
				_Qx(_param.parNumber(), _param.parNumber()) = _sigAmbig * _sigAmbig;
				//cout << "ADD AMB IF " << it->sat() << " " << it->epoch().str_ymdhms() << endl;
				if (_log) _log->comment(3, "gppp", "AMB_IF was added! For Sat PRN " + it->sat() + " Epoch: " + _epoch.str_ymdhms());
			}
			else if (it->getlli(gobs1.gobs()) >= 1 || it->getlli(gobs2.gobs()) >= 1 
				/*|| _isClient && _param[idx].amb_ini == true && _param[idx].beg < _epoch*/)
			{
				if (double_eq(LIF, 0.0) || double_eq(PIF, 0.0)) continue;
				it->addslip(true);
				_param[idx].value(LIF - PIF);
				_param[idx].setTime(_epoch, LAST_TIME); // beg -> end  (because of cycle slip)
				_Qx(idx + 1, idx + 1) = _sigAmbig * _sigAmbig;
				_newAMB[it->sat()] = 1;
			}

			// check amb_obs
			amb_obs_identifier = make_tuple(gobs1.gobs(), gobs2.gobs(), X, X);
			if (idx < 0 || it->getlli(gobs1.gobs()) >= 1 || it->getlli(gobs2.gobs()) >= 1)
				_amb_obs[make_pair(it->sat(), t_gpar::AMB_IF)] = amb_obs_identifier;
			else if (_amb_obs.find(make_pair(it->sat(), t_gpar::AMB_IF)) == _amb_obs.end())
			{
				//std::cerr << "Error: amb_obs not correct!" << it->sat() << " " << _epoch.str_hms() << std::endl;
				if (_log) _log->comment(1, "gpvtflt", "Error: amb_obs not correct!" + it->sat() + " " + _epoch.str_hms());
			}
			else if (_amb_obs[make_pair(it->sat(), t_gpar::AMB_IF)] != amb_obs_identifier)
			{
				//std::cerr << "Warning: amb_obs switched silently!" << it->sat() << " " << _epoch.str_hms() << std::endl;
				if (_log) _log->comment(1, "gpvtflt", "Warning: amb_obs switched silently!" + it->sat() + " " + _epoch.str_hms());
				_amb_obs[make_pair(it->sat(), t_gpar::AMB_IF)] = amb_obs_identifier;
			}
		}
		else if (_observ == RAW_SINGLE || _observ == RAW_DOUBLE || _observ == RAW_ALL || _observ == RAW_MIX)
		{
			int newAmb = 0;
			int nf = 5;
			if (_observ == RAW_SINGLE) nf = 1;
			else if (_observ == RAW_DOUBLE) nf = 2;
			else if (_observ == RAW_ALL || _observ == RAW_MIX) nf = 5;
			//nf = band.size();
			for (const auto& ib : _band_index[gs])
			{
			//	t_gpar::t_type amb_type = t_gpar::AMB_IF;
				if (ib.first > _frequency) continue;
				t_gobs gobsi; double Li, Pi; t_gpar::t_type amb_type;

				if      (ib.first == FREQ_1) { amb_type = t_gpar::AMB_L1; Li = L1; Pi = P1; gobsi = gobs1; }
				else if (ib.first == FREQ_2) { amb_type = t_gpar::AMB_L2; Li = L2; Pi = P2; gobsi = gobs2; }
				else if (ib.first == FREQ_3) { amb_type = t_gpar::AMB_L3; Li = L3; Pi = P3; gobsi = gobs3; }
				else if (ib.first == FREQ_4) { amb_type = t_gpar::AMB_L4; Li = L4; Pi = P4; gobsi = gobs4; }
				else if (ib.first == FREQ_5) { amb_type = t_gpar::AMB_L5; Li = L5; Pi = P5; gobsi = gobs5; }
				else {
					throw std::logic_error("Unknown band!!");
				};

				if (double_eq(Li, 0.0) || double_eq(Pi, 0.0)) continue;
				int idx = -1;
				// add L1 amb - everytime when RAW observations
				idx = _param.getParam(_site, amb_type, it->sat());
				if (idx < 0) {
					t_gpar newPar(it->site(), amb_type, _param.parNumber() + 1, it->sat());
					newPar.value(Li - Pi);
					newPar.setTime(_epoch, LAST_TIME); // beg -> end
					_param.addParam(newPar);
					
					Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
					_Qx(_param.parNumber(), _param.parNumber()) = _sigAmbig * _sigAmbig;
					//         cout << "ADD AMB L1 " << it->sat() << " " << it->epoch().str_ymdhms()  << " " << L1 - cmpObs << endl;
					if (_log) _log->comment(3, "gppp", "RAW AMB_L1 was added! For Sat PRN " + it->sat() + " Epoch: " + _epoch.str_ymdhms());
					newAmb = 1;
					//cout << "AMB will be Added! For Sat PRN " << it->sat()
					//	<< " Epoch: " << _epoch.str_ymdhms() << endl;
				}
				else if (it->getlli(gobsi.gobs()) >= 1
					/*|| _isClient && _param[idx].amb_ini == true && _param[idx].beg < _epoch*/) // check cycle slip
				{
					it->addslip(true);
					_param[idx].value(Li - Pi);
					_param[idx].setTime(_epoch, LAST_TIME); // beg -> end  (because of cycle slip)
					_Qx(idx + 1, idx + 1) = _sigAmbig * _sigAmbig;
					_newAMB[it->sat()] = 1;
					//cout << "AMB will be reset! For Sat PRN " << it->sat()
					//	<< " Epoch: " << _epoch.str_ymdhms() << endl;
				}

				// check amb_obs
				amb_obs_identifier = make_tuple(gobsi.gobs(), X, X, X);
				if (idx < 0 || it->getlli(gobsi.gobs()) >= 1)
					_amb_obs[make_pair(it->sat(), amb_type)] = amb_obs_identifier;
				else if (_amb_obs.find(make_pair(it->sat(), amb_type)) == _amb_obs.end())
				{
					//std::cerr << "Error: amb_obs not correct!" << it->sat() << " " << _epoch.str_hms() << std::endl;
					if (_log) _log->comment(1, "gpvtflt", "Error: amb_obs not correct!" + it->sat() + " " + _epoch.str_hms());
				}
				else if (_amb_obs[make_pair(it->sat(), amb_type)] != amb_obs_identifier)
				{
					//std::cerr << "Warning: amb_obs switched silently!" << it->sat() << " " << _epoch.str_hms() << std::endl;
					if (_log) _log->comment(1, "gpvtflt", "Warning: amb_obs switched silently!" + it->sat() + " " + _epoch.str_hms());
					_amb_obs[make_pair(it->sat(), amb_type)] = amb_obs_identifier;
				}
			}

			if (newAmb) _newAMB[it->sat()] = 1;
		}


	}  // end loop over all observations

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
				//cout << "AMB will be removed! For Sat PRN " << _param[i].prn
				//	<< " Epoch: " << _epoch.str_ymdhms() << endl;
				_amb_obs.erase(make_pair(_param[i].prn, _param[i].parType));

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

void great::t_gpvtflt::_post_turbo_syncAmb()
{
	gtrace("t_gpppflt::_post_turbo_syncAmb");

	for (map<string, int>::iterator it = _newAMB.begin(); it != _newAMB.end(); ) {
		if (it->second > 2) _newAMB.erase(it++);
		else ++it;
	}

	t_gupdateparinfo update_info;
	_update_AMB->update_amb_pars(_epoch, _param, _data, update_info);

	vector<int> remove_id;
	vector<t_gpar> newparlist;
	update_info.get(remove_id);
	update_info.get(newparlist);

	set<string> addSats;
	// First. // Add ambiguity parameter and appropriate rows/columns covar. matrix   
	for (auto par : newparlist) {

		par.index = _param.parNumber() + 1;
		_param.addParam(par);
		_newAMB[par.prn] = 1;
		addSats.insert(par.prn);

		Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
		_Qx(_param.parNumber(), _param.parNumber()) = _sigAmbig * _sigAmbig;

		if (_log) _log->comment(3, "gppp", "AMB_IF was added! For Sat PRN " + par.prn + " Epoch: " + _epoch.str_ymdhms());
	}


	// Second. Remove old pars  // Remove ambiguity parameter and appropriate rows/columns covar. matrix
	sort(remove_id.begin(), remove_id.end());
	for (int i = 0; i < remove_id.size(); i++)
	{
		string prn = _param[remove_id[i] - i - 1].prn;
		int index = _param[remove_id[i] - i - 1].index;
		if (addSats.find(prn) == addSats.end())
		{
			_newAMB.erase(prn);
			if (_log)
				_log->comment(3, "gppp", "AMB will be removed! For Sat PRN " + prn + " Epoch: " + _epoch.str_ymdhms());
		}

		Matrix_remRC(_Qx, index, index);
		_param.delParam(remove_id[i] - i - 1);
		_param.reIndex();
	}
	return;
}

int great::t_gpvtflt::_get_spplsq_result(string site, t_gtime current, t_gtriple & crd, double & clk)
{
	t_gspplsq* _spplsq = new t_gspplsq(site, _set, _gallobj, false);
	_spplsq->setDAT(_gobs, _gnav);
	_spplsq->glog(_log);
	if (_spplsq->processBatch(current, current) <= 0) {
		if (_log) {
			_log->comment(1, site + current.str_ymdhms() + "calculate initial crd failed using spplsq");
		}
		delete _spplsq; _spplsq = nullptr;
		return -1;
	};
	crd = _spplsq->getCrd(current);
	clk = _spplsq->getRecClk(current);
	delete _spplsq; _spplsq = nullptr;
	return 1;
}

void great::t_gpvtflt::Add_rho_azel(/*const t_gtime & runEpoch,*/ const string & site_name,  t_gtriple & xyz_s, const t_gtriple & xyz_r, t_gsatdata & obs_sat)
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
	if (_isBase && site_name != _site)
	{
		string tmp = _site;
		_site = site_name;
		_apply_tides(_epoch, xRec);
		_site = tmp;
	}
	else _apply_tides(_epoch, xRec);

	double tmp = (xyz_s.crd_cvect() - xRec.crd_cvect()).norm_Frobenius();

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

void great::t_gpvtflt::_prtOut(t_gtime& epoch, t_gallpar& X, const SymmetricMatrix& Q, vector<t_gsatdata>& data, ostringstream& os, xml_node& node, bool saveProd)
{
	gtrace("t_gpvtflt::_prtOut");

	// get CRD params
	t_gtriple xyz, ell;
	X.getCrdParam(_site, xyz);
	if (_crd_est == FIX) {
		shared_ptr<t_gobj> grec = _gallobj->obj(_site);
		xyz = grec->crd_arp(epoch);
	}
	xyz2ell(xyz, ell, false);

	// CRD using eccentricities
	t_gtriple xyz_ecc = xyz - _grec->eccxyz(epoch); // MARKER + ECC = ARP

													// get CRD rms  (XYZ)
	double Xrms = 0.0, Yrms = 0.0, Zrms = 0.0,
		Vxrms = 0.0, Vyrms = 0.0, Vzrms = 0.0;
	double cov_xy = 0.0, cov_xz = 0.0, cov_yz = 0.0;
	int icrdx = _param.getParam(_site, t_gpar::CRD_X, "");
	int icrdy = _param.getParam(_site, t_gpar::CRD_Y, "");
	int icrdz = _param.getParam(_site, t_gpar::CRD_Z, "");
	if (icrdx>=0 && icrdy>=0 &icrdz>=0)
	{
		if (Q(icrdx + 1, icrdx + 1) < 0) Xrms = -1; else Xrms = sqrt(Q(icrdx + 1, icrdx + 1));
		if (Q(icrdy + 1, icrdy + 1) < 0) Yrms = -1; else Yrms = sqrt(Q(icrdy + 1, icrdy + 1));
		if (Q(icrdz + 1, icrdz + 1) < 0) Zrms = -1; else Zrms = sqrt(Q(icrdz + 1, icrdz + 1));
		cov_xy = Q(icrdy + 1, icrdx + 1);
		cov_xz = Q(icrdz + 1, icrdx + 1);
		cov_yz = Q(icrdz + 1, icrdy + 1);
	}

	t_gtriple crd_rms(Xrms, Yrms, Zrms);
	t_gtriple vRec(0, 0, 0);

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
		ZHD = X[iztd].zhd;//_gModel->tropoModel()->getZHD(ell, epoch);
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
	bool ins = dynamic_cast<t_gsetinp*>(_set)->input_size("imu") > 0 ? true : false;
	// write kml
	if (_kml && !ins) {
		ostringstream out;
		out << fixed << setprecision(11) << " " << setw(0) << ell[1] * R2D << ',' << ell[0] * R2D;
		string val = out.str();
		
		xml_node root = _doc;
		node = this->_default_node(root, _root.c_str());
		xml_node document = node.child("Document");
		xml_node last_child = document.last_child();
		xml_node placemark = document.insert_child_after("Placemark", last_child);
		if(_amb_state)this->_default_node(placemark, "styleUrl", "#P1"); // value: #P1 (according to precision)
		else this->_default_node(placemark, "styleUrl", "#P3");
		this->_default_node(placemark, "time", int2str(_epoch.sow()).c_str());
		xml_node point = this->_default_node(placemark, "Point");
		this->_default_node(point, "coordinates", val.c_str()); // for point
		
		xml_node Placemark = document.child("Placemark");
		xml_node LineString = Placemark.child("LineString");
		this->_default_node(LineString, "coordinates", val.c_str(), false);  // for line
	}

	string str_dsec = dbl2str(epoch.dsec());

	double bl = 0;
	if (_isBase)
	{
		auto crd_base = _gallobj->obj(_site_base)->crd_arp(epoch);
		t_gtriple tmpell, tmpdxyz,tmpneu;
		xyz2ell(crd_base, tmpell, false);
		tmpdxyz = xyz - crd_base;
		xyz2neu(tmpell, tmpdxyz, tmpneu);
		bl = tmpneu.norm();
	}

	os << fixed << setprecision(4) << " "
		// << epoch.str_ymdhms() << str_dsec.substr(2) << setprecision(4)
		<< " " << epoch.sow() + epoch.dsec();
		if (_crd_est != FIX) {
		os << fixed << setprecision(4)
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
			;
		}
		os<< fixed << setprecision(0)
		<< " " << setw(5) << nsat            // nsat
		<< fixed << setprecision(1)
		//<< " " << setw(3) << gdop            // gdop
		<< fixed << setprecision(2)
		<< " " << setw(5) << pdop            // pdop
		<< fixed << setprecision(2)
		<< " " << setw(8) << _sig_unit       // m0
		<< " " << setw(8) << amb;
		if (_fix_mode != t_gsetamb::FIX_MODE::NO)	os << setprecision(2) << " " << fixed << setw(10) << _ambfix->get_ratio();
		if (_isBase) os << " " <<fixed<<setprecision(3)<< setw(15) << bl;//lvhb
		if (_observ == RAW_MIX)
		{
			map<string, int> satsnum = _param.freq_sats_num(2);
			os << setw(8) << satsnum["Single"] << setw(8) << satsnum["Double"];
		}
		os<< endl;
	// write enu
	t_gtriple xyz_ecc1 = xyz_ecc;
	if (_enufile)
	{
		t_gtriple dxyz, ell, neu;

		ostringstream enu_info;
		//double crt_t = epoch.sow() + epoch.dsec();
		double crt_t = epoch.sod() + epoch.dsec();


		if (_crd_est == CONSTRPAR::KIN || _pos_kin)
		{
			// just for now
			enu_info <<  _site << " " << fixed << setw(8) << epoch.mjd() << 
				fixed << setw(10) << setprecision(1) << crt_t;
			//enu_info << fixed << setw(10) << setprecision(1) << crt_t;
			enu_info << setw(16) << setprecision(4) << xyz_ecc1[0]
				<< setw(16) << setprecision(4) << xyz_ecc1[1]
				<< setw(16) << setprecision(4) << xyz_ecc1[2] << "  " << nsat << endl;
		}
		else
		{
			enu_info << fixed << setw(10) << setprecision(1) << crt_t;
			xyz2ell(xyz_standard, ell, false);
			dxyz = xyz_ecc1 - xyz_standard;
			xyz2neu(ell, dxyz, neu);

			// RMS
			neu_sum_nepo.first++;   // nepoch
			neu_sum_nepo.second += t_gtriple(pow(neu.crd(0),2), pow(neu.crd(1),2),pow(neu.crd(2),2)); // sum

			enu_info << fixed << setw(16) << setprecision(4) << neu[1]
				<< setw(16) << setprecision(4) << neu[0]
				<< setw(16) << setprecision(4) << neu[2] 
				<< " " << setw(8) << amb<< endl;
		}
		_enufile->write(enu_info.str().c_str(), enu_info.str().size());
		_enufile->flush();
	}


	return;
}

void great::t_gpvtflt::_prt_port(t_gtime & epoch, t_gallpar & X, const SymmetricMatrix & Q, vector<t_gsatdata>& data)
{
	gtrace("t_gpvtflt::_prt_port");
	ostringstream os;
	// get CRD params
	t_gtriple xyz, ell;
	X.getCrdParam(_site, xyz);
	if (_crd_est == FIX) {
		shared_ptr<t_gobj> grec = _gallobj->obj(_site);
		xyz = grec->crd_arp(epoch);
	}
	xyz2ell(xyz, ell, false);

	// CRD using eccentricities
	t_gtriple xyz_ecc = xyz - _grec->eccxyz(epoch); // MARKER + ECC = ARP

	double Xrms = 0.0, Yrms = 0.0, Zrms = 0.0,
		Vxrms = 0.0, Vyrms = 0.0, Vzrms = 0.0;
	int icrdx = _param.getParam(_site, t_gpar::CRD_X, "");
	int icrdy = _param.getParam(_site, t_gpar::CRD_Y, "");
	int icrdz = _param.getParam(_site, t_gpar::CRD_Z, "");
	if (icrdx >= 0 && icrdy >= 0 & icrdz >= 0)
	{
		if (Q(icrdx + 1, icrdx + 1) < 0) Xrms = -1; else Xrms = sqrt(Q(icrdx + 1, icrdx + 1));
		if (Q(icrdy + 1, icrdy + 1) < 0) Yrms = -1; else Yrms = sqrt(Q(icrdy + 1, icrdy + 1));
		if (Q(icrdz + 1, icrdz + 1) < 0) Zrms = -1; else Zrms = sqrt(Q(icrdz + 1, icrdz + 1));
	}

	t_gtriple crd_rms(Xrms, Yrms, Zrms);
	t_gtriple vRec(0, 0, 0);
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

	t_gtriple blh;
	xyz2ell(xyz_ecc, blh, true);

	string str_dsec = dbl2str(epoch.dsec());

	double bl = 0;
	if (_isBase)
	{
		auto crd_base = _gallobj->obj(_site_base)->crd_arp(epoch);
		t_gtriple tmpell, tmpdxyz, tmpneu;
		xyz2ell(crd_base, tmpell, false);
		tmpdxyz = xyz - crd_base;
		xyz2neu(tmpell, tmpdxyz, tmpneu);
		bl = tmpneu.norm();
	}

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
		<< fixed << setprecision(0)
		<< " " << setw(5) << nsat            // nsat
		<< fixed << setprecision(1)
		<< fixed << setprecision(2)
		<< " " << setw(5) << pdop            // pdop
		<< fixed << setprecision(2)
		<< " " << setw(8) << amb;
		if (_fix_mode != t_gsetamb::FIX_MODE::NO)	os << setprecision(2) << " " << fixed << setw(10) << _ambfix->get_ratio();
		else                                        os << setprecision(2) << " " << fixed << setw(10) << 0.0;
		os << endl;
	bool ins = dynamic_cast<t_gsetinp*>(_set)->input_size("imu") > 0 ? true : false;
	if (!ins && _maptcp.find(FLT_OUT) != _maptcp.end())_maptcp[FLT_OUT]->run_send(os.str());


}

void great::t_gpvtflt::_get_result(t_gtime& epoch, t_gposdata::data_pos & pos)
{
	gtrace("t_gpvtflt::_get_result");

	double crt = epoch.sow() + epoch.dsec();
	t_gallpar X = _param_fixed;
	// get CRD params
	t_gtriple xyz, xyz_ecc;
	// get CRD rms  (XYZ)
	double Xrms = 0.0, Yrms = 0.0, Zrms = 0.0,
		Vxrms = 0.0, Vyrms = 0.0, Vzrms = 0.0;


	if (X.getCrdParam(_site, xyz) > 0)
	{
		// CRD using eccentricities
		xyz_ecc = xyz - _grec->eccxyz(epoch); // MARKER + ECC = ARP

		int icrdx = _param.getParam(_site, t_gpar::CRD_X, "");
		int icrdy = _param.getParam(_site, t_gpar::CRD_Y, "");
		int icrdz = _param.getParam(_site, t_gpar::CRD_Z, "");
		if (_Qx(icrdx + 1, icrdx + 1) < 0) Xrms = -1; else Xrms = sqrt(_Qx(icrdx + 1, icrdx + 1));
		if (_Qx(icrdy + 1, icrdy + 1) < 0) Yrms = -1; else Yrms = sqrt(_Qx(icrdy + 1, icrdy + 1));
		if (_Qx(icrdz + 1, icrdz + 1) < 0) Zrms = -1; else Zrms = sqrt(_Qx(icrdz + 1, icrdz + 1));

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
void great::t_gpvtflt::_prtOutHeader()
{
	gtrace("t_gpvtflt::_prtOutHeader");
	ostringstream os;

	//lvhb added in 20200818
	if (_isBase)
	{
		auto beg = dynamic_cast<t_gsetgen*>(_set)->beg();
		auto crd_base = _gallobj->obj(_site_base)->crd(beg);
		os << "# base -" << _site_base << " Pos(XYZ): " << fixed << setprecision(4) << crd_base[0] << setw(15) << crd_base[1] << setw(15) << crd_base[2] << endl;
		if (_crd_est == FIX) {
			auto crd = _gallobj->obj(_site)->crd(beg);
			os << "# rover-" << _site << " Pos(XYZ): " << fixed << setprecision(4) << crd[0] << setw(15) << crd[1] << setw(15) << crd[2] << endl;
		}
	}

	// first line
	os << "#" <<
		setw(15) << "Seconds of Week";
	if (_crd_est != FIX) {
		os <<
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
			" " << setw(9) << "Vz-RMS";
		}
	os<<" " << setw(5) << "NSat" <<
		" " << setw(5) << "PDOP" <<
		" " << setw(8) << "sigma0" <<
		" " << setw(10) << "AmbStatus";
		if (_fix_mode != t_gsetamb::FIX_MODE::NO) os<<" " << setw(5) << "Ratio";	//lvhb added
		if (_isBase)	os << " " << setw(15) << "BL";								//lvhb added
		os<<endl;

	// second line
		os << "#" <<
			setw(15) << "(s)";
		if (_crd_est != FIX) {
			os <<
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
				" " << setw(9) << "(m/s)"          // [m/s]
				;
		}
		os<<" " << setw(5) << "(#)" <<
			" " << setw(5) << "(#)" <<
			" " << setw(8) << "(m)" <<
			" " << setw(10) << " ";
		if (_fix_mode != t_gsetamb::FIX_MODE::NO) os << setw(5) << " ";	//lvhb added
		if (_isBase)	os << " " << setw(15) << "(m)";					//lvhb added
		os<<endl;

	// Print flt results
	if (_flt) {
		_flt->write(os.str().c_str(), os.str().size());
		_flt->flush();
	}
}

void great::t_gpvtflt::_combineEWL(t_gsatdata & satdata)
{
}

void great::t_gpvtflt::_generateObsIndex(t_glsqEquationMatrix& equ)
{
	_obs_index.clear();
	auto band_index = dynamic_cast<t_gcombmodel*>(&(*_base_model))->get_band_index();
	auto freq_index = dynamic_cast<t_gcombmodel*>(&(*_base_model))->get_freq_index();
	for (int i = 0; i < equ.num_equ(); i++)
	{
		char sys_char = equ.get_satname(i)[0];
		GSYS sys;
		bool sys_valid = true;
		switch (sys_char)
		{
		case 'G':sys = GSYS::GPS; break;
		case 'R':sys = GSYS::GLO; break;
		case 'E':sys = GSYS::GAL; break;
		case 'C':sys = GSYS::BDS; break;
		case 'J':sys = GSYS::QZS; break;
		default:sys_valid=false; break;
		}
		if (!sys_valid) throw exception();
		GOBSBAND b = equ.get_obscombtype(i).getBand();
		FREQ_SEQ f = freq_index[sys][b];
		GOBSTYPE obstype = equ.get_obscombtype(i).is_code() ? TYPE_C : GOBSTYPE::TYPE_L;
		_obs_index.push_back(make_pair(equ.get_satname(i), make_pair(f, obstype)));
	}
}

void great::t_gpvtflt::_equationToMatrix(t_glsqEquationMatrix& equ, Matrix& A, SymmetricMatrix& P,
	ColumnVector& l, unsigned int nPar)
{
	unsigned nObs = equ.num_equ();
	A.resize(nObs, nPar); A = 0.0;
	P.resize(nObs); P = 0.0;
	l.resize(nObs); l = 0.0;
	for (int i = 0; i < equ.B.size(); i++)
	{
		for (int j = 0; j < equ.B[i].size(); j++)
		{
			if (equ.B[i][j].first > nPar) continue;
			A(i+1, equ.B[i][j].first) = equ.B[i][j].second;
		}
		l(i + 1) = equ.l[i];
		P(i + 1, i + 1) = equ.P[i];
	}
}

void great::t_gpvtflt::_extern_ion_constraint(const t_gtime& runEpoch)
{
	_epoch_cout++;
	t_giono_tecgrid iono_tec;
	vector<t_gsatdata>::iterator it = _data.begin();
	while (it != _data.end())
	{
		int ipar_sion = _param.getParam(_site, t_gpar::SION, it->sat());
		if (ipar_sion > 0)
		{
			double ion_L1, ion_rms;
			t_gtriple site_pos = t_gtriple(_vBanc);
			if (iono_tec.getIonoDelay(_gionex_GIM, *it, runEpoch, site_pos, ion_L1, ion_rms))
			{
				_param[ipar_sion].value(ion_L1);
				_Qx(ipar_sion + 1, ipar_sion + 1) = 1.0 + 0.1 * _epoch_cout;
			}
		}
		it++;
	}

}

bool great::t_gpvtflt::_calculate_timediff_baseline(const string& site, const t_gtime & epoch0, const t_gtime & epoch1, t_gtriple & baseline)
{
	std::vector<gnut::t_gsatdata> data_t0 = _gobs->obs(site, epoch0);
	std::vector<gnut::t_gsatdata> data_t1 = _gobs->obs(site, epoch1);

	if (data_t0.size() < 5 || data_t1.size() < 5) return false;

	t_gtriple crd_t0, crd_t1;
	double clk_t0, clk_t1;
	t_gallpar param_t0;
	t_gallpar param_t1;

	// solve bancroft to get approximate coordinates
	ColumnVector spp_t0, spp_t1;
	_prepareData(site, data_t0); spp_t0 = _vBanc;
	_prepareData(site, data_t1); spp_t1 = _vBanc;

	int max_obs = 1000;
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(max_obs, 4);
	Eigen::VectorXd l = Eigen::VectorXd::Zero(max_obs);
	Eigen::VectorXd x;
	int iobs = 0;

	Eigen::Vector3d rr_t0, rr_t1;
	rr_t0 = Eigen::Vector3d(spp_t0[0], spp_t0[1], spp_t0[2]);
	rr_t1 = Eigen::Vector3d(spp_t1[0], spp_t1[1], spp_t1[2]);

	// use only L1 carrier phase data
	// construct a least-square problem
	for (auto& satdata_t0 : data_t0)
	{
		string sat = satdata_t0.sat();
		for (auto& satdata_t1 : data_t1)
		{
			if (sat == satdata_t1.sat())
			{
				vector < GOBSBAND > band = dynamic_cast<t_gsetgnss*>(_set)->band(satdata_t0.gsys());
				GOBSBAND b1 = band[0];
				t_gobs gobsL1 = t_gobs(satdata_t0.select_phase(b1));
				double obs_L1_t0 = satdata_t0.obs_L(gobsL1);
				double obs_L1_t1 = satdata_t1.obs_L(gobsL1);
				if (double_eq(obs_L1_t0, 0.0) || double_eq(obs_L1_t1, 0.0)
					|| satdata_t0.lli()[gobsL1.gobs()] > 0 
					|| satdata_t1.lli()[gobsL1.gobs()] > 0) continue;
				
				Eigen::Vector3d rs_t0, rs_t1;
				rs_t0 = Eigen::Vector3d(satdata_t0.satcrd().crd(0), satdata_t0.satcrd().crd(1), satdata_t0.satcrd().crd(2));
				rs_t1 = Eigen::Vector3d(satdata_t1.satcrd().crd(0), satdata_t1.satcrd().crd(1), satdata_t1.satcrd().crd(2));

				double model = (rs_t1.dot(rs_t1) - rs_t0.dot(rs_t0) - 2 * rr_t0.dot(rs_t1 - rs_t0))
					/ ((rs_t0 - rr_t0).norm() + (rs_t1 - rr_t1).norm())
					//+
					//(-2 * rs_t1 + rr_t1 + rr_t0).dot(rr_t1 - rr_t0)
					/// ((rs_t0 - rr_t0).norm() + (rs_t1 - rr_t1).norm())
					+ satdata_t0.clk()
					- satdata_t1.clk();

				double omc = (obs_L1_t1 - obs_L1_t0) - model;
				
				B.block(iobs, 0, 1, 3) = ((-2 * rs_t1 + rr_t1 + rr_t0) / ((rs_t0 - rr_t0).norm() + (rs_t1 - rr_t1).norm())).transpose();
				B(iobs, 3) = 1;
				l(iobs) = omc;
				iobs++;
			}
		}
	}

	B.conservativeResize(iobs, 4);
	l.conservativeResize(iobs);

	x = (B.transpose()*B).inverse()*B.transpose()*l;

	baseline[0] = x(0);
	baseline[1] = x(1);
	baseline[2] = x(2);
	return true;
}

bool great::t_gpvtflt::_post_turbo_select_obs(const t_gtime& epoch, vector<t_gsatdata>& data)
{
	auto iter = data.begin();
	while (iter != data.end())
	{
		const string& obs_rec = iter->site();
		const string& obs_sat = iter->sat();


		bool obs_valid = _slip12->use_of_obs(obs_rec, obs_sat, epoch);
		if (obs_valid)
		{
			iter->tb12(obs_valid);
			if (_frequency >= 3)
			{
				bool obs_13 = _slip13->use_of_obs(obs_rec, obs_sat, epoch);
				iter->tb13(obs_13);
			}
			if (_frequency >= 5)
			{
				bool obs_14 = _slip14->use_of_obs(obs_rec, obs_sat, epoch);
				bool obs_15 = _slip15->use_of_obs(obs_rec, obs_sat, epoch);
				iter->tb14(obs_14);
				iter->tb15(obs_15);
			}
			iter++;
		}
		else
		{
			iter = data.erase(iter);
		}
	}

	if (data.empty())
	{
		_log->logDebug("t_gpvtflt", "_post_turbo_select_obs", "This epoch have no useful data : " + epoch.str_ymdhms());
		return false;
	}

	return true;

}