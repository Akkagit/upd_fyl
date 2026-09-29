#include "gnpp/gInterpol.h"
#include "gall/gallobs.h"
#include "gset/gsetrec.h"
#include "gproc/gppplsq.h"
#include "gnpp/gDelaunay.h"

great::t_gInterpol::t_gInterpol(t_gsetbase* set, t_gallproc* data, t_glog* log):
	_set(set),
	_alldata(data),
	_log(log),
	_isBase(false),
	_isCompAug(true)
{
	if (dynamic_cast<t_gsetgen*>(_set)->list_base().size()) _isBase = true;
	_allaug = dynamic_cast<t_gaug*>((*_alldata)[t_gdata::AUG]);
	_allobj = dynamic_cast<t_gallobj*>((*_alldata)[t_gdata::ALLOBJ]);
	_site = dynamic_cast<t_gsetgen*>(_set)->list_rover()[0];
	t_gtime begt = dynamic_cast<t_gsetgen*>(_set)->beg();
	t_gtime endt = dynamic_cast<t_gsetgen*> (_set)->end();

	//Step1:set interpolation site list from gset
	auto sites = dynamic_cast<t_gsetgen*>(_set)->recs();
	_list_ipsite.clear();
	int i = 0;
	for (auto it = sites.begin(); it != sites.end(); it++) {
		if (*it == _site) continue;
		_list_ipsite.push_back(*it);
		//auto crd = dynamic_cast<t_gsetrec*>(_set)->get_crd_xyz(*it);
	}
	_nipSite = _list_ipsite.size();

	//Setp2:initialize posproc class
	_site_vrs = "";
	if (_isBase) {
		_site_vrs = dynamic_cast<t_gsetgen*>(_set)->list_base()[0];
		//add vrs-site obj
		if (_allobj->obj(_site_vrs) == nullptr) {
			shared_ptr<t_grec> rec = make_shared<t_grec>();
			rec->name(_site_vrs);
			rec->id(_site_vrs);
			_allobj->add(rec);
		}
	}
}

great::t_gInterpol::~t_gInterpol()
{
}

void great::t_gInterpol::_setCrd(ColumnVector& dx, ColumnVector& dy, ColumnVector& dz)
{
	_crd_dx = dx;
	_crd_dy = dy;
	_crd_dz = dz;
}

double great::t_gInterpol::_MLCM(ColumnVector& AugData)
{
	// double dAug = 0.0;
	int n = AugData.Nrows();

	_coef_a.CleanUp();
	_coef_a.ReSize(n); _coef_a = 0;

	ColumnVector unit(n); unit = 1.0;
	Matrix B; B.ReSize(3, n);
	ColumnVector l(3); l = 0.0; l(1) = 1.0;
	B.Row(1) = unit.t();
	B.Row(2) = _crd_dx.t();
	B.Row(3) = _crd_dy.t();

	//X= (BtB).i() * (BtL)
	//_coef_a << (B.t() * B).i() * (B.t() * l);
	_coef_a << B.t() * (B * B.t()).i() * (l);
	//cout << _coef_a << endl;
	Matrix dAug = _coef_a.t() * AugData;

	return dAug(1, 1);
}

double great::t_gInterpol::_URTKInterpol(ColumnVector& dAugData)
{
	auto XY_12 = _crd_dx(2) * _crd_dy(1) - _crd_dx(1) * _crd_dy(2);
	auto XY_23 = _crd_dx(2) * _crd_dy(3) - _crd_dx(3) * _crd_dy(2);
	auto XY_13 = _crd_dx(3) * _crd_dy(1) - _crd_dx(1) * _crd_dy(3);

	auto a = XY_23 / XY_12;
	auto b = XY_13 / XY_12;

	auto daug = a * dAugData(1) + b * dAugData(2);

	return daug;
}
double great::t_gInterpol::_NRTKInterpol(ColumnVector& dAugData)
{
	ColumnVector v; v.ReSize(2); v = 0.0;
	v(1) = _crd_dx(3);
	v(2) = _crd_dy(3);
	Matrix A;	A.ReSize(2, 2); A = 0.0;
	A(1, 1) = _crd_dx(1); A(1, 2) = _crd_dy(1);
	A(2, 1) = _crd_dx(2); A(2, 2) = _crd_dy(2);

	ColumnVector X(1); X.ReSize(1); X = 0.0;
	ColumnVector l(2); l.ReSize(2); l = 0.0;
	l(1) = dAugData(1);
	l(2) = dAugData(2);
	X << v.t() * (A.i() * l);

	cout << "V: " << endl;
	cout << fixed << setprecision(5) << setw(15) << v << endl;
	cout << "A: " << endl;
	cout << fixed << setprecision(5) <<setw(15)<< A << endl;
	cout << "X: " << endl;
	cout << fixed << setprecision(5) << setw(15) << X << endl;

	return X(1);
}

bool great::t_gInterpol::interpolAug(const t_gtime& nowT,  vector<t_gsatdata>& data_rover, vector<t_gsatdata>& data_base)
{
	t_gtime curT = nowT;
	_data_rover.erase(_data_rover.begin(), _data_rover.end());

	auto gobs = dynamic_cast<t_gallobs*>((*_alldata)[t_gdata::ALLOBS]);
	_data_rover = gobs->obs(_site, nowT);
	if (_data_rover.size() == 0) return false;
	map<string, t_gsatdata> rdata;

	if (!_SppProcOneEpoch(nowT, _site, &_crd_rover, NULL, &rdata)) return false;
	_crd_rover += _allobj->obj(_site)->eccxyz(nowT);

	/* base site setting for rtk */
	map<string, t_gsatdata> bdata;
	if (_isBase) {
		if (!_selMajorSite(nowT)) return false;
		_data_base.erase(_data_base.begin(), _data_base.end());
		_data_base = gobs->obs(_strMajorSite, nowT);
		if (_data_base.size() == 0) return false;
		//if (!_Proc->SppProcOneEpoch(nowT, _strMajorSite, &crd_major, &clk_major, &bdata)) return false;
	}

	int nSat_cor = 0;

	ColumnVector dx, dy, dz;
	vector<t_gtriple> crd_base;
	ColumnVector Aug_l, Aug_p;

	ostringstream os;
	os << "> " << left << setw(35) << nowT.str("%Y %m %d %H %M %p")
		<< setw(5) << "" << endl;

	for (auto iter = _data_rover.begin(); iter != _data_rover.end();)
	{
		auto sat = iter->sat();
		auto gs = iter->gsys();

		t_gsatdata obsdata_major;
		if (_isBase)
		{
			for (int j = 0; j < _data_base.size(); j++) {
				if (_data_base[j].sat() == sat) {
					obsdata_major = (_data_base)[j]; break;
				}
			}
			if (obsdata_major.sat() == "") {
				iter = _data_rover.erase(iter);
				continue;
			}
		}
		os << right << setw(3) << iter->sat();

		bool bcalcor = false;
		t_gobsgnss gnss(sat);
		vector < GOBSBAND > band = dynamic_cast<t_gsetgnss*>(_set)->band(gs);
		for (int i = 0; i < band.size(); i++)
		{
			GOBSBAND b = band[i];

			bool test = false;//自己改自己,lvhb added for PPP-RTK in 202009
			if (test)
			{
				auto aug = _allaug->getAug(_site, nowT, sat);
				if (!aug) continue;

				bool valid = false;
				if (aug->find(make_pair(AUGTYPE::TYPE_P, b)) != aug->end() &&
					aug->find(make_pair(AUGTYPE::TYPE_L, b)) != aug->end())
				{
					if (_isCompAug)	valid = true;
					else if (aug->find(make_pair(AUGTYPE::TYPE_ION, b)) != aug->end() &&
						aug->find(make_pair(AUGTYPE::TYPE_TRP, GOBSBAND::BAND_1)) != aug->end())
					{
						valid = true;
					}
				}
				auto gobsP = iter->select_range(b);
				auto gobsL = iter->select_phase(b);
				if (gobsP == X || gobsL == X) valid = false;
				if (valid) {
					auto Aug_p_rover = (*aug)[make_pair(AUGTYPE::TYPE_P, b)];
					auto Aug_l_rover = (*aug)[make_pair(AUGTYPE::TYPE_L, b)];

					//_allaug->add_data(_site, curT, sat, make_pair(AUGTYPE::TYPE_P, b), Aug_p_rover);
					//_allaug->add_data(_site, curT, sat, make_pair(AUGTYPE::TYPE_L, b), Aug_l_rover);
					//Aug_l_rover = Aug_p_rover = 0.0;

					double corrP = iter->obs_C(gobsP) - Aug_p_rover;
					double corrL = (iter->obs_L(gobsL) - Aug_l_rover) / gnss.wavelength(b);
					iter->addobs(gobsP, corrP);
					iter->addobs(gobsL, corrL);

					bcalcor = true;
				}
				else {
					if (!_isBase) {
						auto gobsP = iter->select_range(b, true);
						auto gobsL = iter->select_phase(b, true);
						iter->addobs(gobsP, 0);
						iter->addobs(gobsL, 0);
					}
				}
				continue;
			}

			_Rectangle.ReSize(4); _Rectangle = 0;

			int iobs = 0;
			dx.ReSize(_nipSite); dy.ReSize(_nipSite); dz.ReSize(_nipSite);;
			Aug_l.ReSize(_nipSite); Aug_p.ReSize(_nipSite);
			dx = 0; dy = 0; dz = 0; Aug_p = 0; Aug_l = 0; crd_base.clear();
			for (int j = 0; j < _nipSite; j++)
			{
				string site = _list_ipsite[j];
				auto aug = _allaug->getAug(site, nowT, sat);
				if (!aug) continue;

				auto crd = dynamic_cast<t_gsetrec*>(_set)->get_crd_xyz(site);
				crd += _allobj->obj(site)->eccxyz(nowT);

				bool valid = false;
				if (aug->find(make_pair(AUGTYPE::TYPE_P, b)) != aug->end() &&
					aug->find(make_pair(AUGTYPE::TYPE_L, b)) != aug->end())
				{
					if (_isCompAug)	valid = true;
					else if (aug->find(make_pair(AUGTYPE::TYPE_ION, b)) != aug->end() &&
						aug->find(make_pair(AUGTYPE::TYPE_TRP, GOBSBAND::BAND_1)) != aug->end())
					{
						valid = true;
					}
				}
				if (valid) {
					++iobs;
					dx(iobs) = -_crd_rover[0] + crd[0];
					dy(iobs) = -_crd_rover[1] + crd[1];
					dz(iobs) = -_crd_rover[2] + crd[2];
					Aug_p(iobs) = (*aug)[make_pair(AUGTYPE::TYPE_P, b)];
					Aug_l(iobs) = (*aug)[make_pair(AUGTYPE::TYPE_L, b)];
					crd_base.push_back(crd);
					//_isInterpolRange(&crd);
				}
			}//end nipsite

			//interpolation
			if (iobs > 2 /*&& _isInterpolRange()*/)
			{
				if (iobs > 3 && 0) {//select 3 sites by Delaunay method
					t_gDelaunay delaunay;
					int index[3] = { 0 };
					delaunay.getIntSite(_crd_rover, crd_base, index);
					ColumnVector dx_new, dy_new, dz_new, Aug_p_new, Aug_l_new;
					dx_new = dx; dy_new = dy; dz_new = dz; Aug_p_new = Aug_p; Aug_l_new = Aug_l;
					for (auto i = 1; i <= 3; i++) {
						dx.Row(i) = dx_new.Row(index[i - 1]);
						dy.Row(i) = dy_new.Row(index[i - 1]);
						dz.Row(i) = dz_new.Row(index[i - 1]);
						Aug_p.Row(i) = Aug_p_new.Row(index[i - 1]);
						Aug_l.Row(i) = Aug_l_new.Row(index[i - 1]);
					}
					iobs = 3;
				}
				dx.Rows(1, iobs); dy.Rows(1, iobs); dz.Rows(1, iobs);
				Aug_p.Rows(1, iobs); Aug_l.Rows(1, iobs);
				/*t_spt_gInterpol interpol = make_shared<t_gInterpol>(dx, dy, dz);*/
				_setCrd(dx, dy, dz);
				auto Aug_p_rover = _MLCM(Aug_p);
				auto Aug_l_rover = _MLCM(Aug_l);

				os << right
					<< setw(12) << fixed << setprecision(4) << Aug_p_rover
					<< setw(12) << fixed << setprecision(4) << Aug_l_rover;

				auto gobsP = iter->select_range(b);
				auto gobsL = iter->select_phase(b);
				if (gobsP == X || gobsL == X) continue;
				//add obs
				if (_isBase)//vrs
				{
					if (obsdata_major.select_phase(b) == X || obsdata_major.select_range(b) == X
						|| double_eq(obsdata_major.obs_L(gobsL), 0.0) || double_eq(obsdata_major.obs_C(gobsP), 0.0)) {
						continue;
					}
					//creating VRS observations- referenced satellite
					if (!bcalcor)
					{
						//major site coordinate
						auto xyz_r = dynamic_cast<t_gsetrec*>(_set)->get_crd_xyz(_strMajorSite);
						xyz_r += _allobj->obj(_strMajorSite)->eccxyz(nowT);
						_cal_sat_inf(obsdata_major, xyz_r);//major site
						_cal_sat_inf(*iter, _crd_rover);//rover site
					}

					if (double_eq(obsdata_major.rho(), 0.0) || double_eq(iter->rho(), 0.0))	continue;
					double rho_AV_sat = iter->rho() - obsdata_major.rho();

					t_spt_gobs obs_vrs = _check_vrs_prn_pt(nowT, sat);

					auto aug_major = _allaug->getAug(_strMajorSite, nowT, sat);
					if (!aug_major) continue;
					if (aug_major->find(make_pair(AUGTYPE::TYPE_P, b)) == aug_major->end() &&
						aug_major->find(make_pair(AUGTYPE::TYPE_L, b)) == aug_major->end())
					{
						continue;
					}
					auto Aug_p_major = (*aug_major)[make_pair(AUGTYPE::TYPE_P, b)];
					auto Aug_l_major = (*aug_major)[make_pair(AUGTYPE::TYPE_L, b)];

					//new method by lvhb,202010
					double trop_satb = _getTropModValue(nowT, obsdata_major);
					double trop_satr = _getTropModValue(nowT, *iter);
					double tmpddtrop = (trop_satr - trop_satb);
					Aug_p_rover += tmpddtrop;
					Aug_l_rover += tmpddtrop;

					double obsP = obsdata_major.obs_C(gobsP) + rho_AV_sat + Aug_p_rover - Aug_p_major /*- clk_major*/;
					double obsL = (obsdata_major.obs_L(gobsL) + rho_AV_sat + Aug_l_rover - Aug_l_major /*- clk_major - Amb_major*/) / obs_vrs->wavelength(b);

					obs_vrs->addobs(gobsP, obsP);
					obs_vrs->addobs(gobsL, obsL);

					gobs->addobs(obs_vrs);
				}
				else {
					_allaug->add_data(_site, curT, sat, make_pair(AUGTYPE::TYPE_P, b), Aug_p_rover);
					_allaug->add_data(_site, curT, sat, make_pair(AUGTYPE::TYPE_L, b), Aug_l_rover);

					double corrP = (iter->obs_C(gobsP) - Aug_p_rover);
					double corrL = (iter->obs_L(gobsL) - Aug_l_rover) / gnss.wavelength(b);
					iter->addobs(gobsP, corrP);
					iter->addobs(gobsL, corrL);
				}
				bcalcor = true;
			}
			else {
				if (!_isBase) {
					auto gobsP = iter->select_range(b, true);
					auto gobsL = iter->select_phase(b, true);
					iter->addobs(gobsP, 0);
					iter->addobs(gobsL, 0);
				}
			}
		}//end band

		os << endl;

		if (bcalcor) nSat_cor++;
		else if (!_isBase) {
			iter = _data_rover.erase(iter);
			continue;
		}
		iter++;
	}//end data one epoch
	if (nSat_cor < 5) return false;
	else
	{
		if (_augfile)
		{
			_augfile->write(os.str().c_str(), os.str().size());
			_augfile->flush();
		}
	}
	//output lvhb
	data_rover.erase(data_rover.begin(), data_rover.end());
	data_rover = _data_rover;
	if (_isBase) {
		_data_base.erase(_data_base.begin(), _data_base.end());
		_data_base = gobs->obs(_site_vrs, nowT);
		if (_data_base.size() < 6) return false;
		//output lvhb
		data_base.erase(data_base.begin(), data_base.end());
		data_base = _data_base;
	}

	return true;
}

void great::t_gInterpol::printAugFile(t_giof* augfile)
{
	if (augfile) _augfile = augfile;
}

bool great::t_gInterpol::_selMajorSite(const t_gtime& nowT)
{
	//t_gtime beg = dynamic_cast<t_gsetgen*>(_set)->beg();
	t_gtime end = dynamic_cast<t_gsetgen*>(_set)->end();
	auto grec_vrs = _allobj->obj(_site_vrs);
	//auto xyz = _crd_rover + _allobj->obj(_site)->eccxyz(nowT);
	grec_vrs->crd(_crd_rover, t_gtriple(1, 1, 1), nowT, end + 30, true);

	double preLengthBL = 0.0; _strMajorSite = "";
	for (int j = 0; j < _nipSite; j++)
	{
		string site = _list_ipsite[j];
		if (!_allaug->getAug(site, nowT)) continue;

		auto crd = dynamic_cast<t_gsetrec*>(_set)->get_crd_xyz(site);
		crd[0] -= _crd_rover[0];
		crd[1] -= _crd_rover[1];
		crd[2] -= _crd_rover[2];
		double lengthBL = crd.norm();
		if (_strMajorSite == "") {
			_strMajorSite = site;
			preLengthBL = lengthBL;
		}
		else if (lengthBL < preLengthBL)
		{
			preLengthBL = lengthBL;
			_strMajorSite = site;
		}
	}
	if (_strMajorSite == "") return false;

	string ant = _allobj->obj(_strMajorSite)->ant(nowT);
	grec_vrs->ant(ant, nowT);
	_allobj->sync_pcvs();

	if (dynamic_cast<t_gsetproc*>(_set)->realtime()) {
		while (!_avaliable(nowT, &_site, &_strMajorSite))
			t_gtime::gmsleep(1);
	}

	return true;
}

t_spt_gobs great::t_gInterpol::_check_vrs_prn_pt(const t_gtime& nowT, string& sat)
{
	auto gobs = dynamic_cast<t_gallobs*>((*_alldata)[t_gdata::ALLOBS]);
	auto obs_temp = gobs->obs_pt(_site_vrs, nowT);
	t_spt_gobs obs_ref(nullptr);
	for (auto iter = obs_temp.begin(); iter != obs_temp.end(); iter++) {
		if ((*iter)->sat() == sat) {
			obs_ref = *iter;
			break;
		}
	}
	if (!obs_ref) {
		obs_ref = make_shared<t_gobsgnss>(_site_vrs, sat, nowT);
	}
	return obs_ref;
}

bool great::t_gInterpol::_SppProcOneEpoch(const t_gtime& nowT, const string& site, t_gtriple* crd, double* reclk, map<string, t_gsatdata>* satdata)
{
	auto gobs = dynamic_cast<t_gallobs*>((*_alldata)[t_gdata::ALLOBS]);
	auto gorb = dynamic_cast<t_gallnav*>((*_alldata)[t_gdata::GRP_EPHEM]);
	auto gobj = dynamic_cast<t_gallobj*>((*_alldata)[t_gdata::ALLOBJ]);
	auto gbias = dynamic_cast<t_gallbias*>((*_alldata)[t_gdata::ALLBIAS]);
	shared_ptr<t_gspplsq> spplsq = make_shared<t_gspplsq>(site, _set, gobj);
	spplsq->setDAT(gobs, gorb);
	spplsq->setDCB(gbias);
	if (spplsq->processBatch(nowT, nowT) <= 0) {
		if (_log) {
			_log->comment(0, site + "spplsq failed");
		}
		return false;
	}
	if (crd) *crd = spplsq->getCrd(nowT);
	if (reclk)*reclk = spplsq->getRecClk(nowT);
	if (satdata) spplsq->outcomsat(*satdata);
	return true;
}

void great::t_gInterpol::cal_rho_azel(const string& site_name, t_gtriple& xyz_s, const t_gtriple& xyz_r, t_gsatdata& obs_sat)
{
	//t_spt_gobs pvtflt = make_shared<t_gpvtflt>(_site, site_base, set, data);
}

void great::t_gInterpol::_cal_sat_inf(t_gsatdata& obs_sat, const t_gtriple& xyz_r)
{
	auto precise_bias = make_shared<t_gprecisebias>(_alldata, _log, _set);
	auto gmodel = dynamic_cast<t_gprecisemodel*>(precise_bias->precisemodel());
	t_gallpar parm;

	int ipar = 0;
	string tmpsite = obs_sat.site();
	//xyz
	t_gpar parx(tmpsite, t_gpar::CRD_X, ++ipar, ""); parx.value(xyz_r[0]);
	t_gpar pary(tmpsite, t_gpar::CRD_Y, ++ipar, ""); pary.value(xyz_r[1]);
	t_gpar parz(tmpsite, t_gpar::CRD_Z, ++ipar, ""); parz.value(xyz_r[2]);
	parm.addParam(parx);
	parm.addParam(pary);
	parm.addParam(parz);
	//clk
	double clk = 0.0;
	//if (!SppProcOneEpoch(obs_sat.epoch(), tmpsite, NULL, &clk)) clk = 0.0;
	t_gpar parclk(tmpsite, t_gpar::CLK, ++ipar, ""); parclk.value(clk);
	parm.addParam(parclk);

	gmodel->_update_obs_info(obs_sat.epoch(), obs_sat, parm);
	gmodel->_prepare_obs(obs_sat.epoch(), parm);
	gmodel->_update_obs_info(obs_sat);
}

double great::t_gInterpol::_getTropModValue(const t_gtime& nowT, t_gsatdata& obsdata)
{
	double delay = 0.0;
	auto strEst = dynamic_cast<t_gsetgen*>(_set)->estimator();
	auto isFLT = (_isBase || strEst == "FLT");

	//saas
	auto trpModStr = dynamic_cast<t_gsetproc*>(_set)->tropo_model();
	shared_ptr<t_gtropo> tropoModel;
	if (trpModStr == SAASTAMOINEN)		 tropoModel = make_shared<t_saast>();
	else if (trpModStr == DAVIS)         tropoModel = make_shared<t_davis>();
	else if (trpModStr == HOPFIELD)      tropoModel = make_shared<t_hopf>();
	//else if (trpModStr == MOPS)          tropoModel = make_shared<t_blindmops>();
	else return 0.0;

	string site = obsdata.site();
	auto grec = dynamic_cast<t_gallobj*>((*_alldata)[t_gdata::ALLOBJ])->obj(site);
	// Cartesian coordinates to ellipsodial coordinates
	t_gtriple xyz, ell;
	xyz = dynamic_cast<t_gsetrec*>(_set)->get_crd_xyz(site);
	if (xyz.zero()) {
		xyz = grec->crd_arp(nowT);
	}
	else {
		xyz = xyz + grec->eccxyz(nowT);
	}
	xyz2ell(xyz, ell, false);

	double zwd = 0.0, zhd = 0.0;
	zwd = tropoModel->getZWD(ell, nowT);
	zhd = tropoModel->getZHD(ell, nowT);

	string sat = obsdata.sat();
	double ele = obsdata.ele();
	if (double_eq(ele, 0.0)) {
		if (!isFLT) {
			xyz = xyz - grec->eccxyz(nowT);
		}
		_cal_sat_inf(obsdata, xyz);
		ele = obsdata.ele();
	}

	auto tropo_mf = dynamic_cast<t_gsetproc*>(_set)->tropo_mf();

	double mfh, mfw, dmfh, dmfw;
	mfh = mfw = dmfh = dmfw = 0.0;
	if (tropo_mf == GMF || tropo_mf == NO_MF)
	{
		t_gmf mf;
		mf.gmf(nowT.mjd(), ell[0], ell[1], ell[2], G_PI / 2.0 - ele,
			mfh, mfw, dmfh, dmfw);
	}
	else if (tropo_mf == COSZ)
	{
		mfh = mfw = 1 / sin(ele);
		dmfh = dmfw = -(cos(ele)) / (sin(ele) * sin(ele));
	}

	delay = mfh * zhd + mfw * zwd;

	return delay;
}

bool great::t_gInterpol::_avaliable(t_gtime nowT, string* site, string* site_base)
{
	auto gobs = dynamic_cast<t_gallobs*>((*_alldata)[t_gdata::ALLOBS]);
	auto gorb = dynamic_cast<t_gallnav*>((*_alldata)[t_gdata::GRP_EPHEM]);
	//auto nppmodel = dynamic_cast<t_gsetproc*>(_set)->npp_model();
	bool isClient = dynamic_cast<t_gsetgen*>(_set)->isClient();

	bool b = nowT < gobs->end_obs(site ? (*site) : _site);
	if (_isBase && b) {
		if (site_base) {
			b = nowT < gobs->end_obs((*site_base));
		}
		else b = false;
	}

	return b;
}

bool great::t_gInterpol::_isInterpolRange(const t_gtriple* crdT)
{
	//_Rectangle: 1 : x_max 2 : x_min; 3:y_max; 4:y_min
		// To jude _crd_rover in rectangle
	if (crdT == NULL)
	{
		cerr << endl;
		cerr << _Rectangle << endl;
		cerr << _crd_rover[0] << " " << _crd_rover[1] << endl;
		if (_crd_rover[0] > _Rectangle(1) || _crd_rover[0] < _Rectangle(2)) return false;	//x
		if (_crd_rover[1] > _Rectangle(3) || _crd_rover[1] < _Rectangle(4)) return false;	//y

		return true;
	}

	auto crd = *crdT;

	if (double_eq(_Rectangle(1), 0.0)) {
		_Rectangle(1) = crd[0];//x_max
		_Rectangle(2) = 0;//x_min
		_Rectangle(3) = crd[1];//y_max
		_Rectangle(4) = 0;//y_min
	}
	else if (double_eq(_Rectangle(2), 0.0)) {
		//X
		if (_Rectangle(1) < crd[0]) {
			_Rectangle(2) = _Rectangle(1);
			_Rectangle(1) = crd[0];
		}
		else _Rectangle(2) = crd[0];
		//Y
		if (_Rectangle(3) < crd[1]) {
			_Rectangle(4) = _Rectangle(3);
			_Rectangle(3) = crd[1];
		}
		else _Rectangle(4) = crd[1];

	}
	else {
		//x
		if (_Rectangle(1) < crd[0]) _Rectangle(1) = crd[0];
		else if (_Rectangle(2) > crd[0]) _Rectangle(2) = crd[0];
		//y
		if (_Rectangle(3) < crd[1]) _Rectangle(3) = crd[1];
		else if (_Rectangle(4) > crd[1]) _Rectangle(4) = crd[1];
	}

	return true;
}

double great::t_gInterpol::_corrTropHeight(const t_gtime& nowT, int sign, const t_gtriple& bcrd, t_gtriple& rcrd, const t_gtriple& crd, t_gsatdata bdata_sat, t_gsatdata bdata_bsat, t_gsatdata rdata_sat, t_gsatdata rdata_bsat)
{
	double delay = 0.0;
	//saas
	auto trpModStr = dynamic_cast<t_gsetproc*>(_set)->tropo_model();
	shared_ptr<t_gtropo> tropoModel;
	if (trpModStr == SAASTAMOINEN)		 tropoModel = make_shared<t_saast>();
	else if (trpModStr == DAVIS)         tropoModel = make_shared<t_davis>();
	else if (trpModStr == HOPFIELD)      tropoModel = make_shared<t_hopf>();
	//else if (trpModStr == MOPS)          tropoModel = make_shared<t_blindmops>();
	else return false;
	auto tropo_mf = dynamic_cast<t_gsetproc*>(_set)->tropo_mf();

	string bsat = bdata_bsat.sat();
	string sat = bdata_sat.sat();


	t_gtriple bell, rell, ell_h0;
	xyz2ell(bcrd, bell, false);
	xyz2ell(rcrd, rell, false);
	xyz2ell(crd, ell_h0, false);

	double zwd = tropoModel->getZWD(rell, nowT);
	double zhd = tropoModel->getZHD(rell, nowT);

	bell[2] = ell_h0[2];
	rell[2] = ell_h0[2];
	double zwd_h0 = tropoModel->getZWD(rell, nowT) - zwd;
	double zhd_h0 = tropoModel->getZHD(rell, nowT) - zhd;


	t_gtriple xyz_r, xyz_b;
	ell2xyz(rell, xyz_r, false);
	ell2xyz(bell, xyz_b, false);
	t_gtriple rb_satcrd = rdata_bsat.satcrd(), r_satcrd = rdata_sat.satcrd(),
		bb_satcrd = bdata_bsat.satcrd(), b_satcrd = bdata_sat.satcrd();
	cal_rho_azel(rdata_bsat.site(), rb_satcrd, xyz_r, rdata_bsat);
	cal_rho_azel(rdata_sat.site(), r_satcrd, xyz_r, rdata_sat);
	cal_rho_azel(bdata_bsat.site(), bb_satcrd, xyz_b, bdata_bsat);
	cal_rho_azel(bdata_sat.site(), b_satcrd, xyz_b, bdata_sat);

	double bbsat_ele = bdata_bsat.ele();
	//double bsat_ele = bdata_sat.ele();
	//double rbsat_ele = rdata_bsat.ele();
	double rsat_ele = rdata_sat.ele();

	//auto tropo_mf = dynamic_cast<t_gsetproc*>(_set)->tropo_mf();

	for (int j = 0; j < 2; j++)
	{
		auto tmpell = bell;
		auto tmpele = bbsat_ele;
		if (j)
		{
			tmpell = rell;
			tmpele = rsat_ele;
		}
		double mfh, mfw, dmfh, dmfw;
		mfh = mfw = dmfh = dmfw = 0.0;
		if (tropo_mf == GMF)
		{
			t_gmf mf;
			mf.gmf(nowT.mjd(), tmpell[0], tmpell[1], tmpell[2], G_PI / 2.0 - tmpele,
				mfh, mfw, dmfh, dmfw);
		}
		else if (tropo_mf == COSZ)
		{
			mfh = mfw = 1 / sin(tmpele);
			dmfh = dmfw = -(cos(tmpele)) / (sin(tmpele) * sin(tmpele));
		}

		delay = (mfh * zhd_h0 + mfw * zwd_h0) - delay;
	}

	return delay;
}
