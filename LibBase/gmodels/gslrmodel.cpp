#include "gmodels/gslrmodel.h"
#include "gutils/gtime.h"
#include "gutils/ginfolog.h"
#include "gutils/gleo.h"
#include <algorithm>

// #define DEBUG_LOG
// #define DEBUG_GLO
namespace great
{
	t_gslrmodel::t_gslrmodel(t_gallproc* data, t_glog* log, t_gsetbase* setting) :
		t_gpppmodel("", log, setting),
		_tide(0),
		_allerp(0),
		_trs2crs_2000(nullptr),
		_scf2crs_leo(nullptr)
	{
		_tide = new t_gtideIERS(log);

		setOTL(dynamic_cast<t_gallotl*>((*data)[t_gdata::ALLOTL]));
		setOBJ(dynamic_cast<t_gallobj*>((*data)[t_gdata::ALLOBJ]));
		_allorb = dynamic_cast<t_gallorb*>((*data)[t_gdata::ALLORB]);
		_allerp = dynamic_cast<t_gpoleut1*>((*data)[t_gdata::ALLPOLEUT1]);
		_allde = dynamic_cast<t_gnavde*>((*data)[t_gdata::ALLDE]);
		_allatt = dynamic_cast<t_gallatt*>((*data)[t_gdata::AllATTITUDE]);
		_alllra= dynamic_cast<t_galllra*>((*data)[t_gdata::ALLLRA]);
		_allslrf= dynamic_cast<t_gslrf*>((*data)[t_gdata::SLRF]);
		_allpreproc= dynamic_cast<t_gallpreproc*>((*data)[t_gdata::ALLPREPROC]);
		_allrb = dynamic_cast<t_gallrb*>((*data)[t_gdata::ALLRB]);
		if (_allrb == nullptr)
		{
			_allrb = new t_gallrb();
		}

		_trs2crs_2000 = new t_gtrs2crs(false, _allerp);

	}

	t_gslrmodel::~t_gslrmodel()
	{
		if (_tide)
		{
			delete _tide;
			_tide = nullptr;
		}
	}

	void t_gslrmodel::setOTL(t_gallotl* otl)
	{
		if (_tide)
		{
			_tide->setOTL(otl);
		}
	}


	bool t_gslrmodel::_cmb_equ(const t_gtime & epoch, t_gallpar& pars, string sat_leo,t_glsqEquationMatrix& equ,string slr_site)
	{

		equ.clear_allequ();

		_sat_leo = sat_leo;

		_crt_rec = _sat_leo;

		vector<pair<int, double> > coef;

		double omc = 0.0;
		double wgt = 0.0;

		//0.01(1-3)  0.005(4-6)  0.001(7-9)
		//0.1(10-12)  1(13-15)
		double sigmaSLR = 0.01;
		wgt = 1.0 / pow(sigmaSLR, 2);

		if (!_omc_obs(epoch, omc) ||
			!_prt_obs(epoch, pars, coef,slr_site,strs2strl(_sat_leo)))
		{
			return false;
		}

		t_gpreprocslr  preproc_slr;
		preproc_slr = _allpreproc->getSLROBS(strs2strl(_sat_leo),epoch);

		//equ.add_equ(coef, wgt, omc, _crt_rec, "", t_gobscombtype());
		t_gobscombtype type("SLR");
		equ.add_equ(coef, wgt, omc, preproc_slr._sta_name, _crt_rec, type, false);
		
		/*equ.add_equ(B, P, l, vector<string>({ "X","Y","Z" }),
			vector<string>(3, *iter),
			vector<t_gobscombtype>(3, t_gobscombtype()));*/

	}

	bool t_gslrmodel::apply_erp(const t_gtime & epoch,t_gallpar& pars, string site)
	{

		ColumnVector xyz_trs(0.0, 3);
		xyz_trs.ReSize(3, 1);
		int site_indexX = pars.getParam(site, t_gpar::CRD_X, "");
		int site_indexY = pars.getParam(site, t_gpar::CRD_Y, "");
		int site_indexZ = pars.getParam(site, t_gpar::CRD_Z, "");
		
		xyz_trs(1) = pars.getPar(site_indexX).value();
		xyz_trs(2) = pars.getPar(site_indexY).value();
		xyz_trs(3) = pars.getPar(site_indexZ).value();

		_trs2crs_2000->calcRotMat(epoch, true, true, true);
	    _siteRotmat = _trs2crs_2000->getRotMat();
		_drdxpole = _trs2crs_2000->getMatDx() * xyz_trs;
		_drdypole = _trs2crs_2000->getMatDy() * xyz_trs;
		_drdut1 = _trs2crs_2000->getMatDu() * xyz_trs;

		int eop_index = pars.getParam("", t_gpar::XPOLE, "");
		if (eop_index >= 0) {

			// get unit_vector in TRS
			ColumnVector unit_rec2sat(3);

			unit_rec2sat = _crs_rec_crd.crd_cvect() - _crs_sta_crd.crd_cvect();
			unit_rec2sat = unit_rec2sat / unit_rec2sat.NormFrobenius();

			// get drotdxpole
			Matrix drdxpole = _drdxpole;
			drdxpole = (drdxpole / RAD2SEC);
			// get drotdypole
			Matrix drdypole = _drdypole;
			drdypole = (drdypole / RAD2SEC);
			// get drotdut1
			Matrix drdut1 = _drdut1;
			drdut1 = drdut1 / RAD2TSEC;
			double diff = (epoch - pars[eop_index].beg) / 86400.0;

			// dxdxpole dxdypole dxdut1
			Matrix dxdpoleut1(3, 6);
			// dx/xpole
			dxdpoleut1.Column(1) = drdxpole;
			// dx/dxpole
			dxdpoleut1.Column(2) = drdxpole * diff;
			// dx/ypole
			dxdpoleut1.Column(3) = drdypole;
			// dx/dypole
			dxdpoleut1.Column(4) = drdypole * diff;
			// dx/ut1
			dxdpoleut1.Column(5) = drdut1;
			// dx/dut1
			dxdpoleut1.Column(6) = drdut1 * diff;

			_drderp.ReSize(1, 6);
			_drderp = unit_rec2sat.t() * dxdpoleut1;
		}
		return true;
	}

	bool t_gslrmodel::_prt_obs(const t_gtime & epoch, t_gallpar& pars, vector<pair<int, double>>& coeff, string slr_site ,string sat)
	{
		/*t_gtriple groundEll;
		xyz2ell(_trs_rec_crd, groundEll, false);*/
		int parNum = pars.parNumber();

		apply_erp(epoch, pars, slr_site);

		vector<t_gpar> orbpars = pars.getOrbPar(strs2strl(_crt_rec));
		Matrix partial(0.0, orbpars.size(), 3);
		// TODO
		_allorb->getLeoPartial(epoch+_two_way/2.0, strs2strl(_crt_rec), partial, orbpars);
		_partialSite = partial;

		/*t_gtriple xyz_leo(0.0, 0.0, 0.0);
		if (_allorb->get_crs_sat_crd(epoch, strs2strl(_site), xyz_leo[0], xyz_leo[1], xyz_leo[2])) {

			_crs_rec_crd = xyz_leo;
		}*/

		
		// set sat parameter PXSAT in lsq
		int site_index = pars.getParam(_crt_rec, t_gpar::PXSAT, "");
		//obsdata.addsatindex(par[sat_index].index);
		if (site_index >= 0)
		{
			_siteIndex = pars.getParIndex(site_index);
		}

		for (int ipar = 0; ipar < parNum; ipar++)
		{
			
			double coeff_value = _Partial(epoch,  pars.getPar(ipar),slr_site,sat);
			if (coeff_value != 0.0)
			{
				//note: This is bad from 1 in coeff vector
				coeff.push_back(make_pair(ipar + 1, coeff_value));
			}
		}

		return true;
	}

	double t_gslrmodel::_Partial(const t_gtime& epoch,  t_gpar& par, string slr_site,string slr_sat)
	{
		double mfw, dmfw, mfh, dmfh;
		mfw = dmfw = mfh = dmfh = 0.0;
		t_gtriple ell(0.0, 0.0, 0.0);
		ColumnVector unit_rec2sat(3);
		unit_rec2sat = -_crs_rec_crd.crd_cvect() + _crs_sta_crd.crd_cvect();
		unit_rec2sat = unit_rec2sat / unit_rec2sat.NormFrobenius();

		switch (par.parType)
		{
		case t_gpar::CRD_X:
		{
			if (slr_site == par.site)
			{
				return (par.value() - _trs_rec_crd.crd(0)) / _slr_rho;
			}
		}
		case t_gpar::CRD_Y:
		{
			if (slr_site == par.site)
			{
				return  (par.value() - _trs_rec_crd.crd(1)) / _slr_rho;
			}
		}
		case t_gpar::CRD_Z:
		{
			if (slr_site == par.site)
			{
				return  (par.value() - _trs_rec_crd.crd(2)) / _slr_rho;
			}
		}
		case t_gpar::CLK:
		{
			return 0.0;
		}
		case t_gpar::RB:
		{
			if (slr_site == par.site && slr_sat==par.prn)
			{
				return 1.0;
			}
			
		}
		}

		if (par.parType >= t_gpar::XPOLE && par.parType <= t_gpar::DUT1) {

			if (epoch > par.beg&&epoch < par.end)
			{
				return _Partial_Erp(par);
			}
		}
		else if (par.parType >= t_gpar::PXSAT && par.parType <= t_gpar::EMP_Br) {

			if (_crt_rec == par.site)
			{
				return _Partial_Orb_Leo(par, epoch);
			}
			

		}

		return 0.0;
		
	}

	double t_gslrmodel::_Partial_Erp(t_gpar& par)
	{

		switch (par.parType)
		{
		case t_gpar::XPOLE:
			return -_drderp(1, 1);
			break;
		case t_gpar::YPOLE:
			return -_drderp(1, 3);
			break;
		case t_gpar::DXPOLE:
			return -_drderp(1, 2);
			break;
		case t_gpar::DYPOLE:
			return -_drderp(1, 4);
			break;
		case t_gpar::UT1:
			return _drderp(1, 5);
			break;
		case t_gpar::DUT1:
			return _drderp(1, 6);
			break;
		default:
			return 0.0;
			break;
		}
	}

	double t_gslrmodel::_Partial_Orb_Leo(t_gpar& par,t_gtime epoch)
	{
		
		if (par.beg > epoch && par.lPWC)
		{
			return 0.0;
		}

		//Matrix funct = obsdata.orbfunct();
		Matrix funct = _partialSite;
		int index_begin = _siteIndex;

		// get unit_vector in crs
		ColumnVector unit_rec2sat(3);

		ColumnVector x_f(0.0, 3), v_f(0.0, 3), x_j(0.0, 3), v_j(0.0, 3);

		// get
		//unit_rec2sat = -obsdata.sat2reccrs().crd_cvect();
		unit_rec2sat =_crs_rec_crd.crd_cvect()- _crs_sta_crd.crd_cvect();
		unit_rec2sat = unit_rec2sat / unit_rec2sat.NormFrobenius();

		/*cout << endl;
		cout << satdata.epoch().mjd()<<"  "<< satdata.epoch().sod() <<"   u2s  "<<setw(30)<<setprecision(15)<<scientific<< unit_rec2sat << endl;*/
		// order: P X0 Y0 Z0 VX0..VZ0 solarPar
		funct = funct * unit_rec2sat;
		funct = funct * 1E3;

		//cout << setw(20) << setprecision(13)<< scientific << funct << endl;
		return funct(par.index - index_begin + 1, 1);

	}


	bool t_gslrmodel::_omc_obs(const t_gtime & epoch, double& omc)
	{
		int jd;
		double sod, sodr, sodb;
		double fjd;
		t_gtime now, rtime, btime, time;
		int ista;
		Matrix rot_l2f;
		Matrix rotmat;
		ColumnVector xyz(0.0, 3), xtmp;
		t_gtriple xsun, vsun, xlun, vlun;
		t_gtriple xsat, vsat;
		ColumnVector x_f(0.0, 3), v_f(0.0, 3), x_j(0.0, 3), v_j(0.0, 3), r1(0.0, 3), dloudx(0.0, 3), u_r1(0.0, 3);
		double _xsat[3], _vsat[3];
		ColumnVector cxsat(0.0, 3), cvsat(0.0, 3);
		ColumnVector cxsat_km(0.0, 3), cvsat_km(0.0, 3);
		ColumnVector cxsun(0.0, 3), cxsun_km(0.0, 3);
		double dlra_f = 0.0, dlra = 0.0;
		double r1length = 0.0;
		double azim_leo, elev_leo;
		double azim, elev;
		double trpdel, ztdpart;
		double temp;
		double press_wvp;
		double fcul_ztd, fcul_zhd, fcul_zwd;
		double fak_map;
		double std;
		double ztd_dif;
		double reldel0, reldel1;
		double sitrad, satrad;
		Matrix xtmp_leo1;
		ColumnVector xtmp_leo(0.0, 3), xtmp_lrr(0.0, 3);
		xyz.ReSize(3, 1);
		xtmp.ReSize(3, 1);

		rotmat.ReSize(3, 3);

		xtmp_leo1.ReSize(3, 1);

		//now.tsys(t_gtime::UTC);
		now = epoch;
		now.tsys(t_gtime::GPS);

		t_gpreprocslr  preproc_slr;
		preproc_slr = _allpreproc->getSLROBS(strs2strl(_sat_leo),now);

		Matrix rot_scf2crsleo(0.0, 3, 3); //rotation matrix of SCF to CRS, just for leo!
		Matrix rot_scf2crs(0.0, 3, 3);
		Matrix rot_trs2crs(0.0, 3, 3);   // rotation matrix of TRS to CRS
		Matrix rot_du(0.0, 3, 3);
		
		_machine_offset(xyz);
		
		_scf2crs_leo = new t_gscf2crsleo(strs2strl(_sat_leo));

		map<string, map<string,double>> RB = _allrb->getrb();
		double range_bias;

		if (RB.empty())
		{
			range_bias = 0.0;
		}
		else
		{
			range_bias = (RB)[strs2strl(_sat_leo)][preproc_slr._sta_name];
		}

		//t_gprecisemodel gmodel=new()

		jd = epoch.mjd();

		//now.tsys(t_gtime::GPS);

		sod = now.sod() + now.dsec();

		sodr = preproc_slr._two_way_time / 2.0 + sod;

		sodb = preproc_slr._two_way_time + sod;

		_two_way = preproc_slr._two_way_time;

		for (ista = 0; ista < _allslrf->getSize(); ista++)
		{
			if (preproc_slr._code == _allslrf->getCode(ista))
			{
				break;
			}
		}
		if (ista >= _allslrf->getSize())
		{
			write_log_info(_log, 1, "NOTE", "The code is not found in SLRF2008");
			return false;
		}

		rot_l2f.ReSize(3, 3);

		enu2trs(_allslrf->getOnegslrf(ista)._blh[0], _allslrf->getOnegslrf(ista)._blh[1], rot_l2f);

		fjd = jd + (sod + 32.1840 + 19.00) / 86400.0;

		_allde->get_pv(fjd, "SUN", xsun, vsun);
		_allde->get_pv(fjd, "MOON", xlun, vlun);

		for (int ij = 0; ij < 3; ij++)
		{
			cxsun(ij + 1) = xsun.crd(ij);

		}

		rtime.from_mjd(jd, (int)sodr, sodr - (int)sodr);

		//m
		_allorb->get_crs_sat_crd(rtime, strs2strl(_sat_leo), _xsat[0], _xsat[1], _xsat[2]);
		xsat.set(_xsat);

		//m/s
		_allorb->get_crs_sat_vel(rtime, strs2strl(_sat_leo), _vsat[0], _vsat[1], _vsat[2]);
		vsat.set(_vsat);
		
		for (int ij = 0; ij < 3; ij++)
		{
			cxsat(ij + 1) = _xsat[ij];
			cvsat(ij + 1) = _vsat[ij];
		}

		t_gtime tdt_grtime = rtime;
		tdt_grtime.tsys(t_gtime::TT);
		double tdt = tdt_grtime.dmjd();

		_trs2crs_2000->calcRotMat(tdt_grtime, true, true, true);
		rot_trs2crs = _trs2crs_2000->getRotMat();

		//scf2crs
		_scf2crs_leo->calcRotMat(rtime, _allatt->attitude(strs2strl(_sat_leo)), cxsat, cvsat, rot_trs2crs, rot_scf2crsleo);

		xtmp = rot_scf2crsleo * xyz;

		for (int m = 0; m < 3; m++)
		{
			cxsat(m + 1) = xsat.crd(m) + xtmp(m + 1);
		}

		_crs_rec_crd.set(cxsat);

		// SET TRS in epoch TR [include earth rotation]
		_update_rot_matrix(epoch);
		_trs_rec_crd = t_gtriple(_trs2crs_2000->getRotMat().t() * _crs_rec_crd.crd_cvect());

		for (int m = 0; m < 3; m++)
		{
			double ttmp = _allslrf->getOnegslrf(ista)._xsta[m] + _allslrf->getOnegslrf(ista)._offset[m] * 0.001;
			x_f(m + 1) = _allslrf->getOnegslrf(ista)._xsta[m] + _allslrf->getOnegslrf(ista)._offset[m] * 0.001;
			x_f(m + 1) = x_f(m + 1)*1000.0;
			v_f(m + 1) = 0.0;
		}

		time.from_mjd(jd, (int)sod, (sod - (int)sod));

		_calculate_sta_pos(time, &x_f, &v_f, &x_j, &v_j);
		
		_crs_sta_crd.set(x_j);
		

		for (int m = 0; m < 3; m++)
		{
			r1(m + 1) = cxsat(m + 1) - x_j(m + 1);
		}


		r1length = sqrt(r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3));
		_slr_rho = r1length;

		//calculate elevation and azimuth for lra
		xtmp_leo1 = r1.t() * rot_scf2crsleo;
		for (int m = 0; m < 3; m++)
		{
			xtmp_leo(m + 1) = xtmp_leo1(1, m + 1);
		}
		for (int m = 0; m < 3; m++)
		{
			xtmp_leo(m + 1) = -xtmp_leo(m + 1);
		}
		string _sat;
		_sat = strs2strl(_sat_leo);
		if (_sat.find("SENTINEL") != string::npos)
		{
			xtmp_leo(1) -= 1.13403;
			xtmp_leo(2) -= 0.637905;
			xtmp_leo(3) -= 0.80118;
			xtmp_lrr(1) = xtmp_leo(2);
			xtmp_lrr(2) = -xtmp_leo(1);
			xtmp_lrr(3) = xtmp_leo(3);
			azim_leo = atan2(xtmp_lrr(2), xtmp_lrr(1));
			elev_leo = atan2(xtmp_lrr(3), sqrt(xtmp_lrr(1)*xtmp_lrr(1) + xtmp_lrr(2)*xtmp_lrr(2)));
		}
		else
		{
			azim_leo = atan2(xtmp_leo(2), xtmp_leo(1));
			elev_leo = atan2(sqrt(xtmp_leo(1)*xtmp_leo(1) + xtmp_leo(2)*xtmp_leo(2)), xtmp_leo(3));
		}
		get_lracor(_sat,elev_leo, azim_leo, &dlra_f);


		for (int m = 0; m < 3; m++)
		{
			double ttmp = _allslrf->getOnegslrf(ista)._xsta[m] + _allslrf->getOnegslrf(ista)._offset[m] * 0.001;
			x_f(m + 1) = _allslrf->getOnegslrf(ista)._xsta[m] + _allslrf->getOnegslrf(ista)._offset[m] * 0.001;
			x_f(m + 1) = x_f(m + 1)*1000.0;
			v_f(m + 1) = 0.0;
		}

		btime.from_mjd(jd, (int)sodb, sodb - (int)sodb);

		_calculate_sta_pos(btime, &x_f, &v_f, &x_j, &v_j);

		for (int m = 0; m < 3; m++)
		{
			r1(m + 1) = cxsat(m + 1) - x_j(m + 1);
		}

		for (int m = 0; m < 3; m++)
		{
			dloudx(m + 1) = r1(m + 1) / sqrt(r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3));
		}

		r1length = r1length + sqrt(r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3));

		xtmp_leo1 = r1.t() * rot_scf2crsleo;

		for (int m = 0; m < 3; m++)
		{
			xtmp_leo(m + 1) = xtmp_leo1(1, m + 1);
		}

		for (int m = 0; m < 3; m++)
		{
			xtmp_leo(m + 1) = -xtmp_leo(m + 1);
		}
		
		_sat = strs2strl(_sat_leo);
		if (_sat.find("SENTINEL") != string::npos)
		{
			xtmp_leo(1) -= 1.13403;
			xtmp_leo(2) -= 0.637905;
			xtmp_leo(3) -= 0.80118;
			xtmp_lrr(1) = xtmp_leo(2);
			xtmp_lrr(2) = -xtmp_leo(1);
			xtmp_lrr(3) = xtmp_leo(3);
			azim_leo = atan2(xtmp_lrr(2), xtmp_lrr(1));
			elev_leo = atan2(xtmp_lrr(3), sqrt(xtmp_lrr(1)*xtmp_lrr(1) + xtmp_lrr(2)*xtmp_lrr(2)));
		}
		else
		{
			azim_leo = atan2(xtmp_leo(2), xtmp_leo(1));
			elev_leo = atan2(sqrt(xtmp_leo(1)*xtmp_leo(1) + xtmp_leo(2)*xtmp_leo(2)), xtmp_leo(3));
		}
	    get_lracor(_sat,elev_leo, azim_leo, &dlra);


		Matrix xtmp1;
		xtmp1.ReSize(1, 3);

		xtmp1 = r1.t()*rot_trs2crs;

		for (int m = 0; m < 3; m++)
		{
			xtmp(m + 1) = xtmp1(1, m + 1);
		}

		xtmp1 = xtmp.t()*rot_l2f;

		for (int m = 0; m < 3; m++)
		{
			xtmp(m + 1) = xtmp1(1, m + 1);
		}

		azim = atan2(xtmp(2), xtmp(1));
		elev = atan2(xtmp(3), sqrt(xtmp(1)*xtmp(1) + xtmp(2)*xtmp(2)));

		//_gnpt->getAllNpt(_sat, i).setElev(elev*R2D);

		troposphere_delay_slr(elev, &trpdel, &ztdpart, now, ista);


		temp = preproc_slr._dry_tem - 273.150;
		press_wvp = preproc_slr._atm_humidity*exp(2.30260*(0.78580 + 7.50*temp / (237.30 + temp)));

		trop_slr_fculzd_hpa(now, ista, elev, &fcul_ztd, &fcul_zhd, &fcul_zwd, press_wvp);

		fak_map = trop_slr_fcul_a(elev, now, ista);

		std = fak_map * fcul_ztd;
		ztd_dif = (trpdel - std * 1e-3)*1e3;
		trpdel = std * 1.e-3;

		reldel0 = -2.0*(x_j(1)*v_j(1) + x_j(2)*v_j(2) + x_j(3)*v_j(3)) / (CLIGHT*1e-3);
		reldel0 = reldel0 + 2.0*(cxsat(1)*1e-3*cvsat(1)*1e-3 + cxsat(2)*1e-3*cvsat(2)*1e-3 + cxsat(3)*1e-3*cvsat(3)*1e-3) / (CLIGHT*1e-3);

		sitrad = sqrt(x_j(1)*x_j(1) + x_j(2)*x_j(2) + x_j(3)*x_j(3));
		satrad = sqrt(cxsat(1)*1e-3*cxsat(1)*1e-3 + cxsat(2)*1e-3*cxsat(2)*1e-3 + cxsat(3)*1e-3*cxsat(3)*1e-3);

		reldel1 = 2.0*EARTH_GM / ((CLIGHT*1e-3)*(CLIGHT*1e-3));
		reldel1 = log((sitrad + satrad + (r1length * 1e-3) / 2.0) / (sitrad + satrad - (r1length * 1e-3) / 2.0));
		reldel1 = 2.0*EARTH_GM / ((CLIGHT*1e-3)*(CLIGHT*1e-3)) * log((sitrad + satrad + (r1length * 1e-3) / 2.0) / (sitrad + satrad - (r1length * 1e-3) / 2.0));

		//double omc = (preproc_slr._two_way_time*(CLIGHT*1e-3) - (r1length*1e-3 + 2.0*trpdel + 2.0*reldel1))*1e3 / 2.0;
		//omc = ((dlra / 1e3 + dlra_f / 1e3) / 2.0);
		omc = (preproc_slr._two_way_time*(CLIGHT*1e-3) - (r1length*1e-3 + 2.0*trpdel + 2.0*reldel1))*1e3 / 2.0 -range_bias+ ((dlra / 1e3 + dlra_f / 1e3) / 2.0);
		ostringstream out;
		out << omc;
		//write_log_info(_log, 1, "OMC", out.str());
		writeLogInfo(_log, 0, "NOTE", out.str());

		return true;

	}

	void t_gslrmodel::enu2trs(double lat, double lon, Matrix& rotmat)
	{
		double coslat, sinlat, coslon, sinlon;

		coslat = cos(lat - G_PI / 2.0);
		sinlat = sin(lat - G_PI / 2.0);

		coslon = cos(-G_PI / 2.0 - lon);
		sinlon = sin(-G_PI / 2.0 - lon);

		rotmat(1, 1) = coslon;
		rotmat(1, 2) = sinlon * coslat;
		rotmat(1, 3) = sinlon * sinlat;
		rotmat(2, 1) = -sinlon;
		rotmat(2, 2) = coslon * coslat;
		rotmat(2, 3) = coslon * sinlat;
		rotmat(3, 1) = 0.0;
		rotmat(3, 2) = -sinlat;
		rotmat(3, 3) = coslat;

	}

	double t_gslrmodel::trop_slr_fcul_a(double elev, t_gtime now, int ista)
	{
		double EPSILON, SINE, T_C, COSPHI;
		double A10, A11, A12, A13;
		double A20, A21, A22, A23;
		double A30, A31, A32, A33;
		double A1, A2, A3;
		double MAP_ZEN;
		double trop_slr_FCUL_A;

		t_gpreprocslr  preproc_slr;
		preproc_slr = _allpreproc->getSLROBS(strs2strl(_sat_leo),now);

		EPSILON = elev;
		SINE = sin(EPSILON);
		T_C = preproc_slr._dry_tem - 273.1500;

		COSPHI = cos(_allslrf->getOnegslrf(ista)._blh[0]);

		A10 = 0.121008e-02;
		A11 = 0.17295e-05;
		A12 = 0.3191e-04;
		A13 = -0.18478e-07;

		A20 = 0.304965e-02;
		A21 = 0.2346e-05;
		A22 = -0.1035e-03;
		A23 = -0.1856e-07;

		A30 = 0.68777e-01;
		A31 = 0.1972e-04;
		A32 = -0.3458e-02;
		A33 = 0.1060e-06;

		A1 = A10 + A11 * T_C + A12 * COSPHI + A13 * _allslrf->getOnegslrf(ista)._blh[2] * 1000.0;
		A2 = A20 + A21 * T_C + A22 * COSPHI + A23 * _allslrf->getOnegslrf(ista)._blh[0] * 1000.0;
		A3 = A30 + A31 * T_C + A32 * COSPHI + A33 * _allslrf->getOnegslrf(ista)._blh[0] * 1000.0;

		MAP_ZEN = (1.000 + A1 / (1.000 + A2 / (1.000 + A3)));

		trop_slr_FCUL_A = MAP_ZEN / (SINE + A1 / (SINE + A2 / (SINE + A3)));

		return trop_slr_FCUL_A;

	}

	void t_gslrmodel::trop_slr_fculzd_hpa(t_gtime now, int ista, double elev, double* FCUL_ZTD, double* FCUL_ZHD, double* FCUL_ZWD, double WVP)
	{
		double k0, k1, k2, k3, XC, CORR;
		double SIGMA, W0, W1, W2, W3, FH, FNH;
		double F;

		t_gpreprocslr  preproc_slr;
		preproc_slr = _allpreproc->getSLROBS(strs2strl(_sat_leo),now);

		XC = 375.00;
		k0 = 238.01850;
		k1 = 19990.9750;
		k2 = 57.3620;
		k3 = 579.551740;

		W0 = 295.2350;
		W1 = 2.64220;
		W2 = -0.0323800;
		W3 = 0.0040280;

		SIGMA = 1000.0 / (preproc_slr._wave*1000.0);

		F = 1 - 0.00266*cos(2.0 * _allslrf->getOnegslrf(ista)._blh[0]) - 0.00028e-3 * _allslrf->getOnegslrf(ista)._blh[2] * 1000.0;

		CORR = 1.00 + 0.534e-6 * (XC - 450);

		FH = 0.010*CORR*((k1*(k0 + SIGMA * SIGMA)) / ((k0 - SIGMA * SIGMA)*(k0 - SIGMA * SIGMA)) + k3 * (k2 + SIGMA * SIGMA) / ((k2 - SIGMA * SIGMA)*(k2 - SIGMA * SIGMA)));

		*FCUL_ZHD = 2.416579e-3 * FH*preproc_slr._atm_pressure / F;

		FNH = 0.0031010*(W0 + 3.0*W1*SIGMA*SIGMA + 5.0*W2*SIGMA*SIGMA*SIGMA*SIGMA + 7.0*W3*SIGMA*SIGMA*SIGMA*SIGMA*SIGMA*SIGMA);

		*FCUL_ZWD = 1e-4 * (5.316e0*FNH - 3.759*FH)*WVP / F;

		*FCUL_ZTD = *FCUL_ZHD + *FCUL_ZWD;
	}

	void t_gslrmodel::troposphere_delay_slr(double elev, double* delay, double* partial,t_gtime now,int ista)
	{
		double T0, fw, es, fpos, fwvl, e0, a, k, b;

		t_gpreprocslr  preproc_slr;
	    preproc_slr = _allpreproc->getSLROBS(strs2strl(_sat_leo),now);
		T0 = preproc_slr._dry_tem;
		fw = 1.000620 + 3.14e-6 * preproc_slr._atm_pressure + 5.6e-7 * (T0 - 273.15)*(T0 - 273.15);
		es = 0.01*exp(1.2378847e-5 * T0*T0 - 1.9121316e-2 * T0 + 33.937110470 - 6.3431645e3 / T0);
		fpos = 1.0 - 0.0026*cos(2.0*_allslrf->getOnegslrf(ista)._blh[0]) - 0.00031*_allslrf->getOnegslrf(ista)._blh[2];
		fwvl = 0.9650 + 0.0164 / preproc_slr._wave / preproc_slr._wave + 0.000228*pow((preproc_slr._wave), -4);
		e0 = preproc_slr._atm_humidity * es*fw * 1e-2;
		a = 0.002357* preproc_slr._atm_pressure + 0.000141*e0;
		k = 1.163 - 0.00968*cos(2 * _allslrf->getOnegslrf(ista)._blh[0]) - 0.00104*T0 + 0.00001435*preproc_slr._atm_pressure;

		b = 1.084e-8 * preproc_slr._atm_pressure*T0*k + 4.734e-8 * preproc_slr._atm_pressure*preproc_slr._atm_pressure / T0 * 2.0 / (3.0 - 1.0 / k);
		*partial = 1.0 / (sin(elev) + b / (a + b) / (sin(elev) + 0.010));
		*delay = fwvl / fpos * (a + b)*(*partial) * 1e-3;

	}

	void t_gslrmodel::get_lracor(string sat,double zeni0, double azim0, double* dlra)
	{
		if (!_alllra)
		{
			*dlra = 0.0;
			return;
		}
		double zeni;
		double azim;
		int ntmp;
		int iazi, izen;
		double x1, x2;
		double alpha;

		zeni = zeni0 * R2D;

		if (zeni > _alllra->zen2(sat))
			zeni = _alllra->zen2(sat) - 1e-5;
		if (zeni < _alllra->zen1(sat))
			zeni = _alllra->zen1(sat);

		if (_alllra->dzen(sat) == 0.0)
		{
			write_log_info(_log, 1, "slr", "ERROR : slr ==> dzen in lra is zero");
			return;
		}

		azim = azim0 * R2D;
		if (azim < 0.0)
			azim = azim + 360.0;
		if (sat.find("SENTINEL") != string::npos)
		{
			if (azim > 60.0)
			{
				ntmp = int(azim / 60.0);
				azim = azim - 60.0*ntmp;
			}
		}
		else
		{
			if (azim > 90.0)
			{
				ntmp = int(azim / 90.0);
				azim = azim - 90.0*ntmp;
			}
			if (azim > 45.0)
				azim = 90.0 - azim;
		}

		iazi = 0;

		if (_alllra->dazi(sat) != 0)
			iazi = int(azim / _alllra->dazi(sat) - 1e-5);

		if (_alllra->dzen(sat) != 0)
			izen = int((zeni - _alllra->zen1(sat)) / _alllra->dzen(sat));

		x1 = _alllra->getlra(sat,izen, iazi);
		x2 = _alllra->getlra(sat,izen + 1, iazi);

		if (iazi != 0)
		{
			alpha = (azim - _alllra->azi1(sat)) / _alllra->dazi(sat) - iazi;
			x1 = x1 + (_alllra->getlra(sat,izen, iazi + 1) - _alllra->getlra(sat,izen, iazi))*alpha;
			x2 = x2 + (_alllra->getlra(sat,izen + 1, iazi + 1) - _alllra->getlra(sat,izen + 1, iazi))*alpha;
		}

		alpha = (zeni - _alllra->zen1(sat)) / _alllra->dzen(sat) - izen;
		*dlra = *dlra + x1 + (x2 - x1)*alpha;

		if (sat.find("SENTINEL") != string::npos)
		{
			*dlra = *dlra*(-1);
		}

	}



	void t_gslrmodel::_machine_offset(ColumnVector& xyz)
	{
		int islr = 0;
		double              _slre;            ///< the term S : the vector from the center of mass to the reference point
		double              _slrn;            ///< slre and slrn are always negative
		double              _slru;

		string _sat;
		_sat = strs2strl(_sat_leo);

		if (_sat.compare("SWARM-A") == 0)
		{
			islr = 1;
			_slre = -0.4633;
			_slrn = 0.51979;
			_slru = 0.30281;
		}
		else if (_sat.compare("SWARM-B") == 0)
		{
			islr = 2;
			_slre = -0.46356;
			_slrn = 0.52097;
			_slru = 0.30239;
		}
		else if (_sat.compare("SWARM-C") == 0)
		{
			islr = 3;
			_slre = -0.46592;
			_slrn = 0.52009;
			_slru = 0.30289;
		}
		else if (_sat.compare("GRACE-A") == 0)
		{
			islr = 4;
			_slre = -0.6000;
			_slrn = -0.3275;
			_slru = 0.2178;
		}
		else if (_sat.compare("GRACE-B") == 0)
		{
			islr = 4;
			_slre = -0.6000;
			_slrn = -0.3275;
			_slru = 0.2178;
		}
		else if (_sat.compare("GRACE-C") == 0)
		{
			islr = 6;
			_slre = -0.6000;
			_slrn = -0.3275;
			_slru = 0.2208;
		}
		else if (_sat.compare("GRACE-D") == 0)
		{
			islr = 7;
			_slre = -0.6000;
			_slrn = -0.3275;
			_slru = 0.2208;
		}
		else if (_sat.compare("SENTINEL-3A") == 0)
		{
			islr = 8;
			_slre = -0.3630;
			_slrn = 0.4279;
			_slru = 0.7922;
		}
		else if (_sat.compare("SENTINEL-3B") == 0)
		{
			islr = 8;
			_slre = -0.3530;
			_slrn = 0.4359;
			_slru = 0.7892;
		}
		else if (_sat.compare("HAIYANG-2A") == 0)
		{
			islr = 9;
			_slre = -1.01840;
			_slrn = -0.21930;
			_slru = 0.98730;
		}
		else if (_sat.compare("TIANGONG-2") == 0)
		{
			islr = 10;
			_slre = 0.0195;
			_slrn = -0.0044;
			_slru = 1.6112;

		}
		else if (_sat.compare("JASON-2") == 0)
		{
			islr = 11;
			_slre = 0.2172;
			_slrn = 0.5981;
			_slru = 0.6829;

		}
		else if (_sat.compare("JASON-3") == 0)
		{
			islr = 12;
			_slre = 0.2172;
			_slrn = 0.5981;
			_slru = 0.6829;

		}
		xyz(1) = _slre;
		xyz(2) = _slrn;
		xyz(3) = _slru;
	}

	void t_gslrmodel::_calculate_sta_pos(t_gtime time, ColumnVector* x_f, ColumnVector* v_f, ColumnVector* x_j, ColumnVector* v_j)
	{
		Matrix rot_scf2crsleo(0.0, 3, 3); //rotation matrix of SCF to CRS, just for leo!
		Matrix rot_trs2crs(0.0, 3, 3);   // rotation matrix of TRS to CRS
		Matrix rot_du(0.0, 3, 3);
		// TRS2CRS 
		t_gtime tdt_gtime = time;
		//ChangeTimeSys(GPS2TDT, time, tdt_gtime);
		tdt_gtime.tsys(t_gtime::TT);

		//t_gtime tdt_grtime = rtime;
		//tdt_grtime.tsys(t_gtime::TT);
		//double tdt = tdt_gtime.dmjd(false);

		_trs2crs_2000->calcRotMat(tdt_gtime, true, true, true);
		rot_trs2crs = _trs2crs_2000->getRotMat();
		rot_du = _trs2crs_2000->getMatDu();
		rot_du = rot_du * OMGE_DOT_GAL;

		_rotate_posvel(*x_f, *v_f, *x_j, v_j, rot_trs2crs, rot_du);

		t_gtriple tx_f(*x_f);

		bool tide_valid = _apply_rec_tides(time, tx_f);
		if (!tide_valid)
		{
			write_log_info(_log, 1, "NOTE", "apply tide failed");
		}

		//tx_f = tx_f + _crt_obj->eccxyz(time);
		for (int m = 0; m < 3; m++)
		{
			(*x_f)(m + 1) = tx_f.crd(m);
		}

		(*x_j) = rot_trs2crs * (*x_f);
	}

	void t_gslrmodel::_rotate_posvel(ColumnVector& x_f, ColumnVector& v_f, ColumnVector& x_j, ColumnVector* v_j, Matrix& rotmat, Matrix& drotmat)
	{
		ColumnVector xhelp(0.0, 3);

		x_j = rotmat * x_f;
		*v_j = rotmat * v_f;
		xhelp = drotmat * x_f;
		*v_j = *v_j + xhelp;
	}

	bool t_gslrmodel::_apply_rec_tides(const t_gtime& epoch, t_gtriple& rec)
	{
		if (!_tide)
		{
			write_log_info(_log, 1, "ERROR", "The tide ptr is nullptr.");
			return false;
		}

		_update_rot_matrix(epoch);
		double xpole = _trs2crs_2000->getXpole();
		double ypole = _trs2crs_2000->getYpole();
		double gast = _trs2crs_2000->getGmst();
		Matrix rot_trs2crs = _trs2crs_2000->getRotMat();
		string _site = "";

		t_gtriple tide(0.0, 0.0, 0.0);
		try
		{
			// solid tide
			t_gtideIERS* solid_ptr = dynamic_cast<t_gtideIERS*>(_tide);
			if (!solid_ptr)
			{
				write_log_info(_log, 1, "ERROR", "can not get the solid tide ptr.");
				return false;
			}
			t_gtriple solid_earth = solid_ptr->tide_solid(epoch, rec, rot_trs2crs, _allde);

			// ocean load
			t_gtriple load_ocean = _tide->load_ocean(epoch, _site, rec);

			//pole tide
			t_gtideIERS* pole_ptr = dynamic_cast<t_gtideIERS*>(_tide);
			if (!pole_ptr)
			{
				write_log_info(_log, 1, "ERROR", "can not get the pole tide ptr.");
				return false;
			}
			t_gtriple tide_pole = pole_ptr->tide_pole_pod(epoch, xpole, ypole, rec);

			//atmosph
			t_gtriple load_atmosph = _tide->load_atmosph();

			// freq tide
			t_gtideIERS* freq_ptr = dynamic_cast<t_gtideIERS*>(_tide);
			if (!freq_ptr)
			{
				write_log_info(_log, 1, "ERROR", "can not get the freq tide ptr.");
				return false;
			}
			t_gtriple tide_freq = freq_ptr->tide_freq(_site, rec, gast);

			// jdhuang  for test
			tide = solid_earth +
				load_ocean +
				tide_pole +
				load_atmosph +
				tide_freq;

#ifdef DEBUG_PRECISEMODEL
			cout << "Tide" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd(false)
				<< " " << epoch.sod(false)
				<< " " << tide[0]
				<< " " << tide[1]
				<< " " << tide[2] << endl;

			cout << "solid_earth" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd(false)
				<< " " << epoch.sod(false)
				<< " " << solid_earth[0]
				<< " " << solid_earth[1]
				<< " " << solid_earth[2] << endl;

			cout << "load_ocean" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd(false)
				<< " " << epoch.sod(false)
				<< " " << load_ocean[0]
				<< " " << load_ocean[1]
				<< " " << load_ocean[2] << endl;

			cout << "tide_pole" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd(false)
				<< " " << epoch.sod(false)
				<< " " << tide_pole[0]
				<< " " << tide_pole[1]
				<< " " << tide_pole[2] << endl;

			cout << "load_atmosph" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd(false)
				<< " " << epoch.sod(false)
				<< " " << load_atmosph[0]
				<< " " << load_atmosph[1]
				<< " " << load_atmosph[2] << endl;

			cout << "tide_freq" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd(false)
				<< " " << epoch.sod(false)
				<< " " << tide_freq[0]
				<< " " << tide_freq[1]
				<< " " << tide_freq[2] << endl;
#endif // DEBUG

		}
		catch (...)
		{
			write_log_info(_log, 1, "ERROR", "can not get the tide ptr.");
			return false;
		}

		//unit to m
		rec = rec + tide * 1.e3;
		return true;

	}

	void t_gslrmodel::_update_rot_matrix(const t_gtime& epoch)
	{
		if (!_trs2crs_2000 && !_allerp)
		{
			write_log_info(_log, 1, "t_gorbmodel::_update_rotmatrix ", "have no poleut1 data,cant' compute trs2crs!!");
			throw exception();
		}

		if (!_trs2crs_2000)
		{
			_trs2crs_2000 = new t_gtrs2crs(false, _allerp);
		}

		t_gtime tdt = epoch;
		tdt.tsys(t_gtime::TT);

		if (tdt != _trs2crs_2000->getCurtEpoch())
		{
			_trs2crs_2000->calcRotMat(tdt, true, true, true);
		}

#ifdef DEBUG_TRS
		//jdhuang : check for the trs2crs_2000
		cout << "Mark = trs2crs  : "
			<< setw(25) << tdt.mjd(false)
			<< setw(25) << tdt.sod(false)
			<< setw(25) << tdt.dsec(false) << endl
			<< setw(25) << fixed << setprecision(15)
			<< trs2crs_2000->rotmat() << endl
			<< trs2crs_2000->rotdx() << endl
			<< trs2crs_2000->rotdy() << endl
			<< trs2crs_2000->rotdu() << endl;
#endif // DEBUG
		return;

	}

}
