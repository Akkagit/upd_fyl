// history :
//  jdhuang 2019-05-12 add the case of RAW_DOUBLE in Weight_Obs


#include "gmodels/gprecisemodel.h"

#include "gutils/ginfolog.h"
#include "gutils/gleo.h"
#include "gutils/gstring.h"
#include <algorithm>
#include <assert.h>

// #define DEBUG_LOG
// #define DEBUG_GLO
namespace great
{
	t_gprecisemodel::t_gprecisemodel(t_gallproc* data, t_glog* log, t_gsetbase* setting) :
		t_gpppmodel("", log, setting),
		_minElev(7),
		_minElevLeo(1),
		_gdata_erp(0),
		_scf2crs_leo(t_gscf2crsleo("")),
		_rec_obj_flag(make_pair("", t_gtime()))
	{
		// Get Settings
		_tide        = shared_ptr<t_gtide>(new t_gtideIERS(log));
		_minElev     = dynamic_cast<t_gsetproc*>(setting)->minimum_elev();
		_minElevLeo  = dynamic_cast<t_gsetproc*>(setting)->minimum_elev_leo(); // for LEO by zhangwei
		_weight      = dynamic_cast<t_gsetproc*>(setting)->weighting();

		_mf_ztd      = dynamic_cast<t_gsetproc*>(setting)->tropo_mf();
		_mf_grd      = dynamic_cast<t_gsetproc*>(setting)->grad_mf();

		_frequency  = dynamic_cast<t_gsetproc*>(setting)->frequency();
		_ion_model  = dynamic_cast<t_gsetproc*>(setting)->ion_model();
		_realtime   = dynamic_cast<t_gsetproc*>(setting)->realtime(); //add by xiongyun
		_ultrasp3   = dynamic_cast<t_gsetproc*>(setting)->ultrasp3(); //add by xiongyun

		_crd_est = dynamic_cast<t_gsetproc*>(setting)->crd_est();

		_isTriple    = false;

		log->logInfo("t_gprecisemodel", "min Elev",     format("%16.8f", _minElev));
		log->logInfo("t_gprecisemodel", "min Elev leo", format("%16.8f", _minElevLeo));
		log->logInfo("t_gprecisemodel", "weight",       format("%16.8f", _weight));
		log->logInfo("t_gprecisemodel", "mf ztd", dynamic_cast<t_gsetproc*>(setting)->ztdmpfunc2str(_mf_ztd));
		log->logInfo("t_gprecisemodel", "mf grd", dynamic_cast<t_gsetproc*>(setting)->grdmpfunc2str(_mf_grd));

		_sigCodeGPS = dynamic_cast<t_gsetgnss*>(setting)->sigma_C(GPS);
		_sigCodeGLO = dynamic_cast<t_gsetgnss*>(setting)->sigma_C(GLO);
		_sigCodeGAL = dynamic_cast<t_gsetgnss*>(setting)->sigma_C(GAL);
		_sigCodeBDS = dynamic_cast<t_gsetgnss*>(setting)->sigma_C(BDS);
		_sigCodeQZS = dynamic_cast<t_gsetgnss*>(setting)->sigma_C(QZS);
		_sigPhaseGPS = dynamic_cast<t_gsetgnss*>(setting)->sigma_L(GPS);
		_sigPhaseGLO = dynamic_cast<t_gsetgnss*>(setting)->sigma_L(GLO);
		_sigPhaseGAL = dynamic_cast<t_gsetgnss*>(setting)->sigma_L(GAL);
		_sigPhaseBDS = dynamic_cast<t_gsetgnss*>(setting)->sigma_L(BDS);
		_sigPhaseQZS = dynamic_cast<t_gsetgnss*>(setting)->sigma_L(QZS);
		
		log->logInfo("t_gprecisemodel", "sigCodeGPS ", format("%16.8f", _sigCodeGPS));
		log->logInfo("t_gprecisemodel", "sigCodeGLO ", format("%16.8f", _sigCodeGLO));
		log->logInfo("t_gprecisemodel", "sigCodeGAL ", format("%16.8f", _sigCodeGAL));
		log->logInfo("t_gprecisemodel", "sigCodeBDS ", format("%16.8f", _sigCodeBDS));
		log->logInfo("t_gprecisemodel", "sigCodeQZS ", format("%16.8f", _sigCodeQZS));
		log->logInfo("t_gprecisemodel", "sigPhaseGPS", format("%16.8f", _sigPhaseGPS));
		log->logInfo("t_gprecisemodel", "sigPhaseGLO", format("%16.8f", _sigPhaseGLO));
		log->logInfo("t_gprecisemodel", "sigPhaseGAL", format("%16.8f", _sigPhaseGAL));
		log->logInfo("t_gprecisemodel", "sigPhaseBDS", format("%16.8f", _sigPhaseBDS));
		log->logInfo("t_gprecisemodel", "sigPhaseQZS", format("%16.8f", _sigPhaseQZS));

		//for LEO processing by zhangwei
		_sigCodeGPSLEO = dynamic_cast<t_gsetgnss*>(setting)->sigma_C_LEO(GPS);
		_sigCodeGLOLEO = dynamic_cast<t_gsetgnss*>(setting)->sigma_C_LEO(GLO);
		_sigCodeGALLEO = dynamic_cast<t_gsetgnss*>(setting)->sigma_C_LEO(GAL);
		_sigCodeBDSLEO = dynamic_cast<t_gsetgnss*>(setting)->sigma_C_LEO(BDS);
		_sigCodeQZSLEO = dynamic_cast<t_gsetgnss*>(setting)->sigma_C_LEO(QZS);
		_sigPhaseGPSLEO = dynamic_cast<t_gsetgnss*>(setting)->sigma_L_LEO(GPS);
		_sigPhaseGLOLEO = dynamic_cast<t_gsetgnss*>(setting)->sigma_L_LEO(GLO);
		_sigPhaseGALLEO = dynamic_cast<t_gsetgnss*>(setting)->sigma_L_LEO(GAL);
		_sigPhaseBDSLEO = dynamic_cast<t_gsetgnss*>(setting)->sigma_L_LEO(BDS);
		_sigPhaseQZSLEO = dynamic_cast<t_gsetgnss*>(setting)->sigma_L_LEO(QZS);

		log->logInfo("t_gprecisemodel", "sigCodeGPSLEO  ", format("%16.8f", _sigCodeGPSLEO));
		log->logInfo("t_gprecisemodel", "sigCodeGLOLEO  ", format("%16.8f", _sigCodeGLOLEO));
		log->logInfo("t_gprecisemodel", "sigCodeGALLEO  ", format("%16.8f", _sigCodeGALLEO));
		log->logInfo("t_gprecisemodel", "sigCodeBDSLEO  ", format("%16.8f", _sigCodeBDSLEO));
		log->logInfo("t_gprecisemodel", "sigCodeQZSLEO  ", format("%16.8f", _sigCodeQZSLEO));
		log->logInfo("t_gprecisemodel", "sigPhaseGPSLEO ", format("%16.8f", _sigPhaseGPSLEO));
		log->logInfo("t_gprecisemodel", "sigPhaseGLOLEO ", format("%16.8f", _sigPhaseGLOLEO));
		log->logInfo("t_gprecisemodel", "sigPhaseGALLEO ", format("%16.8f", _sigPhaseGALLEO));
		log->logInfo("t_gprecisemodel", "sigPhaseBDSLEO ", format("%16.8f", _sigPhaseBDSLEO));
		log->logInfo("t_gprecisemodel", "sigPhaseQZSLEO ", format("%16.8f", _sigPhaseQZSLEO));


		_band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::GPS);
		_band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::GAL);
		_band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::GLO);
		_band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::BDS);
		_band_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::QZS);

		_freq_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::GPS);
		_freq_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::GAL);
		_freq_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::GLO);
		_freq_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::BDS);
		_freq_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::QZS);

		// Get Data
		setOTL(dynamic_cast<t_gallotl*>((*data)[t_gdata::ALLOTL]));
		setOBJ(dynamic_cast<t_gallobj*>((*data)[t_gdata::ALLOBJ]));
		_gall_nav    = dynamic_cast<t_gallnav*>((*data)[t_gdata::GRP_EPHEM]);
		_gall_orb    = dynamic_cast<t_gallorb*>((*data)[t_gdata::ALLORB]);
		_gdata_erp   = dynamic_cast<t_gpoleut1*>((*data)[t_gdata::ALLPOLEUT1]);
		_gdata_navde = dynamic_cast<t_gnavde*>((*data)[t_gdata::ALLDE]);
		_gall_att    = dynamic_cast<t_gallatt*>((*data)[t_gdata::AllATTITUDE]);
		_gdata_ion   = dynamic_cast<t_gion*>((*data)[t_gdata::ION]);                   //< add by jdhuang,used for getting UCUD amb value

		// add by glfeng
		_gallbias = dynamic_cast<t_gallbias*>((*data)[t_gdata::ALLBIAS]);
		_gifcb    = dynamic_cast<t_gifcb*>((*data)[t_gdata::IFCB]);
	}

	t_gprecisemodel::~t_gprecisemodel()
	{


#ifdef DEBUG_LOG
		if (_debug)
		{
			delete _debug;
			_debug = nullptr;
		}

		if (_B_log)
		{
			delete _B_log;
			_B_log = nullptr;
		}

		if (_P_log)
		{
			delete _P_log;
			_P_log = nullptr;
		}

		if (_L_log)
		{
			delete _L_log;
			_L_log = nullptr;
		}

		if (_dx_log)
		{
			delete _dx_log;
			_dx_log = nullptr;
		}
#endif // DEBUG_LOG
	}

	void t_gprecisemodel::setOTL(t_gallotl* otl)
	{
		_tide->setOTL(otl);
	}
	double t_gprecisemodel::windUp(t_gsatdata& satdata, const ColumnVector& rRec)
	{
		gtrace("t_gprecisemodel::windUp");

		t_gtime epoch = satdata.epoch();
		string  prn   = satdata.sat();
		ColumnVector rSat = _trs_sat_crd.crd_cvect();

		// glfeng not change time
		double Mjd = epoch.mjd() + epoch.sod() / 86400.0;

		// First time - initialize to zero
		// -------------------------------
		if (_phase_windup.find(_site) == _phase_windup.end() ||
			_phase_windup[_site].find(prn) == _phase_windup[_site].end() ||
			_phase_windup[_site][prn].size() == 0)
		{
			_phase_windup[_site][prn][epoch] = 0.0;
		}

		// Compute the correction for new time
		// -----------------------------------
		if (_phase_windup[_site][prn].find(epoch) == _phase_windup[_site][prn].end() ||
			_phase_windup[_site][prn].size() == 1) {

			// the last epoch
			double dphi0 = _phase_windup[_site][prn].rbegin()->second;
			ColumnVector rho = rRec - rSat;  //glfeng
			rho /= rho.norm_Frobenius();

			ColumnVector i(3);
			ColumnVector j(3);
			ColumnVector k(3);

			// attitude model
			string antype = "";
			if (_gallobj != 0) {
				shared_ptr<t_gobj>  sat_obj = _gallobj->obj(satdata.sat());
				shared_ptr<t_gpcv>  sat_pcv;
				if (sat_obj != 0)  sat_pcv = sat_obj->pcv(satdata.epoch());
				if (sat_pcv != 0)  antype = sat_pcv->anten();
			}
			if (_attitudes == YAW_NOMI)		   attitude(satdata, "", i, j, k);				//nominal modeling
			else if (_attitudes == YAW_RTCM) attitude(satdata, satdata.yaw(), i, j, k);	//value from RTCM used
			else							   attitude(satdata, antype, i, j, k);			//default
																							//
			if (antype.find("BLOCK IIR") != string::npos)
			{
				i = (-1.0)*i;
				j = (-1.0)*j;
			}
			double rlength = DotProduct(rho, i);
			ColumnVector dipSat = i - rlength * rho - crossproduct(rho, j);

			// Receiver unit Vectors rx, ry
			// ----------------------------
			ColumnVector rx(3);
			ColumnVector ry(3);

			double recEll[3]; xyz2ell(rRec.data(), recEll, false);
			double neu[3];

			if (_grec->id_type() == t_gdata::REC_LEO)
			{
				update_rotmat_scf(_site, _rec_epo, _crs_rec_crd.crd_cvect(), _crs_rec_vel.crd_cvect());
				// for LEO
				for (int irow = 0; irow < 3; irow++)
				{
					rx[irow] = _rot_scf2trs_leo[irow][0];
					ry[irow] = -_rot_scf2trs_leo[irow][1];
				}
			}
			else
			{
				neu[0] = 1.0;
				neu[1] = 0.0;
				neu[2] = 0.0;
				neu2xyz(recEll, neu, rx.data());

				neu[0] = 0.0;
				neu[1] = -1.0;
				neu[2] = 0.0;
				neu2xyz(recEll, neu, ry.data());
			}

			// Effective Dipole of the Receiver Antenna
			// ----------------------------------------
			rlength = DotProduct(rho, rx);
			ColumnVector dipRec = rx - rlength * rho + crossproduct(rho, ry);

			// Resulting Effect
			// ----------------
			double alpha = DotProduct(dipSat, dipRec) / (dipSat.norm_Frobenius() * dipRec.norm_Frobenius());

			if (alpha > 1.0) alpha = 1.0;
			if (alpha < -1.0) alpha = -1.0;

			double dphi = acos(alpha) / 2.0 / G_PI;  // in cycles

			if (DotProduct(rho, crossproduct(dipSat, dipRec)) < 0.0) dphi = -dphi;

			_phase_windup[_site][prn][epoch] = floor(dphi0 - dphi + 0.5) + dphi;
		}

#ifdef DEBUG_CORR
		cout << "Wind up " << epoch.str_ymdhms() << " " << prn << " " << _windUpSum[prn] << " " << satdata.orb_angle()*R2D << endl;
		//    int ooo; cin >> ooo;
#endif   

		satdata.addwind(_phase_windup[_site][prn][epoch]);
		return _phase_windup[_site][prn][epoch];
	}

	bool t_gprecisemodel::_prepare_obs(const t_gtime& crt_epo, t_gallpar& pars)
	{
		_factorP = 1;
		if (!this->_gall_nav || !this->_gallobj)
		{
			write_log_info(_log, 1, "NOTE", "no navgation data or atx data for epoch " + crt_epo.str_ymdhms());
			return false;
		}

		// compute reciver time 
		_rec_epo = crt_epo - _crt_rec_clk;
		_crt_obs.addrecTime(_rec_epo);

		// get Rec crd
		bool apply_obj_valid;
		if (_crt_obj->id_type() == t_gdata::REC)
		{
			apply_obj_valid = _apply_rec(crt_epo, _rec_epo, pars);
		}
		else
		{
			apply_obj_valid = _apply_leo(crt_epo, _rec_epo, pars);
		}

		if (!apply_obj_valid)
		{
			write_log_info(_log, 1, "NOTE", "can not apply site in " + crt_epo.str_ymdhms());
			return false;
		}

		// get sat crd
		bool apply_sat_valid = _apply_sat(_rec_epo, _sat_epo);
		if (!apply_sat_valid)
		{
			write_log_info(_log, 1, "NOTE", "can not apply sat in " + crt_epo.str_ymdhms());
			return false;
		}
		_crt_obs.addsatTime(_sat_epo);

#ifdef DEBUG
		cout << "site crs" << scientific << setw(25) << setprecision(15) << _crdSiteCrs[0] << "  " << _crdSiteCrs[1] << "  " << _crdSiteCrs[2] << "  " << endl;
		cout << "sat crs" << scientific << setw(25) << setprecision(15) << _crdSatCrs[0] << "  " << _crdSatCrs[1] << "  " << _crdSatCrs[2] << "  " << endl;
#endif // DEBUG21

		// get f01 PCO
		this->_crs_rec_pco = this->_crs_rec_crd;
		this->_crs_sat_pco = this->_crs_sat_crd;

		shared_ptr<t_gobj>  sat_obj = this->_gallobj->obj(_crt_sat);
		shared_ptr<t_gobj>  rec_obj = this->_gallobj->obj(_crt_rec);

		shared_ptr<t_gpcv>  sat_pcv = (sat_obj != 0) ? sat_obj->pcv(crt_epo) : nullptr;
		shared_ptr<t_gpcv>  rec_pcv = (rec_obj != 0) ? rec_obj->pcv(crt_epo) : nullptr;

		if (!_isCalPCO)
		{
			sat_pcv = nullptr;
			rec_pcv = nullptr;
		}


		if (sat_pcv && !_realtime) //add _realtime by xiongyun 
		{
			// Satellite phase center offset
			t_gtriple pco(0, 0, 0);
			if (sat_pcv->pcoS(_crt_obs, pco, _observ, _band_index[_crt_sys][FREQ_1], _band_index[_crt_sys][FREQ_2]) > 0)
			{
				Matrix _rot_matrix = _RotMatrix_Ant(_crt_obs, _sat_epo, sat_obj, true);
				this->_crs_sat_pco += t_gtriple(_rot_matrix*pco.crd_cvect());
			}

		}

		if (rec_pcv)
		{
			// Receiver phase center offset
			t_gtriple pco(0.0, 0.0, 0.0);
			if (rec_pcv->pcoR(_crt_obs, pco, _observ, _band_index[_crt_sys][FREQ_1], _band_index[_crt_sys][FREQ_2]) > 0)
			{
				Matrix _rot_matrix = _RotMatrix_Ant(_crt_obs, _rec_epo, rec_obj, true);
				this->_crs_rec_pco += t_gtriple(_rot_matrix*pco.crd_cvect());
			}
		}
		
		// compute satclk[s]
		double sat_clk = _crt_sat_clk;

		// compute reldelay[m]
		double reldelay = 0.0;
		bool rel_valid = apply_reldelay(this->_crs_rec_pco, this->_crs_rec_vel, this->_crs_sat_pco, this->_crs_sat_vel, reldelay);
		if (!rel_valid)
		{
			return false;
		}
		_crt_obs.addclk(sat_clk * CLIGHT - reldelay); // include reldelay
		_crt_obs.addreldelay(reldelay);               // for debug

		// set sat parameter PXSAT in lsq
		int sat_index = pars.getParam("", t_gpar::PXSAT, _crt_sat);
		//obsdata.addsatindex(par[sat_index].index);
		if (sat_index >= 0) 
		{
			_sat_index = pars.getParIndex(sat_index);

			// get unit_vector in TRS
			ColumnVector unit_rec2sat(3);

			// change TRS to CRS
			unit_rec2sat = _crs_sat_pco.crd_cvect() - _crs_rec_pco.crd_cvect();
			unit_rec2sat = unit_rec2sat / unit_rec2sat.NormFrobenius();

			// order: P X0 Y0 Z0 VX0..VZ0 solarPar
			_sat_partial *= unit_rec2sat * 1E3;
		}

		int eop_index = pars.getParam("", t_gpar::XPOLE, "");
		if (_grec->id_type() != t_gdata::REC_LEO && eop_index >= 0 ) {

			// get unit_vector in TRS
			ColumnVector unit_rec2sat(3);

			unit_rec2sat = _crs_sat_pco.crd_cvect() - _crs_rec_pco.crd_cvect();
			unit_rec2sat = unit_rec2sat / unit_rec2sat.NormFrobenius();

			// get drotdxpole
			Matrix drdxpole = _dr_dxpole;
			drdxpole = drdxpole / RAD2SEC;
			// get drotdypole
			Matrix drdypole = _dr_dypole;
			drdypole = drdypole / RAD2SEC;
			// get drotdut1
			Matrix drdut1 = _dr_dut1;
			drdut1 = drdut1 / RAD2TSEC;
			double diff = (_rec_epo - pars[eop_index].beg) / 86400.0;

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

			_dr_derp.ReSize(1, 6);
			_dr_derp = unit_rec2sat.t() * dxdpoleut1;
		}

		// addrho
		double tmp = (_crs_sat_crd - _crs_rec_crd).norm();
		_crt_obs.addrho(tmp);

		// add drate
		_crt_obs.adddrate((DotProduct((_crs_sat_vel - _crs_rec_vel).crd_cvect(),
			(_crs_sat_pco - _crs_rec_pco).crd_cvect())) / (CLIGHT * tmp));//yjqin

		// add azim && elev
		//t_gtriple xyz_rho = _crdSatTrs - _crdSiteTrs;
		t_gtriple xyz_rho = _crs_sat_pco - _crs_rec_pco;
		t_gtriple ell_r, neu_s;

		Matrix _rot_matrix = _RotMatrix_Ant(_crt_obs, _rec_epo, _crt_obj, true);
		neu_s = t_gtriple(_rot_matrix.t()*xyz_rho.crd_cvect());

		double NE2 = neu_s[0] * neu_s[0] + neu_s[1] * neu_s[1];
		double ele = acos(sqrt(NE2) / _crt_obs.rho());
		/*if (neu_s[2] < 0.0)
		{
			ele *= -1.0;
		}*/

		if (sqrt(NE2) / _crt_obs.rho() > 1.0)
		{
			_crt_obs.addele(0.0);
		}
		else
		{
			_crt_obs.addele(ele);
		}

		double azi = atan2(neu_s[1], neu_s[0]);
		if (azi < 0)
		{
			azi += 2 * G_PI;
		}
		_crt_obs.addazi(azi);


		///add for another elev and azi used in calculating weight matric
		t_gtriple xyz_rh = _trs_sat_crd - _trs_rec_crd;
		t_gtriple ell_(0, 0, 0), neu_sa(0, 0, 0), xRec(0, 0, 0), xyz_s(0, 0, 0);
		xyz2ell(_trs_rec_crd, ell_, false);
		xyz2neu(ell_, xyz_rh, neu_sa);
		double rho0 = sqrt(pow(_trs_rec_crd[0] - xyz_s[0], 2) + pow(_trs_rec_crd[1] - xyz_s[1], 2) + pow(_trs_rec_crd[2] - xyz_s[2], 2));
		double dPhi = OMEGA * rho0 / CLIGHT;
		xRec[0] = _trs_rec_crd[0] * cos(dPhi) - _trs_rec_crd[1] * sin(dPhi);
		xRec[1] = _trs_rec_crd[1] * cos(dPhi) + _trs_rec_crd[0] * sin(dPhi);
		xRec[2] = _trs_rec_crd[2];
		double NE2_ = neu_sa[0] * neu_sa[0] + neu_sa[1] * neu_sa[1];
		double ele_ = acos(sqrt(NE2_) / _crt_obs.rho());
		/*if (neu_s[2] < 0.0)
		{
			ele_ *= -1.0;
		}*/ //commented by zhangwei
		_crt_obs.addele_leo(ele_);

		// check elevation cut-off
		if (_crt_obj->id_type() == t_gdata::REC && _crt_obs.ele_deg() < _minElev)
		{
			if (_log)
			{
				_log->comment(1, "t_gprecisemodel::prepareObs", "Prepare fail! the elevation is too small");
			}
			return false;
		}
		if (_crt_obj->id_type() == t_gdata::REC_LEO && ele_ * 180.0 / G_PI < _minElevLeo) // for LEO satellite
		{
			_factorP = 1.0E12;
		}
		return true;
	}

	bool t_gprecisemodel::_update_obs_info(const t_gtime & epoch, t_gsatdata & obsdata, t_gallpar& pars)
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
		_trs_sat_crd = obsdata.satcrd();//lvhb added in 20200912

		bool rec_clk_valid = _update_obj_clk("rec"+_crt_rec, epoch, pars, _crt_rec_clk, _obj_clk);
		bool sat_clk_valid = _update_obj_clk("sat"+_crt_sat, epoch, pars, _crt_sat_clk, _obj_clk);
		if (!rec_clk_valid || !sat_clk_valid)
		{
			write_log_info(_log, 1, "NOTE", "no rec or sat clk for epoch " + epoch.str_ymdhms());
			return false;
		}
		return true;
	}

	bool t_gprecisemodel::_update_obs_info(t_gsatdata & obsdata)
	{
		// smgcz
		obsdata = _crt_obs;
		//obsdata.addwind(_crt_obs.wind());
		return true;
	}

	void t_gprecisemodel::update_obj_clk(const string & obj, const t_gtime & epo, double clk)
	{
		_obj_clk[obj].first = epo;
		_obj_clk[obj].second = clk;
		_rec_clk[obj] = clk;
	}

	bool t_gprecisemodel::_update_obj_clk(const string& obj, const t_gtime& crt_epoch,
		t_gallpar& pars, double& clk, map<string, pair<t_gtime, double>>& obj_clk)
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
		auto nppmodel = dynamic_cast<t_gsetproc*>(_settings)->npp_model();
		int base_size = dynamic_cast<t_gsetgen*>(_settings)->list_base().size();
		// update clk
		if (obj_clk[name].first != crt_epoch)
		{
			int clk_valid = _gall_nav->clk(name, crt_epoch, &clk, &clk_rms, &dclk);
			if (clk_valid < 0) clk = 0.0;
			if (_realtime && 
				type == "sat" && 
				(base_size == 0) && (nppmodel != NPP_MODEL::URTK)&&
				!_get_sat_clk_corr(crt_epoch, name, clk, dclk))
			{
				write_log_info(_log, 1, "NOTE", "no sat" + name + "  clk correction for " + crt_epoch.str_ymdhms());
				return false;
			} // zhshen

			obj_clk[name].first = crt_epoch;
			obj_clk[name].second = clk;
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
			idx_clk = pars.getParam(name, t_gpar::CLK, "");
		}
		else if (type == "sat")
		{
			idx_clk = pars.getParam("", t_gpar::CLK_SAT, name);
		}
		else
		{
			idx_clk = -1;
		}

		if (double_eq(clk, 0.0) && idx_clk < 0)
		{
			write_log_info(_log, 1, "NOTE", "no rec clk for " + crt_epoch.str_ymdhms());
			return false;
		}
		//-- no clk file but clk par exist
		if (double_eq(clk, 0.0) && idx_clk >= 0)
		{
			clk = pars[idx_clk].value() / CLIGHT;
		}
		//-- update clk par
		if (idx_clk >= 0)
		{
			pars[idx_clk].value(clk * CLIGHT);
		}

		return true;
	}

	bool t_gprecisemodel::_omc_obs_ALL(const t_gtime& crt_epo, t_gallpar& pars, t_gobs& gobs, double& omc)
	{
		double obs = 0.0;
		if (gobs.is_code())
		{
			obs = _crt_obs.obs_C(gobs);
			//apply_dcb(_gallbias, _crt_obs, obs, _band_index[_crt_sys], _observ, _freq_model, _dcb_model, &gobs, &gobs);
			//correct_DCB("CAS", _gallbias, _crt_obs, obs, _observ, &gobs);
			//obs = obs + ApplyFreqDcb("CAS", _observ, _gallbias, _crt_obs, &gobs, &gobs);
			// 暂时只考虑gps的频内偏差
			//double raw_dcb = intraFreqDcb("CAS", _gallbias, _crt_obs, &gobs);
			//obs = obs + raw_dcb;
			//obs = obs - raw_dcb;
			// already correct dcb

		}
		else if (gobs.is_phase())
		{
			obs = _crt_obs.obs_L(gobs);
		}
		else
		{
			throw logic_error("only support code and phase observation.");
		}

		// jdhuang : if obs = 0.0, return false
		if (double_eq(obs, 0.0))
		{
			_log->logDebug("t_gprecisemodel", "_omc_obs_ALL", "obs type is empty, the obs is zero!");
			return false;
		}

		t_gtime epo = crt_epo;
		double ModelObs = cmpObs(epo, pars, _crt_obs, gobs);
		if (ModelObs < 0)
		{
			_log->logDebug("t_gprecisemodel", "_omc_obs_ALL", "Compute ModelObs error!");
			return false;  //glfeng
		}

		double Amb = 0.0;
		omc = (obs - ModelObs - Amb);

		return true;
	}

	// jdhuang fix
	bool t_gprecisemodel::_wgt_obs_ALL(t_gobs& gobs1, double factorP, double& wgt)
	{
		GOBSTYPE type = gobs1.type();
		if (type != TYPE_C &&
			type != TYPE_P &&
			type != TYPE_L)
		{
			_log->logError("t_gprecisemodel", "_wgt_obs", " type should be TYPE_C/TYPE_P/TYPE_L");
			return false;
		}

		double factor = 0.0;
		double sigRange = 0.0;
		double sigPhase = 0.0;

		// jdhuang : fix the bug
		factor = factorP;
        if (_crt_sys == GSYS::BDS) factor = 2.0;	
		//lijie:if geo
		//if (_crt_sys == GSYS::BDS && _crt_sat <= "C05") factor = 5.0;          // same with panda glfeng
		if (_crt_sys == GSYS::BDS && t_gsys::bds_geo(_crt_sat)) factor = 5.0;          

		// get the sys sigRange
		if (type == TYPE_C || type == TYPE_P)
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
		auto b1 = gobs1.band();
		if (type == TYPE_C || type == TYPE_P)  
		{
			sigRange = sigRange /** 2*/;
		}
		else if (type == TYPE_L)
		{
			sigPhase = sigPhase * _crt_obs.wavelength(b1) /** 2*/;
		}
		else
		{
			sigRange = 0.0;
			sigPhase = 0.0;
		}

		if (sigPhase == 0.0 && sigRange == 0.0) return false;

		switch (type)
		{
		case TYPE_L:
			if (sigPhase == 0.0) return false;
			wgt = 1.0 / pow(factor*sigPhase, 2);
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
		return true;
	}

	bool t_gprecisemodel::_wgt_obs_ALL(t_gobs& gobs1, double factorP, double& wgt, double snr_value)
	{
		GOBSTYPE type = gobs1.type();
		if (type != TYPE_C &&
			type != TYPE_P &&
			type != TYPE_L)
		{
			_log->logError("t_gprecisemodel", "_wgt_obs", " type should be TYPE_C/TYPE_P/TYPE_L");
			return false;
		}

		double factor = 0.0;
		double sigRange = 0.0;
		double sigPhase = 0.0;

		// jdhuang : fix the bug
		factor = factorP;
		if (_crt_sys == GSYS::BDS) factor = 2.0;
		//lijie:if geo
		//if (_crt_sys == GSYS::BDS && _crt_sat <= "C05") factor = 5.0;          // same with panda glfeng
		if (_crt_sys == GSYS::BDS && t_gsys::bds_geo(_crt_sat)) factor = 5.0;

		// get the sys sigRange
		if (type == TYPE_C || type == TYPE_P)
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
			case SNR:
				factor *= sqrt(134.02 * pow(10, -(snr_value / 17.91))); 
				break;
			default:
				cerr << "gspplsq: we can't deal with this WeightObs method!";
				return false;
			}
		}

		// get the combine obs factor
		// TODO finish other weight method 
		auto b1 = gobs1.band();
		if (type == TYPE_C || type == TYPE_P)
		{
			sigRange = sigRange /** 2*/;
		}
		else if (type == TYPE_L)
		{
			sigPhase = sigPhase * _crt_obs.wavelength(b1) /** 2*/;
		}
		else
		{
			sigRange = 0.0;
			sigPhase = 0.0;
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
		return true;
	}



	bool t_gprecisemodel::_prt_obs(const t_gtime & epoch, t_gallpar& pars, t_gobs & gobs, vector<pair<int, double>>& coeff)
	{
		t_gtriple groundEll;
		xyz2ell(_trs_rec_crd, groundEll, false);
		int parNum = pars.parNumber();

		auto par_list = pars.getPartialIndex(_crt_rec, _crt_sat);
		//for (int ipar = 0; ipar < parNum; ipar++)
		for (int ipar: par_list)
		{
			//double coeff_value = _Partial(epoch, _crt_obs, gobs1, gobs2, pars.getPar(ipar));
			double coeff_value = _Partial(epoch, _crt_obs, gobs, pars.getPar(ipar));
			if (coeff_value != 0.0)
			{
				//note: This is bad from 1 in coeff vector
				coeff.push_back(make_pair(ipar+1, coeff_value));
			}
		}
		return true;
	}




	int t_gprecisemodel::outlierDetect(vector<t_gsatdata>& data)
	{
		vector<t_gsatdata>::iterator itMaxVcodeNORM = data.end();
		vector<t_gsatdata>::iterator itMaxVphaseNORM = data.end();

		//vector<t_gsatdata>::iterator itMaxVcodeORIG = data.end();
		//vector<t_gsatdata>::iterator itMaxVphaseORIG = data.end();

		vector<t_gsatdata>::iterator itDataErase = data.end();

		double maxVcodeNORM = 0.0;
		double maxVphaseNORM = 0.0;

		double maxVcodeORIG = 0.0;
		double maxVphaseORIG = 0.0;

		// find maximal code/phase residuals
		maxVcodeNORM = _maxres(false, data, itMaxVcodeNORM, RES_NORM);
		maxVphaseNORM = _maxres(true, data, itMaxVphaseNORM, RES_NORM);

		int maxVphaseNORM_loc = itMaxVphaseNORM - data.begin() + 1;
		int maxVcodeNORM_loc = itMaxVcodeNORM - data.begin() + 1;


		//maxVcodeORIG = _maxres(false, data, itMaxVcodeORIG, RES_ORIG);
		//maxVphaseORIG = _maxres(true, data, itMaxVphaseORIG, RES_ORIG);


		if (_check_outl(true, maxVphaseNORM, itMaxVphaseNORM, maxVphaseORIG, itMaxVphaseNORM, itDataErase, data)) { data.erase(itDataErase); return maxVphaseNORM_loc; }
		if (_check_outl(false, maxVcodeNORM, itMaxVcodeNORM, maxVcodeORIG, itMaxVcodeNORM, itDataErase, data)) { data.erase(itDataErase);  return maxVcodeNORM_loc; }

		return 0;

	}


//	double t_gprecisemodel::cmpObs(t_gtime& epoch, t_gallpar& param, t_gsatdata& obsdata, t_gobs& gobs, bool com)
//	{
//
//		//return (t_gpppmodel::cmpObs(epoch, param, obsdata, gobs, com));
//
//		gtrace("t_gpppmodel::cmpObs");
//
//		double spp_model = t_gsppmodel::cmpObs(epoch, param, obsdata, gobs);
//		if (spp_model < 0) return -1;
//
//		// Cartesian coordinates to ellipsodial coordinates
//		t_gtriple ell;
//		xyz2ell(_trs_rec_crd, ell, false);
//
//		// Wind up correction 
//		double wind = 0.0;
//		if (gobs.is_phase())
//		{
//			double wavelength = 0.0;
//			if (_observ == IONO_FREE)
//			{
//				GSYS gs = obsdata.gsys();
//				if (!_isTriple)
//				{
//					GOBSBAND b1 = t_gsys::band_priority(gs, FREQ_1);
//					GOBSBAND b2 = t_gsys::band_priority(gs, FREQ_2);
//					wavelength = obsdata.wavelength_L3(b1, b2);
//				}
//				if (_frequency == 3 && _isTriple)
//				{
//					GOBSBAND b1 = t_gsys::band_priority(gs, FREQ_1);
//					GOBSBAND b3 = t_gsys::band_priority(gs, FREQ_3);
//					wavelength = obsdata.wavelength_L3(b1, b3);
//				}
//			}
//			else
//			{
//				wavelength = obsdata.wavelength(gobs.band());
//			}
//
//			if (fabs(obsdata.wind()) > 0) wind = obsdata.wind() * wavelength;
//			else wind = windUp(obsdata, _trs_rec_crd.crd_cvect()) * wavelength;
//		}
//
//		// BDS2/QZSS Beta < 4deg, not add to equation!
//		if (_crt_sys == GSYS::BDS && obsdata.sat() <= "C15" || _crt_sys == GSYS::QZS)
//		{
//			if(fabs(obsdata.beta()) <= 4.0*D2R) return -1;
//		}
//
//		// Phase center variation correction
//		double pcv_R = 0.0;
//		double pcv_S = 0.0;
//		double pco_R = 0.0;
//		double pco_S = 0.0;
//		if (_gallobj)
//		{
//			shared_ptr<t_gobj>  sat_obj = _gallobj->obj(_crt_sat);
//			shared_ptr<t_gobj>  rec_obj = _gallobj->obj(_crt_rec);
//
//			shared_ptr<t_gpcv>  sat_pcv;
//			shared_ptr<t_gpcv>  rec_pcv;
//
//			if (sat_obj != 0)  sat_pcv = sat_obj->pcv(epoch);
//			if (rec_obj != 0)  rec_pcv = rec_obj->pcv(epoch);
//
//			GOBS_LC lc = LC_IF;
//			if (_observ != IONO_FREE)
//			{
//				if (t_gsys::band2freq(obsdata.gsys(), gobs.band()) == FREQ_1) lc = LC_L1;
//				if (t_gsys::band2freq(obsdata.gsys(), gobs.band()) == FREQ_2) lc = LC_L2;
//				if (t_gsys::band2freq(obsdata.gsys(), gobs.band()) == FREQ_3) lc = LC_L3;
//				if (t_gsys::band2freq(obsdata.gsys(), gobs.band()) == FREQ_4) lc = LC_L4;
//				if (t_gsys::band2freq(obsdata.gsys(), gobs.band()) == FREQ_5) lc = LC_L5;
//			}
//
//			GOBSBAND b1 = t_gsys::band_priority(obsdata.gsys(), FREQ_1);
//			GOBSBAND b2 = t_gsys::band_priority(obsdata.gsys(), FREQ_2);
//			if (sat_pcv && com)
//			{
//				// Satellite phase center variation
//				// -- Satellite phase center offset
//				t_gtriple pco(0, 0, 0);
//				if (sat_pcv->pcoS(obsdata, pco, lc, b1, b2) > 0)
//				{
//					string antenna = sat_pcv->anten();
//					t_gtriple dx(0, 0, 0);
//					ColumnVector i(3);
//					ColumnVector j(3);
//					ColumnVector k(3);
//					if (_attitudes == YAW_NOMI)	   attitude(obsdata, "", i, j, k);				//nominal modeling
//					else if (_attitudes == YAW_RTCM) attitude(obsdata, obsdata.yaw(), i, j, k);	//value from RTCM used
//					else							   attitude(obsdata, antenna, i, j, k);		//default
//					dx[0] = pco[0] * i(1) + pco[1] * j(1) + pco[2] * k(1);
//					dx[1] = pco[0] * i(2) + pco[1] * j(2) + pco[2] * k(2);
//					dx[2] = pco[0] * i(3) + pco[1] * j(3) + pco[2] * k(3);
//					sat_pcv->pco_proj(pco_S, obsdata, _trs_rec_crd, dx);
//					obsdata.addpco(dx);
//				}
//
//			}
//
//			if (rec_pcv != 0)
//			{
//				// Receiver phase center offset
//				t_gtriple pco(0.0, 0.0, 0.0);
//				if (rec_pcv->pcoR(obsdata, pco, lc, b1, b2) > 0)
//				{
//					Matrix _rot_matrix = _RotMatrix_Ant(obsdata, _rec_epo, rec_obj, false);
//					t_gtriple dx(_rot_matrix*(pco.crd_cvect()));
//					rec_pcv->pco_proj(pco_R, obsdata, _trs_rec_crd, dx);
//				}
//				pco_R *= -1;
//			}
//
//			if (gobs.is_phase())
//			{
//				if (sat_pcv != 0 && com)
//				{
//					// Satellite phase center variation
//					//sat_pcv->pcvS(pcv_S, obsdata, _trs_rec_crd, LC_IF);
//					sat_pcv->pcvS(pcv_S, obsdata, _trs_rec_crd, lc, b1, b2);
//				}
//				if (rec_pcv != 0)
//				{
//					// Receiver phase center variation
//					//rec_pcv->is_noazi(true);
//					rec_pcv->pcvR(pcv_R, obsdata, lc, b1, b2);
//				}
//			}
//		}
//
//		// 用来调试的，别删了
//		//cout << setw(20) << " EPOCH is : " << epoch.str_mjdsod() << " " << setw(6) << obsdata.sat() << setw(6) << obsdata.site() 
//		//	 << setw(20) << " Band is : " << setw(6) << gobs.band() 
//		//	 << setw(20) << " spp_model  " << fixed << left << setw(20) << setprecision(6) << spp_model
//		//	 << setw(20) << " wind       " << fixed << left << setw(20) << setprecision(6) << wind
//		//	 << setw(20) << " pco_R      " << fixed << left << setw(20) << setprecision(6) << pco_R
//		//	 << setw(20) << " pco_S      " << fixed << left << setw(20) << setprecision(6) << pco_S
//		//	 << setw(20) << " pcv_R      " << fixed << left << setw(20) << setprecision(6) << pcv_R
//		//	 << setw(20) << " pcv_S      " << fixed << left << setw(20) << setprecision(6) << pcv_S
//		//	 << endl;
//
//#ifdef DEBUG_GLO
//		cout << fixed << setprecision(5) << setw(10) << pco_R << setw(10) << pco_S << setw(10) << pcv_R << setw(10) << pcv_S
//			<< setw(10) << wind;
//#endif
//		return spp_model +
//			wind +
//			pco_R +
//			pco_S +
//			pcv_R +
//			pcv_S;
//	}

	double t_gprecisemodel::cmpObs(t_gtime& epoch, t_gallpar& param, t_gsatdata& obsdata, t_gobs& gobs)
	{
		gtrace("t_gpppmodel::cmpObs");
		if (_gallobj == nullptr)
		{
			throw logic_error("gallobj is nullptr in precisemodel cmpObs");
		}

		double spp_model = t_gsppmodel::cmpObs(epoch, param, obsdata, gobs);
		if (spp_model < 0) return -1;

		// Cartesian coordinates to ellipsodial coordinates
		t_gtriple ell;
		xyz2ell(_trs_rec_crd, ell, false);

		// band used
		GOBSBAND band = gobs.band();

		// Wind up correction 
		double wind = 0.0;
		if (gobs.is_phase())
		{
			double wavelength = obsdata.wavelength(gobs.band());
			if (fabs(obsdata.wind()) > 0)
			{
				wind = obsdata.wind() * wavelength;
			}
			else
			{
				wind = windUp(obsdata, _trs_rec_crd.crd_cvect()) * wavelength;
			}
		}

		// ion correction
		double ion = 0.0;
		auto band_1 = _band_index[obsdata.gsys()][FREQ_1];
		ion = ionDelay(epoch, param, obsdata, _ion_model, band_1, gobs);

		// Phase center variation correction
		double pcv_R = 0.0;
		double pcv_S = 0.0;
		double pco_R = 0.0;
		double pco_S = 0.0;

		if (_isCalPCO) {//lvhb modified in 20200918

			shared_ptr<t_gobj>  sat_obj = _gallobj->obj(_crt_sat);
			shared_ptr<t_gobj>  rec_obj = _gallobj->obj(_crt_rec);

			shared_ptr<t_gpcv>  sat_pcv = (sat_obj != nullptr) ? sat_obj->pcv(epoch) : nullptr;
			shared_ptr<t_gpcv>  rec_pcv = (rec_obj != nullptr) ? rec_obj->pcv(epoch) : nullptr;

			if (sat_pcv != nullptr && !_realtime)//add _realtime by xiongyun
			{
				// Satellite phase center variation
				// -- Satellite phase center offset
				t_gtriple pco(0, 0, 0);
				if (sat_pcv->pcoS_raw(obsdata, pco, band) > 0)
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
					sat_pcv->pco_proj(pco_S, obsdata, _trs_rec_crd, dx);
					obsdata.addpco(dx);
				}
			}

			if (rec_pcv != nullptr)
			{
				// Receiver phase center offset
				t_gtriple pco(0.0, 0.0, 0.0);
				if (rec_pcv->pcoR_raw(obsdata, pco, band) > 0)
				{
					Matrix _rot_matrix = _RotMatrix_Ant(obsdata, _rec_epo, rec_obj, false);
					t_gtriple dx(_rot_matrix * (pco.crd_cvect()));
					rec_pcv->pco_proj(pco_R, obsdata, _trs_rec_crd, dx);
				}
				pco_R *= -1;
			}

			if (gobs.is_phase())
			{
				if (sat_pcv != 0 && (!_realtime || _ultrasp3))//add _realtime by xiongyun
				{
					// Satellite phase center variation
					sat_pcv->pcvS_raw(pcv_S, obsdata, band, _trs_rec_crd);
				}
				if (rec_pcv != 0)
				{
					// Receiver phase center variation
					rec_pcv->pcvR_raw(pcv_R, obsdata, band);
				}
			}
		}

		// 用来调试的，别删了
	//	cout << setw(20) << " EPOCH is : " << epoch.str_mjdsod() << " " << setw(6) << obsdata.sat() << setw(6) << obsdata.site() 
	//		 << setw(20) << " Band is : " << setw(6) << gobs.band() 
	//		 << setw(20) << " spp_model  " << fixed << left << setw(20) << setprecision(6) << spp_model
	//		 << setw(20) << " wind       " << fixed << left << setw(20) << setprecision(6) << wind
	//		 << setw(20) << " pco_R      " << fixed << left << setw(20) << setprecision(6) << pco_R
	//		 << setw(20) << " pco_S      " << fixed << left << setw(20) << setprecision(6) << pco_S
	//		 << setw(20) << " pcv_R      " << fixed << left << setw(20) << setprecision(6) << pcv_R
	//		 << setw(20) << " pcv_S      " << fixed << left << setw(20) << setprecision(6) << pcv_S
	//		 << endl;

#ifdef DEBUG_GLO
		cout << fixed << setprecision(5) << setw(10) << pco_R << setw(10) << pco_S << setw(10) << pcv_R << setw(10) << pcv_S
			<< setw(10) << wind;
#endif
		return spp_model +
			wind +
			pco_R +
			pco_S +
			pcv_R +
			pcv_S + 
			ion;
	}

	double t_gprecisemodel::getZHD(const string& site, const t_gtime& epo)
	{
		t_gtriple xyz, ell;
		_grec = _gallobj->obj(site);
		xyz = _grec->crd(epo);
		xyz2ell(xyz, ell, false);
		double ZHD = this->tropoModel()->getZHD(ell, epo);
		return ZHD;

	}
	double t_gprecisemodel::getZWD(const string& site, const t_gtime& epo)
	{
		t_gtriple xyz, ell;
		_grec = _gallobj->obj(site);
		xyz = _grec->crd(epo);
		xyz2ell(xyz, ell, false);
		double ZWD = this->tropoModel()->getZWD(ell, epo);
		return ZWD;
	}

	double t_gprecisemodel::tropoDelay(t_gtime& epoch, t_gallpar& param, t_gtriple site_ell, t_gsatdata& satdata)
	{
		gtrace("t_gprecisemodel::tropoDelay");

		if (site_ell[2] > 1E4)
		{
			return 0.0;
		}
		else
		{
			if (_tropoModel == 0)
			{
				if (_log)
				{
					_log->comment(0, "gppp", "Tropo Model setting is not correct. Default used! Check config.");
				}
				else
				{
					cerr << "gppp - Tropo Model setting is not correct. Default used! Check config.\n";
				}
				_tropoModel = make_shared<t_saast>();
			}

			double ele = satdata.ele();
			double azi = satdata.azi();

			double delay = 0.0;
			double zwd = 0.0;
			double zhd = 0.0;
			int i, j, k;
			t_gtriple ell;

			xyz2ell(_trs_rec_crd, ell, false);

			if (ell[2] > 1E4) {
				return 100.0;
			}

			i = param.getParam(_site, t_gpar::TRP, "");
			j = param.getParam(_site, t_gpar::GRD_N, "");
			k = param.getParam(_site, t_gpar::GRD_E, "");

			if (i >= 0)
			{
				zwd = param[i].value();
				
				//zhd = param[i].apriori();
				//zhd = param[i].zhd;
				// update zwd if zwd par is new
				// change by wangbo for PPP-RTK
				if (param[i].apriori()>1E-4 && (zwd == 0.0 || epoch == param[i].beg))
				{
					zwd = _tropoModel->getZWD(ell, epoch);
					param[i].value(zwd);
					
				}
				zhd = _tropoModel->getZHD(ell, epoch);
					param[i].zhd = zhd;
				//zhd = param[i].apriori();
			}
			else
			{
				if (_tropoModel != 0)
				{
					zwd = _tropoModel->getZWD(ell, epoch);
					zhd = _tropoModel->getZHD(ell, epoch);
				}
			}

			double mfh, mfw, dmfh, dmfw;
			mfh = mfw = dmfh = dmfw = 0.0;
			if (_tropo_mf == GMF)
			{
				t_gmf mf;
				mf.gmf(epoch.mjd(), ell[0], ell[1], ell[2], G_PI / 2.0 - ele,
					mfh, mfw, dmfh, dmfw);
			}
			else if (_tropo_mf == COSZ)
			{
				mfh = mfw = 1 / sin(ele);
				dmfh = dmfw = -(cos(ele)) / (sin(ele) * sin(ele));
			}
			else return 0.0;

			satdata.addmfH(mfh);
			satdata.addmfW(mfw);

			delay = mfh * zhd + mfw * zwd;

#ifdef DEBUG
			//debug by ZHJ
			cout << fixed << setprecision(5)
				<< "ele" << satdata.ele_deg() << endl
				<< "dry delay " << zhd << endl
				<< "dry map" << mfh << endl
				<< "wet delay" << zwd << endl
				<< "wet map" << mfw << endl;
#endif // !DEBUG

			double grdN, grdE;
			grdN = grdE = 0.0;

			if (j >= 0 && k >= 0) {
				if (_grad_mf == TILTING) {
					grdN = param[j].value() * dmfw * cos(azi);
					grdE = param[k].value() * dmfw * sin(azi);
					satdata.addmfG(dmfw);
				}
				else if (_grad_mf == CHEN_HERRING) {
					double mfg = 1.0 / (sin(ele)*tan(ele) + 0.0032);
					grdN = param[j].value() * 1000.0 * mfg * cos(azi);  grdN /= 1000.0;
					grdE = param[k].value() * 1000.0 * mfg * sin(azi);  grdE /= 1000.0;
					satdata.addmfG(mfg);
				}
				else if (_grad_mf == BAR_SEVER) {
					double mfg = mfw * (1 / tan(ele));
					grdN = param[j].value() * mfg * cos(azi);
					grdE = param[k].value() * mfg * sin(azi);
					satdata.addmfG(mfg);
				}

				delay += grdN + grdE;
			}


#ifdef DEBUG_CORR
			cout << epoch.str("EPOCH: %H:%M:%S") << endl << fixed << setprecision(3);
			cout << satdata.sat() << " Ell:" << ell[0] << " " << ell[1] << " " << ell[2]
				<< " Hydrostatic part: " << zhd //_tropoModel->getZHD(ell, epoch)
				<< " Wet part: " << zwd
				<< " mfh: " << mfh
				<< " mfw: " << mfw
				<< " Delay: " << delay << endl << endl;
#endif
			return delay;
		}
	}

	double t_gprecisemodel::ionDelay(t_gtime & epoch, t_gallpar & param, t_gsatdata& satdata, IONMODEL& ion_model, GOBSBAND& band_1, t_gobs& gobs)
	{
		gtrace("t_gsppmodel::ionoDelay");
		if (band_1 == BAND || band_1 == BAND_A || band_1 == BAND_B || band_1 == BAND_C || band_1 == BAND_D)
		{
			return 0.0;
		}

		if (gobs.band() == BAND || gobs.band() == BAND_A || gobs.band() == BAND_B || gobs.band() == BAND_C || gobs.band() == BAND_D)
		{
			return 0.0;
		}

		if (ion_model == IONMODEL::VION)
		{
			throw logic_error(" not support VION model, use SION in XML");
		}

		double iono_delay = 0.0;

		// jdhuang : from G01_F
		double f1 = satdata.frequency(band_1);;
		double fk = satdata.frequency(gobs.band());
		double alfa = 0.0;

		if (gobs.is_phase())
		{
			alfa = -(f1 * f1) / (fk * fk);
		}
		if (gobs.is_code())
		{
			alfa = (f1 * f1) / (fk * fk);
		}

		// ionosphere slant delay parameter
		int i = param.getParam(_site, t_gpar::SION, satdata.sat());
		if (i >= 0)
		{
			iono_delay = alfa * param[i].value();
		}

#ifdef DEBUG
		cout << "IONO delay: " << gsatdata.sat() << " " << epoch.str_ymdhms() << " iono delay = " << iono_delay << endl;
		//    int ooo; cin >> ooo;
#endif

		return iono_delay;
	}

	double t_gprecisemodel::ifbCorrection(t_gallpar & param, string & sat, string & rec, FREQ_MODEL & freq_model, IFB_MODEL & ifb_model, FREQ_SEQ& freq, t_gobs & gobs)
	{
		// Inter frequency clocks bias
		double ifb_offset = 0.0;
		int i = -1;
		switch (ifb_model)
		{
		case IFB_MODEL::EST_IFB:
		{
			if (freq == FREQ_1) return 0.0;
			if (freq == FREQ_2) return 0.0;
			if (freq == FREQ_3) i = param.getParam(rec, t_gpar::IFB_C3, sat);
			if (freq == FREQ_4) i = param.getParam(rec, t_gpar::IFB_C4, sat);
			if (freq == FREQ_5) i = param.getParam(rec, t_gpar::IFB_C5, sat);
			break;
		}
		case IFB_MODEL::EST_REC_IFB:
		{
			auto gsys = t_gsys::sat2gsys(sat);
			if (freq == FREQ_1) return 0.0;
			if (freq == FREQ_2) return 0.0;
			if (freq == FREQ_3)
			{
				if(gsys == GPS) i = param.getParam(rec, t_gpar::GPS_IFB_C3, "");
				if(gsys == GAL) i = param.getParam(rec, t_gpar::GAL_IFB_C3, "");
				if(gsys == BDS) i = param.getParam(rec, t_gpar::BDS_IFB_C3, "");
				if(gsys == QZS) i = param.getParam(rec, t_gpar::QZS_IFB_C3, "");
			}
			if (freq == FREQ_4)
			{
				if (gsys == GPS) i = param.getParam(rec, t_gpar::GPS_IFB_C4, "");
				if (gsys == GAL) i = param.getParam(rec, t_gpar::GAL_IFB_C4, "");
				if (gsys == BDS) i = param.getParam(rec, t_gpar::BDS_IFB_C4, "");
				if (gsys == QZS) i = param.getParam(rec, t_gpar::QZS_IFB_C4, "");
			}
			if (freq == FREQ_5)
			{
				if (gsys == GPS) i = param.getParam(rec, t_gpar::GPS_IFB_C5, "");
				if (gsys == GAL) i = param.getParam(rec, t_gpar::GAL_IFB_C5, "");
				if (gsys == BDS) i = param.getParam(rec, t_gpar::BDS_IFB_C5, "");
				if (gsys == QZS) i = param.getParam(rec, t_gpar::QZS_IFB_C5, "");
			}
			break;
		}
		case IFB_MODEL::EST_SAT_IFB:
		{
			if (freq == FREQ_1) return 0.0;
			if (freq == FREQ_2) return 0.0;
			if (freq == FREQ_3) i = param.getParam("", t_gpar::SAT_IFB_C3, sat);
			if (freq == FREQ_4) i = param.getParam("", t_gpar::SAT_IFB_C4, sat);
			if (freq == FREQ_5) i = param.getParam("", t_gpar::SAT_IFB_C5, sat);
			break;
		}
		default:
			return 0.0;
		}

		if (i >= 0)
		{
			ifb_offset = param[i].value();
		}

		return ifb_offset;
	}

	//double t_gprecisemodel::interFreqDcb(const string& dcb, t_gallbias* allbias, t_gsatdata& satdata, t_gobs* gobs1, t_gobs* gobs2)
	//{

	//	if (allbias == nullptr)                                throw logic_error("dcb data is empty!      ");
	//	if (gobs1 == nullptr || gobs2 == nullptr)              throw logic_error("gob1 or gobs2 is nullptr");

	//	auto gsys = satdata.gsys();
	//	auto gsat = satdata.sat();
	//	auto grec = satdata.site();
	//	auto gepo = satdata.epoch();

	//	t_gobs g1 = (gobs1->attr() == ATTR) ? satdata.id_range(gobs1->band()) : gobs1->gobs();
	//	t_gobs g2 = (gobs2->attr() == ATTR) ? satdata.id_range(gobs2->band()) : gobs2->gobs();

	//	g1.gobs2to3(satdata.gsys());
	//	g2.gobs2to3(satdata.gsys());

	//	double dcb_value = 0.0;

	//	GOBS obs1 = t_gsys::gobs_priority(gsys, g1.gobs());
	//	GOBS obs2 = t_gsys::gobs_priority(gsys, g2.gobs());

	//	if (obs1 == GOBS::X || obs2 == GOBS::X)
	//	{
	//		cout << gobs2str(obs1) + " " + gobs2str(obs2) + " " << " " + gsat + " " << " " + grec + " " << gepo.str_mjdsod("wrong type of GNSS sys, only support GPS and GAL.") << endl;
	//		cerr << gobs2str(obs1) + " " + gobs2str(obs2) + " " << " " + gsat + " " << " " + grec + " " << gepo.str_mjdsod("wrong type of GNSS sys, only support GPS and GAL.") << endl;
	//		_log->logError("t_gprecisemodel", "interFreqDcb", gobs2str(obs1) + " " + gobs2str(obs2) + " " + " " + gsat + " " + " " + grec + " ");
	//		_log->logError("t_gprecisemodel", "interFreqDcb", gepo.str_mjdsod("wrong type of GNSS sys, only support GPS and GAL."));
	//		return 0.0;
	//	}


	//	dcb_value = allbias->get(gepo, gsat, obs1, obs2, dcb);
	//	if (double_eq(dcb_value, 0.0))
	//	{
	//		cout << gobs2str(g1.gobs()) + " " + gobs2str(g2.gobs()) + " " << " " + gsat + " " << " " + grec + " " << endl;
	//		//throw logic_error("error type for code in apply dcb : " + gobs2str(g1.gobs()) + " " + gobs2str(g2.gobs()));
	//		_log->logError("t_gprecisemodel", "interFreqDcb", gobs2str(obs1) + " " + gobs2str(obs2) + " " + " " + gsat + " " + " " + grec + " ");
	//		_log->logError("t_gprecisemodel", "interFreqDcb", gepo.str_mjdsod("wrong type of GNSS sys, only support GPS and GAL."));
	//		return 0.0;
	//	}

	//	return dcb_value;
	//}

	//double t_gprecisemodel::intraFreqDcb(const string& dcb, t_gallbias* allbias, t_gsatdata& satdata, t_gobs* gobs)
	//{
	//	if (allbias == nullptr) throw logic_error("dcb data is empty  !      ");
	//	if (gobs == nullptr) throw logic_error("gob      is nullptr!      ");

	//	t_gobs g1 = (gobs->attr() == ATTR) ? satdata.id_range(gobs->band()) : gobs->gobs();
	//	g1.gobs2to3(satdata.gsys());

	//	auto gsat = satdata.sat();
	//	auto gepo = satdata.epoch();

	//	double dcb_value = 0.0;
	//	if (satdata.gsys() == GSYS::GPS)
	//	{
	//		switch (g1.gobs())
	//		{
	//		case C1C: dcb_value = allbias->get(gepo, gsat, C1W, C1C, dcb); break;
	//		case C1P: dcb_value = 0.0;
	//		case C1Y: dcb_value = 0.0;
	//		case C1W: dcb_value = 0.0;                                     break;
	//		case C2C: dcb_value = allbias->get(gepo, gsat, C2W, C2C, dcb); break;
	//		case C2L: dcb_value = allbias->get(gepo, gsat, C2W, C2L, dcb); break;
	//		case C2X: dcb_value = allbias->get(gepo, gsat, C2W, C2X, dcb); break;
	//		case C2P: dcb_value = 0.0;
	//		case C2Y: dcb_value = 0.0;
	//		case C2W: dcb_value = 0.0;                                     break;
	//		case C5Q: dcb_value = allbias->get(gepo, gsat, C5X, C5Q, dcb); break;
	//		case C5X: dcb_value = 0.0;                                     break;
	//		default:  dcb_value = 0.0;                                     break;
	//		}
	//	}
	//	return dcb_value;
	//}

	//double t_gprecisemodel::cmbIntraFreqDcb(const string& dcb, OBSCOMBIN& obscombin, t_gallbias* allbias, t_gsatdata& satdata, t_gobs* gobs1, t_gobs* gobs2)
	//{
	//	double dcb_value = 0.0;
	//	double intra_dcb_1 = 0.0;
	//	double intra_dcb_2 = 0.0;

	//	auto   gsys = _crt_obs.gsys();
	//	auto   freq_1 = _freq_index[gsys][gobs1->band()];
	//	auto   freq_2 = _freq_index[gsys][gobs2->band()];

	//	// jdhuang :
	//	// add intra freq dcb
	//	intra_dcb_1 = intraFreqDcb(dcb, _gallbias, _crt_obs, gobs1);
	//	intra_dcb_2 = intraFreqDcb(dcb, _gallbias, _crt_obs, gobs2);

	//	//
	//	if (obscombin == MW_COMBIN)
	//	{
	//		// modified for correct P1 code bias glfeng
	//		double l1 = satdata.wavelength(gobs1->band());
	//		double l2 = satdata.wavelength(gobs2->band());

	//		if (double_eq(l2, 0.0))
	//		{
	//			_log->logError("t_gprecisemodel", "ApplyCmbTripleDcb", "wrong type of GNSS FREQ, only support FEQ 1 2 3.");
	//			return 0.0;
	//		}
	//		double fact = l1 / l2;
	//		dcb_value = -(intra_dcb_1 / l1 + intra_dcb_2 / l2) * (1.0 - fact) / (1.0 + fact);

	//		return dcb_value;
	//	}

	//	double alfa12 = 0.0;
	//	double beta12 = 0.0;
	//	satdata.coef_ionofree(_band_index[gsys][FREQ_1], alfa12, _band_index[gsys][FREQ_2], beta12);
	//	dcb_value = intra_dcb_1 * alfa12 + intra_dcb_2 * beta12;

	//	return dcb_value;
	//}

	//double t_gprecisemodel::cmbInterFreqDcb(const string& dcb, OBSCOMBIN& obscombin, t_gallbias* allbias, t_gsatdata& satdata, t_gobs* gobs1, t_gobs* gobs2)
	//{
	//	double dcb_value = 0.0;
	//	double inter_dcb_12 = 0.0;
	//	double inter_dcb_13 = 0.0;
	//	double inter_dcb_23 = 0.0;
	//	double inter_dcb_14 = 0.0;
	//	double inter_dcb_24 = 0.0;
	//	double inter_dcb_15 = 0.0;
	//	double inter_dcb_25 = 0.0;
	//	auto   gsys = _crt_obs.gsys();
	//	auto   freq1 = _freq_index[gsys][gobs1->band()];
	//	auto   freq2 = _freq_index[gsys][gobs2->band()];

	//	t_gobs obs1 = _crt_obs.select_range(_band_index[gsys][FREQ_1], true);
	//	t_gobs obs2 = _crt_obs.select_range(_band_index[gsys][FREQ_2], true);
	//	t_gobs obs3 = _crt_obs.select_range(_band_index[gsys][FREQ_3], true);
	//	t_gobs obs4 = _crt_obs.select_range(_band_index[gsys][FREQ_4], true);
	//	t_gobs obs5 = _crt_obs.select_range(_band_index[gsys][FREQ_5], true);

	//	double alfa12 = 0.0;
	//	double alfa13 = 0.0;
	//	double beta12 = 0.0;
	//	double beta13 = 0.0;
	//	satdata.coef_ionofree(_band_index[gsys][FREQ_1], alfa12, _band_index[gsys][FREQ_2], beta12);
	//	satdata.coef_ionofree(_band_index[gsys][FREQ_1], alfa13, _band_index[gsys][FREQ_3], beta13);


	//	if (obscombin == RAW_ALL)
	//	{
	//		switch (freq1)
	//		{
	//		case FREQ_1:
	//			inter_dcb_12 = interFreqDcb(dcb, _gallbias, _crt_obs, &obs1, &obs2);
	//			dcb_value = + beta12 * inter_dcb_12;
	//			break;
	//		case FREQ_2:
	//			inter_dcb_12 = interFreqDcb(dcb, _gallbias, _crt_obs, &obs1, &obs2);
	//			dcb_value = - alfa12 * inter_dcb_12;
	//			break;
	//		case FREQ_3:
	//			inter_dcb_13 = interFreqDcb(dcb, _gallbias, _crt_obs, &obs1, &obs3);
	//			inter_dcb_23 = interFreqDcb(dcb, _gallbias, _crt_obs, &obs2, &obs3);
	//			dcb_value = -alfa12 * inter_dcb_13 - beta12 * inter_dcb_23;
	//			break;
	//		case FREQ_4:
	//			inter_dcb_14 = interFreqDcb(dcb, _gallbias, _crt_obs, &obs1, &obs4);
	//			inter_dcb_24 = interFreqDcb(dcb, _gallbias, _crt_obs, &obs2, &obs4);
	//			dcb_value = -alfa12 * inter_dcb_14 - beta12 * inter_dcb_24;
	//			break;
	//		case FREQ_5:
	//			inter_dcb_15 = interFreqDcb(dcb, _gallbias, _crt_obs, &obs1, &obs5);
	//			inter_dcb_25 = interFreqDcb(dcb, _gallbias, _crt_obs, &obs2, &obs5);
	//			dcb_value = -alfa12 * inter_dcb_15 - beta12 * inter_dcb_25;
	//			break;
	//		default:
	//			dcb_value = 0.0;
	//			_log->logError("t_gprecisemodel", "ApplyTripleDcb", "wrong type of GNSS FREQ, only support FEQ 1 2 3.");
	//			break;
	//		}
	//	}
	//	else
	//	{
	//		if (freq1 == FREQ_1 && freq2 == FREQ_3)
	//		{
	//			inter_dcb_12 = interFreqDcb(dcb, _gallbias, _crt_obs, &obs1, &obs2);
	//			inter_dcb_13 = interFreqDcb(dcb, _gallbias, _crt_obs, &obs1, &obs3);
	//			dcb_value =  beta12 * inter_dcb_12 - beta13 * inter_dcb_13;
	//		}
	//		else
	//		{
	//			dcb_value = 0.0;
	//		}
	//	}
	//	return dcb_value;
	//}

	//double t_gprecisemodel::ApplyFreqDcb(const string& dcb, OBSCOMBIN& obscombin, t_gallbias* allbias, t_gsatdata& satdata, t_gobs* gobs1, t_gobs* gobs2)
	//{
	//	if (obscombin == RAW_ALL)
	//	{
	//		return ApplyRawTripleDcb(dcb, allbias, satdata, gobs1);
	//	}
	//	else
	//	{
	//		return ApplyCmbTripleDcb(dcb, obscombin, allbias, satdata, gobs1, gobs2);
	//	}
	//}

	//double t_gprecisemodel::ApplyRawTripleDcb(const string& dcb, t_gallbias* allbias, t_gsatdata& satdata, t_gobs* gobs)
	//{
	//	double dcb_value = 0.0;
	//	double intra_dcb = 0.0;
	//	double inter_dcb_12 = 0.0;
	//	double inter_dcb_13 = 0.0;
	//	double inter_dcb_23 = 0.0;

	//	auto   gsys = _crt_obs.gsys();
	//	auto   freq = _freq_index[gsys][gobs->band()];

	//	// jdhuang :
	//	// add intra freq dcb
	//	intra_dcb = intraFreqDcb(dcb, _gallbias, _crt_obs, gobs);

	//	if (freq == FREQ_1 && _frequency == 3)
	//	{
	//		t_gobs gobs1 = _crt_obs.select_range(_band_index[gsys][FREQ_1], true);
	//		t_gobs gobs2 = _crt_obs.select_range(_band_index[gsys][FREQ_2], true);
	//		inter_dcb_12 = interFreqDcb(dcb, _gallbias, _crt_obs, &gobs1, &gobs2);
	//	}
	//	else if (freq == FREQ_2 && _frequency == 3)
	//	{
	//		t_gobs gobs1 = _crt_obs.select_range(_band_index[gsys][FREQ_1], true);
	//		t_gobs gobs2 = _crt_obs.select_range(_band_index[gsys][FREQ_2], true);
	//		inter_dcb_12 = interFreqDcb(dcb, _gallbias, _crt_obs, &gobs1, &gobs2);
	//	}
	//	else if (freq == FREQ_3 && _frequency == 3)
	//	{
	//		t_gobs gobs1 = _crt_obs.select_range(_band_index[gsys][FREQ_1], true);
	//		t_gobs gobs2 = _crt_obs.select_range(_band_index[gsys][FREQ_2], true);
	//		t_gobs gobs3 = _crt_obs.select_range(_band_index[gsys][FREQ_3], true);
	//		inter_dcb_12 = interFreqDcb(dcb, _gallbias, _crt_obs, &gobs1, &gobs2);
	//		inter_dcb_13 = interFreqDcb(dcb, _gallbias, _crt_obs, &gobs1, &gobs3);
	//		inter_dcb_23 = interFreqDcb(dcb, _gallbias, _crt_obs, &gobs2, &gobs3);
	//	}
	//	else
	//	{
	//		inter_dcb_12 = 0.0;
	//		inter_dcb_13 = 0.0;
	//		inter_dcb_23 = 0.0;
	//	}

	//	double alfa = 0.0;
	//	double beta = 0.0;
	//	satdata.coef_ionofree(_band_index[gsys][FREQ_1], alfa, _band_index[gsys][FREQ_2], beta);

	//	switch (freq)
	//	{
	//	case FREQ_1:
	//		dcb_value = intra_dcb + beta * inter_dcb_12;
	//		break;
	//	case FREQ_2:
	//		dcb_value = intra_dcb - alfa * inter_dcb_12;
	//		break;
	//	case FREQ_3:
	//		dcb_value = intra_dcb - alfa * inter_dcb_13 - beta * inter_dcb_23;
	//		break;
	//	default:
	//		dcb_value = 0.0;
	//		_log->logError("t_gprecisemodel", "ApplyTripleDcb", "wrong type of GNSS FREQ, only support FEQ 1 2 3.");
	//		break;
	//	}

	//	return dcb_value;
	//}

	//double t_gprecisemodel::ApplyCmbTripleDcb(const string& dcb, OBSCOMBIN& obscombin, t_gallbias* allbias, t_gsatdata& satdata, t_gobs* gobs1, t_gobs* gobs2)
	//{
	//	double dcb_value = 0.0;
	//	double intra_dcb_1 = 0.0;
	//	double intra_dcb_2 = 0.0;
	//	double inter_dcb_12 = 0.0;
	//	double inter_dcb_13 = 0.0;

	//	auto   gsys = _crt_obs.gsys();
	//	auto   freq_1 = _freq_index[gsys][gobs1->band()];
	//	auto   freq_2 = _freq_index[gsys][gobs2->band()];

	//	// jdhuang :
	//	// add intra freq dcb
	//	intra_dcb_1 = intraFreqDcb(dcb, _gallbias, _crt_obs, gobs1);
	//	intra_dcb_2 = intraFreqDcb(dcb, _gallbias, _crt_obs, gobs2);

	//	//
	//	if (obscombin == MW_COMBIN)
	//	{
	//		// modified for correct P1 code bias glfeng
	//		double l1 = satdata.wavelength(gobs1->band());
	//		double l2 = satdata.wavelength(gobs2->band());

	//		if (double_eq(l2, 0.0))
	//		{
	//			_log->logError("t_gprecisemodel", "ApplyCmbTripleDcb", "wrong type of GNSS FREQ, only support FEQ 1 2 3.");
	//			return 0.0;
	//		}
	//		double fact = l1 / l2;
	//		dcb_value = -(intra_dcb_1 / l1 + intra_dcb_2 / l2) * (1.0 - fact) / (1.0 + fact);

	//		return dcb_value;
	//	}


	//	if (freq_1 == FREQ_1 && freq_2 == FREQ_2 && _frequency == 3)
	//	{
	//		inter_dcb_12 = interFreqDcb(dcb, _gallbias, _crt_obs, gobs1, gobs2);
	//	}
	//	else if (freq_1 == FREQ_1 && freq_2 == FREQ_3 && _frequency == 3)
	//	{
	//		t_gobs gobs1 = _crt_obs.select_range(_band_index[gsys][FREQ_1], false);
	//		t_gobs gobs2 = _crt_obs.select_range(_band_index[gsys][FREQ_2], false);
	//		t_gobs gobs3 = _crt_obs.select_range(_band_index[gsys][FREQ_3], false);
	//		inter_dcb_12 = interFreqDcb(dcb, _gallbias, _crt_obs, &gobs1, &gobs2);
	//		inter_dcb_13 = interFreqDcb(dcb, _gallbias, _crt_obs, &gobs1, &gobs3);
	//	}
	//	else
	//	{
	//		inter_dcb_12 = 0.0;
	//		inter_dcb_13 = 0.0;
	//	}

	//	double alfa12 = 0.0;
	//	double alfa13 = 0.0;
	//	double beta12 = 0.0;
	//	double beta13 = 0.0;
	//	satdata.coef_ionofree(_band_index[gsys][FREQ_1], alfa12, _band_index[gsys][FREQ_2], beta12);
	//	satdata.coef_ionofree(_band_index[gsys][FREQ_1], alfa13, _band_index[gsys][FREQ_3], beta13);

	//	if (freq_1 == FREQ_1 && freq_2 == FREQ_2)
	//	{
	//		dcb_value = intra_dcb_1 * alfa12 + intra_dcb_2 * beta12;
	//	}
	//	else if (freq_1 == FREQ_1 && freq_2 == FREQ_3)
	//	{
	//		dcb_value = intra_dcb_1 * alfa13 + intra_dcb_2 * beta13;
	//		dcb_value = dcb_value + beta12 * inter_dcb_12 - beta13 * inter_dcb_13;
	//	}
	//	else
	//	{
	//		dcb_value = 0.0;
	//		_log->logError("t_gprecisemodel", "ApplyCmbTripleDcb", "wrong type of GNSS FREQ, only support FEQ 1 2 3.");
	//	}

	//	return dcb_value;
	//}

	//void t_gprecisemodel::addCycleSlip(shared_ptr<t_gcycleslip> cycleslip_12, shared_ptr<t_gcycleslip> cycleslip_13)
	//{
	//	_cycleslip_12 = cycleslip_12;
	//	_cycleslip_13 = cycleslip_13;
	//}

	void t_gprecisemodel::_update_rot_matrix(const t_gtime& epoch)
	{
		//lvhb modified in 20200730
		if (_gdata_erp->isEmpty()) return;
		if (!_gdata_erp)
		{
			if (_log)
				_log->comment(1, "t_gorbmodel::_update_rotmatrix ", "have no poleut1 data,cant' compute trs2crs!!");
			throw exception();
		}

		t_gtime tdt = epoch;
		tdt.tsys(t_gtime::TT);

		auto find_iter = _trs2crs_list.find(tdt);
		if (find_iter == _trs2crs_list.end())
		{
			_trs2crs_2000 = t_gtrs2crs(false, _gdata_erp);
			_trs2crs_2000.calcRotMat(tdt, true, true, true);
			_trs2crs_list.insert(make_pair(tdt, _trs2crs_2000));

			auto before_iter = _trs2crs_list.lower_bound(tdt - 300.0);
			if (before_iter != _trs2crs_list.begin())
			{
				_trs2crs_list.erase(_trs2crs_list.begin(), --before_iter);
			}
		}
		else {
			_trs2crs_2000 = find_iter->second;
		}

#ifdef DEBUG_TRS
		//jdhuang : check for the trs2crs_2000
		cout << "Mark = trs2crs  : "
			<< setw(25) << tdt.mjd()
			<< setw(25) << tdt.sod()
			<< setw(25) << tdt.dsec() << endl
			<< setw(25) << fixed << setprecision(15)
			<< trs2crs_2000->rotmat() << endl
			<< trs2crs_2000->rotdx() << endl
			<< trs2crs_2000->rotdy() << endl
			<< trs2crs_2000->rotdu() << endl;
#endif // DEBUG
		return;

	}

	bool t_gprecisemodel::_apply_rec(const t_gtime& crt_epo, const t_gtime& rec_epo, t_gallpar& pars)
	{
		// TODO : improve 
		if (_rec_obj_flag == make_pair(_crt_rec, rec_epo))
		{
			return true;
		}
		else
		{
			_rec_obj_flag = make_pair(_crt_rec, rec_epo);
		}

		int ix = pars.getParam(_crt_rec, t_gpar::CRD_X, "");
		int iy = pars.getParam(_crt_rec, t_gpar::CRD_Y, "");
		int iz = pars.getParam(_crt_rec, t_gpar::CRD_Z, "");

		t_gtriple trs_rec_xyz(0.0, 0.0, 0.0);
		if (ix >= 0 && iy >= 0 && iz >= 0)
		{
			trs_rec_xyz = t_gtriple(pars.getParValue(ix), pars.getParValue(iy), pars.getParValue(iz));
			if (trs_rec_xyz.zero())
			{
				trs_rec_xyz = _crt_obj->crd(rec_epo);
				pars.setParValue(ix, trs_rec_xyz[0]);
				pars.setParValue(iy, trs_rec_xyz[1]);
				pars.setParValue(iz, trs_rec_xyz[2]);
			}
		}
		else
		{
			if (_crd_est != CONSTRPAR::FIX)
			{
				write_log_info(_log, 1, "NOTE", "can't get the crd of" + _site + " in time" + rec_epo.str_ymdhms());
				return false;
			}
			trs_rec_xyz = _crt_obj->crd(rec_epo);
		}

		if (trs_rec_xyz.zero())
		{
			write_log_info(_log, 1, "NOTE", "can't get the crd of" + _site + " in time" + rec_epo.str_ymdhms());
			return false;
		}

		// JDHUANG : 备注了，因为不知道为什么崩溃
		bool tide_valid = _apply_rec_tides(rec_epo, trs_rec_xyz);
		if (!tide_valid)
		{
			write_log_info(_log, 1, "NOTE", "apply tide failed for " + _site + " in time" + rec_epo.str_ymdhms());
			return false;
		}

		//ARP Correction: modified by lvhb in 20200916
		//bool emptyBase = dynamic_cast<t_gsetgen*>(_settings)->list_base().empty();
		string strEst = dynamic_cast<t_gsetgen*>(_settings)->estimator();
		bool isFLT = (/*!emptyBase ||*/ strEst == "FLT");
		// ARP Correction
		if(!isFLT || _crd_est == CONSTRPAR::FIX) trs_rec_xyz = trs_rec_xyz + _crt_obj->eccxyz(rec_epo);
		_trs_rec_crd = trs_rec_xyz;

		// TRS2CRS 
		_update_rot_matrix(rec_epo);
		Matrix trs2crs = _trs2crs_2000.getRotMat();
		ColumnVector crs_rec_xyz = trs2crs * trs_rec_xyz.crd_cvect();
		_crs_rec_crd = t_gtriple(crs_rec_xyz);

		// jdhuang 
		// get vel
		Matrix dtrs2crs = _trs2crs_2000.getMatDu() * OMGE_DOT;
		ColumnVector trs_rec_vel = dtrs2crs * trs_rec_xyz.crd_cvect();
		_crs_rec_vel = t_gtriple(trs_rec_vel);

		if (_gall_orb)
		{
			// set rot  partial matrix in obsdata
			ColumnVector xyz_trs = _crt_obj->crd(rec_epo).crd_cvect();
			_rec_rotmat = _trs2crs_2000.getRotMat();
			_dr_dxpole = _trs2crs_2000.getMatDx() * xyz_trs;
			_dr_dypole = _trs2crs_2000.getMatDy() * xyz_trs;
			_dr_dut1 = _trs2crs_2000.getMatDu() * xyz_trs;
		}
		return true;

	}

	bool t_gprecisemodel::_apply_leo(const t_gtime & crt_epo, const t_gtime & rec_epo, t_gallpar& pars)
	{
		// TODO : improve 
		if (_rec_obj_flag == make_pair(_crt_rec, rec_epo))
		{
			return true;
		}
		else
		{
			_rec_obj_flag = make_pair(_crt_rec, rec_epo);
		}
		t_gtriple xyz_leo(0.0, 0.0, 0.0);
		if (_gall_orb->get_crs_sat_crd(rec_epo, strs2strl(_site), xyz_leo[0], xyz_leo[1], xyz_leo[2])) {

			_crs_rec_crd = xyz_leo;
			_update_rot_matrix(rec_epo);
			_trs_rec_crd = t_gtriple(_trs2crs_2000.getRotMat().t()*xyz_leo.crd_cvect());

			// get crs vel
			t_gtriple vel_leo(0.0, 0.0, 0.0);
			if (!_gall_orb->get_crs_sat_vel(rec_epo, strs2strl(_site), vel_leo[0], vel_leo[1], vel_leo[2]))
			{
				return false;
			}
			_crs_rec_vel = vel_leo;

			vector<t_gpar> orbpars = pars.getOrbPar(strs2strl(_site));
			Matrix partial(0.0, orbpars.size(), 3);
			// TODO
			_gall_orb->getLeoPartial(rec_epo, strs2strl(_crt_rec), partial, orbpars);
			_rec_partial = partial;

			// set sat parameter PXSAT in lsq
			int site_index = pars.getParam(_crt_rec, t_gpar::PXSAT, "");
			//obsdata.addsatindex(par[sat_index].index);
			if (site_index >= 0)
			{
				_rec_index = pars.getParIndex(site_index);
			}

		}
		else
		{
			// geomatric leo estimate

			// use update par value
			int _idx, _idy, _idz;

			//string _name_short = strl2strs((_allkin->get_sat_name()).substr(0, 15));
			_idx = pars.getParam(_crt_rec, t_gpar::CRD_X, "", _crt_epo, _crt_epo);
			_idy = pars.getParam(_crt_rec, t_gpar::CRD_Y, "", _crt_epo, _crt_epo);
			_idz = pars.getParam(_crt_rec, t_gpar::CRD_Z, "", _crt_epo, _crt_epo);

			if (_idx == -1 || _idy == -1 || _idz == -1) {
				write_log_info(_log, 1,"Can't find the Crd par in LSQ pars in LEO geometrc!!");
				return false;
			}

			xyz_leo = t_gtriple(pars[_idx].value(), pars[_idy].value(), pars[_idz].value());


			// get trs&&crs crd
			_trs_rec_crd = xyz_leo;

			// TRS2CRS 
			_update_rot_matrix(rec_epo);
			Matrix trs2crs = _trs2crs_2000.getRotMat();
			ColumnVector crs_rec_xyz = trs2crs * _trs_rec_crd.crd_cvect();
			_crs_rec_crd = t_gtriple(crs_rec_xyz);

			// get crs vel
			t_gtriple xyz_r_next_time;
			_gall_nav->lagrange_pos(strs2strl(_crt_rec), rec_epo + 0.001, xyz_r_next_time);
			//xyz_r_next_time[0] = _allkin->get_value(rec_epo + 0.001, "PX");
			//xyz_r_next_time[1] = _allkin->get_value(rec_epo + 0.001, "PY");
			//xyz_r_next_time[2] = _allkin->get_value(rec_epo + 0.001, "PZ");
			xyz_r_next_time = (xyz_r_next_time - _trs_rec_crd);
			ColumnVector trs_rec_vel = xyz_r_next_time.crd_cvect();
			_crs_rec_vel = t_gtriple(trs2crs * trs_rec_vel);

		}
		
		// update the position of LEO(TRS)
		t_gtriple pos_std(1.0, 1.0, 1.0);
		_grec->crd(_trs_rec_crd, pos_std, rec_epo, LAST_TIME);

		return true;
	}

	bool t_gprecisemodel::_apply_sat(const t_gtime& rec_epo, t_gtime& sat_epo)
	{
		// ITERATION
		// compute sat coord(CRS)  clk(estimated) (rewrite in orb model)
		double delay = 0.0;
		while (true)
		{
			sat_epo = rec_epo - delay;
			if (_realtime && !_ultrasp3)//add by xiongyun for realtime
			{
				bool sat_posvel_valid = _get_crs_sat_crdvel(sat_epo, _crt_sat, _crs_sat_crd, _crs_sat_vel);
				if (!sat_posvel_valid)
				{
					write_log_info(_log, 1, "NOTE", "can not get sat pos or vel for " + _crt_sat);
					return false;
				}
				_update_rot_matrix(rec_epo);
				_trs_sat_crd = t_gtriple(_trs2crs_2000.getRotMat().t() * _crs_sat_crd.crd_cvect());
			}
			else
			{
				// Get CRS
				bool sat_pos_valid = _get_crs_sat_crd(sat_epo, _crt_sat, _crs_sat_crd);
				if (!sat_pos_valid)
				{
					write_log_info(_log, 1, "NOTE", "can not get sat pos for " + _crt_sat);
					return false;
				}
				if (double_eq(_crs_sat_crd[0] * _crs_sat_crd[1] * _crs_sat_crd[2], 0.0) || abs(_crs_sat_crd[0]) >= 1E18)
				{
					write_log_info(_log, 1, "NOTE", "can not get sat pos for " + _crt_sat);
					return false;
				}
				// SET TRS in epoch TR [include earth rotation]
				_update_rot_matrix(rec_epo);
				_trs_sat_crd = t_gtriple(_trs2crs_2000.getRotMat().t() * _crs_sat_crd.crd_cvect());

				bool sat_vel_valid = _get_crs_sat_vel(sat_epo, _crt_sat, _crs_sat_vel);
				if (!sat_vel_valid)
				{
					write_log_info(_log, 1, "NOTE", "can not get sat vel for " + _crt_sat);
					return false;
				}
			}

			// PCO corr sat
			double pco_R = 0.0, pco_S = 0.0;
			if (_gallobj != 0)
			{
				_update_rot_matrix(rec_epo);
				_crt_obs.addcrd(_trs_sat_crd);
				ColumnVector x_earth = (_crs_sat_crd.crd_cvect().t()*_trs2crs_2000.getMatDu() / RAD2TSEC).t();
				_crt_obs.addvel(t_gtriple(_trs2crs_2000.getRotMat().t()*_crs_sat_vel.crd_cvect() - x_earth));
				shared_ptr<t_gobj>  sat_obj = _gallobj->obj(_crt_sat);
				shared_ptr<t_gobj>  rec_obj = _gallobj->obj(_crt_rec);

				shared_ptr<t_gpcv>  sat_pcv;
				shared_ptr<t_gpcv>  rec_pcv;

				if (sat_obj != 0)  sat_pcv = sat_obj->pcv(_crt_epo);
				if (rec_obj != 0)  rec_pcv = rec_obj->pcv(_crt_epo);

				GOBS_LC lc = LC_L1;
				if (sat_pcv && (!_realtime || _ultrasp3))//add _realtime by xiongyun(do not correct pco in realtime when use bnccor)
				{
					// Satellite phase center offset
					t_gtriple pco(0, 0, 0);

					if (sat_pcv->pcoS(_crt_obs, pco, lc, _band_index[_crt_sys][FREQ_1], _band_index[_crt_sys][FREQ_2]) > 0)
					{
						Matrix _rot_matrix = _RotMatrix_Ant(_crt_obs, sat_epo, sat_obj, false);
						t_gtriple dx(_rot_matrix*(pco.crd_cvect()));
						sat_pcv->pco_proj(pco_S, _crt_obs, _trs_rec_crd, dx);
					}

				}

				if (rec_pcv)
				{
					// Receiver phase center offset
					t_gtriple pco(0.0, 0.0, 0.0);
					if (rec_pcv->pcoR(_crt_obs, pco, lc, _band_index[_crt_sys][FREQ_1], _band_index[_crt_sys][FREQ_2]) > 0)
					{
						Matrix _rot_matrix = _RotMatrix_Ant(_crt_obs, rec_epo, rec_obj, false);
						t_gtriple dx(_rot_matrix*(pco.crd_cvect()));
						rec_pcv->pco_proj(pco_R, _crt_obs, _trs_rec_crd, dx);
					}
					pco_R *= -1;
				}
			}

			double delay_temp;
			delay_temp = (_crs_sat_crd - _crs_rec_crd).norm();
			delay_temp += pco_R + pco_S;

			// jdhuang
			if (abs(delay_temp / CLIGHT - delay) < 1E-9) break;
			delay = delay_temp / CLIGHT;
		}

		if (_gall_orb)
		{
			// set funct
			Matrix partial;
			bool part_valid = _gall_orb->getGnsPartial(sat_epo, _crt_sat, partial);
			if (!part_valid)
			{
				write_log_info(_log, 1, "NOTE", "orb file have sat:" + _crt_sat);
				return false;
			}
			_sat_partial = partial;

		}


		// add sat crd vel 
		// Notice this sat crd should be in the TRS
		_crt_obs.addcrd(_trs_sat_crd);
		// modefied by zhshen
		ColumnVector x_earth = (_crs_sat_crd.crd_cvect().t()*_trs2crs_2000.getMatDu() / RAD2TSEC).t();
		_crt_obs.addvel(t_gtriple(_trs2crs_2000.getRotMat().t()*_crs_sat_vel.crd_cvect() - x_earth));
		return true;
	}


	bool t_gprecisemodel::_get_crs_sat_crd(const t_gtime & sat_epoch, const string & sat, t_gtriple & crs_sat_crd)
	{
		if (_gall_orb && _gall_nav)
		{
			write_log_info(_log,2, "NOTE", "there are orb and nav, use orb first");
		}

		if (_gall_orb)	return _get_crs_sat_crd(sat_epoch, sat, _gall_orb, crs_sat_crd);
		if (_gall_nav)	return _get_crs_sat_crd(sat_epoch, sat, _gall_nav, crs_sat_crd);
		write_log_info(_log, 2, "NOTE", "there are no orb and nav");
		return false;
	}

	bool t_gprecisemodel::_get_crs_sat_crd(const t_gtime & sat_epoch, const string & sat, t_gallorb * orb, t_gtriple & crs_sat_crd)
	{
		bool pos_valid = false;
		if (orb)
		{
			pos_valid = orb->get_crs_sat_crd(sat_epoch, sat, crs_sat_crd[0], crs_sat_crd[1], crs_sat_crd[2]);
		}
		if (pos_valid) return pos_valid;
		else write_log_info(_log, 1, "NOTE", "can not get sat crd in epoch" + sat_epoch.str_ymdhms());
		return false;
	}

	bool t_gprecisemodel::_get_crs_sat_crd(const t_gtime & sat_epoch, const string & sat, t_gallnav * nav, t_gtriple & crs_sat_crd)
	{
		bool pos_valid = false;
		if (nav)
		{
			double xyz_sat[3];
			pos_valid = (nav->pos(sat, sat_epoch, xyz_sat) >= 0) ? true : false;

			// TRS2CRS
			_update_rot_matrix(sat_epoch);
			crs_sat_crd = t_gtriple(xyz_sat);
			crs_sat_crd = t_gtriple(_trs2crs_2000.getRotMat()*crs_sat_crd.crd_cvect());
		}

		if (pos_valid) return pos_valid;
		else write_log_info(_log, 2, "NOTE", sat_epoch.str_ymdhms("can not get sat crd in epoch",false));
		return false;
	}

	bool t_gprecisemodel::_get_crs_sat_crdvel(const t_gtime& sat_epoch, const string& sat, t_gtriple& crs_sat_crd, t_gtriple& crs_sat_vel)
	{
		if (!_gall_nav)
		{
			write_log_info(_log, 2, "NOTE", "there are no nav");
			return false;
		}

		bool pos_valid = false;
		bool cor_valid = false;
		double xyz_corr[3] = { 0.0, 0.0, 0.0 };
		double vel[3] = { 0.0, 0.0, 0.0 };
		double vel_corr[3] = { 0.0, 0.0, 0.0 };
		double var[3] = { 0.0, 0.0, 0.0 };
		double clk = 0.0, dclk = 0.0;
		double clk_corr = 0.0;
		double dclk_corr = 0.0;

		double xyz_sat[3] = { 0.0, 0.0, 0.0 };
		pos_valid = (_gall_nav->pos(sat, sat_epoch, xyz_sat, var, vel, true) >= 0) ? true : false;

		// get the eph iod issue
		int iod = _gall_nav->get_iod(sat, sat_epoch);
		if (iod < 0) return false;

		cor_valid = _gall_nav->get_pos_clk_correction(sat, sat_epoch, iod, xyz_corr, vel_corr, clk_corr, dclk_corr);
		if (!cor_valid) pos_valid = false;
		_correction_orb(xyz_sat, vel, xyz_corr, vel_corr);

		// TRS2CRS
		_update_rot_matrix(sat_epoch);
		crs_sat_crd = t_gtriple(xyz_sat);
		crs_sat_crd = t_gtriple(_trs2crs_2000.getRotMat() * crs_sat_crd.crd_cvect());
		crs_sat_vel = t_gtriple(vel);
		crs_sat_vel = t_gtriple(_trs2crs_2000.getRotMat() * crs_sat_vel.crd_cvect());

		if (pos_valid) return pos_valid;
		else write_log_info(_log, 2, "NOTE", sat_epoch.str_ymdhms("can not get sat crd or vel in epoch", false));
		return false;
	}

	bool t_gprecisemodel::_get_sat_clk_corr(const t_gtime & sat_epoch, const string & sat, double& clk, double& dclk)
	{
		if (!_gall_nav)
		{
			write_log_info(_log, 2, "NOTE", "there are no nav");
			return false;
		}

		bool clk_valid = true;
		int cor_valid = false;
		double xyz_corr[3] = { 0.0, 0.0, 0.0 };
		double vel[3] = { 0.0, 0.0, 0.0 };
		double vel_corr[3] = { 0.0, 0.0, 0.0 };
		double var[3] = { 0.0, 0.0, 0.0 };
		double clk_corr = 0.0;
		double dclk_corr = 0.0;

		// get the eph iod issue
		int iod = _gall_nav->get_iod(sat, sat_epoch);
		if (iod < 0) return false;

		cor_valid = _gall_nav->get_pos_clk_correction(sat, sat_epoch, iod, xyz_corr, vel_corr, clk_corr, dclk_corr);
		if (cor_valid < 0 || double_eq(clk_corr, 0.0)) clk_valid = false;
		_correction_clk(clk, dclk, clk_corr, dclk_corr);

		if (clk_valid) return clk_valid;
		else write_log_info(_log, 2, "NOTE", sat_epoch.str_ymdhms("can not get sat clk in epoch", false));
		return false;
	}

	int t_gprecisemodel::_correction_orb(double* xyz, double* vxyz, double* xyz_corr, double* vxyz_corr)
	{
		ColumnVector pos(xyz, 3), vel(vxyz, 3), rao_pos(xyz_corr, 3), rao_vel(vxyz_corr, 3), out(0.0, 3);
		rao2xyz(pos, vel, rao_pos, out);
		for (size_t i = 0; i < 3; i++)
		{
			xyz[i] -= out(i + 1);
		}
		rao2xyz(pos, vel, rao_vel, out);
		for (size_t i = 0; i < 3; i++)
		{
			vxyz[i] -= out(i + 1);
		}
		return 1;
	}

	int t_gprecisemodel::_correction_clk(double & clk, double & dclk, double & clk_corr, double & dclk_corr)
	{
		clk += clk_corr / CLIGHT;
		dclk += dclk_corr / CLIGHT;
		return 1;
	}

	bool t_gprecisemodel::_get_crs_sat_vel(const t_gtime& sat_epoch, string sat, t_gtriple& sat_vel)
	{
		if (!_gall_orb && !_gall_nav)
		{
			write_log_info(_log, 1, "NOTE", "cannot get vel " + sat);
			return false;
		}

		if ((_gall_orb && _gall_nav) ||
			(_gall_orb && !_gall_nav))
		{
			return _gall_orb->get_crs_sat_vel(sat_epoch, sat, sat_vel[0], sat_vel[1], sat_vel[2]);
		}

		if (_gall_nav && !_gall_orb)
		{
			// get vel in CRS
			t_gtriple satcrd_next(0.0, 0.0, 0.0);
			bool crd_valid = _get_crs_sat_crd(sat_epoch + 0.00001, sat, _gall_nav, satcrd_next);
			if (!crd_valid) return false;
			sat_vel = (satcrd_next - _crs_sat_crd) * 100000;
			return true;
		}

		return true;
	}

	// jdhuang : fix
	// use _freq_index replace band_priority
	//double t_gprecisemodel::_Partial(const t_gtime& epoch, t_gsatdata& obsdata, const t_gobs & gobs1, const t_gobs& gobs2, t_gpar& par)
	double t_gprecisemodel::_Partial(const t_gtime& epoch, t_gsatdata& obsdata, const t_gobs & gobs, t_gpar& par)
	{
		double mfw, dmfw, mfh, dmfh;
		mfw = dmfw = mfh = dmfh = 0.0;
		t_gtriple ell(0.0, 0.0, 0.0);

		switch (par.parType)
		{
		case t_gpar::CRD_X: if (obsdata.site() == par.site)
		{
			return (par.value() - obsdata.satcrd().crd(0)) / obsdata.rho();
		}
							else return 0.0;
		case t_gpar::CRD_Y: if (obsdata.site() == par.site)
		{
			return (par.value() - obsdata.satcrd().crd(1)) / obsdata.rho();
		}
							else return 0.0;
		case t_gpar::CRD_Z: if (obsdata.site() == par.site)
		{
			return (par.value() - obsdata.satcrd().crd(2)) / obsdata.rho();
		}
							else return 0.0;
		case t_gpar::CLK:
		{
			if (obsdata.site() == par.site) return 1.0 - obsdata.drate();
			else return 0.0;
		}
		case t_gpar::CLK_G:
		{
			if (obsdata.site() == par.site && obsdata.gsys() == GSYS::GPS) return 1.0 - obsdata.drate();
			else return 0.0;
		}
		case t_gpar::CLK_R:
		{
			if (obsdata.site() == par.site && obsdata.gsys() == GSYS::GLO) return 1.0 - obsdata.drate();
			else return 0.0;
		}
		case t_gpar::CLK_C:
		{
			if (obsdata.site() == par.site && obsdata.gsys() == GSYS::BDS) return 1.0 - obsdata.drate();
			else return 0.0;
		}
		case t_gpar::CLK_E:
		{
			if (obsdata.site() == par.site && obsdata.gsys() == GSYS::GAL) return 1.0 - obsdata.drate();
			else return 0.0;
		}
		case t_gpar::CLK_J:
		{
			if (obsdata.site() == par.site && obsdata.gsys() == GSYS::QZS) return 1.0 - obsdata.drate();
			else return 0.0;
		}
		case t_gpar::IFCB_F3: if (t_gsys::band2freq(obsdata.gsys(), gobs.band()) == FREQ_3
			&& par.prn == obsdata.sat() && _isTriple) return -1.0; else return 0.0;
		case t_gpar::IFCB_F4: if (t_gsys::band2freq(obsdata.gsys(), gobs.band()) == FREQ_4
			&& par.prn == obsdata.sat()) return -1.0; else return 0.0;
		case t_gpar::IFCB_F5: if (t_gsys::band2freq(obsdata.gsys(), gobs.band()) == FREQ_5
			&& par.prn == obsdata.sat()) return -1.0; else return 0.0;
		case t_gpar::IFB_C3: if (_freq_index[obsdata.gsys()][gobs.band()] == FREQ_3 && par.prn == obsdata.sat() && par.site == obsdata.site()) return 1.0; else return 0.0;		
		case t_gpar::IFB_C4: if (_freq_index[obsdata.gsys()][gobs.band()] == FREQ_4 && par.prn == obsdata.sat() && par.site == obsdata.site()) return 1.0; else return 0.0;
		case t_gpar::IFB_C5: if (_freq_index[obsdata.gsys()][gobs.band()] == FREQ_5 && par.prn == obsdata.sat() && par.site == obsdata.site()) return 1.0; else return 0.0;
		case t_gpar::SAT_IFB_C3: if (_freq_index[obsdata.gsys()][gobs.band()] == FREQ_3 && par.prn == obsdata.sat()) return 1.0; else return 0.0;
		case t_gpar::SAT_IFB_C4: if (_freq_index[obsdata.gsys()][gobs.band()] == FREQ_4 && par.prn == obsdata.sat()) return 1.0; else return 0.0;
		case t_gpar::SAT_IFB_C5: if (_freq_index[obsdata.gsys()][gobs.band()] == FREQ_5 && par.prn == obsdata.sat()) return 1.0; else return 0.0;
		case t_gpar::GPS_REC_IFB_C3: if (obsdata.gsys() == GPS && _freq_index[obsdata.gsys()][gobs.band()] == FREQ_3 && par.site == obsdata.site()) return 1.0; else return 0.0;
		case t_gpar::GPS_REC_IFB_C4: if (obsdata.gsys() == GPS && _freq_index[obsdata.gsys()][gobs.band()] == FREQ_4 && par.site == obsdata.site()) return 1.0; else return 0.0;
		case t_gpar::GPS_REC_IFB_C5: if (obsdata.gsys() == GPS && _freq_index[obsdata.gsys()][gobs.band()] == FREQ_5 && par.site == obsdata.site()) return 1.0; else return 0.0;
		
		case t_gpar::GAL_REC_IFB_C3: if (obsdata.gsys() == GAL && _observ == RAW_ALL &&_freq_index[obsdata.gsys()][gobs.band()] == FREQ_3 && par.site == obsdata.site()) return 1.0; else return 0.0;
		case t_gpar::GAL_REC_IFB_C4: if (obsdata.gsys() == GAL && _observ == RAW_ALL &&_freq_index[obsdata.gsys()][gobs.band()] == FREQ_4 && par.site == obsdata.site()) return 1.0; else return 0.0;
		case t_gpar::GAL_REC_IFB_C5: if (obsdata.gsys() == GAL && _observ == RAW_ALL &&_freq_index[obsdata.gsys()][gobs.band()] == FREQ_5 && par.site == obsdata.site()) return 1.0; else return 0.0;

		case t_gpar::CLK_SAT:    if (obsdata.sat() == par.prn)              return -1.0; else return 0.0;
		case t_gpar::G15CLK_SAT: if (obsdata.sat() == par.prn && _isTriple) return -1.0; else return 0.0;
		case t_gpar::E15CLK_SAT: if (obsdata.sat() == par.prn && _isTriple) return -1.0; else return 0.0;
		case t_gpar::C15CLK_SAT: if (obsdata.sat() == par.prn && _isTriple) return -1.0; else return 0.0;
		case t_gpar::J15CLK_SAT: if (obsdata.sat() == par.prn && _isTriple) return -1.0; else return 0.0;
		case t_gpar::TRP:
			if (obsdata.site() == par.site) {
				xyz2ell(_trs_rec_crd, ell, false);
				_getmf(par, obsdata, ell, epoch, mfw, mfh, dmfw, dmfh);
				return mfw;
			}
			else return 0.0;
		case t_gpar::SION:
		{
			if (obsdata.site() == par.site) 
			{
				// jdhuang : form G01_F
				auto gsys = obsdata.gsys();
				double f1 = obsdata.frequency(_band_index[gsys][FREQ_1]);
				double fk = obsdata.frequency(gobs.band());
				double alfa = 0.0;
				if (gobs.is_phase() && par.prn == obsdata.sat()) 
				{ 
					alfa = -(f1*f1) / (fk*fk); 
				}
				if (gobs.is_code() && par.prn == obsdata.sat()) 
				{
					alfa = (f1*f1) / (fk*fk); 
				}
				return alfa;
			}
			else return 0.0;
		}
		case t_gpar::VION:
		{
			if (obsdata.site() == par.site)
			{
				// jdhuang : change from G01_F
				auto gsys = obsdata.gsys();
				double f1 = obsdata.frequency(_band_index[gsys][FREQ_1]);
				double fk = obsdata.frequency(gobs.band());
				double mf = 1.0 / sqrt(1.0 - pow(R_SPHERE / (R_SPHERE + 450000.0) * sin(G_PI / 2.0 - obsdata.ele()), 2));
				double alfa = 0.0;
				if (gobs.is_phase() && par.prn == obsdata.sat()) { alfa = -(f1*f1) / (fk*fk); }
				if (gobs.is_code() && par.prn == obsdata.sat()) { alfa = (f1*f1) / (fk*fk); }
				return alfa * mf;
			}
			else return 0.0;
		}
		case t_gpar::P1P2G_REC:
		{
			if (obsdata.site() == par.site)
			{
				double f1 = G01_F;
				double fk = obsdata.frequency(gobs.band());
				double alfa = (f1*f1) / (fk*fk);
				double beta = (G02_F*G02_F) / (G01_F*G01_F - G02_F * G02_F);
				FREQ_SEQ freq = t_gsys::band2freq(obsdata.gsys(), gobs.band());
				if (obsdata.gsys() == GPS && gobs.is_code() && (freq == FREQ_1 || freq == FREQ_2)) { return -alfa * beta; }
				else return 0.0;
			}
			else return 0.0;
		}
		case t_gpar::P1P2E_REC:
		{
			if (obsdata.site() == par.site) {
				double f1 = E01_F;
				double fk = obsdata.frequency(gobs.band());
				double alfa = (f1*f1) / (fk*fk);
				double beta = (E05_F*E05_F) / (E01_F*E01_F - E05_F * E05_F);
				FREQ_SEQ freq = t_gsys::band2freq(obsdata.gsys(), gobs.band());
				if (obsdata.gsys() == GAL && gobs.is_code() && (freq == FREQ_1 || freq == FREQ_2)) { return -alfa * beta; }
				else return 0.0;
			}
			else return 0.0;
		}
		case t_gpar::GRD_N:
			if (obsdata.site() == par.site) {
				xyz2ell(_trs_rec_crd, ell, false);
				_getmf(par, obsdata, ell, epoch, mfw, mfh, dmfw, dmfh);
				if (_mf_grd == CHEN_HERRING) {
					double sinel = sin(obsdata.ele());
					double tanel = tan(obsdata.ele());
					double cosaz = cos(obsdata.azi());
					return (1.0 / (sinel*tanel + 0.0032)) * cosaz;
				}
				else if (_mf_grd == TILTING) {
					double cosaz = cos(obsdata.azi());
					return dmfw * cosaz;
				}
				else if (_mf_grd == BAR_SEVER) {
					double tanel = tan(obsdata.ele());
					double cosaz = cos(obsdata.azi());
					return mfw * (1.0 / tanel) * cosaz;
				}
				else cerr << "Grad N mapping function is not set up correctly!!!" << endl;
			}
			else return 0.0;
		case t_gpar::GRD_E:
			if (obsdata.site() == par.site) {
				_getmf(par, obsdata, ell, epoch, mfw, mfh, dmfw, dmfh);
				if (_mf_grd == CHEN_HERRING) {
					double sinel = sin(obsdata.ele());
					double tanel = tan(obsdata.ele());
					double sinaz = sin(obsdata.azi());
					return (1.0 / (sinel*tanel + 0.0032)) * sinaz;
				}
				else if (_mf_grd == TILTING) {
					double sinaz = sin(obsdata.azi());
					return dmfw * sinaz;
				}
				else if (_mf_grd == BAR_SEVER) {
					double tanel = tan(obsdata.ele());
					double sinaz = sin(obsdata.azi());
					return mfw * (1 / tanel) * sinaz;
				}
				else cerr << "Grad E mapping function is not set up correctly!!!" << endl;
			}
			else return 0.0;

		case t_gpar::GLO_ISB:
			if (obsdata.site() == par.site && obsdata.gsys() == GLO/*&& gobs.is_code()*/) return 1.0;
			else return 0.0;
		case t_gpar::GLO_IFCB:
			if (obsdata.site() == par.site && !gobs.is_phase() && obsdata.gsys() == GLO && par.prn == obsdata.sat()) return 1.0;
			else return 0.0;
		case t_gpar::GLO_IFPB:
			if (obsdata.site() == par.site && gobs.is_phase() && obsdata.gsys() == GLO && par.prn == obsdata.sat()) return 1.0;
			else return 0.0;
		case t_gpar::GLO_IFB:
			if (obsdata.site() == par.site && obsdata.gsys() == GLO && par.prn == obsdata.sat() && gobs.is_code()) return 1.0;
			else return 0.0;
		case t_gpar::GAL_ISB:
			if (obsdata.site() == par.site && obsdata.gsys() == GAL /*&& gobs.is_code()*/) return 1.0;
			else return 0.0;
		case t_gpar::BDS_ISB:
			if (obsdata.site() == par.site && obsdata.gsys() == BDS /*&& gobs.is_code()*/) return 1.0;
			else return 0.0;
		case t_gpar::QZS_ISB:
			if (obsdata.site() == par.site && obsdata.gsys() == QZS /*&& gobs.is_code()*/) return 1.0;
			else return 0.0;
		case t_gpar::IFB_QZS:
			if (obsdata.site() == par.site && obsdata.gsys() == QZS && _freq_index[obsdata.gsys()][gobs.band()] == FREQ_3 && gobs.is_code()) return 1.0;
			else return 0.0;
		case t_gpar::IFB_GPS:
			if (obsdata.site() == par.site && obsdata.gsys() == GPS && _freq_index[obsdata.gsys()][gobs.band()] == FREQ_3&& gobs.is_code()) return 1.0;
			else return 0.0;
		case t_gpar::IFB_GAL:
			if (obsdata.site() == par.site && obsdata.gsys() == GAL && _freq_index[obsdata.gsys()][gobs.band()] == FREQ_3 && gobs.is_code()) return 1.0;
			else return 0.0;
		case t_gpar::IFB_GAL_2:
			if (obsdata.site() == par.site && obsdata.gsys() == GAL && _freq_index[obsdata.gsys()][gobs.band()] == FREQ_4 && gobs.is_code()) return 1.0;
			else return 0.0;
		case t_gpar::IFB_GAL_3:
			if (obsdata.site() == par.site && obsdata.gsys() == GAL && _freq_index[obsdata.gsys()][gobs.band()] == FREQ_5 && gobs.is_code()) return 1.0;
			else return 0.0;
		case t_gpar::IFB_BDS:
			if (obsdata.site() == par.site && obsdata.gsys() == BDS && _freq_index[obsdata.gsys()][gobs.band()] == FREQ_3 && gobs.is_code()) return 1.0;
			else return 0.0;
		case t_gpar::IFB_BDS_2:
			if (obsdata.site() == par.site && obsdata.gsys() == BDS && _freq_index[obsdata.gsys()][gobs.band()] == FREQ_4 && gobs.is_code()) return 1.0;
			else return 0.0;
		case t_gpar::IFB_BDS_3:
			if (obsdata.site() == par.site && obsdata.gsys() == BDS && _freq_index[obsdata.gsys()][gobs.band()] == FREQ_5 && gobs.is_code()) return 1.0;
			else return 0.0;
		}

		if (par.parType >= t_gpar::XPOLE && par.parType <= t_gpar::DUT1)
		{
			if (_grec->id_type() == t_gdata::REC_LEO) 
			{
				return 0.0;
			}
			return _Partial_Eop(par, obsdata);
		}
		else if (par.parType >= t_gpar::PXSAT && par.parType <= t_gpar::ECOM2_X1s)
		{
			if (obsdata.sat() == par.prn)
			{
				return _Partial_Orb(par, obsdata);
			}
			else if (obsdata.site() == par.site) 
			{
				return _Partial_Orb_Leo(par, obsdata);
			}
			else
			{
				return 0.0;
			}
		}

		return 0.0;
	}


	double t_gprecisemodel::_Partial_Eop(t_gpar& par, t_gsatdata& obsdata)
	{
		switch (par.parType)
		{
		case t_gpar::XPOLE:
			return -_dr_derp(1, 1);
			break;
		case t_gpar::YPOLE:
			return -_dr_derp(1, 3);
			break;
		case t_gpar::DXPOLE:
			return -_dr_derp(1, 2);
			break;
		case t_gpar::DYPOLE:
			return -_dr_derp(1, 4);
			break;
		case t_gpar::UT1:
			return _dr_derp(1, 5);
			break;
		case t_gpar::DUT1:
			return _dr_derp(1, 6);
			break;
		default:
			return 0.0;
			break;
		}
	}


	double t_gprecisemodel::_Partial_Orb(t_gpar& par, t_gsatdata& obsdata)
	{
		double tmp = 0.0;
		try
		{
			tmp = _sat_partial(par.index - _sat_index + 1, 1);
		}
		catch(...)
		{
			throw logic_error("check your ics and force setting in xml.");
		}
		return tmp;
	}

	double t_gprecisemodel::_Partial_Orb_Leo(t_gpar& par, t_gsatdata& obsdata)
	{

		if (par.beg > obsdata.epoch() && par.lPWC)
		{
			return 0.0;
		}

		//Matrix funct = obsdata.orbfunct();
		Matrix funct = _rec_partial;
		int index_begin = _rec_index;

		// get unit_vector in crs
		ColumnVector unit_rec2sat(3);
		unit_rec2sat = _crs_rec_pco.crd_cvect() - _crs_sat_pco.crd_cvect();
		unit_rec2sat = unit_rec2sat / unit_rec2sat.NormFrobenius();

		// order: P X0 Y0 Z0 VX0..VZ0 solarPar
		funct = funct * unit_rec2sat;
		funct = funct * 1E3;

		return funct(par.index - index_begin + 1, 1);

	}


	void t_gprecisemodel::_getmf(t_gpar& par, t_gsatdata& satData, const t_gtriple& crd, const t_gtime& epoch, double& mfw, double& mfh, double& dmfw, double& dmfh)
	{
		if (par.parType != t_gpar::TRP && par.parType != t_gpar::GRD_N && par.parType != t_gpar::GRD_E) return;

		double ele = satData.ele();

		if (_mf_ztd == COSZ) 
		{
			mfw = mfh = 1.0 / sin(ele);
		}
		else if (_mf_ztd == GMF) 
		{
			t_gmf mf;
			mf.gmf(epoch.mjd(), crd[0], crd[1], crd[2], G_PI / 2.0 - ele, mfh, mfw, dmfh, dmfw);
		}
		else cerr << "ZTD mapping function is not set up correctly!!!" << endl;
	}

	Matrix t_gprecisemodel::_RotMatrix_Ant(t_gsatdata& obsdata, const t_gtime& epoch, shared_ptr<t_gobj> obj, bool isCRS)
	{
		t_gdata::ID_TYPE type = obj->id_type();
		Matrix rotmatrix(3, 3);
		t_gtriple ell(0.0, 0.0, 0.0);
		double sinPhi, cosPhi, sinLam, cosLam;
		if (type == t_gdata::TRN) {
			string antenna = (obj->pcv(_crt_epo))->anten();
			t_gtriple dx(0, 0, 0);
			ColumnVector i(3);
			ColumnVector j(3);
			ColumnVector k(3);
			if (_attitudes == YAW_NOMI)	   attitude(obsdata, "", i, j, k);				//nominal modeling
			else if (_attitudes == YAW_RTCM) attitude(obsdata, obsdata.yaw(), i, j, k);	//value from RTCM used
			else							   attitude(obsdata, antenna, i, j, k);		//default
			rotmatrix.Column(1) = i;
			rotmatrix.Column(2) = j;
			rotmatrix.Column(3) = k;
		}
		else if (type == t_gdata::REC) {
			xyz2ell(_trs_rec_crd, ell, false);
			sinPhi = sin(ell[0]);
			cosPhi = cos(ell[0]);
			sinLam = sin(ell[1]);
			cosLam = cos(ell[1]);
			rotmatrix << -sinPhi * cosLam << -sinLam << +cosPhi * cosLam
				<< -sinPhi * sinLam << +cosLam << +cosPhi * sinLam
				<< +cosPhi << 0.0 << +sinPhi;
		}
		else if (type == t_gdata::REC_LEO) {
			// rot scf2trs
			update_rotmat_scf(obsdata.site(), epoch, _crs_rec_crd.crd_cvect(), _crs_rec_vel.crd_cvect());
			rotmatrix.Column(1) = _rot_scf2trs_leo.Column(2);
			rotmatrix.Column(2) = _rot_scf2trs_leo.Column(1);
			rotmatrix.Column(3) = _rot_scf2trs_leo.Column(3);
		}
		else {
			throw exception();
		}

		if (isCRS) {
			// change trs2crs
			_update_rot_matrix(epoch);
			rotmatrix = _trs2crs_2000.getRotMat()*rotmatrix;
		}

		return rotmatrix;
	}

	void t_gprecisemodel::update_rotmat_scf(string leosat, const t_gtime& epoch, const ColumnVector& pos, const ColumnVector& vel)
	{

		string leoname = strs2strl(leosat);
		_update_rot_matrix(epoch);
		Matrix rot_trs2crs = _trs2crs_2000.getRotMat();

		if (_scf2crs_leo.getSatName()=="")
		{
			_scf2crs_leo =t_gscf2crsleo(leoname);
		}

		if (leoname != _scf2crs_leo.getSatName())
		{
			_scf2crs_leo = t_gscf2crsleo(leoname);
		}

		if (epoch != _scf2crs_leo.getCurtEpoch())
		{
			if (_gall_att != nullptr)
			{
				_scf2crs_leo.calcRotMat(epoch, _gall_att->attitude(leoname), pos, vel, rot_trs2crs, _rot_scf2crs_leo);

			}
			else
			{
				_scf2crs_leo.calcRotMat(epoch, NULL, pos, vel, rot_trs2crs, _rot_scf2crs_leo);

			}
			_rot_scf2trs_leo = rot_trs2crs.t() * _rot_scf2crs_leo;
		}
		return;
	}


	bool t_gprecisemodel::_apply_rec_tides(const t_gtime& epoch, t_gtriple& rec)
	{
		if (!_tide)
		{
			_log->comment(0, "ERROR", "The tide ptr is nullptr.");
			return false;
		}

		_update_rot_matrix(epoch);
		double xpole = _trs2crs_2000.getXpole();
		double ypole = _trs2crs_2000.getYpole();
		double gast = _trs2crs_2000.getGmst();
		Matrix rot_trs2crs = _trs2crs_2000.getRotMat();

		t_gtriple tide(0.0, 0.0, 0.0);
		try
		{
			// solid tide
			t_gtideIERS* solid_ptr = dynamic_cast<t_gtideIERS*>(_tide.get());
			if (!solid_ptr)
			{
				_log->comment(0, "ERROR", "can not get the solid tide ptr.");
				return false;
			}
			t_gtriple solid_earth = solid_ptr->tide_solid(epoch, rec, rot_trs2crs, _gdata_navde);

			// ocean load
			t_gtriple load_ocean = _tide->load_ocean(epoch, _site, rec);

			//pole tide
			t_gtideIERS* pole_ptr = dynamic_cast<t_gtideIERS*>(_tide.get());
			if (!pole_ptr)
			{
				_log->comment(0, "ERROR", "can not get the pole tide ptr.");
				return false;
			}
			t_gtriple tide_pole = pole_ptr->tide_pole_pod(epoch, xpole, ypole, rec);

			//atmosph
			t_gtriple load_atmosph = _tide->load_atmosph();

			// freq tide
			t_gtideIERS* freq_ptr = dynamic_cast<t_gtideIERS*>(_tide.get());
			if (!freq_ptr)
			{
				_log->comment(0, "ERROR", "can not get the freq tide ptr.");
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
				<< " " << epoch.mjd()
				<< " " << epoch.sod()
				<< " " << tide[0]
				<< " " << tide[1]
				<< " " << tide[2] << endl;

			cout << "solid_earth" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd()
				<< " " << epoch.sod()
				<< " " << solid_earth[0]
				<< " " << solid_earth[1]
				<< " " << solid_earth[2] << endl;

			cout << "load_ocean" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd()
				<< " " << epoch.sod()
				<< " " << load_ocean[0]
				<< " " << load_ocean[1]
				<< " " << load_ocean[2] << endl;

			cout << "tide_pole" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd()
				<< " " << epoch.sod()
				<< " " << tide_pole[0]
				<< " " << tide_pole[1]
				<< " " << tide_pole[2] << endl;

			cout << "load_atmosph" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd()
				<< " " << epoch.sod()
				<< " " << load_atmosph[0]
				<< " " << load_atmosph[1]
				<< " " << load_atmosph[2] << endl;

			cout << "tide_freq" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd()
				<< " " << epoch.sod()
				<< " " << tide_freq[0]
				<< " " << tide_freq[1]
				<< " " << tide_freq[2] << endl;
#endif // DEBUG

		}
		catch (...)
		{
			_log->comment(0, "ERROR", "can not get the tide ptr.");
			return false;
		}

		//unit to m
		rec = rec + tide * 1.e3;
		return true;

	}

	bool t_gprecisemodel::_check_obs_valid(t_gobsgnss& obsdata, const OBSCOMBIN& obscom, const GOBSBAND& b1, const GOBSBAND& b2)
	{
		// jdhuang :
		// todo, i am not sure we need to check following ptrs in this function ?
		if (!_gall_nav || !_gallobj || !_gallbias)
		{
			write_log_info(_log, 1, "NOTE", "no nav or obj or bias");
			return false;
		}

		// jdhuang
		// follow is only used for IF model
		// todo, we need to change it to support freq 1 + 3 IF model
		GSYS gs = obsdata.gsys();
		t_gobs gobs1;
		t_gobs gobs2;

		double ObsP = 0.0, ObsL = 0.0;
		switch (obscom)
		{
		case gnut::IONO_FREE:
		{
			gobs1 = t_gobs(obsdata.select_range(b1));
			gobs2 = t_gobs(obsdata.select_range(b2));
			ObsP = obsdata.P3(gobs1, gobs2);
			gobs1 = t_gobs(obsdata.select_phase(b1));
			gobs2 = t_gobs(obsdata.select_phase(b2));
			ObsL = obsdata.L3(gobs1, gobs2);
			if (double_eq(ObsP, 0.0) || double_eq(ObsL, 0.0))
			{
				write_log_info(_log, 1, "NOTE", "no range or phase data");
				return false;
			}
			break;
		}
		case gnut::RAW_SINGLE:
		{
			gobs1 = t_gobs(obsdata.select_range(b1));
			ObsP = obsdata.getobs(gobs1.gobs());
			gobs1 = t_gobs(obsdata.select_phase(b1));
			ObsL = obsdata.getobs(gobs1.gobs());
			if (double_eq(ObsP, 0.0) || double_eq(ObsL, 0.0))
			{
				write_log_info(_log, 1, "NOTE", "no range or phase data");
				return false;
			}
			break;
		}
		case gnut::RAW_DOUBLE:
		{
			gobs1 = t_gobs(obsdata.select_range(b1));
			gobs2 = t_gobs(obsdata.select_phase(b1));
			ObsP = obsdata.getobs(gobs1.gobs());
			ObsL = obsdata.getobs(gobs2.gobs());
			if (double_eq(ObsP, 0.0) || double_eq(ObsL, 0.0))
			{
				write_log_info(_log, 1, "NOTE", "no range or phase data");
				return false;
			}
			gobs1 = t_gobs(obsdata.select_range(b2));
			gobs2 = t_gobs(obsdata.select_phase(b2));
			ObsP = obsdata.getobs(gobs1.gobs());
			ObsL = obsdata.getobs(gobs2.gobs());
			if (double_eq(ObsP, 0.0) || double_eq(ObsL, 0.0))
			{
				write_log_info(_log, 1, "NOTE", "no range or phase data");
				return false;
			}
			break;
		}
		case gnut::RAW_ALL:
		{
			std::map<FREQ_SEQ, GOBSBAND>  sys_bands = _band_index[gs];
			for (const auto& band : sys_bands)
			{
				FREQ_SEQ freqX = band.first;
				GOBSBAND bandX = band.second;
				gobs1 = t_gobs(obsdata.select_range(bandX, true));
				gobs2 = t_gobs(obsdata.select_phase(bandX, true));
				ObsP = 0.0;
				ObsL = 0.0;
				ObsP = obsdata.getobs(gobs1.gobs());
				ObsL = obsdata.getobs(gobs2.gobs());
				if ((double_eq(ObsP, 0.0) || double_eq(ObsL, 0.0)))
				{
					write_log_info(_log, 1, "NOTE", "Missing observations");
					return false;
				}
			}
			break;

			//// jdhuang : fix  ==========================================================================
			//// case 1 : 1+2, 1+3, 
			//// case 2 : 1+2+3, but without 2 or 3
			//// case 3 :
			//std::map<FREQ_SEQ, GOBSBAND>  sys_bands = _band_index[gs];
			//GOBSBAND  band1 = _band_index[gs][FREQ_1];
			//GOBSBAND  band2 = _band_index[gs][FREQ_2];

			//// jdhuang 
			//// add to check dif case
			//bool isDoubleFreq = false;
			//bool isTripleFreq = false;
			//bool isMultipFreq = false;

			//// currunetly, L1 is neccessary, we do not want to process freq 2 + 3 or something else
			//
			//// jdhuang
			//// todo
			//// freq > 3, can not process , tpl return false
			//if (sys_bands.size() > 3) return false;

			//if (sys_bands.size() == 2) isDoubleFreq = true;
			//if (sys_bands.size() == 3) isTripleFreq = true;
			//if (sys_bands.size() > 3)  isMultipFreq = true;		
			//if (_freq_model == FREQ_MODEL::TRIPLE && !isTripleFreq) return false;
			//if (_freq_model == FREQ_MODEL::DOUBLE && !isDoubleFreq) return false;

			//// jdhuang fix:
			//// add more case
			//// for (std::vector<GOBSBAND>::const_iterator iter = iter_beg; iter != iter_end; ++iter)
			//set<GOBSBAND> valid_band;
			//for (const auto& band : sys_bands)
			//{
			//	FREQ_SEQ freqX = band.first;
			//	GOBSBAND bandX = band.second;
			//	gobs1 = t_gobs(obsdata.select_range(bandX, true));
			//	gobs2 = t_gobs(obsdata.select_phase(bandX, true));
			//	ObsP = 0.0;
			//	ObsL = 0.0;
			//	ObsP = obsdata.getobs(gobs1.gobs());
			//	ObsL = obsdata.getobs(gobs2.gobs());

			//	// double freq, they should both have ObsL and ObsP
			//	if (isDoubleFreq && (double_eq(ObsP, 0.0) || double_eq(ObsL, 0.0)))
			//	{
			//		write_log_info(_log, 1, "NOTE", "Missing observations for double freq");
			//		return false;
			//	}
			//	// triple freq, band1 should have ObsL and ObsP
			//	if (isTripleFreq && (bandX == band1 || bandX == band2) && (double_eq(ObsP, 0.0) || double_eq(ObsL, 0.0)))
			//	{
			//		write_log_info(_log, 1, "NOTE", "Missing Band observations for Triple freq");
			//		return false;
			//	}
			//	// triple freq, add valid
			//	if (isTripleFreq && !double_eq(ObsP, 0.0) && !double_eq(ObsL, 0.0))
			//	{
			//		valid_band.insert(band.second);
			//	}
			//}

			//// triple freq, but noly band1 valid, return false
			//if (isTripleFreq && (valid_band.size() < 2))
			//{
			//	write_log_info(_log, 1, "NOTE", "Missing observations for Triple freq");
			//	return false;
			//}
			//break;

			//// jdhuang : fix  ==========================================================================
		}
		case gnut::WL_COMBIN: // add by BoWong 2019/11/13
		{
			gobs1 = t_gobs(obsdata.select_range(b1));
			gobs2 = t_gobs(obsdata.select_range(b2));
			ObsP = obsdata.PWL(gobs1, gobs2);
			gobs1 = t_gobs(obsdata.select_phase(b1));
			gobs2 = t_gobs(obsdata.select_phase(b2));
			ObsL = obsdata.LWL(gobs1, gobs2);
			if (double_eq(ObsP, 0.0) || double_eq(ObsL, 0.0))
			{
				write_log_info(_log, 1, "NOTE", "no range or phase data");
				return false;
			}
			break;
		}
		default:
		{
			write_log_info(_log, 1, "NOTE", "no obs mode valid");
			return false;
		}
		}

		return true;
	}


	bool apply_reldelay(t_gtriple crd_site, t_gtriple vel_site, t_gtriple crd_sat, t_gtriple vel_sat,
		double& reldelay)
	{
		//reldelay = 2.0*(DotProduct(crd_sat.crd_cvect(), vel_sat.crd_cvect()) - DotProduct(crd_site.crd_cvect(), vel_site.crd_cvect())) / CLIGHT;
		reldelay = 2.0*(DotProduct(crd_sat.crd_cvect(), vel_sat.crd_cvect()))/ CLIGHT;
		ColumnVector xsat = crd_sat.crd_cvect();
		ColumnVector xsite = crd_site.crd_cvect();

		double r = xsite.norm_Frobenius() + xsat.norm_Frobenius();
		ColumnVector xsat2site = (xsite - xsat);
		double r_site2sat = xsat2site.norm_Frobenius();

		reldelay += 2.0*GM_CGCS / CLIGHT / CLIGHT * log((r + r_site2sat) / (r - r_site2sat));

		return true;

	}


	bool apply_IFCB_GPS(t_gsatdata& satdata, t_gifcb* ifcb, double& L, OBSCOMBIN obscombin)
	{
		t_gtime epoch = satdata.epoch();
		string sat = satdata.sat();
		double ifcb_GPS;


		GOBSBAND b1 = t_gsys::band_priority(satdata.gsys(), FREQ_1);
		GOBSBAND b2 = t_gsys::band_priority(satdata.gsys(), FREQ_2);
		GOBSBAND b3 = t_gsys::band_priority(satdata.gsys(), FREQ_3);

		double l1 = satdata.wavelength(b1);
		double l2 = satdata.wavelength(b2);
		double l3 = satdata.wavelength(b3);

		double c2 = l1 / l3;
		double c1 = 1.0 / (1.0 - c2 * c2);
		c2 = c2 * c2 * c1;

		if (satdata.gsys() == GSYS::GPS)
		{
			// liugege
			one_epoch_ifcb epoch_ifcb = ifcb->get_epo_ifcb(epoch);
			if (epoch_ifcb.find(sat) == epoch_ifcb.end()) return false;

			if (ifcb->get_ifcb_sigma(epoch, sat) > 0.2 ||
				ifcb->get_ifcb_npoint(epoch, sat) <= 2)
			{
				return false;
			}
			ifcb_GPS = ifcb->get_ifcb_value(epoch, sat);
			if (obscombin == IONO_FREE)
			{
				L += ifcb_GPS;
			}
			else if (obscombin == RAW_ALL)
			{
				L += ifcb_GPS / (-c2);
			}

		}

		return true;
	}



	bool correct_DCB(string ac, t_gallbias* allbias, t_gsatdata& satdata, double& P, OBSCOMBIN obscombin, t_gobs* gobs1, t_gobs* gobs2, t_gobs* gobs3)
	{
		if (!allbias || !gobs1 || obscombin != RAW_ALL && !gobs2)
		{
			return false;
		}

		GOBSBAND band1, band2, band3;
		t_gobs fobs1, fobs2, fobs3;            // local gobs
		double wavelen1, wavelen2, wavelen3;
		double coef1, coef2;                   // coefficient

		// constant
		t_gtime epoch = satdata.epoch();
		string sat = satdata.sat();
		GSYS gsys = satdata.gsys();

		band1 = t_gsys::band_priority(gsys, FREQ_1);
		band2 = t_gsys::band_priority(gsys, FREQ_2);
		band3 = t_gsys::band_priority(gsys, FREQ_3);

		// Judge gobs
		if (gobs1->attr() == ATTR) fobs1 = satdata.id_range(gobs1->band());      // automatic GOBS selection
		else fobs1 = gobs1->gobs();                                              // specific GOBS
		fobs1.gobs2to3(gsys);                                                    // 2.xx to 3.xx

		if (gobs2 != 0) {
			if (gobs2->attr() == ATTR) fobs2 = satdata.id_range(gobs2->band());  // automatic GOBS selection
			else fobs2 = gobs2->gobs();                                          // specific GOBS
			fobs2.gobs2to3(gsys);                                                // 2.xx to 3.xx
		}
		if (gobs3 != 0) {
			if (gobs3->attr() == ATTR) fobs3 = satdata.id_range(gobs3->band());  // automatic GOBS selection
			else fobs3 = gobs3->gobs();                                          // specific GOBS
			fobs3.gobs2to3(gsys);                                                // 2.xx to 3.xx
		}

		// coefficients
		wavelen1 = satdata.wavelength(band1);
		wavelen2 = satdata.wavelength(band2);
		coef2 = wavelen1 / wavelen2;
		coef1 = 1.0 / (1.0 - coef2 * coef2);
		coef2 = coef2 * coef2 * coef1;

		// UDUC
		if (obscombin == RAW_ALL) 
		{
			FREQ_SEQ freq = t_gsys::band2freq(gsys, fobs1.band());
			if (freq == FREQ_1) {
				fobs2 = t_gobs(satdata.select_range(band2));
				fobs2.gobs2to3(gsys);   // 2.xx to 3.xx
				double corr = allbias->get(epoch, sat, fobs1.gobs(), fobs2.gobs(), ac);
				P -= -coef2 * corr;
			}
			else if (freq == FREQ_2) {
				fobs2 = fobs1;
				fobs1 = t_gobs(satdata.select_range(band1));
				fobs1.gobs2to3(gsys);   // 2.xx to 3.xx
				double corr = allbias->get(epoch, sat, fobs1.gobs(), fobs2.gobs(), ac);
				P -= -coef1 * corr;
			}
			else if (freq == FREQ_3) {
				fobs3 = fobs1;
				fobs1 = t_gobs(satdata.select_range(band1));
				fobs2 = t_gobs(satdata.select_range(band2));
				fobs1.gobs2to3(gsys);   // 2.xx to 3.xx
				fobs2.gobs2to3(gsys);   // 2.xx to 3.xx
				double corr1 = allbias->get(epoch, sat, fobs1.gobs(), fobs2.gobs(), ac);
				double corr2 = allbias->get(epoch, sat, fobs1.gobs(), fobs3.gobs(), ac);
				P -= -coef1 * corr2 + coef2 * (corr2 - corr1);
			}
			else if (freq == FREQ_4) {
				// To do ....
			}
		}
		// Combination
		else
		{
			// Triple IF DCB correction
			if (gobs3 != 0) {
				if (obscombin == IONO_FREE && gsys != GLO) {
					double corr1 = allbias->get(epoch, sat, fobs1.gobs(), fobs2.gobs(), ac);
					double corr2 = allbias->get(epoch, sat, fobs1.gobs(), fobs3.gobs(), ac);
					wavelen3 = satdata.wavelength(band3);
					double coef3, coef4;
					coef4 = wavelen1 / wavelen3;
					coef3 = 1.0 / (1.0 - coef4 * coef4);
					coef4 = coef4 * coef4 * coef3;
					P -= -coef2 * corr1 + coef4 * corr2;
				}
			}
			// dual 12/13
			else
			{
				double corr = allbias->get(epoch, sat, fobs1.gobs(), fobs2.gobs(), ac);
				double corr1 = coef2 * corr;
				double corr2 = coef1 * corr;

				if (obscombin == MW_COMBIN)
				{
					double fact = wavelen1 / wavelen2;
					P += -(corr1 / wavelen1 + corr2 / wavelen2) * (1.0 - fact) / (1.0 + fact);
				}
				else if (obscombin == EWL_COMBIN)
				{
					wavelen3 = satdata.wavelength(band3);
					double alpha = (wavelen2 * wavelen3 - wavelen1 * wavelen1) / (wavelen2 * wavelen2 - wavelen1 * wavelen1);
					P += -(corr1 * (1 - alpha) + corr2 * alpha) * (wavelen3 - wavelen2) / (wavelen2 * wavelen3);
				}
				else if (obscombin == IONO_FREE)
				{
					double alfa, beta;
					satdata.coef_ionofree(band1, alfa, band2, beta);
					P += alfa * corr1 + beta * corr2;
				}

			}
		}
		return true;
	}


	int Check_Range(vector<double>& omc, double& recclk, double& sigma)
	{
		double mean = 0.0;
		double sig = 0.0;
		vector<double> res;

		for (auto it = omc.begin(); it != omc.end();) {
			if (fabs(*it) < 10e-8) {
				it = omc.erase(it);
				continue;
			}
			mean += *it;
			it++;
		}

		if (omc.size() == 1 || (omc.size() == 2 && fabs(omc[0] - omc[1]) > 1000)) {
			recclk = 0.0;
			return 1;
		}

		mean /= omc.size();

		for (auto it = omc.begin(); it != omc.end(); it++) {
			sig += (*it - mean) * (*it - mean);
			res.push_back(fabs(*it - mean));
		}
		sig = sqrt(sig / (omc.size() - 1));

		double sig1 = 2000;
		double mean1 = 0.0;
		bool isFound = true;

		while (isFound && sig1 > 1000 && omc.size() >= 3)
		{
			auto maxRes = max_element(res.begin(), res.end());
			auto maxPos = omc.begin() + distance(res.begin(), maxRes);
			double MaxValue = *maxPos;

			mean1 = 0.0;
			sig1 = 0.0;
			res.clear();
			for (auto it = omc.begin(); it != omc.end(); it++) {
				if (it == maxPos) continue;
				mean1 += *it / (omc.size() - 1);
			}

			for (auto it = omc.begin(); it != omc.end(); it++) {
				if (it == maxPos) continue;
				sig1 += (*it - mean1) * (*it - mean1);
				res.push_back(fabs(*it - mean1));
			}
			sig1 = sqrt(sig1 / (omc.size() - 2));

			isFound = fabs(MaxValue - mean1) > 3 * (sig1 > 1000 / 3.0 ? sig1 : 1000 / 3.0);
			if (!isFound) continue;
			mean = mean1;
			sig = sig1;
			omc.erase(maxPos);
		}

		recclk = mean;
		sigma = sig;
		return omc.size();
	}
}
