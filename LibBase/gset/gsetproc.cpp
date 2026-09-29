/*
*
* @verbatim
	 (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  @endverbatim
*
* @file		gsetproc.cpp
* @brief	implements process setting class
* @author   Jan Dousa
* @version	1.0.0
* @date		2012-10-23
*
*/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>

#include "gset/gsetproc.h"
#include "gutils/gstring.h"
#include "gproc/gsppflt.h"

#include <omp.h>

using namespace std;
using namespace pugi;

#define XMLKEY_PROC_IFCB "ifcb_model"
#define XMLKEY_PROC_ION_ORDER  "iono_order"
#define XMLKEY_PROC_IFB  "ifb_model"
#define XMLKEY_PROC_DCB  "dcb_model"
#define XMLKEY_PROC_FREQ "freq_model"
#define XMLKEY_PROC_READ_OFILE "read_ofile_mode"
using namespace great;

namespace gnut 
{
	// Constructor
	// ----------
	t_gsetproc::t_gsetproc()
		: t_gsetbase()
	{
		_set.insert(XMLKEY_PROC);
		_phase = true;
		_tropo = true;
		_iono = true;
		_tropo_grad = false;
		_tropo_slant = false;
		_tropo_model = SAASTAMOINEN;
		_ztd_model = PWC;
		_tropo_mf = GMF;
		_iono_mf = ICOSZ;
		_grad_mf = TILTING;
		_obs_weight = SINEL2;
		_res_type = RES_NORM;
		_obs_combin = IONO_FREE;
		_cbiaschar = ORIG;

		_sig_init_ztd = 0.1;
		_sig_init_vion = 10;
		_sig_init_grd = 0.0005;
		_sig_init_amb = 1000.0;
		_sig_init_crd = 100.0;
		_sig_init_vel = 10.0;
		_sig_init_glo = 1000.0;
		_sig_init_gal = 1000.0;
		_sig_init_bds = 1000.0;
		_sig_init_qzs = 1000.0;
		_minimum_elev = 10;
		_max_res_norm = 10;
		_crd_est = "EST";
		_pos_kin = false;
		_use_eclipsed = true;
		_auto_band = false;
		_frequency = 2;
		_realtime = false;
		_sd_sat = false;
		_simulation = false;
		//rtk process
		_basepos = SPP;
		_nppmodel = NPP_MODEL::NO;
		_minsat = static_cast<size_t>(SPP_MINSAT);

	}


	// Destructor
	// ----------
	t_gsetproc::~t_gsetproc()
	{}

	// Return value
    // ----------
	bool   t_gsetproc::isleo()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("isleo").as_bool();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	bool   t_gsetproc::tropo()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("tropo").as_bool();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	bool   t_gsetproc::iono()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("iono").as_bool();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	bool   t_gsetproc::tropo_slant()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("tropo_slant").as_bool();

		_gmutex.unlock(); return tmp;
	}


	// Return value
	// ----------
	bool   t_gsetproc::tropo_grad()
	{
		_gmutex.lock();

		// temporary for legacy mode with keyword "gradient"
		bool tmp = false;
		bool tmp1 = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("tropo_grad").as_bool();
		bool tmp2 = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("gradient").as_bool();   // legacy mode
		if (tmp1 || tmp2) tmp = true;

		_gmutex.unlock(); return tmp;
	}


	// Return value
	// ----------
	bool   t_gsetproc::phase()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("phase").as_bool();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	bool   t_gsetproc::doppler()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("doppler").as_bool();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	bool   t_gsetproc::pos_kin()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("pos_kin").as_bool();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	bool   t_gsetproc::use_eclipsed()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("use_eclipsed").as_bool();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	bool   t_gsetproc::auto_band()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("auto_band").as_bool();

		_gmutex.unlock(); return tmp;
	}

	bool t_gsetproc::sd_sat()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sd_sat").as_bool();

		_gmutex.unlock(); return tmp;
	}
	// Return value
	// ----------
	bool   t_gsetproc::IRC()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("IRC").as_bool();

		_gmutex.unlock(); return tmp;
	}


	// Return value
    // ----------
	bool   t_gsetproc::pre_fix_residuals()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("pre_fix_residuals").as_bool();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	bool   t_gsetproc::carrier_range()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("carrier_range").as_bool();

		_gmutex.unlock(); return tmp;
	}

	bool t_gsetproc::apply_carrier_range()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("apply_carrier_range").as_bool();

		_gmutex.unlock(); return tmp;
	}


	bool   t_gsetproc::ambfix()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("ambfix").as_bool();

		_gmutex.unlock(); return tmp;
	}
	bool t_gsetproc::realtime()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("realtime").as_bool();

		_gmutex.unlock(); return tmp;
	}
	bool t_gsetproc::ultrasp3()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("ultrasp3").as_bool();

		_gmutex.unlock(); return tmp;
	}
	int t_gsetproc::dynamics()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("dynamics").as_bool();

		_gmutex.unlock(); return tmp;
	}
	BASEPOS t_gsetproc::basepos()
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("basepos").value();
		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
		BASEPOS type = str2basepos(tmp);
		_gmutex.unlock(); return type;
	}
	BASEPOS t_gsetproc::str2basepos(string str)
	{
		BASEPOS mode;
		if (str.compare("RINEXH") == 0) mode = RINEXH;
		else if (str.compare("SPP") == 0) mode = SPP;
		else if (str.compare("CFILE") == 0) mode = CFILE;
		else {
			mode = SPP;
		}
		return mode;
	}
	// Return value
	// ----------
	double t_gsetproc::sig_init_ztd()
	{
		_gmutex.lock();

		double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_ztd").as_double();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	double t_gsetproc::sig_init_vion()
	{
		_gmutex.lock();

		double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_vion").as_double();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	double t_gsetproc::sig_init_grd()
	{
		_gmutex.lock();

		double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_grd").as_double();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	double t_gsetproc::sig_init_crd()
	{
		_gmutex.lock();

		double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_crd").as_double();

		_gmutex.unlock(); return tmp;
	}

	double t_gsetproc::sig_init_vel()
	{
		_gmutex.lock();

		double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_vel").as_double();

		_gmutex.unlock(); return tmp;
	}


	// Return value
	// ----------
	double t_gsetproc::sig_init_amb()
	{
		_gmutex.lock();

		double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_amb").as_double();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	double t_gsetproc::sig_init_glo()
	{
		_gmutex.lock();

		double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_glo").as_double();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	double t_gsetproc::sig_init_gal()
	{
		_gmutex.lock();

		double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_gal").as_double();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	double t_gsetproc::sig_init_bds()
	{
		_gmutex.lock();

		double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_bds").as_double();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	double t_gsetproc::sig_init_qzs()
	{
		_gmutex.lock();

		double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_qzs").as_double();

		_gmutex.unlock(); return tmp;
	}

	double t_gsetproc::sig_ref_clk()
	{
		double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_ref_clk").as_double();

		return tmp;
	}

	// Return value
	// ----------
	double t_gsetproc::minimum_elev()
	{
		_gmutex.lock();

		double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("minimum_elev").as_double();

		_gmutex.unlock(); return tmp;
	}

	// Return value for LEO
	// ----------
	double t_gsetproc::minimum_elev_leo()
	{
		_gmutex.lock();

		double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("minimum_elev_leo").as_double();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	double t_gsetproc::max_res_norm()
	{
		_gmutex.lock();

		double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("max_res_norm").as_double();

		_gmutex.unlock(); return tmp;
	}

	int t_gsetproc::read_ofile_step()
	{
		_gmutex.lock();
		int tmp;
		if (_doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).child(XMLKEY_ROM).attribute("step").empty())
		{
			_gmutex.unlock();
			return 5000;
		}
		else
		{
			tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).child(XMLKEY_ROM).attribute("step").as_int();
		}
		_gmutex.unlock(); 
		return tmp;
	}

	int t_gsetproc::minsat()
	{
		_gmutex.lock();

		double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("min_sat").as_int();

		_gmutex.unlock(); return tmp;
	}

	string t_gsetproc::ref_clk()
	{
		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("ref_clk").value();
		return tmp;
	}

	string t_gsetproc::ref_rec()
	{
		return string();
	}

	LSQMODE t_gsetproc::lsq_mode()
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("lsq_mode").value();

		_gmutex.unlock(); 
		if (tmp == "LSQ")
		{
			return LSQ;
		}
		else
		{
			return EPO;
		}
	}

	SLIPMODEL t_gsetproc::slip_model()
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("slip_model").value();

		_gmutex.unlock();
		if (tmp == "turboedit")
		{
			return SLIPMODEL::TURBO_EDIT;
		}
		else
		{
			return SLIPMODEL::DEF_DETECT_MODEL;
		}
	}

	IONMODEL t_gsetproc::ion_model()
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("ion_model").value();

		_gmutex.unlock();
		if (tmp == "SION")
		{
			return IONMODEL::SION;
		}
		else if(tmp == "VION")
		{
			return IONMODEL::VION;
		}
		else
		{
			return IONMODEL::DEF_ION;
		}
	}

	// Return value
	// ----------
	TROPMODEL t_gsetproc::tropo_model()
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("tropo_model").value();

		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

		TROPMODEL TM = str2tropmodel(tmp);

		if (TM == DEF_TROPMODEL) {
			xml_node parent = _doc.child(XMLKEY_ROOT);
			xml_node node = _default_node(parent, XMLKEY_PROC);
			tmp = tropmodel2str(_tropo_model);
			_default_attr(node, "tropo_model", tmp);

			TM = _tropo_model;
		}

		_gmutex.unlock(); return TM;
	}

	ZTDMODEL t_gsetproc::ztd_model(double* dt)
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("ztd_model").value();

		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

		*dt = 0.0;
		ZTDMODEL TM = str2ztdmodel(tmp);
		if (TM == PWC) {
			int loc = tmp.find(":");
			*dt = str2dbl(tmp.substr(loc + 1));
		}

		if (TM == DEF_ZTDMODEL) {
			xml_node parent = _doc.child(XMLKEY_ROOT);
			xml_node node = _default_node(parent, XMLKEY_PROC);
			tmp = ztdmodel2str(_ztd_model);
			_default_attr(node, "tropo_model", tmp);

			TM = _ztd_model;
			*dt = 60;
		}

		_gmutex.unlock(); return TM;
	}

	// Return value
	// ----------
	string t_gsetproc::trop_corr()
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("trop_corr").value();

		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	ATTITUDES t_gsetproc::attitudes()
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("attitudes").value();

		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

		ATTITUDES at = str2attitudes(tmp);

		_gmutex.unlock(); return at;
	}

	// Return value
	// ----------
	CBIASCHAR t_gsetproc::cbiaschar()
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("cbiaschar").value();

		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

		CBIASCHAR cb = str2cbiaschar(tmp);

		_gmutex.unlock(); return cb;
	}


	// Return value
	// ----------
	SYSBIASMODEL t_gsetproc::sysbias_model()
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sysbias_model").value();

		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

		SYSBIASMODEL model = str2sysbiasmodel(trim(tmp));

		_gmutex.unlock(); return model;

	}


	// Return value
	// ----------
	SYSBIASMODEL t_gsetproc::ifbbias_model()
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("ifbbias_model").value();

		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

		SYSBIASMODEL model = str2sysbiasmodel(trim(tmp));

		_gmutex.unlock(); return model;

	}


	// Return value
	// ----------
	GRDMPFUNC t_gsetproc::grad_mf()
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("grad_mf").value();

		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

		GRDMPFUNC MF = str2grdmpfunc(tmp);

		if (MF == DEF_GRDMPFUNC) {
			xml_node parent = _doc.child(XMLKEY_ROOT);
			xml_node node = _default_node(parent, XMLKEY_PROC);
			tmp = grdmpfunc2str(_grad_mf);
			_default_attr(node, "grad_mf", tmp);

			MF = _grad_mf;
		}

		_gmutex.unlock(); return MF;
	}

	// Return value
	// ----------   
	OBSWEIGHT t_gsetproc::weighting()
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("obs_weight").value();

		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

		OBSWEIGHT WG = str2obsweight(tmp);

		if (WG == DEF_OBSWEIGHT) {
			xml_node parent = _doc.child(XMLKEY_ROOT);
			xml_node node = _default_node(parent, XMLKEY_PROC);
			tmp = obsweight2str(_obs_weight);
			_default_attr(node, "obs_weight", tmp);

			WG = _obs_weight;
		}

		_gmutex.unlock(); return WG;
	}

	// Return value
	// ----------   
	RESIDTYPE t_gsetproc::residuals()
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("residuals").value();

		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

		RESIDTYPE RS = str2residtype(tmp);

		if (RS == DEF_RESIDTYPE) {
			xml_node parent = _doc.child(XMLKEY_ROOT);
			xml_node node = _default_node(parent, XMLKEY_PROC);
			tmp = residtype2str(_res_type);
			_default_attr(node, "residuals", tmp);

			RS = _res_type;
		}

		_gmutex.unlock(); return RS;
	}

	// Return value
	// ----------   
	OBSCOMBIN t_gsetproc::obs_combin()
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("obs_combination").value();

		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

		OBSCOMBIN OC = str2obscombin(tmp);

		if (OC == DEF_OBSCOMBIN) {
			xml_node parent = _doc.child(XMLKEY_ROOT);
			xml_node node = _default_node(parent, XMLKEY_PROC);
			tmp = obscombin2str(_obs_combin);
			_default_attr(node, "obs_combination", tmp);

			OC = _obs_combin;
		}

		_gmutex.unlock(); return OC;
	}

	// Return value
   // ----------
	ZTDMPFUNC t_gsetproc::tropo_mf()
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("tropo_mf").value();

		_gmutex.unlock();

		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

		ZTDMPFUNC MF = str2ztdmpfunc(tmp);

		if (MF == DEF_ZTDMPFUNC) {
			xml_node parent = _doc.child(XMLKEY_ROOT);
			xml_node node = _default_node(parent, XMLKEY_PROC);
			tmp = ztdmpfunc2str(_tropo_mf);
			_default_attr(node, "tropo_mf", tmp);

			MF = _tropo_mf;
		}

		return MF;
	}

	// Return value
	// ----------
	IONMPFUNC t_gsetproc::iono_mf()
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("iono_mf").value();

		_gmutex.unlock();

		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

		IONMPFUNC MF = str2ionmpfunc(tmp);

		if (MF == DEF_IONMPFUNC) {
			xml_node parent = _doc.child(XMLKEY_ROOT);
			xml_node node = _default_node(parent, XMLKEY_PROC);
			tmp = ionmpfunc2str(_iono_mf);
			_default_attr(node, "iono_mf", tmp);

			MF = _iono_mf;
		}

		 return MF;
	}

	// Return value
	// ----------
	CONSTRPAR t_gsetproc::crd_est()
	{
		_gmutex.lock();

		CONSTRPAR constr = EST;
		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("crd_constr").value();

		_gmutex.unlock();

		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

		if (tmp.empty()) tmp = _crd_est; // default

		if (     tmp.compare("FIX") == 0) constr = FIX;
		else if (tmp.compare("EST") == 0) constr = EST;
		else if (tmp.compare("KIN") == 0) constr = KIN;
		else if (tmp.compare("SIMU_KIN") == 0) constr = SIMU_KIN;
		else constr = FIX;


		return constr;
	}

	// settings check
	// ----------
	void t_gsetproc::check()
	{
		_gmutex.lock();

		// check existence of nodes/attributes
		xml_node parent = _doc.child(XMLKEY_ROOT);
		xml_node node = _default_node(parent, XMLKEY_PROC);

		// check existence of attributes
		_default_attr(node, "phase", _phase);
		_default_attr(node, "tropo", _tropo);
		_default_attr(node, "iono", _iono);
		_default_attr(node, "tropo_grad", _tropo_grad);
		_default_attr(node, "tropo_model", _tropo_model);
		_default_attr(node, "tropo_slant", _tropo_slant);
		_default_attr(node, "sig_init_crd", _sig_init_crd);
		_default_attr(node, "sig_init_vel", _sig_init_vel);
		_default_attr(node, "sig_init_ztd", _sig_init_ztd);
		_default_attr(node, "sig_init_grd", _sig_init_grd);
		_default_attr(node, "sig_init_vion", _sig_init_vion);
		_default_attr(node, "sig_init_amb", _sig_init_amb);
		_default_attr(node, "sig_init_glo", _sig_init_glo);
		_default_attr(node, "sig_init_gal", _sig_init_gal);
		_default_attr(node, "sig_init_bds", _sig_init_bds);
		_default_attr(node, "sig_init_qzs", _sig_init_qzs);
		_default_attr(node, "minimum_elev", _minimum_elev);
		_default_attr(node, "max_res_norm", _max_res_norm);
		_default_attr(node, "crd_constr", _crd_est);
		_default_attr(node, "pos_kin", _pos_kin);
		_default_attr(node, "sd_sat", _sd_sat);
		_default_attr(node, "min_sat", _minsat);
		_default_attr(node, "use_eclipsed", _use_eclipsed);
		_default_attr(node, "auto_band", _auto_band);
		_default_attr(node, "realtime", _realtime);
		_default_attr(node, "simulation", _simulation);
		//rtk process
		_default_attr(node, "basepos", _basepos);
		//_default_attr(node, "nppmodel", _nppmodel);

		_gmutex.unlock(); return;
	}


	// help body
	// ----------
	void t_gsetproc::help()
	{
		_gmutex.lock();

		cerr << " <process \n"
			<< "   phase=\"" << _phase << "\" \n"
			<< "   tropo=\"" << _tropo << "\" \n"
			<< "   iono=\"" << _iono << "\" \n"
			<< "   tropo_grad=\"" << _tropo_grad << "\" \n"
			<< "   tropo_model=\"" << tropmodel2str(_tropo_model) << "\" \n"
			<< "   tropo_slant=\"" << _tropo_slant << "\" \n"
			<< "   tropo_mf=\"" << ztdmpfunc2str(_tropo_mf) << "\" \n"
			<< "   iono_mf=\"" << ionmpfunc2str(_iono_mf) << "\" \n"
			<< "   grad_mf=\"" << grdmpfunc2str(_grad_mf) << "\" \n"
			<< "   obs_weight=\"" << obsweight2str(_obs_weight) << "\" \n"
			<< "   residuals=\"" << residtype2str(_res_type) << "\" \n"
			<< "   obs_combination=\"" << obscombin2str(_obs_combin) << "\" \n"
			<< "   sig_init_crd=\"" << _sig_init_crd << "\" \n"
			<< "   sig_init_ztd=\"" << _sig_init_ztd << "\" \n"
			<< "   sig_init_vion=\"" << _sig_init_vion << "\" \n"
			<< "   minimum_elev=\"" << _minimum_elev << "\" \n"
			<< "   max_res_norm=\"" << _max_res_norm << "\" \n"
			<< "   basepos=\"" << _basepos << "\" \n"
			//<< "   nppmodel=\"" << _nppmodel << "\" \n"
			<< " />\n";

		cerr << "\t<!-- process description:\n"
			<< "\t phase [0,1]   .. use carrier phase data\n"
			<< "\t tropo [0,1]   .. estimate troposphere\n"
			<< "\t tropo_grad    .. tropospheric horizontal gradient models\n"
			<< "\t tropo_model   .. tropospheric model (SAAS, GPT, ...)\n"
			<< "\t tropo_slant   .. tropospheric slant delays produced\n"
			<< "\t tropo_mf      .. tropospheric mapping function (COSZ, GMF, ...)\n"
			<< "\t grad_mf       .. tropo gradient mapping function (TILTING,CHEN_HERRING, BAR_SEVER ...)\n"
			<< "\t obs_weight    .. observation elevation dependant weighting (EQUAL, SINEL, SINEL2, SINEL4, CODPHA, MLTPTH)\n"
			<< "\t residuals     .. type of residuals (RES_ORIG, RES_NORM, RES_ALL)\n"
			<< "\t sig_init_crd  .. accuracy of initial coordinates [m]\n"
			<< "\t sig_init_ztd  .. accuracy of initial zenith path delay [m]\n"
			<< "\t sig_init_vion .. accuracy of initial vertical iono path delay [m]\n"
			<< "\t minimum_elev  .. elevation angle cut-off [degree]\n"
			<< "\t max_res_norm  .. maximal normalized residuals\n"
			<< "\t basepos  .. base site coordinate\n"
			<< "\t nppmodel  .. npp model \n"
			<< "\t -->\n\n";

		_gmutex.unlock(); return;
	}


	// convert str to ATTITUDES enum
	// ----------   
	ATTITUDES t_gsetproc::str2attitudes(string ati)
	{
		ATTITUDES at;

		if (ati.compare("NOMINAL") == 0)  at = YAW_NOMI;
		else if (ati.compare("RTCM") == 0) at = YAW_RTCM;
		else {
			at = DEF_YAWMODEL;
			stringstream ostr;
			ostr << "Unsupported attitude model (" << at << ")! Used default value ( YAW_MODEL )";
			_add_log("gsetproc", ostr.str());
		}
		return at;
	}

	// convert str to GRDMPFUNC enum
	// ----------   
	GRDMPFUNC t_gsetproc::str2grdmpfunc(string mf)
	{
		GRDMPFUNC MF;

		if (mf.compare("TILTING") == 0) MF = TILTING;
		else if (mf.compare("CHEN_HERRING") == 0) MF = CHEN_HERRING;
		else if (mf.compare("BAR_SEVER") == 0) MF = BAR_SEVER;
		else {
			MF = DEF_GRDMPFUNC;
			stringstream ostr;
			ostr << "Unsupported GRD mapping function (" << mf << ")! Used default value (" << grdmpfunc2str(_grad_mf) << ")";
			_add_log("gsetproc", ostr.str());
		}
		return MF;
	}

	// convert str to CBIASCHAR enum
	// ----------   
	CBIASCHAR t_gsetproc::str2cbiaschar(string cb)
	{
		CBIASCHAR CB;

		if (cb.compare("2CHAR") == 0) CB = CHAR2;
		else if (cb.compare("3CHAR") == 0) CB = CHAR3;
		else if (cb.compare("ORIG") == 0) CB = ORIG;
		else {
			CB = DEF_CBIASCHAR;
			stringstream ostr;
			ostr << "Unsupported forcing code bias signals! Used default value (" << cbiaschar2str(_cbiaschar) << ")";
			_add_log("gsetproc", ostr.str());
		}
		return CB;
	}


	// convert str to SYSBIASMODEL enum
	// ----------   
	SYSBIASMODEL t_gsetproc::str2sysbiasmodel(string sys)
	{
		SYSBIASMODEL model;

		if      (sys.compare("AUTO+CON") == 0) model = SYSBIASMODEL::AUTO_CON;
		else if (sys.compare("AUTO+RWK") == 0) model = SYSBIASMODEL::AUTO_RWK;
		else if (sys.compare("ISB+CON")  == 0) model = SYSBIASMODEL::ISB_CON;
		else if (sys.compare("ISB+RWK")  == 0) model = SYSBIASMODEL::ISB_RWK;
		else if (sys.compare("AUTO+WHIT") == 0) model = SYSBIASMODEL::AUTO_WHIT;
		//else if (sys.compare("IFB+CON")  == 0) model = SYSBIASMODEL::IFB_CON;
		//else if (sys.compare("IFB+RWK")  == 0) model = SYSBIASMODEL::IFB_RWK;
		else {
			//	model = SYSBIASMODEL::AUTO_CON;
			model = SYSBIASMODEL::AUTO_WHIT;
			stringstream ostr;
			ostr << "Unsupported forcing code bias signals! Used default value (AUTO+CON)";
			_add_log("gsetproc", ostr.str());
		}
		return model;
	}

	// convert str to ZTDMPFUNC enum
	// ----------   
	ZTDMPFUNC t_gsetproc::str2ztdmpfunc(string mf)
	{
		ZTDMPFUNC MF;

		if (mf.compare("COSZ") == 0) MF = COSZ;
		else if (mf.compare("GMF") == 0) MF = GMF;
		else if (mf.compare("NO_MF") == 0) MF = NO_MF;
		else {
			MF = DEF_ZTDMPFUNC;
			stringstream ostr;
			ostr << "Unsupported ZTD mapping function (" << mf << ")! Used default value (" << ztdmpfunc2str(_tropo_mf) << ")";
			_add_log("gsetproc", ostr.str());
		}
		return MF;
	}

	// convert str to IONMPFUNC enum
	// ----------   
	IONMPFUNC t_gsetproc::str2ionmpfunc(string mf)
	{
		IONMPFUNC MF;

		if (mf.compare("COSZ") == 0) MF = ICOSZ;
		else if (mf.compare("QFAC") == 0) MF = QFAC;
		else if (mf.compare("NONE") == 0) MF = NONE;
		else {
			MF = DEF_IONMPFUNC;
			stringstream ostr;
			ostr << "Unsupported ION mapping function (" << mf << ")! Used default value (" << ionmpfunc2str(_iono_mf) << ")";
			_add_log("gsetproc", ostr.str());
		}
		return MF;
	}

	// convert str to OBSWEIGHT enum
	// ----------   
	OBSWEIGHT t_gsetproc::str2obsweight(string wg)
	{
		OBSWEIGHT WG;

		if (wg.compare("EQUAL") == 0) WG = EQUAL;
		else if (wg.compare("SINEL") == 0) WG = SINEL;
		else if (wg.compare("SINEL1") == 0) WG = SINEL;
		else if (wg.compare("SINEL2") == 0) WG = SINEL2;
		else if (wg.compare("SINEL4") == 0) WG = SINEL4;
		else if (wg.compare("CODPHA") == 0) WG = CODPHA;
		else if (wg.compare("MLTPTH") == 0) WG = MLTPTH;
		else if (wg.compare("PARTELE") == 0) WG = PARTELE;
		else if (wg.compare("SNR") == 0) WG = SNR;
		else {
			WG = DEF_OBSWEIGHT;
			stringstream ostr;
			ostr << "Unsupported observation weighting model (" << wg << ")! Used default value (" << obsweight2str(_obs_weight) << ")";
			_add_log("gsetproc", ostr.str());
		}
		return WG;
	}

	// convert str to TROPMODEL enum
	// ----------   
	TROPMODEL t_gsetproc::str2tropmodel(string tm)
	{
		TROPMODEL TM;

		if (tm.compare("SAASTAMOINEN") == 0) TM = SAASTAMOINEN;
		else if (tm.compare("DAVIS") == 0) TM = DAVIS;
		else if (tm.compare("HOPFIELD") == 0) TM = HOPFIELD;
		else if (tm.compare("MOPS") == 0) TM = MOPS;
		else if (tm.compare("GPTW") == 0) TM = GPTW;
		else if (tm.compare("GPT2W") == 0) TM = GPT2W;
		else if (tm.compare("GAL27") == 0) TM = GAL27;
		else if (tm.compare("GALTROPO27") == 0) TM = GALTROPO27;
		else if (tm.compare("EXTERN") == 0) TM = EXTERN;
		else {
			TM = DEF_TROPMODEL;
			stringstream ostr;
			ostr << "Unsupported tropospheric model (" << tm << ")! Used default value (" << tropmodel2str(_tropo_model) << ")";
			_add_log("gsetproc", ostr.str());
		}
		return TM;
	}

	ZTDMODEL t_gsetproc::str2ztdmodel(string ztd)
	{
		ZTDMODEL ZTD;
		if (ztd.substr(0, 3).compare("PWC") == 0) ZTD = PWC;
		else if (ztd.substr(0, 3).compare("STO") == 0) ZTD = STO;
		else {
			ZTD = DEF_ZTDMODEL;
		}
		return ZTD;
	}

	// convert str to RESIDTYPE enum
	// ----------   
	RESIDTYPE t_gsetproc::str2residtype(string rs)
	{
		RESIDTYPE RS;

		if (rs.compare("RES_ORIG") == 0) RS = RES_ORIG;
		else if (rs.compare("RES_NORM") == 0) RS = RES_NORM;
		else if (rs.compare("RES_ALL") == 0) RS = RES_ALL;
		else {
			RS = DEF_RESIDTYPE;
			stringstream ostr;
			ostr << "Unsupported type of residuals (" << rs << ")! Used default value (" << residtype2str(_res_type) << ")";
			_add_log("gsetproc", ostr.str());
		}
		return RS;
	}

	// convert str to OBSCOMB enum
	// ----------   
	OBSCOMBIN t_gsetproc::str2obscombin(string oc)
	{
		OBSCOMBIN OC;

		     if (oc.compare("IONO_FREE") == 0)        OC = IONO_FREE;
		else if (oc.compare("RAW_SINGLE") == 0)       OC = RAW_SINGLE;
		else if (oc.compare("RAW_DOUBLE") == 0)       OC = RAW_DOUBLE;
		else if (oc.compare("RAW_ALL") == 0)          OC = RAW_ALL;
		else if (oc.compare("MW_COMBIN") == 0)        OC = MW_COMBIN;
		else if (oc.compare("WL_COMBIN") == 0)        OC = WL_COMBIN;
		else if (oc.compare("RAW_MIX") == 0)          OC = RAW_MIX;
		else if (oc.compare("IF_P1") == 0)	          OC = IF_P1;
		else {
			OC = DEF_OBSCOMBIN;
			stringstream ostr;
			ostr << "Unsupported observations combination (" << oc << ")! Used default value (" << obscombin2str(_obs_combin) << ")";
			_add_log("gsetproc", ostr.str());
		}
		return OC;
	}

	// convert GRDMPFUNC enum to str
	// ---------
	string t_gsetproc::grdmpfunc2str(GRDMPFUNC MF)
	{
		switch (MF) {
		case TILTING:       return "TILTING";
		case CHEN_HERRING:  return "CHEN_HERRING";
		case BAR_SEVER:     return "BAR_SEVER";
		case DEF_GRDMPFUNC: return "NOT DEFINED";
		}
		return "";
	}

	// convert ZTDMPFUNC enum to str
	// ---------
	string t_gsetproc::ztdmpfunc2str(ZTDMPFUNC MF)
	{
		switch (MF) {
		case COSZ:          return "COSZ";
		case GMF:           return "GMF";
		case NO_MF:           return "NO_MF";
		case DEF_ZTDMPFUNC: return "NOT DEFINED";
		}
		return "";
	}

	// convert IONMPFUNC enum to str
	// ---------
	string t_gsetproc::ionmpfunc2str(IONMPFUNC MF)
	{
		switch (MF) {
		case COSZ:          return "COSZ";
		case QFAC:          return "QFAC";
		case NONE:          return "NONE";
		case DEF_IONMPFUNC: return "NOT DEFINED";
		}
		return "";
	}

	// convert OBSWEIGHT enum to str
	// ---------   
	string t_gsetproc::obsweight2str(OBSWEIGHT WG)
	{
		switch (WG) {
		case EQUAL:         return "EQUAL";
		case SINEL:         return "SINEL";
		case SINEL2:        return "SINEL2";
		case SINEL4:        return "SINEL4";
		case CODPHA:        return "CODPHA";
		case MLTPTH:        return "MLTPTH";
		case PARTELE:		return "PARTELE";
		case DEF_OBSWEIGHT: return "NOT DEFINED";
		}
		return "";
	}

	// convert TROPMODEL enum to str
	// ---------   
	string t_gsetproc::tropmodel2str(TROPMODEL TM)
	{
		switch (TM) {
		case SAASTAMOINEN:  return "SAASTAMOINEN";
		case DAVIS:         return "DAVIS";
		case HOPFIELD:      return "HOPFIELD";
		case MOPS:          return "MOPS";
		case GPTW:          return "GPTW";
		case GPT2W:         return "GPT2W";
		case GAL27:         return "GAL27";
		case GALTROPO27:    return "GALTROPO27";
		case EXTERN:        return "EXTERN";
		case DEF_TROPMODEL: return "NOT DEFINED";
		}
		return "";
	}

	// convert ZTD MODEL enum to str
	// ---------
	string t_gsetproc::ztdmodel2str(ZTDMODEL ZTD)
	{
		switch (ZTD)
		{
		case gnut::DEF_ZTDMODEL:
			return "NOT DEFINED";
		case gnut::PWC:
			return "PWC";
		case gnut::STO:
			return "STO";
		}
		return "";
	}

	// convert RESIDTYPE enum to str
	// ---------
	string t_gsetproc::residtype2str(RESIDTYPE RS)
	{
		switch (RS) {
		case RES_ORIG:      return "RES_ORIG";
		case RES_NORM:      return "RES_NORM";
		case RES_ALL:       return "RES_ALL";
		case DEF_RESIDTYPE: return "NOT DEFINED";
		}
		return "";
	}

	// convert OBSCOMB enum to str
	// ---------
	string t_gsetproc::obscombin2str(OBSCOMBIN OC)
	{
		switch (OC) {
		case IONO_FREE:             return "IONO_FREE";
		case RAW_SINGLE:            return "RAW_SINGLE";
		case RAW_DOUBLE:            return "RAW_DOUBLE";
		case RAW_ALL:               return "RAW_ALL";
		case DEF_OBSCOMBIN:         return "NOT DEFINED";
		case MW_COMBIN:             return "MW_COMBIN";
		case RAW_MIX:               return "RAW_MIX";
		case IF_P1:                 return "IF_P1";
		}
		return "";
	}

	// convert CBIASCHAR enum to str
	// ---------
	string t_gsetproc::cbiaschar2str(CBIASCHAR CB)
	{
		switch (CB) {
		case CHAR2:         return "2CHAR";
		case CHAR3:         return "3CHAR";
		case ORIG:          return "ORIG";
		case DEF_CBIASCHAR: return "NOT DEFINED";
		}
		return "";
	}

	t_gobscombtype::t_gobscombtype():
		_obs_type(GOBSTYPE::TYPE),
		_obs_band(GOBSBAND::BAND),
		_obs_combine(OBSCOMBIN::DEF_OBSCOMBIN)
	{
	}

	t_gobscombtype::t_gobscombtype(const t_gobscombtype& other) :
		_obs_type(other._obs_type),
		_obs_band(other._obs_band),
		_obs_combine(other._obs_combine)
	{
	}

	t_gobscombtype::t_gobscombtype(string obscombtype):
		_obs_type(GOBSTYPE::TYPE),
		_obs_band(GOBSBAND::BAND),
		_obs_combine(OBSCOMBIN::DEF_OBSCOMBIN)
	{
		if (obscombtype.length() < 2) {
			// TODO: need throwing error in constructor
			// throw ...
			return;
		}
		string str_type = obscombtype.substr(0, 1);
		if (str_type == "P") {
			str_type == "C";
		}
		_obs_type = str2gobstype(str_type);

		if (obscombtype == "SLR")
		{
			_obs_type = str2gobstype(obscombtype.substr(2,1));
		}

		string str_band = obscombtype.substr(1);
		if (str_band == "C") {
			_obs_combine = IONO_FREE;
		}
		else {
			_obs_band = str2gobsband(str_band);
		}
	}

	t_gobscombtype::t_gobscombtype(t_gobs obstype, OBSCOMBIN combtype):
		_obs_type(obstype.type()),
		_obs_band(obstype.band()),
		_obs_combine(combtype)
	{
	}

	t_gobscombtype::t_gobscombtype(GOBSTYPE t, GOBSBAND b, OBSCOMBIN obscomb):
		_obs_type(t),
		_obs_band(b),
		_obs_combine(obscomb)
	{
	}

	string t_gobscombtype::convert2str() const
	{
		string ans("");
		ans += gobstype2str(_obs_type==TYPE_C?TYPE_P:_obs_type);
		if (_obs_combine == IONO_FREE) {
			ans += "C";
		}
		else {
			ans += gobsband2str(_obs_band);
		}
		return ans;
	}

	bool t_gobscombtype::operator==(const t_gobscombtype & g) const
	{
		if (this->convert2str() == g.convert2str()) {
			return true;
		}
		else {
			return false;
		}

	}

	bool t_gobscombtype::operator<(const t_gobscombtype & g) const
	{
		if (this->convert2str() < g.convert2str()) {
			return true;
		}
		else {
			return false;
		}
		
	}

	bool t_gobscombtype::is_phase() const
	{
		if (_obs_type == TYPE_L) {
			return true;
		}
		else {
			return false;
		}
	}

	bool t_gobscombtype::is_SLR() const
	{
		if (_obs_type == TYPE_SLR) {
			return true;
		}
		else {
			return false;
		}
	}

	bool t_gobscombtype::is_code() const
	{
		if (_obs_type == TYPE_C || _obs_type == TYPE_P) {
			return true;
		}
		else {
			return false;
		}
	}

	///< judge whether simulation data
	bool t_gsetproc::simulation()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("simulation").as_bool();

		_gmutex.unlock(); return tmp;
	}

	// jdhuang : add
	IFCB_MODEL t_gsetproc::ifcb_model()
	{
		//get rkf model node
		string ifcb = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).child_value(XMLKEY_PROC_IFCB);

		if (ifcb.empty()) return IFCB_MODEL::COR;

		// delete spaces
		ifcb.erase(remove(ifcb.begin(), ifcb.end(), ' '), ifcb.end());
		ifcb = str2upper(ifcb);

		if (ifcb == "EST") return IFCB_MODEL::EST;
		if (ifcb == "COR") return IFCB_MODEL::COR;

		return IFCB_MODEL::DEF;
	}

	IONO_ORDER t_gsetproc::iono_order()
	{
		//get rkf model node
		string order = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).child_value(XMLKEY_PROC_ION_ORDER);

		if (order.empty()) return IONO_ORDER::NONE;

		// delete spaces
		order.erase(remove(order.begin(), order.end(), ' '), order.end());
		order = str2upper(order);

		if (order == "FIRST")   return IONO_ORDER::FIRST_ORDER;
		if (order == "SECOND")  return IONO_ORDER::SECOND_ORDER;
		if (order == "THIRD")   return IONO_ORDER::THIRD_ORDER;
		if (order == "NONE")    return IONO_ORDER::NONE;

		throw logic_error("ERROR type of iono order, you should use FIRST/SECOND/NONE");
	}

	NPP_MODEL t_gsetproc::npp_model() {
		//get rkf model node
		string npp_model = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("npp_model").value();

		if (npp_model.empty()) return NPP_MODEL::NO;

		// delete spaces
		npp_model.erase(remove(npp_model.begin(), npp_model.end(), ' '), npp_model.end());
		npp_model = str2upper(npp_model);

		if (npp_model == "NO")   return NPP_MODEL::NO;
		if (npp_model == "VRS")  return NPP_MODEL::VRS;
		if (npp_model == "URTK")   return NPP_MODEL::URTK;
		if (npp_model == "PPP_RTK")    return NPP_MODEL::PPP_RTK;

		throw logic_error("ERROR type of iono order, you should use FIRST/SECOND/NONE");
	}

	string t_gsetproc::range_smooth_mode(int* smt_windows)
	{
		*smt_windows = 20;
		xml_node tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).child("smooth_range");

		string tmp_mode = tmp_set.attribute("mode").value();
		if (tmp_mode.empty()) return "NONE";
		
		// delete spaces
		tmp_mode.erase(remove(tmp_mode.begin(), tmp_mode.end(), ' '), tmp_mode.end());
		tmp_mode = str2upper(tmp_mode);

		if (tmp_mode != "DOPPLER" && tmp_mode != "PHASE") return "NONE";

		string tmp_window = tmp_set.attribute("window").value();
		if (tmp_window.empty()) return tmp_mode;

		// delete spaces
		tmp_window.erase(remove(tmp_window.begin(), tmp_window.end(), ' '), tmp_window.end());
		int win_int = str2int(tmp_window);
		*smt_windows = win_int <= 0 ? 20 : win_int;

		return tmp_mode;
	}

	bool t_gsetproc::bds_code_bias_correction()
	{
		_gmutex.lock();

		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("bds_code_bias_corr").as_bool();

		_gmutex.unlock(); return tmp;
	}

	FREQ_MODEL t_gsetproc::freq_model()
	{
		//get rkf model node
		string ucud_model = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).child_value(XMLKEY_PROC_FREQ);

		if (ucud_model.empty()) throw logic_error("freq model is empty.");

		// delete spaces
		ucud_model = great::trim(ucud_model);
		ucud_model = str2upper(ucud_model);

		if (ucud_model == "DOUBLE") return FREQ_MODEL::DOUBLE;
		if (ucud_model == "TRIPLE") return FREQ_MODEL::TRIPLE;

		return FREQ_MODEL::DOUBLE;
	}

	IFB_MODEL t_gsetproc::ifb_model()
	{
		//get rkf model node
		string ifb = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).child_value(XMLKEY_PROC_IFB);

		if (ifb.empty()) return IFB_MODEL::NONE;

		// delete spaces
		ifb.erase(remove(ifb.begin(), ifb.end(), ' '), ifb.end());
		ifb = str2upper(ifb);

		if (ifb == "EST_REC_IFB") return IFB_MODEL::EST_REC_IFB;
		if (ifb == "EST_SAT_IFB") return IFB_MODEL::EST_SAT_IFB;
		if (ifb == "EST_IFB")     return IFB_MODEL::EST_IFB;
		if (ifb == "COR_REC_IFB") return IFB_MODEL::COR_REC_IFB;
		if (ifb == "COR_SAT_IFB") return IFB_MODEL::COR_SAT_IFB;
		if (ifb == "COR_IFB")     return IFB_MODEL::COR_IFB;
		if (ifb == "NONE")        return IFB_MODEL::NONE;
		
		throw logic_error("ERROR type of ifb model, you should use EST_REC_IFB/EST_SAT_IFB/EST_IFB/COR_REC_IFB/COR_SAT_IFB/COR_IFB/NONE");
	}

	//DCB_MODEL t_gsetproc::dcb_model()
	//{
	//	//get rkf model node
	//	string dcb = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).child_value(XMLKEY_PROC_DCB);

	//	if (dcb.empty()) throw logic_error("ifb model is empty.");

	//	// delete spaces
	//	dcb.erase(remove(dcb.begin(), dcb.end(), ' '), dcb.end());
	//	dcb = str2upper(dcb);

	//	if (dcb == "CODE") return DCB_MODEL::CODE;
	//	if (dcb == "CAS")  return DCB_MODEL::CAS;
	//	if (dcb == "NONE") return DCB_MODEL::NONE;

	//	throw logic_error("ERROR type of dcb model, you should use CODE/CAS/NONE");
	//}

	// Hwang Shih: add
	READ_OFILE_MODE t_gsetproc::read_ofile_mode()
	{
		//get read O files mode node
		string read_ofile = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).child_value(XMLKEY_PROC_READ_OFILE);
		if (read_ofile.empty()) return READ_OFILE_MODE::READALL;

		// delete spaces
		read_ofile.erase(remove(read_ofile.begin(), read_ofile.end(), ' '), read_ofile.end());
		read_ofile = str2upper(read_ofile);

		if (read_ofile == "REALTIME") return READ_OFILE_MODE::REALTIME;
		if (read_ofile == "READALL")  return READ_OFILE_MODE::READALL;

		throw logic_error("ERROR type of read o file mode, you should use REALTIME/READALL");
	}


	int t_gsetproc::frequency()
	{
		_gmutex.lock();

		int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("frequency").as_int();

		_gmutex.unlock();
		if (tmp == 0)
		{
			return 2;
		}
		else
		{
			return tmp;
		}
	}

	int t_gsetproc::num_threads()
	{
		_gmutex.lock();

		int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("num_threads").as_int();

		_gmutex.unlock();

#ifdef USE_OPENMP
		if (tmp == 0) {
			// default set 1
			tmp = 1;
		}
		else if (tmp==-1) {
			// set number of cpu cores
			tmp = omp_get_num_procs();
		}
#else
		// if no use openmp number of threads is 1 
		tmp = 1;
#endif
		return tmp;
	}

	GRDMODEL t_gsetproc::grd_model(double* dt)
	{
		_gmutex.lock();

		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("grd_model").value();

		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

		*dt = 0.0;
		GRDMODEL TM = str2grdmodel(tmp);
		if (TM == GRD_PWC) {
			int loc = tmp.find(":");
			*dt = str2dbl(tmp.substr(loc + 1));
		}

		if (TM == DEF_GRDMODEL) {

			TM = GRD_PWC;
			*dt = 120;
		}

		_gmutex.unlock(); return TM;
	}

	GRDMODEL t_gsetproc::str2grdmodel(string grd)
	{
		GRDMODEL GRD;
		string tmp = trim(grd);
		if (tmp.substr(0, 3).compare("PWC") == 0) GRD = GRD_PWC;
		else if (tmp.substr(0, 3).compare("STO") == 0) GRD = GRD_STO;
		else {
			GRD = DEF_GRDMODEL;
		}
		return GRD;
	}
	
} // namespace
