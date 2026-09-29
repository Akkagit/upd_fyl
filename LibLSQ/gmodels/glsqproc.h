/**
*
* @verbatim
	History
	 -1.0 ZhengHJ      2019-04-04  created
	 -1.1 glfeng       2019-09-05  add function _batch_turboedit and variable _turboedit_software
	 -1.2 Huang Jiande 2020-04-12  merge UC and IF model
	 -1.3 Huang Jiande 2020-06-21  Remove glsqprocIF and glsqprocUC
	 -1.4 Huang Jiande 2020-09-11  Format the functions
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    glsqproc.h
* @brief		virtual class for process
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-18
*
*/


#ifndef GPROCMODEL_H
#define GPROCMODEL_H

#include "gutils/ExportLSQ.h"
#include "gdata/gdata.h"
#include "gall/gallproc.h"
#include "gall/gallprod.h"
#include "gall/gallrb.h"
#include "gall/gallsta.h"
#include "gall/gallobj.h"
#include "gall/gallobs.h"
#include "gall/gallics_nav.h"
#include "gall/gallics_leo.h"
#include "gall/gallsatparam.h"
#include "gutils/gleo.h"
#include "gproc/glsq.h"
#include "gmodels/gmodel.h"
#include "gmodels/gprecisemodel.h"
#include "gmodels/gpceprecisemodel.h"
#include "gmodels/gpppleomodel.h"
#include "gmodels/gcombmodel.h"
#include "gproc/glsqmatrix.h"
#include "gdata/gsatdata.h"
#include "gutils/gcycleslip.h"
#include "gset/gsetamb.h"
#include "gambfix/gambiguity.h"
#include "gdata/gupd.h"
#include "gdata/gion.h"
#include "gdata/gambcon.h"
#include "gall/gallambupd.h"
#include "gmodels/gslrmodel.h"
#include "gdata/gslrf.h"

namespace great
{
	typedef vector< tuple < string, string, string > > t_lsq_equ_info;

	/**
	* @brief virtual class for process
	*/
	class LSQ_LIBRARY_EXPORT t_glsqproc
	{
	public:
		t_glsqproc();
		t_glsqproc(t_gsetbase* set, t_gallproc* data, t_glog* log);
		virtual ~t_glsqproc();

		void add_coder(const vector<t_gcoder*>& coder); //added by liyuhao

		virtual bool ProcessBatch(t_gallproc* data, const t_gtime&beg, const t_gtime& end);
		bool processOneEpoch(const t_gtime& crt_epoch, std::vector<t_gsatdata>& crt_obs) { return _processOneEpoch(crt_epoch, crt_obs); };

		virtual bool GenerateProduct();

		/**
		* @brief batch processing .
		* @param[in] useIon   whether use ion file.
		*/
		void setUseIon(bool useIon) { _useIon_cmd = useIon; };
		
		/**
		* @brief batch processing .
		* @param[in] ambfix   whether fix amb.
		*/
		void setUseAmbcon(bool ambfix) { _ambfix_cmd = ambfix; };


		set<string>& sat_list() { return _sat_list; }
		set<string>& rec_list() { return _rec_list; }
		set<string>& sys_list() { return _sys_list; }

		bool                 get_slr()      { return _slr; };
		double	             get_obs_intv() { return _obs_intv; };
		set<string>          get_rec_list() { return _rec_list; };
		t_glog*              get_glog()     { return _glog; };
		t_gtime              get_beg_time() { return _beg_time; };
		t_gtime              get_end_time() { return _end_time; };
		t_gtime              get_crt_time() { return _crt_time; };
		t_gallobs*           get_gall_obs() { return _gall_obs; };
		t_glsq*              get_lsq()      { return _lsq; };
		t_gallrecover*       get_gall_rcv() { return _gall_rcv; };
		t_gslrmodel*         get_slr_model(){ return _slr_model; };
		t_glsqEquationMatrix get_crt_equ()  { return _crt_equ; };

		void set_crt_time(t_gtime time) { _crt_time = time; };
		void set_end_time(t_gtime time) { _end_time = time; };
		void set_beg_time(t_gtime time) { _beg_time = time; };
		void set_gall_obs(t_gallobs* obs) { _gall_obs = obs; };

	protected:

		t_glog*     _glog = nullptr;	///< log file for process
		t_gsetbase* _gset = nullptr;	///< setting from xml file

		/// The following can get from xml file
		set<string>             _rec_list;				  ///< all sites of proc
		set<string>             _sat_list;				  ///< all sats of proc
		set<string>             _sys_list;				  ///< sat systems of proc
		set<string>             _leo_list;                ///< short name of leos
		vector<string>          _rover_list;              ///< rover sites for PPPRTK
		map<string, t_gtriple>  _rec_crds;                ///< site coordinates
		map<string, t_gtriple>  _rec_stds;                ///< site coordinates std

		int     _frequency    = 0;
		double	_obs_intv     = 0.0;					  ///< sampling of proc
		double  _pwc_intv     = 0.0;                      ///< sampling of pwc
		double  _grd_pwc_intv = 0.0;                      ///< sampling of pwc for gradients, glfeng 
		t_gtime _beg_time;					            ///< beg time of proc
		t_gtime _end_time;				                ///< end time of proc
		t_gtime _crt_time;					            ///< current epoch of proc

		map< GSYS, map<FREQ_SEQ, GOBSBAND> > _band_index;
		map< GSYS, map<GOBSBAND, FREQ_SEQ> > _freq_index;

		// jdhuang : 全部加上了初值，随便加的

		LSQMODE              _lsq_mode = LSQMODE::LSQMODE_UNDEF;
		OBSCOMBIN            _obs_mode = OBSCOMBIN::DEF_OBSCOMBIN;
		CONSTRPAR            _crd_est = CONSTRPAR::CONSTRPAR_UNDEF;
		t_gsetamb::FIX_MODE  _fix_mode = t_gsetamb::FIX_MODE::NO;
		t_gsetamb::UPD_MODE  _upd_mode = t_gsetamb::UPD_MODE::OSB;
		SLIPMODEL            _slip_model = SLIPMODEL::SLIPMODEL_UNDEF;
		ZTDMODEL             _ztd_model = ZTDMODEL::DEF_ZTDMODEL;
		SYSBIASMODEL         _sysbias_model = SYSBIASMODEL::AUTO_WHIT;
		IFCB_MODEL           _ifcb_model = IFCB_MODEL::IFCB_MODEL_UNDEF;
		IONO_ORDER           _iono_order = IONO_ORDER::NONE;
		IFB_MODEL            _ifb_model = IFB_MODEL::NONE;
		DD_MODEL             _dd_model = DD_MODEL::NONE;                        ///< dd fix mode
		GRDMODEL             _grd_model = GRDMODEL::DEF_GRDMODEL;                //glfeng
		bool                 _tropo_grad    = false;    //glfeng
		bool                 _lite_turboedit= false;    //glfeng     true:MW+GF+Gap, false:read log file

		/// The following is public for all lsq
		std::string	_crt_rec;					///< current site of proc
		std::string	_crt_sat;					///< current sat  of proc

		t_gallobs*       _gall_obs     = nullptr;			    ///< all observ info
		t_gallobj*       _gall_obj     = nullptr;			    ///< all obj info include pcv and intial crd
		t_gallrecover*   _gall_rcv     = nullptr;			    ///< all recover info
		t_gallbias*      _gall_bias    = nullptr;

		t_gallproc*      _gall_proc    = nullptr;			    ///< all     the data used in proc
		t_gallics_nav*   _gall_gns     = nullptr;			    ///< all gns ics data used in proc 
		t_gallics_leo*   _gall_leo     = nullptr;			    ///< all leo ics data used in proc
		t_gpoleut1*      _allerp = nullptr;
		t_gallnav*       _gall_nav     = nullptr;               ///< used in kin ppp
		t_gallprec*      _gall_prec    = nullptr;
		t_gallpreproc*   _gall_preproc = nullptr;               //SLR parts
		t_gallrb*        _gall_rb      = nullptr;               //SLR parts

		t_gpoleut1*      _gdata_erp    = nullptr;			    ///< all poleut1 data for init ERP parameters
		t_gion*          _gdata_ion    = nullptr;
		t_gambcon*       _gdata_ambcon = nullptr;
		t_gslrf*         _gdata_slrf   = nullptr;               //SLR parts
		t_gallsatinfo*   _gdata_satinfo = nullptr;              ///< sat information for glonass fid

		t_glsq*	                    _lsq        = nullptr;		///< least square estimator
		shared_ptr<t_gcycleslip>    _slip12     = nullptr;      ///< turboedit class for IF
		shared_ptr<t_gcycleslip>    _slip13     = nullptr;	    ///< turboedit class for 13 if
		shared_ptr<t_gcycleslip>    _slip14     = nullptr;
		shared_ptr<t_gcycleslip>    _slip15     = nullptr;
		shared_ptr<t_gbasemodel> _base_model    = nullptr;
		shared_ptr<t_gbiasmodel> _bias_model	= nullptr;
		t_gprecisemodel*         _bias_pcemodel = nullptr;
		t_gslrmodel*             _slr_model     = nullptr; //SLR parts

		// record now process data
		t_glsqEquationMatrix _crt_equ;
		map<string, map<t_gtime, t_gtriple> > _kin_xyz;   ///< xyz coord reference at each epoch, for update parameters

		 
		map<int, string>   _slr_list;
		map<string, vector<string>>  _sat_sta;
		bool               _slr = false;   ///< use SLR observations or not

		set<string> _mode;  ///< POD mode(S/D/K/SD)

		double             _recclk_threshold = 0.0;

		int _num_threads = 1;

	protected:
		double _maxres_norm = 0.0;                    //added by xiongyun
		t_glog* _clk_log      = nullptr;
		t_giof* _satclkfile   = nullptr;
		t_giof* _rtnetclkfile = nullptr;
		t_giof* _satclk13file = nullptr;
		bool   _IRC = false;                          ///< estimate integer recovery clock [true=1, false=0] (add by xiongyun)
		bool   _carrier_range = false;                ///< use carrier range data [true=1, false=0] (add by xiongyun)
		bool   _use_carrier_range = false;            // added by jqwu
		bool   _realtime;
		bool   _ultrasp3 = false;                     //added by xiongyun(whether use sp3 file in realtime)
		map<GSYS, double>   _maxres_C;                //added by xiongyun
		map<GSYS, double>   _maxres_L;                //added by xiongyun
		set<GSYS> sys_support;                        //added by xiongyun
		map<string, vector<string>> _outlier_preres;  //added by xiongyun
		map<string, t_gallrecover*> _siteres;         //added by xiongyun for site res
		vector<t_gcoder*>           _gcoder;          //added by liyuhao
		bool _isCreatLog = false;  ///< whether the process batch creat a log file, only when the glog is nullptr
		bool _ambfix_cmd = false;  ///< ambfix  or not, set in cmd pars
		bool _ambfix_xml = false;
		bool _useIon_cmd = false;  ///< use ion or not, set in cmd pars


		typedef tuple<int, string, string, t_gtime, t_gtime> t_amb_index;
		map<tuple<t_amb_index, t_amb_index, t_amb_index, t_amb_index>, double> _Dd_amb_L12_used; ///< L1 ambcon while fix with L1 and L2, used in checking L1 and L3 ambcon.
		map<tuple<t_amb_index, t_amb_index, t_amb_index, t_amb_index>, double> _Dd_amb_L13_used; ///< L1 ambcon while fix with L1 and L3, used in checking L1 and L2 ambcon.

	protected:
		/// The lsq process module for override
		virtual bool _init_crs_erp_pars(t_glsq* lsq);
		virtual bool _init_rec_trop_pars(t_glsq* lsq);
		virtual bool _init_gns_vion_pars(t_glsq* lsq);
		virtual bool _init_gns_sion_pars(t_glsq* lsq);
		virtual bool _init_rec_crd_pars(t_glsq* lsq, bool isPPP = true);
		virtual bool _init_rec_clk_pars(t_glsq* lsq);
		virtual bool _init_rec_clk_pars_irc(t_glsq* lsq);//add by xiongyun
		virtual bool _init_sat_clk_pars(t_glsq* lsq);
		virtual bool _init_sat_ifcb_pars(t_glsq* lsq);//add by xiongyun
		virtual bool _init_gns_ics_pars(t_glsq* lsq);
		virtual bool _init_leo_ics_pars(t_glsq* lsq);
		virtual bool _init_sys_bia_pars(t_glsq* lsq);
		virtual bool _init_sys_rec_ifb_pars(t_glsq* lsq);
		virtual bool _init_sys_sat_ifb_pars(t_glsq* lsq);
		virtual bool _init_sys_ifb_pars(t_glsq* lsq);
		virtual bool _init_IF_amb_pars(t_glsq* lsq);
		virtual bool _init_sys_bia_pars_pce(t_glsq* lsq);//add by xiongyun
		void _set_isb_par_setting(t_glsq* lsq, t_gpar& par, SYSBIASMODEL model, double value, double apriori);

		// for triple-frequency IF PPP
		virtual bool _init_rec_clk13_pars(t_glsq* lsq);
		virtual bool _init_IF_amb13_pars(t_glsq* lsq);

		// -- for whole
		virtual bool _initLsqProdData(t_gallprod* data);
		virtual bool _initLsqProcPars(t_glsq* lsq);
		virtual bool _initLsqProcData(t_gallproc* data);

		// for SLR observations
		bool _InitSLRFData(t_gslrf* _gslrf);
		bool _initSLRProcPars(t_glsq* lsq);
		bool _init_slr_pars(t_glsq* lsq);
		bool _init_slr_rb(t_glsq* lsq);
		bool _init_crs_erp_pars_slr(t_glsq* lsq);

		//virtual bool _update_lsq_par(const t_gtime& epoch, vector<t_gsatdata>& obsdata);
		bool _select_rec_obs(const string& rec, const t_gtime& epoch, vector<t_gsatdata>& all_obs, vector<t_gsatdata>& rec_obs);
		bool _apply_carrier_range(const string& rec, const string& sat, const t_gtime& epoch, t_gsatdata& obs);
		bool _init_sat_clk_IF(vector<t_gsatdata>& crt_obs, ColumnVector& l, t_glsq* lsq);
		bool _init_sat_clk_ALL(vector<t_gsatdata>& crt_obs, ColumnVector& l, t_glsq* lsq);
		bool _init_sat_clk();

		// jdhuang
		bool _init_xml_settings(t_gsetbase* set);
		bool _init_lsq_settings(t_gsetbase* set, t_gallproc* data, t_glog* log);

		/**
		* @brief extract clkfile from recover data
		* @return success or fail
		*/
		bool _extract_clkfile(t_gallrecover& recover_data, t_gallprec::clk_type type);
		bool _extract_clkfile(t_gallrecover& recover_data, t_gallprec::clk_type type, bool isTriple);
		bool _extract_resfile(t_gallrecover& recover_data, string resfile="");
		bool _extract_rbfile(t_gallrecover& recover_data);
		bool _extract_stafile(t_gallrecover& recover_data);
		bool _extract_ambupd(t_gallambupd* allambupd);
		bool _get_rb_std(t_gallrb& rb_data);
		bool _extract_eopfile(t_gallrecover& recover_data);
		bool _extract_polefile(t_gallrecover& recover_data);
		bool _extract_xmlcrd(t_gallrecover& recover_data);

		/**
		* @brief read kin file
		* @return void
		*/
		virtual void _read_kinfile();

		// data

		string _recclk_file;	///< outputs recclk file name

		bool _myLog = false;
		bool _create_log(t_glog* log);
		bool _delete_log(t_glog* log);
		bool _add_sion_pars_value(const t_gtime& beg, const t_gtime& end, t_glsq* lsq, t_gion* _gion);
		void _extract_ion_file(t_gallrecover& recover_data);
		bool _get_site_crd(string site, t_gtriple& site_crd); //added by xiongyun to get site crd

		/**
		* @brief add single combine dd ambcon for raw all processs.
		*   L_NL(1) - LNL(2) - L_NL(3) + L_NL(4)
		*   L_WL(1) - LWL(2) - L_WL(3) + L_WL(4)
		* @param[in] crt_epoch   curruent epoch.
		* @param[in] lsq         lsq class.
		* @param[in] ambcon      ambcon data class.
		* @param[in] isL12       add for L1 and L2, L1 and L3.
		* @return	  true - successfully;
		*            false - failure
		*/		
		bool _pod_add_Dd_ambcon(t_glsq* lsq, t_gambcon* ambcon);

		/**
		* @brief add single combine dd ambcon for raw all processs.
		*   L_NL(1) - LNL(2) - L_NL(3) + L_NL(4)
		*   L_WL(1) - LWL(2) - L_WL(3) + L_WL(4)
		* @param[in] crt_epoch   curruent epoch.
		* @param[in] lsq         lsq class.
		* @param[in] ambcon      ambcon data class.
		* @param[in] isL12       add for L1 and L2, L1 and L3.
		* @return	  true - successfully;
		*            false - failure
		*/	
		bool _pod_add_RAW_CB_Dd_ambcon(const t_gtime& crt_epoch, t_glsq* lsq, t_gambcon* ambcon, const bool& isL12 = false);

		//single difference
		/**
		* @brief add single combine dd ambcon for raw all processs.
		*   L_NL(1) - LNL(2)
		*   L_WL(1) - LWL(2)
		* @param[in] lsq         lsq class.
		* @param[in] Dd_amb      ambcon data class.
		* @return	  true - successfully;
		*            false - failure
		*/
		bool _pod_add_one_Sd_ambcon(t_glsq* lsq, const t_gamb_Dd* Dd_amb);

		virtual bool _processOneEpoch(const t_gtime& crt_epoch, std::vector<t_gsatdata>& crt_obs);
		virtual bool _processOneRec(const t_gtime& crt_epoch, const std::string& crt_rec, std::vector<t_gsatdata>& crt_obs);
		virtual bool _processOneRecOneSat(const t_gtime& crt_epoch, const std::string& rec, const std::string& sat, t_gsatdata& crt_obs);


		//pre-fix residuals
		bool _prefix_residuals(vector<t_gsatdata>& vec_obs);
		bool _is_prefix_res;
	};


}
#endif