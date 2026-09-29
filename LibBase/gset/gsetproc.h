/**
*
* @verbatim
	 (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  @endverbatim
*
* @file		gsetproc.h
* @brief	implements process setting class
* @author   Jan Dousa
* @version	1.0.0
* @date		2012-10-23
*
*/

#ifndef GSETPROC_H
#define GSETPROC_H

#define XMLKEY_PROC "process" ///< The defination of process node in XML
#define XMLKEY_ROM "read_ofile_mode"

#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gutils/gtypeconv.h"
#include "gutils/gobs.h"
#include "gset/gsetbase.h"

using namespace std;
using namespace pugi;

namespace gnut
{

	enum CONSTRPAR { EST, FIX, KIN, SIMU_KIN, CONSTRPAR_UNDEF};

	enum GRDMPFUNC { DEF_GRDMPFUNC, TILTING, CHEN_HERRING, BAR_SEVER };
	enum ZTDMPFUNC { DEF_ZTDMPFUNC, COSZ, GMF,NO_MF };//lvhb added "NO_MF" for npp, 20200917
	enum IONMPFUNC { DEF_IONMPFUNC, ICOSZ, QFAC, NONE };
	enum OBSWEIGHT { DEF_OBSWEIGHT, EQUAL, SINEL, SINEL2, SINEL4, CODPHA, MLTPTH, PARTELE,SNR };
	enum TROPMODEL { DEF_TROPMODEL, SAASTAMOINEN, DAVIS, HOPFIELD, MOPS, GPTW, GPT2W, GAL27, GALTROPO27, EXTERN };
	enum ZTDMODEL  { DEF_ZTDMODEL, PWC, STO };
	enum RESIDTYPE { DEF_RESIDTYPE, RES_ORIG, RES_NORM, RES_ALL };
	enum OBSCOMBIN { DEF_OBSCOMBIN, IONO_FREE, RAW_SINGLE, RAW_DOUBLE, RAW_ALL, MW_COMBIN, EWL_COMBIN, WL_COMBIN, RAW_MIX, IF_P1};
	enum PHASEBIAS { DEF_POSTPRD, RTCM_PRD, CNES_PRD, ESTI_PRD };
	enum ATTITUDES { DEF_YAWMODEL, YAW_NOMI, YAW_RTCM };
	enum CBIASCHAR { DEF_CBIASCHAR, ORIG, CHAR2, CHAR3 };
	enum SYSBIASMODEL { AUTO_CON, AUTO_RWK, ISB_CON, ISB_RWK, AUTO_WHIT };   //glfeng
	enum GRDMODEL { DEF_GRDMODEL, GRD_PWC, GRD_STO };             //glfeng
	enum BASEPOS { RINEXH, SPP, CFILE };

	/** @brief  enum for lsq model */
	enum LSQMODE
	{
		EPO,	///< EPO model:solve every epoch
		LSQ,    ///< LSQ model:Overall solution
		LSQMODE_UNDEF
	};

	// 周跳探测方法
	enum SLIPMODEL
	{
		DEF_DETECT_MODEL, ///< 
		TURBO_EDIT,
		SLIPMODEL_UNDEF
	};

	// 电离层模型
	enum IONMODEL
	{
		VION,
		SION,
		DEF_ION,
		IONMODEL_UNDEF
	};

	enum class IFCB_MODEL
	{
		EST,
		COR,
		DEF,
		IFCB_MODEL_UNDEF
	};

	enum class IFB_MODEL
	{
		EST_REC_IFB,
		EST_SAT_IFB,
		EST_IFB,
		COR_REC_IFB,
		COR_SAT_IFB,
		COR_IFB,
		NONE,
		IFB_MODEL_UNDEF
	};
	enum class NPP_MODEL
	{
		NO,
		VRS,
		URTK,
		PPP_RTK,
		NPP_MODEL_UNDEF
	};
	enum class DCB_MODEL
	{
		CODE,
		CAS,
		NONE,
		DCB_MODEL_UNDEF
	};

	enum class FREQ_MODEL
	{
		DOUBLE, // double freq
		TRIPLE, // triple freq
		QUADRUPLE, // 
		QUINTUPLE,
		FREQ_MODEL_UNDEF
	};

	enum class IONO_ORDER
	{
		FIRST_ORDER,
		SECOND_ORDER,
		THIRD_ORDER,
		NONE,
		IONO_ORDER_UNDEF
	};

	// Hwang Shih added
	enum class READ_OFILE_MODE
	{
		REALTIME, // Processing while reading O files
		READALL,  // Processing after reading all of the O files which is defaulted
		READ_OFILE_MODE_UNDEF
	};

	// The Define for LC PC
	class BASE_LIBRARY_EXPORT t_gobscombtype
	{

	public:
		t_gobscombtype();
		t_gobscombtype(const t_gobscombtype& other);
		t_gobscombtype(string obscombtype);
		t_gobscombtype(t_gobs obstype, OBSCOMBIN combtype);
		t_gobscombtype(GOBSTYPE t, GOBSBAND b, OBSCOMBIN obscomb);

		string convert2str() const;

		bool operator ==(const t_gobscombtype& g) const;
		bool operator <(const t_gobscombtype& g) const;

		bool is_phase() const;
		bool is_code() const;
		bool is_SLR() const;

		GOBSBAND getBand() { return _obs_band; };// Added by lvhb 20200428

	protected:

		GOBSTYPE  _obs_type;
		GOBSBAND  _obs_band;
		OBSCOMBIN _obs_combine;

	};


	/// The class for setting of process modelue in XML file
	class BASE_LIBRARY_EXPORT t_gsetproc : public virtual t_gsetbase
	{
	public:
		/// constructor
		t_gsetproc();
		/// destructor
		~t_gsetproc();

		/// settings check
		void check(); 
		/// settings help
		void help();                                  

		/**
		 * @brief         get the decimals for sampling interval (for high-rate) in process
		 * @return bool : decimals for sampling interval (for high-rate) in process
		 */
		bool   isleo();
		bool   tropo();
		bool   iono();
		bool   tropo_slant();
		bool   tropo_grad();
		bool   phase();
		bool   doppler();
		bool   pos_kin();
		bool   use_eclipsed();
		bool   auto_band();
		bool   sd_sat();
		bool   IRC();
		bool   carrier_range(); //xy
		bool   apply_carrier_range(); // jqwu 
		bool   pre_fix_residuals();

		bool   ambfix();

		bool   realtime();   // added by zhShen
		bool   ultrasp3();   // added by xiongyun
		int    dynamics();
		int    frequency();
		int	   num_threads(); 
		double sig_init_ztd();
		double sig_init_vion();
		double sig_init_grd();
		double sig_init_crd();
		double sig_init_vel();
		double sig_init_amb();
		double sig_init_glo();
		double sig_init_gal();
		double sig_init_bds();
		double sig_init_qzs();

		double sig_ref_clk();
		double minimum_elev();
		double minimum_elev_leo(); // add by zhangwei for leo processing
		double max_res_norm();
		int read_ofile_step();
		int minsat();

		string ref_clk();
		string ref_rec();
		LSQMODE lsq_mode();
		SLIPMODEL slip_model();
		IONMODEL ion_model();

		TROPMODEL tropo_model();
		ZTDMODEL  ztd_model(double* dt = 0);
		ZTDMPFUNC tropo_mf();
		IONMPFUNC iono_mf();
		GRDMPFUNC grad_mf();
		CONSTRPAR crd_est();
		OBSWEIGHT weighting();
		RESIDTYPE residuals();
		OBSCOMBIN obs_combin();
		ATTITUDES attitudes();
		string    trop_corr();
		CBIASCHAR cbiaschar();

		SYSBIASMODEL sysbias_model();        //glfeng
		SYSBIASMODEL ifbbias_model();        //LX
		GRDMODEL grd_model(double* dt = 0);  //glfeng

		GRDMPFUNC str2grdmpfunc(string mf);
		ZTDMPFUNC str2ztdmpfunc(string mf);
		IONMPFUNC str2ionmpfunc(string mf);
		OBSWEIGHT str2obsweight(string wg);
		TROPMODEL str2tropmodel(string tm);
		ZTDMODEL   str2ztdmodel(string ztd);
		RESIDTYPE str2residtype(string rs);
		OBSCOMBIN str2obscombin(string oc);
		ATTITUDES str2attitudes(string at);
		CBIASCHAR str2cbiaschar(string cb);
		SYSBIASMODEL str2sysbiasmodel(string sys);  //glfeng
		GRDMODEL str2grdmodel(string grd);          //glfeng

		string grdmpfunc2str(GRDMPFUNC MF);
		string ztdmpfunc2str(ZTDMPFUNC MF);
		string ionmpfunc2str(IONMPFUNC MF);
		string obsweight2str(OBSWEIGHT WG);
		string tropmodel2str(TROPMODEL TM);
		string ztdmodel2str(ZTDMODEL ZTD);
		string residtype2str(RESIDTYPE RS);
		string obscombin2str(OBSCOMBIN OC);
		string cbiaschar2str(CBIASCHAR CB);

		bool simulation();

		IFCB_MODEL                          ifcb_model(); // jdhuang : fcb
		FREQ_MODEL                          freq_model(); // jdhuang : used in raw all
		IFB_MODEL                           ifb_model();

		IONO_ORDER                          iono_order(); // jdhuang : used for ion order
		//rtk
		BASEPOS basepos();		///< obtain base coordinate
		BASEPOS str2basepos(string str);
		NPP_MODEL npp_model();
		READ_OFILE_MODE						read_ofile_mode();// Hwang Shih : used for reading O files
		// smooth range
		string range_smooth_mode(int* smt_windows = 0);
		// correct bds code bias
		bool bds_code_bias_correction();

	protected:

		bool   _phase;                                ///< use phase data [true=1, false=0]
		bool   _tropo;                                ///< estimate troposphere [true=1, false=0]   
		bool   _iono;                                 ///< estimate ionosphere [true=1, false=0]   
		bool   _tropo_grad;                           ///< estimate troposphere gradinet [true=1, false=0]
		bool   _tropo_slant;                          ///< estimate tropo slant delays
		TROPMODEL _tropo_model;                       ///< tropospheric model [SAASTAMOINEN, DAVIS, HOPFIELD, ...]
		ZTDMODEL  _ztd_model;						  ///< ztdmodle [PWC STO]
		ZTDMPFUNC _tropo_mf;                          ///< tropo mapping function [COSZ, NIELL, GMF, VMF1, ... ]
		IONMPFUNC _iono_mf;                           ///< iono mapping function [COSZ, QFAC, NONE, ... ]   
		GRDMPFUNC _grad_mf;                           ///< grad mapping function [tilt, CHH]
		OBSWEIGHT _obs_weight;                        ///< observations weighting
		RESIDTYPE _res_type;                          ///< types of residuals
		OBSCOMBIN _obs_combin;                        ///< observation combination 
		ATTITUDES _attitudes;						  ///< satellite attitude model
		CBIASCHAR _cbiaschar;						  ///< forcing code bias signals

		double _sig_init_ztd;                         ///< accuracy of initial zenith path delays [m]
		double _sig_init_vion;                        ///< accuracy of initial vertical iono path delays [m]
		double _sig_init_grd;                         ///< accuracy of initial tropo gradients [m]   
		double _sig_init_crd;                         ///< accuracy of initial coordinates [m] 
		double _sig_init_vel;                         ///< accuracy of initial velocities [m/s]
		double _sig_init_amb;                         ///< accuracy of initial ambiguity [m]
		double _sig_init_glo;                         ///< accuracy of initial GLONASS system time difference
		double _sig_init_gal;                         ///< accuracy of initial Galileo system time difference
		double _sig_init_bds;                         ///< accuracy of initial BeiDou system time difference
		double _sig_init_qzs;                         ///< accuracy of initial QZSS system time difference
		double _minimum_elev;                         ///< elevation angle cut-off [degree]
		double _max_res_norm;                         ///< normalized residuals threshold
		string _crd_est;                              ///< FIX or estimate CRD
		bool   _pos_kin;                              ///< static/kinematic receiver (true == kinematic)
		bool   _use_eclipsed;                         ///< used eclipsed satellites
		bool   _auto_band;                            ///< automatic band selection (mainly for Anubis purpose)
		int   _dynamics;
		int   _frequency;
		bool  _realtime;
		bool  _simulation;                            ///< judge whether simulation data

		int _minsat;                                 ///<  minimum satellite number
		BASEPOS _basepos;
		bool   _sd_sat;                              ///< single differented between sat and sat_ref
		NPP_MODEL _nppmodel;
	private:
	};

} // namespace

#endif
