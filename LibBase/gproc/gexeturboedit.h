/**
*
* @verbatim
	History
	 -1.0 Guolong Feng  2020-03-20 created
  @endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file		    gexeturboedit.h
* @brief		main control class for turboedit
* @author       Guolong Feng, Wuhan University
* @date		    2020-03-20
*
*/


#ifndef GEXETURBOEDIT_H
#define GEXETURBOEDIT_H

#include "gutils/ExportBase.h"
#include "gutils/gfileconv.h"
#include "gutils/gturboedit.h"

#include "gset/gsetgen.h"
#include "gset/gsetgnss.h"
#include "gset/gsetrec.h"
#include "gset/gsetbase.h"
#include "gset/gsetinp.h"
#include "gset/gsetturboedit.h"

#include "gio/glog.h"
#include "gio/gfile.h"

#include "gall/gallambflag.h"
#include "gall/gallobs.h"
#include "gall/gallnav.h"
#include "gall/gallproc.h"

#include "gproc/gsppflt.h"
#include "gproc/gqualitycontrol.h"

#include "gcoders/satparameters.h"  // For LEO GLO
#include "gcoders/ambflag.h"
#include <algorithm>
#include <stdint.h>
#include <string>
#if  defined _WIN32 || defined _WIN64
#include <io.h>
#include <direct.h> 
#define ACCESS(fileName,accessMode) _access(fileName,accessMode)
#define MKDIR(path) _mkdir(path)
#else
#include <unistd.h>
#include <sys/stat.h>
#define ACCESS(fileName,accessMode) access(fileName,accessMode)
#define MKDIR(path) mkdir(path,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif

using namespace gnut;

namespace great
{
	// ---------------------
	enum DATAFLAG {

		GOOD = 0,
		GFJUMP = 2,   // ambiguity inserted due to jump in ionosphere-lane
		MWJUMP = 3,   // ambiguity inserted due to jump in wide-lane		
		GAP = 5,   // ambiguity inserted after long gap
		LLI = 8,   // ambiguity inserted accoording to LLI in rinex-file

		PCBAD = 16,  // not used because of range is bad

		NOBRD = 22,  // not used because not broadcast ephemeris is available
		GFBAD = 24,  // not used because ionosphere - lane is bad
		MWBAD = 25,  // not used because wide - lane is bad

		ASHADOW = 27,  // not used because II / IIA after shadow
		SHRT = 28,  // not used because data piece is too short
		LOWELE = 29,  // not used because of lower elevation
		NO4 = 30,  // not used because the data set is not complete
		NODATA = 31,  // No data in rinex-file
	};


	/** @brief Specified sat data for post-processing turboedit. */
	class BASE_LIBRARY_EXPORT t_gturbo_post
	{
	public:
		/** @brief default constructor. */
		t_gturbo_post() {};
		/** @brief default destructor. */
		virtual ~t_gturbo_post() {};
		
		map<int, int>                iflag;         ///< epoch, flag_index
		map<int, double>             elev;          ///< epoch, elevation
		map<int, map<int, double>>   obs;           ///< epoch, index, obs
		int                          last_epo;      ///< last epoch
	};

	/** @brief Smooth information for real-time turboedit. */
	class BASE_LIBRARY_EXPORT t_gturbo_rt
	{
	public:
		/** @brief default constructor. */
		t_gturbo_rt() {};
		/** @brief default destructor. */
		virtual ~t_gturbo_rt() {};
		map<GOBS, bool>       mwslip;
		map<GOBS, bool>       gfslip;
		map<GOBS, t_gtime>    last_time;
		map<GOBS, int>        epo_num;                      ///< Store the number of consecutive epochs
		map<GOBS, double>     GF;                           ///< store GF value
		map<GOBS, pair<double, double> >   smooth_MW;       ///< store MW smooth value, std
		map<GOBS, vector<double>>          origin_MWs;      ///< sort according _gap
	};

	/**
	* @brief main control class for caculate clk
	*/
	class BASE_LIBRARY_EXPORT t_gexeturboedit
	{
	public:
		
		t_gexeturboedit(t_gsetbase* set, t_glog* glog, t_gallobs* gobs, t_gallnav* gnav);
		t_gexeturboedit(t_gsetbase* set, t_glog* glog, t_gallproc* gdata);
		virtual ~t_gexeturboedit();

		bool ProcessBatch(string site, const t_gtime& begT, const t_gtime& endT, double sampling);
		void setOBJ(t_gallobj*  gobj) { _allobj = gobj; };
		void setOBS(t_gallobs*  gobs) { _allobs = gobs; };
		void setGLOFreqID(t_gallsatinfo* gsatpars);  // For simulation data
		void set_site_crd(string site, t_gtriple crd);  // set site crd added by xiongyun
		void setSingleEpoData(vector<t_gsatdata>* input_data);
		bool checkSlipWithDoppler(string site, const t_gtime& begT, const t_gtime& endT, double sampling); // For npp Client

		shared_ptr<t_gcycleslip> getSiteCycleSlip(string site, int idx) { return _cycleslip[site][idx]; };
	protected:		
		t_gsetbase*                    _set;              ///< set
		t_glog*                        _log;              ///< log file
		///< Global Settings	       
		int                            _epo_num;
		int                            _interval;
		t_gtime                        _beg_time;
		t_gtime                        _end_time;
		set<string>                    _sys_list;
		set<string>                    _sat_list;
		set<string>                    _rec_list;
		map<string, t_gtriple>         _rec_crds;        ///< site coordinates
		

		///< input data
		t_gallobs*       _allobs;
		t_gallnav*       _allnav;
        t_gallobj*       _allobj;                         ///< For sppflt

		///< Post-processing Settings
		map<int, map<GSYS, pair<GOBSBAND, GOBSBAND>>>        _post_all_band;    ///< set post-processing bands
		map<GSYS, pair<GOBSBAND, GOBSBAND>>                  _post_crt_band;    ///< set current post-processing bands
		map<GSYS, map<GOBSBAND, FREQ_SEQ> >                  _freq_index;
		map<GSYS, vector<GOBS>>                              _post_phase_obs;   ///< set phase observation type
		map<string, t_gturbo_post>                           _post_data;
		map<int, t_gtime>                                    _post_time;
		map<string, vector<shared_ptr<t_gcycleslip>>>        _cycleslip;
		t_gallambflag*                                       _allambflag;       ///< save_ambflag
		map<string, t_gsppflt*>                              _sppflt;		    ///< sppflt, provide init crd
		map<string, double>                                  _scaling_factor;   ///< For check large cycle, post-processing
		map<string, int>                                     _glonass_freq_id;  ///< glonass frequency id
		

		///< Real-time processing mode
		int                                                 _smooth_window;
		int                                                 _frequency;
		map<string, map<string, map<GOBS, t_gturbo_rt>>>    _rt_data;
		map<string, map<string, shared_ptr<t_gobsgnss>>>    _epoDataPre;
		t_goutliers_process                                 _outliers_proc;
		OBSCOMBIN                                           _observ;
		t_giof*                                             _debug_turbo = nullptr;

		// for npp
		vector<shared_ptr<t_gobsgnss>>                      _inputEpoData;
		vector<t_gsatdata>*                                 _satdata_npp = nullptr;
		map<string, map<string, int>>                       _map_nppdata_idx; // site sat idx
		map<string, map<string, shared_ptr<t_gobsgnss>>>    _epoDataPre_npp;

		///< currenct processing site
		t_gtime		                         _crt_time;       
		int                                  _crt_epo;
		string                               _crt_rec;                         ///< For Log file


		///< turboedit settings
		bool   _realtime;           ///< 判断是否实时模式
		bool   _liteMode;           ///< lite mode or not, same with real time mode
		bool   _amb_output;         ///< 判断是否输出到log文件
		bool   _ephemeris;          ///< whether use ephemeris, default true
		bool   _simulation;         ///< whether simulation data, default False;
		bool   _site_kinematic;     ///< kinematic station is True, default False，判断是否动态测站
		bool   _check_pc;           ///< check PC, default 
		bool   _check_mw;           ///< check MW conbination (TURBOEDIT), default True
		bool   _check_gf;           ///< check GF conbination (TURBOEDIT), default True
		bool   _check_sf;           ///< check single frequency data, default is not.
		bool   _check_gap;          ///< check data gap, ambiguity is inserted if data missing longer
		bool   _check_short;        ///< check short arc, data piece shorter considered as short piece
		//bool   _check_with_doppler; ///< check single frequency phase with doppler observations

		//double _delta_PL_limit;     ///< epoch-difference range - epoch-difference phase
		//double _delta_LD_limit;     ///< epoch-difference phase - integration of doppler observation
        double _pc_limit;           ///< default 250m		
		double _mw_limit;           ///< widelane limit, default 4 cycles
		double _gf_limit;           ///< geometry_free_limit, default 1 cycles       ///////// in real-time/lite mode , unit is meters
		double _gf_rms_limit;       ///< default 2 cycles
		double _sf_limit;           ///< default ????
		double _min_elev;           ///< cutoff elevation, default is 0
		int    _gap_limit;          ///< epoch number, ambiguity is inserted if data missing longer, default 20 epoch
		int    _short_limit;        ///< epoch number, data piece shorter considered as short piece, default 10 epoch
		bool   _check_statistics;   ///< For Network Processing, PPP should Close; Whether consider following threshold, default not=false
		double _min_percent;        ///< The percentage of epochs(>4 sats) in all useful epochs
		int    _min_mean_nprn;      ///< Minimum mean sat numbers (sum of obs_number(>4 sats) / epoch_number(> 4 sats) )
		int    _max_mean_namb;      ///< Ambiguities Number / sat_number 
	
		///< Functions
		/**
		* @brief  output control informations
		*/
		void _outputControlInfo();

		/**
		* @brief  initial function of processing turboedit
		*/
		void _initCommonProcess();

		/**
		* @brief  initial function of post-processing turboedit
		*/
		void _initRealTimeProcess();

		/**
		* @brief  main function of real-time processing turboedit
		*/
		bool _realtimeProcess(const t_gtime& begT, const t_gtime& endT);

		/**
		* @brief  Turboedit Check[MW+GF] of real-time processing turboedit
		*/
		void _realTimeMWGF(shared_ptr<t_gobsgnss> oneObs, GOBSBAND b1, GOBSBAND b2);
		
        void _realTimeMWGF_test(shared_ptr<t_gobsgnss> oneObs, GOBSBAND b1, GOBSBAND b2); //xiongyun test

		void _realTimeMWGF_new(shared_ptr<t_gobsgnss> oneObs, GOBSBAND b1, GOBSBAND b2);

		void _checkSlipWithDoppler(shared_ptr<t_gobsgnss> ObsPre, shared_ptr<t_gobsgnss> Obs);

		/**
		* @brief  initial function of post-processing turboedit
		*/
		void _initPostProcess();

		/**
		* @brief  main function of post-processing turboedit
		*/
		bool _postProcess();

		/**
		* @brief  initial _post_data/_post_time for each station
		*/
		void _initPostTurboInfo();

		/**
	    * @brief  init SPPflt for kinematic stations
	    */
		bool _siteSPP(const string& site, const t_gtime& now);
	
		/**
		* @brief  basic check is performed. Receiver clock bias is computed
        *         and corrected to DB%ti() if BRD is used.
		*/
		void _postBasicCheck();

		/**
		* @brief  check observations vaid or not, as well as combine observation
		* @param[in] obsdata   single observation
		* @return              valid true; not false
		*/
		bool _postPrepareObs(t_gobsgnss& obsdata);

		/**
		* @brief  compute distance between satellite and station, satellite
		*         elevation, satellite clock bias,  using the specified brd
        *         record, or the nearest record.
		* @param[in] sat        sat name
		* @param[in] distance   distance between satellite and station, meters 
		* @param[in] elev       elevation, degree
		* @param[in] sat_clk    satellite clock bias, meters
		* @return               Ok is true; not is false
		*/
		bool _getSatInfo(string sat, double& distance, double& elev, double& sat_clk);
		
		/**
		* @brief  Combination of observations
		* @param[in] sat   sat name
		*/
		void _postCheckLargeSlips(string sat);

		/**
		* @brief  Check range residuals at the epoch for outliers
		*/
		void _postCheckRange();
		
		/**
		* @brief  remove data piece with data no longer than len_short and 
        *         insert ambiguities if data gap longer than len_gap
		* @param[in] length_short  length of short piece, seconds
		* @param[in] sat           sat need to be remove short
		* @return                  false: removeShort work have been done; true: need further removement        
		*/
		bool _postCheckGapShort(double length_short, string sat);

		/**
		* @brief  single-station GPS data editing based on the turboedit
		*         algorithm by J. Blewitt. 
		* @param[in] sat       sat name
		*/
		void _postCheckMWGF(string sat);
	
		/**
		* @brief   edit M-W combined observations (Nw )
		* @param[in] sat       sat name
		*/		
		void _postEditWidelane(string sat);

		/**
		* @brief  check jumps in ionospheric observaitons
		* @param[in] sat       sat name
		*/
		void _postCheckIonosphere(string sat);

		/**
		* @brief  check for jump using polynomial-fitting
		* @param[in]     strflag        name of time series
		* @param[in]     ti             time series (seconds)
		* @param[in]     values         values corresponding to time series
		* @param[in&out] flag           data flag
		* @param[out]    res            vector with residuals (diff values)
		* @param[out]    rms            unit middle error of residuals
		*/
		bool _checkForJump(string strflag, const vector<double>& ti, const vector<double>& values,
			vector<int>& flag, vector<double>& res, double& rms);

		/**
		* @brief  approximation of a set of discrete functional values
		*         by a polynomial and determining the integral of the
		*         polynomial using differences of consecutive values
		* @param[in]     degree        degree of polynomial
		* @param[in]     X             vector with independent arguments
		* @param[in]     Y             vector with functional values
		* @param[in]     L             l(i)=0: value # i is used for approx.
		*                                  =1: value # i is not used for approx.
		*                          =2,3,4,...: value # i is not used at all                        
		* @param[in&out] C             vector with coefficients of polynomial
		* @param[out]    V             vector with residuals (diff approx-y)
		* @param[out]    ipt           used data flagged as its position
		* @param[out]    rms           unit middle error of residuals
		*/
		bool _polynomialFit(int degree, const vector<double>& X, const vector<double>& Y, const vector<int>& L, 
			vector<double>& C, vector<double>& V, vector<int>& ipt, double& rms);

		/**
		* @brief   write the report in log files
		* @param[out] decision   decide whether remove short again
		*/	
		void _postPrintReport(string& decision);

		/**
		* @brief   add ambflag information of a site to _all_ambflag (ambflag decoder)
		* @return     success or not
		*/
		bool _postSetAmbflag();

		/**
		* @brief  set ambflag data to _all_ambflag
		* @param[in]  sat_epo_beg_end       ambflag infomation[sat/valid beg epo/valid end epo]
		* @return     success or not
		*/
		bool _postSetAmbflagData(map<string, map<int, int>>& sat_epo_beg_end);

		/**
		* @brief  get flag description(written in ambflag files)
		* @param[in]  sat         sat name
		* @param[in]  str3        AMB/DEL/BAD
		* @param[in]  epo         current epoch
		* @return                 descriptions
		*/
		string _strFlagDescription(string sat, string str3, int epo);

		/**
		* @brief  Encoder ambflag files
		* @param[in]  logsuffix   suffix of ambflag(e.g., log) file
		* @return                 success or not
		*/
		bool _postEncoderProduct(string logsuffix);

		///< Tool Functions
		/**
		* @brief  set data flag to the kind what the string described.
		* @param[in] iflag  original iflag
		* @param[in] type   string description
		* @return           new iflag
		*/
		int  _setFlag(int iflag, string type);

		/**
		* @brief  find the next epoch with the flag means `flag_type`
		* @param[in] istart     start epoch
		* @param[in] istop      end epoch. Search is done from istart to istop
		* @param[in] flg        flags
		* @param[in] flag_type  string flag
		* @param[in] direction  "+" means forward, while "-" means backward
		* @return               return found epoch, if can't find, return -999
		*/
		int _findFlag(int istart, int istop, const vector<int>& flg, string flag_type, string direction);

		/**
		* @brief  check the data flag and string-decribtion are consistant. If yes, return true.
		* @param[in] iflag  input iflag
		* @param[in] type   string description
		* @return           true are consistent, otherwise is not
		*/
		bool _isTrue(int iflag, string type);

		/**
		* @brief  search no more than 5 points are needs forwards and backwards to check LG Jump
		* @param[in] flg        searched flags
		* @param[in] crt_epo    current epoch
		* @param[in] limit_epo  set Search range (if "-",use_epo > limit_epo ; if "+", use_epo < limit_epo)
		* @param[in] mode       "+"  means forward, while "-" means backward
		* @param[in] use_epo    result
		* @return               number of eligible epochs
		*/
		int _searchEpoch(const vector<int>& flg, int crt_epo, int limit_epo, string mode, int& use_epo);

		/**
		* @brief     transform
		* @param[in] data_flag       string format of data_flag
		* @return                    DATAFLAG type of data_flag
		*/
		int _str2DataFlag(string data_flag);
	};
	
}

#endif /*GTURBOEDIT_H*/