/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2019-04-02  created
	 -1.1 ZhengHJ  2019-04-18  add Doxygen format comment
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gppplsq.h
* @brief		main control class for ppp estimate
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-18
*
*/


#ifndef GPPPLSQ_H
#define GPPPLSQ_H

#include <algorithm>
#include "gio/gfile.h"
#include "gio/gxml.h"
#include "gall/gallambupd.h"
#include "gcoders/ambupd.h"
#include "gambfix/gambiguity.h"
#include "gmodels/glsqproc.h"
#include "gdata/gupd.h"
//#include "gdata/gupd_epoch.h"
#include "gproc/gspplsq.h"
#include "gproc/gqualitycontrol.h"
#include "gutils/ExportLSQ.h"
namespace great
{

	/**
	* @brief main control class for PPP
	*/
	class LSQ_LIBRARY_EXPORT t_gppplsq :public t_glsqproc, public t_gxml
	{
	public:

		t_gppplsq(string site, t_gsetbase* set, t_gallproc* data, t_glog* log);
		virtual ~t_gppplsq();

		bool GenerateProduct() override;
		bool ProcessBatch(t_gallproc* data, const t_gtime& beg, const t_gtime& end) override;

		//Wangbo add: for real time NRTK/PPP-RTK
		bool InitProc(t_gallproc* data,const t_gtime& beg, const t_gtime& end);
		virtual int  ProcessOneEpoch(const t_gtime& now, vector<t_gsatdata>* obs=NULL);
		bool ExtractResfile();

		void setObsMode(OBSCOMBIN obsmode) { _obs_mode = obsmode; };
		void setClient(bool isClient) { _isClient = isClient; };
		void setWLParams(t_gallpar& params);
		void resetLsqUpdate(t_gsetbase* set, t_gallproc* data, t_glog* log);
		t_gallpar& getAllPar();

	protected:
		virtual void _setOut();
		virtual int _SolveEpoch();

		bool _init_ppp_setings();
		t_gallbias*      _gallbias = nullptr;
		bool _generateAmbupdIF(t_gallambupd* allambupd);
		bool _generateAmbupdALL(t_gallambupd* allambupd);
		//void _combineMW(t_gsatdata& satdata, string ac ="CAS");
		void _combineMW(t_gsatdata& satdata);
		void _combineEWL(t_gsatdata& satdata);
		void _combineEWL24(t_gsatdata& satdata);
		void _combineEWL25(t_gsatdata& satdata);
		void _combineLW(t_gsatdata& satdata);
		void _combineLE(t_gsatdata& satdata);
		void _combineLE24(t_gsatdata& satdata);
		void _combineLE25(t_gsatdata& satdata);
		void _combineLWL(t_gsatdata& satdata);
		void _combineLWL124(t_gsatdata& satdata);
		void _combineLWL125(t_gsatdata& satdata);
		map<GSYS, double> _recWidelaneUPD(string site, std::vector<t_gsatdata>& crt_obs);
		map<GSYS, double> _recExtraWidelaneUPD(string site, std::vector<t_gsatdata>& crt_obs, string mode);
		map<GSYS, double> _recWidelaneUPD_RAW(string site, std::vector<t_gsatdata>& crt_obs);
		map<GSYS, double> _recExtraWidelaneUPD_RAW(string site, std::vector<t_gsatdata>& crt_obs, string mode);
		int _addFixConstraintAFIF(const t_gtime& crt_epoch, std::vector<t_gsatdata>& crt_obs, string mode);
		int _addFixConstraintAFIF_RAW(const t_gtime& crt_epoch, t_glsq*  glsq, string mode);
		
		void _printOut(const t_glsq& lsqestimator, string mode);
		void _printOutHeader();
		void _extractTropDelay(const t_glsq& lsqestimator);

		virtual bool _initLsqProdData(t_gallprod* data) override;
		virtual bool _initLsqProcPars(t_glsq* lsq)      override;
		virtual bool _init_sys_rec_ifb_pars(t_glsq* lsq) override;
		virtual bool _initLsqProcData(t_gallproc* data) override;

		bool _processOneRecOneSat(const t_gtime& crt_epoch, const std::string& rec, const std::string& sat, t_gsatdata& crd_obs) override;
		bool _processOneRec(const t_gtime& crt_epoch, const std::string& crt_rec, std::vector<t_gsatdata>& crt_obs) override;

		bool _getSatRef();
		void _IFConstrainedByWL(t_glsq* lsq);
	protected:
		// data/processing
		t_gupd*       _gupd = nullptr;
		//t_gupd_epoch* _gupd_epoch = nullptr;
		t_giof*       _fltfile = nullptr;		///< outputs every epoch result file
		t_giof*       _enufile = nullptr;		///< outputs enu coord file 
		t_giof*       _tropfile = nullptr;
		t_gambiguity* _ambfix = nullptr;
		t_gspplsq*    _spplsq = nullptr;		///< spplsq, provide init crd
		vector<t_gsatdata> _crt_obs_list;


		// Fixing Ambiguities
		bool _isFirstFix;
		t_map_MW _MW, _EWL, _EWL24, _EWL25, _LW, _LE, _LE24, _LE25, _LWL, _LWL124, _LWL125, _rMW, _rEWL, _rEWL24, _rEWL25, _IMW, _IEWL, _PW, _AFIF;
		t_map_amb_time _AmbValidTime;


		string           _recclk_file;	       ///< outputs recclk file name
		t_gtriple         xyz_standard;	       ///< enu coord reference
		pair<int, t_gtriple>  neu_sum_nepo;
		map<string, int> _glofrq_num;          ///< glofrq number, set GLO frequency number

		// For PPPLEO settings
		bool                            _simulation;     ///< judge PPPLSQ/PPPLEO, wether simulation	

		// add for NPP-RTK
		map<string, double> _crt_ele;            ///< current satellites elevation
		bool _isClient;///< rover in NPP-RTK
		set<string> _sat_ref;
		t_gallpar _wl_params;

		// Quality_control
		t_gqualitycontrol  _quality_control;

		// add for kml
		bool               _kml;
		string             _kml_name;
		xml_node           _kml_line;

		// add for RAW_MIX
		set<string>        _single_freq_sats;

		// add for realtime
		t_gtime _wl_Upd_time;
		t_gtime _ewl_Upd_time;
		t_gtime _ewl24_Upd_time;
		t_gtime _ewl25_Upd_time;

	};
}

#endif/*  GPPPLSQ_H */