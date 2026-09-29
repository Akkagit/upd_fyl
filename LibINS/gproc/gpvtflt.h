#ifndef GPVTFLT_H
#define GPVTFLT_H

/**
*
* @verbatim
History
-1.0  2020-06-11  created
@endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file		    gpvtflt.cpp
* @brief		main control class for RTK/PPP Processing
*
* @author       lvhb,zhshen,zyx Wuhan University
* @version		1.0.0
* @date		    2020-06-11
				2020-10-25 lvhb simplified old codes,RTK/PPP function
*
*/


#include "gproc/gpppflt.h"
#include "gproc/gqualitycontrol.h"
#include "gmodels/gpppmodel.h"
#include "gambfix/gambiguity.h"
#include "gdata/gposdata.h"
#include "gmodels/gcombmodel.h"
#include "gutils/ExportINS.h"
#include "gproc/gpreproc.h"
#include "gproc/gexeturboedit.h"
#include "gmodels/giono.h"
#include "gproc/gupdateparIF.h"
#include "gproc/gupdateparALL.h"
#include "gnpp/gInterpol.h" //lvhb added in 20201105

using namespace gsins;

namespace great
{
	class INS_LIBRARY_EXPORT t_gpvtflt :public t_gpppflt
	{
	public:

		/** @brief default constructor. */
		//t_gpvtflt();

		/**
		* @brief Constructor
		* @note set parameter value by setS
		*/
		t_gpvtflt(string mark, string mark_base, t_gsetbase* set, t_gallproc* allproc = nullptr);

		/** @brief default destructor. */
		virtual ~t_gpvtflt();

		/** @brief processBatch. */
		virtual int processBatch(const t_gtime& beg, const t_gtime& end, bool prtOut);

		/** @brief Initializing some settings. */
		virtual bool InitProc(const t_gtime& begT, const t_gtime& endT, double* subint = NULL);
		
		/**
		* @brief One epoch processing
		* @note virtual function,lvhb Created for NPP
		* @param[in] Now	current time
		* @param[in] data_rover
		* @param[in] data_base
		* @return
			   -1,failed; 
				0,float;
				1,fixed.
		*/
		virtual int ProcessOneEpoch(const t_gtime& now, vector<t_gsatdata>* data_rover=NULL, vector<t_gsatdata>* data_base=NULL);
		
		/**
		* @brief Get Uncombined Parameters
		* @note lvhb added for NPP
		* @param[out] data_base
		* @return Uncombined parameter
		*/
		t_gallpar& getUCparam(set<string>* satref = NULL); // For NPP,lvhb

		/**
		* @brief Forward and Backward process
		* @note virtual function
		* @param[in] beg	begin time
		* @param[in] end	end time
		* @return
			retval > 0.   success
			retval < 0.   fail
		*/
		virtual int processBatchFB(const t_gtime& beg, const t_gtime& end);

		virtual void Add_UPD(t_gupd* gupd);

		/** @brief calculate azel and rho
		* @Note,lvhb modified for PVT/NPP
		* @param[in] current site name
		* @param[in] current satellite position
		* @param[in] current site  position
		* @param[in/out] current site  data
		* @return void */
		void Add_rho_azel(/*const t_gtime &runEpoch, */const string& site_name,  t_gtriple& xyz_s, const t_gtriple& xyz_r, t_gsatdata& obs_sat);

		/** @brief Setting given reference satellites.
		* @Note, lvhb added for NPP in 202010 
		* @param[in] current reference satellites
		*/
		void setRefsat(map<GSYS, string> satref) { _ipSatRep = satref; }
	protected:

		/*************** RTK process function,created by lvhb *****************/
		
		/**  @brief initiate RTK parameters */
		int _rtkinit();

		/**  @brief Setting base or rover coordinates in SPP/RINEXH/FIX */
		int _setCrd();

		/**  @brief update RTK SD ambiguity parameters */
		void _udsdAmb();

		/** @brief lvhb added to adapting baseline changing in 20200821. */
		virtual void _delPar(const t_gpar::t_type par);

		/**
		* @brief select common sats between base and rover sites,created by lvhb
		* @param[in/out] base-site satellite data
		* @param[in/out] rover-site satellite data
		* @return number of common sats
		*/
		int _selcomsat(vector<t_gsatdata>& data_base, vector<t_gsatdata>& data_rover);

		/**  @brief To jude whether exist UD Phase/Range residuals in current satellite 
		* return true: exist
				 false: no exist
		*/
		bool _valid_residual(bool phase_process, string sat_name, enum FREQ_SEQ& f, map<string, map<FREQ_SEQ, pair<int, int>>>& index_l);
		
		/**
		* @brief add pseudorange and carrier-phase measurement to model. Rewritten for RTK application.
		* @note,Lvhb created  for RTK/SD-PPP
		* @param[in/out] A			Coff Matrix
		* @param[in/out] l			measurement vector
		* @param[in/out] P			measurement variance
		* @return
			retval > 0.   success
			retval < 0.   fail
		*/
		int _combineDD(Matrix& A, SymmetricMatrix& P, ColumnVector& l);

		/*************** rtk process /end *****************/

		/** @brief kalman filter time update. */

		void _syncAmb();
		void _post_turbo_syncAmb();

		int _get_spplsq_result(string site, t_gtime current, t_gtriple& crd, double& clk);

		/*obtain satcrd,satclk,satele,satazi,rho, it must be two bands*/
		virtual int _prepareData(const string& ssite, vector<t_gsatdata>& sdata);
		
		/** @brief prepare data.	*/
		virtual int _prepareData();

		virtual int _processEpoch(const t_gtime& runEpoch);

		virtual int _amb_resolution();

		bool _getSatRef();

		virtual int _outlierDetect(const ColumnVector& v, const SymmetricMatrix& Qsav, string & sat);
		virtual int _outlierDetect(const Eigen::VectorXd& v, const SymmetricMatrix& Qsav, string & sat);

		virtual string _getsat_from_obs(int idx);

		/**
		* @brief add doppler measurement to model
		*
		* @param[in] satdata	satellite data
		* @param[in] iobs		index
		* @param[in] XYZ		ECEF coordinate
		* @param[out] A			Coff Matrix
		* @param[out] l			measurement vector
		* @param[out] P			measurement variance
		* @return
			retval > 0.   success
			retval < 0.   fail
		*/
		int _addObsD(t_gsatdata& satdata, unsigned int& iobs, t_gallpar& param, t_gtriple& XYZ, Matrix& A, ColumnVector& l, DiagonalMatrix& P);

		/** @brief lvhb added post residual for rtk/ins */
		int _postRes(const Matrix& A, const SymmetricMatrix& P, const ColumnVector& l, Matrix& pA,  SymmetricMatrix& pP,  ColumnVector& pl);
		
		/** @brief lvhb modified for PPP/RTK processing */
		virtual int _combineMW(t_gsatdata& satdata);

		virtual bool _external_pos(const t_gtriple& xyz_r, const t_gtriple& rms);
		virtual bool _pos_virtual_obs(const t_gtriple& xyz_r, const t_gtriple& rms, Matrix& A, ColumnVector& l, SymmetricMatrix& P);

		/** @brief process one epoch Velocity. */
		virtual int _processEpochVel();

		void _print_ele();

		virtual int _prt_info(const t_gtime& runEpoch);

		/** @brief kalman filter time update. */
		virtual void _predict();

		/** @brief print  the result. */
		void _prtOut(t_gtime& epo, t_gallpar& X, const SymmetricMatrix& Q, vector<t_gsatdata>& data, ostringstream& os, xml_node& node, bool saveProd = true);

		void _prt_port(t_gtime& epo, t_gallpar& X, const SymmetricMatrix& Q, vector<t_gsatdata>& data);

		void _get_result(t_gtime& epo, t_gposdata::data_pos& pos);

		/** @brief print  the result. */
		void _prtOutHeader();

		void _combineEWL(t_gsatdata& satdata);

		void _equationToMatrix(t_glsqEquationMatrix& equ, Matrix& A, SymmetricMatrix& P, ColumnVector& l,
			unsigned int nPar);
		void _generateObsIndex(t_glsqEquationMatrix& equ);

		void _extern_ion_constraint(const t_gtime& runEpoch);

		bool _calculate_timediff_baseline(const string& site, const t_gtime& epoch0, const t_gtime& epoch1, t_gtriple& baseline);

		bool _post_turbo_select_obs(const t_gtime& epoch, vector<t_gsatdata>& data);

	protected:

		t_whitenoise*       _dclkStoModel;		/// clock drift model 
		t_whitenoise*       _velStoModel;		/// velocity model

		t_gtime _beg_time;					            ///< beg time of proc
		t_gtime _end_time;				                ///< end time of proc
		bool _isFirstFix;
		t_gupd* _gupd = nullptr;
		t_gambiguity* _ambfix = nullptr;
		t_gsetamb::FIX_MODE  _fix_mode;
		t_gsetamb::UPD_MODE  _upd_mode = t_gsetamb::UPD_MODE::UPD;
		t_map_MW _MW, _EWL, _LW, _LE, _LWL, _rMW, _rEWL, _IMW, _IEWL, _PW, _AFIF;
		bool _amb_state;
		int _fixed_pct;
		t_gallpar           _param_fixed;
		t_gtriple _vel;
		SymmetricMatrix     _Qx_vel;
		map<string, double> _crt_ele;

		double _max_res_norm;			// maximum normalize res
		bool   _isBase;					/**<if using base,true; other: fasle */
		bool   _sd_sat;					// single differented between sat and sat_ref
		bool _pos_constrain;
		t_gtriple _extn_pos, _extn_rms; // external position and its rms
		Matrix _post_A; SymmetricMatrix _post_P; ColumnVector _post_l;  // post A P L for RTK

		string _site_base;				///< base site name;
		vector<t_gsatdata> _data_base;	///< base data in every epoch;
		vector<pair<string, pair<FREQ_SEQ, GOBSTYPE>>> _obs_index;
		t_gpppmodel* _gModel_base; 
		set<string> _sat_ref;
		map<GSYS, string> _ipSatRep;   ///< input reference satellite,lvhb
		// for RAW_ALL PPP   : G01, AMB_L1, L1C   X   X   X
		// for IONO_FREE PPP : G01, AMB_IF, L1C L2W   X   X
		// for RAW_ALL RTK   : G01, AMB_L1, L1C   X L1C   X
		// for IONO_FREE RTK : G01, AMB_IF, L1C L2W L1C L2W
		map<pair<string, t_gpar::t_type>, tuple<GOBS,GOBS,GOBS,GOBS>> _amb_obs;
		
		ColumnVector _vBanc_base;
		map< GSYS, map<FREQ_SEQ, GOBSBAND> > _band_index;
		map< GSYS, map<GOBSBAND, FREQ_SEQ> > _freq_index;
		int _frequency;
		// PRECISEMODEL STUFF
		t_gallproc*      _allproc = nullptr;
        shared_ptr<t_gbasemodel> _base_model = nullptr;
		shared_ptr<t_gupdatepar> _update_AMB = nullptr;
		shared_ptr<t_gcycleslip>    _slip12 = nullptr;      ///< turboedit class for IF
		shared_ptr<t_gcycleslip>    _slip13 = nullptr;	    ///< turboedit class for 13 if
		shared_ptr<t_gcycleslip>    _slip14 = nullptr;
		shared_ptr<t_gcycleslip>    _slip15 = nullptr;

		//preprocess
		SLIPMODEL _slip_model;
		bool      _turbo_liteMode;
		shared_ptr<t_gexeturboedit> _gturboedit;
		shared_ptr<t_gpreproc> _gpre;

		// glfeng
		t_gqualitycontrol  _gquality_control;
		t_gionex*          _gionex_GIM;
		bool               _iono_constraint;
		int                _epoch_cout;

		//for npp,lvhb 20200917
		bool _isClient = false;
		t_spt_gInterpol _gInterpol;

		// add for realtime
		t_gtime _wl_Upd_time;
		t_gtime _ewl_Upd_time;

		//by qzy
		bool _realtime;
	};

	/** @brief Test by lvhb 201906 */
	template<class T1, class T2>
	void t_out(T1 const& name, T2 const& matrix) {
		cout << name << endl;
		cout << fixed << setprecision(6) << setw(15) << matrix << endl;
		return;
	}
}






#endif