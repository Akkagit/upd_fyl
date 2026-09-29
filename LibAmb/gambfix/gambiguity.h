/**
*
* @verbatim
History
-1.0 BoWong  2019-03-05  creat the file.
-1.1 glfeng  2019-03-10  Adjusted some functions and variables.
-1.2 BoWong  2019-04-09  Adding Doxygen Style Code Remarks.
-1.3 lhb	 2019-04-20	 Adding default constructor and lambda_solve() for rtk
-1.4 zhshen	 2020-04-01	 Adding the amb function for kalman filter
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambguity.cpp
* @brief	ambiguity fixing by using LAMBDA.
*
* @author   glfeng BoWong , Wuhan University
* @version	1.0.0
* @date		2019-04-09
*
*/

#ifndef GAMBIGUITY_H
#define GAMBIGUITY_H 

#include <string>
#include <map>
#include <memory>
#include <iostream>
#include "gdata/gupd.h"
//#include "gdata/gupd_epoch.h"
#include "gset/gsetamb.h"
#include "gset/gsetout.h"
#include "gambfix/glambda.h"
#include "gambfix/gbdeci.h"
#include "gambfix/gambcommon.h"
#include "gutils/gmatrixconv.h"
#include "gutils/gcycleslip.h"
#include "gproc/glsq.h"
#include "gproc/gflt.h"
#include "gutils/gsys.h"
#include "gutils/ginfolog.h"


//#define DEBUG 1

using namespace std;
using namespace gnut;

namespace great 
{
	/**
	* @brief class for storing ambiguity resolution common value.
	*/
	class AMB_LIBRARY_EXPORT t_gamb_cmn
	{
	public:
		t_gamb_cmn();
		t_gamb_cmn(const t_gtime& t, t_glsq* lsq);
		t_gamb_cmn(const t_gtime& t, t_gflt* flt);

		~t_gamb_cmn();


		void active_amb(map<string, int> active_amb);
		int active_amb(string site);
		void amb_fixed(bool b);
		bool amb_fixed();
		void set_ratio(double r);
		double get_ratio();
		void set_boot(double b);
		double get_boot();

		t_gtime now()const;
		double sigma0()const;

		t_gallpar param()const;
		ColumnVector dx()const;
		SymmetricMatrix Qx()const;

	private:
		t_gtime _now;
		int _nobs_total, _npar_number;
		double _sigma0;
		double _vtpv;
		t_gallpar _param;
		ColumnVector _dx;
		ColumnVector _stdx;
		SymmetricMatrix _Qx; 
		map<string, int> _active_amb;
		bool _amb_fixed;
		double _ratio;
		double _boot;
	};
	/**
	* @brief class for fix ambiguities epoch-wisely.
	*/
	class AMB_LIBRARY_EXPORT t_gambiguity
	{
	public:

		/** @brief default constructor. */
		t_gambiguity();

		/*@brief constructor, init some internal variables.
		* @param[in]  site	site name.
		* @param[in]  gset	ambiguity fixing setting of xml.
		*/
		t_gambiguity(string site, t_gsetbase* gset);

		/** @brief default destructor. */
		virtual ~t_gambiguity();

		/**
		* @brief batch processing loop over epoch.
		* @param[in] cycle_slip Initialize total number of the checked ambiguity set.
		* @param[in] glsq       least squares estimations include parameters and Qx dx etc.
		* @param[in] t          current epoch.
		* @return	 0 - compute successfully
		*            other - failure
		*/
		virtual int processBatch(const t_gtime& t, t_glsq*  glsq, t_gcycleslip* cycle_slip, string mode);
		virtual int processBatch(const t_gtime& t, t_gflt*  gflt, t_gcycleslip* cycle_slip,string mode);

		bool amb_fixed();
		/** 
		*@brief set upd file. 
		* @param[in] gupd upd data
		*/
		void setUPD(t_gupd* gupd) { _gupd = gupd; }
		//void setUPD_epoch(t_gupd_epoch* gupd_epoch) { _gupd_epoch = gupd_epoch; }

		/**
		*@brief set MW.
		* @param[in] t_map_MW& MW
		*/
		void setMW(map<string, map<int, double>>& MW) { _MW = MW; }

		/**
		*@brief set elevation in raw_all mode.
		* @param[in] t_map_MW& MW
		*/
		void setELE( map<string, double> &_crt_ele ) { _ELE = _crt_ele; }
		/** 
		*@brief set general log file. 
		* @param[in] log log file
		*/
		void setLOG(t_glog* log) { _glog = log; }

		/**
		*@brief set IF/ALL.
		* @param[in] OBSCOMBIN type
		*/
		void setObsType(OBSCOMBIN type) { _obstype = type; }

		/**
		*@brief set glonass freq ID.
		* @param[in] t_gallnav* allnav
		*/
		void setWaveLength(map<string, int>& glofrq_num);

		/**
		* @brief Set ambiguity resolution reference satellite, 
		*        used in AR with ref sat in PPP-RTK NRTK etc.
		* @param[in] satRef  Reference satellites of each system or all system.
		* @return	 True  - set reference satellite successfully
		*            False - failure
		* lvhb,wb added or modified for npp 202006
		*/
		void setSatRef(set<string>& satRef) { _sat_refs = satRef; };
		void updateFixParam(t_gallpar& param,ColumnVector& dx,ColumnVector* stdx = NULL); 
		t_gallpar& getFinalParams();	
		t_gallpar& updateUCAmb(set<string>* satref = NULL);
		t_DD_amb& getDD() { return _DD;}	//lvhb added for rtk/ins, 20200826
		double get_ratio() { return _outRatio; };

		/** 
		* @brief set ratio threshold and lambda reduction and search 
		*/
		double lambdaSolve(double & ratio, const Matrix & anor, const ColumnVector & fltpar, ColumnVector& ibias, bool parlamb = false);
		//double lambdaSolve(double& ratio,const Matrix& anor, const vector<double>& fltpar, vector<int>& ibias,bool parlamb = false);

		ColumnVector out_firstambD() { return _mDia; }

		void setActiveAmb(int max) { _max_active_amb_one_epo = max; }

	protected:
		CONSTRPAR _crd_est;
		OBSCOMBIN _obstype;  ///< UDUC/IF ambiguity fixing
		t_gtime _beg;        ///< begin epoch
		t_gtime _end;        ///< end epoch
		t_gtime _crt_time;   ///< current epoch
		double _interval;    ///< sampling interval /s
		set<string> _sys;    ///< system

		string _site;         ///< current site 
		set<string> _sat_rm;  ///< satellites being removed
		set<string> _sat_refs; ///< reference satellites

		int                                  _frequency;
		t_DD_amb                             _DD;            ///< DD over all baslines
		t_DD_amb                             _DD_save;            ///< DD over all baslines
		map<string, map<int, double> >       _MW;            ///< MW for IF ambiguity fixing of current time
		map<string, double>                  _ELE;
		map<GSYS, map<FREQ_SEQ, GOBSBAND> >  _band_index;
		map<string, map<string, double> >    _sys_wavelen;   ///< wave length
		map<string, map<string, bool> > _WL_flag;//  <DD WL flag: true or false
		map<string, map<string, bool> > _EWL_flag;//  <DD EWL flag: true or false
		map<string, map<string, bool> > _EWL24_flag;//  <DD EWL24 flag: true or false
		map<string, map<string, bool> > _EWL25_flag;//  <DD EWL25 flag: true or false
		map<string, map<string, int> > _IWL;//  < Fixed DD WL ambiguity
		map<string, map<string, int> > _IEWL;//  < Fixed DD EWL ambiguity
		map<string, map<string, int> > _IEWL24;//  < Fixed DD EWL24 ambiguity
		map<string, map<string, int> > _IEWL25;//  < Fixed DD EWL25 ambiguity
		t_gupd*       _gupd;        ///< upd data
		//t_gupd_epoch* _gupd_epoch;  ///< upd data
		t_gsetbase*   _gset;        ///< set from xml
		t_glog*       _glog;        ///< log file
		
		t_gsetamb::FIX_MODE  _fix_mode;                          ///< set ambiguity fixing mode
		t_gsetamb::UPD_MODE  _upd_mode;                          ///< set wl and nl upd mode
		bool                 _part_fix;                          ///< set whether take partial ambiguity fixed mode, default is false
		int                  _part_fix_num;                      ///< if value's size less than num, stop part fix
        double               _ratio;                             ///< threshold in LAMBDA method
		double               _boot;                             ///< threshold of bootstrapping rate in amb fix
        double               _min_common_time;                   ///< the Minimum common time of two observation arc

		map<string, double> _map_EWL_decision;                    ///< deriation, sigma in WL/NL-cycle
        map<string, double> _map_WL_decision;                    ///< deriation, sigma in WL/NL-cycle
		map<string, double> _map_NL_decision;

		bool _is_first;                                          ///< whether it is first epoch to be fixed
		bool _is_first_wl;
		bool _is_first_ewl, _is_first_ewl24, _is_first_ewl25;
		bool _is_first_nl;
		ColumnVector _mDia;

		double _outRatio;										 ///< outRatio
		t_giof* _ratiofile = nullptr;                            ///< ratio file
		t_giof* _bootfile = nullptr;                             ///< BootStrapping file

		double* _pdE = nullptr;
		double* _pdC = nullptr;

		ostringstream _os_ratio;
		ostringstream _os_boot;

		int _max_active_amb_one_epo;

		bool _amb_fixed;
		int _total_amb_num, _fixed_amb_num;

		t_gallpar _param;
		map<string, string>  _single_freq_sats;

		// for real-time
		bool    _realtime;
		t_gtime _ewl_Upd_time, _ewl24_Upd_time, _ewl25_Upd_time;
		t_gtime _wl_Upd_time;
   
	protected:
		
		/**
		* @brief get UPD correction for sat
		* @param[in] mode - UPD mode[optional choose: NL WL]
		* @return	 true - compute successfully , false - failure
		*/
		bool _getSingleUpd(string mode, t_gtime t, string sat, double& value, double& sigma);

		/**
		* @brief apply UPD correction for sat
		* @param[in] t    - target epoch
		* @return	 true - compute successfully , false - failure
		*/
		bool _applyUpd(t_gtime t);
		bool _applyWLUpd(t_gtime t, string mode);

		/**
		* @brief check if a new ambiguity is dependant of the already selected set expressed by a set of orth. unit vector Ei.
		* @param[in] isFirst            judge whether is the first.
		* @param[in] iNamb              number of related one-way ambiguties.
		* @param[in] iNdef              number of already selected independent ones if independent +1 as output.
		* @param[in] iN_oneway          number of oneway ambiguities used in the to be checked input .
		* @param[in] arriIpt2ow         position of the correponding ow-ambiguities.
		* @param[in] iMaxamb_ow         number of total number of ow-ambiguties for the differenced ambiguties to be checked, not zero, for allocation, zero for check.
		* @param[in] iMaxamb_for_check  number of the total independent ones of the checked set, only used for the allocaiton .
		* @param[out] iNdef
		* @param[out] arriIpt2ow
		* @return	 true - compute successfully , false - failure
		*/
		bool _checkAmbDepend(bool isFirst, int iNamb, int *iNdef, int iN_oneway, int *arriIpt2ow, int iMaxamb_ow, int iMaxamb_for_check);


		/**
		/**
		* @brief define double-difference ambiguities over one baseline.
		* @param[in] amb_cmn    estimations include parameters and Qx dx etc.
		* @return	 true - compute successfully , false - failure
		*/
		int _defineDDAmb(t_gamb_cmn* amb_cmn);
       /*
		* @brief define double-difference widelane ambiguities over one baseline.
        * @return	 true - compute successfully , false - failure
        */
		bool _calDDAmbWL();
		bool _calDDAmbWLWL();//Wangbo: add for PPP-RTK Client WL
		bool _calDDAmbWLALL(t_gamb_cmn* amb_cmn, string mode);

		/**
		* @brief fix widelane and narrowlane ambiguities.
		* @return	 true - compute successfully , false - failure
		*/
		bool _fixAmbIF();
		bool _fixAmbUDUC();
		bool _fixAmbWL();
		bool _fixAmbWL(string mode);
		

		/**
		* @brief select from (or reorder) a set of DD-ambiguities with their widelane and narrowlane ambiguities.
		* First, with fixed widelane, and their narrowlane near to an integer
		* Second, with fixed widelane
		* Last, form an full range and indepandent set
		* @param[in]  korder  =1, both fixed, =2 also not fixed.
		* @param[in]  namb    number of ambiguity
		* @param[out] ndef    number of fixed ambiguity
		* @return	 true - compute successfully , false - failure
		*/
		int _selectAmb(int korder, int namb);

		
		/**
		* @brief get covariance-matrix for or DD-ambiguities based on the defined DD-ambiguities and the Qx with oneway.
		* @param[in] filter      kalman filter estimations include parameters and Qx dx etc.
		* @param[in] covariance  Covariance matrix of DD.
		* @param[int] value      rnl column matrix of DD.
		* @param[out] covariance
		* @param[out] value
		* @return	 true - compute successfully , false - failure
		*/
		bool _prepareCovariance(t_gamb_cmn* amb_cmn, SymmetricMatrix& covariance, vector<double>& value);
		bool _prepareCovarianceWL(t_gamb_cmn* amb_cmn, SymmetricMatrix& covariance, vector<double>& value, string mode);

		/**
		* @brief resolve integer ambiguities using LAMBDA-method.
		* @param[in] anor    inverted N-matrix, full-matrix stored in one dim.
		* @param[in] fltpar  parameter of free solution, ambiguities at the end.
		* @param[int] ibias  fixed solution.
		* @param[out] ibias
		* @return	  ratio value
		*/
		double _lambdaSearch(const Matrix& anor, const vector<double>& fltpar, vector<int>& ibias,double* boot);

		/**
		/**
		* @brief resolve integer ambiguities using LAMBDA-method [combination of _prepareCovariance+_lambdaSearch].
		* @param[in]  amb_cmn    estimations include parameters and Qx dx etc.
		* @param[int] fixed_amb      fixed solution.
		* @return	 true - compute successfully , false - failure
		*/
		bool _ambSolve(t_gamb_cmn* amb_cmn, vector<int>& fixed_amb, string mode);

		/**
		* @brief applied fixing constraints into Qx to get the fixed solutiond.
		* @param[in] glsq  least squares estimations include parameters and Qx dx etc.
		* @return	 true - compute successfully , false - failure
		*/
		bool _addFixConstraint(t_glsq* glsq);
		bool _addFixConstraint(t_gflt* gflt);
		bool _addFixConstraint_before(t_gflt* gflt);
		bool _addFixConstraintWL(t_glsq* glsq, string mode);
		bool _addFixConstraintWL(t_gflt* gflt,string mode);

		/**
		* @brief init ratio-file.
		* @return	void
		*/
		void _initRatiofile();
		void _initBootfile();
		/**
		* @brief write ratio-file.
		* @return	void
		*/
		void _writeRatio(double ratio);
		void _writeBoot(double BootStrapping);

		/**
		* @brief Calculate GLONASS receiver upd <Single Difference can't eliminate receiver upd in FDMA>.
		* @return	void
		*/
		double _glonassRecUpd(string site);

		/**
		* @brief Find Single-Difference include reference satellite
		* @return  True  - Successful
		*          False - Failure
		*/
		bool _findRefSD();

		/**
		* @brief Recover No-difference ambiguity both frequency.
		* @param[in] satref Reference Satellites.
		* @param[in] glsq   Least squares estimations include parameters and Qx dx etc.
		* @param[in] MW     MW obs.
		* @return	 True  - No-Difference  ambiguity recovery successfully
		*            False - failure
		*/
		bool _NDRecovery(string satref, t_gallpar& params, map<string, map<FREQ_SEQ, double>>& singleND);

	};

	bool _ddCompare(const t_dd_ambiguity& dd1, const t_dd_ambiguity& dd2);
}

#endif
