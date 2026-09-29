#ifndef GPVTFLT_OLD_H
#define GPVTFLT_OLD_H


/**
* @file gpvtflt.h
* @brief implements pvt client

* @verbatim
		History
		-1.0	ShenZhiheng	       2019-01-01 creat the file.
		-1.1	ShenZhiheng	       2019-03-20 add some comment.
		-1.2	ShenZhiheng	       2020-04-01 add ambiguity fixed function
  @endverbatim

* @author ShenZhiheng
* @version 1.0
* @date 2019-01-01
* @license
*/
//#include "gproc/gspplsq.h"
#include "gmodels/gpppmodel.h"
#include "gproc/gpppflt.h"
#include "gambfix/gambiguity.h"
#include "gdata/gposdata.h"
#include "gutils/ExportINS.h"
using namespace great;
//using namespace gsins;

namespace gsins
{
	class INS_LIBRARY_EXPORT t_gpvtflt_old :public t_gpppflt
	{
	public:

		/** @brief default constructor. */
		t_gpvtflt_old();

		/**
		* @brief Constructor
		* @note set parameter value by setS
		*/
		t_gpvtflt_old(string mark, string mark_base, t_gsetbase* set);


		/** @brief default destructor. */
		virtual ~t_gpvtflt_old();

		/** @brief processBatch. */
		virtual int processBatch(const t_gtime& beg, const t_gtime& end, bool prtOut);

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
		* @param[in] current epoch
		* @param[in] current site
		* @param[in] current site satellite data
		* @param[in] prepare process or not
		* @return void */
		void Add_rho_azel(/*const t_gtime &runEpoch, */const string& site_name, const t_gtriple& xyz_s, const t_gtriple& xyz_r, t_gsatdata& obs_sat);

	protected:

		/**
		* @brief initiate rtk parameters
		*/
		int _rtkinit();
		int _setCrd();
		/**
		* @brief select common sats between base and rover sites
		* @param[in] base-site satellite data
		* @param[in] rover-site satellite data
		* @return number of common sats
		*/
		int _selcomsat(vector<t_gsatdata>& data_base, vector<t_gsatdata>& data_rover);
		/** @brief kalman filter time update. */
		void _udsdAmb();
		/*obtain satcrd,satclk,satele,satazi,rho, it must be two bands*/
		virtual int _prepareData(const string& ssite, vector<t_gsatdata>& sdata);

		virtual int _processEpoch_new(const t_gtime& runEpoch);

		virtual int _outlierDetect(const ColumnVector& v, const SymmetricMatrix& Qsav, string & sat);
		virtual int _outlierDetect(const Eigen::VectorXd& v, const SymmetricMatrix& Qsav, string & sat);

		virtual string _getsat_from_obs(int idx);

		bool _valid_residual(bool phase_process, string sat_name, enum FREQ_SEQ &f, map<string, map<FREQ_SEQ, pair<int, int>>>& index_l);
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

		/**
		* @brief add pseudorange and carrier-phase measurement to model. Rewritten for RTK application.
		*
		* @param[in] satdata	satellite data
		* @param[in] iobs		index
		* @param[in] param      undifferenced params
		* @param[in] XYZ		ECEF coordinate
		* @param[out] A			Coff Matrix
		* @param[out] l			measurement vector
		* @param[out] P			measurement variance
		* @return
			retval > 0.   success
			retval < 0.   fail
		*/

		int _combineDD(Matrix& A, SymmetricMatrix& P, ColumnVector& l);
		virtual int _addObsP(t_gsatdata& satdata, unsigned int& iobs, Matrix& A, ColumnVector& l, SymmetricMatrix& P);

		virtual int _addObsL(t_gsatdata& satdata, unsigned int& iobs, Matrix& A, ColumnVector& l, SymmetricMatrix& P);

		virtual int _combineMW(t_gsatdata& satdata);

		virtual bool _external_pos(const t_gtriple& xyz_r, const t_gtriple& rms);

		/** @brief prepare data.	*/
		virtual int _prepareData();

		/** @brief process one epoch Velocity. */
		virtual int _processEpochVel();

		void _print_ele();

		/** @brief kalman filter time update. */
		virtual void _predict();

		/** @brief print  the result. */
		void _prtOut(t_gtime& epo, t_gallpar& X, const SymmetricMatrix& Q, vector<t_gsatdata>& data, ostringstream& os, xml_node& node, bool saveProd = true);

		void _get_result(t_gtime& epo, t_gposdata::data_pos& pos);

		/** @brief print  the result. */
		void _prtOutHeader();

		void _combineEWL(t_gsatdata& satdata);

		t_whitenoise*       _dclkStoModel;		/// clock drift model 
		t_whitenoise*       _velStoModel;		/// velocity model

		t_gtime _beg_time;					            ///< beg time of proc
		t_gtime _end_time;				                ///< end time of proc
		bool _isFirstFix;
		t_gupd* _gupd = nullptr;
		t_gambiguity* _ambfix = nullptr;
		t_gsetamb::FIX_MODE  _fix_mode;
		t_map_MW _MW, _EWL, _LW, _LE, _LWL, _rMW, _rEWL, _IMW, _IEWL, _PW, _AFIF;
		bool _amb_state;
		int _fixed_pct;
		t_gtriple _vel;
		SymmetricMatrix     _Qx_vel;

		double _max_res_norm;			// maximum normalize res
		bool   _isBase;					/**<if using base,true; other: fasle */
		bool   _sd_sat;					// single differented between sat and sat_ref
		bool _pos_constrain;
		t_gtriple _extn_pos, _extn_rms; // external position and its rms
		string _site_base;				// <base site name;
		vector<t_gsatdata> _data_base;	// < base data in every epoch;
		vector<pair<string, pair<FREQ_SEQ, GOBSTYPE>>> _obs_index;
		t_gpppmodel* _gModel_base;
	};
}






#endif