#ifndef GINTEGRATION_H
#define GINTEGRATION_H

/**
* @file gintegration.h
* @brief loosely/tightly coupled integration key algorithm
* @note
		GNSS PPP/RTK and INS integration,including LCI, TCI
		most algorithm is transplanted from the t_glci, t_gtci, t_grtci.

* @verbatim
		History
		-1.0	ShenZhiheng	       2020-05-27 creat the file.
  @endverbatim

* @author ShenZhiheng
* @date 2020-05-27
* @license
*/

#include "gins/gins.h"
#include "gutils/gmutex.h"
#include "gset/gsetbase.h"
#include "gdata/gimudata.h"
#include "gio/glog.h"
#include "gio/giof.h"
#include "gset/gsetign.h"
#include "gset/gcfg_ign.h"
#include "gproc/gpvtflt.h"
#include "gdata/gposdata.h"
#include "gutils/ExportINS.h"
using namespace gnut;

namespace gsins
{

	/**
	* @class t_gintegration
	*
	* t_gintegration Class for loosely-coupled integration
	*
	* loosely-coupled integration based on Adaptive and Sequential kalman filter
	* Inherited from t_gsinskf
	*
	*/
	class INS_LIBRARY_EXPORT t_gintegration : virtual public t_gsinskf,
											  virtual public t_gpvtflt
	{
	public:
		/**
		* @brief Constructor
		* @note set parameter value by gset
		*/
		t_gintegration(string site, string site_base, t_gsetbase* gset, t_gallproc *allproc);


		virtual void Add_IMU(t_gimudata* imu);

		/**
		* @brief Initial Function
		*
		* virual function,inherited from t_gsinskf
		*/
		virtual void init() {}

		/**
		* @brief set measurement matrix Hk
		* virual function,and inherited from t_gsinskf
		*/
		virtual void set_Hk();

		
		/** @brief processBatch. */
		virtual int processBatch(const t_gtime& beg, const t_gtime& end, bool prtOut)override;
		virtual int processBatch(const t_gtime& beg, const t_gtime& end);
		virtual int processBatch2(const t_gtime& beg, const t_gtime& end, bool prtOut,bool beg_end);
		virtual int processBatch_test(const t_gtime& beg, const t_gtime& end);

		//set the proportion of processing
		virtual void Set_pro(t_gtime *fwd_crt, t_gtime *bwd_crt, double* fwd_per, double *bwd_per);

	protected:

		virtual int _init();

		virtual bool _ins_init();
		virtual bool _gnss_init();

		virtual bool _valid_ins_constraint();

		/** @brief prepare data.	*/
		virtual int _prepareData() override;
		/*obtain satcrd,satclk,satele,satazi,rho, it must be two bands*/
		virtual int _prepareData(const string& ssite, vector<t_gsatdata>& sdata)override;
		/** @brief process one epoch.	*/
		virtual int _processEpoch(const t_gtime& runEpoch) override;

		virtual int _merge(Matrix& A);
		virtual int _merge_init();

		virtual int _GNSS_Update();
		virtual int _ZUPT_Update();
		virtual int _NHC_Update();
		virtual int _ATT_Update();  ///add wh 
		virtual int _ZUPT_POS_Update();
		virtual MEAS_TYPE _getPOS(t_gposdata::data_pos& pos);
		virtual MEAS_TYPE _get_DD_Fixed_Obs();
		virtual int _getMeas();
		virtual int _getMeas_test();

		virtual bool _time_valid(t_gtime gt, t_gtime inst);

		virtual bool _avaliable(t_gtime crt);

		virtual void time_update(double kfts, int fback = 1, double inflation = 1.0);
		virtual int _meas_update(double fading = 1.0);
		
		virtual bool _set_meas();
		virtual bool _resize();

		virtual int outlier_detect();

		/** @brief Kalman Feedback. */
		virtual void feedback(double kfts);

		virtual int _gnss_feedback(const ColumnVector& dx);

		virtual int _reset();

		virtual int _prt_ins_kml();

		void _prt_port();

		///wh
		virtual bool _ZUPT_Check();


	protected:

		t_gtime _ins_beg, _gnss_beg;
		t_gtime _ins_end, _gnss_end;
		t_gtime _ins_crt, _gnss_crt;
		IGN_TYPE _ign_type;
		Eigen::MatrixXd _gnss_design_mat;
		Eigen::MatrixXd _gnss_meas_cov;
		Eigen::VectorXd _gnss_meas_vector;
		t_gimudata* _imudata;
		bool _initial_merge;
		//bool _sync;                           // GNSS and INS processing are Synchronous;

		///  for foot mounted ZUTP 
		vector<double> _Tz; 
		vector<double> cal_Tz();
		Eigen::Vector3d _mf;
		vector<Eigen::Vector3d> _f;
		double _thresh;
		t_gtime *_fwd_crt, *_bwd_crt;   //forward,backward time 
		double *_fwd_per, *_bwd_per;   //forward,backward processing percentage 

	};

}




#endif