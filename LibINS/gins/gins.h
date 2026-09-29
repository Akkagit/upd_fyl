#ifndef GINS_H
#define GINS_H

/**
* @file gins.h
* @brief ins core algorithm
* @details
		SINS update algorithm,including aligncoarse,mechanical arrangement.
		SINS initial alignment on moving base.
		SINS base integration class inskf,provide some interface to lci and tci.

* @verbatim
		History
		-0.1	ShenZhiheng	       2018-10-15 creat the file.
		-0.5	ShenZhiheng	       2019-03-20 add transition class t_gsinskf to inheriting for lci and tci.
		-0.7	ShenZhiheng	       2019-12-26 restruction.
  @endverbatim

* @author ShenZhiheng
* @email
* @version 1.0
* @date 2018-10-15
* @license
*/


#include "gbase.h"
#include "gquat.h"
#include "gearth.h"
#include "gavar.h"
#include "gkalman.h"
#include "gimu.h"
#include "gdata/gimudata.h"
#include "gutils/gmutex.h"
#include "gio/gxml.h"
using namespace gnut;

namespace gsins
{


	/**
	* @class t_gsins
	* @brief t_gsins Class is ins core class for ins updating.
	*
	* t_gsins is used for sins core algorithms,including align_coarse,error tranform,mechanical arrangement.
	* and the main function is ins info updating.
	*/
	class INS_LIBRARY_EXPORT t_gsins
	{
	public:


		/**
		* @brief Constructor
		* @note set parameter value
		*/
		t_gsins(const t_gquat& qnb0 = qI, const Eigen::Vector3d& vn0 = Eigen::Vector3d::Zero(), const Eigen::Vector3d& pos0 = Eigen::Vector3d::Zero(), double t0 = 0.0);


		/**
		* @brief Constructor
		* @note set parameter value by gset
		*/
		t_gsins(gnut::t_gsetbase* set);


		/**
		* @brief Constructor
		* @note set parameter value by gset
		*/
		void set_posvel(const Eigen::Vector3d& pos0 = Eigen::Vector3d::Zero(), const Eigen::Vector3d& vn0 = Eigen::Vector3d::Zero());

		/**
		* @brief set Markov time
		*
		* @param[in] tauG		Gyro Markov time
		* @param[in] tauA		Acce Markov time
		*
		*/
		void set_tau(const Eigen::Vector3d& tauG, const Eigen::Vector3d& tauA);

		/**
		* @brief coarse alignment
		*
		* @attention data must be static
		* @param[in] wmm	angle increment
		* @param[in] vmm	velocity increment
		* @return
			initial attitude
		*/
		Eigen::Vector3d align_coarse(const Eigen::Vector3d& wmm, const Eigen::Vector3d& vmm);


		/**
		* @brief coarse alignment in horizon by accelerator
		*
		* @param[in] vmm	velocity increment
		* @return
			initial horizon attitude
		*/
		Eigen::Vector3d acc_align(const Eigen::Vector3d& vmm);


		/**
		* @brief sins mechanical arrangement
		*
		* @param[in] wm		angle increment
		* @param[in] vm		velocity increment
		* @param[in] scm	scheme of processing
		*
		*/
		void Update(const vector<Eigen::Vector3d>& wm, const vector<Eigen::Vector3d>& vm, const t_scheme& scm);


		/**
		* @brief sins error transfer matrix
		*/
		void err_mat();


		/**
		* @brief print ins file header into otringstream
		*/
		void prt_header(ostringstream& os);


		/**
		* @brief print ins parameter into otringstream
		*
		* @verbatim
			position is exported to (X,Y,Z),(B,L,H),....
			velocity is exported to (Vx,Vy,Vz),(Ve,Vn,VU),...
			attitude is exported to (pitch,roll,yaw[-180-180]),...
			eb is exported to (ebx,eby,ebz)[deg/h]
			db is exported to (dbx,dby,dbz)[mg]
		  @endverbatim
		*
		* @param[in] os   out os
		*
		*/
		void prt_sins(ostringstream& os);

		void int_sec(bool _beg_end=true);


		void _debug_ins_info();

	public:
		double nts, t;								/// ultimate interval,now
		t_gimu imu;									/// imu data
		t_gearth eth;								/// earth info
		Eigen::Vector3d att, vn, pos, an;			/// attitude,velocity,position,accelarate
		Eigen::Vector3d fb, fn, wib, web, wnb, vb;	/// special force,angular
		Eigen::Vector3d eb, db;						/// gyro and acce bias
		Eigen::Vector3d _tauG, _tauA;				/// Markov time
		Eigen::Matrix3d Cnb, Cbn, Kg, Ka;			/// DCM,Scale
		t_gquat qnb;								/// quat
		Eigen::Matrix3d Faa, Fav, Fap, Fva, Fvv, Fvp, Fpv, Fpp;   /// Error Matrix
		Eigen::Matrix3d Fng;						/// transfer Matrix from d_enu to d_blh

		Eigen::Vector3d ve, ae, fe, pos_ecef;
		Eigen::Matrix3d Ceb, Cbe;			/// DCM
		t_gquat qeb;

		double output_freq;                /// imu output frepency. wh
		bool output_intsec;                /// imu output in integer second. wh
		map<double, t_gquat> vec_qnb;             /// map qnb for interpolation

		double pure_ins_time;

	protected:

		ODR_TYPE _order;  /// output ins data format
	};



	/**
	* @class t_gsinskf
	* @brief t_gsinskf Class is SINS integration base class for coupled integration
	*
	* t_gsinskf inherits the class t_gkalman1,
	* and implement the purely function of t_gkalman1,such as set_Hk.
	* Moreover,add some constraint for integration,ZUPT,MC,YAW and so on.
	* forward and backward filter.
	*
	*/
	class INS_LIBRARY_EXPORT t_gsinskf : public t_gkalman
	{
	public:

		/** @brief default constructor. */
		t_gsinskf();

		/**
		* @brief Constructor
		* @note set parameter value by gset
		*/
		t_gsinskf(gnut::t_gsetbase *gset,string* filename = NULL);

		/** @brief default destructor. */
		virtual ~t_gsinskf();

		/**
		* @brief Initial Function
		*
		* purely virual function,and should be inherited
		*/
		virtual void init(void) = 0;

		/**
		* @brief Initial Alignment
		*
		* purely virual function,and should be inherited
		*/
		virtual bool align_coarse(const vector<Eigen::Vector3d>& wm, const vector<Eigen::Vector3d>& vm, const t_scheme& scm);

		/**
		* @brief Initial Alignment by position vector
		*
		* purely virual function,and should be inherited
		*/
		virtual bool align_pva(const Eigen::Vector3d& pos, const t_scheme& scm);

		/**
		* @brief Initial Alignment by velocity vector
		*
		* purely virual function,and should be inherited
		*/
		virtual bool align_vva(const Eigen::Vector3d& vel, const t_scheme& scm);

		/**
		* @brief Initial Alignment.
		* @details fast cascaded initial alignment method, including sa, ma, pva, vva.
		* purely virual function,and should be inherited
		*/
		virtual bool cascaded_align(Eigen::Vector3d& pos, Eigen::Vector3d& vel);

		/**
		* @brief Initial Alignment Using Multi-antenna/Moving-antenna Device.
		* @param[in] pos			approximate position(BLH)
		* @param[in] lb				delta-lever vector in body frame
		* @param[in] le				delta-lever vector in ECEF frame
		* @param[in] scm			scheme
		*/
		virtual bool align_moving_antenna(const Eigen::Vector3d& pos,
			const vector<Eigen::Vector3d>& lb,
			const vector<Eigen::Vector3d>& le,
			const t_scheme& scm);

		/**
		* @brief set Ft
		*
		* virtual function
		*/
		virtual void set_Ft(void);

		/**
		* @brief set Ft
		*
		* virtual function
		*/
		virtual void set_Ft1();

		/**
		* @brief set measurement matrix Hk
		*
		* purely virual function,and should be inherited
		*
		*/
		virtual void set_Hk(void) = 0;


		/**
		* @brief set measurement Zk
		*
		* purely virual function,and should be inherited
		*
		*/
		virtual void set_Meas(void) {}

		/**
		* @brief set measurement
		*
		* transmit Zero info for ZUPT
		*
		*/
		virtual void set_meas_ZUPT(double t);

		/**
		* @brief set measurement
		*
		* transmit motion info for motion constraint
		*
		*/
		virtual void set_meas_NHC(double t);

		/**
		* @brief set measurement
		*
		* @note transmit Yaw info for yaw constraint
		*
		*/
		virtual void set_yaw(const double& yaw, double t);

		///add wh
		virtual void set_meas_ATT(double t);
		virtual void set_meas_ZUPT_POS(double t);

		/** @brief body motion state. */
		MOTION_TYPE motion_state();

		/** @brief body measurement state. */
		MEAS_TYPE meas_state();

		/**
		* @brief feedback
		*
		* virtual function
		*
		*/
		virtual void feedback(double kfts);

		/**
		* @brief write info into file
		*/
		virtual void write();

		/**
		* @brief write info into file according to the settings
		*/
		virtual void write_tmp(bool _beg_end = true);
		virtual void write_tmp2(double t, const Eigen::VectorXd & Xk, const Eigen::VectorXd & Pk);

		virtual int _prt_ins_kml() = 0;

		/**
		* @brief print Xk into otringstream
		*
		* @param[in] t				time
		* @param[in] Xk				vector of towrite
		* @param[in] IsError		Xk format
		* @param[out] os			ostringstream
		*
		*/
		virtual void prt_Xk(double t, const Eigen::VectorXd& Xk, ostringstream& os, bool IsError);

		/**
		* @brief print Pk into otringstream
		*
		* @param[in] t				time
		* @param[in] Pk				vector of towrite
		* @param[out] os			ostringstream
		*
		*/
		virtual void prt_Pk(double t, const Eigen::MatrixXd& Pk, ostringstream& os);

		/**
		* @brief fusion algorithm
		*
		* @param[in] Xk1			forward result
		* @param[in] Xk2			backward result
		* @param[in] Pk1			forward variance
		* @param[in] Pk2			backward variance
		* @param[out] Xk_Optimal	optimal result
		* @param[out] Pk_Optimal	optimal variance
		*
		*/
		virtual void fusion(const Eigen::VectorXd& Xk1, const Eigen::VectorXd& Pk1, const Eigen::VectorXd& Xk2, const Eigen::VectorXd& Pk2, Eigen::VectorXd& Xk_Optimal, Eigen::VectorXd& Pk_Optimal);

		/**
		* @brief attitude fusion algorithm
		*
		* @param[in] att1			forward attitude
		* @param[in] att2			backward attitude
		* @param[in] p1				forward variance
		* @param[in] p2				backward variance
		* @return
			optimal attitude
		*/
		virtual Eigen::Vector3d att_fusion(const Eigen::Vector3d& att1, const Eigen::Vector3d& att2, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2);


		/**
		* @brief Kalman filter smoother
		*
		* forward and backward smoother
		*
		*/
		virtual void smoother();


		/** @brief setting. */
		void gset(gnut::t_gsetbase* set) { _setkf = set; };


		/** @brief output setting. */
		void set_out();

		virtual void _restore(const Eigen::VectorXd& X_store, const Eigen::MatrixXd& P_store);
		virtual void _recover(Eigen::VectorXd& X_recover, Eigen::MatrixXd& P_recover);

	public:
		MEAS_TYPE Flag;					            /// Measurement Flag,see also MEAS_TYPE
		FLT_TYPE FilterFlag;			            /// Filter Flag,see also FLT_TYPE
		int fltflag;					            /// Direction Flag
		double tmeas, MeasYaw, MeasVf, MeasHgt;	    /// meas time,meas yaw,meas velocity foward,meas height
		Eigen::Vector3d lever;			            /// arm lever
		Eigen::Vector3d	MeasVel, MeasPos, MeasAtt;  /// meas velosity and position  & attitude addwh
		Eigen::Vector3d	_Cov_MeasVn, _Cov_MeasPos, _Cov_MeasAtt;  /// meas velosity and position convariance
		t_gsins sins;					            /// sins info

		map<double, Eigen::VectorXd> Xk_forward, Xk_backward;		/// forward and backward state
		map<double, Eigen::MatrixXd> Pk_forward, Pk_backward;		/// forward and backward varianc
		map<double, Eigen::MatrixXd> Xk_optimal, Pk_optimal;
		map<double, t_gsins> sins_forward, sins_backward;			/// forward and backward sins info

	protected:

		t_gavar _avar;				// avar
		t_scheme _shm;				// processing scheme (important!)
		int _align_count;			// aligncoarse data number
		bool _aligned;				// aligned
		Eigen::Vector3d wmm, vmm;	// total wm/vm increment (in alignment)
		vector<Eigen::Vector3d> _wm, _vm,_mm;//  wm/vm increment (every epoch)  mm for mag  add wh
		gnut::t_gsetbase* _setkf;	// set
		gnut::t_giof* _fins;		// ins t_giof
		gnut::t_giof* _fkf;			// xk t_giof
		gnut::t_giof* _fpk;			// pk t_giof		
		/*map<double, Eigen::VectorXd> Xk_forward, Xk_backward;		/// forward and backward state
		map<double, Eigen::MatrixXd> Pk_forward, Pk_backward;		/// forward and backward varianc*/
		/*map<double, t_gsins> sins_forward, sins_backward;			/// forward and backward sins info*/

		Eigen::VectorXd _Xk_store;
		Eigen::MatrixXd _Pk_store;
		set<MEAS_TYPE> _Meas_Type;

		Eigen::Vector3d _first_pos;
		bool _first_align;
		string _ofilename = "";

#ifdef BMUTEX
		boost::mutex  _mutex;
#endif
		t_gmutex      _mutex;
	};


	/**
	* @class t_gmove_align
	*
	* SINS alignment on moving base
	*
	*/
	class INS_LIBRARY_EXPORT t_gmove_align
	{
	public:

		/**
		* @brief Constructor
		*
		* set parameter
		*
		*/
		t_gmove_align(const t_gsins& sins0, const Eigen::Vector3d& lever0, double T0);

		/**
		* @brief set external pos and vn
		*
		* @param[in] vn			velocity
		* @param[in] pos		position
		* @param[in] t			time
		*
		*/
		void set_Meas(const Eigen::Vector3d& vn, const Eigen::Vector3d& pos, double t);

		/**
		* @brief alignment
		*
		* @param[in] wm			angle increment
		* @param[in] vm			velocity increment
		* @param[in] scm		scheme of processing
		*
		*/
		void align(const vector<Eigen::Vector3d>& wm, const vector<Eigen::Vector3d>& vm, const t_scheme& scm);

	public:
		double tmeas, nt, T;		/// meas time,total time,interval of meas
		t_gsins sins;				/// ins info
		Eigen::Vector3d pos0, lever, MeasVel, MeasPos, LastPos, LastVn, dvk, T_Alpha, T_Beta;
		Eigen::Matrix3d Ce0n0, Cne, Cee0, Cnn0, Cn0b0, Cb0b;		/// DCM
		t_gquat qb0n0;												/// initial quat
		Eigen::Matrix4d K;											/// cost matrix

	};
}



#endif
