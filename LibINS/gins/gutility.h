#ifndef UTILITY_H
#define UTILITY_H

/**
* @file gutility.h
* @brief common algorithm
* @details implement common algorithm and define some enums

* @verbatim
		History
		-1.0	ShenZhiheng	       2018-10-15 creat the file.
		-1.1	ShenZhiheng	       2019-03-20 add some comments.
  @endverbatim

* @author ShenZhiheng
* @version 1.0
* @date 2018-10-15
* @license
*/


#include <string>
#include "Eigen/Eigen"
#include "gutils/gsysconv.h"
#include "gutils/ExportINS.h"

namespace gsins
{

	#define UNICORE_ACC_SCALE          400 / pow(2,31)    //m/s/LSB
	#define UNICORE_GYRO_SCALE         2160 / pow(2,31)   //deg/LSB
	#define STARNETO_ACC_SCALE         0.05 / pow(2, 15)
	#define STARNETO_GYRO_SCALE        0.1 / (3600 * 256.0)
	#define STARNETO_G                 9.80144145
	#define IMR_FSAS_GYRO_SCALE        0.1 / pow(2,8)     // arcsec/LSB
	#define IMR_FSAS_ACCE_SCALE        0.05 / pow(2,15)   // m/s/LSB

	/**
	* @enum FLT_TYPE
	* @brief filter mode,included forward,backward,forward and backward,RTS
	*/
	enum FLT_TYPE { FORWARD, BACKWARD,FBS, RTS };

	/**
	* @enum CPS_TYPE
	* @brief Compensation mode,included Cone,Poly,One plus previous
	*/
	enum CPS_TYPE { CONE, ONE_PLUS_PRE, POLY };

	/**
	* @enum MEAS_TYPE
	* @brief Measurement mode,included GNSS,ZUPT,Motion Constrain,Yaw et al
	*/
	enum MEAS_TYPE {
		GNSS_MEAS = 1,
		POS_MEAS, VEL_MEAS, POS_VEL_MEAS,
		DD_UC_OBS_MEAS, DD_IF_OBS_MEAS,
		FIXED_DD_UC_OBS_MEAS,

		MOTION_MEAS = 20,
		ZUPT_MEAS, NHC_MEAS, YAW_MEAS, ODO_MEAS, Hgt_MEAS,
		ATT_MEAS, ZUPT_POS_MEAS,

		NO_MEAS = 100,
		VIS_MEAS,DEFAULT_MEAS = 999,
		LIDAR_MEAS=200
	};

	/**
	* @enum MOTION_TYPE
	* @brief motion mode
	* @note be used for imu simulation,includes different motion mode
	*/
	enum MOTION_TYPE {
		m_keep, m_accelerate, m_decelerate,
		m_yawleft, m_yawright, m_pitchup, m_pitchdown,
		m_rollleft, m_rollright, m_turnleft, m_turnright,
		m_climb, m_descent, m_s, m_8, 
		m_straight, m_static, m_default
	};

	/**
	* @enum UNIT_TYPE
	* @brief imu data unit
	*/
	enum UNIT_TYPE { RAD, DEG, RPS, DPS, RPH, DPH, MPS, MPS2 };


	/**
	* @enum ODR_TYPE
	* @brief INS output format
	*/
	enum ODR_TYPE { XYZ_XYZ_PRY_NW, BLH_ENU_PRY_NW, HOLO_ODO };

	/**
	* @enum LCI_TYPE
	* @brief GNSS position mode for loosely conpled integration mode 
	*/
	enum LCI_TYPE { SPP, PPP, RTK, PPPRTK };

	/**
	* @enum ALIGN_TYPE
	* @brief INS align method, SA is static base alignment, MA is moving base alignment
	* and the PVA is position vector alignment, and the VVA is velocity vector alignment
	*/
	enum ALIGN_TYPE { AUTO, SA, MA, PVA, VVA, ALIGN_DEFAULT = 999 };

	/**
	* @enum IGN_TYPE
	* @brief GNSS/INS integration option.
	* PURE_INS is ins process only;
	* LCI is loosely coupled integration;
	* TCI is tightly coupled intergration;
	* STCI is semi-tightly coupled integration.
	* MULTIGN_LC is GNSS/INS/VISION loosely coupled integrartion
	* MULTIGN_TC is GNSS/INS/VISION tightly coupled integrartion
	* MULTIGN_STC is GNSS/INS/VISION semi-tightly coupled integrartion
	* MULTIGN_LIDAR_LC is GNSS/INS/LIDAR loosely coupled integrartion
	* MULTIGN_LIDAR_TC is GNSS/INS/LIDAR tightly coupled integrartion
	* MULTIGN_LIDAR_STC is GNSS/INS/LIDAR semi-tightly coupled integrartion
	* LIDAR_INS is INS/LIDAR tightly coupled integration
	* VIO is INS/VISION tightly coupled integration
	*/
	enum IGN_TYPE { PURE_INS, LCI, TCI, STCI, MULTIGN_LC, MULTIGN_TC, MULTIGN_STC, VIO , MULTIGN_LIDAR_LC, MULTIGN_LIDAR_TC, MULTIGN_LIDAR_STC,LIDAR_INS, IGN_DEFAULT };



	class t_gglv;

	INS_LIBRARY_EXPORT extern const t_gglv glv;
	
	/**
	* @class t_gglv
	* @brief define some global const values,such as PI,deg et al
	* @author ShenZhiheng
	*/
	class INS_LIBRARY_EXPORT t_gglv
	{

	public:
		t_gglv(double Re = 6378137.0, double f = (1.0 / 298.257), double wie0 = 7.2921151467e-5, double g0 = 9.7803267714);
	
		double Re, f, wie, g0, mg, ug;

		const double PI = 3.14159265358979;
		const double PI_2 = PI / 2.0;
		const double PI_4 = PI / 4.0;
		const double _2PI = 2 * PI;
		const double EPS = 1e-10;
		const double INF = 1e100;
		const double deg = PI / 180.0;
		const double min = deg / 60.0;
		const double sec = min / 60.0;
		const double ppm = 1.0e-6;
		const double hur = 3600.0;
		const double dps = deg / 1.0;
		const double dph = deg / hur;
		const double dpsh = deg / sqrt(hur);
		const double dphpsh = dph / sqrt(hur);
		const double mpsh = 1 / sqrt(hur);
		const double mpspsh = 1 / 1 / sqrt(hur);
		const double ppmpsh = ppm / sqrt(hur);
		const double secpsh = sec / sqrt(hur);

		double mgpsHz, ugpsHz, mgpsh, ugpsh;

	};


	/**
	* @struct t_scheme
	* @brief  processing scheme
	*/
	struct INS_LIBRARY_EXPORT t_scheme
	{
		double t;						/// add:now time for Unequal intervals
		double ts, T;
		int freq;						/// IMU Freq
		double start, end;				/// start time and end time.
		double delay;					/// the GNSS misalignment time diff with imu.
		double max_pdop;				/// the max PDOP threshold for LCI.
		int min_sat;		     		/// the minimum satellite number threshold.
		int align_time;			    	/// coarse align time.
		int nSamples;
		int Cps;
		bool align;
		int GM; 
		

		t_scheme() {}
		/**
		* @brief Constructor
		* set parameter value by set
		*/
		t_scheme(gnut::t_gsetbase* set);
	};


	/** @brief lagrange interpolation. */
	double Lagrange(double t[], double x[], int n, double t1);  

	int sign(double d);

	/** @brief string to UNIT_TYPE. */
	UNIT_TYPE str2Unit(const string& s);


	/** @brief string to FLT_TYPE. */
	FLT_TYPE str2Flt(const string& s);

	/** @brief string to FLT_TYPE. */
	IGN_TYPE str2ign(const string& s);

	/** @brief string to ODR_TYPE. */
	ODR_TYPE str2odr(const string& s);

	/** @brief string to ODR_TYPE. */
	LCI_TYPE str2lci(const string& s);


	/** @brief string to ALIGN_TYPE. */
	ALIGN_TYPE str2align(const string& s);

	MEAS_TYPE ins2meas(const int& i);

	CPS_TYPE str2cps(const string& s);

	/**
	* @brief Cart to Geod
	* overwrite gnut's transform function
	*/
	INS_LIBRARY_EXPORT Eigen::Vector3d Cart2Geod(const Eigen::Vector3d& CartPos, bool b);	
	
	/**
	* @brief Geod to Cart
	* overwrite gnut's transform function
	*/
	INS_LIBRARY_EXPORT Eigen::Vector3d Geod2Cart(const Eigen::Vector3d& GeodPos, bool b);

	/**
	* @brief Geod to Cart
	* overwrite gnut's transform function
	*/
	INS_LIBRARY_EXPORT Eigen::Vector3d XYZ2ENU(const Eigen::Vector3d& xyz, const Eigen::Vector3d& xyz_ref);

	INS_LIBRARY_EXPORT Eigen::MatrixXd NewMat2Eigen(const Matrix& newmat);
	INS_LIBRARY_EXPORT Matrix Eigen2newMat(const Eigen::MatrixXd& eigen);
	INS_LIBRARY_EXPORT SymmetricMatrix Eigen2BaseMatrix(const Eigen::MatrixXd& eigen);
	INS_LIBRARY_EXPORT Eigen::MatrixXd BaseMatrix2Eigen(const SymmetricMatrix& newmat);
	INS_LIBRARY_EXPORT Eigen::VectorXd Columns2VectorXd(const ColumnVector& newmat);
	INS_LIBRARY_EXPORT ColumnVector VectorXd2Columns(const Eigen::VectorXd& eigen);

}



#endif
