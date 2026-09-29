#ifndef SIMULATION_H
#define SIMULATION_H


/**
* @file gsimulation.h
* @brief inertial measurement unit simulation generator

* @verbatim
		The IMU sensor simulation can be seen as an inverse data processing problem of traditional SINS updating algorithm.
		The class t_gsim is IMU data simulation,including Trajectory Profile Simulation, SIMU Data Simulation.
		step1:Design motion states;
		step2:IMU Data inversion according to motion states;
		step3:add noise and error into IMU  Data;
		step4:simulate GNSS result(Position and Velocity) according to true trajectory
  @endverbatim

* @verbatim
		History
		-1.0	ShenZhiheng	       2018-10-30 creat the file.
		-1.1	ShenZhiheng	       2019-06-10 rename and comment
  @endverbatim

* @author ShenZhiheng
* @version 1.0
* @date 2018-10-30
* @license
*/


#include "gins.h"
#include "gutility.h"
#include "gset/gcfg_simuins.h"
#include <map>
#include "gutils/ExportINS.h"

namespace gsins
{
	/**
	* @struct t_state
	* @brief Motion state parameter at certain time,included V-forward,angle velocity,acceleration
	*
	*/
	struct INS_LIBRARY_EXPORT t_state
	{
		double T, Vy;
		Eigen::Vector3d AngVel;
		Eigen::Vector3d Accel;
	};

	/**
	* @struct t_imu_data
	* @brief one IMU data at certain time
	*
	*/
	struct INS_LIBRARY_EXPORT t_imu_data
	{
		double t;
		Eigen::Vector3d wm;
		Eigen::Vector3d vm;
	};

	/**
	* @struct t_path
	* @brief carrier's physical properties at certain time,included attitude,velocity,position
	*
	*/
	struct INS_LIBRARY_EXPORT t_path
	{
		double t;
		Eigen::Vector3d att;
		Eigen::Vector3d vn;
		Eigen::Vector3d pos;
		Eigen::Vector3d eb;
		Eigen::Vector3d db;
	};

	/**
	* @class t_gsim
	* @brief INS data generator
	*
	* The class t_gsim is IMU data simulation,including Trajectory Profile Simulation, 
	* SIMU Data Simulation.
	*
	*/
	class INS_LIBRARY_EXPORT t_gsim
	{

	public:

		/**
		* @brief Constructor
		* @note set parameter value
		*
		* @param[in] att		initial attitude
		* @param[in] vel		initial velocity
		* @param[in] pos		initial position
		* @param[in] ts			interval
		*
		*/
		t_gsim(const Eigen::Vector3d& att, const Eigen::Vector3d& vel, const Eigen::Vector3d& pos, double ts);
		t_gsim(string mark, t_gsetbase* set);

		~t_gsim();

		/**
		* @brief get motion config file
		*
		* @param[in] filename
		*
		*/
		void trj_config(const string& filename);
		void trj_config(const vector<string>& motion_strvec);


		/**
		* @brief add t_gsim::t_state into motion vector
		*
		* @param[in] mode		motion mode
		* @param[in] w			angualr
		* @param[in] a			accelarate
		* @param[in] T1			Time
		*
		*/
		void gen_mot(const MOTION_TYPE& mode, const double& T, const double& w, const double& a, const double& T1); 
		

		/**
		* @brief Inversion algorithm
		* @note generate imu data according to motion vector
		*/
		void trj_imu();


		/**
		* @brief Simulation main function
		* @note generate ultimate path and write file
		*
		* @param[in] fin		input configure file
		* @param[in] f1			output imu data file
		* @param[in] f2			output pos data file
		* @ f1 is imu data file,f2 is true trajectory file
		*
		*/
		void gen_path(t_gimudata *imu, t_gposdata *pos);

		void gen_path_from_groundtruth_test(const string& fgt, const vector<string>& fothers, t_gimudata *imu, t_gposdata *pos, int fs_gt, int fs_imu);

		void set_out();
	
	private:
		double gen_whitenoise(double mean, double stdvar);
		void gen_whitenoise(double mean, double stdvar, int len, vector<double> &data);
		void gen_markov1(double sigma, double tau, double ts, int len, vector<double>&data);

	private:
		vector<t_state> _mot;			/// motion vector
		vector<t_path> _path;			/// path
		vector<t_imu_data> _data;		/// imu data point
		map<string, MOTION_TYPE> _MODE;	/// MODE
		t_gearth _eth;					/// earth information
		t_gimu _imu;						/// imu info
		double _Vy, _g, _ts, _T;

		t_giof *_fins;
		t_gsetbase * _set;
		string _site;
	};


}



#endif