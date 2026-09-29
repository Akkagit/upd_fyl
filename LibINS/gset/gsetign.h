#ifndef GSETIGN_H
#define GSETIGN_H

/**
* @file gsetign.h
* @brief implements settings for all coupled integration, including lci, tci, rtci.
* @author ShenZhiheng
* @email
* @version 1.0
* @date 2019-12-27
* @license
*/


#define XMLKEY_IGN "integration"

#include "gset/gsetbase.h"
#include "gins/gutility.h"
using namespace gnut;

namespace gsins
{
	/**
	* @class t_gsetign
	* @brief coupled integration settings
	* @author ShenZhiheng
	*/
	class INS_LIBRARY_EXPORT t_gsetign : public virtual t_gsetbase
	{
	public:
		t_gsetign();
		~t_gsetign();

		void check();
		void help();

		Eigen::Vector3d initial_misalignment_std();
		Eigen::Vector3d initial_vel_std();
		Eigen::Vector3d initial_pos_std();
		Eigen::Vector3d initial_gyro_std();
		Eigen::Vector3d initial_acce_std();

		Eigen::Vector3d min_misalignment_std();
		Eigen::Vector3d min_vel_std();
		Eigen::Vector3d min_pos_std();
		Eigen::Vector3d min_gyro_std();
		Eigen::Vector3d min_acce_std();

		Eigen::Vector3d misalignment_psd();
		Eigen::Vector3d vel_psd();
		Eigen::Vector3d pos_psd();
		Eigen::Vector3d gyro_psd();
		Eigen::Vector3d acce_psd();

		Eigen::Vector3d meas_vel_noise();
		Eigen::Vector3d meas_pos_noise();

		string order();

		int nq();
		int nr();
		double TS();
		double delay_t();
		double max_pdop();
		int    min_sat();
		Eigen::Vector3d lever();
		IGN_TYPE ign_type();
		gsins::FLT_TYPE fltmode();
		gsins::LCI_TYPE gnssmode();
	protected:
		Eigen::Vector3d _initial_misalignment_std;       /**< Standard deviation of initial misalignment angle [deg]. */
		Eigen::Vector3d _initial_vel_std;				 /**< Standard deviation of initial velocity [m/s]. */
		Eigen::Vector3d _initial_pos_std;				 /**< Standard deviation of initial position [m]. */
		Eigen::Vector3d _initial_gyro_std;				 /**< Standard deviation of initial gyro drift [deg/h]. */
		Eigen::Vector3d _initial_acce_std;				 /**< Standard deviation of initial accelerator bias [mg]. */

		Eigen::Vector3d _min_misalignment_std;			 /**< Standard deviation lower bound of misalignment angle [deg]. */
		Eigen::Vector3d _min_vel_std;					 /**< Standard deviation lower bound of velocity [m/s]. */
		Eigen::Vector3d _min_pos_std;					 /**< Standard deviation lower bound of position [m]. */
		Eigen::Vector3d _min_gyro_std;					 /**< Standard deviation lower bound of gyro drift [deg/h]. */
		Eigen::Vector3d _min_acce_std;					 /**< Standard deviation lower bound of accelerator bias [mg]. */

		Eigen::Vector3d _misalignment_psd;				 /**< Spectral density of misalignment angle [deg/sqrt(h)]. */
		Eigen::Vector3d _vel_psd;						 /**< Spectral density of velocity [mg/sqrt(Hz)]. */
		Eigen::Vector3d _pos_psd;						  
		Eigen::Vector3d _gyro_psd;
		Eigen::Vector3d _acce_psd;

		Eigen::Vector3d _meas_vel_noise;				 /// Measurement noise of velocity [m/s]
		Eigen::Vector3d _meas_pos_noise;			     /// Measurement noise of position [m]

		Eigen::Vector3d _lever;							 /// GNSS Receiver center's lever relative to imu center on b Coor
		int _nq, _nr;									 /// parameter dimension and measurement dimension
		gsins::FLT_TYPE  _fltmode;						 /// filter mode (forward, forward and backward smooth, RTS)
		double _TS;										 /// measurement data interval
		string _order;									 /// the order for pos file (size is 6)
		double _delay_t;								 /// the GNSS misalignment time diff with imu
		double _max_pdop;								 /// the max PDOP threshold for LCI
		int    _min_sat;								 /// the minimum sat num threshold.
		gsins::LCI_TYPE _gnssmode;	
		IGN_TYPE _ign_type;

	};
}


#endif