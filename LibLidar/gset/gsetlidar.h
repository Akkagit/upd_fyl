#ifdef USE_LIDAR#ifdef USE_LIDAR#ifdef USE_LIDAR
#ifndef GSETLIDAR_H
#define GSETLIDAR_H

/**
* @file gsetlidar.h
* @brief implements settings for lidar
* @verbatim
		History
		-1.0	WangHuidan	       2020-09-29 create the file.
  @endverbatim

* @author WangHuidan
* @email
* @version 1.0
* @date 2020-09-29
* @license
*/

#define XMLKEY_LIDAR "lidar"
#include "glidar/gutility.h"
#include "gset/gsetbase.h"
#include "gutils/ExportLIDAR.h"
using namespace gnut;

namespace glidar
{
	/**
	* @class t_gsetlidar
	*
	* @brief t_gsetlidar Class provides some parameters used in back end.
	*
	*
	*/
	class LIDAR_LIBRARY_EXPORT t_gsetlidar : public virtual t_gsetbase
	{
	public:
		/** @brief constructor. */
		t_gsetlidar();

		/** @brief default constructor. */
		~t_gsetlidar();

		/**
		* @brief empty function
		*/
		void check();

		/**
		* @brief empty function
		*/
		void help();

		/**
		* @brief get lidardata interval
		* @return double
		*/
		double ts();

		/**
		* @brief get imudata interval
		* @return double
		*/
		double imu_ts();

		/**
		* @brief get lidardata frequency
		* @return int
		*/
		int freq();

		/**
		* @brief get imudata frequency
		* @return int
		*/
		int imu_freq();

		/**
		* @brief get rotation from lidar to imu
		* @return Eigen::Matrix3d
		*/
		Eigen::Matrix3d R_lidar_imu();

		/**
		* @brief get translation from lidar to imu
		* @return Eigen::Vector3d
		*/
		Eigen::Vector3d t_lidar_imu();

		/**
		* @brief get transformation from lidar to imu
		* @return Eigen::Isometry3d
		*/
		Eigen::Isometry3d T_lidar_imu();
		
		/*Lidar Parameters*/
		int n_scans();
		double scan_period();
		double map_resolution_sharp();
		double map_resolution_surf();
		int max_iteration();
		bool use_scan();
		bool use_map();
		bool distortion();
		bool proj_lastmap();
		int map_gap();
		/*Lidar Parameters*/

		/*Estimation Parameters*/
		
		int max_lidar_state_size();
		double close_threshold();
		double surf_dis_threshold();
		double sharp_dis_threshold();
		double screenout_threshold();
		double lidar_observation_noise();
		double huber_parameter();
		/*Estimation Parameters*/

		string output_path();

		/*Optimization Parameters*/
		
		//double estimation_precision();
		/*Optimization Parameters*/

	protected:

	
	};
}
#endif
#endif // USE_LIDAR