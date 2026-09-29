#ifndef GLIDAR_H
#define GLIDAR_H

/**
* @file glidar.h
* @brief lidar core algorithm
* @details
		Starting from the receiving pointcloud data, including classification, 
		lidar odometry, lidar mapping, find the corresponding observations in 
		each part for measurement update
* @verbatim
		History
		-1.0	WangHuidan	       2020-10-07 create the file.
  @endverbatim

* @author WangHuidan
* @email
* @version 1.0
* @date 2020-10-07
* @license
*/

#include "gutils/gmutex.h"
#include "gutils/ExportLIDAR.h"
#include "gutility.h"
#include "glidarOdometry.h"
#include "glidarMapping.h"
#include "glidarproc.h"
#include "glidarframe.h"
#include "gpublish.h"

using namespace gnut;
using namespace std;

namespace glidar
{
	
	class LIDAR_LIBRARY_EXPORT t_glidar
	{
	public:
		/**
		* @brief constructor.
		*
		* @param[in] set setbase control
		*/
		t_glidar(t_gsetbase* set);
		
		/**
		* @brief destructor
		*/
		~t_glidar();
		
		/**
		* @brief screen out the observations whose residual is larger than n times of the MSE
		* @para[in]	H	jacobian matrix of the lidar residual function
		* @para[in] r	residual vector of the lidar residual function
		* @para[in] n	to identify the dimension of observations
		*/
		void screenOut(Eigen::MatrixXd &H, Eigen::VectorXd &r, int n);
		
		/**
		* @brief screen out the observations whose residual is larger than n times of the MSE
		* @para[in]	H	jacobian matrix of the lidar residual function
		* @para[in] r	residual vector of the lidar residual function
		* @para[in] n	to identify the threshold of this function
		*/
		void Huber(Eigen::MatrixXd &H, Eigen::VectorXd &r, float ther);
		
		/**
		* @brief calculate the measurement jacobian for the lidar observations
		* @para[in] obs	observation of the lidar measurement
		* @para[in]	H	jacobian matrix of the lidar residual function
		* @para[in] r	residual vector of the lidar residual function
		* @para[in] use_3d	true: H of corner observation will consider 3D directions.
		* @para[in] use_3d	false: H of corner observation will consider only one direction.
		* @para[in] pursue	a value for Huber function as its' thershold "ther"
		*/
		void lidarMeasurementJacobian(lidarMappingObs obs, Eigen::MatrixXd &H, Eigen::VectorXd &r,bool use_3d, float ther);
		void lidarMeasurementJacobian(lidarOdometryObs obs, Eigen::MatrixXd &H, Eigen::VectorXd &r, bool use_3d, float ther);

		void lidarMeasurementJacobian_test(lidarMappingObs obs, Eigen::MatrixXd &H, Eigen::VectorXd &r, bool use_3d, float ther);



		
	public:
		Eigen::Isometry3d T_lidar_imu;
		Eigen::Vector3d t_lidar_imu;
		Eigen::Matrix3d R_lidar_imu;

		bool use_scan;
		bool use_map;
		bool proj_lastmap;
		bool use_corrdistort;
		int map_gap;
		
		Eigen::MatrixXd Phi;
		Eigen::MatrixXd P_LL;
		Eigen::MatrixXd P_IL;
		Eigen::MatrixXd P_IL_sav;
		Eigen::MatrixXd P_LL_sav;
		Eigen::MatrixXd Dx_sav;
		int Dx_row;

		vector<StateIDType> FrameBuffer;			// store lidarframes
		
		LidarStateServer lidar_states;
		IMUState imu_state;
		StateIDType imu_state_id;
		StateIDType next_id;

		//transformation parameters in lidarmapping from w frame to e frame
		Eigen::Matrix3d R_w_e;
		Eigen::Vector3d t_w_e;
		
		int max_lidarstate_size;
		double observation_noise;
		int num_of_lidar;
		bool residual_threshold;
		string output_path;

		bool first_odo=false;
		bool first_map = true;
		
	protected:
		t_glidarproc *lidarproc;
		t_gpublish publisher;
	
	};
	
}

#endif