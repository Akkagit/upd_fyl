#ifndef GPUBLISH_H
#define GPUBLISH_H

#pragma once
#include "gviewer.h"
#include "gutility.h"
#include "gutils/ExportLIDAR.h"
#include <Eigen/Eigen>
#include <fstream>

namespace glidar
{
	/**
	* @class t_gpublish
	* @brief t_gpublish Class is used for trajectory display
	*/
	class LIDAR_LIBRARY_EXPORT t_gpublish
	{
	public:
		/** @brief default constructor. */
		t_gpublish();

		/** @brief default destructor. */
		~t_gpublish();

		/**
		* @brief initialize the viewer for display
		*/
		void Initialize(string path);

		/**
		* @brief update the parameters for trajectory display
		*
		* @param[in] imu_state		store the imu information to show on the viewer
		* @param[in] lidar_states	store the lidar information to show on the viewer
		* @param[in] init_campos	the initial value of the camera position
		*/
		void UpdateNewState(const IMUState &imu_state, const LidarStateServer &lidar_states, Eigen::Vector3d init_campos=Eigen::Vector3d::Zero(), vector<Eigen::Vector3d> lidarPointcloud= vector<Eigen::Vector3d>());

	private:
		ofstream outputFile;
		gviewer *viewer;				// the core class of trajectory display
		bool firstFlag_imu = true;		// A indicator to determine the value of _init_imupos
		bool firstFlag_lidar = true;	// A indicator to determine the value of _init_campos,
										// maybe imgdata is later than imudata
		Eigen::Vector3d _init_lidarpos;	// transfer the zero of the IMU coordinate
		Eigen::Vector3d _init_imupos;	// transfer the zero of the IMU coordinate
		Eigen::Matrix3d R_e_n;          // rotation from ECEF frame to ENU frame for first imudata
	};

}

#endif