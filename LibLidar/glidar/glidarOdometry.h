#pragma once
#ifndef LOAM_LASERODOMETRY_H
#define LOAM_LASERODOMETRY_H

/**
* @file glidarOdometry.h
* @brief lidar scan to scan algorithm
* @details
		define the t_glidarOdometry class used to find the observation
		between the adjacent lidar frames.
* @verbatim
		History
		-1.0	WangHuidan	       2020-10-09 create the file.
  @endverbatim

* @author Wanghuidan
* @email
* @version 1.0
* @date 2020-10-09
* @license
*/

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/features/normal_3d.h>
#include <pcl/registration/icp.h>
#include <pcl/registration/transformation_estimation_point_to_plane.h>
#include <fstream>
#include <time.h>
#include "gutility.h"
#include "gset/gsetbase.h"
#include "gutils/gmutex.h"
#include "gutils/ExportLIDAR.h"

namespace glidar
{
	/**
	 * @brief Implementation of the LOAM laser odometry component.
	 */
	class LIDAR_LIBRARY_EXPORT t_glidarOdometry
	{
	public:
		t_glidarOdometry();

		~t_glidarOdometry() {};

		/**
		* @brief main function to find correspond "scan to scan" observation
		*
		* @param[in] frame1		last lidar frame
		* @param[in] frame2		current lidar frame
		*/
		void process(LidarFrame& frame1, LidarFrame& frame2);

		/**
		* @brief remove the distortion of pointcloud in current lidar frame
		*
		* @param[in] frame1		last lidar frame
		* @param[in] frame2		current lidar frame
		*/
		void removeDistortion(LidarFrame& frame1, LidarFrame& frame2);

		//ceres optimize,need to find how to add ceres lib to this solution
		//void ceresOptimize();

		///< tools:change the type to store the pointcloud position, transform the position to e frame at the same time
		vector<Eigen::Vector3d> PointXYZI2Vector3d(pcl::PointCloud<pcl::PointXYZI>& xyzi,Eigen::Matrix3d R_l_e,Eigen::Vector3d t_l_e);

	protected:
		/** @brief	remove the history information, etc. */
		void reset();

		/** 
		 * @brief Transform the given point to the start of the sweep.
		 * 
		 * @param pi the point to transform
		 * @param po the point instance for storing the result
		 */
		void transformToStart(const pcl::PointXYZI pi, pcl::PointXYZI & po,bool distortion=false);

		/**
		 * @brief Transform the given point to the end of the sweep.
		 *
		 * @param pi the point to transform
		 * @param po the point instance for storing the result
		 */
		void transformToEnd(const pcl::PointXYZI pi, pcl::PointXYZI & po,bool distortion=false);

	public:
		///< transformation from last frame to current frame
		Eigen::Matrix3d last_curr_rot;
		Eigen::Vector3d last_curr_trans;
		///< transfprmation from current frame to last frame
		Eigen::Matrix3d curr_last_rot;
		Eigen::Vector3d curr_last_trans;

		float scanPeriod_;       ///< time per scan
		long frameCount_;        ///< number of processed frames
		size_t maxIterations_;   ///< maximum number of iterations
		bool systemInited_;      ///< initialization flag

		///< for the publish of pointcloud
		vector<Eigen::Vector3d> test_lastcorner;
		vector<Eigen::Vector3d> test_lastsurf;

		pcl::KdTreeFLANN<pcl::PointXYZI> lastCornerKDTree_;   ///< last corner cloud KD-tree
		pcl::KdTreeFLANN<pcl::PointXYZI> lastSurfaceKDTree_;  ///< last surface cloud KD-tree

		CorrespondCornerFeature correspondCornerFeature_;
		CorrespondSurfFeature   correspondSurfFeature_;
		lidarOdometryObs lidarOdoObs;

	};


}///< end of namespace glidar

#endif 
