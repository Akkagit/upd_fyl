#pragma once
#ifndef LOAM_LaserMapping_H
#define LOAM_LaserMapping_H

/**
* @file glidarMapping.h
* @brief lidar scan to map algorithm
* @details
		define the t_glidarMapping class used for build lidar pointcloud map, 
		find the observation between lidar frame and the map.
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
#include <pcl/filters/voxel_grid.h>
#include "gutility.h"
#include "gset/gsetbase.h"
#include "gutils/gmutex.h"
#include "gutils/ExportLIDAR.h"



namespace glidar
{

	class LIDAR_LIBRARY_EXPORT t_glidarMapping
	{
	public:
		t_glidarMapping();

		~t_glidarMapping() {};

		/**
		* @brief add pointcloud in lidar frame to the lidarmap 
		*
		* @param[in] frame	lidar frame(with the pointcloud after distortion correction)
		*/
		void addPointcloudToMap(LidarFrame frame);

		/**
		* @brief lidar processing main function
		*
		* @param[in] frame	lidar frame(with the pointcloud after distortion correction)
		*/
		void process(LidarFrame frame,bool lastmap);

		//use ceres to optimize the problem
		//void ceresOptimize();
		
		//clear the history information
		void reset();

		/**
		* @brief transform point from w frame to e frame
		*
		* @param[in] point_w	coordinate of point in w frame
		*/
		Eigen::Vector3d translateToE(Eigen::Vector3d point_w);

		/**
		* @brief renew the transfrom from w frame to e frame
		*
		* @param[in] R t	new transfrom from w frame to e frame
		*/
		void renewWframe(Eigen::Matrix3d R, Eigen::Vector3d t);

	protected:
		///< position id of last lidarframe center
		int lastCenterCubeI;
		int lastCenterCubeJ;
		int lastCenterCubeK;

		/**
		* @brief transform point from w frame to e frame
		*
		* @param[in] point_w	coordinate of point in w frame
		*/
		void transformAssociateToMap();

		/**
		* @brief transform point from lidar frame to w frame
		*
		* @param[in] pi		coordinate of point in lidar frame
		* @param[out] po	coordinate of point in w frame
		* @param[in] lastmap	if we choose to transform the map to last lidar frame
		*/
		void pointAssociateToMap(pcl::PointXYZI & pi, pcl::PointXYZI & po,bool lastmap);

		/**
		* @brief transform the map points from old w frame to new w frame
		*
		* @param[in] R_tonewW	rotation from old to new W frame
		* @param[in] t_tonewW	translation from old to new W frame 
		*/
		void toNewW(Eigen::Matrix3d R_tonewW, Eigen::Vector3d t_tonewW);

		/**
		* @brief add corner points to the map
		*
		* @param[in] cornerPoints	coordinates of points in lidar frame
		*/
		void AddCornerPointToMap(pcl::PointCloud<pcl::PointXYZI> cornerPoints);

		/**
		* @brief add surf points to the map
		*
		* @param[in] surfPoints		coordinates of points in lidar frame
		*/
		void AddSurfPointToMap(pcl::PointCloud<pcl::PointXYZI> surfPoints);

		/**
		* @brief downsize the lidar map
		*/
		void downSizeMap();

		/**
		* @brief zoom the lidar map while the first lidar frame arrive
		*
		* @param[in] centerCubeIJK	position of center point of the first lidar frame
		*/
		void firstZoom(int centerCubeI, int centerCubeJ, int centerCubeK);

		/**
		* @brief zoom the lidar map(not the first frame) to asure there has 5*5*5 
		*  adjacent cubes around the lidar frame center.
		*
		* @param[in] centerCubeIJK	position of center point of the lidar frame
		*/
		void Zoom(int centerCubeI, int centerCubeJ, int centerCubeK);

		/**
		* @brief find correspond submap for current lidar frame
		*
		* @param[in] centerCubeIJK	position of center point of the lidar frame
		*/
		void FindSubmap(int centerCubeI, int centerCubeJ, int centerCubeK);

		///< tools:change the type to store the pointcloud position
		vector<Eigen::Vector3d> PointXYZI2Vector3d(pcl::PointCloud<pcl::PointXYZI>::Ptr xyzi);
		
		///< tools:transform the PointXYZI
		void transPointXYZI(pcl::PointXYZI &p, Eigen::Matrix3d& R, Eigen::Vector3d& t);

		///< tools:transform the PointXYZI
		void transPointXYZI(pcl::PointCloud<pcl::PointXYZI>::Ptr xyzi, Eigen::Matrix3d& R, Eigen::Vector3d& t);

	public:

		bool newCornerPointsSharp_;       ///< flag if a new sharp corner cloud has been received
		bool newCornerPointsLessSharp_;   ///< flag if a new less sharp corner cloud has been received
		bool newSurfPointsFlat_;          ///< flag if a new flat surface cloud has been received
		bool newSurfPointsLessFlat_;      ///< flag if a new less flat surface cloud has been received
		bool newLaserCloudFullRes_;       ///< flag if a new full resolution cloud has been received
		bool newImuTrans_;                ///< flag if a new IMU transformation information cloud has been received

		float scanPeriod_;       ///< time per scan
		bool systemInited_;      ///< initialization flag

		///< store the pose of first lidar frame(used as the transform between w frame and e frame)
		Eigen::Matrix3d first_R_l_e;
		Eigen::Vector3d first_t_l_e;
		///< pose of current lidar frame
		Eigen::Matrix3d curr_R_l_e;
	    Eigen::Vector3d curr_t_l_e;
		Eigen::Matrix3d curr_R_l_w;
		Eigen::Vector3d curr_t_l_w;
		Eigen::Matrix3d curr_R_w_l;
		Eigen::Vector3d curr_t_w_l;

		///< pose of last lidar frame
		Eigen::Matrix3d last_R_l_e;
		Eigen::Vector3d last_t_l_e;
		Eigen::Matrix3d last_R_w_l;
		Eigen::Vector3d last_t_w_l;
		
		///< transfrom from current frame to last frame
		Eigen::Matrix3d curr_last_R;
		Eigen::Vector3d curr_last_t;

		///< store the submap points used in last frame(used for visibility)
		vector<Eigen::Vector3d> submap_corner;
		vector<Eigen::Vector3d> submap_surf;

		///< store the observation of lidar mapping
		CorrespondCornerFeature correspondCornerFeature_;
		CorrespondSurfFeature   correspondSurfFeature_;
		lidarMappingObs lidarMapObs;

	///< @parameters:used for construct the lidar pointcloud map
	public:
		///< divide the spatical grid
		int frameCount = 0;

		int laserCloudCenWidth;
		int laserCloudCenHeight;
		int laserCloudCenDepth;

		int laserCloudWidth;
		int laserCloudHeight;
		int laserCloudDepth;

		int laserCloudNum; ///< 4851

		int laserCloudValidNum;
		int laserCloudSurroundNum;
		int laserCloudValidInd[125];
		int laserCloudSurroundInd[125];


		///< surround points in submap to build kdtree
		pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudCornerFromMap;
		pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudSurfFromMap;

		///< points in every cube
		pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudCornerArray[4851];
		pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudSurfArray[4851];

		///< kd-tree
		pcl::KdTreeFLANN<pcl::PointXYZI>::Ptr kdtreeCornerFromMap;
		pcl::KdTreeFLANN<pcl::PointXYZI>::Ptr kdtreeSurfFromMap;

		///< use for downsize the pointcloud
		pcl::VoxelGrid<pcl::PointXYZI> downSizeFilterCorner;
		pcl::VoxelGrid<pcl::PointXYZI> downSizeFilterSurf;
		pcl::VoxelGrid<pcl::PointXYZI> downSizeMapCorner;
		pcl::VoxelGrid<pcl::PointXYZI> downSizeMapSurf;

		///< temp paramter
		std::vector<int> pointSearchInd;
		std::vector<float> pointSearchSqDis;

		pcl::PointXYZI pointOri, pointSel;

	};


}//end of namespace loam

#endif
