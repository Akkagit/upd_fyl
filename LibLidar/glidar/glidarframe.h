

#ifndef GLIDARFRAME_H
#define GLIDARFRAME_H

/**
* @file glidarframe.h
* @brief frame core algorithm
* @details
		LIDAR define the t_glidarframe class used for remove some poor points and classified the pointcloud
		by their curvature and scan ID.
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
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/filter.h>
#include "gutility.h"
#include "gset/gsetlidar.h"
#include "gutils/gmutex.h"
#include "glidarproc.h"
#include "gutils/ExportLIDAR.h"
#include <Eigen/Eigen>
using namespace gnut;


typedef std::pair<size_t, size_t> IndexRange;

namespace glidar
{ 
	//@brief Class for registering point clouds received from multi-laser lidars.
	class LIDAR_LIBRARY_EXPORT t_glidarframe:public t_glidarproc
	{
	public:
		t_glidarframe(t_gsetbase* _set);
		~t_glidarframe() {};

		/**
		* @brief lidar frame main function
		*
		* @param[in] t		 cur lidar time
		*/
		virtual LidarFrame ProcessBatch(const double& t);

		/**
		* @brief read lidar pointcloud from bin file
		*
		* @param[in]in_file		lidar file path
		*/
		pcl::PointCloud<pcl::PointXYZI> readKittiBinData(std::string &in_file);

		/**
		* @brief read lidar pointcloud from txt file
		*
		* @param[in] in_file	lidar file path
		*/
		pcl::PointCloud<pcl::PointXYZI> readKittiTxtData(std::string &filepath);

		/**
		* @brief classification of the lidar pointcloud
		*
		* @param[in] laserCloudIn	lidar pointcloud
		*/
		void handleCloud(pcl::PointCloud<pcl::PointXYZI>& laserCloudIn);

		/**
		* @brief remove the close points in lidar pointcloud
		*
		* @param[in] cloud_in	lidar pointcloud
		* @param[in] thres		distance thershold of the close point
		* @param[out] cloud_out	lidar pointcloud
		*/
		void removeClosedPointCloud(pcl::PointCloud<pcl::PointXYZI> &cloud_in,
			pcl::PointCloud<pcl::PointXYZI> &cloud_out, float thres);

		//clear the history information
		void clear();

	public:
		void process(pcl::PointCloud<pcl::PointXYZI>& laserCloudIn/*, const Time& scanTime*/);

		std::vector<pcl::PointCloud<pcl::PointXYZI> > _laserCloudScans;

		pcl::PointCloud<pcl::PointXYZI> _laserCloud;   ///< full resolution input cloud
		std::vector<IndexRange> _scanIndices;          ///< start and end indices of the individual scans withing the full resolution cloud

		pcl::PointCloud<pcl::PointXYZI> _cornerPointsSharp;      ///< sharp corner points cloud
		pcl::PointCloud<pcl::PointXYZI> _cornerPointsLessSharp;  ///< less sharp corner points cloud
		pcl::PointCloud<pcl::PointXYZI> _surfacePointsFlat;      ///< flat surface points cloud
		pcl::PointCloud<pcl::PointXYZI> _surfacePointsLessFlat;

		int n_scans;
		double close_threshold;
		double scanPeriod;
		double time;							// timestamp of the lidar
		dataLIDAR lidar_msg;					// original observation data

	};

}//end of namespace glidar

#endif

