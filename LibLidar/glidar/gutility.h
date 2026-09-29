
#ifndef GUTILITY_H
#define GUTILITY_H

/**
* @file gutility.h
* @brief common algorithm
* @details implement common algorithm and define some enums

* @verbatim
		History
		-1.0	WangHuidan	       2020-10-07 create the file.
  @endverbatim

* @author WangHuidan
* @version 1.0
* @date 2020-10-07
* @license
*/


#include <string>
#include "Eigen/Eigen"
#include "opencv2/opencv.hpp"
#include <gutils/ExportLIDAR.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#ifdef USE_CERES
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#endif // USE_CERES

namespace glidar
{
	/**
	* @typedef StateIDType
	* @brief StateID
	*/
	typedef long long int StateIDType;

	/**
	* @enum Meas_TYPE
	* @brief filter update mode,included GNSS/VIS,GNSS,VIS and NO
	*/
	enum Meas_Type {All_Meas,GNSS_Meas,Lidar_Meas,No_Meas};

	/**
	* @enum Distortion_Type
	* @brief type of distortion,included radtan and equidistant
	*/
	//enum Distortion_Type { radtan, equidistant};


	/**
	* @struct lidarPath
	* @brief  save lidar time and path
	*/
	struct lidarPath
	{
		double t;           /**< time */
		string lidar_path;
	};

	/**
	* @struct dataLIDAR
	* @brief  save lidar time and its pointer object
	*/
	struct dataLIDAR
	{
		double t;
		pcl::PointCloud<pcl::PointXYZI> lidar;
	};

	/**
	* @struct PointCloud
	* @brief timestamp and featurepoints used in back end optimization
	*/
	//struct PointCloud
	//{
	//	double time;					// timestamp of current image
	//	vector<FeaturePoint> features;	// feature points set of current image
	//};

	/**
	* @struct IMUState
	* @brief some information only used in visualization
	*/
	struct IMUState
	{
		StateIDType id;						// id of imu state,not use
		double time;						// timstamp of imu state,not use
		Eigen::Quaterniond orientation;		// attitude of imu state
		Eigen::Vector3d position;			// position of imu state

		IMUState() : id(0), time(0),
			orientation(Eigen::Quaterniond::Identity()),
			position(Eigen::Vector3d::Zero()) {}

	};

	/**
	* @struct LIDARState
	* @brief store lidar information used in back end optimization
	*/
	struct LIDARState
	{
		StateIDType id;						// id of lidar state
		double time;						// time when the lidar is recorded
		Eigen::Quaterniond orientation;		// attitude of lidar state
		Eigen::Vector3d position;			// position of lidar state
		bool isKeyFrame;					// identification of current lidar

		/**
		* @brief Constructor
		* set initial parameter
		*/
		LIDARState() : id(0), time(0),
			orientation(Eigen::Quaterniond::Identity()),
			position(Eigen::Vector3d::Zero()),
			isKeyFrame(false) {}

		/**
		* @brief Constructor
		* set initial parameter by StateID
		*/
		LIDARState(const StateIDType& new_id) : id(new_id), time(0),
			orientation(Eigen::Quaterniond::Identity()),
			position(Eigen::Vector3d::Zero()),
			isKeyFrame(false) {}

	};

	/**
	* @struct Correspond corner points between different lidar frame
	* @brief  observation used in back end updating
	*/
	struct CorrespondCornerFeature
	{
		std::vector<Eigen::Vector3d> currCornerPointCloud_inlast;
		std::vector<Eigen::Vector3d> currCornerPointCloud;
		std::vector<Eigen::Vector3d> correspondCornerPointCloudA;
		std::vector<Eigen::Vector3d> correspondCornerPointCloudB;
	};

	/**
	* @struct Correspond surf points between different lidar frame
	* @brief  observation used in back end updating
	*/
	struct CorrespondSurfFeature
	{
		std::vector<Eigen::Vector3d> currSurfPointCloud_inlast;
		std::vector<Eigen::Vector3d> currSurfPointCloud;
		std::vector<Eigen::Vector3d> correspondSurfPointCloudA;
		std::vector<Eigen::Vector3d> correspondSurfPointCloudB;
		std::vector<Eigen::Vector3d> correspondSurfPointCloudC;
		std::vector<Eigen::Vector3d> norm;
		std::vector <double> negative_OA_dot_norm;
	};

	/**
	* @struct observation of lidar odometry
	* @brief  correspond points and pose between two lidar frames
	*/
	struct lidarOdometryObs
	{
		///< pose of last lidarframe
		Eigen::Matrix3d last_R_l_e;
		Eigen::Vector3d last_t_l_e;
		///< pose of current lidarframe
		Eigen::Matrix3d curr_R_l_e;
		Eigen::Vector3d curr_t_l_e;

		CorrespondCornerFeature cornerFeature;
		CorrespondSurfFeature surfFeature;

		lidarOdometryObs()
		{
			last_R_l_e = Eigen::Matrix3d::Identity();
			last_t_l_e = Eigen::Vector3d::Zero();
			curr_R_l_e = Eigen::Matrix3d::Identity();
			curr_t_l_e = Eigen::Vector3d::Zero();
		}
	};

	/**
	* @struct observation of lidar mapping
	* @brief  corespond points and pose between lidar frame and lidar map
	*/
	struct lidarMappingObs
	{
		///< pose of current lidarframe
		Eigen::Matrix3d curr_R_l_e;
		Eigen::Vector3d curr_t_l_e;
		Eigen::Matrix3d last_R_l_e;
		Eigen::Vector3d last_t_l_e;

		CorrespondCornerFeature cornerFeature;
		CorrespondSurfFeature surfFeature;

		lidarMappingObs()
		{
			curr_R_l_e = Eigen::Matrix3d::Identity();
			curr_t_l_e = Eigen::Vector3d::Zero();
			last_R_l_e = Eigen::Matrix3d::Identity();
			last_t_l_e = Eigen::Vector3d::Zero();
		}
	};

	/**
	* @struct to store information in one lidarframe
	* @brief  contain the pointclouds after classification and current pose
	*/
	struct LidarFrame
	{
		bool empty;

		pcl::PointCloud<pcl::PointXYZI> LessSurf;
		pcl::PointCloud<pcl::PointXYZI> LessSharp;
		pcl::PointCloud<pcl::PointXYZI> Surf;
		pcl::PointCloud<pcl::PointXYZI> Sharp;
		pcl::PointCloud<pcl::PointXYZI> fullCloud;

		///< store the lidar frame pose
		Eigen::Matrix3d R_l_e;
		Eigen::Vector3d t_l_e;

		LidarFrame()
		{
			empty = true;

			LessSurf.clear();
			LessSharp.clear();
			Surf.clear();
			Sharp.clear();
			fullCloud.clear();

			R_l_e = Eigen::Matrix3d::Identity();
			t_l_e = Eigen::Vector3d::Zero();
		}
	};

	/**
	* @typedef LidarStateServer
	* @brief store historical lidar information
	*/
	typedef std::map<StateIDType, LIDARState, std::less<int>,
		Eigen::aligned_allocator<std::pair<const StateIDType, LIDARState>>> LidarStateServer;

	/**
	* @brief transform a vector to skew-symmetric matrix
	* @param[in] v
	* @return Eigen::Matrix3d
	*/
	 Eigen::Matrix3d skew(const Eigen::Vector3d& v);

	 /**
	 * @brief get the rotation from navigation frame to ECEF frame 
	 * @param[in] BLH
	* @return Eigen::Matrix3d XYZ
	*/
	 Eigen::Matrix3d R_ENU_ECEF(const Eigen::Vector3d &BLH);

	 /**
	 * @brief transform 
	 * @param[in] X
	 * @return Eigen::Vector3d BLH
	 */
	 Eigen::Vector3d XYZ2BLH(const Eigen::Vector3d &X);

	 Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d& w);

	 /** \brief Calculate the squared difference of the given two points.
	 *
	 * @param a The first point.
	 * @param b The second point.
	 * @return The squared difference between point a and b.
	 */
	 template <typename PointT>
	 inline float calcSquaredDiff(const PointT& a, const PointT& b)
	 {
		 float diffX = a.x - b.x;
		 float diffY = a.y - b.y;
		 float diffZ = a.z - b.z;

		 return diffX * diffX + diffY * diffY + diffZ * diffZ;
	 }

	 /**
	 * @brief split the string by "seperator"
	 * @param[in] s		string need to be splited
	 * @param[in] seperator	symbol used for split the string.
	 */
	 vector<std::string> split(const string &s, const string &seperator);
}



#endif
