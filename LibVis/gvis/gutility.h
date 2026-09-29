
#ifndef GUTILITY_H
#define GUTILITY_H

/**
* @file gutility.h
* @brief common algorithm
* @details implement common algorithm and define some enums

* @verbatim
		History
		-1.0	LiShengyu	       2020-05-20 create the file.
  @endverbatim

* @author LiShengyu
* @version 1.0
* @date 2020-05-20
* @license
*/


#include <string>
#include "Eigen/Eigen"
#include <gutils/ExportVIS.h>
#include "opencv2/opencv.hpp"
#ifdef USE_CERES
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#endif // USE_CERES

namespace gvis
{
	/**
	* @typedef StateIDType
	* @brief StateID
	*/
	typedef long long int StateIDType;

	/**
	* @typedef FeatureIDType
	* @brief FeatureID
	*/
	typedef long long int FeatureIDType;


	/**
	* @enum Meas_TYPE
	* @brief filter update mode,included GNSS/VIS,GNSS,VIS and NO
	*/
	enum Meas_Type {All_Meas,GNSS_Meas,Vis_Meas,No_Meas};

	/**
	* @enum Distortion_Type
	* @brief type of distortion,included radtan and equidistant
	*/
	enum Distortion_Type { radtan, equidistant};


	/**
	* @struct imgPath
	* @brief  save image time and path
	*/
	struct imgPath
	{
		double t;           /**< time */
		string img0_path;
		string img1_path;
	};

	/**
	* @struct dataIMG
	* @brief  save image time and its pointer object
	*/
	struct dataIMG
	{
		double t;
		shared_ptr<cv::Mat> img0;
		shared_ptr<cv::Mat> img1;
	};

	/**
	* @struct FeaturePoint
	* @brief  featurepoint used in front end tracking
	*/
	struct FeaturePoint
	{
		long long int id;			// feature id             
		cv::Point2f cam0_point;		// feature pixel coordinate in left image 
		cv::Point2f cam1_point;		// feature pixel coordinate in rgiht image
		double depth;               // feature initial depth in left image,only used in stereo_orb model
		float response;             // feature response
		int lifetime;               // feature lifetime
	};

	/**
	* @struct PointCloud
	* @brief timestamp and featurepoints used in back end optimization
	*/
	struct PointCloud
	{
		double time;					// timestamp of current image
		vector<FeaturePoint> features;	// feature points set of current image
	};

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
			orientation(Eigen:: Quaterniond::Identity()),
			position(Eigen::Vector3d::Zero()) {}
		
	};

	/**
	* @struct CAMState
	* @brief store camera information used in back end optimization
	*/
	struct CAMState
	{
		StateIDType id;						// id of camera state
		double time;						// time when the image is recorded
		Eigen::Quaterniond orientation;		// attitude of camera state
		Eigen::Vector3d position;			// position of camera state
		bool isKeyFrame;					// identification of current image
		double mtracking_rate;				// feature tracking rate of current image

		/**
		* @brief Constructor
		* set initial parameter
		*/
		CAMState() : id(0), time(0),
			orientation(Eigen::Quaterniond::Identity()),
			position(Eigen::Vector3d::Zero()),
			mtracking_rate(0),
			isKeyFrame(false) {}

		/**
		* @brief Constructor
		* set initial parameter by StateID
		*/
		CAMState(const StateIDType& new_id) : id(new_id), time(0),
			orientation(Eigen::Quaterniond::Identity()),
			position(Eigen::Vector3d::Zero()),
			mtracking_rate(0),
			isKeyFrame(false) {}

	};

	/**
	* @typedef CamStateServer
	* @brief store historical camera information
	*/
	typedef std::map<StateIDType, CAMState, std::less<int>,
		Eigen::aligned_allocator<std::pair<const StateIDType, CAMState>>> CamStateServer;

	/**
	* @struct Feature_param
	* @brief store some parameters used in depth optimization,only used in the initialization of t_gfeature
	*/
	struct Feature_param
	{
		double translation_threshold;			// transformation matrix threshold
		double huber_epsilon;					// kernel function
		double estimation_precision;			// threshold of whether the update quantity is iterative
		double	initial_damping;				// initial lambda of Levenberg-Marquart
		int outler_loop_max_iteration;			// iterations
		int inner_loop_max_iteration;			// iterations of depth
		Eigen::Isometry3d T_cam0_cam1;			// transformation matrix between cam0 and cam1
		bool stereo;							// stereo or mono

	};

	#ifdef USE_CERES
	struct ReprojectionError
	{
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

			//  // ReprojectionError(const cv::Point2f &measure_point_uv,
			//  //                   const Eigen::Isometry3d &cam0_cami_transform) : image_uv(measure_point_uv),
			//  //                                                                   cam0_cami(cam0_cami_transform) {}
			//  // template <typename T>
			//  // bool operator()(const T *const feature_point_world, T *residuals) const
			//  // {
			//  //     Eigen::Matrix<T, 3, 1> point;
			//  //     point << feature_point_world[0], feature_point_world[1], feature_point_world[2];//

			//  //     point = cam0_cami.inverse() * point;
			//  //     T point_reproject_image[2] = {point(0) / point(2), point(1) / point(2)};
			//  //     residuals[0] = point_reproject_image[0] - T(image_uv.x);
			//  //     residuals[1] = point_reproject_image[1] - T(image_uv.y);
			//  //     return true;
			//  // }

			ReprojectionError(const cv::Point2f &measure_point_uv,
				const Eigen::Isometry3d &cam0_cami_transform) : image_uv(measure_point_uv)
		{
			//Eigen::Isometry3d cami_cam0 = cam0_cami_transform.inverse();
			
			quat = Eigen::Quaterniond(cam0_cami_transform.rotation());
				//utiltool::attitude::RotationMartix2Quaternion(cami_cam0.rotation());
			cam0_cami_translation[0] = cam0_cami_transform.translation()(0);
			cam0_cami_translation[1] = cam0_cami_transform.translation()(1);
			cam0_cami_translation[2] = cam0_cami_transform.translation()(2);
		}
		template <typename T>
		bool operator()(const T *const feature_point_world, T *residuals) const
		{
			T point[3] = { T(feature_point_world[0]), T(feature_point_world[1]), T(1.0) };
			T point_result[3];
			T cam0_cami_quat[4] = { T(quat.w()),
								   T(quat.x()),
								   T(quat.y()),
								   T(quat.z()) };
			T point_translation[3] = { T(feature_point_world[2] * cam0_cami_translation[0]),
									  T(feature_point_world[2] * cam0_cami_translation[1]),
									  T(feature_point_world[2] * cam0_cami_translation[2]) };
			ceres::QuaternionRotatePoint(cam0_cami_quat, point, point_result);
			point_result[0] += point_translation[0];
			point_result[1] += point_translation[1];
			point_result[2] += point_translation[2];

			T point_reproject_image[2] = { T(point_result[0] / point_result[2]),
										  T(point_result[1] / point_result[2]) };
			residuals[0] = -point_reproject_image[0] + T(image_uv.x);
			residuals[1] = -point_reproject_image[1] + T(image_uv.y);
			//cout << "residual:" << residuals[0] << "," << residuals[1] << endl;
			return true;
		}

		static ceres::CostFunction *Create(const cv::Point2f &measure_point_uv,
			const Eigen::Isometry3d &cam0_cami_transform)
		{
			return new ceres::AutoDiffCostFunction<ReprojectionError, 2, 3>(new ReprojectionError(measure_point_uv, cam0_cami_transform));
		}

		const cv::Point2f image_uv;
		Eigen::Quaterniond quat;
		double cam0_cami_translation[3];
		// const Eigen::Isometry3d cam0_cami;
	};
	#endif // USE_CERES

	/** @brief string to Distortion_Type. */
	Distortion_Type str2distortion(const string& s);
	
	/** @brief Distortion_Type to string. */
	string distortion2str(const Distortion_Type& s);


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

	

}



#endif
