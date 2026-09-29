#ifdef USE_VIS#ifdef USE_VIS#ifdef USE_VIS
#ifndef GSETVIS_H
#define GSETVIS_H

/**
* @file gsetvis.h
* @brief implements settings for vis
* @verbatim
		History
		-1.0	LiShengyu	       2020-05-20 create the file.
  @endverbatim

* @author LiShengyu
* @email
* @version 1.0
* @date 2020-05-20
* @license
*/

#define XMLKEY_VIS "vis"
#include "gvis/gutility.h"
#include "gset/gsetbase.h"
#include "gutils/ExportVIS.h"
using namespace gnut;

namespace gvis
{
	/**
	* @class t_gsetvis
	*
	* @brief t_gsetvis Class provides some parameters used in back end.
	*
	*
	*/
	class VIS_LIBRARY_EXPORT t_gsetvis : public virtual t_gsetbase
	{
	public:
		/** @brief constructor. */
		t_gsetvis();

		/** @brief default constructor. */
		~t_gsetvis();

		/**
		* @brief empty function
		*/
		void check();

		/**
		* @brief empty function
		*/
		void help();

		/**
		* @brief get imgdata interval
		* @return double
		*/
		double ts();

		/**
		* @brief get imudata interval
		* @return double
		*/
		double imu_ts();

		/**
		* @brief get imgdata frequency
		* @return int
		*/
		int freq();

		/**
		* @brief get imudata frequency
		* @return int
		*/
		int imu_freq();

		/**
		* @brief get img0 distortion model,ratan or equidistant
		* @return Distortion_Type
		*/
		gvis::Distortion_Type cam0_distortion_model();

		/**
		* @brief get img1 distortion model,ratan or equidistant
		* @return Distortion_Type
		*/
		gvis::Distortion_Type cam1_distortion_model();

		/**
		* @brief get cam0 resolution
		* @return Eigen::Vector2d
		*/
		Eigen::Vector2d cam0_resolution();

		/**
		* @brief get cam1 resolution
		* @return Eigen::Vector2d
		*/
		Eigen::Vector2d cam1_resolution();

		/**
		* @brief get cam0 intrinsics
		* @return Eigen::Vector4d
		*/
		Eigen::Vector4d cam0_intrinsics();

		/**
		* @brief get cam1 intrinsics
		* @return Eigen::Vector4d
		*/
		Eigen::Vector4d cam1_intrinsics();

		/**
		* @brief get cam0 distortion coeffs
		* @return Eigen::Vector4d
		*/
		Eigen::Vector4d cam0_distortion_coeffs();

		/**
		* @brief get cam1 distortion coeffs
		* @return Eigen::Vector4d
		*/
		Eigen::Vector4d cam1_distortion_coeffs();

		/**
		* @brief get rotation from cam0 to cam1
		* @return Eigen::Matrix3d
		*/
		Eigen::Matrix3d R_cam0_cam1();

		/**
		* @brief get translation from cam0 to cam1
		* @return Eigen::Vector3d
		*/
		Eigen::Vector3d t_cam0_cam1();

		/**
		* @brief get transformation from cam0 to cam1
		* @return Eigen::Isometry3d
		*/
		Eigen::Isometry3d T_cam0_cam1();

		/**
		* @brief get rotation from cam0 to imu
		* @return Eigen::Matrix3d
		*/
		Eigen::Matrix3d R_cam0_imu();

		/**
		* @brief get translation from cam0 to imu
		* @return Eigen::Vector3d
		*/
		Eigen::Vector3d t_cam0_imu();

		/**
		* @brief get transformation from cam0 to imu
		* @return Eigen::Isometry3d
		*/
		Eigen::Isometry3d T_cam0_imu();

		/**
		* @brief get number of camera 1 or 2
		* @return camera number
		*/
		int num_of_cam();
		
		/*Stereo Optical Flow Parameters*/
		int grid_row();
		int grid_col();
		int grid_min_feature_num();
		int grid_max_feature_num();
		int pyramid_levels();
		int patch_size();
		int fast_threshold();
		int fscalefactor();
		double ransac_threshold();
		double stereo_threshold();
		int max_iteration();
		double track_precision();
		double mbf();
		/*Stereo Optical Flow Parameters*/

		/*Stereo ORB Optical Flow Parameters*/
		bool usingstereorecify();
		Eigen::Matrix3d K_L();
		Eigen::MatrixXd D_L();
		Eigen::Matrix3d R_L();
		Eigen::MatrixXd P_L();


		Eigen::Matrix3d K_R();
		Eigen::MatrixXd D_R();
		Eigen::Matrix3d R_R();
		Eigen::MatrixXd P_R();

		int nfeatures();
		int max_fast();
		int min_fast();
		int max_cnt();
		/*Stereo ORB Optical Flow Parameters*/


		/*Estimation Parameters*/
		int max_cam_state_size();
		double position_std_threshold();
		double rotation_threshold();
		double translation_threshold();
		double feature_observation_noise();
		/*Estimation Parameters*/

		/*Depth Optimization Parameters*/
		double huber_epsilon();
		double estimation_precision();
		double initial_damping();
		int outler_loop_max_iteration();
		int inner_loop_max_iteration();
		/*Depth Optimization Parameters*/

		/*Estimated Paramters*/
		bool estimate_extrinsic();

		Eigen::Vector3d initial_extrinsic_rotation_cov();

		Eigen::Vector3d initial_extrinsic_translation_cov();

	protected:

	
	};
}
#endif
#endif // USE_VIS