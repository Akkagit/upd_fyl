
#ifndef GIMGPROC_H
#define GIMGPROC_H


/**
* @file gimgproc.h
* @details 
		define some common functions and all parameters that may be used in subclasses
* @endverbatim
* @verbatim
	History
	-1.0	LiShengyu	      2020-05-20 create the file.
  @endverbatim

* @author LiShengyu
* @version 1.0
* @date 2020-05-20
* @license
*/


#include "gvis/gutility.h"
#include "gutils/gstring.h"
#include <Eigen/Eigen>
#include "gutils/ExportVIS.h"
#include "gset/gsetbase.h"


using namespace std;
using namespace gnut;

namespace gvis {

	/**
	* @class t_imgproc
	* @brief t_imgproc Class for reading parameters and virtual functions
	*/
	class VIS_LIBRARY_EXPORT t_gimgproc
	{
	public:

		/**
		* @brief constructor.
		*
		* @param[in]  _set     setbase control
		*/
		t_gimgproc(t_gsetbase* _set);

		/** @brief destructor. */
		~t_gimgproc() {}

		/**
		* @brief save imu observation in _vecimu
		* all imu observations' timestamps are between cur frame and pre frame
		*
		* @param[in] t		 imu observation time
		* @param[in] gyro	 corresponding angular velocity observation
		* @param[in] acc     corresponding acceleration observation
		*/
		void load_imuobs(const double &t, const vector<Eigen::Vector3d> &gyro, const vector<Eigen::Vector3d> &acc);

		/**
		* @brief get cur image observation
		* @param[in] t			 img observation time
		* @param[in] img_path	 corresponding image path
		*/
		void load_imgobs(const double &t,const imgPath& img_path);

		/**
		* @brief image processing main function
		*
		* @param[in] t		 cur image time
		*/
		virtual PointCloud ProcessBatch(const double& t);


		/**
		* @struct imumsg
		* @brief  imudata vector only used in gimgproc class
		*/
		struct imumsg
		{
			double t;
			Eigen::Vector3d linear_acceleration;
			Eigen::Vector3d angular_velocity;
		};

	public:
		//new add;
		imgPath cur_img_path;

		vector<imumsg> _vecimu;			// vector store imu data between pre frame and cur frame
		PointCloud _pointCloud;			// vector store pointcloude

		double ts;						// camera sample interval. //
		int freq;						// camera sample freq. //
		double imu_ts;					// imu sample interval. //
		int imu_freq;					// imu sample freq //

		/* Basic parameters */
		Distortion_Type cam0_distortion_model;				
		Distortion_Type cam1_distortion_model;
		Eigen::Vector2d cam0_resolution;
		Eigen::Vector2d cam1_resolution;
		Eigen::Vector4d cam0_intrinsics;
		Eigen::Vector4d cam1_intrinsics;
		Eigen::Vector4d cam0_distortion_coeffs;
		Eigen::Vector4d cam1_distortion_coeffs;
		Eigen::Matrix3d R_cam0_cam1;
		Eigen::Vector3d t_cam0_cam1;
		Eigen::Isometry3d T_cam0_cam1;
		Eigen::Matrix3d R_cam0_imu;
		Eigen::Vector3d t_cam0_imu;
		Eigen::Isometry3d T_cam0_imu;
		Eigen::Isometry3d T_cam1_imu;
		Eigen::Vector3d t_cam1_imu;
		Eigen::Matrix3d R_cam1_imu;
		/* Basic parameters */
		
		/* Stereo Optical flow */
		int grid_row;				
		int grid_col;
		int grid_min_feature_num;	// min feature number in a grid
		int grid_max_feature_num;	// max feature number in a gird
		int pyramid_levels;			
		int patch_size;				// optical flow tracking window size
		int fast_threshold;			// feature detection threshold
		double ransac_threshold;	// RANSAC threshold
		double stereo_threshold;	// stereo match threshold
		int max_iteration;			// max iterations
		double track_precision;		//
		/* Stereo Optical flow */

		/* Stereo Orb Optical Flow */
		int nfeatures;
		int max_fast;
		int min_fast;
		int max_cnt;				// use in st storb mono
		int fscalefactor;
		double mbf;
		bool usingstereorecify;
		Eigen::Matrix3d K_L;
		Eigen::MatrixXd D_L;
		Eigen::MatrixX3d R_L;
		Eigen::MatrixXd P_L;
		Eigen::Matrix3d K_R;
		Eigen::MatrixXd D_R;
		Eigen::MatrixX3d R_R;
		Eigen::MatrixXd P_R;
		/* Stereo Orb Optical Flow */
		
		/* Estimator parameter */
		int max_cam_state_size;
		double position_std_threshold;
		double rotation_threshold;
		double translation_threshold;
		double feature_observation_noise;
		/* Estimator parameter */

		/* Depth optimation parameter */
		double huber_epsilon;
		double estimation_precision;
		double initial_damping;
		int outler_loop_max_iteration;
		int inner_loop_max_iteration;
		/* Depth optimation parameter */

		bool estimate_extrinsic;
		Eigen::Vector3d initial_extrinsic_rotation_cov;
		Eigen::Vector3d initial_extrinsic_translation_cov;
	};
}





#endif
