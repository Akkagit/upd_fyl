
#ifndef GVIS_H
#define GVIS_H

/**
* @file gvis.h
* @brief vis core algorithm
* @details
		Starting from the receiving pointcloud data, including adding 
		observation values, searching for lost feature points and the 
		camera state window is full for measurement update 
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

#include "gutils/gmutex.h"
#include "gutils/ExportVIS.h"
#include "gutility.h"
#include "gfeature.h"
#include "gimgproc.h"
#include "gmonolk.h"
#include "gstereolk.h"
#include "gstereoorblk.h"
#include "gpublish.h"

using namespace gnut;
using namespace std;

namespace gvis
{
	/**
	* @typedef MapServer
	* @brief collection of all pointclouds
	*/
	typedef map<FeatureIDType, t_gfeature> MapServer;

	/**
	* @typedef EstimatedPoint
	* @brief 
	*/
	typedef vector<pair<FeatureIDType, t_gfeature>> EstimatedPoint;

	/**
	* @class t_gvis
	* @brief t_gvis Class is vis core class in dataprocessing
	*
	* t_gvis is used for vis add feature observations,camera meaurement update
	* 
	*/

	class VIS_LIBRARY_EXPORT t_gvis
	{
	public:
		/**
		* @brief constructor.
		*
		* @param[in]  set        setbase control
		*/
		t_gvis(t_gsetbase* set);
		

		/**
		* @brief destructor
		*
		*/
		~t_gvis();

	

		/**
		* @brief read Chi square distribution data
		*/
		void readChisquare_test();

		/**
		* @brief determine key frame or not
		* @param[in]  msg     cur pointcloud data
		* @return bool
		*/
		bool keyframeCheck(const PointCloud& msg);

		/**
		* @brief add cur pointcloud into mapserver
		* @param[in]  msg     cur pointcloud data
		*/
		void addFeatureObservations(const PointCloud &msg);

		/**
		* @brief process camera measurement update
		*/
		//void gns_time_update();




		/**
		* @brief get one feature's Jacobian matrix in one observation 
		* and calculate residual
		*
		* @param[in] cam_state_id		camera state number flag
		* @param[in] feature_id			feature id
		* @param[out] H_x				Jacobian about camera pose
		* @param[out] H_f				Jacobian about feature point in the world frame
		* @param[out] r					reprojection error
		*/
		void measurementJacobian(
			const StateIDType& cam_state_id,
			const FeatureIDType& feature_id,
			Eigen::Matrix<double, 4, 6>& H_x, Eigen::Matrix<double, 4, 3>& H_f, Eigen::Vector4d& r);

		/**
		* @brief get one feature's Jacobian matrix in one observation for mono
		* and calculate residual
		*
		* @param[in] cam_state_id		camera state number flag
		* @param[in] feature_id			feature id
		* @param[out] H_x				Jacobian about camera pose
		* @param[out] H_f				Jacobian about feature point in the world frame
		* @param[out] r					reprojection error
		*/
		void measurementJacobian(
			const StateIDType& cam_state_id,
			const FeatureIDType& feature_id,
			Eigen::Matrix<double, 2, 6>& H_x, Eigen::Matrix<double, 2, 3>& H_f, Eigen::Vector2d& r);

		/**
		* @brief get one feature's Jacobian matrix in all observations
		* and calculate residual
		*
		* @param[in] feature_id			feature id
		* @param[in] cam_state_ids		camera states number flag
		* @param[out] H_x				Jacobian about all camera poses		
		* @param[out] r					reprojection errors
		*/
		bool featureJacobian(const FeatureIDType& feature_id,
			const std::vector<StateIDType>& cam_state_ids,
			Eigen::MatrixXd& H_x, Eigen::VectorXd& r);

		/**
		* @brief find the two camera state id to be removed
		*
		* @param[out] rm_cam_state_ids		camera state ids to be removed
		*/
		void findRedundantCamStates(vector<StateIDType>& rm_cam_state_ids);

		/**
		* @brief find the two camera state id to be removed in nonkeyframes
		*
		* @param[out] rm_cam_state_ids		camera state ids to be removed
		*/
		void findNonKeyCamStates(vector<StateIDType>& rm_cam_state_ids);


		// these functions just for test
		void test_cam_time_update();
		void test_featureJacobian(const FeatureIDType& feature_id,
			const std::vector<StateIDType>& cam_state_ids,
			Eigen::MatrixXd& H_x, Eigen::VectorXd& r);

	public:
		
	    map<int, double> chi_squared_test_table;    // store chi square test table

		CamStateServer cam_states;					// store all camera state poses
		CAMState extrinsic;                         // store extrinsic for estimate

		IMUState imu_state;							// used for trajectory display
		MapServer map_server;						// store all feature observations
		StateIDType imu_state_id;
		StateIDType next_id;						// used for recursive imu_state_id
		//EstimatedPoint es_point;					// save feature which are estimated		
		//FeatureIDType es_id;						// save feature ids which are estimated	
		bool isKeyFrame = false;					// A indicator to show if the frame is keyframe or not

		Eigen::MatrixXd Phi;						// only used for test,not used in main function
		Eigen::MatrixXd P_IC;
		Eigen::MatrixXd P_CC;
		Eigen::MatrixXd P_IC_sav;
		Eigen::MatrixXd P_CC_sav;
		Eigen::VectorXd Dx_sav;
		int Dx_row;

		int nonFrame_num = 0;                       // num of nonkeyframe not used 
		vector<StateIDType> nonKeyFrameBuffer;		// store nonkeyframes	
		vector<StateIDType> KeyFrameBuffer;			// store keyframes
		
		double tracking_rate;						// feature tracking rate between cur frame and pre frame
		Feature_param f_param;						// parameters used in depth optimization
	
		Eigen::Matrix3d R_cam0_imu;                 // rotation from cam0 to imu
		Eigen::Vector3d t_cam0_imu;					// translation from cam0 to imu
		Eigen::Isometry3d T_cam0_imu;				// transformation from cam0 to imu

		Eigen::Matrix3d R_cam0_cam1;				// rotation from cam0 to cam1
		Eigen::Vector3d t_cam0_cam1;				// translation from cam0 to cam1
		Eigen::Isometry3d T_cam0_cam1;				// transformation from cam0 to cam1

		Eigen::Vector4d cam0_intrinsics;			// cam0 intrinsics used in keyframecheck()
		int max_camstate_size;						// max num of camera states	
		double observation_noise;					// observation noise of feature
		int num_of_cam;								// camera number 1 or 2
		bool stereo;								// false-number 1  true-number 2
		bool usingstereorecity;						// using orb or not

	protected:
		//new add for data get
		t_gimgproc *imgproc;
		t_gpublish publisher;

	};



}



#endif