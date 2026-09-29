
#ifndef GFEATURE_H
#define GFEATURE_H

/**
* @file gfeature.h
* @brief feature core algorithm
* @details
		VIS initial depth calculation
* @verbatim
		History
		-1.0	LiShengyu       2020-05-20 creat the file.
  @endverbatim

* @author LiShengyu
* @email
* @version 1.0
* @date 2020-05-20
* @license
*/


#include "gset/gsetbase.h"
#include "gutils/gmutex.h"
#include "gutility.h"
#include "gutils/ExportVIS.h"
#include <Eigen/Eigen>
#ifdef USE_CERES
#include <ceres/ceres.h>
#endif // USE_CERES
using namespace gnut;

namespace gvis
{

	/**
	* @class t_gfeature
	* @brief t_gfeature Class is vis core class for depth optimization.
	*
	* t_gfeature is used for vis core algorithms,including triangulation and depth optimization
	*/
	class VIS_LIBRARY_EXPORT t_gfeature
	{
	public:
		/** @brief constructor. */
		t_gfeature();

		/**
		* @brief constructor.
		*
		* @param[in]  _param    depth optimization parameters
		*/
		t_gfeature(const Feature_param &_param);

		/**
		* @brief constructor.
		*
		* @param[in]  _param     depth optimization parameters
		* @param[in]  new_id     feature_id

		*/
		t_gfeature(const Feature_param &_param,const FeatureIDType& new_id);

		/**
		* @brief judge the feature can be triangulated or not
		*
		* @param[in] cam_states		 all historical camera information
		* @return
			retval = true.   can be triangulated
			retval = false.	 can not be triangulated
		*/
		bool triangulatePoint(const CamStateServer& cam_states) const;

		/**
		* @brief used for mono,judge whether the feature is static or not
		*
		* @param[in] cam_states		 all historical camera information
		* @return
		retval = true.   obvious movement
		retval = false.  nonobvious movement
		*/
		bool checkMotion(const CamStateServer& cam_states) const;

		/**
		* @brief triangulate the initial depth of the feature
		*
		* @param[in] T_c1_c2	transformation matrix between frame_1 and frame_2 [not in the same time]
		* @param[in] z1			observations(X/Z,Y/Z) in frame_1
		* @param[in] z2			observations(X/Z,Y/Z) in frame_2
		* @param[out] p			feature point(X,Y,Z) in frame_1
		*/
		void generateInitialGuess(
			const Eigen::Isometry3d& T_c1_c2, const Eigen::Vector2d& z1,
			const Eigen::Vector2d& z2, Eigen::Vector3d& p) const;

		/**
		* @brief calculate reprojection error
		*
		* @param[in] T_c1_ci	transformation matrix between frame_1 and frame_i
		* @param[in] x			feature in frame_1(X/Z,Y/Z,1/Z)
		* @param[in] z			observations(X/Z,Y/Z) in frame_i
		* @param[out] e			reprojection error
		*/
		void cost(const Eigen::Isometry3d& T_c0_ci,
			const Eigen::Vector3d& x, const Eigen::Vector2d& z,
			double& e) const;

		/**
		* @brief calculate jacobian of depth
		*
		* @param[in] T_c1_ci	transformation matrix between frame_0 and frame_i
		* @param[in] x			feature in frame_0(X/Z,Y/Z,1/Z)
		* @param[in] z			observations(X/Z,Y/Z) in frame_i
		* @param[out] J			jacobian
		* @param[out] r			residual of reprojection
		* @param[out] w			huber weight
		*/
		void jacobian(const Eigen::Isometry3d& T_c0_ci,
			const Eigen::Vector3d& x, const Eigen::Vector2d& z,
			Eigen::Matrix<double, 2, 3>& J, Eigen::Vector2d& r,
			double& w) const;

		/**
		* @brief calculate the optimized position of the feature in the world frame.
		*
		* @param[in] cam_states		all historical camera information
		* @return
		retval = true.   depth was solved successfully
		retval = false.  depth was not solved successfully
		*/
		bool initializePosition(const CamStateServer& cam_states);

		#ifdef USE_CERES
		/**
		* @brief construct ceres problem
		*
		* @param[in] cam_poses		T_cami_cam0
		* @param[in] measurements	all measurements corresponding cam poses
		* @param[out] pos			x,y,z in cam frame not x y 1/z
		* @return
		retval = true.   depth was solved successfully
		retval = false.  depth was not solved successfully
		*/
		bool LMOptimizatePosition(const vector<Eigen::Isometry3d>&cam_poses, const vector<Eigen::Vector2d>&measurements, Eigen::Vector3d& pos);

		/**
		* @brief calculate the optimized position of the feature in the world frame by ceres.
		*		 
		* @param[in] cam_states		all historical camera information
		* @return
		retval = true.   depth was solved successfully
		retval = false.  depth was not solved successfully
		*/
		bool LMOptimizatePositionAndCheck(const CamStateServer& cam_states);
		#endif // USE_CERES
		

	public:


		std::map<StateIDType, Eigen::Vector4d, std::less<StateIDType>,
			Eigen::aligned_allocator<
			std::pair<const StateIDType, Eigen::Vector4d> > > observations;		//store all the observations of the feature 

		FeatureIDType id;				// cur feature id
		Eigen::Vector3d position;		// 3d postion of the feature in  world frame.
		bool is_initialized;			// A indicator to show if the 3d postion of the feature has been initialized or not.
		double initial_depth;			// initial depth from stero orb font end 
		bool is_KeyFrame;				// A indicator to show if feature is generated on the keyframe
		bool is_initialized_NonKey;     // A indicator to show if feature is initialized by nonkeyframe
		bool isLost = false;			// A indicator to show if feature track lost
		

		/**
		* @brief  some parameters in Feature_param,only used in the initialization of t_gfeature
		*/
		double translation_threshold;			// transformation matrix threshold
		double huber_epsilon;					// kernel function
		double estimation_precision;			// threshold of whether the update quantity is iterative
		double	initial_damping;				// initial lambda of Levenberg-Marquart
		int outler_loop_max_iteration;			// iterations
		int inner_loop_max_iteration;			// iterations of depth
		Eigen::Isometry3d T_cam0_cam1;			// transformation matrix between cam0 and cam1
		bool stereo;							// stereo or num
		
	};



}



#endif
