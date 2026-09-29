
#ifndef GMONOLK_H
#define GMONOLK_H


/**
* @file gmonolk.h
* @brief optical flow method comes from vins
* @detials
	 mono optical flow can be performed

* @verbatim
	History
	-1.0	LiShengyu	      2020-05-20 create the file.
  @endverbatim

* @author LiShengyu
* @version 1.0
* @date 2020-05-20
* @license
*/


#include "gutils/ExportVIS.h"
#include "gset/gsetbase.h"
#include "gutility.h"
#include "gimgproc.h"

using namespace std;
using namespace gnut;

namespace gvis {

	/**
	* @class t_gmonolk
	* @brief t_gmonolk Class for completing the tracking of feature points
	*
	*/
	class VIS_LIBRARY_EXPORT t_gmonolk : public t_gimgproc
	{
	public:

		/**
		* @brief constructor.
		*
		* @param[in]  s        setbase control
		* @param[in]  version  version of the gcoder
		* @param[in]  sz       size of the buffer
		*/
		t_gmonolk(gnut::t_gsetbase* _set);

		/** @brief destructor. */
		~t_gmonolk() {}

		/**
		* @brief image processing main function
		*
		* @param[in] t		 cur image time
		*/
		virtual PointCloud ProcessBatch(const double& t);

		/**
		* @brief Initialization of feature extraction
		* parameters converted from eigen form to CV form
		*/
		bool Initialize();

		/**
		* @brief Initialization of feature extraction
		*
		* @param[in] pt		 cur image time
		*/
		bool inBorder(const cv::Point2f &pt);

		/**
		* @brief Remove the unmarked elements within a vector.
		*
		*/
		void _reduceVector(vector<cv::Point2f> &v, vector<uchar> status);

		/**
		* @brief Remove the unmarked elements within a vector.
		*
		*/
		void _reduceVector(vector<int> &v, vector<uchar> status);

		/**
		* @brief set mask in img to remove feature that stay too close
		*
		*
		*/
		void setMask();

		/**
		* @brief ransac to remove outliers
		*
		*
		*/
		void rejectWithF();

		/**
		* @brief calculate the distance between pt1 and pt2
		*
		*
		*/
		double distance(cv::Point2f &pt1, cv::Point2f &pt2);

		/**
		* @brief imu observations predict camera rotation
		*
		* @param[out] cam0_R_p_c		cam0 rotation from pre to cur
		* @param[out] cam1_R_p_c		cam1 rotation from pre to cur
		*
		*/
		void integrateImuData(cv::Matx33f& cam0_R_p_c);

		/**
		* @brief use imu observations to predict feature in cur frame
		*
		* @param[in] input_pts			feature points in cam0
		* @param[in] R_p_c				cam0 rotation from pre to cur
		* @param[in] intrinsics			cam0 intrinsics
		* @param[in] compensated_pts	feature points after predicting
		*
		*/
		void predictFeatureTracking(
			const vector<cv::Point2f>& input_pts,
			const cv::Matx33f& R_p_c,
			const cv::Vec4d& intrinsics,
			vector<cv::Point2f>& compensated_pts);

		/**
		* @brief optical flow tracking is performed on the non-first image
		*
		*/
		void trackFeatures();

		/**
		* @brief publish cur feature tracking result
		*
		*/
		void publish();

		/**
		* @brief add new feature if cur feature number is less than threshold
		*
		*/
		void addnewFeatures();

		/**
		* @brief feature tracking effect of visual display
		*
		*/
		void drawFeatures();

		/**
		* @brief dedistorts the feature points
		*
		*/
		void rescalePoints(
			vector<cv::Point2f>& pts1, vector<cv::Point2f>& pts2,
			float& scaling_factor);

		/**
		* @brief  two point ransac to remove outliers
		*
		* @param[in] pts1					feature points in cam0
		* @param[in] pts2					feature points in cam1
		* @param[in] R_p_c					cam0 rotation from pre to cur
		* @param[in] distortion_model		cam0 distortion type ,ratan or
		* @param[in] distortion_coeff		cam0 distortion coeff
		* @param[in] inlier_error			threshold
		* @param[in] success_probability	inlier success probability
		* @param[out] inlier_markers		inlier -1, outlier 0
		*
		*/
		void twoPointRansac(
			const vector<cv::Point2f>& pts1, const vector<cv::Point2f>& pts2,
			const cv::Matx33f& R_p_c, const cv::Vec4d& intrinsics,
			const std::string& distortion_model,
			const cv::Vec4d& distortion_coeffs,
			const double& inlier_error,
			const double& success_probability,
			vector<int>& inlier_markers);

		/**
		* @brief dedistorts the feature points
		*
		*/
		void undistortPoints(
			const vector<cv::Point2f>& pts_in,
			const cv::Vec4d& intrinsics,
			const string& distortion_model,
			const cv::Vec4d& distortion_coeffs,
			vector<cv::Point2f>& pts_out,
			const cv::Matx33d &rectification_matrix = cv::Matx33d::eye(),
			const cv::Vec4d &new_intrinsics = cv::Vec4d(1, 1, 0, 0));

	public:
		
		int row, col;		 
		cv::Mat mask;							// to remove feature points that stay close
		cv::Mat prev_img, cur_img;				// mat form used for KLT
		vector<cv::Point2f> n_pts;				// used in add new features
		vector<cv::Point2f> prev_pts, cur_pts;	// record pts in cur and prev feature point
		vector<int> prev_ids,cur_ids;			// record id in cur and prev img
		vector<int> track_cnt;					// track number

		double cur_time;
		double prev_time;
		int n_id;

		bool IsInitialized = false;				// used for class parameter initial
		bool flow_back = true;					// KLT back to check
	
		int before_tracking;		// cur feature number before tracking
		int after_tracking;			// cur feature number after tracking	
		int after_setmask;			// cur feature number after set mask
		int after_ransac;			// cur feature number after ransac
		int after_twopoint_ransac;	// cur feature number after two point ransac

		string _cam0_distortion_model;
		cv::Vec4d _cam0_intrinsics;
		cv::Vec4d _cam0_distortion_coeffs;

		cv::Matx33d _R_cam0_imu;
		cv::Vec3d _t_cam0_imu;
	};
}





#endif