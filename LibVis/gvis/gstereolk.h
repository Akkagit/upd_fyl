
#ifndef GSTEREOLK_H
#define GSTEREOLK_H


/**
* @file gstereolk.h
* @brief optical flow tracking from MSCKF
* @details
		VIS the feature points are divided into the grid and optical flow tracking 
		is carried out,this method runs fastest
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



#include "gutils/ExportVIS.h"
#include "gset/gsetbase.h"
#include "gutility.h"
#include "gimgproc.h"

using namespace std;
using namespace gnut;



namespace gvis {

	/**
	* @class t_gstereolk
	* @brief t_gstereolk Class for completing the tracking of feature points
	*
	*/
	class VIS_LIBRARY_EXPORT t_gstereolk: public t_gimgproc
	{
	public:
		
		/**
		* @brief constructor.
		*
		* @param[in]  s        setbase control
		*/
		t_gstereolk(gnut::t_gsetbase* _set);

		/** @brief destructor. */
		~t_gstereolk() {}

		/**
		* @brief Initialization of feature extraction
		* parameters converted from eigen form to CV form
		*/
		bool Initialize();

		/**
		* @brief image processing main function
		*
		* @param[in] t		 cur image time
		*/
		virtual PointCloud ProcessBatch(const double& t);

		/**
		* @brief create cur image pyramid
		*
		*/
		void createImagePyramids();
		

		/**
		* @brief extract the feature points on the first image
		*
		*/
		void initializeFirstFrame();

		/**
		* @brief The left and right feature points are matched, 
		* and the wrong point pairs are deleted
		*
		* @param[in] cam0_points		feature points in cam0
		* @param[in] cam1_points		feature points in cam1
		* @param[in] inlier_markers     an identifier to indicate whether it is a correct match
		*/
		void stereoMatch(
			const vector<cv::Point2f>& cam0_points,
			vector<cv::Point2f>& cam1_points,
			vector<unsigned char>& inlier_markers);

		/**
		* @brief feature tracking effect of visual display
		*
		*/
		void drawFeaturesStereo();

		/**
		* @brief publish cur feature tracking result
		*
		*/
		void publish();
		
		/**
		* @brief optical flow tracking is performed on the non-first image
		*
		*/
		void trackFeatures();

		/**
		* @brief imu observations predict camera rotation
		*
		* @param[out] cam0_R_p_c		cam0 rotation from pre to cur
		* @param[out] cam1_R_p_c		cam1 rotation from pre to cur
		*
		*/
		void integrateImuData(cv::Matx33f& cam0_R_p_c, cv::Matx33f& cam1_R_p_c);

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
		* @brief add new feature if cur feature number is less than threshold
		*
		*/
		void addNewFeatures();

		/**
		* @brief put feauture in grids,delete if one grid'feature is more than threshold
		*
		*/
		void pruneGridFeatures();

		/**
		* @brief dedistorts the feature points
		*
		*/
		vector<cv::Point2f> distortPoints(
			const vector<cv::Point2f>& pts_in,
			const cv::Vec4d& intrinsics,
			const string& distortion_model,
			const cv::Vec4d& distortion_coeffs);

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
		
		/**
		* @brief sort by response
		*
		*/
		static bool featureCompareByResponse(
			const FeaturePoint& f1,
			const FeaturePoint& f2);

		/**
		* @brief Remove the unmarked elements within a vector.
		*
		*/
		template <typename T>
		void removeUnmarkedElements(
			const std::vector<T>& raw_vec,
			const std::vector<unsigned char>& markers,
			std::vector<T>& refined_vec);

		/**
		* @brief normalized feature point,taken scale factor into consideration
		*
		*/
		void rescalePoints(
			vector<cv::Point2f>& pts1, vector<cv::Point2f>& pts2,
			float& scaling_factor);

		/**
		* @brief sort by response for keypoint
		*
		*/
		static bool keyPointCompareByResponse(
			const cv::KeyPoint& pt1,
			const cv::KeyPoint& pt2);

		/**
		* @brief sort by lifetime
		*
		*/
		static bool featureCompareByLifetime(
			const FeaturePoint& f1,
			const FeaturePoint& f2);

		

	public:

		/**
		* @typedef CamStateServer
		* @brief store all feature observations
		*/
		typedef std::map<int, std::vector<FeaturePoint> > GridFeatures;

		
		cv::Ptr<cv::Feature2D> detector_ptr;				// Fast feature point detector

		dataIMG curr_img_msg;								// store information about cur img
		dataIMG prev_img_msg;								// store information about pre img,used for tracking in next frame
		cv::Mat cam0_curr_img;								// cur img0 clone
		cv::Mat cam1_curr_img;								// cur img1 clone
		std::vector<cv::Mat> curr_cam0_pyramid;				// cur img0 pyramid
		std::vector<cv::Mat> curr_cam1_pyramid;				// cur img1 pyramid
		std::vector<cv::Mat> prev_cam0_pyramid;				// pre img0 pyramid used for optical flow
		long long int next_feature_id;						// ID for the next new feautre
		std::shared_ptr<GridFeatures> prev_features_ptr;    // store all feature observations in pre frame
		std::shared_ptr<GridFeatures> curr_features_ptr;	// store all feature observations in cur frame


		int before_tracking;		// cur feature number before tracking
		int after_tracking;			// cur feature number after tracking	
		int after_matching;			// cur feature number after stereomatch
		int after_ransac;			// cur feature number after two point ransac
		map<long long int , int> track_cnt;		// not used in stereo optical flow
		bool isFirstImg;		    // A indicator to show if it is the first frame
		bool IsInitialized;			// A indicator to show if the system is initialized

		/* some parameter in parent class ,but form from eigen to cv */
		cv::Matx33d _R_cam0_cam1;
		cv::Vec3d _t_cam0_cam1;
		cv::Matx33d _R_cam0_imu;
		cv::Vec3d _t_cam0_imu;
		cv::Matx33d _R_cam1_imu;

		string _cam0_distortion_model;
		cv::Vec4d _cam0_intrinsics;
		cv::Vec4d _cam0_distortion_coeffs;

		string _cam1_distortion_model;
		cv::Vec4d _cam1_intrinsics;
		cv::Vec4d _cam1_distortion_coeffs;
		/* some parameter in parent class ,but form from eigen to cv */
	
	};
}






#endif