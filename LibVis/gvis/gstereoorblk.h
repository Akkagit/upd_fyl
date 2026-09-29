
#ifndef GSTEREOORBLK_H
#define GSTEREOORBLK_H

/**
* @file gstereoorblk.h
* @brief optical flow tracking from ORB-SLAM2
* @details
		VIS the font end are carried out through optical flow tracking combined 
		with ORB descriptor matching£¬this method run the slowest
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


#include "string"
#include <vector>
#include <map>
#include <memory>
#include <thread>
#include <fstream>
#include <opencv2/opencv.hpp>
#include "opencv2/core/core.hpp"  
#include "opencv2/core/utility.hpp"  
#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "gutils/ExportVIS.h"
#include "gset/gsetbase.h"
#include "gutility.h"
#include "gimgproc.h"
#include "gfeature.h"
#include "gframe.h"
#include "ORBextractor.h"
#include "ORBmatcher.h"

using namespace std;
using namespace ORB_SLAM2;

namespace gvis
{

	/**
	* @class t_gstereoorblk
	* @brief t_gstereoorblk Class for completing the tracking of feature points 
	* with ORB descriptor in
	*
	*/

	class VIS_LIBRARY_EXPORT t_gstereoorblk:public t_gimgproc
	{
	public:
		/*
		*@brief constructor.
		*
		* @param[in]  s        setbase control
		*/
		t_gstereoorblk(gnut::t_gsetbase* _set);

		/** @brief destructor. */
		~t_gstereoorblk();

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
		* @brief calibrate the image according to the parameters
		*
		* @param[in] img				input image
		* @param[out] UndistortImg		ouput image after undistorting
		* @param[in] K0	
		* @param[in] D0
		* @param[in] R0
		* @param[in] P0
		*/
		void RectWholeImg(const cv::Mat &img, cv::Mat &UndistortImg, const cv::Mat &K0,
			const cv::Mat &D0, const cv::Mat &R0, const cv::Mat &P0);

		/**
		* @brief encapsulate the function RectWholeImg()
		*
		*/
		void RectifyImg();

		/**
		* @brief extract the feature points on the first image
		*
		*/
		void initializeFirstFrame();

		/**
		* @brief optical flow tracking is performed on the non-first image
		*
		*/
		void FeatureTracker();

		/**
		* @brief add new feature if cur feature number is less than threshold
		*
		*/
		void addNewFeatures();


		/**
		* @brief publish the front-end tracking results in _pointCloud
		*
		*/
		void publish();

		/**
		* @brief feature tracking effect of visual display
		*
		*/
		void drawImage();
		//void drawoneImage();


		/**
		* @brief feature extraction and corresponding descriptor calculation 
		* were performed on the left and right images
		*
		*/
		void featureDetect();

		/**
		* @brief The left and right feature points are matched,
		* and the wrong point pairs are deleted
		*
		* @param[in] KeyPtLeft				feature points in cam0
		* @param[in] KeyPtRight				feature points in cam1
		* @param[in] PtDescriptorsLeft		cam0 feature descriptor
		* @param[in] PtDescriptorsRight		cam1 feature descriptor
		* @param[out] uPtRight				not used
		* @param[out] ptDepth				stereo depth 
		*/
		void StereoMatch(vector<cv::KeyPoint> &KeyPtLeft,
			vector<cv::KeyPoint> &KeyPtRight, cv::Mat PtDescriptorsLeft,
			cv::Mat PtDescriptorsRight, vector<float> &uPtRight,
			vector<float> &ptDepth);

		/**
		* @brief the feature in cam0 and cam1 are constrained by polar line
		*
		* @param[in] cam0_track_points		feature points in cam0
		* @param[in] cam1_track_points		feature points in cam1
		* @param[in] status					an identifier to indicate whether it is a correct match
		*/
		void Epipolar_Constraint(vector<cv::Point2f> &cam0_track_points,
			vector<cv::Point2f> &cam1_track_points, vector<uchar> &status);

		/**
		* @brief create cur image pyramid
		*
		*/
		void createImagePyramids();

		/**
		* @brief circle the image to ensure sparse feature points
		*
		* @param[in] point		feature points in cam0
		*/
		
		cv::Mat setmask(std::vector<cv::Point2f> &point);

		/**
		* @brief delete the error according to status
		*
		* @param[in] status		an identifier to indicate whether it is a correct match
		*/
		template <typename vecType>
		void reduceVector(vector<vecType> &v, vector<uchar>  &status);

		/**
		* @brief determine if the feature point is within the boundary
		*
		* @param[in] pt		pixel coordinate 
		*/
		bool inBorder(const cv::Point2f &pt);

		/**
		* @brief normalized feature point
		*
		* @param[in] intrinsics		camera intrinsics
		* @param[in] p				pixel coordinate
		* @param[out] P				normalized coordinate
		*
		*/
		void Normalised(cv::Vec4d intrinsics, const Eigen::Vector2d& p, Eigen::Vector3d& P);

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
			vector<uchar>& inlier_markers);

		/**
		* @brief normalized feature point,taken scale factor into consideration
		*
		*/
		void rescalePoints(
			vector<cv::Point2f>& pts1, vector<cv::Point2f>& pts2,
			float& scaling_factor);

		/**
		* @brief create the orb extract thread
		*
		* @param[in] im				img
		* @param[in] flag			flag=0,create for img0
		*							flag=1,create for img1
		*/
		void ExtractORB(const cv::Mat &im, int flag);

		/**
		* @brief sort by tracking number
		*
		*/
		static bool featureCompareByTrack(const t_gframe::PointsAttri & f1, const t_gframe::PointsAttri & f2);


	public:
		/* some parameter in parent class ,but form from eigen to cv */
		 std::string  _cam0_distortion_model;
		 cv::Vec2i  _cam0_resolution;
		 cv::Vec4d  _cam0_intrinsics;
		 cv::Vec4d  _cam0_distortion_coeffs;
		 
		 std::string  _cam1_distortion_model;
		 cv::Vec2i  _cam1_resolution;
		 cv::Vec4d  _cam1_intrinsics;
		 cv::Vec4d  _cam1_distortion_coeffs;

		cv::Matx33d  _R_cam0_imu;
		cv::Vec3d  _t_cam0_imu;
		cv::Matx33d  _R_cam1_imu;
		cv::Vec3d  _t_cam1_imu;
		/* some parameter in parent class ,but form from eigen to cv */
		
		/* stereo-rectify parameters */
		cv::Mat  _K_L;
		cv::Mat  _D_L;
		cv::Mat  _R_L;
		cv::Mat  _P_L;
		cv::Mat  _K_R;
		cv::Mat  _D_R;
		cv::Mat  _R_R;
		cv::Mat  _P_R;
		float  mb;		// Baseline length
		/* stereo-rectify parameters */

		ORBextractor* mpORBextractorLeft, *mpORBextractorRight;  //ORB feature extraction for img0 and img1
		bool isFirstImg;		// A indicator to show if it is the first frame
		bool IsInitialized;		// A indicator to show if the system is initialized

		t_gframe *prev_frame;	// store information about cur img
		t_gframe *cur_frame;	// store information about pre img,used for tracking in next frame

		vector<cv::KeyPoint>  mvKeysLeft, mvKeysRight;	// store the feature points extracted from the img0 and img1
		cv::Mat mDescriptorsLeft, mDescriptorsRight;	// descriptors, corresponding to feature, are stored by column

		std::vector<float> mvuRight;					// the x-coordinate of the right feature point of the final interpolation
		std::vector<float> mvDepth;						// stereo depth
		long long int npointID;							// cur feature id

	};
}



#endif // !gstereoorblk_H
