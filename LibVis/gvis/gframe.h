
#ifndef GFRAME_H
#define GFRAME_H

/**
* @file gframe.h
* @brief frame core algorithm
* @details
		VIS define the t_gframe class only used in stereo orb font end 
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
#include <opencv2/opencv.hpp>
#include "ORBmatcher.h"
#include "gutility.h"
#include "gset/gsetbase.h"
#include "gutils/gmutex.h"
#include "gutils/ExportVIS.h"
#include <Eigen/Eigen>
using namespace gnut;


using namespace std;
using namespace ORB_SLAM2;

namespace gvis
{

	/**
	* @class t_gframe
	* @brief t_gfeature Class is core class for stereo orb font end.
	*
	* t_gfeature is used for vis core algorithms,including pyramid,keypoint,timestamp and descriptor 
	*/
	class VIS_LIBRARY_EXPORT t_gframe
	{

	public:
		/**
		* @struct PointsAttri
		* @brief img,pyramid,keypoint,timestamp and descriptor,only used in stereo orb font end
		*/
		struct PointsAttri
		{
			cv::Point2f  cam0_pt;
			cv::Point2f  cam1_pt;
			
			cv::KeyPoint cam0_keypt;
			cv::KeyPoint cam1_keypt;

			cv::Mat cam0_descriptor;
			cv::Mat cam1_descriptor;

			float  depth = 0.0;
			int id = 0;
			int track_cnt = 0;
			double response = 0.0;
		};
		double time;							// timestamp of the img
		cv::Mat cam0_img;						// A copy of the original observation data
		cv::Mat cam1_img;
		dataIMG img_msg;						// original observation data
		std::vector<cv::Mat> cam0_pyramid;		// pyramid of cam0
		std::vector<cv::Mat> cam1_pyramid;		// pyramid of cam1
		vector<PointsAttri> vpointInFrame;		// the feature in cur frame
	};
}
#endif // !GFRAME_H
