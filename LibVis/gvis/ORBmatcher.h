
#ifndef ORBMATCHER_H
#define ORBMATCHER_H

/**
* @file ORBMATCHER.h
* @brief ORB feature matcher core algorithm
* @verbatim
		History
		-1.0	LiShengyu	       2020-05-20 creat the file.
  @endverbatim

* @author LiShengyu
* @email
* @version 1.0
* @date 2020-05-20
* @license
*/

#include<vector>
#include<opencv2/core/core.hpp>
#include<opencv2/features2d/features2d.hpp>

namespace ORB_SLAM2
{

	class ORBmatcher
	{
	public:

		ORBmatcher(float nnratio = 0.6, bool checkOri = true);

		// Computes the Hamming distance between two ORB descriptors
		static int DescriptorDistance(const cv::Mat &a, const cv::Mat &b);

		// Search matches between Frame keypoints and projected MapPoints. Returns number of matches
		// Used to track the local map (Tracking)
		// int SearchByProjection(Frame &F, const std::vector<MapPoint*> &vpMapPoints, const float th = 3);

	public:

		static const int TH_LOW;
		static const int TH_HIGH;
		static const int HISTO_LENGTH;

	protected:
		float mfNNratio;
		bool mbCheckOrientation;
	};
}


#endif // !1