
#include "gstereoorblk.h"
#include "gset/gsetvis.h"
#include "opencv2/core/eigen.hpp"

using namespace cv;
using namespace Eigen;
using namespace gvis;


gvis::t_gstereoorblk::t_gstereoorblk(gnut::t_gsetbase* _set) :t_gimgproc(_set),
isFirstImg(true), npointID(0), IsInitialized(false)
{
	Initialize();//这个放在构造函数
	return;
}

gvis::t_gstereoorblk::~t_gstereoorblk()
{
	destroyAllWindows();
	return;
}


bool gvis::t_gstereoorblk::Initialize()
{
	mpORBextractorLeft = new ORBextractor(nfeatures, fscalefactor, pyramid_levels, max_fast, min_fast);
	mpORBextractorRight = new ORBextractor(nfeatures, fscalefactor, pyramid_levels, max_fast, min_fast);

	
	eigen2cv(R_cam1_imu, _R_cam1_imu);
	eigen2cv(t_cam1_imu, _t_cam1_imu);
	eigen2cv(R_cam0_imu, _R_cam0_imu);
	eigen2cv(t_cam0_imu, _t_cam0_imu);
	
	_cam0_distortion_model = distortion2str(cam0_distortion_model);
	_cam0_intrinsics << cam0_intrinsics(0), cam0_intrinsics(1), cam0_intrinsics(2), cam0_intrinsics(3);
	_cam0_distortion_coeffs << cam0_distortion_coeffs(0), cam0_distortion_coeffs(1), cam0_distortion_coeffs(2), cam0_distortion_coeffs(3);

	_cam1_distortion_model = distortion2str(cam1_distortion_model);
	_cam1_intrinsics << cam1_intrinsics(0), cam1_intrinsics(1), cam1_intrinsics(2), cam1_intrinsics(3);
	_cam1_distortion_coeffs << cam1_distortion_coeffs(0), cam1_distortion_coeffs(1), cam1_distortion_coeffs(2), cam1_distortion_coeffs(3);

	eigen2cv(K_L, _K_L);
	eigen2cv(D_L, _D_L);
	eigen2cv(R_L, _R_L);
	eigen2cv(P_L, _P_L);
	eigen2cv(K_R, _K_R);
	eigen2cv(D_R, _D_R);
	eigen2cv(R_R, _R_R);
	eigen2cv(P_R, _P_R);

	
	_cam0_intrinsics= cv::Vec4d(_P_L.at<double>(0, 0),_P_L.at<double>(1, 1), _P_L.at<double>(0, 2), _P_L.at<double>(1, 2));

	_cam1_intrinsics = cv::Vec4d(_P_R.at<double>(0, 0),_P_R.at<double>(1, 1),_P_R.at<double>(0, 2), _P_R.at<double>(1, 2));


	Isometry3d T_l_re = Isometry3d::Identity();
	Isometry3d T_r_re = Isometry3d::Identity();
	Matrix3d R_l_re, R_r_re;
	cv::Mat cvT_cam0_imu, cvT_cam1_imu;
	cv::cv2eigen(_R_L, R_l_re);
	cv::cv2eigen(_R_R, R_r_re);
	T_l_re.linear() = R_l_re.transpose();
	T_r_re.linear() = R_r_re.transpose();
	
	T_cam0_imu = T_cam0_imu * T_l_re;
	T_cam1_imu = T_cam1_imu * T_r_re;
	T_cam0_cam1 = T_cam1_imu.inverse() * T_cam0_imu;
	
	R_cam0_imu = T_cam0_imu.linear();
	t_cam0_imu = T_cam0_imu.translation();
	R_cam1_imu = T_cam1_imu.linear();
	t_cam1_imu = T_cam1_imu.translation();
	R_cam0_cam1 = T_cam0_cam1.linear();
	t_cam0_cam1 = T_cam0_cam1.translation();
	
	cv::eigen2cv(R_cam0_imu, _R_cam0_imu);
	cv::eigen2cv(t_cam0_imu, _t_cam0_imu);
	cv::eigen2cv(R_cam1_imu, _R_cam1_imu);
	cv::eigen2cv(t_cam1_imu, _t_cam1_imu);

	mb = mbf / _cam0_intrinsics(0);
	return true;
}

PointCloud gvis::t_gstereoorblk::ProcessBatch(const double& t)
{
	cout << "==================================" << endl;
	cout << "start process image" << endl;

	//for (int i = 0; i < _vecimg.size(); i++)
	//{
	//	imgPath &img_msg = _vecimg.at(i);
	//	if (double_eq(img_msg.t, t))
	//	{
	//		//在这里生成新的share_ptr

	//		cv::Mat _img0 = cv::imread(img_msg.img0_path, 0);
	//		cv::Mat _img1 = cv::imread(img_msg.img1_path, 0);

	//		dataIMG _img_msg;
	//		_img_msg.t = img_msg.t;
	//		_img_msg.img0 = make_shared<cv::Mat>(_img0);
	//		_img_msg.img1 = make_shared<cv::Mat>(_img1);
	//		cur_frame = new t_gframe;
	//		cur_frame->time = img_msg.t;
	//		cur_frame->cam0_img = _img_msg.img0->clone();
	//		cur_frame->cam1_img = _img_msg.img1->clone();
	//		cur_frame->img_msg = _img_msg;
	//		
	//		break;
	//	}
	//}

	cv::Mat _img0 = cv::imread(cur_img_path.img0_path, 0);
	cv::Mat _img1 = cv::imread(cur_img_path.img1_path, 0);

	dataIMG _img_msg;
	_img_msg.t = cur_img_path.t;
	_img_msg.img0 = make_shared<cv::Mat>(_img0);
	_img_msg.img1 = make_shared<cv::Mat>(_img1);
	cur_frame = new t_gframe;
	cur_frame->time = cur_img_path.t;
	cur_frame->cam0_img = _img_msg.img0->clone();
	cur_frame->cam1_img = _img_msg.img1->clone();
	cur_frame->img_msg = _img_msg;

	// 图像均衡化
	//if (EQUALIZE)
	//{
	//	cv::Mat img_l, img_r;
	//	cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE(3.0, cv::Size(8, 8));

	//	clahe->apply(cur_frame->cam0_img, img_l);
	//	clahe->apply(cur_frame->cam1_img, img_r);
	//	
	//	cur_frame->cam0_img = img_l.clone();
	//	cur_frame->cam1_img = img_r.clone();
	//}
	if (!IsInitialized)
	{
		//Initialize();
		IsInitialized = true;
		
	}


    // 图像校正
	if (usingstereorecify)
		RectifyImg();

	// 构建图像金字塔
	createImagePyramids();
	
	if (isFirstImg)
	{
		initializeFirstFrame();
		isFirstImg = false; 
		//只保留最后一个imu
		imumsg last_imu = _vecimu.at(_vecimu.size() - 1);
		_vecimu.clear();
		_vecimu.push_back(last_imu);
		prev_frame = new t_gframe;
	}
	else
	{
		FeatureTracker();

		addNewFeatures();

	}
		
	publish();

	// 图像显示
	drawImage();
	//drawoneImage();
	// 数据交换
	delete prev_frame;
	prev_frame = cur_frame;
	cur_frame = NULL;
	return _pointCloud;
}



void gvis::t_gstereoorblk::initializeFirstFrame()
{

	featureDetect();
	

	std::vector<float> vuRightTemp;
	std::vector<float> vDepthTemp;

	StereoMatch(mvKeysLeft, mvKeysRight, mDescriptorsLeft, mDescriptorsRight,
		vuRightTemp, vDepthTemp);


	cv::Mat add_mask = cv::Mat(cur_frame->cam0_img.rows, cur_frame->cam0_img.cols, CV_8UC1, cv::Scalar(255));

	for (int i = 0; i < mvKeysLeft.size(); i++)
	{


		if (add_mask.at<uchar>(mvKeysLeft[i].pt) == 255)
		{
			t_gframe::PointsAttri newPoint;
			newPoint.cam0_pt = mvKeysLeft[i].pt;
			newPoint.cam1_pt = mvKeysRight[i].pt;
			newPoint.depth = vDepthTemp[i];
			newPoint.id = npointID++;
			newPoint.track_cnt = 0;

			cur_frame->vpointInFrame.push_back(newPoint);
			cv::circle(add_mask, mvKeysLeft[i].pt, 30, 0, -1);
		}
	}

}


void gvis::t_gstereoorblk::featureDetect()
{
	mvKeysLeft.clear();
	mvKeysRight.clear();
	mDescriptorsLeft.release();
	mDescriptorsRight.release();


	thread threadLeft(&t_gstereoorblk::ExtractORB, this, cur_frame->cam0_img, 0);
	thread threadRight(&t_gstereoorblk::ExtractORB, this, cur_frame->cam1_img, 1);

	threadLeft.join();
	threadRight.join();

	if ((mvKeysLeft.empty()) || (mvKeysRight.empty()))
	{
		cerr << "ERROR: no Detected Points." << endl;
		exit(-1);
	}

}


void gvis::t_gstereoorblk::StereoMatch(vector<cv::KeyPoint> &KeyPtLeft,
	vector<cv::KeyPoint> &KeyPtRight, cv::Mat PtDescriptorsLeft,
	cv::Mat PtDescriptorsRight,vector<float> &uPtRight,
    vector<float> &ptDepth)
{


	int N = KeyPtLeft.size();


	uPtRight = vector<float>(N, -1.0f);
	ptDepth = vector<float>(N, -1.0f);
	vector<int> mvindexRight = vector<int>(N, -1);



	int mnScaleLevels = mpORBextractorLeft->GetLevels();
	float mfScaleFactor = mpORBextractorLeft->GetScaleFactor();
	float mfLogScaleFactor = log(mfScaleFactor);
	vector<float> mvScaleFactors = mpORBextractorLeft->GetScaleFactors();
	vector<float> mvInvScaleFactors = mpORBextractorLeft->GetInverseScaleFactors();
	vector<float> mvLevelSigma2 = mpORBextractorLeft->GetScaleSigmaSquares();
	vector<float> mvInvLevelSigma2 = mpORBextractorLeft->GetInverseScaleSigmaSquares();




	const int thOrbDist = (ORBmatcher::TH_HIGH + ORBmatcher::TH_LOW) / 2;

	const int nRows = mpORBextractorLeft->mvImagePyramid[0].rows;




	vector<vector<size_t>> vRowIndices(nRows, vector<size_t>());

	for (int i = 0; i < nRows; i++)
		vRowIndices[i].reserve(200);


	const int Nr = KeyPtRight.size();
	for (int iR = 0; iR < Nr; iR++)
	{

		const cv::KeyPoint &kp = KeyPtRight[iR];
		const float &kpY = kp.pt.y;


		const float r = 2.0f*mvScaleFactors[KeyPtRight[iR].octave];


		const int maxr = ceil(kpY + r);
		const int minr = floor(kpY - r);


		for (int yi = minr; yi <= maxr; yi++)
			vRowIndices[yi].push_back(iR);
	}

	
	//lsy
	
	const float minZ = mb;
	//const float minZ = 5;//认为10m内的点是不可信的
	const float minD = 0;                                 

	const float maxD = mbf / minZ;
	//const float maxD = 70;  

	
	//cout << "maxD:" << maxD << endl;

	vector<pair<int, int> > vDistIdx;
	vDistIdx.reserve(N);

	int count = 0;

	for (int iL = 0; iL < N; iL++)
	{
		//cout << "!!N!!:" << iL << endl;
		const cv::KeyPoint &kpL = KeyPtLeft[iL];
		const int &levelL = kpL.octave;

		const float &vL = kpL.pt.y;
		const float &uL = kpL.pt.x;


		const vector<size_t> &vCandidates = vRowIndices[vL];

		if (vCandidates.empty())
			continue;


		const float minU = uL - maxD;
		const float maxU = uL - minD;

		if (maxU < 0)
			continue;

		int bestDist = ORBmatcher::TH_HIGH;
		//cout << "bestDist:" << bestDist << endl;
	//	cout << "theorbdist:" << thOrbDist << endl;
		size_t bestIdxR = 0;

	
		const cv::Mat &dL = PtDescriptorsLeft.row(iL);

		//cout << "vCandidates:" << vCandidates.size() << endl;
		for (size_t iC = 0; iC < vCandidates.size(); iC++)
		{
			const size_t iR = vCandidates[iC];
			const cv::KeyPoint &kpR = KeyPtRight[iR];

			
			if (kpR.octave<levelL - 1 || kpR.octave>levelL + 1)
				continue;

			const float &uR = kpR.pt.x;


		
			if (uR >= minU && uR <= maxU)
			{
				const cv::Mat &dR = PtDescriptorsRight.row(iR);
			
				const int dist = ORBmatcher::DescriptorDistance(dL, dR);

				if (dist < bestDist)
				{
					bestDist = dist;
					bestIdxR = iR;
				}
			}
		}



		if (bestDist < thOrbDist)
		{
			
			mvindexRight[iL] = bestIdxR;


			const float uR0 = KeyPtRight[bestIdxR].pt.x;

		
			const float scaleFactor = mvInvScaleFactors[kpL.octave];
			const float scaleduL = round(kpL.pt.x*scaleFactor);
			const float scaledvL = round(kpL.pt.y*scaleFactor);

	
			const float scaleduR0 = round(uR0*scaleFactor);



			const int w = 5;

			
			cv::Mat IL = mpORBextractorLeft->mvImagePyramid[kpL.octave].rowRange(scaledvL - w, scaledvL + w + 1).colRange(scaleduL - w, scaleduL + w + 1);

			
			IL.convertTo(IL, CV_32F);
			IL = IL - IL.at<float>(w, w) *cv::Mat::ones(IL.rows, IL.cols, CV_32F);

			int bestDist = INT_MAX;
			int bestincR = 0;
		
			const int L = 5;
			vector<float> vDists;
			vDists.resize(2 * L + 1);



			const float iniu = scaleduR0 - L - w;
			const float endu = scaleduR0 + L + w + 1;
			if (iniu < 0 || endu >= mpORBextractorRight->mvImagePyramid[kpL.octave].cols)
				continue;

			
			for (int incR = -L; incR <= +L; incR++)
			{
				cv::Mat IR = mpORBextractorRight->mvImagePyramid[kpL.octave].rowRange(scaledvL - w, scaledvL + w + 1).colRange(scaleduR0 + incR - w, scaleduR0 + incR + w + 1);
				IR.convertTo(IR, CV_32F);


				IR = IR - IR.at<float>(w, w) *cv::Mat::ones(IR.rows, IR.cols, CV_32F);

			
				float dist = cv::norm(IL, IR, cv::NORM_L1);


				if (dist < bestDist)
				{
			
					bestDist = dist;
				
					bestincR = incR;
				}


				vDists[L + incR] = dist;
			}

		
			if (bestincR == -L || bestincR == L)
				continue;

			const float dist1 = vDists[L + bestincR - 1];
			const float dist2 = vDists[L + bestincR];
			const float dist3 = vDists[L + bestincR + 1];

			const float deltaR = (dist1 - dist3) / (2.0f*(dist1 + dist3 - 2.0f*dist2));

			
			if (deltaR < -1 || deltaR>1)
				continue;


			float bestuR = mvScaleFactors[kpL.octave] * ((float)scaleduR0 + (float)bestincR + deltaR);

	
			float disparity = (uL - bestuR);
			//cout << "disparity:" << disparity << endl;
		
			if (disparity >= minD && disparity < maxD)
			{
				if (disparity <= 0)
				{
					disparity = 0.01;
					bestuR = uL - 0.01;
				}

		
				ptDepth[iL] = mbf / disparity;
	
				uPtRight[iL] = bestuR;
				//cout << "bestuR:"<<bestuR << endl;

				vDistIdx.push_back(pair<int, int>(bestDist, iL));
			}
		}
	}



	sort(vDistIdx.begin(), vDistIdx.end());
	const float median = vDistIdx[vDistIdx.size() / 2].first;

	const float thDist = 1.5f*1.4f*median;
//	cout << "median:" << median << endl;
	//cout << "thedist:" << thDist << endl;

	for (int i = vDistIdx.size() - 1; i >= 0; i--)
	{
		if (vDistIdx[i].first < thDist)
			break;
		else
		{
			uPtRight[vDistIdx[i].second] = -1;
			ptDepth[vDistIdx[i].second] = -1;
		}
	}


	vector<uchar> status;
	status = vector<uchar>(N, 1);
	for (int i =0; i < ptDepth.size(); i++)
	{
		if (ptDepth[i] == -1)
		{
			status[i] = 0;
		}
	}


	reduceVector(KeyPtLeft, status);
	reduceVector(mvindexRight, status);
	reduceVector(uPtRight, status);
	reduceVector(ptDepth, status);

	vector<cv::KeyPoint> vKeysRightTemp;
	for (int i = 0; i < mvindexRight.size(); i++)
	{
		vKeysRightTemp.push_back(KeyPtRight[mvindexRight[i]]);
	}
	KeyPtRight = vKeysRightTemp;

}


void gvis::t_gstereoorblk::FeatureTracker()
{
	//cout << "pre_img feature Num: " << prev_frame->vpointInFrame.size() << endl;
	int before_tracking, after_tracking, after_ransac;
	before_tracking = prev_frame->vpointInFrame.size();

	Matx33f cam0_R_p_c;
	Matx33f cam1_R_p_c;
	integrateImuData(cam0_R_p_c, cam1_R_p_c);

	if (prev_frame->vpointInFrame.size() == 0)
		return;


	vector<Point2f> cam0_pre_pts(0);
	vector<Point2f> cam1_pre_pts(0);
	for (int i = 0; i < prev_frame->vpointInFrame.size(); i++)
	{
		cam0_pre_pts.push_back(prev_frame->vpointInFrame[i].cam0_pt);
		cam1_pre_pts.push_back(prev_frame->vpointInFrame[i].cam1_pt);

	}
	vector<Point2f> cam0_track_points(0);
	vector<Point2f> cam1_track_points(0);
	vector<uchar> cam0_LK_status(0);
	vector<uchar> cam1_LK_status(0);


	predictFeatureTracking(cam0_pre_pts,
		cam0_R_p_c, _cam0_intrinsics, cam0_track_points);

	predictFeatureTracking(cam1_pre_pts,
		cam1_R_p_c, _cam1_intrinsics, cam1_track_points);



	calcOpticalFlowPyrLK(
		prev_frame->cam0_pyramid, cur_frame->cam0_pyramid,
		cam0_pre_pts, cam0_track_points,
		cam0_LK_status, noArray(),
		Size(patch_size, patch_size),
		pyramid_levels,
		TermCriteria(TermCriteria::COUNT + TermCriteria::EPS,
			max_iteration,
			track_precision),
		cv::OPTFLOW_USE_INITIAL_FLOW);


	calcOpticalFlowPyrLK(
		prev_frame->cam1_pyramid, cur_frame->cam1_pyramid,
		cam1_pre_pts, cam1_track_points,
		cam1_LK_status, noArray(),
		Size(patch_size, patch_size),
		pyramid_levels,
		TermCriteria(TermCriteria::COUNT + TermCriteria::EPS,
			max_iteration,
			track_precision),
		cv::OPTFLOW_USE_INITIAL_FLOW);




	for (int i = 0; i < cam0_track_points.size(); ++i)
	{
		if ((cam0_LK_status[i] == 0) || (cam1_LK_status[i] == 0))
		{
			cam0_LK_status[i] = cam1_LK_status[i] = 0;
			continue;
		}
		//删除边界点	
		//lsy 
		if (cam0_track_points[i].y < 1 ||
			cam0_track_points[i].y > cur_frame->cam0_img.rows - 1 ||
			cam0_track_points[i].x < 1 ||
			cam0_track_points[i].x > cur_frame->cam0_img.cols - 1)
			cam0_LK_status[i] = cam1_LK_status[i] = 0;

		if (cam1_track_points[i].y < 1 ||
			cam1_track_points[i].y > cur_frame->cam1_img.rows - 1 ||
			cam1_track_points[i].x < 1 ||
			cam1_track_points[i].x > cur_frame->cam1_img.cols - 1)
			cam0_LK_status[i] = cam1_LK_status[i] = 0;
	}


	reduceVector(cam0_track_points, cam0_LK_status);
	reduceVector(cam1_track_points, cam0_LK_status);
	reduceVector(cam0_pre_pts, cam0_LK_status);
	reduceVector(cam1_pre_pts, cam0_LK_status);
	reduceVector(prev_frame->vpointInFrame, cam0_LK_status);


	

	vector<uchar> cam0_Fundamental_inliers(0);
//	cv::findFundamentalMat(cam0_track_points, cam1_track_points, cv::FM_RANSAC, 1, 0.99, cam0_Fundamental_inliers);
	//cout << cam0_track_points.at(0).x << " " << cam0_track_points.at(0).y << endl;
	for (int i = 0; i < cam0_track_points.size(); i++)
	{
		//左右 纵向差6 就剔除？
		//lsy:
		if(fabs(cam0_track_points[i].y-cam1_track_points[i].y)<6)
		//if (fabs(cam0_track_points[i].y - cam1_track_points[i].y) < 6)
			cam0_Fundamental_inliers.push_back(1);
		else
			cam0_Fundamental_inliers.push_back(0);
	}

	reduceVector(cam0_track_points, cam0_Fundamental_inliers);
	reduceVector(cam1_track_points, cam0_Fundamental_inliers);
	reduceVector(cam0_pre_pts, cam0_Fundamental_inliers);
	reduceVector(cam1_pre_pts, cam0_Fundamental_inliers);
	reduceVector(prev_frame->vpointInFrame, cam0_Fundamental_inliers);


	after_tracking = cam0_track_points.size();
	//cout << "before ransac points Num: " << cam0_track_points.size() << endl;

	

	//9.执行双目ransac剔除前后两帧跟踪错误的点
	vector<uchar> cam0_ransac_inliers(0);
	twoPointRansac(cam0_pre_pts, cam0_track_points,
		cam0_R_p_c, _cam0_intrinsics, _cam0_distortion_model,
		_cam0_distortion_coeffs, ransac_threshold,
		0.99, cam0_ransac_inliers);

	vector<uchar> cam1_ransac_inliers(0);
	twoPointRansac(cam1_pre_pts, cam1_track_points,
		cam1_R_p_c, _cam1_intrinsics, _cam1_distortion_model,
		_cam1_distortion_coeffs, ransac_threshold,
		0.99, cam1_ransac_inliers);

	// 9.移除所有ransac后的外点
	for (int i = 0; i < cam0_track_points.size(); ++i)
	{
		if ((cam0_ransac_inliers[i] == 0) || (cam1_ransac_inliers[i] == 0))
		{
			cam0_ransac_inliers[i] = cam1_ransac_inliers[i] = 0;
			continue;
		}
	}

	reduceVector(cam0_track_points, cam0_ransac_inliers);
	reduceVector(cam1_track_points, cam0_ransac_inliers);
	reduceVector(cam0_pre_pts, cam0_Fundamental_inliers);
	reduceVector(cam1_pre_pts, cam0_Fundamental_inliers);
	reduceVector(prev_frame->vpointInFrame, cam0_ransac_inliers);

	after_ransac = cam0_track_points.size();
	printf("candidates: %d; track: %d; ransac: %d/%d=%f; after_ransac：%d \n",
		before_tracking, after_tracking, after_ransac,after_tracking,
		static_cast<double>(after_ransac) /
		(static_cast<double>(before_tracking) + 1e-5), after_ransac);
	cout << "==================================" << endl;

	//这里对cam0_track_points 和cam1_track_points进行取整操作
	/*for (auto& cam0_track_point : cam0_track_points)
	{
		int a= int(cam0_track_point.x + 0.5);
		int b= int(cam0_track_point.y + 0.5);
		cam0_track_point.x = a;
		cam0_track_point.y = b;
		
	}
	for (auto& cam1_track_point : cam1_track_points)
	{
		int a = int(cam1_track_point.x + 0.5);
		int b = int(cam1_track_point.y + 0.5);
		cam1_track_point.x = a;
		cam1_track_point.y = b;
	}*/

	//重新计算描述子
	
	//vector<cv::KeyPoint> cam0_kpts, cam1_kpts;
	////把cam0_tracked_points转为cam0_pt_left
	//for (auto cam0_track_point : cam0_track_points)
	//{
	//	cv::KeyPoint temp_cam0_kpt;
	//	//cout << cam0_track_point.x << "," << cam0_track_point.y << endl;
	//	temp_cam0_kpt.pt.x = cam0_track_point.x;
	//	temp_cam0_kpt.pt.y = cam0_track_point.y;
	//	cam0_kpts.emplace_back(temp_cam0_kpt);
	//}
	//for (auto cam1_track_point : cam1_track_points)
	//{
	//	cv::KeyPoint temp_cam1_kpt;
	//	//cout << cam1_track_point.x << "," << cam1_track_point.y << endl;
	//	temp_cam1_kpt.pt.x = cam1_track_point.x;
	//	temp_cam1_kpt.pt.y = cam1_track_point.y;
	//	cam1_kpts.emplace_back(temp_cam1_kpt);
	//}

	//cv::Ptr<cv::ORB> orb = cv::ORB::create();
	//vector<float> uRight_temp, u_ptDepth_temp;
	//cv::Mat img0_descriptor, img1_descriptor;
	
	//cout<<"计算描述子之前："<< cam0_kpts.size() << "," << cam1_kpts.size() << endl;

	//temp
	//vector<cv::KeyPoint> temp_cam0, temp_cam1,d_temp,d_temp_1;
	//temp_cam0 = cam0_kpts; temp_cam1 = cam1_kpts;
	//bool fla = false;

	//计算描述子
	//orb->compute(cur_frame->cam0_img, cam0_kpts, img0_descriptor);
	//orb->compute(cur_frame->cam1_img, cam1_kpts, img1_descriptor);
	//cout <<"before cal descriptor:"<< cam0_kpts.size()<<","<<cam1_kpts.size() << endl;
	//if (cam0_kpts.size() != cam1_kpts.size())
	//{
	//	cout << "边界错误：" << endl;
	//	cin.get();
	//}
	//for (int i = 0; i < temp_cam0.size(); i++)
	//{
	//	fla = false;
	//	//cout << "match 前：" << cam0_kpts_before.at(i).pt.x << endl;
	//	for (int j = 0; j < cam0_kpts.size(); j++)
	//	{
	//		if (fabs(temp_cam0.at(i).pt.x - cam0_kpts.at(j).pt.x) < 1e-6
	//			&&fabs(temp_cam0.at(i).pt.y - cam0_kpts.at(j).pt.y) < 1e-6)
	//		{
	//			fla= true;
	//		}
	//	}
	//	if (fla == false)
	//	{
	//		//cout << "无:" << cam0_kpts_before.at(i).pt.x << endl;
	//		d_temp.emplace_back(temp_cam0.at(i));
	//		d_temp_1.emplace_back(temp_cam1.at(i));
	//		//d_cam1_kpts.emplace_back(cam1_kpts_before.at(i));
	//	}
	//}

	//for (int i = 0; i < d_temp.size(); i++)
	//{
	//	cout << d_temp.at(i).pt.x << "," << d_temp.at(i).pt.y <<","<< d_temp_1.at(i).pt.x <<","<< d_temp_1.at(i).pt.y <<endl;
	//}

	//vector<cv::KeyPoint> cam0_kpts_before, cam1_kpts_before,d_cam0_kpts,d_cam1_kpts;
	//cam0_kpts_before = cam0_kpts;
	//cam1_kpts_before = cam1_kpts;
	//StereoMatch(cam0_kpts, cam1_kpts,img0_descriptor,img1_descriptor,uRight_temp, u_ptDepth_temp);
	//cout << "after_stereoMatch:" << cam0_kpts.size() << "," << cam1_kpts.size() << endl;
	//bool flag = false;
	//for (int i = 0; i < cam0_kpts_before.size(); i++)
	//{
	//	flag = false;
	//	//cout << "match 前：" << cam0_kpts_before.at(i).pt.x << endl;
	//	for (int j = 0; j < cam0_kpts.size(); j++)
	//	{
	//		if (fabs(cam0_kpts_before.at(i).pt.x - cam0_kpts.at(j).pt.x) < 1e-6
	//			&&fabs(cam0_kpts_before.at(i).pt.y - cam0_kpts.at(j).pt.y) < 1e-6)
	//		{
	//			flag = true;
	//		}
	//	}
	//	if (flag == false)
	//	{
	//		//cout << "无:" << cam0_kpts_before.at(i).pt.x << endl;
	//		d_cam0_kpts.emplace_back(cam0_kpts_before.at(i));
	//		d_cam1_kpts.emplace_back(cam1_kpts_before.at(i));
	//	}
	//}

	/*for (int i = 0; i < d_cam0_kpts.size(); i++)
	{
		cout << d_cam0_kpts.at(i).pt.x << "," << d_cam0_kpts.at(i).pt.y<<","<<d_cam1_kpts.at(i).pt.x<<","<<d_cam1_kpts.at(i).pt.y << endl;
	}
	cout << endl;

	cout << cam0_kpts_before.size() << endl;
	cout << cam0_kpts.size() << endl;
	cout << d_cam0_kpts.size() << endl;*/
	//画图
	//int rows = cur_frame->cam0_img.rows;
	//int cols = cur_frame->cam0_img.cols;
	//Mat out_img(rows, cols * 2, CV_8UC3);
	//cvtColor(cur_frame->cam0_img, out_img.colRange(0, cols), cv::COLOR_GRAY2BGR);
	//cvtColor(cur_frame->cam1_img, out_img.colRange(cols, cols * 2), cv::COLOR_GRAY2BGR);
	//Mat out_img_after(rows, cols * 2, CV_8UC3);
	//cvtColor(cur_frame->cam0_img, out_img_after.colRange(0, cols), cv::COLOR_GRAY2BGR);
	//cvtColor(cur_frame->cam1_img, out_img_after.colRange(cols, cols * 2), cv::COLOR_GRAY2BGR);
	//// 画栅格线
	//int height = rows / 30;

	//for (int i = 1; i < 30; i++)
	//{
	//	line(out_img, Point(0, height*i), Point(cols * 2, height*i), Scalar(0, 255, 0), 1.8);
	//	line(out_img_after, Point(0, height*i), Point(cols * 2, height*i), Scalar(0, 255, 0), 1.8);
	//}

	//for (int i = 0; i < d_cam0_kpts.size(); i++)
	//{
	//	
	//	cv::Point2f pt0 =d_cam0_kpts.at(i).pt;
	//	circle(out_img, pt0, 4, Scalar(0, 0, 255), 2);
	//	string name;
	//	//lsy
	//	name = to_string(i + 1);
	//	cv::putText(out_img, name, pt0, cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255, 0, 0));
	//	cv::Point2f pt1 =d_cam1_kpts.at(i).pt + Point2f(cols, 0.0);
	//	circle(out_img, pt1, 4, Scalar(255, 0, 255), 2);
	//	cv::putText(out_img, name, pt1, cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255, 0, 0));
	//}
	//imwrite("d_img.jpg", out_img);
	//cin.get();

	cur_frame->vpointInFrame = prev_frame->vpointInFrame;
	for (int i = 0; i < cur_frame->vpointInFrame.size(); i++)
	{
		//whd
		//cur_frame->vpointInFrame[i].cam0_pt = cam0_kpts[i].pt;
		//cur_frame->vpointInFrame[i].cam1_pt = cam1_kpts[i].pt;
		cur_frame->vpointInFrame[i].cam0_pt = cam0_track_points[i];
		cur_frame->vpointInFrame[i].cam1_pt = cam1_track_points[i];
		cur_frame->vpointInFrame[i].track_cnt++;
	}

	
		    
}



void gvis::t_gstereoorblk::addNewFeatures()
{
	int addNum = max_cnt - static_cast<int>(cur_frame->vpointInFrame.size());

	if (addNum > 0)
	{


		featureDetect();

		//// 设置mask防止跟踪点被重复检测,只设置左侧影像即可
		//cv::Mat mask = setmask(cam0_cur_pts);
		//vector<uchar> status = vector<uchar>(mvKeysLeft.size(), 1);
		//for (int i = 0; i < mvKeysLeft.size(); i++)
		//{
		//	if (mask.at<uchar>(mvKeysLeft[i].pt) == 0)
		//	{
		//		status[i] = 0;
		//	}
		//}
		//reduceVector(mvKeysLeft, status);
		//// 删除描述子
		//cv::Mat DescriptorsLeftTemp;
		//for (int i = 0; i < status.size(); i++)
		//{
		//	if (status[i])
		//	{
		//		DescriptorsLeftTemp.push_back(mDescriptorsLeft.rowRange(i,i+1));
		//	}
		//}
		//mDescriptorsLeft = DescriptorsLeftTemp.clone();



		std::vector<float> vuRightTemp;
		std::vector<float> vDepthTemp;

		StereoMatch(mvKeysLeft, mvKeysRight, mDescriptorsLeft, mDescriptorsRight,
			vuRightTemp, vDepthTemp);


		vector<t_gframe::PointsAttri> vnewPoints;


		for (int i = 0; i < mvKeysLeft.size(); i++)
		{
			t_gframe::PointsAttri  newPoint;
			newPoint.cam0_pt = mvKeysLeft[i].pt;
			newPoint.cam1_pt = mvKeysRight[i].pt;
			newPoint.depth = vDepthTemp[i];
			newPoint.id = npointID++;
			newPoint.track_cnt = 0;

			vnewPoints.push_back(newPoint);
		}


		vnewPoints.insert(vnewPoints.end(), cur_frame->vpointInFrame.begin(), cur_frame->vpointInFrame.end());


		sort(vnewPoints.begin(), vnewPoints.end(), featureCompareByTrack);



		cv::Mat add_mask = cv::Mat(cur_frame->cam0_img.rows, cur_frame->cam0_img.cols, CV_8UC1, cv::Scalar(255));

		int addCnt = 0;

		cur_frame->vpointInFrame.clear();
		for (int i = 0; i < vnewPoints.size(); i++)
		{
			if (addCnt > addNum)
				break;


			if (add_mask.at<uchar>(vnewPoints[i].cam0_pt) == 255)
			{
				cur_frame->vpointInFrame.push_back(vnewPoints[i]);
				addCnt++;


				cv::circle(add_mask, vnewPoints[i].cam0_pt, 30, 0, -1);
			}
		}

	}
	

}



void gvis::t_gstereoorblk::publish()
{
	_pointCloud.features.clear();
	FeaturePoint feature_msg;

	vector<cv::Point2f> norm_cur0_pts(cur_frame->vpointInFrame.size());
	vector<cv::Point2f> norm_cur1_pts(cur_frame->vpointInFrame.size());



	for (int i = 0; i < cur_frame->vpointInFrame.size(); i++)
	{
		Eigen::Vector3d tmp_p1, tmp_p2;

		Normalised(_cam0_intrinsics, Eigen::Vector2d(cur_frame->vpointInFrame[i].cam0_pt.x, cur_frame->vpointInFrame[i].cam0_pt.y), tmp_p1);
		Eigen::Vector3d p1_new = Eigen::Vector3d(tmp_p1.x(), tmp_p1.y(), 1);
		norm_cur0_pts[i] = cv::Point2f(p1_new(0), p1_new(1));



		Normalised(_cam1_intrinsics, Eigen::Vector2d(cur_frame->vpointInFrame[i].cam1_pt.x, cur_frame->vpointInFrame[i].cam1_pt.y), tmp_p2);
		Eigen::Vector3d p2_new = Eigen::Vector3d(tmp_p2.x(), tmp_p2.y(), 1);
		norm_cur1_pts[i] = cv::Point2f(p2_new(0), p2_new(1));
	}


	_pointCloud.time = cur_frame->time;


	for (int i = 0; i < norm_cur0_pts.size(); i++)
	{

		feature_msg.id = cur_frame->vpointInFrame[i].id;
		feature_msg.cam0_point = Point2f(norm_cur0_pts[i].x,norm_cur0_pts[i].y);
		feature_msg.cam1_point = Point2f(norm_cur1_pts[i].x, norm_cur1_pts[i].y);
		feature_msg.depth = cur_frame->vpointInFrame[i].depth;

		_pointCloud.features.push_back(feature_msg);
	}

}

//void gvis::t_gstereoorblk::drawoneImage()
//{
//	int rows = cur_frame->cam0_img.rows;
//	int cols = cur_frame->cam0_img.cols;
//	Mat out_img(rows, cols, CV_8UC3);
//	cvtColor(cur_frame->cam0_img, out_img.colRange(0, cols), cv::COLOR_GRAY2BGR);
//	for (int i = 0; i < cur_frame->vpointInFrame.size(); i++)
//	{
//
//		cv::Point2f pt0 = cur_frame->vpointInFrame[i].cam0_pt;
//		circle(out_img, pt0, 4, Scalar(0, 0, 255), 2);
//		//string name;
//		//lsy
//		//name = to_string(cur_frame->vpointInFrame[i].depth);
//		//name = to_string(cur_frame->vpointInFrame[i].id) + ", " + to_string(cur_frame->vpointInFrame[i].depth);
//		//name = to_string(track_cnt[i]) + ", " + to_string(mvDepth[i]);
//		//name = to_string(cur_frame->vpointInFrame[i].track_cnt);
//		//cv::putText(out_img, name, pt0, cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255, 0, 0));
//
//
//		//cv::Point2f pt1 = cur_frame->vpointInFrame[i].cam1_pt + Point2f(cols, 0.0);
//		//circle(out_img, pt1, 4, Scalar(255, 0, 255), 2);
//		//cv::putText(out_img, name, pt1, cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255, 0, 0));
//
//
//		//line(out_img, pt0, pt1, Scalar(255, 0, 0), 1);
//	}
//	resize(out_img, out_img, cv::Size(out_img.cols, out_img.rows));
//	//lsy
//	//imwrite("Feature.jpg", out_img);
//	//cin.get();
//	//cv::resize(out_img, out_img, cv::Size(out_img.cols / 2.0, out_img.rows / 2.0));
//	imshow("Feature", out_img);
//	waitKey(5);
//	return;
//}


void gvis::t_gstereoorblk::drawImage()
{
	int rows = cur_frame->cam0_img.rows;
	int cols = cur_frame->cam0_img.cols;
	Mat out_img(rows, cols * 2, CV_8UC3);
	cvtColor(cur_frame->cam0_img, out_img.colRange(0, cols), cv::COLOR_GRAY2BGR);
	cvtColor(cur_frame->cam1_img, out_img.colRange(cols, cols * 2), cv::COLOR_GRAY2BGR);


	// 画格网
	//for (int i = 1; i < grid_row; i++)
	//{
	//	line(out_img, Point(0,grid_height*i ), Point(cols*2, grid_height*i), Scalar(255, 0, 0), 1.8);
	//}

	//for (int i = 1; i < grid_col*2; i++)
	//{
	//	line(out_img, Point(grid_width*i, 0), Point(grid_width*i,rows), Scalar(255, 0, 0), 1.8);
	//}

	// 画栅格线
	int height = rows / 30;

	for (int i = 1; i < 30; i++)
	{
		line(out_img, Point(0, height*i), Point(cols * 2, height*i), Scalar(0, 255, 0), 1.8);
	}


	for (int i=0; i < cur_frame->vpointInFrame.size(); i++)
	{

		cv::Point2f pt0 = cur_frame->vpointInFrame[i].cam0_pt;
		circle(out_img, pt0, 4, Scalar(0, 0, 255), 2);
		string name;
		//lsy
		//name = to_string(cur_frame->vpointInFrame[i].depth);
		name = to_string(cur_frame->vpointInFrame[i].id) + ", " +to_string(cur_frame->vpointInFrame[i].depth);
		//name = to_string(track_cnt[i]) + ", " + to_string(mvDepth[i]);
		//name = to_string(cur_frame->vpointInFrame[i].track_cnt);
		cv::putText(out_img, name, pt0, cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255, 0, 0));


		cv::Point2f pt1 = cur_frame->vpointInFrame[i].cam1_pt + Point2f(cols, 0.0);
		circle(out_img, pt1, 4, Scalar(255, 0, 255), 2);
		cv::putText(out_img, name, pt1, cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255, 0, 0));


		//line(out_img, pt0, pt1, Scalar(255, 0, 0), 1);
	}




	resize(out_img, out_img, cv::Size(out_img.cols, out_img.rows));
	//lsy
	//imwrite("Feature.jpg", out_img);
	//cin.get();
	cv::resize(out_img, out_img, cv::Size(out_img.cols / 2.0, out_img.rows / 2.0));
	imshow("feature", out_img);
	waitKey(5);
	return;

}




//***************************************************************************************
cv::Mat gvis::t_gstereoorblk::setmask(std::vector<cv::Point2f> &point)
{
	// mask初值255
	cv::Mat mask = cv::Mat(cur_frame->cam0_img.rows, cur_frame->cam0_img.cols, CV_8UC1, cv::Scalar(255));


	for (int i = 0; i < point.size(); i++)
	{
		if (mask.at<uchar>(point[i]) == 255)
		{
			// 已有像素点圆内画为0
			// 现有点周围15个像素之内不能取得特征点，这个阈值放小使得有更多点参与匹配
			cv::circle(mask, point[i], 15, 0, -1);
		}
	}

	return mask;
}



template <typename vecType>
void gvis::t_gstereoorblk::reduceVector(vector<vecType> &v, vector<uchar>  &status)
{
	int j = 0;
	for (int i = 0; i < int(v.size()); i++)
	{
		if (status[i])
			v[j++] = v[i];
	}

	v.resize(j);
}



void gvis::t_gstereoorblk::Normalised(cv::Vec4d intrinsics,const Eigen::Vector2d& p, Eigen::Vector3d& P)
{
	
	double mx_d, my_d, mx2_d, mxy_d, my2_d, mx_u, my_u;
	double rho2_d, rho4_d, radDist_d, Dx_d, Dy_d, inv_denom_d;
	double m_inv_K11, m_inv_K13, m_inv_K22, m_inv_K23;

	m_inv_K11 = 1.0 / intrinsics(0);
	m_inv_K13 = -intrinsics(2) / intrinsics(0);
	m_inv_K22 = 1.0 / intrinsics(1);
	m_inv_K23 = -intrinsics(3) / intrinsics(1);


	// Lift points to normalised plane
	mx_d = m_inv_K11 * p(0) + m_inv_K13;
	my_d = m_inv_K22 * p(1) + m_inv_K23;


	mx_u = mx_d;
	my_u = my_d;

	P << mx_u, my_u, 1.0;
}


bool gvis::t_gstereoorblk::inBorder(const cv::Point2f &pt)
{
	const int BORDER_SIZE = 1;
	int pt_x = cvRound(pt.x);
	int pt_y = cvRound(pt.y);
	return (BORDER_SIZE <= pt_x ) && (pt_x < cur_frame->cam0_img.cols - BORDER_SIZE)
		&& (BORDER_SIZE <= pt_y) && (pt_y < cur_frame->cam0_img.rows - BORDER_SIZE);
}



void gvis::t_gstereoorblk::RectifyImg()
{
	cv::Mat imLeftRect, imRightRect;
	
	RectWholeImg(cur_frame->cam0_img, imLeftRect, _K_L, _D_L, _R_L, _P_L);

	RectWholeImg(cur_frame->cam1_img, imRightRect, _K_R, _D_R, _R_R, _P_R);
	
	cur_frame->cam0_img = imLeftRect.clone();
	cur_frame->cam1_img = imRightRect.clone();
}


void gvis::t_gstereoorblk::RectWholeImg(const cv::Mat &img, cv::Mat &UndistortImg, const cv::Mat &K0, const cv::Mat &D0, const cv::Mat &R0, const cv::Mat &P0)
{
	cv::Mat TempM1, TempM2;
	
	cv::initUndistortRectifyMap(K0, D0, R0, P0.rowRange(0, 3).colRange(0, 3), img.size(), CV_32F, TempM1, TempM2);

	cv::remap(img, UndistortImg, TempM1, TempM2, cv::INTER_LINEAR);
}



void gvis::t_gstereoorblk::Epipolar_Constraint(vector<Point2f> &cam0_track_points,
	vector<Point2f> &cam1_track_points, vector<uchar> &status)
{


	const cv::Matx33d R_cam0_cam1 =_R_cam1_imu.t() * _R_cam0_imu;
	const cv::Vec3d t_cam0_cam1 = _R_cam1_imu.t() * (_t_cam0_imu - _t_cam1_imu);


	const cv::Matx33d t_cam0_cam1_hat(
		0.0, -t_cam0_cam1[2], t_cam0_cam1[1],
		t_cam0_cam1[2], 0.0, -t_cam0_cam1[0],
		-t_cam0_cam1[1], t_cam0_cam1[0], 0.0);
	const cv::Matx33d E = t_cam0_cam1_hat * R_cam0_cam1;

	double norm_pixel_unit = 4.0 / (
		_cam0_intrinsics[0] + _cam0_intrinsics[1] +
		_cam1_intrinsics[0] + _cam1_intrinsics[1]);

	cv::Vec3d pt0, pt1;

	for (unsigned int i = 0; i < cam0_track_points.size(); i++)
	{
		if (status[i] == 1)
		{
		
			Eigen::Vector3d tmp_p1, tmp_p2;
			Normalised(_cam0_intrinsics, Eigen::Vector2d(cam0_track_points[i].x, cam0_track_points[i].y), tmp_p1);
			Normalised(_cam1_intrinsics, Eigen::Vector2d(cam1_track_points[i].x, cam1_track_points[i].y), tmp_p2);
			pt0 = cv::Vec3d(tmp_p1.x(), tmp_p1.y(), 1);
			pt1 = cv::Vec3d(tmp_p2.x(), tmp_p2.y(), 1);

		
			cv::Vec3d epipolar_line = E * pt0;

	
			double error = fabs((pt1.t() * epipolar_line)[0]) / sqrt(
				epipolar_line[0] * epipolar_line[0] +
				epipolar_line[1] * epipolar_line[1]);

			if (error > stereo_threshold * norm_pixel_unit)
				status[i] = 0;

		}
	}
}





void gvis::t_gstereoorblk::integrateImuData(
	cv::Matx33f& cam0_R_p_c, cv::Matx33f& cam1_R_p_c)
{

	//
	//Vec3f mean_ang_vel(0.0, 0.0, 0.0);

	//for (auto &iter : vecimu)
	//{

	//	mean_ang_vel += Vec3f(iter.angular_velocity.x(),
	//		iter.angular_velocity.y(), iter.angular_velocity.z());

	//}
	//mean_ang_vel *= 1.0f / vecimu.size();

	//Vec3f cam0_mean_ang_vel = _R_cam0_imu.t() * mean_ang_vel;
	//Vec3f cam1_mean_ang_vel = _R_cam1_imu.t() * mean_ang_vel;

	//double dtime = (cur_frame->time -prev_frame->time);
	//Rodrigues(cam0_mean_ang_vel*dtime, cam0_R_p_c);
	//Rodrigues(cam1_mean_ang_vel*dtime, cam1_R_p_c);
	//cam0_R_p_c = cam0_R_p_c.t();
	//cam1_R_p_c = cam1_R_p_c.t();


	//vecimu.clear();
	//return;


	auto begin_iter = _vecimu.begin();
	while (begin_iter != _vecimu.end())
	{
		if (begin_iter->t -
			prev_frame->time < -0.01)
			++begin_iter;
		else
			break;
	}

	auto end_iter = begin_iter;
	while (end_iter != _vecimu.end())
	{
		/*if (end_iter->t -
			curr_img_msg.t < 0.005)*/
			//修改
		if (end_iter->t -
			cur_frame->time < 0.002)
			++end_iter;
		else
			break;
	}
	//cout << "imu_size:" << _vecimu.size() << endl;
	//cout << "prev_img:" << prev_frame->time << endl;
	//cout << "cur_img:" << cur_frame->time << endl;
	//for (auto imu : _vecimu)
	//{
	//	cout << "imu_t:" << imu.t << endl;
	//}
	//cout << "use_imu_size:" << end_iter - begin_iter << endl;
	//cout << "prev_img_t:" <<setprecision(10)<< prev_img_msg.t << endl;
	Vec3f mean_ang_vel(0.0, 0.0, 0.0);
	for (auto iter = begin_iter; iter < end_iter; ++iter)

	{
		mean_ang_vel += Vec3f(iter->angular_velocity.x(),
			iter->angular_velocity.y(), iter->angular_velocity.z());
		//cout << "imu_t:" << setprecision(10) << iter->t << "  " << iter->angular_velocity.transpose() << endl;
	}
	//cout << "cur_img_t:" << setprecision(10) << curr_img_msg.t << endl;

	if (end_iter - begin_iter > 0)
		mean_ang_vel *= 1.0f / (end_iter - begin_iter);
	//cout << "end_iter-begin_iter:" << end_iter - begin_iter << endl;

	Vec3f cam0_mean_ang_vel = _R_cam0_imu.t() * mean_ang_vel;
	Vec3f cam1_mean_ang_vel = _R_cam1_imu.t() * mean_ang_vel;


	double dtime = (cur_frame->time - prev_frame->time);
	//cout << "dtime:" << setprecision(10)<<dtime << endl;
	Rodrigues(cam0_mean_ang_vel*dtime, cam0_R_p_c);
	Rodrigues(cam1_mean_ang_vel*dtime, cam1_R_p_c);
	cam0_R_p_c = cam0_R_p_c.t();
	cam1_R_p_c = cam1_R_p_c.t();
	/*for (auto imu : vecimu)
		cout << "befroe_imu_time:" << setprecision(10) << imu.t << endl;*/


	_vecimu.erase(_vecimu.begin(), end_iter);
	//for (auto imu : vecimu)
	//	cout << "after_imu_time:" << setprecision(10) << imu.t << endl;

	return;
}




void gvis::t_gstereoorblk::predictFeatureTracking(
	const vector<cv::Point2f>& input_pts,
	const cv::Matx33f& R_p_c,
	const cv::Vec4d& intrinsics,
	vector<cv::Point2f>& compensated_pts)
{

	
	if (input_pts.size() == 0)
	{
		compensated_pts.clear();
		return;
	}
	compensated_pts.resize(input_pts.size());


	cv::Matx33f K(
		intrinsics[0], 0.0, intrinsics[2],
		0.0, intrinsics[1], intrinsics[3],
		0.0, 0.0, 1.0);



	cv::Matx33f H = K * R_p_c * K.inv();

	for (int i = 0; i < input_pts.size(); ++i)
	{
		cv::Vec3f p1(input_pts[i].x, input_pts[i].y, 1.0f);

		cv::Vec3f p2 = H * p1;
		compensated_pts[i].x = p2[0] / p2[2];
		compensated_pts[i].y = p2[1] / p2[2];
	}

	return;
}


//创建当前观测影像的金字塔用于光流跟踪
void gvis::t_gstereoorblk::createImagePyramids()
{
	

	buildOpticalFlowPyramid(
		cur_frame->cam0_img, cur_frame->cam0_pyramid,
		cv::Size(patch_size, patch_size),
		pyramid_levels, true, cv::BORDER_REFLECT_101,
		cv::BORDER_CONSTANT, false);

	buildOpticalFlowPyramid(
		cur_frame->cam1_img, cur_frame->cam1_pyramid,
		cv::Size(patch_size, patch_size),
		pyramid_levels, true, cv::BORDER_REFLECT_101,
		cv::BORDER_CONSTANT, false);


	mpORBextractorLeft->mPyramid = cur_frame->cam0_pyramid;
	mpORBextractorRight->mPyramid = cur_frame->cam1_pyramid;
}





void gvis::t_gstereoorblk::twoPointRansac(
	const vector<cv::Point2f>& pts1, const vector<cv::Point2f>& pts2,
	const cv::Matx33f& R_p_c, const cv::Vec4d& intrinsics,
	const std::string& distortion_model,
	const cv::Vec4d& distortion_coeffs,
	const double& inlier_error,
	const double& success_probability,
	vector<uchar>& inlier_markers)
{

	// Check the size of input point size.
	if (pts1.size() != pts2.size())
		printf("Sets of different size (%d and %d) are used...\n",
			pts1.size(), pts2.size());


	double norm_pixel_unit = 2.0 / (intrinsics[0] + intrinsics[1]);


	int iter_num = static_cast<int>(
		ceil(log(1 - success_probability) / log(1 - 0.7*0.7)));


	inlier_markers.clear();
	inlier_markers.resize(pts1.size(), 1);


	
	vector<Point2f> pts1_undistorted;
	vector<Point2f> pts2_undistorted;

	for (int i = 0; i < pts1.size(); i++)
	{
		Eigen::Vector3d pts1_norm;
		Eigen::Vector3d pts2_norm;

		Normalised(_cam0_intrinsics, Eigen::Vector2d(pts1[i].x, pts1[i].y), pts1_norm);
		Normalised(_cam1_intrinsics, Eigen::Vector2d(pts2[i].x, pts2[i].y), pts2_norm);

		pts1_undistorted.push_back(Point2f(pts1_norm.x(), pts1_norm.y()));
		pts2_undistorted.push_back(Point2f(pts2_norm.x(), pts2_norm.y()));
	}


	//undistortPoints(
	//	pts1, intrinsics, distortion_model,
	//	distortion_coeffs, pts1_undistorted);
	//undistortPoints(
	//	pts2, intrinsics, distortion_model,
	//	distortion_coeffs, pts2_undistorted);



	for (auto& pt : pts1_undistorted)
	{
		Vec3f pt_h(pt.x, pt.y, 1.0f);
		Vec3f pt_hc = R_p_c * pt_h;
		pt.x = pt_hc[0];
		pt.y = pt_hc[1];
	}


	float scaling_factor = 0.0f;
	rescalePoints(pts1_undistorted, pts2_undistorted, scaling_factor);


	norm_pixel_unit *= scaling_factor;




	vector<Point2d> pts_diff(pts1_undistorted.size());
	for (int i = 0; i < pts1_undistorted.size(); ++i)
		pts_diff[i] = pts1_undistorted[i] - pts2_undistorted[i];





	double mean_pt_distance = 0.0;
	int raw_inlier_cntr = 0;
	for (int i = 0; i < pts_diff.size(); ++i)
	{
		double distance = sqrt(pts_diff[i].dot(pts_diff[i]));


		if (distance > 50.0*norm_pixel_unit)
		{
			inlier_markers[i] = 0;
		}
		else
		{
			mean_pt_distance += distance;
			++raw_inlier_cntr;
		}
	}

	mean_pt_distance /= raw_inlier_cntr;


	if (raw_inlier_cntr < 3)
	{
		for (auto& marker : inlier_markers)
			marker = 0;
		return;
	}



	if (mean_pt_distance < norm_pixel_unit)
	{
		printf("Degenerated motion...\n");
		for (int i = 0; i < pts_diff.size(); ++i)
		{
			if (inlier_markers[i] == 0)
				continue;


			if (sqrt(pts_diff[i].dot(pts_diff[i])) >
				inlier_error*norm_pixel_unit)
				inlier_markers[i] = 0;
		}
		return;
	}


	MatrixXd coeff_t(pts_diff.size(), 3);
	for (int i = 0; i < pts_diff.size(); ++i)
	{
		coeff_t(i, 0) = pts_diff[i].y;
		coeff_t(i, 1) = -pts_diff[i].x;
		coeff_t(i, 2) = pts1_undistorted[i].x*pts2_undistorted[i].y -
			pts1_undistorted[i].y*pts2_undistorted[i].x;
	}



	vector<int> raw_inlier_idx;
	for (int i = 0; i < inlier_markers.size(); ++i)
	{
		if (inlier_markers[i] != 0)
	
			raw_inlier_idx.push_back(i);
	}


	vector<int> best_inlier_set;
	double best_error = 1e10;

	cv::RNG randomNumber;


	for (int iter_idx = 0; iter_idx < iter_num; ++iter_idx)
	{

		int pair_idx1 = raw_inlier_idx[randomNumber.uniform(
			0, raw_inlier_idx.size())];

		int idx_diff = randomNumber.uniform(
			1, raw_inlier_idx.size());


		int pair_idx2 = pair_idx1 + idx_diff < raw_inlier_idx.size() ?
			pair_idx1 + idx_diff : pair_idx1 + idx_diff - raw_inlier_idx.size();


		Vector2d coeff_tx(coeff_t(pair_idx1, 0), coeff_t(pair_idx2, 0));
		Vector2d coeff_ty(coeff_t(pair_idx1, 1), coeff_t(pair_idx2, 1));
		Vector2d coeff_tz(coeff_t(pair_idx1, 2), coeff_t(pair_idx2, 2));


		vector<double> coeff_l1_norm(3);
		coeff_l1_norm[0] = coeff_tx.lpNorm<1>();
		coeff_l1_norm[1] = coeff_ty.lpNorm<1>();
		coeff_l1_norm[2] = coeff_tz.lpNorm<1>();
		int base_indicator = min_element(coeff_l1_norm.begin(),
			coeff_l1_norm.end()) - coeff_l1_norm.begin();


		Vector3d model(0.0, 0.0, 0.0);
		if (base_indicator == 0)
		{
			Matrix2d A;
			A << coeff_ty, coeff_tz;
			Vector2d solution = A.inverse() * (-coeff_tx);
			model(0) = 1.0;
			model(1) = solution(0);
			model(2) = solution(1);
		}
		else if (base_indicator == 1)
		{
			Matrix2d A;
			A << coeff_tx, coeff_tz;
			Vector2d solution = A.inverse() * (-coeff_ty);
			model(0) = solution(0);
			model(1) = 1.0;
			model(2) = solution(1);
		}
		else {
			Matrix2d A;
			A << coeff_tx, coeff_ty;
			Vector2d solution = A.inverse() * (-coeff_tz);
			model(0) = solution(0);
			model(1) = solution(1);
			model(2) = 1.0;
		}


		VectorXd error = coeff_t * model;

		vector<int> inlier_set;
		for (int i = 0; i < error.rows(); ++i)
		{
			if (inlier_markers[i] == 0)
				continue;
			if (std::abs(error(i)) < inlier_error*norm_pixel_unit)
				inlier_set.push_back(i);
		}


		if (inlier_set.size() < 0.2*pts1_undistorted.size())
			continue;


		VectorXd coeff_tx_better(inlier_set.size());
		VectorXd coeff_ty_better(inlier_set.size());
		VectorXd coeff_tz_better(inlier_set.size());
		for (int i = 0; i < inlier_set.size(); ++i) {
			coeff_tx_better(i) = coeff_t(inlier_set[i], 0);
			coeff_ty_better(i) = coeff_t(inlier_set[i], 1);
			coeff_tz_better(i) = coeff_t(inlier_set[i], 2);
		}

		Vector3d model_better(0.0, 0.0, 0.0);
		if (base_indicator == 0)
		{
			MatrixXd A(inlier_set.size(), 2);
			A << coeff_ty_better, coeff_tz_better;
			Vector2d solution =
				(A.transpose() * A).inverse() * A.transpose() * (-coeff_tx_better);
			model_better(0) = 1.0;
			model_better(1) = solution(0);
			model_better(2) = solution(1);
		}
		else if (base_indicator == 1)
		{
			MatrixXd A(inlier_set.size(), 2);
			A << coeff_tx_better, coeff_tz_better;
			Vector2d solution =
				(A.transpose() * A).inverse() * A.transpose() * (-coeff_ty_better);
			model_better(0) = solution(0);
			model_better(1) = 1.0;
			model_better(2) = solution(1);
		}
		else
		{
			MatrixXd A(inlier_set.size(), 2);
			A << coeff_tx_better, coeff_ty_better;
			Vector2d solution =
				(A.transpose() * A).inverse() * A.transpose() * (-coeff_tz_better);
			model_better(0) = solution(0);
			model_better(1) = solution(1);
			model_better(2) = 1.0;
		}

	
		VectorXd new_error = coeff_t * model_better;

		double this_error = 0.0;
		for (const auto& inlier_idx : inlier_set)
			this_error += std::abs(new_error(inlier_idx));
		this_error /= inlier_set.size();


		if (inlier_set.size() > best_inlier_set.size())
		{
			best_error = this_error;
			best_inlier_set = inlier_set;
		}
	}


	inlier_markers.clear();
	inlier_markers.resize(pts1.size(), 0);
	for (const auto& inlier_idx : best_inlier_set)
		inlier_markers[inlier_idx] = 1;



	return;
}




void gvis::t_gstereoorblk::rescalePoints(
	vector<cv::Point2f>& pts1, vector<cv::Point2f>& pts2,
	float& scaling_factor)
{
	scaling_factor = 0.0f;

	
	for (int i = 0; i < pts1.size(); ++i)
	{
		scaling_factor += sqrt(pts1[i].dot(pts1[i]));
		scaling_factor += sqrt(pts2[i].dot(pts2[i]));
	}


	scaling_factor = (pts1.size() + pts2.size()) /
		scaling_factor * sqrt(2.0f);


	for (int i = 0; i < pts1.size(); ++i)
	{
		pts1[i] *= scaling_factor;
		pts2[i] *= scaling_factor;
	}

	return;
}



void gvis::t_gstereoorblk::ExtractORB(const cv::Mat &im, int flag)
{
	if (flag == 0)
		(*mpORBextractorLeft)(im, cv::Mat(), mvKeysLeft, mDescriptorsLeft);
	else
		(*mpORBextractorRight)(im, cv::Mat(), mvKeysRight, mDescriptorsRight);
}



bool gvis::t_gstereoorblk::featureCompareByTrack(const t_gframe::PointsAttri & f1, const t_gframe::PointsAttri & f2)
{
	return f1.track_cnt > f2.track_cnt;
}


