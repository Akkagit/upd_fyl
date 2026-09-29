#ifdef USE_VIS
#include "gstereolk.h"
#include "gset/gsetvis.h"
#include "opencv2/core/eigen.hpp"
using namespace Eigen;
using namespace gvis;
using namespace cv;
//using namespace gnut;

gvis::t_gstereolk::t_gstereolk(gnut::t_gsetbase* _set) :t_gimgproc(_set),
prev_features_ptr(new GridFeatures()),
curr_features_ptr(new GridFeatures()),
isFirstImg(true),
IsInitialized(false)
{
	next_feature_id = 0;

}

PointCloud gvis::t_gstereolk::ProcessBatch(const double& t)
{
	cout << "==================================" << endl;
	cout << "start process image" << endl;
	

	cv::Mat _img0 = cv::imread(cur_img_path.img0_path, 0);
	cv::Mat _img1 = cv::imread(cur_img_path.img1_path, 0);
	dataIMG _img_msg;
	_img_msg.t = cur_img_path.t;
	_img_msg.img0 = make_shared<cv::Mat>(_img0);

	_img_msg.img1 = make_shared<cv::Mat>(_img1);
	cam0_curr_img = *(_img_msg.img0);
	cam1_curr_img = *(_img_msg.img1);
	//cout << cam0_curr_img.cols << ", " << cam1_curr_img.cols << endl;
	curr_img_msg = _img_msg;

	// 图像均衡化
	// if (EQUALIZE)
	// {
	// cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE(3.0, cv::Size(8, 8));

	// clahe->apply(cam0_curr_img, cam0_curr_img);
	// clahe->apply(cam1_curr_img, cam1_curr_img);
	// }


	if (!IsInitialized)
	{
		Initialize();
		IsInitialized = true;
	}

	createImagePyramids();
	
	if (isFirstImg)
	{

		initializeFirstFrame();

		std::cout << "detect first image" << std::endl;
		isFirstImg = false;
	
		//只保留最后一个imu
		imumsg last_imu = _vecimu.at(_vecimu.size() - 1);
		_vecimu.clear();
		_vecimu.push_back(last_imu);
		
	
	}
	else
	{

		trackFeatures();

		addNewFeatures();

		pruneGridFeatures();

	}



	drawFeaturesStereo();
	//cout << "img_t:" << setprecision(8)<<curr_img_msg.t << endl;
	publish();
	//cin.get();
	prev_img_msg = curr_img_msg;
	prev_features_ptr = curr_features_ptr;
	std::swap(prev_cam0_pyramid, curr_cam0_pyramid);



	curr_features_ptr.reset(new GridFeatures());
	for (int code = 0; code < grid_row * grid_col; ++code)
	{
		(*curr_features_ptr)[code] = vector<FeaturePoint>(0);
	}
	
	return _pointCloud;
}

bool gvis::t_gstereolk::Initialize()
{
	// Create feature detector（fast）
	detector_ptr = cv::FastFeatureDetector::create(fast_threshold);

	eigen2cv(R_cam0_cam1, _R_cam0_cam1);
	eigen2cv(t_cam0_cam1, _t_cam0_cam1);
	eigen2cv(R_cam0_imu, _R_cam0_imu);
	eigen2cv(t_cam0_imu, _t_cam0_imu);
	_R_cam1_imu = _R_cam0_imu * _R_cam0_cam1.t();
	_cam0_distortion_model = distortion2str(cam0_distortion_model);
	_cam0_intrinsics << cam0_intrinsics(0), cam0_intrinsics(1), cam0_intrinsics(2), cam0_intrinsics(3);
	_cam0_distortion_coeffs << cam0_distortion_coeffs(0), cam0_distortion_coeffs(1), cam0_distortion_coeffs(2), cam0_distortion_coeffs(3);

	_cam1_distortion_model = distortion2str(cam1_distortion_model);
	_cam1_intrinsics << cam1_intrinsics(0), cam1_intrinsics(1), cam1_intrinsics(2), cam1_intrinsics(3);
	_cam1_distortion_coeffs << cam1_distortion_coeffs(0), cam1_distortion_coeffs(1), cam1_distortion_coeffs(2), cam1_distortion_coeffs(3);
	return true;
}

void gvis::t_gstereolk::createImagePyramids()
{
	buildOpticalFlowPyramid(
		cam0_curr_img, curr_cam0_pyramid,
		cv::Size(patch_size, patch_size),
		pyramid_levels, true, cv::BORDER_REFLECT_101,
		cv::BORDER_CONSTANT, false);

	buildOpticalFlowPyramid(
		cam1_curr_img, curr_cam1_pyramid,
		cv::Size(patch_size, patch_size),
		pyramid_levels, true, cv::BORDER_REFLECT_101,
		cv::BORDER_CONSTANT, false);
}

void gvis::t_gstereolk::initializeFirstFrame()
{

	const Mat& img = cam0_curr_img;


	static int grid_height = img.rows / grid_row;
	static int grid_width = img.cols / grid_col;


	
	vector<KeyPoint> new_features(0);
	detector_ptr->detect(img, new_features);

	

	vector<cv::Point2f> cam0_points(new_features.size());
	for (int i = 0; i < new_features.size(); ++i)
		cam0_points[i] = new_features[i].pt;

	vector<cv::Point2f> cam1_points(0);
	vector<unsigned char> inlier_markers(0);
	stereoMatch(cam0_points, cam1_points, inlier_markers);


	vector<cv::Point2f> cam0_inliers(0);
	vector<cv::Point2f> cam1_inliers(0);
	vector<float> response_inliers(0);
	for (int i = 0; i < inlier_markers.size(); ++i)
	{
		if (inlier_markers[i] == 0)
			continue;
		cam0_inliers.push_back(cam0_points[i]);
		cam1_inliers.push_back(cam1_points[i]);


		response_inliers.push_back(new_features[i].response);
	}
	
	GridFeatures grid_new_features;
	for (int code = 0; code < grid_row * grid_col; ++code)
		grid_new_features[code] = vector<FeaturePoint>(0);

	for (int i = 0; i < cam0_inliers.size(); ++i)
	{
		const cv::Point2f& cam0_point = cam0_inliers[i];
		const cv::Point2f& cam1_point = cam1_inliers[i];
		const float& response = response_inliers[i];

		int row = static_cast<int>(cam0_point.y / grid_height);
		int col = static_cast<int>(cam0_point.x / grid_width);
		int code = row * grid_col + col;


		FeaturePoint new_feature;
		new_feature.response = response;
		new_feature.cam0_point = cam0_point;
		new_feature.cam1_point = cam1_point;
		grid_new_features[code].push_back(new_feature);
	}


	for (auto& item : grid_new_features)
		std::sort(item.second.begin(), item.second.end(),
			&t_gstereolk::featureCompareByResponse);


	for (int code = 0; code < grid_row * grid_col; ++code)
	{

		vector<FeaturePoint>& features_this_grid = (*curr_features_ptr)[code];
		vector<FeaturePoint>& new_features_this_grid = grid_new_features[code];


		for (int k = 0; k < grid_min_feature_num &&
			k < new_features_this_grid.size(); ++k)
		{
			features_this_grid.push_back(new_features_this_grid[k]);
			features_this_grid.back().id = next_feature_id++;
			features_this_grid.back().lifetime = 1;
		}
	}
	
	return;

}

void gvis::t_gstereolk::stereoMatch(const vector<cv::Point2f>& cam0_points,
	vector<cv::Point2f>& cam1_points,
	vector<unsigned char>& inlier_markers)
{
	if (cam0_points.size() == 0)
		return;


	if (cam1_points.size() == 0)
	{
		
		vector<cv::Point2f> cam0_points_undistorted;
		string _cam0_distortion_model = distortion2str(cam0_distortion_model);

		undistortPoints(cam0_points, _cam0_intrinsics, _cam0_distortion_model,
			_cam0_distortion_coeffs, cam0_points_undistorted,
			_R_cam0_cam1);


		cam1_points = distortPoints(cam0_points_undistorted, _cam1_intrinsics,
			_cam1_distortion_model, _cam1_distortion_coeffs);
	}




	calcOpticalFlowPyrLK(curr_cam0_pyramid, curr_cam1_pyramid,
		cam0_points, cam1_points,
		inlier_markers, noArray(),
		Size(patch_size, patch_size),
		pyramid_levels,
		TermCriteria(TermCriteria::COUNT + TermCriteria::EPS,
			max_iteration, track_precision),
		cv::OPTFLOW_USE_INITIAL_FLOW);



	for (int i = 0; i < cam1_points.size(); ++i)
	{
		if (inlier_markers[i] == 0)
			continue;
		if (cam1_points[i].y < 0 ||
			cam1_points[i].y > cam1_curr_img.rows - 1 ||
			cam1_points[i].x < 0 ||
			cam1_points[i].x > cam1_curr_img.cols - 1)
		{
			inlier_markers[i] = 0;
		}

	}

	const cv::Matx33d t_cam0_cam1_hat(
		0.0, -t_cam0_cam1[2], t_cam0_cam1[1],
		t_cam0_cam1[2], 0.0, -t_cam0_cam1[0],
		-t_cam0_cam1[1], t_cam0_cam1[0], 0.0);
	const cv::Matx33d E = t_cam0_cam1_hat * _R_cam0_cam1;



	vector<cv::Point2f> cam0_points_undistorted(0);
	vector<cv::Point2f> cam1_points_undistorted(0);

	undistortPoints(
		cam0_points, _cam0_intrinsics, _cam0_distortion_model,
		_cam0_distortion_coeffs, cam0_points_undistorted);
	undistortPoints(
		cam1_points, _cam1_intrinsics, _cam1_distortion_model,
		_cam0_distortion_coeffs, cam1_points_undistorted);



	double norm_pixel_unit = 4.0 / (
		_cam0_intrinsics[0] + _cam0_intrinsics[1] +
		_cam1_intrinsics[0] + _cam1_intrinsics[1]);


	for (int i = 0; i < cam0_points_undistorted.size(); ++i)
	{
		if (inlier_markers[i] == 0)
			continue;


		cv::Vec3d pt0(cam0_points_undistorted[i].x,
			cam0_points_undistorted[i].y, 1.0);
		cv::Vec3d pt1(cam1_points_undistorted[i].x,
			cam1_points_undistorted[i].y, 1.0);


		cv::Vec3d epipolar_line = E * pt0;


		double error = fabs((pt1.t() * epipolar_line)[0]) / sqrt(
			epipolar_line[0] * epipolar_line[0] +
			epipolar_line[1] * epipolar_line[1]);


		if (error > stereo_threshold * norm_pixel_unit)
			inlier_markers[i] = 0;
	}

	return;
}

void gvis::t_gstereolk::drawFeaturesStereo()
{

	Scalar tracked(0, 255, 0);
	Scalar new_feature(0, 255, 255);


	static int grid_height =
		cam0_curr_img.rows / grid_row;
	static int grid_width =
		cam0_curr_img.cols / grid_col;


	int img_height = cam0_curr_img.rows;
	int img_width = cam0_curr_img.cols;
	Mat out_img(img_height, img_width * 2, CV_8UC3);

	cvtColor(cam0_curr_img, out_img.colRange(0, img_width), cv::COLOR_GRAY2BGR);
	cvtColor(cam1_curr_img, out_img.colRange(img_width, img_width * 2), cv::COLOR_GRAY2BGR);



	for (int i = 1; i < grid_row; ++i)
	{
		Point pt1(0, i * grid_height);
		Point pt2(img_width * 2, i * grid_height);
		line(out_img, pt1, pt2, Scalar(255, 0, 0));
	}
	for (int i = 1; i < grid_col; ++i)
	{
		Point pt1(i * grid_width, 0);
		Point pt2(i * grid_width, img_height);
		line(out_img, pt1, pt2, Scalar(255, 0, 0));
	}
	for (int i = 1; i < grid_col; ++i)
	{
		Point pt1(i * grid_width + img_width, 0);
		Point pt2(i * grid_width + img_width, img_height);
		line(out_img, pt1, pt2, Scalar(255, 0, 0));
	}





	vector<long long int> prev_ids(0);
	for (const auto &grid_features : *prev_features_ptr)
		for (const auto &feature : grid_features.second)
			prev_ids.push_back(feature.id);



	map<long long int, Point2f> prev_cam0_points;
	map<long long int, Point2f> prev_cam1_points;
	for (const auto &grid_features : *prev_features_ptr)
		for (const auto &feature : grid_features.second)
		{
			prev_cam0_points[feature.id] = feature.cam0_point;
			prev_cam1_points[feature.id] = feature.cam1_point;
		}


	map<long long int, Point2f> curr_cam0_points;
	map<long long int, Point2f> curr_cam1_points;
	for (const auto &grid_features : *curr_features_ptr)
		for (const auto &feature : grid_features.second)
		{
			curr_cam0_points[feature.id] = feature.cam0_point;
			curr_cam1_points[feature.id] = feature.cam1_point;
		}



	for (const auto &id : prev_ids)
	{

		if (prev_cam0_points.find(id) != prev_cam0_points.end() &&
			curr_cam0_points.find(id) != curr_cam0_points.end())
		{
			cv::Point2f prev_pt0 = prev_cam0_points[id];
			cv::Point2f prev_pt1 = prev_cam1_points[id] + Point2f(img_width, 0.0);
			cv::Point2f curr_pt0 = curr_cam0_points[id];
			cv::Point2f curr_pt1 = curr_cam1_points[id] + Point2f(img_width, 0.0);

			circle(out_img, curr_pt0, 3, tracked, -1);
			circle(out_img, curr_pt1, 3, tracked, -1);
			line(out_img, prev_pt0, curr_pt0, tracked, 1);
			line(out_img, prev_pt1, curr_pt1, tracked, 1);

			prev_cam0_points.erase(id);
			prev_cam1_points.erase(id);
			curr_cam0_points.erase(id);
			curr_cam1_points.erase(id);
		}
	}


	for (const auto &new_cam0_point : curr_cam0_points)
	{
		cv::Point2f pt0 = new_cam0_point.second;
		cv::Point2f pt1 = curr_cam1_points[new_cam0_point.first] +
			Point2f(img_width, 0.0);

		circle(out_img, pt0, 3, new_feature, -1);
		circle(out_img, pt1, 3, new_feature, -1);
	}

	cv::resize(out_img, out_img, out_img.size());
	//imwrite("stereo.jpg", out_img);
    cv::resize(out_img, out_img, cv::Size(out_img.cols / 2.0, out_img.rows / 2.0));
	imshow("Feature", out_img);
	waitKey(5);
	return;
}

void gvis::t_gstereolk::publish()
{
	//sensor_msgs::FeaturePoint feature_msg;
	_pointCloud.features.clear();
	FeaturePoint feature_msg;
	//pointCloud = sensor_msgs::PointCloud();
	_pointCloud.time = curr_img_msg.t;

	vector<long long int> curr_ids(0);
	vector<Point2f> curr_cam0_points(0);
	vector<Point2f> curr_cam1_points(0);


	for (const auto& grid_features : (*curr_features_ptr))
	{
		for (const auto& feature : grid_features.second)
		{
			curr_ids.push_back(feature.id);
			curr_cam0_points.push_back(feature.cam0_point);
			curr_cam1_points.push_back(feature.cam1_point);
		}
	}
	vector<Point2f> curr_cam0_points_undistorted(0);
	vector<Point2f> curr_cam1_points_undistorted(0);

	undistortPoints(
		curr_cam0_points, _cam0_intrinsics, _cam0_distortion_model,
		_cam0_distortion_coeffs, curr_cam0_points_undistorted);
	undistortPoints(
		curr_cam1_points, _cam1_intrinsics, _cam1_distortion_model,
		_cam1_distortion_coeffs, curr_cam1_points_undistorted);


	for (int i = 0; i < curr_ids.size(); ++i)
	{
		//cout << "imgproc_id" << feature_msg.id << endl;
		feature_msg.id = curr_ids[i];
		feature_msg.cam0_point = curr_cam0_points_undistorted[i];
		
		feature_msg.cam1_point = curr_cam1_points_undistorted[i];

		_pointCloud.features.push_back(feature_msg);
	}
	return;
}

void gvis::t_gstereolk::trackFeatures()
{

	static int grid_height =
		cam0_curr_img.rows / grid_row;
	static int grid_width =
		cam0_curr_img.cols / grid_col;



	Matx33f cam0_R_p_c;
	Matx33f cam1_R_p_c;
	integrateImuData(cam0_R_p_c, cam1_R_p_c);



	vector<long long int> prev_ids(0);
	vector<int> prev_lifetime(0);
	vector<Point2f> prev_cam0_points(0);
	vector<Point2f> prev_cam1_points(0);

	for (const auto& item : *prev_features_ptr)
	{
		for (const auto& prev_feature : item.second)
		{
			prev_ids.push_back(prev_feature.id);
			prev_lifetime.push_back(prev_feature.lifetime);
			prev_cam0_points.push_back(prev_feature.cam0_point);
			//cout << prev_feature.cam0_point.x <<" ,"<<prev_feature.cam0_point.y << endl;
			prev_cam1_points.push_back(prev_feature.cam1_point);
		}
	}


	before_tracking = prev_cam0_points.size();


	if (prev_ids.size() == 0)
		return;



	vector<Point2f> curr_cam0_points(0);
	vector<unsigned char> track_inliers(0);


	predictFeatureTracking(prev_cam0_points,
		cam0_R_p_c, _cam0_intrinsics, curr_cam0_points);


	calcOpticalFlowPyrLK(
		prev_cam0_pyramid, curr_cam0_pyramid,
		prev_cam0_points, curr_cam0_points,
		track_inliers, noArray(),
		Size(patch_size, patch_size),
		pyramid_levels,
		TermCriteria(TermCriteria::COUNT + TermCriteria::EPS,
			max_iteration,
			track_precision),
		cv::OPTFLOW_USE_INITIAL_FLOW);



	for (int i = 0; i < curr_cam0_points.size(); ++i) {
		if (track_inliers[i] == 0) continue;
		if (curr_cam0_points[i].y < 0 ||
			curr_cam0_points[i].y > cam0_curr_img.rows - 1 ||
			curr_cam0_points[i].x < 0 ||
			curr_cam0_points[i].x > cam0_curr_img.cols - 1)
			track_inliers[i] = 0;
	}


	vector<long long int> prev_tracked_ids(0);
	vector<int> prev_tracked_lifetime(0);
	vector<Point2f> prev_tracked_cam0_points(0);
	vector<Point2f> prev_tracked_cam1_points(0);
	vector<Point2f> curr_tracked_cam0_points(0);


	removeUnmarkedElements(
		prev_ids, track_inliers, prev_tracked_ids);
	removeUnmarkedElements(
		prev_lifetime, track_inliers, prev_tracked_lifetime);
	removeUnmarkedElements(
		prev_cam0_points, track_inliers, prev_tracked_cam0_points);
	removeUnmarkedElements(
		prev_cam1_points, track_inliers, prev_tracked_cam1_points);
	removeUnmarkedElements(
		curr_cam0_points, track_inliers, curr_tracked_cam0_points);


	after_tracking = curr_tracked_cam0_points.size();



	// Outlier removal involves three steps, which forms a close
	// loop between the previous and current frames of cam0 (left)
	// and cam1 (right). Assuming the stereo matching between the
	// previous cam0 and cam1 images are correct, the three steps are:
	//
	// prev frames cam0 ----------> cam1
	//              |                |
	//              |ransac          |ransac
	//              |   stereo match |
	// curr frames cam0 ----------> cam1
	//
	// 1) Stereo matching between current images of cam0 and cam1.
	// 2) RANSAC between previous and current images of cam0.
	// 3) RANSAC between previous and current images of cam1.
	//
	// For Step 3, tracking between the images is no longer needed.
	// The stereo matching results are directly used in the RANSAC.


	vector<Point2f> curr_cam1_points(0);
	vector<unsigned char> match_inliers(0);
	stereoMatch(curr_tracked_cam0_points, curr_cam1_points, match_inliers);

	vector<long long int > prev_matched_ids(0);
	vector<int> prev_matched_lifetime(0);
	vector<Point2f> prev_matched_cam0_points(0);
	vector<Point2f> prev_matched_cam1_points(0);
	vector<Point2f> curr_matched_cam0_points(0);
	vector<Point2f> curr_matched_cam1_points(0);


	removeUnmarkedElements(
		prev_tracked_ids, match_inliers, prev_matched_ids);
	removeUnmarkedElements(
		prev_tracked_lifetime, match_inliers, prev_matched_lifetime);
	removeUnmarkedElements(
		prev_tracked_cam0_points, match_inliers, prev_matched_cam0_points);
	removeUnmarkedElements(
		prev_tracked_cam1_points, match_inliers, prev_matched_cam1_points);
	removeUnmarkedElements(
		curr_tracked_cam0_points, match_inliers, curr_matched_cam0_points);
	removeUnmarkedElements(
		curr_cam1_points, match_inliers, curr_matched_cam1_points);


	after_matching = curr_matched_cam0_points.size();




	vector<int> cam0_ransac_inliers(0);
	twoPointRansac(prev_matched_cam0_points, curr_matched_cam0_points,
		cam0_R_p_c, _cam0_intrinsics, _cam0_distortion_model,
		_cam0_distortion_coeffs, ransac_threshold,
		0.99, cam0_ransac_inliers);

	vector<int> cam1_ransac_inliers(0);
	twoPointRansac(prev_matched_cam1_points, curr_matched_cam1_points,
		cam1_R_p_c, _cam1_intrinsics, _cam1_distortion_model,
		_cam1_distortion_coeffs, ransac_threshold,
		0.99, cam1_ransac_inliers);


	after_ransac = 0;

	for (int i = 0; i < cam0_ransac_inliers.size(); ++i)
	{

		if (cam0_ransac_inliers[i] == 0 || cam1_ransac_inliers[i] == 0)
			continue;

		int row = static_cast<int>(curr_matched_cam0_points[i].y / grid_height);
		int col = static_cast<int>(curr_matched_cam0_points[i].x / grid_width);

		int code = row * grid_col + col;


		(*curr_features_ptr)[code].push_back(FeaturePoint());
		FeaturePoint& grid_new_feature = (*curr_features_ptr)[code].back();

		grid_new_feature.id = prev_matched_ids[i];
		grid_new_feature.lifetime = ++prev_matched_lifetime[i];
		grid_new_feature.cam0_point = curr_matched_cam0_points[i];
		grid_new_feature.cam1_point = curr_matched_cam1_points[i];

		++after_ransac;
	}


	int prev_feature_num = 0;
	for (const auto& item : *prev_features_ptr)
		prev_feature_num += item.second.size();

	int curr_feature_num = 0;
	for (const auto& item : *curr_features_ptr)
		curr_feature_num += item.second.size();

	printf("candidates: %d; track: %d; match: %d; ransac: %d/%d=%f; after_ransac：%d \n",
		before_tracking, after_tracking, after_matching,
		curr_feature_num, prev_feature_num,
		static_cast<double>(curr_feature_num) /
		(static_cast<double>(prev_feature_num) + 1e-5), after_ransac);
	cout << "==================================" << endl;

	return;
}



void gvis::t_gstereolk::integrateImuData(
	cv::Matx33f& cam0_R_p_c, cv::Matx33f& cam1_R_p_c)
{

	auto begin_iter = _vecimu.begin();
	while (begin_iter != _vecimu.end())
	{
		if (begin_iter->t -
			prev_img_msg.t < -0.01)
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
			curr_img_msg.t < 0.002)
			++end_iter;
		else
			break;
	}
	//cout << "imu_size:" << _vecimu.size() << endl;
	//cout << "prev_img:" <<setprecision(10)<< prev_img_msg.t << endl;
	//cout << "cur_img:" <<setprecision(10)<< curr_img_msg.t << endl;
	//for (auto imu : _vecimu)
	//{
	////	cout << "imu_t:" << imu.t << endl;
	//}
	//cout << "use_imu_size:" << end_iter - begin_iter << endl;

	Vec3f mean_ang_vel(0.0, 0.0, 0.0);
	for (auto iter = begin_iter; iter < end_iter; ++iter)

	{
		mean_ang_vel += Vec3f(iter->angular_velocity.x(),
			iter->angular_velocity.y(), iter->angular_velocity.z());
		
	}

	if (end_iter - begin_iter > 0)
		mean_ang_vel *= 1.0f / (end_iter - begin_iter);


	Vec3f cam0_mean_ang_vel = _R_cam0_imu.t() * mean_ang_vel;
	Vec3f cam1_mean_ang_vel = _R_cam1_imu.t() * mean_ang_vel;


	double dtime = (curr_img_msg.t -
		prev_img_msg.t);
	Rodrigues(cam0_mean_ang_vel*dtime, cam0_R_p_c);
	Rodrigues(cam1_mean_ang_vel*dtime, cam1_R_p_c);
	cam0_R_p_c = cam0_R_p_c.t();
	cam1_R_p_c = cam1_R_p_c.t();


	_vecimu.erase(_vecimu.begin(), end_iter);
	return;
	
}



void gvis::t_gstereolk::predictFeatureTracking(
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





void gvis::t_gstereolk::twoPointRansac(
	const vector<cv::Point2f>& pts1, const vector<cv::Point2f>& pts2,
	const cv::Matx33f& R_p_c, const cv::Vec4d& intrinsics,
	const std::string& distortion_model,
	const cv::Vec4d& distortion_coeffs,
	const double& inlier_error,
	const double& success_probability,
	vector<int>& inlier_markers)
{


	if (pts1.size() != pts2.size())
		printf("Sets of different size (%d and %d) are used...\n",
			pts1.size(), pts2.size());


	double norm_pixel_unit = 2.0 / (intrinsics[0] + intrinsics[1]);


	int iter_num = static_cast<int>(
		ceil(log(1 - success_probability) / log(1 - 0.7*0.7)));


	inlier_markers.clear();
	inlier_markers.resize(pts1.size(), 1);


	vector<Point2f> pts1_undistorted(pts1.size());
	vector<Point2f> pts2_undistorted(pts2.size());
	undistortPoints(
		pts1, intrinsics, distortion_model,
		distortion_coeffs, pts1_undistorted);
	undistortPoints(
		pts2, intrinsics, distortion_model,
		distortion_coeffs, pts2_undistorted);



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

		// Construct the model;
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

		// Find all the inliers among point pairs.
		VectorXd error = coeff_t * model;

		vector<int> inlier_set;
		for (int i = 0; i < error.rows(); ++i)
		{
			if (inlier_markers[i] == 0)
				continue;
			if (std::abs(error(i)) < inlier_error*norm_pixel_unit)
				inlier_set.push_back(i);
		}

		// If the number of inliers is small, the current
		// model is probably wrong.
		if (inlier_set.size() < 0.2*pts1_undistorted.size())
			continue;

		// Refit the model using all of the possible inliers.
		VectorXd coeff_tx_better(inlier_set.size());
		VectorXd coeff_ty_better(inlier_set.size());
		VectorXd coeff_tz_better(inlier_set.size());
		for (int i = 0; i < inlier_set.size(); ++i) {
			coeff_tx_better(i) = coeff_t(inlier_set[i], 0);
			coeff_ty_better(i) = coeff_t(inlier_set[i], 1);
			coeff_tz_better(i) = coeff_t(inlier_set[i], 2);
		}

		Vector3d model_better(0.0, 0.0, 0.0);
		if (base_indicator == 0) {
			MatrixXd A(inlier_set.size(), 2);
			A << coeff_ty_better, coeff_tz_better;
			Vector2d solution =
				(A.transpose() * A).inverse() * A.transpose() * (-coeff_tx_better);
			model_better(0) = 1.0;
			model_better(1) = solution(0);
			model_better(2) = solution(1);
		}
		else if (base_indicator == 1) {
			MatrixXd A(inlier_set.size(), 2);
			A << coeff_tx_better, coeff_tz_better;
			Vector2d solution =
				(A.transpose() * A).inverse() * A.transpose() * (-coeff_ty_better);
			model_better(0) = solution(0);
			model_better(1) = 1.0;
			model_better(2) = solution(1);
		}
		else {
			MatrixXd A(inlier_set.size(), 2);
			A << coeff_tx_better, coeff_ty_better;
			Vector2d solution =
				(A.transpose() * A).inverse() * A.transpose() * (-coeff_tz_better);
			model_better(0) = solution(0);
			model_better(1) = solution(1);
			model_better(2) = 1.0;
		}

		// Compute the error and upate the best model if possible.
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

	// Fill in the markers.
	inlier_markers.clear();
	inlier_markers.resize(pts1.size(), 0);
	for (const auto& inlier_idx : best_inlier_set)
		inlier_markers[inlier_idx] = 1;

	//printf("inlier ratio: %lu/%lu\n",
	//    best_inlier_set.size(), inlier_markers.size());

	return;
}



void gvis::t_gstereolk::addNewFeatures()
{
	const Mat& curr_img = cam0_curr_img;

	// Size of each grid.
	static int grid_height =
		cam0_curr_img.rows / grid_row;
	static int grid_width =
		cam0_curr_img.cols / grid_col;



	Mat mask(curr_img.rows, curr_img.cols, CV_8U, Scalar(1));

	for (const auto& features : *curr_features_ptr)
	{
		for (const auto& feature : features.second)
		{
			const int y = static_cast<int>(feature.cam0_point.y);
			const int x = static_cast<int>(feature.cam0_point.x);

			int up_lim = y - 2, bottom_lim = y + 3,
				left_lim = x - 2, right_lim = x + 3;
			if (up_lim < 0)
				up_lim = 0;
			if (bottom_lim > curr_img.rows)
				bottom_lim = curr_img.rows;
			if (left_lim < 0)
				left_lim = 0;
			if (right_lim > curr_img.cols)
				right_lim = curr_img.cols;

			Range row_range(up_lim, bottom_lim);
			Range col_range(left_lim, right_lim);
			mask(row_range, col_range) = 0;
		}
	}



	vector<KeyPoint> new_features(0);
	detector_ptr->detect(curr_img, new_features, mask);



	vector<vector<KeyPoint> > new_feature_sieve(grid_row * grid_col);
	for (const auto& feature : new_features)
	{
		int row = static_cast<int>(feature.pt.y / grid_height); row = (row < grid_row ? row: grid_row - 1);
		int col = static_cast<int>(feature.pt.x / grid_width); col = (col < grid_col ? col : grid_col - 1);
		new_feature_sieve[row*grid_col + col].push_back(feature);
	}
	new_features.clear();



	for (auto& item : new_feature_sieve)
	{

		if (item.size() > grid_max_feature_num)
		{
			std::sort(item.begin(), item.end(),
				&t_gstereolk::keyPointCompareByResponse);
			item.erase(
				item.begin() + grid_max_feature_num, item.end());
		}
		new_features.insert(new_features.end(), item.begin(), item.end());
	}
	int detected_new_features = new_features.size();



	vector<cv::Point2f> cam0_points(new_features.size());
	for (int i = 0; i < new_features.size(); ++i)
		cam0_points[i] = new_features[i].pt;
	vector<cv::Point2f> cam1_points(0);
	vector<unsigned char> inlier_markers(0);

	stereoMatch(cam0_points, cam1_points, inlier_markers);


	vector<cv::Point2f> cam0_inliers(0);
	vector<cv::Point2f> cam1_inliers(0);
	vector<float> response_inliers(0);
	for (int i = 0; i < inlier_markers.size(); ++i)
	{
		if (inlier_markers[i] == 0) continue;
		cam0_inliers.push_back(cam0_points[i]);
		cam1_inliers.push_back(cam1_points[i]);
		response_inliers.push_back(new_features[i].response);
	}


	int matched_new_features = cam0_inliers.size();
	if (matched_new_features < 5 &&
		static_cast<double>(matched_new_features) /
		static_cast<double>(detected_new_features) < 0.1)
	{
		printf("Images at [%f] seems unsynced...\n",
			curr_img_msg.t);
	}




	GridFeatures grid_new_features;
	for (int code = 0; code < grid_row * grid_col; ++code)
		grid_new_features[code] = vector<FeaturePoint>(0);

	for (int i = 0; i < cam0_inliers.size(); ++i)
	{
		const cv::Point2f& cam0_point = cam0_inliers[i];
		const cv::Point2f& cam1_point = cam1_inliers[i];
		const float& response = response_inliers[i];

		int row = static_cast<int>(cam0_point.y / grid_height);
		int col = static_cast<int>(cam0_point.x / grid_width);
		int code = row * grid_col + col;

		FeaturePoint new_feature;
		new_feature.response = response;
		new_feature.cam0_point = cam0_point;
		new_feature.cam1_point = cam1_point;
		grid_new_features[code].push_back(new_feature);
	}


	for (auto& item : grid_new_features)
		std::sort(item.second.begin(), item.second.end(),
			&t_gstereolk::featureCompareByResponse);

	int new_added_feature_num = 0;



	for (int code = 0; code < grid_row * grid_col; ++code)
	{
		vector<FeaturePoint>& features_this_grid = (*curr_features_ptr)[code];
		vector<FeaturePoint>& new_features_this_grid = grid_new_features[code];


		if (features_this_grid.size() >= grid_min_feature_num)
			continue;


		int vacancy_num = grid_min_feature_num -
			features_this_grid.size();


		for (int k = 0; k < vacancy_num && k < new_features_this_grid.size(); ++k)
		{
			features_this_grid.push_back(new_features_this_grid[k]);
			features_this_grid.back().id = next_feature_id++;
			features_this_grid.back().lifetime = 1;

			++new_added_feature_num;
		}
		
	}
	//cout << "add_new_feature:" << new_added_feature_num << endl;
	return;
}



void gvis::t_gstereolk::pruneGridFeatures()
{
	for (auto& item : *curr_features_ptr)
	{
		auto& grid_features = item.second;


		if (grid_features.size() <= grid_max_feature_num)
			continue;

		std::sort(grid_features.begin(), grid_features.end(),
			&t_gstereolk::featureCompareByLifetime);

		grid_features.erase(grid_features.begin() +
			grid_max_feature_num, grid_features.end());
	}
	return;
}




vector<cv::Point2f> gvis::t_gstereolk::distortPoints(
	const vector<cv::Point2f>& pts_in,
	const cv::Vec4d& intrinsics,
	const string& distortion_model,
	const cv::Vec4d& distortion_coeffs)
{
	const cv::Matx33d K(intrinsics[0], 0.0, intrinsics[2],
		0.0, intrinsics[1], intrinsics[3],
		0.0, 0.0, 1.0);

	vector<cv::Point2f> pts_out;

	if (distortion_model == "radtan")
	{
		vector<cv::Point3f> homogenous_pts;
		cv::convertPointsToHomogeneous(pts_in, homogenous_pts);
		cv::projectPoints(homogenous_pts, cv::Vec3d::zeros(), cv::Vec3d::zeros(), K,
			distortion_coeffs, pts_out);
	}

	else if (distortion_model == "equidistant")
	{
		cv::fisheye::distortPoints(pts_in, pts_out, K, distortion_coeffs);
	}

	else
	{
		printf("The model %s is unrecognized, using radtan instead...",
			distortion_model.c_str());
		vector<cv::Point3f> homogenous_pts;
		cv::convertPointsToHomogeneous(pts_in, homogenous_pts);
		cv::projectPoints(homogenous_pts, cv::Vec3d::zeros(), cv::Vec3d::zeros(), K,
			distortion_coeffs, pts_out);
	}

	return pts_out;
}

void gvis::t_gstereolk::undistortPoints(
	const vector<cv::Point2f>& pts_in,
	const cv::Vec4d& intrinsics,
	const string& distortion_model,
	const cv::Vec4d& distortion_coeffs,
	vector<cv::Point2f>& pts_out,
	const cv::Matx33d &rectification_matrix,
	const cv::Vec4d &new_intrinsics)
{

	if (pts_in.size() == 0)
		return;

	const cv::Matx33d K(
		intrinsics[0], 0.0, intrinsics[2],
		0.0, intrinsics[1], intrinsics[3],
		0.0, 0.0, 1.0);

	const cv::Matx33d K_new(
		new_intrinsics[0], 0.0, new_intrinsics[2],
		0.0, new_intrinsics[1], new_intrinsics[3],
		0.0, 0.0, 1.0);


	if (distortion_model == "radtan")
	{
		cv::undistortPoints(pts_in, pts_out, K, distortion_coeffs,
			rectification_matrix, K_new);
	}
	else if (distortion_model == "equidistant")
	{
		cv::fisheye::undistortPoints(pts_in, pts_out, K, distortion_coeffs,
			rectification_matrix, K_new);
	}
	else
	{
		printf("The model %s is unrecognized, use radtan instead...",
			distortion_model.c_str());
		cv::undistortPoints(pts_in, pts_out, K, distortion_coeffs,
			rectification_matrix, K_new);
	}
	return;
}

bool gvis::t_gstereolk::featureCompareByResponse(
	const FeaturePoint& f1,
	const FeaturePoint& f2)
{

	return f1.response > f2.response;
}

template <typename T>
void gvis::t_gstereolk::removeUnmarkedElements(
	const std::vector<T>& raw_vec,
	const std::vector<unsigned char>& markers,
	std::vector<T>& refined_vec)
{
	if (raw_vec.size() != markers.size())
	{
		printf("The input size of raw_vec(%d) and markers(%d) does not match...\n",
			raw_vec.size(), markers.size());
	}
	for (int i = 0; i < markers.size(); ++i)
	{
		if (markers[i] == 0) continue;
		refined_vec.push_back(raw_vec[i]);
	}
	return;
}

void gvis::t_gstereolk::rescalePoints(
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

bool gvis::t_gstereolk::keyPointCompareByResponse(
	const cv::KeyPoint& pt1,
	const cv::KeyPoint& pt2)
{
	return pt1.response > pt2.response;
}





bool gvis::t_gstereolk::featureCompareByLifetime(
	const FeaturePoint& f1,
	const FeaturePoint& f2)
{
	return f1.lifetime > f2.lifetime;
}

//double gvis::t_gstereolk::load_simu_img(const double &t)
//{
//
//}

#endif // USE_VIS