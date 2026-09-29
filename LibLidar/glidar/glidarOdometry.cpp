
#include <pcl/filters/filter.h>
#include <pcl/io/pcd_io.h>
#include <iostream>
#include <pcl/visualization/cloud_viewer.h>


#include "glidarOdometry.h"


namespace glidar
{

	using std::sin;
	using std::cos;
	using std::asin;
	using std::atan2;
	using std::sqrt;
	using std::fabs;
	using std::pow;

	t_glidarOdometry::t_glidarOdometry()
	{
		
	}

	void t_glidarOdometry::reset()
	{
		correspondCornerFeature_.currCornerPointCloud.clear();
		correspondCornerFeature_.correspondCornerPointCloudA.clear();
		correspondCornerFeature_.correspondCornerPointCloudB.clear();
		correspondCornerFeature_.currCornerPointCloud_inlast.clear();
		correspondSurfFeature_.currSurfPointCloud.clear();
		correspondSurfFeature_.correspondSurfPointCloudA.clear();
		correspondSurfFeature_.correspondSurfPointCloudB.clear();
		correspondSurfFeature_.correspondSurfPointCloudC.clear();
		correspondSurfFeature_.currSurfPointCloud_inlast.clear();
		correspondSurfFeature_.norm.clear();
		lidarOdoObs = lidarOdometryObs();

		last_curr_rot = Eigen::Matrix3d::Identity();
		last_curr_trans = Eigen::Vector3d::Zero();
	}

	void t_glidarOdometry::transformToStart(const pcl::PointXYZI pi, pcl::PointXYZI & po,bool distortion)
	{
		double s=1.0;
		if (distortion)
		    s = (pi.intensity - int(pi.intensity)) / scanPeriod_;
		else
		    s = 1.0;
		Eigen::Quaterniond curr_last_q(curr_last_rot);
		Eigen::Quaterniond q_point_last = Eigen::Quaterniond::Identity().slerp(s, curr_last_q);
		Eigen::Matrix3d R_point_last = q_point_last.toRotationMatrix();
		Eigen::Vector3d t_point_last = s * curr_last_trans;

		Eigen::Vector3d point(pi.x, pi.y, pi.z);
		Eigen::Vector3d un_point;
		un_point = R_point_last * point + t_point_last;

		po.x = un_point.x();
		po.y = un_point.y();
		po.z = un_point.z();
		po.intensity = pi.intensity;
	}

	// transform all lidar points to the start of the next frame
	void t_glidarOdometry::transformToEnd(const pcl::PointXYZI pi, pcl::PointXYZI & po,bool distortion)
	{
		// undistort point first
		pcl::PointXYZI un_point_tmp;
		transformToStart(pi, un_point_tmp,distortion);

		Eigen::Vector3d un_point(un_point_tmp.x, un_point_tmp.y, un_point_tmp.z);
		Eigen::Vector3d point_end = curr_last_rot.transpose()*(un_point - curr_last_trans);

		po.x = point_end.x();
		po.y = point_end.y();
		po.z = point_end.z();

		//Remove distortion time info
		po.intensity = int(pi.intensity);
	}

	void t_glidarOdometry::removeDistortion(LidarFrame& frame1, LidarFrame& frame2)
	{
		/*
		 * @brief:correct the distortion of pointcloud in frame2
		 * !!!: now doesn't correct the distortion of the "fullCloud"
		 */
		if (frame1.empty || frame2.empty)
		{
			std::cout << "ERROR:the pose of lidar frames need to be initialized!" << endl;
			getchar();
		}

		curr_last_rot = frame1.R_l_e.transpose()*frame2.R_l_e;
		curr_last_trans = frame1.R_l_e.transpose()*(frame2.t_l_e - frame1.t_l_e);
		if (1)
		{
			int cornerPointsLessSharpNum = frame2.LessSharp.points.size();
			for (int i = 0; i < cornerPointsLessSharpNum; i++)
			{
				transformToEnd(frame2.LessSharp.points[i], frame2.LessSharp.points[i],true);
			}

			int surfPointsLessFlatNum = frame2.LessSurf.points.size();
			for (int i = 0; i < surfPointsLessFlatNum; i++)
			{
				transformToEnd(frame2.LessSurf.points[i], frame2.LessSurf.points[i],true);
			}

			int surfPointsFlatNum = frame2.Surf.points.size();
			for (int i = 0; i < surfPointsFlatNum; i++)
			{
				transformToEnd(frame2.Surf.points[i], frame2.Surf.points[i], true);
			}

			int cornerPointsSharpNum = frame2.Sharp.points.size();
			for (int i = 0; i < cornerPointsSharpNum; i++)
			{
				transformToEnd(frame2.Sharp.points[i], frame2.Sharp.points[i], true);
			}
		}
	}

	void t_glidarOdometry::process(LidarFrame& frame1, LidarFrame& frame2)
	{
		if (frame1.empty || frame2.empty)
		{
			std::cout << "ERROR:the pose of lidar frames need to be initialized!" << endl;
			getchar();
		}
		///< clear the history information
		reset();
		///< set KDTree
		pcl::PointCloud<pcl::PointXYZI>::Ptr kdtreeSharp(new pcl::PointCloud<pcl::PointXYZI>());
		pcl::PointCloud<pcl::PointXYZI>::Ptr kdtreeSurf(new pcl::PointCloud<pcl::PointXYZI>());

		std::vector<int> indices;
		pcl::removeNaNFromPointCloud(frame1.LessSharp, *kdtreeSharp, indices);
		pcl::removeNaNFromPointCloud(frame1.LessSurf, *kdtreeSurf, indices);

		//cout << "frame1 size:" << frame1.LessSharp.size() << "," << frame1.LessSurf.size() << endl;
		//cout << "kdtree size:" << kdtreeSharp->size() << "," << kdtreeSurf->size() << endl;
		lastCornerKDTree_.setInputCloud(kdtreeSharp);
		lastSurfaceKDTree_.setInputCloud(kdtreeSurf);
		//cout << "kdtree size:" << lastCornerKDTree_.getInputCloud()->size()<< "," << lastSurfaceKDTree_.getInputCloud()->size() << endl;
		///< transformation from frame1 to frame2
		curr_last_rot = frame1.R_l_e.transpose()*frame2.R_l_e;
		curr_last_trans = frame1.R_l_e.transpose()*(frame2.t_l_e - frame1.t_l_e);

		cout << "releative pose:" << curr_last_rot.eulerAngles(0, 1, 2).transpose()*180.0 / M_PI << endl;
		cout << "releative transform:" << curr_last_trans.transpose() << endl;

		size_t lastCornerCloudSize = kdtreeSharp->points.size();
		size_t lastSurfaceCloudSize = kdtreeSurf->points.size();

		//cout << "kdtree size:" << lastCornerKDTree_.getInputCloud()->size() << "," << lastSurfaceKDTree_.getInputCloud()->size() << endl;
		if (lastCornerCloudSize > 10 && lastSurfaceCloudSize > 100)
		{
			std::vector<int> pointSearchInd(1);
			std::vector<float> pointSearchSqDis(1);

			size_t cornerPointsSharpNum = frame2.Sharp.points.size();
			size_t surfPointsFlatNum = frame2.Surf.points.size();

			///< position of feature which is projected to the last lidar frame
			pcl::PointXYZI pointSel;

			for (int i = 0; i < cornerPointsSharpNum; i++)
			{
				pointSel = pcl::PointXYZI();
				transformToStart(frame2.Sharp.points[i], pointSel);///< transform to the start time of the scan
				lastCornerKDTree_.nearestKSearch(pointSel, 1, pointSearchInd, pointSearchSqDis);///< find the closest point to pointSel in lastCornerCloud_


				int closestPointInd = -1;
				int minPointInd2 = -1;
				if (pointSearchSqDis[0] < 25)
				{

					closestPointInd = pointSearchInd[0];
					int closestPointScan = int(kdtreeSharp->points[closestPointInd].intensity);

					float pointSqDis, minPointSqDis2 = 25;
					for (int j = closestPointInd + 1; j < kdtreeSharp->points.size(); j++)
					{
						///< if in the same scan line, continue
						if (int(kdtreeSharp->points[j].intensity) <= closestPointScan)
							continue;

						///< if not in nearby scans, end the loop
						if (int(kdtreeSharp->points[j].intensity) > (closestPointScan + 2.5))
							break;

						pointSqDis = calcSquaredDiff(kdtreeSharp->points[j], pointSel);

						if (int(kdtreeSharp->points[j].intensity) > closestPointScan)
						{
							if (pointSqDis < minPointSqDis2)
							{
								minPointSqDis2 = pointSqDis;
								minPointInd2 = j;
							}
						}
					}
					for (int j = closestPointInd - 1; j >= 0; j--)
					{
						///< if in the same scan line, continue
						if (int(kdtreeSharp->points[j].intensity) >= closestPointScan)
							continue;

						if (int(kdtreeSharp->points[j].intensity) < closestPointScan - 2.5)
							break;

						pointSqDis = calcSquaredDiff(kdtreeSharp->points[j], pointSel);

						if (int(kdtreeSharp->points[j].intensity) < closestPointScan)
						{
							if (pointSqDis < minPointSqDis2)
							{
								minPointSqDis2 = pointSqDis;
								minPointInd2 = j;
							}
						}


					}
				}


				if (closestPointInd >= 0 && minPointInd2 >= 0) ///< both closestPointInd and minPointInd2 is valid
				{

					Eigen::Vector3d currCornerPoint(frame2.Sharp.points[i].x,
						frame2.Sharp.points[i].y,
						frame2.Sharp.points[i].z);
					Eigen::Vector3d correspondCornerPointA(kdtreeSharp->points[closestPointInd].x,
						kdtreeSharp->points[closestPointInd].y,
						kdtreeSharp->points[closestPointInd].z);
					Eigen::Vector3d correspondCornerPointB(kdtreeSharp->points[minPointInd2].x,
						kdtreeSharp->points[minPointInd2].y,
						kdtreeSharp->points[minPointInd2].z);

					pcl::PointXYZI currCornerPoint_inlast;
					transformToStart(frame2.Sharp.points[i], currCornerPoint_inlast);

					Eigen::Vector3d currCornerPoint_tolast(currCornerPoint_inlast.x,
						currCornerPoint_inlast.y,
						currCornerPoint_inlast.z);

					Eigen::Vector3d v_lj = correspondCornerPointA - correspondCornerPointB;
					Eigen::Vector3d v_li = currCornerPoint_tolast - correspondCornerPointB;
					Eigen::Vector3d v_lji = v_lj.cross(v_li);
					double length_lj = sqrt(v_lj.x()*v_lj.x() + v_lj.y()*v_lj.y() + v_lj.z()*v_lj.z());
					double length_ljxi = sqrt(v_lji.x()*v_lji.x() + v_lji.y()*v_lji.y() + v_lji.z()*v_lji.z());
					double d_e = length_ljxi / length_lj;
					//lsy
					if (d_e > 0.2)
					{
						continue;
					}
					correspondCornerFeature_.currCornerPointCloud.push_back(currCornerPoint);
					correspondCornerFeature_.currCornerPointCloud_inlast.push_back(currCornerPoint_tolast);
					correspondCornerFeature_.correspondCornerPointCloudA.push_back(correspondCornerPointA);
					correspondCornerFeature_.correspondCornerPointCloudB.push_back(correspondCornerPointB);

				}

			}

			for (int i = 0; i < surfPointsFlatNum; i++)
			{
				pointSearchInd.clear();
				pointSearchSqDis.clear();
				pointSel = pcl::PointXYZI();
				transformToStart(frame2.Surf.points[i], pointSel);

				lastSurfaceKDTree_.nearestKSearch(pointSel, 1, pointSearchInd, pointSearchSqDis);


				int closestPointInd = -1, minPointInd2 = -1, minPointInd3 = -1;
				if (pointSearchSqDis[0] < 25)
				{
					closestPointInd = pointSearchInd[0];
					int closestPointScan = int(kdtreeSurf->points[closestPointInd].intensity);

					float pointSqDis, minPointSqDis2 = 25, minPointSqDis3 = 25;
					for (int j = closestPointInd + 1; j < kdtreeSurf->points.size(); j++)
					{
						if (int(kdtreeSurf->points[j].intensity) > closestPointScan + 2.5)
						{
							break;
						}

						pointSqDis = calcSquaredDiff(kdtreeSurf->points[j], pointSel);

						if (int(kdtreeSurf->points[j].intensity) <= closestPointScan && pointSqDis < minPointSqDis2)
						{

							minPointSqDis2 = pointSqDis;
							minPointInd2 = j;

						}
						else if (int(kdtreeSurf->points[j].intensity) > closestPointScan && pointSqDis < minPointSqDis3)
						{

							minPointSqDis3 = pointSqDis;
							minPointInd3 = j;

						}
					}
					for (int j = closestPointInd - 1; j >= 0; j--)
					{
						if (int(kdtreeSurf->points[j].intensity) < closestPointScan - 2.5)
						{
							break;
						}

						pointSqDis = calcSquaredDiff(kdtreeSurf->points[j], pointSel);

						if (int(kdtreeSurf->points[j].intensity) >= closestPointScan && pointSqDis < minPointSqDis2)
						{

							minPointSqDis2 = pointSqDis;
							minPointInd2 = j;

						}
						else if (int(kdtreeSurf->points[j].intensity) < closestPointScan && pointSqDis < minPointSqDis3)
						{
							minPointSqDis3 = pointSqDis;
							minPointInd3 = j;
						}
					}

					if (minPointInd2 >= 0 && minPointInd3 >= 0)
					{

						Eigen::Vector3d currSurfPoint(frame2.Surf.points[i].x,
							frame2.Surf.points[i].y,
							frame2.Surf.points[i].z);
						Eigen::Vector3d correspondSurfPointA(kdtreeSurf->points[closestPointInd].x,
							kdtreeSurf->points[closestPointInd].y,
							kdtreeSurf->points[closestPointInd].z);
						Eigen::Vector3d correspondSurfPointB(kdtreeSurf->points[minPointInd2].x,
							kdtreeSurf->points[minPointInd2].y,
							kdtreeSurf->points[minPointInd2].z);
						Eigen::Vector3d correspondSurfPointC(kdtreeSurf->points[minPointInd3].x,
							kdtreeSurf->points[minPointInd3].y,
							kdtreeSurf->points[minPointInd3].z);

						pcl::PointXYZI currSurfPoint_inlast;
						transformToStart(frame2.Surf.points[i], currSurfPoint_inlast);
						Eigen::Vector3d currSurfPoint_tolast(currSurfPoint_inlast.x,
							currSurfPoint_inlast.y,
							currSurfPoint_inlast.z);

						Eigen::Vector3d v_jl = correspondSurfPointB - correspondSurfPointA;
						Eigen::Vector3d v_jm = correspondSurfPointC - correspondSurfPointA;
						Eigen::Vector3d v_ji = currSurfPoint_tolast - correspondSurfPointA;
						Eigen::Vector3d vj_lm = v_jl.cross(v_jm);
						Eigen::Vector3d up = Eigen::Vector3d(vj_lm(0)*v_ji(0), vj_lm(1)*v_ji(1), vj_lm(2)*v_ji(2));
						double dj_lm = sqrt(vj_lm(0)*vj_lm(0) + vj_lm(1)*vj_lm(1) + vj_lm(2)*vj_lm(2));
						double d_p = sqrt(up(0)*up(0) + up(1)*up(1) + up(2)*up(2)) / dj_lm;

						//lsy
						if (d_p > 0.1)
						{
							continue;
						}
						/*if (correspondSurfFeature_.currSurfPointCloud.size() > 0)
						{
							continue;
						}*/
						correspondSurfFeature_.currSurfPointCloud.push_back(currSurfPoint);
						correspondSurfFeature_.currSurfPointCloud_inlast.push_back(currSurfPoint_tolast);
						correspondSurfFeature_.correspondSurfPointCloudA.push_back(correspondSurfPointA);
						correspondSurfFeature_.correspondSurfPointCloudB.push_back(correspondSurfPointB);
						correspondSurfFeature_.correspondSurfPointCloudC.push_back(correspondSurfPointC);

					}
				}
			}
		}

		test_lastcorner = PointXYZI2Vector3d(frame1.LessSharp,frame1.R_l_e,frame1.t_l_e);
		test_lastsurf = PointXYZI2Vector3d(frame1.LessSurf, frame1.R_l_e, frame1.t_l_e);
		
		///< store the observation of lidar odometry
		lidarOdoObs.cornerFeature = correspondCornerFeature_;
		lidarOdoObs.surfFeature = correspondSurfFeature_;
		lidarOdoObs.last_R_l_e = frame1.R_l_e;
		lidarOdoObs.last_t_l_e = frame1.t_l_e;
		lidarOdoObs.curr_R_l_e = frame2.R_l_e;
		lidarOdoObs.curr_t_l_e = frame2.t_l_e;
	}

	vector<Eigen::Vector3d> t_glidarOdometry::PointXYZI2Vector3d(pcl::PointCloud<pcl::PointXYZI>& xyzi, Eigen::Matrix3d R_l_e, Eigen::Vector3d t_l_e)
	{
		vector<Eigen::Vector3d> xyz;
		for (int i = 0; i < xyzi.points.size(); i++)
		{
			xyz.push_back(R_l_e*Eigen::Vector3d(xyzi.points[i].x, xyzi.points[i].y, xyzi.points[i].z)+ t_l_e);
		}
		return xyz;
	}


	//void t_glidarOdometry::ceresOptimize()
	//{
	//	cout << "before ceres optimize" << endl;
	//	cout << "relative att:" << curr_last_rot.eulerAngles(0, 1, 2).transpose()*180/M_PI << endl;
	//	cout << "relative pos:" << curr_last_trans.transpose() << endl;
	//	//这里的优化对象实际上是两帧之间的相对位姿
	//	int surf_size = lidarOdoObs.surfFeature.currSurfPointCloud.size();
	//	int corner_size = lidarOdoObs.cornerFeature.currCornerPointCloud.size();
	//	/*cout << "surf_size:" << surf_size << endl;
	//	cout << "corner_size:" << corner_size << endl;
	//	getchar();*/

	//	CorrespondCornerFeature& cornerf = lidarOdoObs.cornerFeature;
	//	CorrespondSurfFeature& surff = lidarOdoObs.surfFeature;

	//	Eigen::Matrix3d last_R = lidarOdoObs.last_R_l_e;
	//	Eigen::Vector3d last_t = lidarOdoObs.last_t_l_e;

	//	double para_q[4] = { 0,0,0,1 };
	//	double para_t[3] = { 0,0,0};

	//	Eigen::Map<Eigen::Quaterniond> q_curr_last(para_q);
	//	Eigen::Map<Eigen::Vector3d> t_curr_last(para_t);

	//	q_curr_last= Eigen::Quaterniond(curr_last_rot);
	//	t_curr_last = curr_last_trans;

	//	clock_t start, stop;

	//	//build optimize problem
	//	ceres::LossFunction *loss_function = new ceres::HuberLoss(0.1);
	//	ceres::LocalParameterization *q_parameterization =
	//		new ceres::EigenQuaternionParameterization();
	//	ceres::Problem::Options problem_options;

	//	ceres::Problem problem(problem_options);
	//	problem.AddParameterBlock(para_q, 4, q_parameterization);
	//	problem.AddParameterBlock(para_t, 3);

	//	if (corner_size <= 0 || surf_size <= 0)
	//	{
	//		cout << "There are some problem in the lidar odometry observation!" << endl;
	//		getchar();
	//		return;
	//	}

	//	for (int i = 0; i < surf_size; i++)
	//	{
	//		//最后一个参数，经过雷达畸变矫正后设为1.0
	//		ceres::CostFunction *cost_function = LidarPlaneFactor::Create(surff.currSurfPointCloud[i], surff.correspondSurfPointCloudA[i], surff.correspondSurfPointCloudB[i], surff.correspondSurfPointCloudC[i], 1.0);
	//		problem.AddResidualBlock(cost_function, loss_function, para_q, para_t);
	//	}

	//	for (int i = 0; i < corner_size; i++)
	//	{
	//		ceres::CostFunction *cost_function = LidarEdgeFactor::Create(cornerf.currCornerPointCloud[i], cornerf.correspondCornerPointCloudA[i], cornerf.correspondCornerPointCloudB[i], 1.0);
	//		problem.AddResidualBlock(cost_function, loss_function, para_q, para_t);
	//	}

	//	start = clock();
	//	ceres::Solver::Options options;
	//	options.linear_solver_type = ceres::DENSE_QR;
	//	options.max_num_iterations = 4;
	//	options.minimizer_progress_to_stdout = false;
	//	ceres::Solver::Summary summary;
	//	ceres::Solve(options, &problem, &summary);
	//	stop = clock();
	//	cout<<summary.BriefReport()<<endl;
	//	cout << "optimize cost time:" << (double)(stop - start) / CLK_TCK << endl;

	//	lidarOdoObs.curr_R_l_e= last_R * q_curr_last.toRotationMatrix();
	//	lidarOdoObs.curr_t_l_e= last_t + last_R * t_curr_last;

	//	curr_last_rot = q_curr_last;
	//	curr_last_trans = t_curr_last;

	//	cout << "after ceres optimize" << endl;
	//	cout << "relative att:" <<curr_last_rot.eulerAngles(0, 1, 2).transpose() * 180 / M_PI << endl;
	//	cout << "relative pos:" << curr_last_trans.transpose() << endl;


	//}


}

