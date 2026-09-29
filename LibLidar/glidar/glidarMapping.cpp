#include "glidarMapping.h"

namespace glidar
{

	t_glidarMapping::t_glidarMapping()
	{
		laserCloudCornerFromMap.reset(new pcl::PointCloud<pcl::PointXYZI>());
		laserCloudSurfFromMap.reset(new pcl::PointCloud<pcl::PointXYZI>());

		///< kd-tree
		kdtreeCornerFromMap.reset(new pcl::KdTreeFLANN<pcl::PointXYZI>());
		kdtreeSurfFromMap.reset(new pcl::KdTreeFLANN<pcl::PointXYZI>());

		frameCount = 0;

		laserCloudCenWidth = 10;
		laserCloudCenHeight = 10;
		laserCloudCenDepth = 5;
		laserCloudWidth = 21;
		laserCloudHeight = 21;
		laserCloudDepth = 11;

		laserCloudNum = laserCloudWidth * laserCloudHeight * laserCloudDepth; ///< 4851

		systemInited_ = false;

		for (int i = 0; i < laserCloudNum; i++)
		{
			laserCloudSurfArray[i].reset(new pcl::PointCloud<pcl::PointXYZI>());
			laserCloudCornerArray[i].reset(new pcl::PointCloud<pcl::PointXYZI>());
		}

		laserCloudValidNum = 0;
		laserCloudSurroundNum = 0;
		for (int i = laserCloudCenWidth - 2; i <= laserCloudCenWidth + 2; i++)
		{
			for (int j = laserCloudCenHeight - 2; j <= laserCloudCenHeight + 2; j++)
			{
				for (int k = laserCloudCenDepth - 1; k <= laserCloudCenDepth + 1; k++)
				{
					///< save the id of cubes in the submap
					if (i >= 0 && i < laserCloudWidth &&
						j >= 0 && j < laserCloudHeight &&
						k >= 0 && k < laserCloudDepth)
					{
						laserCloudValidInd[laserCloudValidNum] = i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k;
						laserCloudValidNum++;
						laserCloudSurroundInd[laserCloudSurroundNum] = i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k;
						laserCloudSurroundNum++;
					}
				}
			}
		}

		///< Set the voxel size parameters used in filter processing
		downSizeFilterCorner.setLeafSize(0.4f, 0.4f, 0.4f);
		downSizeFilterSurf.setLeafSize(0.8f, 0.8f, 0.8f);
		downSizeMapCorner.setLeafSize(0.2, 0.2, 0.2);
		downSizeMapSurf.setLeafSize(0.4, 0.4, 0.4);

		last_R_l_e = Eigen::Matrix3d::Identity();
		last_t_l_e = Eigen::Vector3d::Zero();

		systemInited_ = false;

	}

	void t_glidarMapping::reset()
	{
		///< Clear the corresponding point set of the previous frame
		correspondCornerFeature_.currCornerPointCloud.clear();
		correspondCornerFeature_.correspondCornerPointCloudA.clear();
		correspondCornerFeature_.correspondCornerPointCloudB.clear();
		correspondCornerFeature_.currCornerPointCloud_inlast.clear();

		correspondSurfFeature_.currSurfPointCloud.clear();
		correspondSurfFeature_.correspondSurfPointCloudA.clear();
		correspondSurfFeature_.correspondSurfPointCloudB.clear();
		correspondSurfFeature_.correspondSurfPointCloudC.clear();
		correspondSurfFeature_.currSurfPointCloud_inlast.clear();
		correspondSurfFeature_.negative_OA_dot_norm.clear();
		correspondSurfFeature_.norm.clear();

		lidarMapObs = lidarMappingObs();
	}

	void t_glidarMapping::addPointcloudToMap(LidarFrame frame)
	{
		if (frame.empty)
			return;///< waiting for new data to arrive...

		if (!systemInited_)
		{
			///< store the first pose of lidar frame
			first_R_l_e = frame.R_l_e;
			first_t_l_e = frame.t_l_e;
			curr_R_l_e = frame.R_l_e;
			curr_t_l_e = frame.t_l_e;
			///< get the pose of current frame
			transformAssociateToMap();

			///< get the cube position of current frame center
			int centerCubeI = int((curr_t_l_w.x() + 25.0) / 50.0) + laserCloudCenWidth;
			int centerCubeJ = int((curr_t_l_w.y() + 25.0) / 50.0) + laserCloudCenHeight;
			int centerCubeK = int((curr_t_l_w.z() + 25.0) / 50.0) + laserCloudCenDepth;

			if (curr_t_l_w.x() + 25.0 < 0)
				centerCubeI--;
			if (curr_t_l_w.y() + 25.0 < 0)
				centerCubeJ--;
			if (curr_t_l_w.z() + 25.0 < 0)
				centerCubeK--;

			///< adjust the range of map
			firstZoom(centerCubeI, centerCubeJ, centerCubeK);

			///< add points to correspond map cube
			AddCornerPointToMap(frame.LessSharp);
			AddSurfPointToMap(frame.LessSurf);

			last_R_l_e = curr_R_l_e;
			last_t_l_e = curr_t_l_e;

			systemInited_ = true;
			return;

		}

		curr_R_l_e = frame.R_l_e;
		curr_t_l_e = frame.t_l_e;
		///< get the pose of current frame
		transformAssociateToMap();
		///< get the cube position of current frame center
		int centerCubeI = int((curr_t_l_w.x() + 25.0) / 50.0) + laserCloudCenWidth;
		int centerCubeJ = int((curr_t_l_w.y() + 25.0) / 50.0) + laserCloudCenHeight;
		int centerCubeK = int((curr_t_l_w.z() + 25.0) / 50.0) + laserCloudCenDepth;

		if (curr_t_l_w.x() + 25.0 < 0)
			centerCubeI--;
		if (curr_t_l_w.y() + 25.0 < 0)
			centerCubeJ--;
		if (curr_t_l_w.z() + 25.0 < 0)
			centerCubeK--;

		///< adjust the range of current frame center
		Zoom(centerCubeI, centerCubeJ, centerCubeK);
		AddCornerPointToMap(frame.LessSharp);
		AddSurfPointToMap(frame.LessSurf);
		downSizeMap();
	}

	void t_glidarMapping::renewWframe(Eigen::Matrix3d R, Eigen::Vector3d t)
	{
		Eigen::Matrix3d R_tonewW = R.transpose()*first_R_l_e;
		Eigen::Vector3d t_tonewW = R.transpose()*(first_t_l_e - t);

		///< transfrom map points to new W frame
		toNewW(R_tonewW, t_tonewW);

		first_R_l_e = R;
		first_t_l_e = t;
	}

	void t_glidarMapping::transPointXYZI(pcl::PointXYZI& p, Eigen::Matrix3d& R, Eigen::Vector3d& t)
	{
		Eigen::Vector3d temp_p = R * Eigen::Vector3d(p.x, p.y, p.z) + t;
		p.x = temp_p.x();
		p.y = temp_p.y();
		p.z = temp_p.z();
	}

	void t_glidarMapping::transPointXYZI(pcl::PointCloud<pcl::PointXYZI>::Ptr xyzi, Eigen::Matrix3d& R, Eigen::Vector3d& t)
	{
		for (int i = 0; i < xyzi->points.size(); i++)
		{
			transPointXYZI(xyzi->points[i], R, t);
		}
	}

	void t_glidarMapping::toNewW(Eigen::Matrix3d R_tonewW, Eigen::Vector3d t_tonewW)
	{
		for (int i = 0; i < laserCloudValidNum; i++)
		{
			int id = laserCloudValidInd[i];
			//cout << "id:" << id << endl;
			for (int j = 0; j < laserCloudSurfArray[id]->points.size(); j++)
			{
				transPointXYZI(laserCloudSurfArray[id]->points[j], R_tonewW, t_tonewW);
			}
			for (int k = 0; k < laserCloudCornerArray[id]->points.size(); k++)
			{
				transPointXYZI(laserCloudCornerArray[id]->points[k], R_tonewW, t_tonewW);
			}
		}
	}

	void t_glidarMapping::process(LidarFrame frame,bool lastmap)
	{
		if (frame.empty)
		{
			std::cout << "ERROR:the pose of lidar frame need to be initialized!" << endl;
			getchar();
		}
		reset();

		if (!systemInited_)
		{
			///< store the first frame pose
			first_R_l_e = frame.R_l_e;
			first_t_l_e = frame.t_l_e;
			curr_R_l_e = frame.R_l_e;
			curr_t_l_e = frame.t_l_e;
			///< get the pose of current lidar frame
			transformAssociateToMap();

			///< get the cube position of current frame center
			int centerCubeI = int((curr_t_l_w.x() + 25.0) / 50.0) + laserCloudCenWidth;
			int centerCubeJ = int((curr_t_l_w.y() + 25.0) / 50.0) + laserCloudCenHeight;
			int centerCubeK = int((curr_t_l_w.z() + 25.0) / 50.0) + laserCloudCenDepth;

			if (curr_t_l_w.x() + 25.0 < 0)
				centerCubeI--;
			if (curr_t_l_w.y() + 25.0 < 0)
				centerCubeJ--;
			if (curr_t_l_w.z() + 25.0 < 0)
				centerCubeK--;

			///< adjust the range of lidar map
			firstZoom(centerCubeI, centerCubeJ, centerCubeK);

			///< add points to correspond cube of map
			AddCornerPointToMap(frame.LessSharp);
			AddSurfPointToMap(frame.LessSurf);

			last_R_l_e = curr_R_l_e;
			last_t_l_e = curr_t_l_e;
			systemInited_ = true;
			return;

		}

		curr_R_l_e = frame.R_l_e;
		curr_t_l_e = frame.t_l_e;
		///< get the pose of current frame
		transformAssociateToMap();
	    ///< get the cube position of current frame center
		int centerCubeI = int((curr_t_l_w.x() + 25.0) / 50.0) + laserCloudCenWidth;
		int centerCubeJ = int((curr_t_l_w.y() + 25.0) / 50.0) + laserCloudCenHeight;
		int centerCubeK = int((curr_t_l_w.z() + 25.0) / 50.0) + laserCloudCenDepth;

		if (curr_t_l_w.x() + 25.0 < 0)
			centerCubeI--;
		if (curr_t_l_w.y() + 25.0 < 0)
			centerCubeJ--;
		if (curr_t_l_w.z() + 25.0 < 0)
			centerCubeK--;

		///< adjust the range of lidar map
		Zoom(centerCubeI, centerCubeJ, centerCubeK);
		///< find the correspond submap of lidar map
		FindSubmap(centerCubeI, centerCubeJ, centerCubeK);

		///< points in the choosed submap cubes
		laserCloudCornerFromMap->clear();
		laserCloudSurfFromMap->clear();
		for (int i = 0; i < laserCloudValidNum; i++)
		{
			*laserCloudCornerFromMap += *laserCloudCornerArray[laserCloudValidInd[i]];
			*laserCloudSurfFromMap += *laserCloudSurfArray[laserCloudValidInd[i]];
		}

		///< if use "lastmap", we transform the submap to the coordinate of last frame
		if (lastmap)
		{
			transPointXYZI(laserCloudCornerFromMap,last_R_w_l,last_t_w_l);
			transPointXYZI(laserCloudSurfFromMap, last_R_w_l, last_t_w_l);
		}


		int laserCloudCornerFromMapNum = laserCloudCornerFromMap->points.size();
		int laserCloudSurfFromMapNum = laserCloudSurfFromMap->points.size();

		std::vector<int> index;
		pcl::PointCloud<pcl::PointXYZI>::Ptr inputSurf(new pcl::PointCloud<pcl::PointXYZI>());
		pcl::PointCloud<pcl::PointXYZI>::Ptr inputSharp(new pcl::PointCloud<pcl::PointXYZI>());
		pcl::removeNaNFromPointCloud(frame.Sharp, *inputSharp, index);
		pcl::removeNaNFromPointCloud(frame.Surf, *inputSurf, index);

		///< downsize the feature points
		pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudCornerStack(new pcl::PointCloud<pcl::PointXYZI>());
		downSizeFilterCorner.setInputCloud(inputSharp);
		downSizeFilterCorner.filter(*laserCloudCornerStack);
		int laserCloudCornerStackNum = laserCloudCornerStack->points.size();

		pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudSurfStack(new pcl::PointCloud<pcl::PointXYZI>());
		downSizeFilterSurf.setInputCloud(inputSurf);
		downSizeFilterSurf.filter(*laserCloudSurfStack);
		int laserCloudSurfStackNum = laserCloudSurfStack->points.size();

		///< construct the 
		if (laserCloudCornerFromMapNum > 10 && laserCloudSurfFromMapNum > 50)
		{
			cout << "map size(surf/corner):" << laserCloudSurfFromMapNum << "," << laserCloudCornerFromMapNum << endl;
			kdtreeCornerFromMap->setInputCloud(laserCloudCornerFromMap);
			kdtreeSurfFromMap->setInputCloud(laserCloudSurfFromMap);
		}
		else
		{
			cout << "not enough points in laser map!" << endl;
			getchar();
		}

		///< find the corresponding points for corner features
		for (int i = 0; i < laserCloudCornerStackNum; i++)
		{
			///< transform the point to w frame
			pointOri = laserCloudCornerStack->points[i];
			if (lastmap)
				pointAssociateToMap(pointOri, pointSel,true);
			else
				pointAssociateToMap(pointOri, pointSel,false);

			kdtreeCornerFromMap->nearestKSearch(pointSel, 5, pointSearchInd, pointSearchSqDis);

			if (pointSearchSqDis[4] < 0.5)
			{
				std::vector<Eigen::Vector3d> nearCorners;
				Eigen::Vector3d center(0, 0, 0);
				for (int j = 0; j < 5; j++)
				{
					Eigen::Vector3d tmp(laserCloudCornerFromMap->points[pointSearchInd[j]].x,
						laserCloudCornerFromMap->points[pointSearchInd[j]].y,
						laserCloudCornerFromMap->points[pointSearchInd[j]].z);
					center = center + tmp;
					nearCorners.push_back(tmp);
				}
				center = center / 5.0;

				Eigen::Matrix3d covMat = Eigen::Matrix3d::Zero();
				for (int j = 0; j < 5; j++)
				{
					Eigen::Matrix<double, 3, 1> tmpZeroMean = nearCorners[j] - center;
					covMat = covMat + tmpZeroMean * tmpZeroMean.transpose();
				}

				Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> saes(covMat);

				// if is indeed line feature
				// note Eigen library sort eigenvalues in increasing order
				Eigen::Vector3d unit_direction = saes.eigenvectors().col(2);
				Eigen::Vector3d curr_point(pointOri.x, pointOri.y, pointOri.z);
				Eigen::Vector3d curr_point_w(pointSel.x, pointSel.y, pointSel.z);
				if (saes.eigenvalues()[2] > 3 * saes.eigenvalues()[1])
				{
					Eigen::Vector3d point_on_line = center;
					Eigen::Vector3d point_a, point_b;
					point_a = 0.02 * unit_direction + point_on_line;
					point_b = -0.02 * unit_direction + point_on_line;

					Eigen::Vector3d v_lj = point_a - point_b;
					Eigen::Vector3d v_li = curr_point_w - point_b;
					Eigen::Vector3d v_lji = v_lj.cross(v_li);
					double length_lj = sqrt(v_lj.x()*v_lj.x() + v_lj.y()*v_lj.y() + v_lj.z()*v_lj.z());
					double length_ljxi = sqrt(v_lji.x()*v_lji.x() + v_lji.y()*v_lji.y() + v_lji.z()*v_lji.z());
					double d_e = length_ljxi / length_lj;

					if (d_e > 0.2)
					{
						continue;
					}

					correspondCornerFeature_.currCornerPointCloud.push_back(curr_point);
					correspondCornerFeature_.currCornerPointCloud_inlast.push_back(curr_point_w);
					correspondCornerFeature_.correspondCornerPointCloudA.push_back(point_a);
					correspondCornerFeature_.correspondCornerPointCloudB.push_back(point_b);
				}
			}
			/*
			else if(pointSearchSqDis[4] < 0.01 * sqrtDis)
			{
				Eigen::Vector3d center(0, 0, 0);
				for (int j = 0; j < 5; j++)
				{
					Eigen::Vector3d tmp(laserCloudCornerFromMap->points[pointSearchInd[j]].x,
										laserCloudCornerFromMap->points[pointSearchInd[j]].y,
										laserCloudCornerFromMap->points[pointSearchInd[j]].z);
					center = center + tmp;
				}
				center = center / 5.0;
				Eigen::Vector3d curr_point(pointOri.x, pointOri.y, pointOri.z);
				ceres::CostFunction *cost_function = LidarDistanceFactor::Create(curr_point, center);
				problem.AddResidualBlock(cost_function, loss_function, parameters, parameters + 4);
			}
			*/
		}

		///< find the correspond points for surf features
		for (int i = 0; i < laserCloudSurfStackNum; i++)
		{
			pointOri = laserCloudSurfStack->points[i];
			if(lastmap)
				pointAssociateToMap(pointOri, pointSel,true);
			else
				pointAssociateToMap(pointOri, pointSel,false);

			kdtreeSurfFromMap->nearestKSearch(pointSel, 5, pointSearchInd, pointSearchSqDis);

			Eigen::Matrix<double, 5, 3> matA0;
			Eigen::Matrix<double, 5, 1> matB0 = -1 * Eigen::Matrix<double, 5, 1>::Ones();
			if (pointSearchSqDis[4] < 1.0)
			{

				for (int j = 0; j < 5; j++)
				{
					matA0(j, 0) = laserCloudSurfFromMap->points[pointSearchInd[j]].x;
					matA0(j, 1) = laserCloudSurfFromMap->points[pointSearchInd[j]].y;
					matA0(j, 2) = laserCloudSurfFromMap->points[pointSearchInd[j]].z;
				}
				///< find the norm of plane
				Eigen::Vector3d norm = matA0.colPivHouseholderQr().solve(matB0);
				double negative_OA_dot_norm = 1 / norm.norm();
				norm.normalize();

				///< Here n(pa, pb, pc) is unit norm of plane
				bool planeValid = true;
				for (int j = 0; j < 5; j++)
				{
					///< if OX * n larger than 0.2, then plane is not fit well
					if (fabs(norm(0) * laserCloudSurfFromMap->points[pointSearchInd[j]].x +
						norm(1) * laserCloudSurfFromMap->points[pointSearchInd[j]].y +
						norm(2) * laserCloudSurfFromMap->points[pointSearchInd[j]].z + negative_OA_dot_norm) > 0.2)
					{
						planeValid = false;
						break;
					}
				}
				Eigen::Vector3d curr_point(pointOri.x, pointOri.y, pointOri.z);
				Eigen::Vector3d curr_point_w(pointSel.x, pointSel.y, pointSel.z);
				if (planeValid)
				{
					double residual = norm.dot(curr_point_w) + negative_OA_dot_norm;
					if (residual > 0.1)
						continue;

					correspondSurfFeature_.currSurfPointCloud.push_back(curr_point);
					correspondSurfFeature_.currSurfPointCloud_inlast.push_back(curr_point_w);
					correspondSurfFeature_.norm.push_back(norm);
					correspondSurfFeature_.negative_OA_dot_norm.push_back(negative_OA_dot_norm);
				}
			}

		}
		
		lidarMapObs.cornerFeature = correspondCornerFeature_;
		lidarMapObs.surfFeature = correspondSurfFeature_;
		lidarMapObs.curr_R_l_e = frame.R_l_e;
		lidarMapObs.curr_t_l_e = frame.t_l_e;
		lidarMapObs.last_R_l_e = last_R_l_e;
		lidarMapObs.last_t_l_e = last_t_l_e;

		last_R_l_e = curr_R_l_e;
		last_t_l_e = curr_t_l_e;

		///< used for visibility
		laserCloudCornerFromMap->clear();
		laserCloudSurfFromMap->clear();
		for (int i = 0; i < laserCloudValidNum; i++)
		{
			*laserCloudCornerFromMap += *laserCloudCornerArray[laserCloudValidInd[i]];
			*laserCloudSurfFromMap += *laserCloudSurfArray[laserCloudValidInd[i]];
		}
		submap_surf = PointXYZI2Vector3d(laserCloudSurfFromMap);
		submap_corner = PointXYZI2Vector3d(laserCloudCornerFromMap);
	}

	vector<Eigen::Vector3d> t_glidarMapping::PointXYZI2Vector3d(pcl::PointCloud<pcl::PointXYZI>::Ptr xyzi)
	{
		vector<Eigen::Vector3d> xyz;
		for (int i = 0; i < xyzi->points.size(); i++)
		{
			xyz.push_back(first_R_l_e*Eigen::Vector3d(xyzi->points[i].x, xyzi->points[i].y, xyzi->points[i].z)+first_t_l_e);
		}
		return xyz;
	}

	void t_glidarMapping::transformAssociateToMap()
	{
		curr_R_l_w = first_R_l_e.transpose()*curr_R_l_e;
		curr_t_l_w = first_R_l_e.transpose()*(curr_t_l_e-first_t_l_e);
		curr_R_w_l = curr_R_l_e.transpose()*first_R_l_e;
		curr_t_w_l = curr_R_l_e.transpose()*(first_t_l_e - curr_t_l_e);
		curr_last_R = last_R_l_e.transpose()*curr_R_l_e;
		curr_last_t = last_R_l_e.transpose()*(curr_t_l_e - last_t_l_e);
		last_R_w_l= last_R_l_e.transpose()*first_R_l_e;
		last_t_w_l = last_R_l_e.transpose()*(first_t_l_e - last_t_l_e);

	}

	void t_glidarMapping::pointAssociateToMap(pcl::PointXYZI & pi, pcl::PointXYZI & po,bool lastmap)
	{
		Eigen::Vector3d point_curr(pi.x, pi.y, pi.z);
		
		Eigen::Vector3d point_w;
		if (lastmap)
		{
			point_w = curr_last_R * point_curr + curr_last_t;
		}
		else
		{
			point_w = curr_R_l_w * point_curr + curr_t_l_w;
		}
		
		po.x = point_w.x();
		po.y = point_w.y();
		po.z = point_w.z();
		po.intensity = pi.intensity;
	}

	void t_glidarMapping::AddCornerPointToMap(pcl::PointCloud<pcl::PointXYZI> cornerPoints)
	{
		for (int i = 0; i < cornerPoints.points.size(); i++)
		{
			pointAssociateToMap(cornerPoints.points[i], pointSel,false);

			int cubeI = int((pointSel.x + 25.0) / 50.0) + laserCloudCenWidth;
			int cubeJ = int((pointSel.y + 25.0) / 50.0) + laserCloudCenHeight;
			int cubeK = int((pointSel.z + 25.0) / 50.0) + laserCloudCenDepth;

			if (pointSel.x + 25.0 < 0)
				cubeI--;
			if (pointSel.y + 25.0 < 0)
				cubeJ--;
			if (pointSel.z + 25.0 < 0)
				cubeK--;

			cubeI = cubeI - lastCenterCubeI + laserCloudCenWidth;
			cubeJ = cubeJ - lastCenterCubeJ + laserCloudCenHeight;
			cubeK = cubeK - lastCenterCubeK + laserCloudCenDepth;

			if (cubeI >= 0 && cubeI < laserCloudWidth &&
				cubeJ >= 0 && cubeJ < laserCloudHeight &&
				cubeK >= 0 && cubeK < laserCloudDepth)
			{
				int cubeInd = cubeI + laserCloudWidth * cubeJ + laserCloudWidth * laserCloudHeight * cubeK;
				laserCloudCornerArray[cubeInd]->push_back(pointSel);
			}
		}
	}

	void t_glidarMapping::AddSurfPointToMap(pcl::PointCloud<pcl::PointXYZI> surfPoints)
	{
		for (int i = 0; i < surfPoints.points.size(); i++)
		{
			pointAssociateToMap(surfPoints.points[i], pointSel,false);

			int cubeI = int((pointSel.x + 25.0) / 50.0) + laserCloudCenWidth;
			int cubeJ = int((pointSel.y + 25.0) / 50.0) + laserCloudCenHeight;
			int cubeK = int((pointSel.z + 25.0) / 50.0) + laserCloudCenDepth;

			if (pointSel.x + 25.0 < 0)
				cubeI--;
			if (pointSel.y + 25.0 < 0)
				cubeJ--;
			if (pointSel.z + 25.0 < 0)
				cubeK--;

			cubeI = cubeI - lastCenterCubeI + laserCloudCenWidth;
			cubeJ = cubeJ - lastCenterCubeJ + laserCloudCenHeight;
			cubeK = cubeK - lastCenterCubeK + laserCloudCenDepth;

			if (cubeI >= 0 && cubeI < laserCloudWidth &&
				cubeJ >= 0 && cubeJ < laserCloudHeight &&
				cubeK >= 0 && cubeK < laserCloudDepth)
			{
				int cubeInd = cubeI + laserCloudWidth * cubeJ + laserCloudWidth * laserCloudHeight * cubeK;
				laserCloudSurfArray[cubeInd]->push_back(pointSel);
			}
		}
	}

	void t_glidarMapping::downSizeMap()
	{
		for (int i = 0; i < laserCloudValidNum; i++)
		{
			int ind = laserCloudValidInd[i];

			if (laserCloudCornerArray[ind]->size() < 1)
				continue;

			pcl::PointCloud<pcl::PointXYZI>::Ptr tmpCorner(new pcl::PointCloud<pcl::PointXYZI>());
			downSizeMapCorner.setInputCloud(laserCloudCornerArray[ind]);
			downSizeMapCorner.filter(*tmpCorner);
			laserCloudCornerArray[ind] = tmpCorner;

			pcl::PointCloud<pcl::PointXYZI>::Ptr tmpSurf(new pcl::PointCloud<pcl::PointXYZI>());
			downSizeMapSurf.setInputCloud(laserCloudSurfArray[ind]);
			downSizeMapSurf.filter(*tmpSurf);
			laserCloudSurfArray[ind] = tmpSurf;
		}
	}

	//void t_glidarMapping::ceresOptimize()
	//{
	//	int surf_size = lidarMapObs.surfFeature.currSurfPointCloud.size();
	//	int corner_size = lidarMapObs.cornerFeature.currCornerPointCloud.size();
	//	CorrespondCornerFeature& cornerf = lidarMapObs.cornerFeature;
	//	CorrespondSurfFeature& surff = lidarMapObs.surfFeature;

	//	//q_curr_last(x, y, z, w), t_curr_last

	//	Eigen::Quaterniond q_c_w(curr_R_l_w);
	//	Eigen::Vector3d t_c_w = curr_t_l_w;
	//	double parameters[7] = { q_c_w.x(), q_c_w.y(), q_c_w.z(), q_c_w.w(), t_c_w.x(), t_c_w.y(), t_c_w.z() };
	//	Eigen::Map<Eigen::Quaterniond> q_w_curr(parameters);
	//	Eigen::Map<Eigen::Vector3d> t_w_curr(parameters + 4);

	//	//build optimize problem
	//	ceres::LossFunction *loss_function = new ceres::HuberLoss(0.1);
	//	ceres::LocalParameterization *q_parameterization =
	//		new ceres::EigenQuaternionParameterization();
	//	ceres::Problem::Options problem_options;

	//	ceres::Problem problem(problem_options);
	//	problem.AddParameterBlock(parameters, 4, q_parameterization);
	//	problem.AddParameterBlock(parameters + 4, 3);

	//	if (corner_size <= 0 || surf_size <= 0)
	//	{
	//		cout << "There are some problem in the lidar odometry observation!" << endl;
	//		return;
	//	}

	//	for (int i = 0; i < surf_size; i++)
	//	{
	//		//最后一个参数，经过雷达畸变矫正后设为1.0
	//		ceres::CostFunction *cost_function = LidarPlaneNormFactor::Create(surff.currSurfPointCloud[i], surff.norm[i],surff.negative_OA_dot_norm[i]);
	//		problem.AddResidualBlock(cost_function, loss_function, parameters, parameters + 4);
	//		//surf_correspondence++;
	//	}

	//	for (int i = 0; i < corner_size; i++)
	//	{
	//		ceres::CostFunction *cost_function = LidarEdgeFactor::Create(cornerf.currCornerPointCloud[i], cornerf.correspondCornerPointCloudA[i], cornerf.correspondCornerPointCloudB[i], 1.0);
	//		problem.AddResidualBlock(cost_function, loss_function, parameters, parameters + 4);
	//		//corner_correspondence++;
	//	}

	//	ceres::Solver::Options options;
	//	options.linear_solver_type = ceres::DENSE_QR;
	//	options.max_num_iterations = 4;
	//	options.minimizer_progress_to_stdout = false;
	//	options.check_gradients = false;
	//	options.gradient_check_relative_precision = 1e-4;
	//	ceres::Solver::Summary summary;
	//	ceres::Solve(options, &problem, &summary);

	//	//renew the pose of current lidar frame
	//	/*Matrix3d curr_R_l_e = last_rot * q_curr_last.toRotationMatrix();
	//	Vector3d curr_t_l_e = last_trans + last_rot * t_curr_last;*/

	//}

	

	void t_glidarMapping::firstZoom(int centerCubeI, int centerCubeJ, int centerCubeK)
	{
		//put the cube which contains the center as the cube center ->(10,10,5)
		lastCenterCubeI = centerCubeI;
		lastCenterCubeJ = centerCubeJ;
		lastCenterCubeK = centerCubeK;
	}
	void t_glidarMapping::Zoom(int centerCubeI, int centerCubeJ, int centerCubeK)
	{
		int diffI, diffJ, diffK;
		diffI = centerCubeI - lastCenterCubeI;
		diffJ = centerCubeJ - lastCenterCubeJ;
		diffK = centerCubeK - lastCenterCubeK;

		if (diffI != 0 || diffJ != 0 || diffK != 0)
		{
			std::cout << "last center:" << lastCenterCubeI << "," << lastCenterCubeJ << "," << lastCenterCubeK << endl;
			std::cout << "now center:" << centerCubeI << "," << centerCubeJ << "," << centerCubeK << endl;
			std::cout << "Zoom diff:" << diffI << "," << diffJ << "," << diffK << endl;
		}

		///< to asure there are 5*5*5 adjacent cubes arround (centerCubeI，centerCubeJ，centerCubeK)
		while (diffI < 0)
		{
			for (int j = 0; j < laserCloudHeight; j++)
			{
				for (int k = 0; k < laserCloudDepth; k++)
				{
					int i = laserCloudWidth - 1;
					pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudCubeCornerPointer =
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
					pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudCubeSurfPointer =
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
					
					for (; i >= 1; i--)
					{
						int surf_size = laserCloudSurfArray[i - 1 + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->points.size();
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCornerArray[i - 1 + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudSurfArray[i - 1 + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						if (surf_size != laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->points.size())
						{
							cout << "before move:" << surf_size << endl;
							cout << "after move:" << laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->points.size() << endl;
							std::cout << "Wrong:diffI<0" << endl;
							getchar();
						}
					}
					
					laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->clear();
					laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->clear();
					
					laserCloudCubeCornerPointer->clear();
					laserCloudCubeSurfPointer->clear();
				}
			}

			diffI++;
		}

		while (diffI > 0)
		{
			for (int j = 0; j < laserCloudHeight; j++)
			{
				for (int k = 0; k < laserCloudDepth; k++)
				{
					int i = 0;
					pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudCubeCornerPointer =
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
					pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudCubeSurfPointer =
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
					
					for (; i < laserCloudWidth - 1; i++)
					{
						int surf_size = laserCloudSurfArray[i + 1 + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->points.size();
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCornerArray[i + 1 + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudSurfArray[i + 1 + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						if (surf_size != laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->points.size())
						{
							cout << "Wrong:diffI>0" << endl;
							getchar();
						}
					}
					
					laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->clear();
					laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->clear();
		
					laserCloudCubeCornerPointer->clear();
					laserCloudCubeSurfPointer->clear();
				}
			}

			diffI--;
		}

		while (diffJ < 0)
		{
			for (int i = 0; i < laserCloudWidth; i++)
			{
				for (int k = 0; k < laserCloudDepth; k++)
				{
					int j = laserCloudHeight - 1;
					pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudCubeCornerPointer =
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
					pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudCubeSurfPointer =
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
					
					for (; j >= 1; j--)
					{
						int surf_size = laserCloudSurfArray[i + laserCloudWidth * (j - 1) + laserCloudWidth * laserCloudHeight * k]->points.size();
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCornerArray[i + laserCloudWidth * (j - 1) + laserCloudWidth * laserCloudHeight * k];
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudSurfArray[i + laserCloudWidth * (j - 1) + laserCloudWidth * laserCloudHeight * k];
						if (surf_size != laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->points.size())
						{
							std::cout << "Wrong:diffJ<0" << endl;
							getchar();
						}
					}
					
					laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->clear();
					laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->clear();
		
					laserCloudCubeCornerPointer->clear();
					laserCloudCubeSurfPointer->clear();
				}
			}

			diffJ++;
		}

		while (diffJ >0)
		{
			for (int i = 0; i < laserCloudWidth; i++)
			{
				for (int k = 0; k < laserCloudDepth; k++)
				{
					int j = 0;
					pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudCubeCornerPointer =
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
					pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudCubeSurfPointer =
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
					
					for (; j < laserCloudHeight - 1; j++)
					{
						int surf_size = laserCloudSurfArray[i + laserCloudWidth * (j + 1) + laserCloudWidth * laserCloudHeight * k]->points.size();
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCornerArray[i + laserCloudWidth * (j + 1) + laserCloudWidth * laserCloudHeight * k];
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudSurfArray[i + laserCloudWidth * (j + 1) + laserCloudWidth * laserCloudHeight * k];
						if (surf_size != laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->points.size())
						{
							std::cout << "Wrong:diffJ>0" << endl;
							getchar();
						}
					}
					
					laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->clear();
					laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->clear();
						
					laserCloudCubeCornerPointer->clear();
					laserCloudCubeSurfPointer->clear();
				}
			}

			diffJ--;
		}

		while (diffK < 0)
		{
			for (int i = 0; i < laserCloudWidth; i++)
			{
				for (int j = 0; j < laserCloudHeight; j++)
				{
					int k = laserCloudDepth - 1;
					pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudCubeCornerPointer =
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
					pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudCubeSurfPointer =
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
					
					for (; k >= 1; k--)
					{
						int surf_size = laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * (k - 1)]->points.size();
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * (k - 1)];
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * (k - 1)];
						if (surf_size != laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->points.size())
						{
							std::cout << "Wrong:diffK<0" << endl;
							getchar();
						}
					}
					
					laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->clear();
					laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->clear();
				
					laserCloudCubeCornerPointer->clear();
					laserCloudCubeSurfPointer->clear();
				}
			}

			diffK++;
		}

		while (diffK >0)
		{
			for (int i = 0; i < laserCloudWidth; i++)
			{
				for (int j = 0; j < laserCloudHeight; j++)
				{
					int k = 0;
					pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudCubeCornerPointer =
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
					pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudCubeSurfPointer =
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];

					
					for (; k < laserCloudDepth - 1; k++)
					{
						int surf_size = laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * (k + 1)]->points.size();
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * (k + 1)];
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * (k + 1)];
						if (surf_size != laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->points.size())
						{
							std::cout << "Wrong:diffK>0" << endl;
							getchar();
						}
					}
					
					laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->clear();
					laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k]->clear();
					
					laserCloudCubeCornerPointer->clear();
					laserCloudCubeSurfPointer->clear();
				}
			}

			diffK--;
		}

		lastCenterCubeI = centerCubeI;
		lastCenterCubeJ = centerCubeJ;
		lastCenterCubeK = centerCubeK;
	}

	void t_glidarMapping::FindSubmap(int centerCubeI, int centerCubeJ, int centerCubeK)
	{
		laserCloudValidNum = 0;
		laserCloudSurroundNum = 0;
		centerCubeI = laserCloudCenWidth;
		centerCubeJ = laserCloudCenHeight;
		centerCubeK = laserCloudCenDepth;
		///< find the submap
		for (int i = centerCubeI - 2; i <= centerCubeI + 2; i++)
		{
			for (int j = centerCubeJ - 2; j <= centerCubeJ + 2; j++)
			{
				for (int k = centerCubeK - 1; k <= centerCubeK + 1; k++)
				{
					///< save the id of cubes in the submap
					if (i >= 0 && i < laserCloudWidth &&
						j >= 0 && j < laserCloudHeight &&
						k >= 0 && k < laserCloudDepth)
					{
						laserCloudValidInd[laserCloudValidNum] = i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k;
						laserCloudValidNum++;
						laserCloudSurroundInd[laserCloudSurroundNum] = i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k;
						laserCloudSurroundNum++;
					}
				}
			}
		}


	}

	Eigen::Vector3d t_glidarMapping::translateToE(Eigen::Vector3d point_w)
	{
		return (first_R_l_e * point_w + first_t_l_e);
	}

}

