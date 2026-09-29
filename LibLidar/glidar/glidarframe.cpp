
#include "glidarframe.h"


namespace glidar
{ 

	float cloudCurvature[400000];
	int cloudSortInd[400000];
	int cloudNeighborPicked[400000];
	int cloudLabel[400000];

	glidar::t_glidarframe::t_glidarframe(gnut::t_gsetbase* _set) :t_glidarproc(_set)
	{
		n_scans = dynamic_cast<t_gsetlidar*>(_set)->n_scans();
		close_threshold = dynamic_cast<t_gsetlidar*>(_set)->close_threshold();
		scanPeriod = dynamic_cast<t_gsetlidar*>(_set)->scan_period();
	}

	bool comp(int i, int j) { return (cloudCurvature[i] < cloudCurvature[j]); }


	void t_glidarframe::handleCloud(pcl::PointCloud<pcl::PointXYZI> &laserCloudIn)
	{
		clear();
		process(laserCloudIn);
	}

	void t_glidarframe::clear()
	{
		_laserCloudScans.clear();
		_laserCloud.clear();
		_scanIndices.clear();
		_cornerPointsSharp.clear();
		_cornerPointsLessSharp.clear();
		_surfacePointsFlat.clear();
		_surfacePointsLessFlat.clear();
	}
	
	void t_glidarframe::removeClosedPointCloud(pcl::PointCloud<pcl::PointXYZI> &cloud_in,
		pcl::PointCloud<pcl::PointXYZI> &cloud_out, float thres)
	{
		if (&cloud_in != &cloud_out)
		{
			cloud_out.header = cloud_in.header;
			cloud_out.points.resize(cloud_in.points.size());
		}

		size_t j = 0;

		for (size_t i = 0; i < cloud_in.points.size(); ++i)
		{
			if (cloud_in.points[i].x * cloud_in.points[i].x + cloud_in.points[i].y * cloud_in.points[i].y + cloud_in.points[i].z * cloud_in.points[i].z < thres * thres)
				continue;
			cloud_out.points[j] = cloud_in.points[i];
			j++;
		}
		if (j != cloud_in.points.size())
		{
			cloud_out.points.resize(j);
		}

		cloud_out.height = 1;
		cloud_out.width = static_cast<uint32_t>(j);
		cloud_out.is_dense = true;
	}

	//read lidar file (format:txt)
	pcl::PointCloud<pcl::PointXYZI> t_glidarframe::readKittiTxtData(std::string &filepath)
	{
		pcl::PointCloud<pcl::PointXYZI> pointCloudPtrTmp;
		// load point cloud
		ifstream input;
		input.open(filepath);
		if (!input.good())
		{
			std::cout << "Could not read file: " << filepath << std::endl;
			exit(EXIT_FAILURE);
		}
		while (!input.eof()) {
			string line;
			getline(input, line);
			if (line == "")
				continue;
			vector<std::string> arr = split(line, " ");
			pcl::PointXYZI mPoint;
			mPoint.x = stod(arr[0]);
			mPoint.y = stod(arr[1]);
			mPoint.z = stod(arr[2]);
			pointCloudPtrTmp.points.push_back(mPoint);
		}

		return pointCloudPtrTmp;
	}

	//read lidar file (format:bin)
	pcl::PointCloud<pcl::PointXYZI> t_glidarframe::readKittiBinData(std::string &in_file)
	{
		//load point cloud
		std::fstream input(in_file.c_str(), std::ios::in | std::ios::binary);
		if (!input.good()) {
			std::cout << "Could not read file: " << in_file << std::endl;
			exit(EXIT_FAILURE);
		}
		input.seekg(0, std::ios::beg);


		int lidar_line_max = 64;

		pcl::PointCloud<pcl::PointXYZI> pointCloudPtrTmp;
		for (int i = 0; input.good() && !input.eof(); i++)
		{
			pcl::PointXYZI mPoint;
			input.read((char *)&mPoint.x, sizeof(float));
			input.read((char *)&mPoint.y, sizeof(float));
			input.read((char *)&mPoint.z, sizeof(float));
			input.read((char *)&mPoint.intensity, sizeof(float));
			pointCloudPtrTmp.points.push_back(mPoint);
		}
		input.close();

		return pointCloudPtrTmp;
	}

	LidarFrame glidar::t_glidarframe::ProcessBatch(const double& t)
	{
		cout << "==================================" << endl;
		cout << "start process lidar" << endl;
		clear();

		pcl::PointCloud<pcl::PointXYZI> lidarPoints = readKittiBinData(cur_lidar_path.lidar_path);///< get pointCloud from bin
		handleCloud(lidarPoints);

		LidarFrame frame =LidarFrame();

		std::vector<int> indices;
		frame.fullCloud.clear();
		pcl::removeNaNFromPointCloud(_laserCloud, frame.fullCloud, indices);
		frame.LessSharp.clear();
		pcl::removeNaNFromPointCloud(_cornerPointsLessSharp, frame.LessSharp, indices);
		frame.LessSurf.clear();
		pcl::removeNaNFromPointCloud(_surfacePointsLessFlat, frame.LessSurf, indices);
		frame.Sharp.clear();
		pcl::removeNaNFromPointCloud(_cornerPointsSharp, frame.Sharp, indices);
		frame.Surf.clear();
		pcl::removeNaNFromPointCloud(_surfacePointsFlat, frame.Surf, indices);

		return frame;
	}

	void t_glidarframe::process(pcl::PointCloud<pcl::PointXYZI> &laserCloudIn)
	{
		clear();

		std::vector<int> indices;
		pcl::removeNaNFromPointCloud(laserCloudIn, laserCloudIn, indices);
		removeClosedPointCloud(laserCloudIn, laserCloudIn, close_threshold);

		///< determine scan start and end orientations
		int cloudSize = laserCloudIn.size();
		float startOri = -std::atan2(laserCloudIn[0].y, laserCloudIn[0].x);
		float endOri = -std::atan2(laserCloudIn[cloudSize - 1].y,
			laserCloudIn[cloudSize - 1].x) + 2 * float(M_PI);
		if (endOri - startOri > 3 * M_PI) {
			endOri -= 2 * M_PI;
		}
		else if (endOri - startOri < M_PI) {
			endOri += 2 * M_PI;
		}


		std::vector<int> scanStartInd(n_scans, 0);
		std::vector<int> scanEndInd(n_scans, 0);
		

		bool halfPassed = false;
		int count = cloudSize;
		pcl::PointXYZI point;
		_laserCloudScans.resize(n_scans);

		for (int i = 0; i < cloudSize; i++)
		{
			point.x = laserCloudIn[i].x;
			point.y = laserCloudIn[i].y;
			point.z = laserCloudIn[i].z;

			///< whd+ skip points whoes distance larger than 100
			/*if (point.x * point.x + point.y * point.y + point.z * point.z >10000) {
				count--;
				continue;
			}*/

			///< calculate vertical point angle and scan ID
			float angle = std::atan(point.z / sqrt(point.x * point.x + point.y * point.y)) * 180 / float(M_PI);
			int scanID = 0;
			if (n_scans == 16)
			{
				scanID = int((angle + 15) / 2 + 0.5);
				if (scanID > (n_scans - 1) || scanID < 0)
				{
					count--;
					continue;
				}
			}
			else if (n_scans == 32)
			{
				scanID = int((angle + 92.0 / 3.0) * 3.0 / 4.0);
				if (scanID > (n_scans - 1) || scanID < 0)
				{
					count--;
					continue;
				}
			}
			else if (n_scans == 64)
			{
				if (angle >= -8.83)
					scanID = int((2 - angle) * 3.0 + 0.5);
				else
					scanID = n_scans / 2 + int((-8.83 - angle) * 2.0 + 0.5);

				if (angle > 2 || angle < -24.33 || scanID > 50 || scanID < 0)
				{
					count--;
					continue;
				}
			}
			else
			{
				printf("wrong scan number\n");
				break;
			}

			///< calculate horizontal point angle
			float ori = -std::atan2(point.y, point.x);
			if (!halfPassed)
			{
				if (ori < startOri - M_PI / 2)
				{
					ori += 2 * M_PI;
				}
				else if (ori > startOri + M_PI * 3 / 2)
				{
					ori -= 2 * M_PI;
				}

				if (ori - startOri > M_PI)
				{
					halfPassed = true;
				}
			}
			else
			{
				ori += 2 * M_PI;

				if (ori < endOri - M_PI * 3 / 2)
				{
					ori += 2 * M_PI;
				}
				else if (ori > endOri + M_PI / 2)
				{
					ori -= 2 * M_PI;
				}
			}

			///< calculate relative scan time based on point orientation
			float relTime = (ori - startOri) / (endOri - startOri);
			point.intensity = scanID + scanPeriod * relTime;

			///< projectPointToStartOfSweep(point, relTime);//没IMU，point不变

			_laserCloudScans[scanID].push_back(point);
		}
		cloudSize = count;
		printf("points size %d \n", cloudSize);

		pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloud(new pcl::PointCloud<pcl::PointXYZI>());
		for (int i = 0; i < n_scans; i++)
		{
			scanStartInd[i] = laserCloud->size() + 5;
			*laserCloud += _laserCloudScans[i];
			scanEndInd[i] = laserCloud->size() - 6;
		}
		_laserCloud = *laserCloud;

		///< extract feature
		for (int i = 5; i < cloudSize - 5; i++)
		{
			float diffX = laserCloud->points[i - 5].x + laserCloud->points[i - 4].x + laserCloud->points[i - 3].x + laserCloud->points[i - 2].x + laserCloud->points[i - 1].x - 10 * laserCloud->points[i].x + laserCloud->points[i + 1].x + laserCloud->points[i + 2].x + laserCloud->points[i + 3].x + laserCloud->points[i + 4].x + laserCloud->points[i + 5].x;
			float diffY = laserCloud->points[i - 5].y + laserCloud->points[i - 4].y + laserCloud->points[i - 3].y + laserCloud->points[i - 2].y + laserCloud->points[i - 1].y - 10 * laserCloud->points[i].y + laserCloud->points[i + 1].y + laserCloud->points[i + 2].y + laserCloud->points[i + 3].y + laserCloud->points[i + 4].y + laserCloud->points[i + 5].y;
			float diffZ = laserCloud->points[i - 5].z + laserCloud->points[i - 4].z + laserCloud->points[i - 3].z + laserCloud->points[i - 2].z + laserCloud->points[i - 1].z - 10 * laserCloud->points[i].z + laserCloud->points[i + 1].z + laserCloud->points[i + 2].z + laserCloud->points[i + 3].z + laserCloud->points[i + 4].z + laserCloud->points[i + 5].z;

			cloudCurvature[i] = diffX * diffX + diffY * diffY + diffZ * diffZ;
			cloudSortInd[i] = i;
			cloudNeighborPicked[i] = 0;
			cloudLabel[i] = 0;
		}


		for (int i = 0; i < n_scans; i++)
		{
			if (scanEndInd[i] - scanStartInd[i] < 6)
				continue;
			pcl::PointCloud<pcl::PointXYZI>::Ptr surfacePointsLessFlatScan(new pcl::PointCloud<pcl::PointXYZI>);
			for (int j = 0; j < 6; j++)
			{

				int sp = scanStartInd[i] + (scanEndInd[i] - scanStartInd[i]) * j / 6;
				int ep = scanStartInd[i] + (scanEndInd[i] - scanStartInd[i]) * (j + 1) / 6 - 1;

				std::sort(cloudSortInd + sp, cloudSortInd + ep + 1, comp);
				
				int largestPickedNum = 0;
				for (int k = ep; k >= sp; k--)
				{
					int ind = cloudSortInd[k];

					if (cloudNeighborPicked[ind] == 0 &&
						cloudCurvature[ind] > 0.1)
					{
						largestPickedNum++;
						if (largestPickedNum <= 2)
						{
							cloudLabel[ind] = 2;
							_cornerPointsSharp.push_back(laserCloud->points[ind]);
							_cornerPointsLessSharp.push_back(laserCloud->points[ind]);
						}
						else if (largestPickedNum <= 20)
						{
							cloudLabel[ind] = 1;
							_cornerPointsLessSharp.push_back(laserCloud->points[ind]);
						}
						else
						{
							break;
						}

						cloudNeighborPicked[ind] = 1;

						for (int l = 1; l <= 5; l++)
						{
							float diffX = laserCloud->points[ind + l].x - laserCloud->points[ind + l - 1].x;
							float diffY = laserCloud->points[ind + l].y - laserCloud->points[ind + l - 1].y;
							float diffZ = laserCloud->points[ind + l].z - laserCloud->points[ind + l - 1].z;
							if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05)
							{
								break;
							}

							cloudNeighborPicked[ind + l] = 1;
						}
						for (int l = -1; l >= -5; l--)
						{
							float diffX = laserCloud->points[ind + l].x - laserCloud->points[ind + l + 1].x;
							float diffY = laserCloud->points[ind + l].y - laserCloud->points[ind + l + 1].y;
							float diffZ = laserCloud->points[ind + l].z - laserCloud->points[ind + l + 1].z;
							if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05)
							{
								break;
							}

							cloudNeighborPicked[ind - l] = 1;
						}
					}
				}

				int smallestPickedNum = 0;
				for (int k = sp; k <= ep; k++)
				{
					int ind = cloudSortInd[k];

					if (cloudNeighborPicked[ind] == 0 &&
						cloudCurvature[ind] < 0.1)
					{

						cloudLabel[ind] = -1;
						_surfacePointsFlat.push_back(laserCloud->points[ind]);

						smallestPickedNum++;
						if (smallestPickedNum >= 4)
						{
							break;
						}

						cloudNeighborPicked[ind] = 1;
						for (int l = 1; l <= 5; l++)
						{
							float diffX = laserCloud->points[ind + l].x - laserCloud->points[ind + l - 1].x;
							float diffY = laserCloud->points[ind + l].y - laserCloud->points[ind + l - 1].y;
							float diffZ = laserCloud->points[ind + l].z - laserCloud->points[ind + l - 1].z;
							if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05)
							{
								break;
							}

							cloudNeighborPicked[ind + l] = 1;
						}
						for (int l = -1; l >= -5; l--)
						{
							float diffX = laserCloud->points[ind + l].x - laserCloud->points[ind + l + 1].x;
							float diffY = laserCloud->points[ind + l].y - laserCloud->points[ind + l + 1].y;
							float diffZ = laserCloud->points[ind + l].z - laserCloud->points[ind + l + 1].z;
							if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05)
							{
								break;
							}

							cloudNeighborPicked[ind - l] = 1;
						}
					}
				}

				for (int k = sp; k <= ep; k++)
				{
					if (cloudLabel[k] <= 0)
					{
						surfacePointsLessFlatScan->push_back(laserCloud->points[k]);
					}
				}
			}

			pcl::PointCloud<pcl::PointXYZI> surfPointsLessFlatScanDS;
			pcl::VoxelGrid<pcl::PointXYZI> downSizeFilter;
			downSizeFilter.setInputCloud(surfacePointsLessFlatScan);
			downSizeFilter.setLeafSize(0.2, 0.2, 0.2);
			downSizeFilter.filter(surfPointsLessFlatScanDS);

			_surfacePointsLessFlat += surfPointsLessFlatScanDS;
		}


	}

}