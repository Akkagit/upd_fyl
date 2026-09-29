
#include "gpublish.h"


glidar::t_gpublish::t_gpublish()
{
	
}

glidar::t_gpublish::~t_gpublish()
{
	if (viewer != nullptr) delete viewer;
}

void glidar::t_gpublish::Initialize(string path)
{
	outputFile.open(path, ios::trunc);
	viewer = new gviewer();
	viewer->Show();
}

void glidar::t_gpublish::UpdateNewState(const IMUState &imu_state, const LidarStateServer &lidar_states, Eigen::Vector3d init_lidarpos, vector<Eigen::Vector3d> lidarPointcloud)
{

	if (firstFlag_imu)
	{
		if (imu_state.position.norm() >= 0.0)
		{
			_init_imupos = imu_state.position;
			Eigen::Vector3d BLH = XYZ2BLH(_init_imupos);
			R_e_n = R_ENU_ECEF(BLH);
			_init_lidarpos = init_lidarpos;
			firstFlag_imu = false;
		}
	}
	if (firstFlag_lidar && init_lidarpos.norm() > 0 && lidar_states.size() > 0)
	{
		firstFlag_lidar = false;
		_init_lidarpos = init_lidarpos;
	}

	viewer->SetFrames(vector<pair<Eigen::Matrix3d, Eigen::Vector3d>>(0));
	vector<std::pair<Eigen::Matrix3d, Eigen::Vector3d>> frames;
	Eigen::Vector3d Position;
	Eigen::Matrix3d atti;
	VPointCloud lpc;

	Position = R_e_n * (imu_state.position - _init_imupos);
	atti = imu_state.orientation.toRotationMatrix();
	frames.push_back(make_pair(atti, Position));
	if (lidar_states.size() > 0)
	{
		for (auto & lidar_state : lidar_states)
		{
			Eigen::Vector3d pos_lidar;
			Eigen::Matrix3d att_lidar;
			pos_lidar = R_e_n * (lidar_state.second.position + _init_lidarpos - _init_imupos);
			att_lidar = R_e_n * lidar_state.second.orientation.toRotationMatrix();
			frames.push_back(make_pair(att_lidar, pos_lidar));
		}

	}

	for (auto & pt : lidarPointcloud)
	{
		pt= R_e_n * (pt - _init_imupos);
	}

	//att	
	viewer->SetFrames(frames);

	//trajectory
	viewer->AddNewPos(Position);

	//pointcloud
	viewer->SetPointCloud(lidarPointcloud);

	Eigen::Vector3d att = imu_state.orientation.toRotationMatrix().eulerAngles(0,1,2);
	outputFile << setw(20) << setprecision(15) << imu_state.time << ","
		<< setw(15) << setprecision(10) << imu_state.position(0) << ","
		<< setw(15) << setprecision(10) << imu_state.position(1) << ","
		<< setw(15) << setprecision(10) << imu_state.position(2) << ","
		<< setw(15) << setprecision(6) << imu_state.orientation.x() << " ,"
		<< setw(15) << setprecision(6) << imu_state.orientation.y() << " ,"
		<< setw(15) << setprecision(6) << imu_state.orientation.z() << " ,"
		<< setw(15) << setprecision(6) << imu_state.orientation.w() << " ,"
		<< setw(15) << setprecision(8) << att(0)*180.0 / M_PI << ","
		<< setw(15) << setprecision(8) << att(1)*180.0 / M_PI << " ,"
		<< setw(15) << setprecision(8) << att(2)*180.0 / M_PI <<endl;

}

