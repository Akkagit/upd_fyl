
#include "gpublish.h"


gvis::t_gpublish::t_gpublish()
{
}

gvis::t_gpublish::~t_gpublish()
{
	if (viewer != nullptr) delete viewer;
}

void gvis::t_gpublish::Initialize()
{
	
	viewer = new gviewer();
	viewer->Show();
}

void gvis::t_gpublish::UpdateNewState(const IMUState &imu_state, const CamStateServer &cam_states, const Eigen::Vector3d & init_campos)
{
	
	if (firstFlag_imu)
	{
		if (imu_state.position.norm() >= 0.0)
		{
			_init_imupos = imu_state.position;
			Eigen::Vector3d BLH = XYZ2BLH(_init_imupos);
			R_e_n = R_ENU_ECEF(BLH);
			_init_campos = init_campos;
			firstFlag_imu = false;
		}
	}
	if (firstFlag_cam && init_campos.norm() > 0 && cam_states.size()>0)
	{
		firstFlag_cam = false;
		_init_campos = init_campos;
	}

	viewer->SetFrames(vector<pair<Eigen::Matrix3d, Eigen::Vector3d>>(0));
	vector<pair<Eigen::Matrix3d, Eigen::Vector3d>> frames;
	Eigen::Vector3d Position;
	Eigen::Matrix3d atti;

	Position = R_e_n * (imu_state.position-_init_imupos);
	atti = imu_state.orientation.toRotationMatrix();
	frames.push_back(make_pair(atti, Position));
	if (cam_states.size() > 0)
	{
		for (auto & cam_state : cam_states)
		{
			Eigen::Vector3d pos_cam;
			Eigen::Matrix3d att_cam;
			pos_cam = R_e_n * (cam_state.second.position + _init_campos - _init_imupos);
			att_cam = R_e_n * cam_state.second.orientation.toRotationMatrix();
			frames.push_back(make_pair(att_cam, pos_cam));
		}

	}

	//att
	viewer->SetFrames(frames);

	//trajectory
	viewer->AddNewPos(Position);


}

