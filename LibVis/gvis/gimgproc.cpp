
#include "gimgproc.h"
#include "gset/gsetvis.h"
using namespace Eigen;
using namespace gvis;

gvis::t_gimgproc::t_gimgproc(t_gsetbase* _set)
{
	/*  ***Basic Parameters*** */
	ts = dynamic_cast<t_gsetvis*>(_set)->ts();
	freq = dynamic_cast<t_gsetvis*>(_set)->freq();
	imu_ts = dynamic_cast<t_gsetvis*>(_set)->imu_ts();
	imu_freq = dynamic_cast<t_gsetvis*>(_set)->imu_freq();
	cam0_resolution = dynamic_cast<t_gsetvis*>(_set)->cam0_resolution();
	cam1_resolution = dynamic_cast<t_gsetvis*>(_set)->cam1_resolution();
	cam0_distortion_model = dynamic_cast<t_gsetvis*>(_set)->cam0_distortion_model();
	cam1_distortion_model = dynamic_cast<t_gsetvis*>(_set)->cam1_distortion_model();
	cam0_intrinsics = dynamic_cast<t_gsetvis*>(_set)->cam0_intrinsics();
	cam1_intrinsics = dynamic_cast<t_gsetvis*>(_set)->cam1_intrinsics();
	cam0_distortion_coeffs = dynamic_cast<t_gsetvis*>(_set)->cam0_distortion_coeffs();
	cam1_distortion_coeffs = dynamic_cast<t_gsetvis*>(_set)->cam1_distortion_coeffs();

	R_cam0_cam1 = dynamic_cast<t_gsetvis*>(_set)->R_cam0_cam1();
	t_cam0_cam1 = dynamic_cast<t_gsetvis*>(_set)->t_cam0_cam1();
	T_cam0_cam1 = dynamic_cast<t_gsetvis*>(_set)->T_cam0_cam1();

	R_cam0_imu = dynamic_cast<t_gsetvis*>(_set)->R_cam0_imu();
	t_cam0_imu = dynamic_cast<t_gsetvis*>(_set)->t_cam0_imu();
	T_cam0_imu = dynamic_cast<t_gsetvis*>(_set)->T_cam0_imu();

	T_cam1_imu = T_cam0_imu * T_cam0_cam1.inverse();
	R_cam1_imu = T_cam1_imu.linear();
	t_cam1_imu = T_cam1_imu.translation();

	//Stereo Optical flow 这个是不是应该放在下层光流中更加合适呢？
	grid_row = dynamic_cast<t_gsetvis*>(_set)->grid_row();
	grid_col = dynamic_cast<t_gsetvis*>(_set)->grid_col();
	grid_min_feature_num = dynamic_cast<t_gsetvis*>(_set)->grid_min_feature_num();
	grid_max_feature_num = dynamic_cast<t_gsetvis*>(_set)->grid_max_feature_num();
	pyramid_levels = dynamic_cast<t_gsetvis*>(_set)->pyramid_levels();
	patch_size = dynamic_cast<t_gsetvis*>(_set)->patch_size();
	fast_threshold = dynamic_cast<t_gsetvis*>(_set)->fast_threshold();
	ransac_threshold = dynamic_cast<t_gsetvis*>(_set)->ransac_threshold();
	stereo_threshold = dynamic_cast<t_gsetvis*>(_set)->stereo_threshold();
	max_iteration = dynamic_cast<t_gsetvis*>(_set)->max_iteration();
	track_precision = dynamic_cast<t_gsetvis*>(_set)->track_precision();

	/*cout << "grid_row:" << grid_row << endl;
	cout << "grid_col:" << grid_col << endl;
	cout << "grid_min_feature_num:" << grid_min_feature_num << endl;
	cout << "grid_max_feature_num:" << grid_max_feature_num << endl;
	cout << "pyramid_levels:" << pyramid_levels << endl;
	cout << "patch_size:" << patch_size << endl;
	cout << "fast_threshold:" << fast_threshold << endl;
	cout << "ransac_threshold:" << ransac_threshold << endl;
	cout << "stereo_threshold:" << stereo_threshold << endl;
	cout << "max_iteration:" << max_iteration << endl;
	cout << "track_precision:" << track_precision << endl;*/

	//Stereo Orb Optical Flow
	nfeatures = dynamic_cast<t_gsetvis*>(_set)->nfeatures();
	max_fast = dynamic_cast<t_gsetvis*>(_set)->max_fast();
	min_fast = dynamic_cast<t_gsetvis*>(_set)->min_fast();
	max_cnt = dynamic_cast<t_gsetvis*>(_set)->max_cnt();
	fscalefactor = dynamic_cast<t_gsetvis*>(_set)->fscalefactor();
	usingstereorecify = dynamic_cast<t_gsetvis*>(_set)->usingstereorecify();
	K_L = dynamic_cast<t_gsetvis*>(_set)->K_L();
	D_L = dynamic_cast<t_gsetvis*>(_set)->D_L();
	R_L = dynamic_cast<t_gsetvis*>(_set)->R_L();
	P_L = dynamic_cast<t_gsetvis*>(_set)->P_L();

	K_R = dynamic_cast<t_gsetvis*>(_set)->K_R();
	D_R = dynamic_cast<t_gsetvis*>(_set)->D_R();
	R_R = dynamic_cast<t_gsetvis*>(_set)->R_R();
	P_R = dynamic_cast<t_gsetvis*>(_set)->P_R();
	mbf = dynamic_cast<t_gsetvis*>(_set)->mbf();



	//估计参数
	max_cam_state_size = dynamic_cast<t_gsetvis*>(_set)->max_cam_state_size();
	position_std_threshold = dynamic_cast<t_gsetvis*>(_set)->position_std_threshold();
	rotation_threshold = dynamic_cast<t_gsetvis*>(_set)->rotation_threshold();
	translation_threshold = dynamic_cast<t_gsetvis*>(_set)->translation_threshold();
	feature_observation_noise = dynamic_cast<t_gsetvis*>(_set)->feature_observation_noise();

	//跟踪点深度优化相关参数 feature
	huber_epsilon = dynamic_cast<t_gsetvis*>(_set)->huber_epsilon();
	estimation_precision = dynamic_cast<t_gsetvis*>(_set)->estimation_precision();
	initial_damping = dynamic_cast<t_gsetvis*>(_set)->initial_damping();
	outler_loop_max_iteration = dynamic_cast<t_gsetvis*>(_set)->outler_loop_max_iteration();
	inner_loop_max_iteration = dynamic_cast<t_gsetvis*>(_set)->inner_loop_max_iteration();

	estimate_extrinsic = dynamic_cast<t_gsetvis*>(_set)->estimate_extrinsic();
	if (estimate_extrinsic)
	{
		initial_extrinsic_rotation_cov = dynamic_cast<t_gsetvis*>(_set)->initial_extrinsic_rotation_cov();
		initial_extrinsic_translation_cov = dynamic_cast<t_gsetvis*>(_set)->initial_extrinsic_translation_cov();
	}
}


void gvis::t_gimgproc::load_imuobs(const double &t,const vector<Eigen::Vector3d> & wm,const vector<Eigen::Vector3d> & vm)
{
	
	assert(wm.size() == vm.size());
	/* only one sample can be processed successfully */
	assert(wm.size() == 1);
	for (int i = 0; i < wm.size(); i++)
	{
		imumsg tmp_imu;
		tmp_imu.t = t;
		tmp_imu.angular_velocity = wm[i]/imu_ts;
		tmp_imu.linear_acceleration = vm[i]/imu_ts;
		_vecimu.push_back(tmp_imu);
	}
	
}

void gvis::t_gimgproc::load_imgobs(const double &t, const imgPath &img_path)
{
	cur_img_path = img_path;
}

PointCloud gvis::t_gimgproc::ProcessBatch(const double& t)
{
	PointCloud res;
	return res;
}

