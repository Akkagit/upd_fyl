
#include "glidar.h"
#include "gset/gsetlidar.h"
#include <sstream>
#include <fstream>

using namespace glidar;
using namespace Eigen;


glidar::t_glidar::t_glidar(t_gsetbase* set):imu_state_id(0),next_id(0)
{
	max_lidarstate_size=dynamic_cast<t_gsetlidar*>(set)->max_lidar_state_size();
	observation_noise=dynamic_cast<t_gsetlidar*>(set)->lidar_observation_noise();
	observation_noise*=observation_noise;
	
	R_lidar_imu= dynamic_cast<t_gsetlidar*>(set)->R_lidar_imu();
	t_lidar_imu = dynamic_cast<t_gsetlidar*>(set)->t_lidar_imu();
	T_lidar_imu = dynamic_cast<t_gsetlidar*>(set)->T_lidar_imu();
	use_scan = dynamic_cast<t_gsetlidar*>(set)->use_scan();
	use_map = dynamic_cast<t_gsetlidar*>(set)->use_map();
	use_corrdistort = dynamic_cast<t_gsetlidar*>(set)->distortion();
	map_gap = dynamic_cast<t_gsetlidar*>(set)->map_gap();
	proj_lastmap = dynamic_cast<t_gsetlidar*>(set)->proj_lastmap();
	
	output_path = dynamic_cast<t_gsetlidar*>(set)->output_path();

	lidarproc = new t_glidarframe(set);

	R_w_e = Eigen::Matrix3d::Identity();
	t_w_e = Eigen::Vector3d::Zero();
	
	P_IL = Eigen::MatrixXd::Zero(15, 0);
	P_LL = Eigen::MatrixXd::Zero(0, 0);
	P_IL_sav = Eigen::MatrixXd::Zero(15, 0);
	P_LL_sav= Eigen::MatrixXd::Zero(0, 0);
    Dx_sav = Eigen::VectorXd::Zero(0);
	Dx_row = 0;
	Phi = Eigen::MatrixXd::Zero(15, 15);

	publisher.Initialize(output_path);
}

glidar::t_glidar::~t_glidar()
{
	if(lidarproc!=nullptr)
		delete lidarproc;
}

/**
* @brief screen out the observations whose residual is larger than n times of the MSE
* @para[in]	H	jacobian matrix of the lidar residual function
* @para[in] r	residual vector of the lidar residual function
* @para[in] n	to identify the dimension of observations
*/
void glidar::t_glidar::screenOut(MatrixXd &H, VectorXd &r,int n)
{
	MatrixXd sort_H;
	VectorXd sort_r;
	VectorXd new_r;
	if (n == 1)
	{
		new_r = r;
	}
	else if (n > 1)
	{
		new_r.resize(r.size() / n);
		for (int i = 0; i < new_r.size(); i++)
		{
			double temp = 0;
			for (int j = 0; j < n; j++)
			{
				temp += pow(r(i*n + j), 2);
			}
			new_r(i) = sqrt(temp);
		}
	}
	else
	{
		cout << "Error:Wrong Observation Residual Number!" << endl;
		return;
	}
	///< calculate the MSE of residuals
	double ave=new_r.sum()/new_r.size();
	double rmse=0;
	for (int i = 0; i < new_r.size(); i++)
	{
		rmse += pow((new_r(i) - ave),2);
	}
	rmse=rmse / new_r.size();

	///< use MSE as threshold to chose H and r
	VectorXi ind = VectorXi::LinSpaced(new_r.size(), 0, new_r.size() - 1);///< [0 1 2 3 ... N-1]
	auto rule = [new_r](int i, int j)->bool {
		return new_r(i) > new_r(j);
	};///< regular expression, used for sort
	std::sort(ind.data(), ind.data() + ind.size(), rule);

	sort_r.resize(r.size());
	sort_H.resize(H.rows(),H.cols());
	int id = 0;
	for (int i = 0; i < new_r.size(); i++) {
		if (pow((new_r(i) - ave), 2) < rmse*2 )
			id++;
		for (int j = 0; j < n; j++)
		{
			sort_r(i*n+j) = r(ind(i)*n+j);
			sort_H.row(i*n + j) = H.row(ind(i)*n + j);
		}
	}
	sort_r.conservativeResize(id*n);
	sort_H.conservativeResize(id*n, H.cols());

	H = sort_H;
	r = sort_r;
}

void glidar::t_glidar::Huber(MatrixXd &H,VectorXd &r, float ther)
{
	///< temporarily useless
	double tmp = ther;
	r = tmp * r;
	H = tmp * H;

	/*if (r.norm() < ther)
	{
		return;
	}
	else if (r.norm() > ther)
	{
		double tmp = ther+2*ther*r.norm() /0.1;
		r = tmp * r;
		H = tmp * H;
	}
	else
	{
		double tmp = sqrt(2 * r.norm()*ther - ther * ther) / r.norm();
		r = tmp*r;
		H = tmp * H;
	}*/
}

void glidar::t_glidar::lidarMeasurementJacobian(lidarOdometryObs obs, Eigen::MatrixXd &H, Eigen::VectorXd &r,bool use_3d,float ther)
{
	cout << "lidarodometry" << endl;
	cout << "corner points:" << obs.cornerFeature.currCornerPointCloud.size() << endl;
	cout << "surf points:" << obs.surfFeature.currSurfPointCloud.size() << endl;
	if (lidar_states.size() != 3)
	{
		cout << "there may be some problem" << endl;
		getchar();
	}

	int corner_obsnum= obs.cornerFeature.currCornerPointCloud.size();
	int surf_obsnum = obs.surfFeature.currSurfPointCloud.size();
	int state_num = lidar_states.size() * 6;
	if (use_3d)
		corner_obsnum = corner_obsnum *3;

	Eigen::MatrixXd corner_H = MatrixXd::Zero(corner_obsnum, state_num);
	Eigen::VectorXd corner_r = VectorXd::Zero(corner_obsnum);
	Eigen::MatrixXd surf_H = MatrixXd::Zero(surf_obsnum, state_num);
	Eigen::VectorXd surf_r = VectorXd::Zero(surf_obsnum);

	for (int i = 0; i < obs.cornerFeature.currCornerPointCloud.size(); i++)
	{
		if (use_3d)
		{
			///< initialize the H matrix
			Eigen::MatrixXd H_dp = MatrixXd::Zero(3, 3);
			Eigen::MatrixXd H_px = MatrixXd::Zero(3, 18);
			///< initialize the r matrix
			Eigen::VectorXd rl = VectorXd::Zero(3);

			///< calculate the residual of observation(distance between point and line)
			Eigen::Vector3d v_lj = obs.cornerFeature.correspondCornerPointCloudA[i] - obs.cornerFeature.correspondCornerPointCloudB[i];
			Eigen::Vector3d v_li = obs.cornerFeature.currCornerPointCloud_inlast[i] - obs.cornerFeature.correspondCornerPointCloudB[i];
			Eigen::Vector3d v_ji = obs.cornerFeature.currCornerPointCloud_inlast[i] - obs.cornerFeature.correspondCornerPointCloudA[i];
			Eigen::Vector3d vjl_i = v_ji.cross(v_li);
			Eigen::Vector3d vl_ji = v_lj.cross(v_li);
			double length_lj = sqrt(v_lj(0)*v_lj(0) + v_lj(1)*v_lj(1) + v_lj(2)*v_lj(2));;

			H_dp(0, 1) = -v_lj(2) / length_lj;
			H_dp(0, 2) = v_lj(1) / length_lj;
			H_dp(1, 0) = v_lj(2) / length_lj;
			H_dp(1, 2) = -v_lj(0) / length_lj;
			H_dp(2, 0) = -v_lj(1) / length_lj;
			H_dp(2, 1) = v_lj(0) / length_lj;


			rl(0) = vl_ji.x() / length_lj;
			rl(1) = vl_ji.y() / length_lj;
			rl(2) = vl_ji.z() / length_lj;

			Eigen::Matrix3d Rk = obs.curr_R_l_e;
			Eigen::Matrix3d Rk_1 = obs.last_R_l_e;
			Eigen::Vector3d pk = obs.curr_t_l_e;
			Eigen::Vector3d pk_1 = obs.last_t_l_e;

			//if (first_odo)
			//{
			//	H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
			//}
			//else
			//{
			//	/*H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();*/
			//	H_px.block<3, 3>(0, 0) = skewSymmetric(Rk_1.transpose()*(Rk*obs.cornerFeature.currCornerPointCloud[i] + pk - pk_1));
			//	H_px.block<3, 3>(0, 3) = -Rk_1.transpose();
			//}

			H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
			if (first_odo)
			{
				H_px.block<3, 3>(0, 6) = Eigen::Matrix3d::Zero();
				H_px.block<3, 3>(0, 9) = Eigen::Matrix3d::Zero();
			}
			else
			{
				H_px.block<3, 3>(0, 6) = skewSymmetric(Rk_1.transpose()*(Rk*obs.cornerFeature.currCornerPointCloud[i] + pk - pk_1));
				H_px.block<3, 3>(0, 9) = -Rk_1.transpose();
			}
			H_px.block<3, 3>(0, 12) = -Rk_1.transpose()*Rk*skewSymmetric(obs.cornerFeature.currCornerPointCloud[i]);
			H_px.block<3, 3>(0, 15) = Rk_1.transpose();

			//Huber(H_px, rl, ther);

			corner_H.block<3, 18>(i * 3, 0) = H_dp * H_px;
			corner_r(i * 3) = rl(0);
			corner_r(i * 3 + 1) = rl(1);
			corner_r(i * 3 + 2) = rl(2);
		}
		else
		{
			///< initialize the H matrix
			Eigen::MatrixXd H_dp = MatrixXd::Zero(1, 3);
			Eigen::MatrixXd H_px = MatrixXd::Zero(3, 18);
			///< initialize the r matrix
			VectorXd rl = VectorXd::Zero(1);

			///< calculate the residual of observation(distance between the point and line)
			Vector3d v_lj = obs.cornerFeature.correspondCornerPointCloudA[i] - obs.cornerFeature.correspondCornerPointCloudB[i];
			Vector3d v_li = obs.cornerFeature.currCornerPointCloud_inlast[i] - obs.cornerFeature.correspondCornerPointCloudB[i];
			Vector3d v_ljxi = v_lj.cross(v_li);
			double length_lj = sqrt(v_lj(0)*v_lj(0) + v_lj(1)*v_lj(1) + v_lj(2)*v_lj(2));
			double length_ljxi = sqrt(v_ljxi(0)*v_ljxi(0) + v_ljxi(1)*v_ljxi(1) + v_ljxi(2)*v_ljxi(2));
			double d_e = length_ljxi / length_lj;

			double a = v_lj(1)*v_li(2) - v_lj(2)*v_li(1);
			double b = v_lj(2)*v_li(0) - v_lj(0)*v_li(2);
			double c = v_lj(0)*v_li(1) - v_lj(1)*v_li(0);

			H_dp(0, 0) = -(v_lj(2)*b - v_lj(1)*c) / length_lj / length_ljxi;
			H_dp(0, 1) = -(v_lj(0)*c - v_lj(2)*a) / length_lj / length_ljxi;
			H_dp(0, 2) = -(v_lj(1)*a - v_lj(0)*b) / length_lj / length_ljxi;

			rl(0) = -d_e;

			Eigen::Matrix3d Rk = obs.curr_R_l_e;
			Eigen::Matrix3d Rk_1 = obs.last_R_l_e;
			Eigen::Vector3d pk = obs.curr_t_l_e;
			Eigen::Vector3d pk_1 = obs.last_t_l_e;

			//if (first_odo)
			//{
			//	H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
			//}
			//else
			//{
			//	H_px.block<3, 3>(0, 0) = skewSymmetric(Rk_1.transpose()*(Rk*obs.cornerFeature.currCornerPointCloud[i] + pk - pk_1));
			//	H_px.block<3, 3>(0, 3) = -Rk_1.transpose();
			//	/*H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();*/
			//}

			H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
			if (first_odo)
			{
				H_px.block<3, 3>(0, 6) = Eigen::Matrix3d::Zero();
				H_px.block<3, 3>(0, 9) = Eigen::Matrix3d::Zero();
			}
			else
			{
				H_px.block<3, 3>(0, 6) = skewSymmetric(Rk_1.transpose()*(Rk*obs.cornerFeature.currCornerPointCloud[i] + pk - pk_1));
				H_px.block<3, 3>(0, 9) = -Rk_1.transpose();
			}
			H_px.block<3, 3>(0, 12) = -Rk_1.transpose()*Rk*skewSymmetric(obs.cornerFeature.currCornerPointCloud[i]);
			H_px.block<3, 3>(0, 15) = Rk_1.transpose();

			//Huber(H_px, rl, ther);

			corner_H.block<1, 18>(i, 0) = H_dp * H_px;
			corner_r(i) = rl(0);
			
		}
		
	}
	///< if need to screen out some outpoints
	/*if (use_3d)
		screenOut(corner_H, corner_r, 3);
	else
		screenOut(corner_H, corner_r, 1);*/

	for (int j = 0; j < obs.surfFeature.currSurfPointCloud.size(); j++)
	{

		///< initialize the H matrix
		Eigen::MatrixXd H_dp = MatrixXd::Zero(1, 3);
		Eigen::MatrixXd H_px = MatrixXd::Zero(3, 18);
		///< initialize the r matrix
		VectorXd rl = VectorXd::Zero(1, 1);

		///< calculate residual of observation(distance between point and plane)
		Vector3d v_jl = obs.surfFeature.correspondSurfPointCloudB[j] - obs.surfFeature.correspondSurfPointCloudA[j];
		Vector3d v_jm = obs.surfFeature.correspondSurfPointCloudC[j] - obs.surfFeature.correspondSurfPointCloudA[j];
		Vector3d v_ji = obs.surfFeature.currSurfPointCloud_inlast[j] - obs.surfFeature.correspondSurfPointCloudA[j];
		Vector3d vj_lm = v_jl.cross(v_jm);
		Vector3d up =Vector3d(vj_lm(0)*v_ji(0), vj_lm(1)*v_ji(1), vj_lm(2)*v_ji(2));
		double dj_lm = sqrt(vj_lm(0)*vj_lm(0)+ vj_lm(1)*vj_lm(1)+ vj_lm(2)*vj_lm(2));
		double d_p = sqrt(up(0)*up(0)+ up(1)*up(1)+ up(2)*up(2))/ dj_lm;

		double a = v_jl(1)*v_jm(2) - v_jl(2)*v_jm(1);
		double b = v_jl(2)*v_jm(0) - v_jl(0)*v_jm(2);
		double c = v_jl(0)*v_jm(1) - v_jl(1)*v_jm(0);

		H_dp(0, 0) = -a * a*v_ji(0) / d_p / dj_lm / dj_lm;
		H_dp(0, 1) = -b * b*v_ji(1) / d_p / dj_lm / dj_lm;
		H_dp(0, 2) = -c * c*v_ji(2) / d_p / dj_lm / dj_lm;

		rl(0) = -d_p;

		Eigen::Matrix3d Rk = obs.curr_R_l_e;
		Eigen::Matrix3d Rk_1 = obs.last_R_l_e;
		Eigen::Vector3d pk = obs.curr_t_l_e;
		Eigen::Vector3d pk_1 = obs.last_t_l_e;

		//if (first_odo)
		//{
		//	H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
		//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
		//}
		//else
		//{
		///*	H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
		//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();*/
		//	H_px.block<3, 3>(0, 0) = skewSymmetric(Rk_1.transpose()*(Rk*obs.surfFeature.currSurfPointCloud[j] + pk - pk_1));
		//	H_px.block<3, 3>(0, 3) = -Rk_1.transpose();
		//}

		H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
		H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
		if (first_odo)
		{
			H_px.block<3, 3>(0, 6) = Eigen::Matrix3d::Zero();
			H_px.block<3, 3>(0, 9) = Eigen::Matrix3d::Zero();
		}
		else
		{
			H_px.block<3, 3>(0, 6) = skewSymmetric(Rk_1.transpose()*(Rk*obs.surfFeature.currSurfPointCloud[j] + pk - pk_1));
			H_px.block<3, 3>(0, 9) = -Rk_1.transpose();
		}
		H_px.block<3, 3>(0, 12) = -Rk_1.transpose()*Rk*skewSymmetric(obs.surfFeature.currSurfPointCloud[j]);
		H_px.block<3, 3>(0, 15) = Rk_1.transpose();

		//Huber(H_px, rl, ther);

		surf_H.block<1, 18>(j, 0) = H_dp * H_px;
		surf_r(j) = rl(0);
	}
	///< if need to screen out some outpoints
	//screenOut(surf_H, surf_r, 1);
	first_odo = false;

	H.resize(corner_H.rows() + surf_H.rows(), corner_H.cols());
	r.resize(corner_r.size() + surf_r.size());
	H << corner_H,
		surf_H;
	r << corner_r,
		surf_r;
	
}

void glidar::t_glidar::lidarMeasurementJacobian(lidarMappingObs obs, Eigen::MatrixXd &H, Eigen::VectorXd &r,bool use_3d,float ther)
{

	int corner_obsnum = obs.cornerFeature.currCornerPointCloud.size();
	int surf_obsnum=obs.surfFeature.currSurfPointCloud.size();
	int state_num = lidar_states.size() * 6;

	if (lidar_states.size() != 3)
	{
		cout << "there is something wrong with lidar states" << endl;
	}

	if (use_3d)
		corner_obsnum = corner_obsnum * 3;

	Eigen::MatrixXd corner_H = MatrixXd::Zero(corner_obsnum, state_num);
	Eigen::VectorXd corner_r = VectorXd::Zero(corner_obsnum);
	Eigen::MatrixXd surf_H = MatrixXd::Zero(surf_obsnum, state_num);
	Eigen::VectorXd surf_r = VectorXd::Zero(surf_obsnum);

	Eigen::Matrix3d Rk = obs.curr_R_l_e;
	Eigen::Vector3d pk = obs.curr_t_l_e;
	Eigen::Matrix3d Rk_1 = obs.last_R_l_e;
	Eigen::Vector3d pk_1 = obs.last_t_l_e;
	Eigen::Matrix3d R_tolast = Rk_1.transpose()*Rk;
	Eigen::Vector3d t_tolast = Rk_1.transpose()*(pk - pk_1);


	for (int i = 0; i < obs.cornerFeature.currCornerPointCloud.size(); i++)
	{
		if (use_3d)
		{
			///< initialize the H matrix
			Eigen::MatrixXd H_dp = MatrixXd::Zero(3, 3);
			Eigen::MatrixXd H_px = MatrixXd::Zero(3, 18);
			///< initialize the r matrix
			Eigen::VectorXd rl = VectorXd::Zero(3);

			///< calculate residual of observation(distance between point and line)
			Eigen::Vector3d v_lj = obs.cornerFeature.correspondCornerPointCloudA[i] - obs.cornerFeature.correspondCornerPointCloudB[i];
			Eigen::Vector3d v_li = obs.cornerFeature.currCornerPointCloud_inlast[i] - obs.cornerFeature.correspondCornerPointCloudB[i];
			Eigen::Vector3d v_ji = obs.cornerFeature.currCornerPointCloud_inlast[i] - obs.cornerFeature.correspondCornerPointCloudA[i];
			Eigen::Vector3d vjl_i = v_ji.cross(v_li);
			Eigen::Vector3d vl_ji = v_lj.cross(v_li);
			double length_lj = sqrt(v_lj(0)*v_lj(0) + v_lj(1)*v_lj(1) + v_lj(2)*v_lj(2));;

			H_dp(0, 1) = -v_lj(2) / length_lj;
			H_dp(0, 2) = v_lj(1) / length_lj;
			H_dp(1, 0) = v_lj(2) / length_lj;
			H_dp(1, 2) = -v_lj(0) / length_lj;
			H_dp(2, 0) = -v_lj(1) / length_lj;
			H_dp(2, 1) = v_lj(0) / length_lj;

			rl(0) = vl_ji.x() / length_lj;
			rl(1) = vl_ji.y() / length_lj;
			rl(2) = vl_ji.z() / length_lj;

			//Eigen::Vector3d p_last = R_tolast*obs.cornerFeature.currCornerPointCloud[i]+t_tolast;
			//if (first_map)
			//{
			//	H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
			//}
			//else
			//{
			//	H_px.block<3, 3>(0, 0) = -R_w_e.transpose()*Rk_1*skewSymmetric(p_last);
			//	H_px.block<3, 3>(0, 3) = R_w_e.transpose();
			//	/*H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();*/
			//}
			if (first_map)
			{
				H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
				H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
			}
			else
			{
				H_px.block<3, 3>(0, 0) = skewSymmetric(R_w_e.transpose()*Rk*obs.cornerFeature.currCornerPointCloud[i]) + skewSymmetric(R_w_e.transpose()*(pk - t_w_e));
				H_px.block<3, 3>(0, 3) = -R_w_e.transpose();
			}
			H_px.block<3, 3>(0, 6) = Eigen::Matrix3d::Zero();
			H_px.block<3, 3>(0, 9) = Eigen::Matrix3d::Zero();
			H_px.block<3, 3>(0, 12) = -R_w_e.transpose()*Rk*skewSymmetric(obs.cornerFeature.currCornerPointCloud[i]);
			H_px.block<3, 3>(0, 15) = R_w_e.transpose();

			//Huber(H_px, rl, ther);

			corner_H.block<3, 18>(i * 3, 0) = H_dp * H_px;
			corner_r(i * 3) = rl(0);
			corner_r(i * 3 + 1) = rl(1);
			corner_r(i * 3 + 2) = rl(2);
			//to be continued
		}
		else
		{
			///< initialize the H matrix
			Eigen::MatrixXd H_dp = MatrixXd::Zero(1, 3);
			Eigen::MatrixXd H_px = MatrixXd::Zero(3, 18);
			///< initialize the r matrix
			Eigen::VectorXd rl = VectorXd::Zero(1);

			///< calculate the residual of observation(distance between point and line)
			Eigen::Vector3d v_lj = obs.cornerFeature.correspondCornerPointCloudA[i] - obs.cornerFeature.correspondCornerPointCloudB[i];
			Eigen::Vector3d v_li = obs.cornerFeature.currCornerPointCloud_inlast[i] - obs.cornerFeature.correspondCornerPointCloudB[i];
			Eigen::Vector3d v_ljxi = v_lj.cross(v_li);
			double length_lj = sqrt(v_lj(0)*v_lj(0) + v_lj(1)*v_lj(1) + v_lj(2)*v_lj(2));
			double length_ljxi = sqrt(v_ljxi(0)*v_ljxi(0) + v_ljxi(1)*v_ljxi(1) + v_ljxi(2)*v_ljxi(2));
			double d_e = length_ljxi / length_lj;

			double a = v_lj(1)*v_li(2) - v_lj(2)*v_li(1);
			double b = v_lj(2)*v_li(0) - v_lj(0)*v_li(2);
			double c = v_lj(0)*v_li(1) - v_lj(1)*v_li(0);

			H_dp(0, 0) = -(v_lj(2)*b - v_lj(1)*c) / length_lj / length_ljxi;
			H_dp(0, 1) = -(v_lj(0)*c - v_lj(2)*a) / length_lj / length_ljxi;
			H_dp(0, 2) = -(v_lj(1)*a - v_lj(0)*b) / length_lj / length_ljxi;

			rl(0) = -d_e;

			//Eigen::Vector3d p_last = R_tolast * obs.cornerFeature.currCornerPointCloud[i] + t_tolast;
			//if (first_map)
			//{
			//	H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
			//}
			//else
			//{
			//	/*H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();*/
			//	H_px.block<3, 3>(0, 0) = -R_w_e.transpose()*Rk_1*skewSymmetric(p_last);
			//	H_px.block<3, 3>(0, 3) = R_w_e.transpose();
			//}
			if (first_map)
			{
				H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
				H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
			}
			else
			{
				H_px.block<3, 3>(0, 0) = skewSymmetric(R_w_e.transpose()*Rk*obs.cornerFeature.currCornerPointCloud[i]) + skewSymmetric(R_w_e.transpose()*(pk - t_w_e));
				H_px.block<3, 3>(0, 3) = -R_w_e.transpose();
			}
			H_px.block<3, 3>(0, 6) = Eigen::Matrix3d::Zero();
			H_px.block<3, 3>(0, 9) = Eigen::Matrix3d::Zero();
			H_px.block<3, 3>(0, 12) = -R_w_e.transpose()*Rk*skewSymmetric(obs.cornerFeature.currCornerPointCloud[i]);
			H_px.block<3, 3>(0, 15) = R_w_e.transpose();

			//Huber(H_px, rl, ther);

			corner_H.block<1, 18>(i, 0) = H_dp * H_px;
			corner_r(i) = rl(0);

		}
	}
	///< if need to screen out some outpoints
	/*if (use_3d)
		screenOut(corner_H, corner_r, 3);
	else
		screenOut(corner_H, corner_r, 1);*/

	for (int j = 0; j < obs.surfFeature.currSurfPointCloud.size(); j++)
	{
		///< initialize the H matrix
		Eigen::MatrixXd H_dp = MatrixXd::Zero(1, 3);
		Eigen::MatrixXd H_px = MatrixXd::Zero(3, 18);
		///< initialize the r matrix
		Eigen::VectorXd rl = VectorXd::Zero(1);
		///< calculate residual of observation(distance between point and plane)
		rl(0) = -(obs.surfFeature.norm[j].dot(obs.surfFeature.currSurfPointCloud_inlast[j]) + obs.surfFeature.negative_OA_dot_norm[j]);

		H_dp(0, 0) = -obs.surfFeature.norm[j](0);
		H_dp(0, 1) = -obs.surfFeature.norm[j](1);
		H_dp(0, 2) = -obs.surfFeature.norm[j](2);

		//Eigen::Vector3d p_last = R_tolast * obs.surfFeature.currSurfPointCloud[j] + t_tolast;
		//if (first_map)
		//{
		//	H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
		//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
		//}
		//else
		//{
		//	H_px.block<3, 3>(0, 0) = -R_w_e.transpose()*Rk_1*skewSymmetric(p_last);
		//	H_px.block<3, 3>(0, 3) = R_w_e.transpose();
		//	/*H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
		//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();*/
		//}
		if (first_map)
		{
			H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
		}
		else
		{
			H_px.block<3, 3>(0, 0) = skewSymmetric(R_w_e.transpose()*Rk*obs.surfFeature.currSurfPointCloud[j]) + skewSymmetric(R_w_e.transpose()*(pk - t_w_e));
			H_px.block<3, 3>(0, 3) = -R_w_e.transpose();
		}
		H_px.block<3, 3>(0, 6) = Eigen::Matrix3d::Zero();
		H_px.block<3, 3>(0, 9) = Eigen::Matrix3d::Zero();
		H_px.block<3, 3>(0, 12) = -R_w_e.transpose()*Rk*skewSymmetric(obs.surfFeature.currSurfPointCloud[j]);
		H_px.block<3, 3>(0, 15) = R_w_e.transpose();

		//Huber(H_px, rl, ther);

		surf_H.block<1, 18>(j, 0) = H_dp * H_px;
		surf_r(j) = rl(0);
	}
	///< if need to screen out some outpoints
	//screenOut(surf_H, surf_r, 1);

	first_map = false;

	H.resize(surf_H.rows() + corner_H.rows(), surf_H.cols());
	r.resize(surf_r.size() + corner_r.size());

	H << surf_H,
		corner_H;
	r << surf_r,
		corner_r;

}


void glidar::t_glidar::lidarMeasurementJacobian_test(lidarMappingObs obs, Eigen::MatrixXd &H, Eigen::VectorXd &r, bool use_3d, float ther)
{

	int corner_obsnum = obs.cornerFeature.currCornerPointCloud.size();
	int surf_obsnum = obs.surfFeature.currSurfPointCloud.size();
	int state_num = lidar_states.size() * 6;

	if (lidar_states.size() != 3)
	{
		cout << "there is something wrong with lidar states" << endl;
	}

	if (use_3d)
		corner_obsnum = corner_obsnum * 3;

	Eigen::MatrixXd corner_H = MatrixXd::Zero(corner_obsnum, state_num);
	Eigen::VectorXd corner_r = VectorXd::Zero(corner_obsnum);
	Eigen::MatrixXd surf_H = MatrixXd::Zero(surf_obsnum, state_num);
	Eigen::VectorXd surf_r = VectorXd::Zero(surf_obsnum);

	Eigen::Matrix3d Rk = obs.curr_R_l_e;
	Eigen::Vector3d pk = obs.curr_t_l_e;
	Eigen::Matrix3d Rk_1 = obs.last_R_l_e;
	Eigen::Vector3d pk_1 = obs.last_t_l_e;
	Eigen::Matrix3d R_tolast = Rk_1.transpose()*Rk;
	Eigen::Vector3d t_tolast = Rk_1.transpose()*(pk - pk_1);


	for (int i = 0; i < obs.cornerFeature.currCornerPointCloud.size(); i++)
	{
		if (use_3d)
		{
			///< initialize the H matrix
			Eigen::MatrixXd H_dp = MatrixXd::Zero(3, 3);
			Eigen::MatrixXd H_px = MatrixXd::Zero(3, 18);
			///< initialize the r matrix
			Eigen::VectorXd rl = VectorXd::Zero(3);

			///< calculate residual of observation(distance between point and line)
			Eigen::Vector3d v_lj = obs.cornerFeature.correspondCornerPointCloudA[i] - obs.cornerFeature.correspondCornerPointCloudB[i];
			Eigen::Vector3d v_li = obs.cornerFeature.currCornerPointCloud_inlast[i] - obs.cornerFeature.correspondCornerPointCloudB[i];
			Eigen::Vector3d v_ji = obs.cornerFeature.currCornerPointCloud_inlast[i] - obs.cornerFeature.correspondCornerPointCloudA[i];
			Eigen::Vector3d vjl_i = v_ji.cross(v_li);
			Eigen::Vector3d vl_ji = v_lj.cross(v_li);
			double length_lj = sqrt(v_lj(0)*v_lj(0) + v_lj(1)*v_lj(1) + v_lj(2)*v_lj(2));;

			H_dp(0, 1) = -v_lj(2) / length_lj;
			H_dp(0, 2) = v_lj(1) / length_lj;
			H_dp(1, 0) = v_lj(2) / length_lj;
			H_dp(1, 2) = -v_lj(0) / length_lj;
			H_dp(2, 0) = -v_lj(1) / length_lj;
			H_dp(2, 1) = v_lj(0) / length_lj;

			rl(0) = vl_ji.x() / length_lj;
			rl(1) = vl_ji.y() / length_lj;
			rl(2) = vl_ji.z() / length_lj;

			//Eigen::Vector3d p_last = R_tolast*obs.cornerFeature.currCornerPointCloud[i]+t_tolast;
			//if (first_map)
			//{
			//	H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
			//}
			//else
			//{
			//	H_px.block<3, 3>(0, 0) = -R_w_e.transpose()*Rk_1*skewSymmetric(p_last);
			//	H_px.block<3, 3>(0, 3) = R_w_e.transpose();
			//	/*H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();*/
			//}
			H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
			if (first_map)
			{
				H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
				H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
			}
			else
			{
				H_px.block<3, 3>(0, 0) = skewSymmetric(Rk_1.transpose()*Rk*obs.cornerFeature.currCornerPointCloud[i]) + skewSymmetric(R_w_e.transpose()*(pk - pk_1));
				H_px.block<3, 3>(0, 3) = -Rk_1.transpose();
			}
			H_px.block<3, 3>(0, 12) = -Rk_1.transpose()*Rk*skewSymmetric(obs.cornerFeature.currCornerPointCloud[i]);
			H_px.block<3, 3>(0, 15) = Rk_1.transpose();

			//Huber(H_px, rl, ther);

			corner_H.block<3, 18>(i * 3, 0) = H_dp * H_px;
			corner_r(i * 3) = rl(0);
			corner_r(i * 3 + 1) = rl(1);
			corner_r(i * 3 + 2) = rl(2);
			//to be continued
		}
		else
		{
			///< initialize the H matrix
			Eigen::MatrixXd H_dp = MatrixXd::Zero(1, 3);
			Eigen::MatrixXd H_px = MatrixXd::Zero(3, 18);
			///< initialize the r matrix
			Eigen::VectorXd rl = VectorXd::Zero(1);

			///< calculate the residual of observation(distance between point and line)
			Eigen::Vector3d v_lj = obs.cornerFeature.correspondCornerPointCloudA[i] - obs.cornerFeature.correspondCornerPointCloudB[i];
			Eigen::Vector3d v_li = obs.cornerFeature.currCornerPointCloud_inlast[i] - obs.cornerFeature.correspondCornerPointCloudB[i];
			Eigen::Vector3d v_ljxi = v_lj.cross(v_li);
			double length_lj = sqrt(v_lj(0)*v_lj(0) + v_lj(1)*v_lj(1) + v_lj(2)*v_lj(2));
			double length_ljxi = sqrt(v_ljxi(0)*v_ljxi(0) + v_ljxi(1)*v_ljxi(1) + v_ljxi(2)*v_ljxi(2));
			double d_e = length_ljxi / length_lj;

			double a = v_lj(1)*v_li(2) - v_lj(2)*v_li(1);
			double b = v_lj(2)*v_li(0) - v_lj(0)*v_li(2);
			double c = v_lj(0)*v_li(1) - v_lj(1)*v_li(0);

			H_dp(0, 0) = -(v_lj(2)*b - v_lj(1)*c) / length_lj / length_ljxi;
			H_dp(0, 1) = -(v_lj(0)*c - v_lj(2)*a) / length_lj / length_ljxi;
			H_dp(0, 2) = -(v_lj(1)*a - v_lj(0)*b) / length_lj / length_ljxi;

			rl(0) = -d_e;

			//Eigen::Vector3d p_last = R_tolast * obs.cornerFeature.currCornerPointCloud[i] + t_tolast;
			//if (first_map)
			//{
			//	H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
			//}
			//else
			//{
			//	/*H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();*/
			//	H_px.block<3, 3>(0, 0) = -R_w_e.transpose()*Rk_1*skewSymmetric(p_last);
			//	H_px.block<3, 3>(0, 3) = R_w_e.transpose();
			//}
			H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
			if (first_map)
			{
				H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
				H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
			}
			else
			{
				H_px.block<3, 3>(0, 0) = skewSymmetric(Rk_1.transpose()*Rk*obs.cornerFeature.currCornerPointCloud[i]) + skewSymmetric(R_w_e.transpose()*(pk - pk_1));
				H_px.block<3, 3>(0, 3) = -Rk_1.transpose();
			}
			H_px.block<3, 3>(0, 12) = -Rk_1.transpose()*Rk*skewSymmetric(obs.cornerFeature.currCornerPointCloud[i]);
			H_px.block<3, 3>(0, 15) = Rk_1.transpose();

			//Huber(H_px, rl, ther);

			corner_H.block<1, 18>(i, 0) = H_dp * H_px;
			corner_r(i) = rl(0);

		}
	}
	///< if need to screen out some outpoints
	/*if (use_3d)
		screenOut(corner_H, corner_r, 3);
	else
		screenOut(corner_H, corner_r, 1);*/

	for (int j = 0; j < obs.surfFeature.currSurfPointCloud.size(); j++)
	{
		///< initialize the H matrix
		Eigen::MatrixXd H_dp = MatrixXd::Zero(1, 3);
		Eigen::MatrixXd H_px = MatrixXd::Zero(3, 18);
		///< initialize the r matrix
		Eigen::VectorXd rl = VectorXd::Zero(1);
		///< calculate residual of observation(distance between point and plane)
		rl(0) = -(obs.surfFeature.norm[j].dot(obs.surfFeature.currSurfPointCloud_inlast[j]) + obs.surfFeature.negative_OA_dot_norm[j]);

		H_dp(0, 0) = -obs.surfFeature.norm[j](0);
		H_dp(0, 1) = -obs.surfFeature.norm[j](1);
		H_dp(0, 2) = -obs.surfFeature.norm[j](2);

		//Eigen::Vector3d p_last = R_tolast * obs.surfFeature.currSurfPointCloud[j] + t_tolast;
		//if (first_map)
		//{
		//	H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
		//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
		//}
		//else
		//{
		//	H_px.block<3, 3>(0, 0) = -R_w_e.transpose()*Rk_1*skewSymmetric(p_last);
		//	H_px.block<3, 3>(0, 3) = R_w_e.transpose();
		//	/*H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
		//	H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();*/
		//}
		H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
		H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
		if (first_map)
		{
			H_px.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
			H_px.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
		}
		else
		{
			H_px.block<3, 3>(0, 0) = skewSymmetric(Rk_1.transpose()*Rk*obs.surfFeature.currSurfPointCloud[j]) + skewSymmetric(R_w_e.transpose()*(pk - pk_1));
			H_px.block<3, 3>(0, 3) = -Rk_1.transpose();
		}
		H_px.block<3, 3>(0, 12) = -Rk_1.transpose()*Rk*skewSymmetric(obs.surfFeature.currSurfPointCloud[j]);
		H_px.block<3, 3>(0, 15) = Rk_1.transpose();

		//Huber(H_px, rl, ther);

		surf_H.block<1, 18>(j, 0) = H_dp * H_px;
		surf_r(j) = rl(0);
	}
	///< if need to screen out some outpoints
	//screenOut(surf_H, surf_r, 1);

	first_map = false;

	H.resize(surf_H.rows() + corner_H.rows(), surf_H.cols());
	r.resize(surf_r.size() + corner_r.size());

	H << surf_H,
		corner_H;
	r << surf_r,
		corner_r;

}
