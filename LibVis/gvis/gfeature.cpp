
#include "gfeature.h"
#include <numeric>


gvis::t_gfeature::t_gfeature() :id(0), position(Eigen::Vector3d::Zero()),
is_initialized(false)
{
}

gvis::t_gfeature::t_gfeature(const Feature_param &_param) : id(0), position(Eigen::Vector3d::Zero()),
is_initialized(false)
{
	translation_threshold = _param.translation_threshold;
	huber_epsilon = _param.huber_epsilon;
	estimation_precision = _param.estimation_precision;
	initial_damping = _param.initial_damping;
	outler_loop_max_iteration = _param.outler_loop_max_iteration;
	inner_loop_max_iteration = _param.inner_loop_max_iteration;
	T_cam0_cam1 = _param.T_cam0_cam1;
	stereo = _param.stereo;
	

	//translation_threshold = dynamic_cast<t_gsetvis*>(_set)->translation_threshold();
}

gvis::t_gfeature::t_gfeature(const Feature_param &_param,const FeatureIDType& new_id) : id(new_id),
position(Eigen::Vector3d::Zero()),
is_initialized(false), is_KeyFrame(false), is_initialized_NonKey(false)
{
	translation_threshold = _param.translation_threshold;
	huber_epsilon = _param.huber_epsilon;
	estimation_precision = _param.estimation_precision;
	initial_damping = _param.initial_damping;
	outler_loop_max_iteration = _param.outler_loop_max_iteration;
	inner_loop_max_iteration = _param.inner_loop_max_iteration;
	T_cam0_cam1 = _param.T_cam0_cam1;
	stereo = _param.stereo;
}

bool gvis::t_gfeature::triangulatePoint(const CamStateServer& cam_states) const
{

	//把单目的三角化判断加在这里
	const FeatureIDType& first_cam_id = observations.begin()->first;
	const FeatureIDType& last_cam_id = (--observations.end())->first;


	Eigen::Isometry3d first_cam_pose; /// R_c_w, camera to world
	first_cam_pose.linear() =
		cam_states.find(first_cam_id)->second.orientation.toRotationMatrix();
	first_cam_pose.translation() =
		cam_states.find(first_cam_id)->second.position;

	Eigen::Isometry3d last_cam_pose; /// R_c_w, camera to world
	last_cam_pose.linear() =
		cam_states.find(last_cam_id)->second.orientation.toRotationMatrix();
	last_cam_pose.translation() =
		cam_states.find(last_cam_id)->second.position;


	Eigen::Vector3d feature_direction(
		observations.begin()->second(0),
		observations.begin()->second(1), 1.0);
	feature_direction = feature_direction / feature_direction.norm();
	feature_direction = first_cam_pose.linear()*feature_direction;



	Eigen::Vector3d translation = last_cam_pose.translation() -
		first_cam_pose.translation();
	double parallel_translation =
		translation.transpose()*feature_direction;
	Eigen::Vector3d orthogonal_translation = translation -
		parallel_translation * feature_direction;



	if (orthogonal_translation.norm() >
		translation_threshold)
		return true;
	else
		return false;
}

bool gvis::t_gfeature::checkMotion(const CamStateServer& cam_states) const
{
	std::vector<double> x_vector, y_vector;
	for (auto iter = observations.begin(); iter != (observations.end()); iter++)
	{
		auto iter_2 = iter;
		iter_2++;
		if (iter_2 == observations.end())
			break;
		x_vector.emplace_back(iter->second(0) - iter_2->second(1));
		y_vector.emplace_back(iter->second(0) - iter_2->second(1));
	}
	int size_vector = x_vector.size();
	if (size_vector < 3)
	{
		return false;
	}
	double x_average = std::accumulate(x_vector.begin(), x_vector.end(), 0) / size_vector;
	double y_average = std::accumulate(y_vector.begin(), y_vector.end(), 0) / size_vector;
	double x_std = 0.0, y_std = 0.0;
	for (size_t i = 0; i < size_vector; i++)
	{
		x_std += pow(x_vector[i] - x_average, 2) / size_vector;
		y_std += pow(y_vector[i] - y_average, 2) / size_vector;
	}
	x_std = sqrt(x_std);
	y_std = sqrt(y_std);

	if (x_std > 10 || y_std > 8)
	{
		return false;
	}
	for (size_t i = 0; i < size_vector; i++)
	{
		if (x_vector[i] > x_average + (x_std * 3) || x_vector[i] < x_average - (x_std * 3))
		{
			return false;
		}
		if (y_vector[i] > y_average + (y_std * 3) || y_vector[i] < y_average - (y_std * 3))
		{
			return false;
		}
	}
	return true;
}

void gvis::t_gfeature::generateInitialGuess(
	const Eigen::Isometry3d& T_c1_c2, const Eigen::Vector2d& z1,
	const Eigen::Vector2d& z2, Eigen::Vector3d& p) const
{

	Eigen::Vector3d m = T_c1_c2.linear() * Eigen::Vector3d(z1(0), z1(1), 1.0);

	Eigen::Vector2d A(0.0, 0.0);
	A(0) = m(0) - z2(0)*m(2);
	A(1) = m(1) - z2(1)*m(2);

	Eigen::Vector2d b(0.0, 0.0);
	b(0) = z2(0)*T_c1_c2.translation()(2) - T_c1_c2.translation()(0);
	b(1) = z2(1)*T_c1_c2.translation()(2) - T_c1_c2.translation()(1);


	double depth = (A.transpose() * A).inverse() * A.transpose() * b;
	p(0) = z1(0) * depth;
	p(1) = z1(1) * depth;
	p(2) = depth;
	return;
}

void gvis::t_gfeature::cost(const Eigen::Isometry3d& T_c0_ci,
	const Eigen::Vector3d& x, const Eigen::Vector2d& z,
	double& e) const
{

	const double& alpha = x(0);
	const double& beta = x(1);
	const double& rho = x(2);

	Eigen::Vector3d h = T_c0_ci.linear()*
		Eigen::Vector3d(alpha, beta, 1.0) + rho * T_c0_ci.translation();
	double& h1 = h(0);
	double& h2 = h(1);
	double& h3 = h(2);


	Eigen::Vector2d z_hat(h1 / h3, h2 / h3);


	e = (z_hat - z).squaredNorm();
	return;
}

void gvis::t_gfeature::jacobian(const Eigen::Isometry3d& T_c0_ci,
	const Eigen::Vector3d& x, const Eigen::Vector2d& z,
	Eigen::Matrix<double, 2, 3>& J, Eigen::Vector2d& r,
	double& w) const
{


	const double& alpha = x(0);
	const double& beta = x(1);
	const double& rho = x(2);

	Eigen::Vector3d h = T_c0_ci.linear()*
		Eigen::Vector3d(alpha, beta, 1.0) + rho * T_c0_ci.translation();
	double& h1 = h(0);
	double& h2 = h(1);
	double& h3 = h(2);


	Eigen::Matrix3d W;
	W.leftCols<2>() = T_c0_ci.linear().leftCols<2>();
	W.rightCols<1>() = T_c0_ci.translation();

	J.row(0) = 1 / h3 * W.row(0) - h1 / (h3*h3)*W.row(2);
	J.row(1) = 1 / h3 * W.row(1) - h2 / (h3*h3)*W.row(2);


	Eigen::Vector2d z_hat(h1 / h3, h2 / h3);
	r = z_hat - z;


	double e = r.norm();
	if (e <= huber_epsilon)
		w = 1.0;
	else
		w = huber_epsilon / (2 * e);

	return;
}



bool gvis::t_gfeature::initializePosition(const CamStateServer& cam_states)
{	
	// for keyframe ,need to keep
	if (observations.size() < 5) return false;

	std::vector<Eigen::Isometry3d,
		Eigen::aligned_allocator<Eigen::Isometry3d> > cam_poses(0);
	std::vector<Eigen::Vector2d,
		Eigen::aligned_allocator<Eigen::Vector2d> > measurements(0);

	for (auto& m : observations)
	{
		auto cam_state_iter = cam_states.find(m.first);
		if (cam_state_iter == cam_states.end())
			continue;


		measurements.push_back(m.second.head<2>());
		if(stereo)
		measurements.push_back(m.second.tail<2>());


		Eigen::Isometry3d cam0_pose;
		cam0_pose.linear() = cam_state_iter->second.orientation.toRotationMatrix();  /// R_c_w, camera to world
		cam0_pose.translation() = cam_state_iter->second.position;

		Eigen::Isometry3d cam1_pose;
		if (stereo)
		cam1_pose = cam0_pose * T_cam0_cam1.inverse();

		cam_poses.push_back(cam0_pose);
		if (stereo)
		cam_poses.push_back(cam1_pose);
	}


	Eigen::Isometry3d T_c0_w = cam_poses[0];
	for (auto& pose : cam_poses)
		pose = pose.inverse() * T_c0_w;
	// c0->ci


	Eigen::Vector3d initial_position(0.0, 0.0, 0.0);
	
	generateInitialGuess(cam_poses[cam_poses.size() - 1], measurements[0],
		measurements[measurements.size() - 1], initial_position);

	//if()
	//Eigen::Vector2d z0 = measurements[0];
	//initial_position = Eigen::Vector3d(z0(0)*initial_depth, z0(1)*initial_depth, initial_depth);


	Eigen::Vector3d solution(
		initial_position(0) / initial_position(2),
		initial_position(1) / initial_position(2),
		1.0 / initial_position(2));



	double lambda = initial_damping;
	int inner_loop_cntr = 0;
	int outer_loop_cntr = 0;
	bool is_cost_reduced = false;
	double delta_norm = 0;


	double total_cost = 0.0;
	for (int i = 0; i < cam_poses.size(); ++i)
	{
		double this_cost = 0.0;
		cost(cam_poses[i], solution, measurements[i], this_cost);
		total_cost += this_cost;
	}


	do
	{
		Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
		Eigen::Vector3d b = Eigen::Vector3d::Zero();

		for (int i = 0; i < cam_poses.size(); ++i)
		{
			Eigen::Matrix<double, 2, 3> J;
			Eigen::Vector2d r;
			double w;


			jacobian(cam_poses[i], solution, measurements[i], J, r, w);


			if (w == 1)
			{
				A += J.transpose() * J;
				b += J.transpose() * r;
			}
			else
			{
				double w_square = w * w;
				A += w_square * J.transpose() * J;
				b += w_square * J.transpose() * r;
			}
		}


		do
		{
			Eigen::Matrix3d damper = lambda * Eigen::Matrix3d::Identity();
			Eigen::Vector3d delta = (A + damper).ldlt().solve(b);
			Eigen::Vector3d new_solution = solution - delta;
			delta_norm = delta.norm();

			double new_cost = 0.0;
			for (int i = 0; i < cam_poses.size(); ++i)
			{
				double this_cost = 0.0;
				cost(cam_poses[i], new_solution, measurements[i], this_cost);
				new_cost += this_cost;
			}

			if (new_cost < total_cost)
			{
				is_cost_reduced = true;
				solution = new_solution;
				total_cost = new_cost;
				lambda = lambda / 10 > 1e-10 ? lambda / 10 : 1e-10;
			}
			else
			{
				is_cost_reduced = false;
				lambda = lambda * 10 < 1e12 ? lambda * 10 : 1e12;
			}

		} while (inner_loop_cntr++ <
			inner_loop_max_iteration && !is_cost_reduced);
		inner_loop_cntr = 0;

	} while (outer_loop_cntr++ <
		outler_loop_max_iteration &&
		delta_norm > estimation_precision);




	Eigen::Vector3d final_position(solution(0) / solution(2),
		solution(1) / solution(2), 1.0 / solution(2));
	//cout << "depth:" << 1.0 / solution(2) << endl;
	//cout << "LevelMar:" << final_position.transpose() << endl;
	bool is_valid_solution = true;
	for (const auto& pose : cam_poses)
	{
		Eigen::Vector3d position =
			pose.linear()*final_position + pose.translation();
		if (position(2) <= 0)
		{
			is_valid_solution = false;
			cout << "depth(2)<0:" << endl;
			break;
		}

	}

	position = T_c0_w.linear()*final_position + T_c0_w.translation();

	
	//if(position(2)<1)
	//{
	//	is_valid_solution = false;
	//}

	if (is_valid_solution)
		is_initialized = true;



	return is_valid_solution;
}

#ifdef USE_CERES
bool gvis::t_gfeature::LMOptimizatePosition(const vector<Eigen::Isometry3d>&cam_poses, const vector<Eigen::Vector2d>&measurements, Eigen::Vector3d& pos)
{
	
	ceres::Problem problem;
	Eigen::Isometry3d cam0_tranformation{ Eigen::Isometry3d::Identity() };
	Eigen::Vector3d initial_position_world(0.0, 0.0, 0.0);
	generateInitialGuess(cam_poses[cam_poses.size() - 1], measurements[0],
		measurements[measurements.size() - 1], initial_position_world);

	double depth_initial = initial_position_world[2];
	double position_world[3] = { initial_position_world[0] / initial_position_world[2],
								initial_position_world[1] / initial_position_world[2],
								1.0 / initial_position_world[2] };
	assert(cam_poses.size() == measurements.size());
	for(int i=0;i<cam_poses.size();i++)
	{

		cv::Point2f observation_point = cv::Point2f(measurements.at(i)(0), measurements.at(i)(1));
		Eigen::Isometry3d cam0_cami = cam_poses.at(i);//cam0->cami
		ceres::CostFunction *cost_function = ReprojectionError::Create(observation_point, cam0_cami);
		problem.AddResidualBlock(cost_function, new ceres::HuberLoss(0.1), position_world);
	}
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	options.max_num_iterations = 20;
	ceres::Solver::Summary summary;
	// options.minimizer_progress_to_stdout = true;
	ceres::Solve(options, &problem, &summary);

	pos << position_world[0], position_world[1], 1.0;
	pos *= (1 / position_world[2]);
	
	double delta_depth = fabs(depth_initial - pos(2)) / pos(2);

	//cout << summary.BriefReport() << endl;

	if(delta_depth > 0.25 || summary.final_cost > 1e-3 || summary.num_successful_steps > 10 || pos(2)<1)
	//if (delta_depth > 0.25 || summary.final_cost > 1e-3 || feature.position_world_(2) > 200 || feature.position_world_(2) < 1 || summary.num_successful_steps > 10)
	{
		
		return false;
	}
	// LOG(INFO) << "summary.final_cost: " << summary.final_cost << std::endl;
	//TODO  // feature.position_world_ = cam0_tranformation * feature.position_world_;

	// auto camera_state_end_iter = feature.raw_uv_.rbegin();
	// const auto &camera_state_tmp = map_state_set_[camera_state_end_iter->first];
	// Eigen::Isometry3d cam_n_tran{Eigen::Isometry3d::Identity()};
	// cam_n_tran.rotate(camera_state_tmp.quat_.toRotationMatrix());
	// cam_n_tran.pretranslate(camera_state_tmp.position_);
	// auto tmp_pos = (cam_n_tran.inverse() * feature.position_world_);
	// LOG(ERROR) <<"depth: "<<tmp_pos(2) << std::endl;
	// if (tmp_pos(2) > 50 || tmp_pos(2) < 0.5)
	//     return false;

	return true;
}

bool gvis::t_gfeature::LMOptimizatePositionAndCheck(const CamStateServer& cam_states)
{
	//先把数据位姿、分好 直接传进去 就okLe 
	// for keyframe ,need to keep
	if (observations.size() < 5) return false;

	std::vector<Eigen::Isometry3d> all_cam_poses,part1_cam_poses,part2_cam_poses;
	std::vector<Eigen::Vector2d> all_measurements,part1_measurements,part2_measurements;
	
	size_t i = 0;
	for (auto& m : observations)
	{
		auto cam_state_iter = cam_states.find(m.first);
		if (cam_state_iter == cam_states.end())
			continue;

		all_measurements.push_back(m.second.head<2>());
		if (stereo)
			all_measurements.push_back(m.second.tail<2>());
		if (i % 3 == 0)
		{
			part1_measurements.push_back(m.second.head<2>());
			if (stereo)
				part1_measurements.push_back(m.second.tail<2>());
			part2_measurements.push_back(m.second.head<2>());
			if (stereo)
				part2_measurements.push_back(m.second.tail<2>());
		}
		else if (i % 3 == 1)
		{
			part1_measurements.push_back(m.second.head<2>());
			if (stereo)
				part1_measurements.push_back(m.second.tail<2>());
		}
		else
		{
			part2_measurements.push_back(m.second.head<2>());
			if (stereo)
				part2_measurements.push_back(m.second.tail<2>());
		}


		Eigen::Isometry3d cam0_pose;
		cam0_pose.linear() = cam_state_iter->second.orientation.toRotationMatrix();  /// R_c_w, camera to world
		cam0_pose.translation() = cam_state_iter->second.position;

		Eigen::Isometry3d cam1_pose;
		if (stereo)
			cam1_pose = cam0_pose * T_cam0_cam1.inverse();

		all_cam_poses.push_back(cam0_pose);
		if (stereo)
			all_cam_poses.push_back(cam1_pose);
		if (i % 3 == 0)
		{
			part1_cam_poses.push_back(cam0_pose);
			if (stereo)
				part1_cam_poses.push_back(cam1_pose);
			part2_cam_poses.push_back(cam0_pose);
			if (stereo)
				part2_cam_poses.push_back(cam1_pose);
		}
		else if (i % 3 == 1)
		{
			part1_cam_poses.push_back(cam0_pose);
			if (stereo)
				part1_cam_poses.push_back(cam1_pose);
		}
		else
		{
			part2_cam_poses.push_back(cam0_pose);
			if (stereo)
				part2_cam_poses.push_back(cam1_pose);
		}

		i++;
	}

	Eigen::Isometry3d T_c0_w = all_cam_poses[0];
	for (auto& pose : all_cam_poses)
		pose = pose.inverse() * T_c0_w;
	for (auto& pose : part1_cam_poses)
		pose = pose.inverse() * T_c0_w;
	for (auto& pose : part2_cam_poses)
		pose = pose.inverse() * T_c0_w;

	Eigen::Vector3d all_final_position, part1_final_position, part2_final_position;
	//if(LMOptimizatePosition(all_cam_poses, all_measurements, all_final_position))
	if (LMOptimizatePosition(all_cam_poses,all_measurements,all_final_position) &&
		LMOptimizatePosition(part1_cam_poses,part1_measurements,part1_final_position) &&
		LMOptimizatePosition(part2_cam_poses,part2_measurements,part2_final_position))
	{
		Eigen::Vector3d diff1 = all_final_position - part1_final_position;
		Eigen::Vector3d diff2 = all_final_position - part2_final_position;
		//cout << "LM:" << endl << std::fixed << std::setw(12) << all_final_position.transpose() << endl;
		//cout<< std::setw(12) << part1_final_position.transpose() << endl;
		//cout  << std::setw(12) << part2_final_position.transpose() <<  endl;
		if (abs(diff1(2)) < 8 && abs(diff2(2)) < 8 &&
			abs(diff1(1)) < 2 && abs(diff2(1)) < 2 &&
			abs(diff1(0)) < 2 && abs(diff2(0)) < 2)
		//if(1)
		{
			//cout << "camera_frame_LM:" << all_final_position.transpose() << endl;
			position = T_c0_w * all_final_position;//待定 可能是错误的
			is_initialized = true;
			return true;
		}
	}

	return false;
	//bool is_valid_solution = true;
	//for (const auto& pose : cam_poses)
	//{
	//	Eigen::Vector3d position =
	//		pose.linear()*final_position + pose.translation();
	//	if (position(2) <= 0)
	//	{
	//		is_valid_solution = false;
	//		cout << "depth(2)<0:" << endl;
	//		break;
	//	}

	//}

	//position = T_c0_w.linear()*final_position + T_c0_w.translation();

	////判断一下深度范围
	////if (position(2) > 200 || position(2) < 1)
	//if (position(2) < 1)
	//{
	//	is_valid_solution = false;
	//}

	//if (is_valid_solution)
	//	is_initialized = true;



	//return is_valid_solution;

	
}

#endif // USE_CERES


