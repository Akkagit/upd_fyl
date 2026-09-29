#include "glci.h"
#include "gproc/gpreproc.h"
#include "gutils/gtimesync.h"
#include <iostream>
using namespace gsins;
using namespace std;


/******************************************* t_glci ***********************************************/
t_glci::t_glci():t_gsinskf()
{
}

t_glci::t_glci(t_gsetbase * gset) : t_gsinskf(gset)
{
	kftk = 0.0;
	nq = dynamic_cast<t_gsetign*>(gset)->nq();
	nr = dynamic_cast<t_gsetign*>(gset)->nr();
	Ft = Pk = Eigen::MatrixXd::Zero(nq, nq);
	Hk = Eigen::MatrixXd::Zero(nr, nq);
	Qt = Pmin = Xk = Eigen::VectorXd::Zero(nq);  Pmax = Eigen::VectorXd::Ones(nq)* glv.INF;
	Rt = Zk = Eigen::VectorXd::Zero(nr);  rts = Eigen::VectorXd::Ones(nr);
	Rk = Eigen::MatrixXd::Zero(nr, nr);
	Rmax = Eigen::VectorXd::Ones(nr)* glv.INF; Rmin = Rb = Eigen::VectorXd::Zero(nr); Rbeta = Eigen::VectorXd::Ones(nr);
	FBTau = FBMax = Eigen::VectorXd::Ones(nq)* glv.INF; FBXk = FBTotal = Eigen::VectorXd::Zero(nq);
	measflag = 0;
	FilterFlag = dynamic_cast<t_gsetign*>(gset)->fltmode();
}

void t_glci::init15()
{
	//sins = sins0;
	kftk = sins.t; MeasVel = MeasPos = Eigen::Vector3d::Zero(), Flag = NO_MEAS, fltflag = 0;

	this->lever = dynamic_cast<t_gsetign*>(_setkf)->lever();

	Eigen::VectorXd v(nq), vmax(nq), vmin(nq);

	Eigen::Vector3d angle_tmp, v_tmp, pos_tmp, eb_tmp, db_tmp;
	angle_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_misalignment_std();
	v_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_vel_std();
	pos_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_pos_std();
	eb_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_gyro_std();
	db_tmp = dynamic_cast<t_gsetign*>(_setkf)->initial_acce_std();
	pos_tmp(0) = pos_tmp(0) / glv.Re; pos_tmp(1) = pos_tmp(1) / glv.Re;
	v << angle_tmp*glv.deg, v_tmp, pos_tmp, eb_tmp*glv.dph, db_tmp*glv.mg;
	v = v.array().abs2();
	Pk = v.array().matrix().asDiagonal();

	angle_tmp = dynamic_cast<t_gsetign*>(_setkf)->min_misalignment_std();
	v_tmp = dynamic_cast<t_gsetign*>(_setkf)->min_vel_std();
	pos_tmp = dynamic_cast<t_gsetign*>(_setkf)->min_pos_std();
	eb_tmp = dynamic_cast<t_gsetign*>(_setkf)->min_gyro_std();
	db_tmp = dynamic_cast<t_gsetign*>(_setkf)->min_acce_std();
	pos_tmp(0) = pos_tmp(0) / glv.Re; pos_tmp(1) = pos_tmp(1) / glv.Re;
	vmin << angle_tmp*glv.min, v_tmp, pos_tmp, eb_tmp*glv.dph, db_tmp*glv.mg;
	Pmin = vmin.array().abs2();

	vmax = vmin * 1e6;
	Pmax = vmax.array().abs2();

	angle_tmp = dynamic_cast<t_gsetign*>(_setkf)->misalignment_psd();
	v_tmp = dynamic_cast<t_gsetign*>(_setkf)->vel_psd();
	pos_tmp = dynamic_cast<t_gsetign*>(_setkf)->pos_psd();
	eb_tmp = dynamic_cast<t_gsetign*>(_setkf)->gyro_psd();
	db_tmp = dynamic_cast<t_gsetign*>(_setkf)->acce_psd();
	Qt << angle_tmp*glv.dpsh, v_tmp*glv.mgpsHz, pos_tmp, eb_tmp*glv.dphpsh, db_tmp*glv.mgpsh;
	Qt = Qt.array().abs2();

	v_tmp = dynamic_cast<t_gsetign*>(_setkf)->meas_vel_noise();
	pos_tmp = dynamic_cast<t_gsetign*>(_setkf)->meas_pos_noise();
	pos_tmp(0) = pos_tmp(0) / glv.Re; pos_tmp(1) = pos_tmp(1) / glv.Re;
	Rt << v_tmp, pos_tmp;
	Rt = Rt.array().abs2();

	Rmax = Rt * 1000; Rmin = Rt*0.01;
	Rb << 0.9, 0.9, 0.9, 0.9, 0.9, 0.9;
	FBTau << 1.0, 1.0, 10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0;
	FBMax << glv.INF, glv.INF, glv.INF, glv.INF, glv.INF, glv.INF, glv.INF, glv.INF, glv.INF,
		3000.0*glv.dph, 3000.0*glv.dph, 3000.0*glv.dph, 50.0*glv.mg, 50.0*glv.mg, 50.0*glv.mg;
}


/************************************************************
** Function: Loose Integration Continuous time measurement vector coefficient matrix
** Ft = [ O33, I33, O33, O33, O33
**		  O33, O33, I33, O33, O33 ];
** See also P159 - P161
***********************************************************/
void t_glci::set_Hk()
{
	switch (Flag)
	{
	case POS_VEL_MEAS:
		Hk.block(0, 3, 3, 3) = Eigen::Matrix3d::Identity();
		Hk.block(3, 6, 3, 3) = Eigen::Matrix3d::Identity();
		break;
	case POS_MEAS:
		Hk.block(3, 6, 3, 3) = Eigen::Matrix3d::Identity();
		break;
	case ZUPT_MEAS:
		Hk.block(0, 3, 3, 3) = Eigen::Matrix3d::Identity();
		break;
	case NHC_MEAS:
		Hk.block(0, 3, 3, 3) = sins.Cbn;
		Hk.block(1, 0, 1, nq) = Eigen::VectorXd::Zero(nq);
		break;
	case NO_MEAS:
		Hk = Eigen::MatrixXd::Zero(nr, nq);
		break;
	default:
		break;
	}
}

/************************************************************
** Function: Loose Integration Continuous time measurement vector
** Z = [ sins.vn - gnss.vn
**		sins,pos - gnss.pos ]
** Taking into account the synchronization of time and space
***********************************************************/
void t_glci::set_Meas()
{
	double dt = sins.t - tmeas;
	if (fabs(dt) > 0.5) return;
	switch (Flag)
	{
	case POS_VEL_MEAS:
		Zk.block(0, 0, 3, 1) = sins.vn - MeasVel - sins.an*dt + sins.Cnb*t_gbase::askew(sins.web)*lever;
		set_meas_flag(007);
		Zk.block(3, 0, 3, 1) = sins.pos - MeasPos - sins.eth.v2dp(sins.vn, dt) + sins.Fpv*sins.Cnb*lever;
		set_meas_flag(070);
		MeasVel = MeasPos = Eigen::Vector3d::Zero();
		break;
	case POS_MEAS:
		Zk.block(3, 0, 3, 1) = sins.pos - MeasPos - sins.eth.v2dp(sins.vn, dt) + sins.Fpv*sins.Cnb*lever;
		set_meas_flag(070);
		MeasVel = MeasPos = Eigen::Vector3d::Zero();
		break;
	case ZUPT_MEAS:
		Zk.block(0, 0, 3, 1) = MeasVel - sins.an*dt;
		set_meas_flag(007);
		MeasVel = Eigen::Vector3d::Zero();
		break;
	case NHC_MEAS:
		//Zk.block(0, 0, 3, 1) = MeasVel - sins.Cbn*t_gbase::askew(sins.wnb)*sins.vn;
		Zk.block(0, 0, 3, 1) = MeasVel;
		MeasVel = Eigen::Vector3d::Zero();
		set_meas_flag(005);
		break;
	case YAW_MEAS:

		MeasYaw = 0;
		break;
	case NO_MEAS:
		set_meas_flag(0);
		break;
	default:
		break;
	}
	if (this->measflag)set_Hk();
}

void gsins::t_glci::set_Rk(const Eigen::Vector3d & Qvn, const Eigen::Vector3d & Qpos)
{
	Rt.block(0, 0, 3, 1) = Qvn;
	Rt.block(3, 0, 3, 1) = Qpos;
}

void t_glci::set_Meas(const Eigen::Vector3d& vn, const Eigen::Vector3d& pos, double t)
{
	if (vn != Eigen::Vector3d::Zero()) MeasVel = vn;
	if (pos != Eigen::Vector3d::Zero()) MeasPos = pos;
	this->tmeas = t;
	Flag = POS_VEL_MEAS;
}

int gsins::t_glci::outlier_detect(double threshold)
{
	_v_post = Zk - Hk*Xk;
	Eigen::MatrixXd Rz = Rt.array().matrix().asDiagonal();
	Eigen::MatrixXd Qz = Hk*Pk*Hk.transpose() + Rz;
	Eigen::VectorXd v_norm = Eigen::VectorXd::Zero(nr);
	for (int i = 0; i < nr; i++)
	{
		v_norm(i) = sqrt(1 / Qz(i, i))*abs(_v_post(i));
	}

	for (int i = 0; i < nr; i++)
	{
		if (v_norm(i) > threshold)
		{
			Pk = Pk_Sav;
			Xk = Eigen::VectorXd::Zero(nq);
			return 1;
		}
	}
	return 0;
}

void t_glci::Update(const vector<Eigen::Vector3d>& wm, const vector<Eigen::Vector3d>& vm, const t_scheme& scm)
{
	sins.Update(wm, vm, scm);
	time_update(sins.nts, 0.0);  kftk = sins.t;
	meas_update();
	if (Flag != NO_MEAS && outlier_detect(3))
	{
		Flag = NO_MEAS;
	}
	feedback(glv.INF);
	//if (fmod(kftk, 1)<0.0076)
	//	t_gsinskf::write();
	if (Flag != NO_MEAS)
	{
		//mc_update();
		t_gsinskf::write();
		if (FilterFlag == FBS)
		{
			if (fltflag == 0)
			{
				Xk_forward.insert(make_pair(kftk, Xk));
				Pk_forward.insert(make_pair(kftk, Pk));
				sins_forward.insert(make_pair(kftk, sins));
			}
			else
			{
				Xk_backward.insert(make_pair(kftk, Xk));
				Pk_backward.insert(make_pair(kftk, Pk));
				sins_backward.insert(make_pair(kftk, sins));
			}
		}
	}
}

void t_glci::feedback(double kfts)
{
	try
	{
		t_gkalman::feedback(kfts);
	}
	catch (out_of_range)
	{
		cout << "t_glci::feedback" + to_string(kftk) + ": " + "FeedBack Failed!" << endl;
	}
	sins.qnb = sins.qnb - FBXk.block(0, 0, 3, 1);
	sins.vn = sins.vn - FBXk.block(3, 0, 3, 1);
	sins.pos = sins.pos - FBXk.block(6, 0, 3, 1);
	sins.eb = sins.eb + FBXk.block(9, 0, 3, 1);
	sins.db = sins.db + FBXk.block(12, 0, 3, 1);
}
