#include "gkalman.h"
#include "gbase.h"
using namespace gsins;
using namespace Eigen;


/***************************  class t_gkalman  *********************************/
t_gkalman::t_gkalman()
{

}

t_gkalman::t_gkalman(gnut::t_gsetbase* gset)
{
	kftk = 0.0;
	nq = 15; nr = 6;
	Ft = Pk = MatrixXd::Zero(nq, nq);
	Hk = MatrixXd::Zero(nr, nq);
	Qt = Pmin = Xk = VectorXd::Zero(nq);  Pmax = VectorXd::Ones(nq)* glv.INF;
	Rt = Zk = VectorXd::Zero(nr);  rts = VectorXd::Ones(nr);
	Rk = MatrixXd::Zero(nr, nr);
	Rmax = VectorXd::Ones(nr)* glv.INF; Rmin = Rb = VectorXd::Zero(nr); Rbeta = VectorXd::Ones(70);
	FBTau = FBMax = VectorXd::Ones(nq)* glv.INF; FBXk = FBTotal = VectorXd::Zero(nq);
	measflag = 0; 
}

gsins::t_gkalman::t_gkalman(int _nq, int _nr)
{
	kftk = 0.0;
	Ft = Pk = MatrixXd::Zero(nq, nq);
	Hk = MatrixXd::Zero(nr, nq);
	Qt = Pmin = Xk = VectorXd::Zero(nq);  Pmax = VectorXd::Ones(nq)* glv.INF;
	Rt = Zk = VectorXd::Zero(nr);  rts = VectorXd::Ones(nr);
	Rk = MatrixXd::Zero(nr, nr);
	Rmax = VectorXd::Ones(nr)* glv.INF; Rmin = Rb = VectorXd::Zero(nr); Rbeta = VectorXd::Ones(70);
	FBTau = FBMax = VectorXd::Ones(nq)* glv.INF; FBXk = FBTotal = VectorXd::Zero(nq);
	measflag = 0;
}

void t_gkalman::time_update(double kfts, int fback)
{
	set_Ft();
	MatrixXd Fk = MatrixXd::Identity(nq, nq) + (Ft*kfts);
	Xk = Fk * Xk;
	MatrixXd Qk = (Qt*kfts).array().matrix().asDiagonal();
	Pk = Fk*Pk*(Fk.transpose());  Pk += Qk;
	if (fback)  feedback(kfts);
}

void t_gkalman::set_meas_flag(int flag)
{
	measflag = (flag == 0) ? 0 : (measflag | flag);
}

/************************************************************
** Function: Kalman Filter Measurement Update
** include: Sequential filtering¡¢Measurement noise adaptation
***********************************************************/
void t_gkalman::meas_update(double fading)
{
	Pk_Sav = Pk;
	VectorXd Pxz(nq), Kk(nq), Hi(nq);
	set_Meas();
	for (int i = 0; i < nr; i++)
	{
		if (measflag&(0x01 << i))
		{
			Hi = Hk.block(i, 0, 1, nq).transpose();
			Pxz = Pk*Hi;
			double Pz0 = (Hi.transpose()*Pxz), r = Zk(i) - (Hi.transpose()*Xk);
			//R_adaptive(i, r, Pz0);
			double Pzz = Pz0 + Rt(i) / rts(i);
			Kk = Pxz*(1.0 / Pzz);
			Xk += Kk*r;
			Pk = Pk - Kk*Pxz.transpose();
		}
	}
	if (fading > 1.0) Pk *= fading;
	//Pk_constrain();
	t_gbase::symmetry(Pk);
	set_meas_flag(0);

}

/************************************************************
** Function: Measurement noise adaptation
** See also P128
***********************************************************/
void t_gkalman::R_adaptive(int i, double r, double Pr)
{
	if (Rb(i) > glv.EPS)
	{
		double rr = r*r - Pr;
		if (rr < Rmin(i))	
			rr = Rmin(i);
		if (rr > Rmax(i))	
			Rt(i) = Rmax(i);
		else
			Rt(i) = (1.0 - Rbeta(i))*Rt(i) + Rbeta(i) * rr;
		Rbeta(i) = Rbeta(i) / (Rbeta(i) + Rb(i));
	}
}

/************************************************************
** Function: Variance constraint
***********************************************************/
void t_gkalman::Pk_constrain()
{
	for (int i = 0; i < nq; i++)
	{
		if (Pk(i, i) < Pmin(i) /*&& Pk(i, i) > glv.EPS*/)
		{
			Pk(i, i) = Pmin(i);
		}
		else if (Pk(i, i) > Pmax(i))
		{
			double sqf = sqrt(Pmax(i) / Pk(i, i));//* 0.5;
			for (int j = 0; j < nq; j++)
			{
				Pk(j, i) *= sqf;
			}
			for (int j = 0; j < nq; j++)
			{
				Pk(i, j) *= sqf;
			}
		}
	}
}

/************************************************************
** Function: Feedback correction
***********************************************************/
void t_gkalman::feedback(double fbts)
{
	for (int i = 0; i < nq; i++)
	{
		if (FBTau(i) < glv.INF / 2.0)
		{
			double afa = fbts<FBTau(i) ? fbts / FBTau(i) : 1.0;
			FBXk(i) = Xk(i)*afa;
			if (FBTotal(i) + FBXk(i)>FBMax(i))
				FBXk(i) = FBMax(i) - FBTotal(i);
			else if (FBTotal(i) + FBXk(i) < -FBMax(i))
				FBXk(i) = -FBMax(i) - FBTotal(i);
			Xk(i) -= FBXk(i);
			FBTotal(i) += FBXk(i);
		}
		else
		{
			FBXk(i) = 0.0;
		}
	}
}

int t_gkalman::outlier_detect(double threshold)
{
	double Vstd = _v_post(3)*_v_post(3) + _v_post(4)*_v_post(4) + _v_post(5)*_v_post(5);
	if (Vstd > threshold)return 0;
	else return 1;
}

int gsins::t_gkalman::_IGG_detect()
{
	double k0 = 3, k1 = 8;
	_v_post = Zk - Hk*Xk;
	Eigen::MatrixXd Qz = Hk*Pk*Hk.transpose() + Rk;
	Eigen::VectorXd v_norm = Eigen::VectorXd::Zero(nr);
	for (int i = 0; i < nr; i++)
	{
		v_norm(i) = sqrt(1.0 / Qz(i, i))*abs(_v_post(i));
	}
	for (int i = 0; i < nr; i++)
	{
		if (v_norm(i) <= k0)Rk(i, i) = Rk(i, i);
		else if (v_norm(i) > k0&&v_norm(i) <= k1)
		{
			double gamma = v_norm(i) / k0 * (k1 - k0) / (k1 - v_norm(i));
			Rk(i, i) = gamma*Rk(i, i);
		}
		else
		{
			Rk(i, i) = glv.INF*Rk(i, i);
		}
	}
	return 0;
}

