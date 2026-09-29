#include "gavar.h"
using namespace gsins;
using namespace Eigen;
using namespace std;

gsins::t_gavar::t_gavar()
{
}

/******************************************* t_gavar ***********************************************/
t_gavar::t_gavar(int n, double ts, const VectorXd& r0, const VectorXd& tau)
{
	init(n, ts, r0, tau);
}

void gsins::t_gavar::init(int n, double ts, const Eigen::VectorXd & r0, const Eigen::VectorXd & tau)
{
	this->n = n; this->ts = ts;
	flag = R = Rmax = Rmin = tstau = this->r0 = VectorXd::Zero(n);
	R = r0.array().square();
	tstau = tau.cwiseInverse()*ts;
	Rmax = 100 * R; Rmin = 0.01*R;
}

/************************************************************
** Function: update Variance
** Input: Sampling data sequence
** See also 《惯性仪器测试和数据分析》P199 - P200
***********************************************************/
void t_gavar::update(const VectorXd& r)
{
	VectorXd dr = r - r0;
	VectorXd dr2 = dr.array().square();
	for (int i = 0; i < n; i++)
	{
		if (dr2[i] > R[i]) R[i] = dr2[i];
		else R[i] = (1 - tstau[i])*R[i] + tstau[i] * dr2[i];
		if (R[i] < Rmin[i]) R[i] = Rmin[i];
		if (R[i] > Rmax[i])
		{
			R[i] = Rmax[i];
			flag[i] = 1;
		}
		else
			flag[i] = 0;
	}
	r0 = r;
}

double t_gavar::operator()(int i)const
{
	return flag[i] ? glv.INF : sqrt(R[i]);
}

