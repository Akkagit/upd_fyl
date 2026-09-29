#include "gutility.h"
#include "gset/gsetins.h"
#include "gset/gsetign.h"
using namespace gsins;
using namespace gnut;
const t_gglv gsins::glv = t_gglv();

t_gglv::t_gglv(double Re, double f, double wie0, double g0)
{
	this->Re = Re;
	this->f = f;
	this->wie = wie0;
	this->g0 = g0;
	mg = 1.0e-3*g0;
	ug = 1.0e-6*glv.g0;
	mgpsHz = mg / sqrt(1.0);
	ugpsHz = ug / sqrt(1.0);
	mgpsh = mg / sqrt(hur);
	ugpsh = ug / sqrt(hur);

}


double gsins::Lagrange(double t[], double x[], int n, double t1)
{
	double res = 0.0;
	double *L = new double[n];

	for (int k = 0; k<n; k++)
	{
		double up = 1.0, down = 1.0;
		for (int m = 0; m < n; m++)
		{
			if (m == k)continue;
			up *= (t1 - t[m]);
			down *= (t[k] - t[m]);
		}
		L[k] = up / down;
	}

	for (int i = 0; i<n; i++)
	{
		res += x[i] * L[i];
	}
	delete[] L;
	return res;
}

int gsins::sign(double d)
{
	if (d > glv.EPS)
		return 1;
	else if (d < -glv.EPS)
		return -1;
	else
		return 0;
}

UNIT_TYPE gsins::str2Unit(const string & s)
{
	string tmp = s;
	transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
	UNIT_TYPE res;
	if (tmp == "RAD")
		res = RAD;
	else if (tmp == "DEG")
		res = DEG;
	else if (tmp == "RPS")
		res = RPS;
	else if (tmp == "DPS")
		res = DPS;
	else if (tmp == "RPH")
		res = RPH;
	else if (tmp == "DPH")
		res = DPH;
	else if (tmp == "MPS")
		res = MPS;
	else if (tmp == "MPS2")
		res = MPS2;	 
	return res;
}

FLT_TYPE gsins::str2Flt(const string & s)
{
	string tmp = s;
	transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
	if (tmp == "FORWARD")return FORWARD;
	if (tmp == "BACKWARD")return BACKWARD;
	if (tmp == "FBS")return FBS;
	if (tmp == "RTS")return RTS;
	return FORWARD;
}

IGN_TYPE gsins::str2ign(const string & s)
{
	string tmp = s;
	transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
	if (tmp == "PURE_INS")return PURE_INS;
	if (tmp == "LCI")return LCI;
	if (tmp == "TCI")return TCI;
	if (tmp == "STCI")return STCI;
	if (tmp == "VIO") return VIO;
	if (tmp == "MULTIGN_LC")return MULTIGN_LC;
	if (tmp == "MULTIGN_TC")return MULTIGN_TC;
	if (tmp == "MULTIGN_STC")return MULTIGN_STC;
	if (tmp == "MULTIGN_LIDAR_TC")return MULTIGN_LIDAR_TC;
	if (tmp == "MULTIGN_LIDAR_LC")return MULTIGN_LIDAR_LC;
	if (tmp == "MULTIGN_LIDAR_STC")return MULTIGN_LIDAR_STC;
	if (tmp == "LIDAR_INS")return LIDAR_INS;
	return PURE_INS;
}

ODR_TYPE gsins::str2odr(const string & s)
{
	string tmp = s;
	transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
	if (tmp == "XYZ_XYZ_PRY_NW")return XYZ_XYZ_PRY_NW;
	if (tmp == "BLH_ENU_PRY_NW")return BLH_ENU_PRY_NW;
	if (tmp == "HOLO_ODO")return HOLO_ODO;
	return XYZ_XYZ_PRY_NW;
}

LCI_TYPE gsins::str2lci(const string & s)
{
	string tmp = s;
	transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
	if (tmp == "SPP")return LCI_TYPE::SPP;
	if (tmp == "PPP")return LCI_TYPE::PPP;
	if (tmp == "RTK")return LCI_TYPE::RTK;
	if (tmp == "PPPRTK")return LCI_TYPE::PPPRTK;
	return LCI_TYPE::SPP;
}

ALIGN_TYPE gsins::str2align(const string & s)
{
	string tmp = s;
	transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
	if (tmp == "AUTO")return ALIGN_TYPE::AUTO;
	if (tmp == "SA")return ALIGN_TYPE::SA;
	if (tmp == "MA")return ALIGN_TYPE::MA;
	if (tmp == "PVA")return ALIGN_TYPE::PVA;
	if (tmp == "VVA")return ALIGN_TYPE::VVA;
	return ALIGN_DEFAULT;
}

MEAS_TYPE gsins::ins2meas(const int & i)
{
	return MEAS_TYPE();
}

CPS_TYPE gsins::str2cps(const string & s)
{
	string tmp = s;
	transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
	if (tmp == "POLY")return POLY;
	if (tmp == "CONE")return CONE;
	if (tmp == "ONE_PLUS_PRE")return ONE_PLUS_PRE;
	return POLY;
}

Eigen::Vector3d gsins::Cart2Geod(const Eigen::Vector3d& CartPos, bool b)
{
	double XYZ[3] = { CartPos(0),CartPos(1) ,CartPos(2) }, ell[3];
	gnut::xyz2ell(XYZ, ell, b);
	return Eigen::Vector3d(ell);
}

Eigen::Vector3d gsins::Geod2Cart(const Eigen::Vector3d& GeodPos, bool b)
{
	double ell[3] = { GeodPos(0),GeodPos(1) ,GeodPos(2) }, XYZ[3];
	gnut::ell2xyz(ell, XYZ, b);
	return Eigen::Vector3d(XYZ);
}

INS_LIBRARY_EXPORT Eigen::Vector3d gsins::XYZ2ENU(const Eigen::Vector3d & xyz, const Eigen::Vector3d & xyz_ref)
{
	double XYZ[3] = { xyz(0),xyz(1) ,xyz(2) }, XYZ_REF[3] = { xyz_ref(0),xyz_ref(1) ,xyz_ref(2) }, neu[3];
	gnut::xyz2neu(XYZ, XYZ_REF, neu);
	// east north up
	return Eigen::Vector3d(neu[1], neu[0], neu[2]);
}

Eigen::MatrixXd gsins::NewMat2Eigen(const Matrix & newmat)
{
	int m = newmat.Nrows(), n = newmat.Ncols();
	Eigen::MatrixXd res(m, n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			res(i, j) = newmat(i + 1, j + 1);
		}

	}
	return  res;
}

Matrix gsins::Eigen2newMat(const Eigen::MatrixXd & eigen)
{
	int m = eigen.rows(), n = eigen.cols();
	Matrix res; res.ReSize(m, n); res = 0.0;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			res(i + 1, j + 1) = eigen(i, j);
		}
	}
	return res;
}

INS_LIBRARY_EXPORT SymmetricMatrix gsins::Eigen2BaseMatrix(const Eigen::MatrixXd & eigen)
{
	int m = eigen.rows(), n = eigen.cols();
	SymmetricMatrix res(m);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			res(i + 1, j + 1) = eigen(i, j);
		}
	}
	return  res;
}

INS_LIBRARY_EXPORT Eigen::MatrixXd gsins::BaseMatrix2Eigen(const SymmetricMatrix & newmat)
{
	int m = newmat.Nrows(), n = newmat.Ncols();
	Eigen::MatrixXd res(m, n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			res(i, j) = newmat(i + 1, j + 1);
		}
	}
	return  res;
}

INS_LIBRARY_EXPORT Eigen::VectorXd gsins::Columns2VectorXd(const ColumnVector & newmat)
{
	int m = newmat.Nrows();
	Eigen::VectorXd res(m);
	for (int i = 0; i < m; i++)
	{
		res(i) = newmat(i + 1);
	}
	return  res;
}

INS_LIBRARY_EXPORT ColumnVector gsins::VectorXd2Columns(const Eigen::VectorXd & eigen)
{
	int m = eigen.size();
	ColumnVector res; res.ReSize(m); res = 0.0;
	for (int i = 0; i < m; i++)
	{
		res(i + 1) = eigen(i);
	}
	return  res;
}

gsins::t_scheme::t_scheme(t_gsetbase * set)
{
	t = 0.0;
	ts = dynamic_cast<t_gsetins*>(set)->ts();
	freq = dynamic_cast<t_gsetins*>(set)->freq();
	nSamples = dynamic_cast<t_gsetins*>(set)->subsample();
	Cps = dynamic_cast<t_gsetins*>(set)->cps();
	align = dynamic_cast<t_gsetins*>(set)->align();
	start = dynamic_cast<t_gsetins*>(set)->start();
	end = dynamic_cast<t_gsetins*>(set)->end();
	align_time = dynamic_cast<t_gsetins*>(set)->align_time();

	delay = dynamic_cast<t_gsetign*>(set)->delay_t();
	max_pdop = dynamic_cast<t_gsetign*>(set)->max_pdop();
	min_sat = dynamic_cast<t_gsetign*>(set)->min_sat();
}
