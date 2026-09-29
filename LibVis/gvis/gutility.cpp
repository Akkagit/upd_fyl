
#include "gutility.h"
#include "gset/gsetvis.h"

using namespace gnut;
using namespace gvis;



Distortion_Type gvis::str2distortion(const string & s)
{
	string tmp = s;
	transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
	if (tmp == "radtan")return radtan;
	if (tmp == "equidistant")return equidistant;

	return radtan;

}

string gvis::distortion2str(const Distortion_Type & s)
{
	//string res;
	if (s == 0)
		return string("radtan");
	else if (s == 1)
		return string("equidistant");
	return string("radtan");
}


Eigen::Matrix3d gvis::skew(const Eigen::Vector3d& v)
{
	Eigen::Matrix3d vnx;
	vnx << 0, -v(2), v(1),
		v(2), 0, -v(0),
		-v(1), v(0), 0;
	return vnx;
}

//e系到n系
Eigen::Matrix3d gvis::R_ENU_ECEF(const Eigen::Vector3d &BLH)
{
	Eigen::Matrix3d ans;
	double sinB = sin(BLH(0));
	double cosB = cos(BLH(0));
	double sinL = sin(BLH(1));
	double cosL = cos(BLH(1));
	ans(0, 0) = -sinL;			ans(0, 1) = cosL;				ans(0, 2) = 0;//E
	ans(1, 0) = -cosL * sinB;		ans(1, 1) = -sinL * sinB;			ans(1, 2) = cosB;//N
	ans(2, 0) = cosL * cosB;		ans(2, 1) = sinL * cosB;			ans(2, 2) = sinB;//U
	return ans;
}

Eigen::Vector3d gvis::XYZ2BLH(const Eigen::Vector3d &X)
{
	double a1 = 6378137;
	double e2 = 0.00669437999013;
	double pi = 4 * atan(1);
	Eigen::Vector3d BLH;
	double B1, B2, N;
	BLH(1) = atan(X(1) / X(0));
	BLH(1) = BLH(1);
	if (X(0) < 0 && X(1) > 0) BLH(1) = BLH(1) + pi;
	if (X(0) < 0 && X(1) < 0) BLH(1) = BLH(1) - pi;
	B1 = atan(X(2) / sqrt(X(0)*X(0) + X(1)*X(1)));
	N = a1 / sqrt(1 - e2 * sin(B1)*sin(B1));
	B2 = atan((X(2) + N * e2*sin(B1)) / sqrt(X(0)*X(0) + X(1)*X(1)));
	//迭代法求大地纬度;
	while (abs(B1 - B2) > (0.001 / 3600 * pi / 180))
	{
		B1 = B2;
		N = a1 / sqrt(1 - e2 * sin(B1)*sin(B1));
		B2 = atan((X(2) + N * e2*sin(B1)) / sqrt(X(0)*X(0) + X(1)*X(1)));
	}
	BLH(0) = B2;
	BLH(2) = X(2) / sin(B2) - N * (1 - e2);
	return (BLH);
}
//vector<string> gvis::split(const string &s, const string &seperator) 
//{
//	vector<string> result;
//	typedef string::size_type string_size;
//	string_size i = 0;
//
//	while (i != s.size()) {
//		int flag = 0;
//		while (i != s.size() && flag == 0) {
//			flag = 1;
//			for (string_size x = 0; x < seperator.size(); ++x)
//				if (s[i] == seperator[x]) {
//					++i;
//					flag = 0;
//					break;
//				}
//		}
//
//		flag = 0;
//		string_size j = i;
//		while (j != s.size() && flag == 0) {
//			for (string_size x = 0; x < seperator.size(); ++x)
//				if (s[j] == seperator[x]) {
//					flag = 1;
//					break;
//				}
//			if (flag == 0)
//				++j;
//		}
//		if (i != j) {
//			result.push_back(s.substr(i, j - i));
//			i = j;
//		}
//	}
//	return result;
//}

//Eigen::Quaterniond gvis::RotationVector2Quaternion(const Eigen::Vector3d &R)
//{
//	Eigen::Vector4d q = Eigen::Vector4d::Zero();
//	double n2 = R(0)*R(0) + R(1)*R(1) + R(2)*R(2);
//	double n = 0;
//	double n_2 = 0;
//	double s = 0;
//	if (n2 < 1.0e-8)  //cos(n / 2) = 1 - n2 / 8 + n4 / 384; sin(n / 2) / n = 1 / 2 - n2 / 48 + n4 / 3840
//	{
//		q(0) = 1 - n2 * (1 / 8.0 - n2 / 384.0); s = 1 / 2.0 - n2 * (1 / 48.0 - n2 / 3840.0);
//	}
//	else
//	{
//		n = sqrt(n2); n_2 = n / 2;
//		q(0) = cos(n_2); s = sin(n_2) / n;
//	}
//
//	q(1) = s * R(0); q(2) = s * R(1); q(3) = s * R(2);
//
//	Eigen::Quaterniond out;
//	out.w() = q(0);
//	out.x() = q(1);
//	out.y() = q(2);
//	out.z() = q(3);
//	return out;
//}

