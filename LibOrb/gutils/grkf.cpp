/**
*
* @verbatim
	History
	 -1.0 Jiande Huang	2019-01-04 creat the file.
	 -1.1 Jie Li		2019-04-02 add doxygen style code
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			grkf.CPP
* @brief		CRS coordinate transform to ACR coordinate
*
* @author       Jiande Huang, Wuhan University
* @version		1.0.0
* @date			2019-01-04
*
*/

#include "gutils/grkf.h"
#include <iostream>
#include <math.h>

namespace great
{
	/** @brief show degree and step of runge kutta. */
	void t_grkf::show()
	{
		std::cout << "degree : " << _degree << " h :" << _h << std::endl;
	}

	/** @brief constructor. */
	t_grkf::t_grkf(double step_h, t_gtime tn, double Xn)
	{
		_h = step_h;
		_Xn = Xn;
		_t = tn;
		_degree = 9;
	}

	/** @brief constructor. */
	t_grkf6::t_grkf6()
	{
		_degree = 9;
	}

	/** @brief constructor. */
	t_grkf6::t_grkf6(double step_h, t_gtime tn, double Xn)
		:t_grkf(step_h, tn, Xn)
	{
		_h = step_h;
		_Xn = Xn;
		_t = tn;
		_degree = 9;
		_fi[0] = _fi[1] = _fi[2] = _fi[3] = _fi[4] = _fi[5] = _fi[6] = _fi[7] = _fi[8] = 0.0;
	}

	/** @brief calculate ti and Xi and storage fi_1.
	* @param[in]   i
	* @param[in]   fi_1		partial of X(i-1)
	* @param[out]  ti		ti
	* @param[out]  Xi		Xi
	*/
	bool t_grkf6::updateRkfPoint(const int& i, const double& fi_1, t_gtime& ti, double & Xi)
	{
		double sod = 0.0;
		int isod = 0;
		switch (i)
		{
		case 0:
		{
			_f0 = fi_1;
			sod = _t.sod() + _t.dsec() + _a1 * _h;
			isod = floor(sod);
			ti.from_mjd(_t.mjd(), isod, sod - isod);
			Xi = _Xn + _b10 * _f0 * _h;		//x1
			return true;
		}
		case 1:
		{
			_fi[0] = fi_1;
			sod = _t.sod() + _t.dsec() + _a2 * _h;
			isod = floor(sod);
			ti.from_mjd(_t.mjd(), isod, sod - isod);
			Xi = _Xn + _b21 * _fi[0] * _h;	//x2
			return true;
		}
		case 2:
		{
			_fi[1] = fi_1;
			sod = _t.sod() + _t.dsec() + _a3 * _h;
			isod = floor(sod);
			ti.from_mjd(_t.mjd(), isod, sod - isod);
			Xi = _Xn + _b30 * _f0*_h + _b32 * _fi[1] * _h;
			return true;
		}
		case 3:
		{
			_fi[2] = fi_1;
			sod = _t.sod() + _t.dsec() + _a4 * _h;
			isod = floor(sod);
			ti.from_mjd(_t.mjd(), isod, sod - isod);
			Xi = _Xn + _b40 * _f0* _h +
				_b42 * _fi[1] * _h +
				_b43 * _fi[2] * _h;
			return true;
		}
		case 4:
		{
			_fi[3] = fi_1;
			sod = _t.sod() + _t.dsec() + _a5 * _h;
			isod = floor(sod);
			ti.from_mjd(_t.mjd(), isod, sod - isod);
			Xi = _Xn + _b50 * _f0* _h +
				_b52 * _fi[1] * _h +
				_b53 * _fi[2] * _h +
				_b54 * _fi[3] * _h;
			return true;
		}
		case 5:
		{
			_fi[4] = fi_1;
			sod = _t.sod() + _t.dsec() + _a6 * _h;
			isod = floor(sod);
			ti.from_mjd(_t.mjd(), isod, sod - isod);
			Xi = _Xn + _b60 * _f0 * _h +
				_b62 * _fi[1] * _h +
				_b63 * _fi[2] * _h +
				_b64 * _fi[3] * _h +
				_b65 * _fi[4] * _h;
			return true;
		}
		case 6:
		{
			_fi[5] = fi_1;
			sod = _t.sod() + _t.dsec() + _a7 * _h;
			isod = floor(sod);
			ti.from_mjd(_t.mjd(), isod, sod - isod);
			Xi = _Xn + _b70 * _f0* _h +
				_b72 * _fi[1] * _h +
				_b73 * _fi[2] * _h +
				_b74 * _fi[3] * _h +
				_b75 * _fi[4] * _h +
				_b76 * _fi[5] * _h;
			return true;
		}
		case 7:
		{
			_fi[6] = fi_1;
			ti = _t;
			Xi = _Xn + _b80 * _f0* _h +
				_b83 * _fi[2] * _h +
				_b84 * _fi[3] * _h +
				_b85 * _fi[4] * _h +
				_b86 * _fi[5] * _h;
			return true;
		}
		case 8:
		{
			_fi[7] = fi_1;
			sod = _t.sod() + _t.dsec() + _a9 * _h;
			isod = floor(sod);
			ti.from_mjd(_t.mjd(), isod, sod - isod);
			Xi = _Xn + _b90 * _f0* _h +
				_b92 * _fi[1] * _h +
				_b93 * _fi[2] * _h +
				_b94 * _fi[3] * _h +
				_b95 * _fi[4] * _h +
				_b96 * _fi[5] * _h +
				_b98 * _fi[7] * _h;
			return true;
		}
		case 9:
		{
			_fi[8] = fi_1;
			sod = _t.sod() + _t.dsec() + _h;
			isod = floor(sod);
			ti.from_mjd(_t.mjd(), isod, sod - isod);
			return true;
		}
		default:
			return false;
		}
	}

	/** @brief zero setting _Xn  f0 and _fi. */
	inline void t_grkf6::clearRkfSettings()
	{
		_Xn = _f0 = _fi[0] = _fi[1] = _fi[2] = _fi[3] = _fi[4] = _fi[5] = _fi[6] = _fi[7] = _fi[8] = 0.0;
	}

	/** @brief reset runge kutta data.
	* @param[in]   Xn		initial value of runge kutta integratio
	* @param[in]   t		begin time of runge kutta
	*/
	void t_grkf6::resetRkfPointValue(double Xn, t_gtime t)
	{
		// f_i
		_f0 = _fi[0] = _fi[1] = _fi[2] = _fi[3] = _fi[4] = _fi[5] = _fi[6] = _fi[7] = _fi[8] = 0.0;
		_Xn = Xn;
		_t = t;
	}

	/** @brief calculate integration result.
	* @return	 RKF integration result
	*/
	inline double t_grkf6::getRkfResult()
	{
		return _Xn + _h * _d0 * _f0 +
			_h * _d3 * _fi[2] +
			_h * _d4 * _fi[3] +
			_h * _d5 * _fi[4] +
			_h * _d6 * _fi[5] +
			_h * _d8 * _fi[7] +
			_h * _d9 * _fi[8];
	}
}

