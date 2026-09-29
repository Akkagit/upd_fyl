/**
*
* @verbatim
	History
	 -1.0 jdhuang	2019-01-04 creat the file.
	 -1.1 lijie		2019-04-02 add doxygen style code
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			grkfnav.CPP
* @brief		the class for satellite state integral by runge kutta.
*
* @author       jdhuang, Wuhan University
* @version		1.0.0
* @date			2019-01-04
*
*/
#include "grkfnav.h"

namespace great
{
	/** @brief default constructor. */
	t_grkf_nav::t_grkf_nav()
	{

	}

	/** @brief constructor. */
	t_grkf_nav::t_grkf_nav(int rkf_intetor_num, int rkf_degree, double step_h, t_gtime tn, const vector<double>& Xn)
	{
		_rkf_intetor_num = rkf_intetor_num;
		_rkf_degree = rkf_degree;
		_rkf_Xn = Xn;
		if (rkf_degree == 7)
		{
			for (int i = 0; i < rkf_intetor_num; ++i)
			{
				t_grkf* rkf = new t_grkf6(step_h,tn, Xn[i]);
				_rkf_list.push_back(rkf);
				_rkf_result.push_back(0.0);
			}
		}
	}

	/** @brief default destructor. */
	t_grkf_nav::~t_grkf_nav()
	{
		for (int i = 0; i < _rkf_intetor_num; ++i)
		{
			delete _rkf_list[i];
			_rkf_list[i] = NULL;
		}
	}

	/** @brief get runge kutta result. 
	* @return runge kutta result
	*/
	vector<double>& t_grkf_nav::getRkfResultVector()
	{
		double result = 0.0;
		_rkf_result.clear();
		for (int i = 0; i < _rkf_intetor_num; ++i)
		{
			result = _rkf_list[i]->getRkfResult();
			_rkf_result.push_back(result);
		}
		return _rkf_result;
	}

	/** @brief get initial value(Xn) of runge kutta.
	* @return initial value(Xn) of runge kutta
	*/
	vector<double>& t_grkf_nav::getRkfXnVector()
	{
		return _rkf_Xn;
	}

	/** @brief reset initial value and time. 
	* @param[in]   Xn	initial value 
	* @param[in]   t	begin time of rkf
	*/
	bool t_grkf_nav::resetRkf(const vector<double>& Xn,const t_gtime& t)
	{
		if (_rkf_list.size() < _rkf_intetor_num)
		{
			return false;
		}
		_rkf_Xn = Xn;
		for (int i = 0; i < _rkf_intetor_num; ++i)
		{
			_rkf_list[i]->resetRkfPointValue(Xn[i], t);
		}

		return true;
	}

	/** @brief reset initial value and degree.
	* @param[in]   Xn		initial value
	* @param[in]   degree	order of rkf
	*/
	void t_grkf_nav::resetRkf(const int& degree, const vector<double>& Xn)
	{
	}

	/** @brief zero setting runge kutta dat. */
	bool t_grkf_nav::clearRkf()
	{
		if (_rkf_result.size() < _rkf_intetor_num || 
			_rkf_Xn.size() < _rkf_intetor_num     || 
			_rkf_list.size() < _rkf_intetor_num)
		{
			return false;
		}

		for (int i = 0; i < _rkf_intetor_num; ++i)
		{
			_rkf_result[i] = _rkf_Xn[i] = 0.0;
			_rkf_list[i]->clearRkfSettings();
		}
		return true;
	}

	/** @brief setting show. */
	void t_grkf_nav::showRkfInfo()
	{
	}

	/** @brief update tn and Xn and storage fn for state. 
	* @param[in]   i	
	* @param[in]   fn	partial of X(i-1)
	* @param[out]  tn	ti
	* @param[out]  Xn	Xi
	*/
	bool t_grkf_nav::updatRkfPoint(const int& i, const vector<double>& fn, t_gtime& tn, vector<double>& Xn)
	{
		if (_rkf_list.size() < _rkf_intetor_num)
		{
			return false;
		}

		for (int n = 0; n < _rkf_intetor_num; ++n)
		{
			if (!_rkf_list[n]->updateRkfPoint(i, fn[n], tn, Xn[n]))
			{
				return false;
			}
		}
		return true;
	}

}

