/**
*
* @verbatim
	History
	 -1.0 Yujie Qin	 2019-01-06 creat the file.
	 -1.1 Jie Li     2019-04-02 add doxygen style code
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gadams.cpp
* @brief		adams integration
*
* @author		Yujie Qin, Wuhan University
* @version		1.0.0
* @date			2019-01-06
*
*/

#include "gutils/gadams.h"
#include <iostream>
#include <iomanip>
#include <newmat/newmatio.h>

using namespace ::std;
namespace great
{
	/** @brief default constructor. */
	t_adams::t_adams()
	{

	}

	/** @brief constructor. */
	t_adams::t_adams(int order_adams, double h)
	{
		_order_adams = order_adams;
		_h = h;
		_calcAdamsCoeff();
	}

	/**
	* @brief reset adams integration order
	* @param[in]  order   adams integration order
	*/
	void t_adams::resetAdamsOrder(int order)
	{
		_order_adams = order;
	}

	/**
	* @brief reset adams integration step
	* @param[in]  h   adams integration step
	*/
	void t_adams::resetAdamsH(double h)
	{
		_h = h;
	}

	/** @brief destructor. */
	t_adams::~t_adams()
	{

	}

	/**
	* @brief calculate the new adams state and update state list.
	* @param[in]  mode		adams integration mode
	* @param[in]  partial	partial list
	* @param[out] state		state list
	*/
	bool t_adams::update(const string& mode, list<double>& state, list<double>& partial)
	{
		
		int i = 0;
		double sum = 0.0;
		
		try
		{
			if (mode == "PRE")
			{
				list<double>::reverse_iterator ritlist = partial.rbegin();
				list<double>::reverse_iterator rend = partial.rend();
				for (ritlist; ritlist != rend; ++ritlist)
				{
					sum += (*ritlist)*_b(i + 1);
					i++;
				}
				double tmp = state.back() + _h * sum;
				state.push_front(tmp);
				return true;
			}
			else if (mode == "COR")
			{
				list<double>::reverse_iterator ritlist = partial.rbegin();
				list<double>::reverse_iterator rend = partial.rend();
				for (ritlist; ritlist != rend; ++ritlist)
				{
					sum += (*ritlist)*_betal(i + 1);
					i++;
				}
				double tmp = state.back() + _h * sum;
				state.push_back(tmp);
				state.pop_front(); //delete pre adams state of 12th
				state.pop_front(); //delete the first adams state
				return true;
			}
			return false;
		}
		catch (...)
		{
			return false;
		}
	}

	/**
	* @brief calculate adams coefficient.
	*/
	bool t_adams::_calcAdamsCoeff()
	{
		//ColumnVector factor, gama, gamal;
		//bino.ReSize(_order_adams + 2, _order_adams + 2);
		//factor.resize(_order_adams + 2);
		//gama.resize(_order_adams + 2);
		//gamal.resize(_order_adams + 2);

		double s = 1.0;

		try
		{
			ColumnVector factor(0.0, _order_adams + 2);
			ColumnVector gama(0.0, _order_adams + 2);
			ColumnVector gamal(0.0, _order_adams + 2);
			Matrix bino(0.0, _order_adams + 2, _order_adams + 2);
			_b.resize(_order_adams);
			_betal.resize(_order_adams);
			factor(1) = 1.0;

			for (int i = 2; i <= _order_adams + 1; i++)
			{
				factor(i) = factor(i - 1)*(i - 1);
			}

			for (int i = 1; i <= _order_adams + 1; i++)
			{
				for (int j = 1; j <= i; j++)
				{
					bino(i, j) = factor(i) / factor(j) / factor(i - j + 1);
				}
			}

			for (int i = 1; i <= _order_adams + 1; i++)
			{
				gama(i) = 1.0;
				gamal(i) = 0.0;
				for (int j = 1; j < i; j++)
				{
					gama(i) = gama(i) - (gama(j) / (i - j + 1.0));
					gamal(i) = gamal(i) - (gamal(j) / (i - j + 1.0));
				}
				if (i == 1)
				{
					gamal(i) = 1.0;
				}
			}

			for (int i = 1; i <= _order_adams; i++)
			{
				_b(i) = 0.0;
				_betal(i) = 0.0;

				for (int j = i; j <= _order_adams; j++)
				{
					_b(i) = _b(i) + bino(j, i)*gama(j);
					_betal(i) = _betal(i) + bino(j, i)*gamal(j);
				}
				_b(i) = _b(i) * s;
				_betal(i) = _betal(i) * s;
				s = -s;
			}

			return true;
		}
		catch (...)
		{
			return false;
		}
	}

}