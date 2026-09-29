/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2020-01-05  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    ginverse.h
* @brief		header files of matrix inverse
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-17
*
*/

#ifdef USE_OPENBLAS
#ifndef GINVERSE_OPENBLAS_H
#define GINVERSE_OPENBLAS_H


#include "cblas.h"
#include "lapack.h"
#include <stdexcept>
#include "gproc/ginverse.h"


namespace great
{
	template <class MyMatrix>
	class t_ginverse_openblas:public t_ginverse<MyMatrix>
	{
	public:

		t_ginverse_openblas(bool is_from_zero) :
			t_ginverse<MyMatrix>(is_from_zero)
		{
		}

		/**
		* @brief solve the NEQ * X = W
		* @param[in] num the number of dimension
		* @param[in] NEQ dimension:num*num
		* @note result save into the dx
		*/
		virtual void sovle_NEQ(int num, MyMatrix& NEQ, const MyMatrix& W, MyMatrix& X)
		{
			char U = 'U';
			int Num = num;
			int One = 1;
			int error_info = 0;

			double* neq = new double[num*num];
			double* w = new double[num];

			// input
			for (int i = 0; i < num; i++) {
				for (int j = 0; j < num; j++) {
					neq[i*num + j] = NEQ(i + this->_beg, j + this->_beg);
				}
				w[i] = W(i + this->_beg,this->_beg);
			}

			
			// Cholosky 
			LAPACK_dpotrf(&U, &Num, neq, &Num, &error_info);
			if (error_info != 0)
			{
				throw std::runtime_error("Lapack Solve Wrong!");
			}

			// Solve X
			LAPACK_dpotrs(&U, &Num, &One, neq, &Num, w, &Num, &error_info);
			if (error_info != 0)
			{
				throw std::runtime_error("Lapack Solve Wrong!");
			}

			// Inverse NEQ
			LAPACK_dpotri(&U, &Num, neq, &Num, &error_info);
			if (error_info != 0)
			{
				throw std::runtime_error("Lapack Solve Wrong!");
			}

			// output
			for (int i = 0; i < num; i++) {
				for (int j = 0; j < num; j++) {
					NEQ(i + this->_beg, j + this->_beg) = neq[i*num + j];
				}
				X(i + this->_beg, this->_beg) = w[i];
			}

			delete[] neq;
			delete[] w;
		}

		virtual void solve_x(int num, const MyMatrix& NEQ, const MyMatrix& W, MyMatrix& X)
		{
			char U = 'U';
			int Num = num;
			int One = 1;
			int error_info = 0;

			double* neq = new double[num*num];
			double* w = new double[num];

			// input
			for (int i = 0; i < num; i++) {
				for (int j = 0; j < num; j++) {
					neq[i*num + j] = NEQ(i + this->_beg, j + this->_beg);
				}
				w[i] = W(i + this->_beg, this->_beg);
			}

			// Cholosky 
			LAPACK_dpotrf(&U, &Num, neq, &Num, &error_info);
			if (error_info != 0)
			{
				throw std::runtime_error("Lapack Solve Wrong!");
			}

			// Solve X
			LAPACK_dpotrs(&U, &Num, &One, neq, &Num, w, &Num, &error_info);
			if (error_info != 0)
			{
				throw std::runtime_error("Lapack Solve Wrong!");
			}

			// output
			for (int i = 0; i < num; i++) {
				X(i + this->_beg, this->_beg) = w[i];
			}

			delete[] neq;
			delete[] w;
		}

		virtual void inverse_NEQ(int num, MyMatrix& NEQ)
		{
			char U = 'U';
			int Num = num;
			int One = 1;
			int error_info = 0;

			double* neq = new double[num*num];

			// input
			for (int i = 0; i < num; i++) {
				for (int j = 0; j < num; j++) {
					neq[i*num + j] = NEQ(i + this->_beg, j + this->_beg);
				}
			}

			// Cholosky 
			LAPACK_dpotrf(&U, &Num, neq, &Num, &error_info);
			if (error_info != 0)
			{
				throw std::runtime_error("Lapack Solve Wrong!");
			}

			// Inverse NEQ
			LAPACK_dpotri(&U, &Num, neq, &Num, &error_info);
			if (error_info != 0)
			{
				throw std::runtime_error("Lapack Solve Wrong!");
			}

			// output
			for (int i = 0; i < num; i++) {
				for (int j = 0; j < num; j++) {
					NEQ(i + this->_beg, j + this->_beg) = neq[i*num + j];
				}
			}

			delete[] neq;
		}


	};
}



#endif // !GINVERSE_LAPACK_H

#endif