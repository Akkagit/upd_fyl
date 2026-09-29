/**
*
* @verbatim
	History
	 -1.0 Jiande Huang	2019-01-04 creat the file.
	 -1.1 Jie Li		2019-03-31 add doxygen style code
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			grkfnav.h
* @brief		the class for satellite state integral by runge kutta. 
*
* @author       Jiande Huang, Wuhan University
* @version		1.0.0
* @date			2019-01-04
*
*/

#ifndef GRKFNAV_H
#define GRKFNAV_H

#include "gutils/grkf.h"
#include <vector>

using namespace std;

namespace great
{
	/**
	*@brief	   Class for satellite state integral by runge kutta
	*/
	class t_grkf_nav
	{
	public:
		/** @brief constructor. */
		t_grkf_nav();

		/**
		* @brief constructor. 
		* @param[in]  rkf_intetor_num   point num of rkf
	    * @param[in]  rkf_degree        integration order of rkf
	    * @param[in]  step_h			integration step of rkf
		* @param[in]  tn				initial time of rkf
		* @param[in]  Xn				initial X of rkf
		*/
		t_grkf_nav(int rkf_intetor_num, int rkf_degree, double step_h, t_gtime tn, const vector<double>& Xn);

		/** @brief destructor. */
		virtual ~t_grkf_nav();

		//todo const£¿£¿
		/** 
		* @brief get runge kutta result.
		* @return runge kutta result
		*/
		vector<double>& getRkfResultVector();

		/** 
		* @brief get initial value(Xn) of runge kutta.
		* @return initial value(Xn) of runge kutta
		*/
		vector<double>& getRkfXnVector();
		
		/** 
		* @brief reset initial value and time.
		* @param[in]   Xn	initial value
		* @param[in]   t	begin time of rkf
		* @return   if reset successfully
		*	@retval   true   reset successfully
		*   @retval   false  reset unsuccessfully
		*/
		bool resetRkf(const vector<double>& Xn,const t_gtime& t);

		/** @brief reset initial value and degree.
		* @param[in]   Xn		initial value
		* @param[in]   degree	order of rkf
		* @return   if reset successfully
		*	@retval   true   reset successfully
		*   @retval   false  reset unsuccessfully
		*/
		void resetRkf(const int& degree, const vector <double>& Xn);

		/**
		* @brief zero setting runge kutta dat. 
		* @return   if clear successfully
		*	@retval   true   clear successfully
		*   @retval   false  clear unsuccessfully
		*/
		bool clearRkf();

		/** @brief setting show. */
		void showRkfInfo();

		/** @brief update tn and Xn and storage fn for state.
		* @param[in]   i	degree
		* @param[in]   fn	partial of X(i-1)
		* @param[out]  tn	ti
		* @param[out]  Xn	Xi
		* @return   if update successfully
		*	@retval   true   update successfully
		*   @retval   false  update unsuccessfully
		*/
		bool updatRkfPoint(const int& i,const vector<double>& fn, t_gtime& tn, vector<double>& Xn);

	protected:
		int  _rkf_intetor_num;			///< number of parameters need runge kutta 
		int  _rkf_degree;				///< runge kuttta degree
		vector<t_grkf*> _rkf_list;		///< all data need runge kutta
		vector<double> _rkf_result;		///< runge kutta integral result
		vector<double> _rkf_Xn;			///< runge kutta initial value
	};
}
#endif // !GRKFNAV_H
