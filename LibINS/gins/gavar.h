#ifndef GAVAR_H
#define GAVAR_H

/**
* @file gavar.h
* @details 
* @note Class for counting sequence variance

* @verbatim
		History
		-0.1	ShenZhiheng	       2018-10-15 creat the file.
		-0.3	ShenZhiheng	       2019-03-20 add some base functions
		-0.5	ShenZhiheng	       2019-06-10 rename class name
		-0.7	ShenZhiheng	       2019-12-02 restruct the file
  @endverbatim

* @author ShenZhiheng
* @version 0.7
* @date 2018-10-15
* @license
*/


#include "gutility.h"
#include "gutils/ExportINS.h"

using namespace std;

namespace gsins
{
	/**
	* @class t_gavar
	* @brief t_gavar Class for counting sequence variance
	*
	* Control result by attenuation time factor tau
	*/
	class INS_LIBRARY_EXPORT t_gavar
	{
	public:

		/** @brief default constructor. */
		t_gavar();


		/**
		* @brief constructor.
		*
		* @param[in]  n     number
		* @param[in]  ts    interval
		* @param[in]  r0    initial	standard deviation
		* @param[in]  tau   decay time
		*/
		t_gavar(int n, double ts, const Eigen::VectorXd& r0, const Eigen::VectorXd& tau);


		/**
		* @brief Init
		*
		* @param[in]  n     number
		* @param[in]  ts    interval
		* @param[in]  r0    initial	standard deviation
		* @param[in]  tau   decay time
		*/
		void init(int n, double ts, const Eigen::VectorXd& r0, const Eigen::VectorXd& tau);


		/**
		* @brief update variance
		* @param[in] r    value of type
		*/
		void update(const Eigen::VectorXd& r);

		/**
		* @brief get variance by i
		* @param[in] i			index
		* @return
		The value of index
		*/
		double operator()(int i)const;

	public:
		int n;			/// number of type
		double ts;		/// interval
		Eigen::VectorXd flag, R, Rmax, Rmin, tstau, r0;

	};


}



#endif