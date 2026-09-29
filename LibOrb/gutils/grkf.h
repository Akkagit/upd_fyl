/**
*
* @verbatim
	History
	 -1.0 Jiande Huang 	2019-01-04 creat the file.
	 -1.1 Jie Li		2019-03-31 add doxygen style code
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			grkf.h
* @brief		runge kutta integration
*
* @author       Jiande Huang, Wuhan University
* @version		1.0.0
* @date			2019-01-04
*
*/

#ifndef GRKF_H
#define GRKF_H
#include "gutils/gtime.h"

using namespace gnut;

namespace great
{
	/**
	*@brief	   Base class for runge kutta
	*/
	class t_grkf
	{
	public:
		/** 
		 * @brief default constructor. 
		 */
		t_grkf() {};

		/**
		* @brief constructor.
		* @param[in]  step_h			integration step of rkf
		* @param[in]  tn				initial time of rkf
		* @param[in]  Xn				initial X of rkf
		*/
		t_grkf(double step_h, t_gtime tn, double Xn);

		/** 
		 * @brief default destructor. 
		 */
		virtual ~t_grkf() {};

		/**
		 * @brief show degree and step of runge kutta. 
		 */
		void show();

		/** 
		 * @brief zero setting _Xn  f0 and _fi. 
		 */
		virtual void	clearRkfSettings() = 0;

		/**
		* @brief calculate integration result.
		* @return RKF integration result
		*/
		virtual double	getRkfResult() = 0;

		/** 
		 * @brief reset runge kutta data. 
		 */
		virtual void resetRkfPointValue(double Xn,t_gtime t) {};

		/** 
		 * @brief update ti and Xi and storage fi_1. 
		 */
		virtual bool updateRkfPoint(const int& i,const double& fi_1, t_gtime& ti, double & Xi) { return true; };
	protected:
		int		_degree = 0;	///< degree of runge kutta
		double	_h = 0.0;		///< step of runge kutta
		t_gtime	_t ;			///< begin time of runge kutta
		double	_Xn = 0.0;		///< initial value for integration
	};

	/**
	*@brief	   Class for runge kutta in order 6
	*/
	class t_grkf6 : public t_grkf
	{
	public :
		/** @brief default constructor. */
		t_grkf6();


		/**
		* @brief constructor.
		* @param[in]  step_h			integration step of rkf
		* @param[in]  tn				initial time of rkf
		* @param[in]  Xn				initial X of rkf
		*/
		t_grkf6(double step_h, t_gtime tn,double Xn);

		/** @brief default destructor. */
		virtual ~t_grkf6() {};

		/** 
		* @brief calculate ti and Xi and storage fi_1.
	    * @param[in]   i
	    * @param[in]   fi_1		partial of X(i-1)
	    * @param[out]  ti		ti
	    * @param[out]  Xi		Xi
		* @return  bool
		*	@retval   true    calculate successfully
		*   @retval   false   calculate unsuccessfully
		*/
		bool updateRkfPoint(const int& i, const double& fi_1, t_gtime& ti, double & Xi) override;

		/** @brief reset runge kutta data.
	    * @param[in]   Xn		initial value of runge kutta integratio
	    * @param[in]   t		begin time of runge kutta
	    */
		void resetRkfPointValue(double Xn,t_gtime t) override;

		/** 
		* @brief calculate integration result.
		* @return RKF integration result
		*/
		double getRkfResult() override;

		/** @brief zero setting _Xn  f0 and _fi. */
		void clearRkfSettings() override;

	protected:

		double _fi[9];						///< derivative in xi
		double _f0 = 0.0;					///< derivative in xn

		

		
		const double _a1 = 2.0 / 33.0;     ///< coef defined
		const double _a2 = 4.0 / 33.0;
		const double _a3 = 2.0 / 11.0;
		const double _a4 = 0.50;
		const double _a5 = 2.0 / 3.0;
		const double _a6 = 6.0 / 7.0;
		const double _a7 = 1.0;
		const double _a9 = 1.0;

		
		const double _b10 = 2.0 / 33.0;    ///< init beta
		const double _b30 = 1.0 / 22.0;
		const double _b40 = 43.0 / 64.0;
		const double _b50 = -2383.0 / 486.0;
		const double _b60 = 10077.0 / 4802.0;
		const double _b70 = -733.0 / 176.0;
		const double _b80 = 15.0 / 352.0;
		const double _b90 = -1833.0 / 352.0;
		const double _b21 = 4.0 / 33.0;
		const double _b32 = 3.0 / 22.0;
		const double _b42 = -165.0 / 64.0;
		const double _b52 = 1067.0 / 54.0;
		const double _b62 = -5643.0 / 686.0;
		const double _b72 = 141.0 / 8.0;
		const double _b92 = 141.0 / 8.0;
		const double _b43 = 77.0 / 32.0;
		const double _b53 = -26312.0 / 1701.0;
		const double _b63 = 116259.0 / 16807.0;
		const double _b73 = -335763.0 / 23296.0;
		const double _b83 = -5445.0 / 46592.0;
		const double _b93 = -51237.0 / 3584.0;
		const double _b54 = 2176.0 / 1701.0;
		const double _b64 = -6240.0 / 16807.0;
		const double _b74 = 216.0 / 77.0;
		const double _b84 = 18.0 / 77.0;
		const double _b94 = 18.0 / 7.0;
		const double _b65 = 1053.0 / 2401.0;
		const double _b75 = -4617.0 / 2816.0;
		const double _b85 = -1215.0 / 5632.0;
		const double _b95 = -729.0 / 512.0;
		const double _b76 = 7203.0 / 9152.0;
		const double _b86 = 1029.0 / 18304.0;
		const double _b96 = 1029.0 / 1408.0;
		const double _b98 = 1.0;

		
		const double _c0 = 77.0 / 1440.0;   ///< init c
		const double _c3 = 1771561.0 / 6289920.0;
		const double _c4 = 32.0 / 105.0;
		const double _c5 = 243.0 / 2560.0;
		const double _c6 = 16807.0 / 74880.0;
		const double _c7 = 11.0 / 270.0;
		
		
		const double _d0 = 11.0 / 864.0;   ///< init d
		const double _d3 = _c3;
		const double _d4 = _c4;
		const double _d5 = _c5;
		const double _d6 = _c6;
		const double _d8 = 11.0 / 270.0;
		const double _d9 = 11.0 / 270.0;

	};
}
#endif // !GRKF_H
