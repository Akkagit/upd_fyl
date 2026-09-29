/**
*
* @verbatim
	History
	 -1.0 Yujie Qin	  2019-01-06 creat the file.
	 -1.1 Jie Li      2019-03-31 add doxygen style code
	 -1.2 Yujie Qin   2019-07-04 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			    gadams.h
* @brief		    adams integration
*
* @author		    Yujie Qin, Wuhan University
* @version		    1.0.0
* @date			    2019-01-06
*
*/

#ifndef GADAMS_H
#define GADAMS_H

#include "newmat/newmat.h"
#include "gutils/gconst.h"
#include "gutils/gtime.h"
#include <list>
namespace great
{
	/**
	*@brief	   Class for adams integration
	*/
	class t_adams
	{
	public:
		/** @brief default constructor. */
		t_adams();

		/**
		* @brief default constructor.
		* @param[in]  order_adams   adams integration order
		* @param[in]  h             adams integration step(unit: sec)
		*/
		t_adams(int order_adams, double h);

		/** @brief default destructor. */
		virtual ~t_adams();

		/**
	    * @brief reset adams integration order
	    * @param[in]  order   adams integration order
	    */
		void resetAdamsOrder(int order);

		/**
	    * @brief reset adams integration step
	    * @param[in]  h   adams integration step
	    */
		void resetAdamsH(double h);

		/**
	    * @brief calculate the new adams state and update state list.
	    * @param[in]  mode		adams integration mode
	    * @param[in]  partial	partial list
	    * @param[out] state		state list
	    * @return  bool
	    *	@retval   true    calculate and update successfully
	    *   @retval   false   calculate and update unsuccessfully
	    */
		bool update(const string& mode, list<double>& Yn_1, list<double>& fright);

	private:
		/** 
		* @brief calculate adams coefficient. 
		* @return  bool
	    *	@retval   true    calculate adams successfully
	    *   @retval   false   calculate adams unsuccessfully
	    */
		bool _calcAdamsCoeff();
		
		double _h;				///< adams step(unit: sec).
		string _mode;			///< model of adams integral.	
		int	_order_adams;		///< adams order.
		ColumnVector _b;		///< coefficient b.
		ColumnVector _betal;	///< coefficient β.
	};
}
#endif // !GNUMINTERG_H
