/**
 *
 * @verbatim
 	History
 	-1.0 Yujie Qin	  2019-01-06 creat the file.
 	-1.1 Wei Zhang    2019-03-27 Adding Doxygen Style Code Remarks
 	-1.2 Yujie Qin    2019-07-03 Adding Doxygen Style Code Remarks
 @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file			gattitude.cpp
 * @brief		    Storage the LEO attitude data for one satellite
 *
 * @author  	    Yujie Qin, Wuhan University
 * @version		    1.0.0
 * @date			2019-03-27
 *
 */


#ifndef GATTITUDE_H
#define GATTITUDE_H

#include <vector>
#include <map>

#include "gdata/gdata.h"
#include "gutils/gtime.h"
#include "gutils/ExportBase.h"

using namespace std;
using namespace gnut;

namespace great
{
	/**
	* @brief	   Class for one record of attitude data(quanternion parameters: q1, q2, q3, q4)
	*/
	class BASE_LIBRARY_EXPORT t_attitude_record
	{
	public:
		/** @brief defalut constructor. */
		t_attitude_record();

		/** @brief default destructor. */
		virtual ~t_attitude_record();

		double q1;  ///< q1 of quanternion
		double q2;  ///< q2 of quanternion
		double q3;  ///< q3 of quanternion
		double q4;  ///< q4 of quanternion
	};

	/**
	* @brief	   Class for attitue data storaging of one LEO satellite
	*/
	class BASE_LIBRARY_EXPORT t_gattitude : public t_gdata
	{
	public:
		/** @brief defalut constructor. */
		t_gattitude();

		/** @brief default destructor. */
		virtual ~t_gattitude();

		/**
	    * @brief push time and _onerec into the map.
	    * @param [in]  time        time of attitude data
	    * @param [in]  onerec      one record of attitude data
	    */
		void set_attitude_record(t_gtime time, t_attitude_record onerec);

		/**
	    * @brief set begin and end time of the attitude data.
	    * @param [in]  begtime       begin time of attitude file
	    * @param [in]  endtime       end time of attitude file
	    */
		void set_time(const t_gtime& begtime, const t_gtime& endtime);

		/**
	    * @brief set interval of the data.
	    * @param [in]  inter       interval of attitude data
	    */
		void set_inter(int inter);

		/**
	     * @brief return begin time of the attitude data
	     * @return    begin time
	     */
		t_gtime  beg();

		/**
	    * @brief return end time of the attitude data
	    * @return    end time
	    */
		t_gtime  end();

		/**
	    * @brief return _inter.
	    * @return   interval of the attitude data
	    */
		int      interval();

		/**
	    * @brief return attitude data map for one LEO satellite
	    * @return  attitude data map for one LEO satellite
	    */
		map<t_gtime, t_attitude_record>* attitude();

	protected:
		map<t_gtime, t_attitude_record>  _attitude;  ///< map of the epoch and the corresponding attitude record for one satellite
		t_gtime                     _beg;       ///< begin time of the attitude data
		t_gtime                     _end;       ///< end time of the attitude data
		int                         _inter;     ///< interval of the data
	};
};// namespace

#endif
