/**
*
* @verbatim
    History
     -1.0 Hongmin Zhang 2019-01-06 creat the file.
     -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
	 -1.2 Yujie Qin     2019-03-27 Adding Doxygen Style Code Remarks
	 -1.3 Huang Jiande  2020-09-12 Standard code style
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gallics_leo.h
* @brief		Storage the LEO ICS data(more than one satellite)
*
* @author       Hongmin Zhang, Wuhan University
* @version		1.0.0
* @date		    2019-03-27
*
*/

#ifndef GALLICSLEO_H
#define GALLICSLEO_H

#include "gdata/gics_leo.h"
#include "gdata/gdata.h"

using namespace std;
using namespace gnut;

namespace great
{
	/**
	*@brief	   Class for all LEO ICS data storaging
	*/
	class OI_LIBRARY_EXPORT t_gallics_leo : public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gallics_leo();

		/** @brief default destructor. */
		virtual ~t_gallics_leo();

		/* -------------
	     * @brief get icsleo for one satellite
	     * @param [in]  sat_name     LEO satellite name
	     * @return  ICS data for one LEO satellite
	     */
		t_gics_leo* get_leo_ics(string sat_name);

		/* ----------------------------
	     * @brief whether the ics data is empty.
	     * @return  bool
	     *	@retval   true    ics data is empty
	     *   @retval   false   ics data is existent
	     */
		bool isEmpty();

		/* ----------------------------
	     * @brief push one satellite with ics into the map
	     * @param[in]  sat_name   LEO satellite name
	     * @param[in]  gicsleo    ICS data for one LEO satellite
	     */
		void add_leo_ics(string sat_name, t_gics_leo& gics_leo);

		/* ----------------------------
	     * @brief set begin and end time of the ICS data.
	     * @param[in]  begtime   begin time
	     * @param[in]  endtime   endtime
	     */
		void set_time(const t_gtime& begtime, const t_gtime& endtime);
		void set_beg_time(const t_gtime& beg_time);
		void set_end_time(const t_gtime& end_time);

		/* ----------------------------
	     * @brief return begin time of the ICS data.
	     * @return    begin time
	     */
		t_gtime get_beg_time();

		/* ----------------------------
	     * @brief return end time of the ICS data.
	     * @return    end time
	     */
		t_gtime get_end_time(); 

		/* ----------------------------
	     * @brief return LEO satellite name list.
	     * @return    name list
	     */
		std::set<std::string>	get_leo_sat_list();	

		/* ----------------------------
	     * @brief return number of satellites.
	     * @return    satellite number
	     */
		int get_leo_sat_num();

	protected:
		std::map<std::string, t_gics_leo> _map_ics_leo;   ///< map for all ICS data
		t_gtime _beg_time;                                ///< begin time of the ICS data
		t_gtime _end_time;                                ///< end time of the ICS data
	};
} // namespace

#endif