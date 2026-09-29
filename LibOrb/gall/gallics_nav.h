/**
*
* @verbatim
	History
	 -1.0 Hongmin Zhang		2019-01-06 creat the file.
	 -1.1 Jie Li		    2019-03-30 adding doxygen style code remarks
	 -1.2 Huang Jiande      2020-09-12 Standard code style
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gallics_nav.h
* @brief		Storage the GNS ICS data(more than one satellite)
*
* @author       Hongmin Zhang, Wuhan University
* @version		1.0.0
* @date		    2019-03-30
*
*/


#ifndef GALLICS_NAV_H
#define GALLICS_NAV_H


#include "gdata/gics_nav.h"
#include "gdata/gdata.h"

using namespace std;
using namespace gnut;

namespace great
{
	/**
	*@brief	   Class for all GNS ICS data storaging
	*/
	class OI_LIBRARY_EXPORT t_gallics_nav : public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gallics_nav();

		/** @brief default destructor. */
		virtual ~t_gallics_nav() {};

		/** 
		* @brief whether the ics data is empty.
		* @return  bool
			@retval   true    ics data is empty
		    @retval   false   ics data is existent
		*/
		bool isEmpty();

		/** 
		* @brief get ICS data for one GNS satellite.
		* @param [in]  sat_name     GNS satellite name
		* @return  ICS data for one GNS satellite
		*/
		t_gics_nav* get_nav_ics(std::string sat_name);
		
		/** 
		* @brief push ICS data for one GNS satellite into the data map.
		* @param[in]  sat_name   GND satellite name
		* @param[in]  gics		 ICS data for one GNS satellite
		*/
		void add_nav_ics(std::string sat_name, t_gics_nav gics);
		/**
		* @brief delete ICS data for one GNS satellite in the data map.
		* @param[in]  sat_name   GND satellite name
		*/
		void del_nav_ics(string sat_name);

		t_gtime beg_time;					///< begin time of the ICS data
		t_gtime end_time;					///< end time of the ICS data

		int                   sat_num;		///< number of satellites
		std::string           sys_type;		///< satellite type(LEO/GNS)
		std::set<std::string> sat_list;		///< PRN of satellites
		std::set<std::string> sys_list;		///< abbreviation of systems

		

	protected:
		std::map<std::string, t_gics_nav> _map_nav_ics;    ///< map for all ICS data
	};

}  //namespace

#endif