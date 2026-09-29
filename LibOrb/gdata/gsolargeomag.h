/**
*
* @verbatim
	History
	 -1.0 Yujie Qin		 2019-01-06 creat the file.
	 -1.1 Wei Zhang      2019-03-27 Adding Doxygen Style Code Remarks
	 -1.2 Yujie Qin      2019-07-03 Adding Doxygen Style Code Remarks
	 -1.3 Jiande Huang   2020-09-30 Standard code style and Optimize code logic
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gsolargeomag.h
* @brief		Storage the LEO solar_flux and geomag_kp data for one satellite
*
* @author       Yujie Qin, Wuhan University
* @version		1.0.0
* @date		    2019-03-27
*
*/

#ifndef GSOLARGEOMAG_H
#define GSOLARGEOMAG_H


#include <vector>
#include <map>
#include <stdio.h>
#include <string>
#include "gdata/gdata.h"
#include "gutils/gtime.h"
#include "gutils/ExportOI.h"

using namespace std;
using namespace gnut;
namespace great
{
	/**
	*@brief	   Class for solar_flux/geomag_kp data
	*/
	class OI_LIBRARY_EXPORT t_gsolargeomag : public t_gdata {

	public:
		/** @brief default constructor. */
		t_gsolargeomag();

		/** @brief default destructor. */
		virtual ~t_gsolargeomag() {};

		/** @brief the map stored the record in solar_flux file. */
		typedef map<double, double[2]> t_mapvalue;

		/**
		* @brief set one record with time and corresponding value(solar_flux/geomag_kp).
		* @param[in]    datatype      data type(solar/geomag)
		* @param[in]    mjdtime       modified juliay day
		* @param[in]    value         value in the requested time
		*/
		void set_data(string data_type, double mjd_time, double* value);

		/**
		* @brief set the interval of data(solar_flux/geomag_kp).
		* @param[in]  datatype   data type(solar_flux/geomag_kp)
		* @param[in]  interval   data interval
		*/
		void set_intv(string data_type, double interval);

		/**
		* @brief set begin and end time of data(solar_flux/geomag_kp)
		* @param[in]   datatype  data type(solar/geomag)
		* @param[in]   begin     begin time
		* @param[in]   end       end time
		*/
		void set_time(string data_type, double beg_time, double end_time);

		/**
		* @brief whether the solar_flux and geomag_kp data is empty.
		* @return  bool
		*	@retval   true    solar_flux and geomag_kp data is empty
		*   @retval   false   solar_flux and geomag_kp data is existent
		*/
		bool is_empty();

		/**
		* @brief return solar_flux data.
		* @return  solar_flux data
		*/
		t_mapvalue& get_solar_data() { return _solar_data; };

		/**
		* @brief return geomag_kp data.
		* @return  geomag_kp data
		*/
		map<double, double>& get_geomag_data() { return _geomag_gdata; };//get kp

		/**
		* @brief return interval according to the file type.
		* @param[in]   datatype  data type(solar_flux/geomag_kp)
		* @return      interval of the data
		*/
		double get_intv(string data_type);

		/**
		* @brief return begin time according to the file type.
		* @param[in]   datatype  data type(solar_flux/geomag_kp)
		* @return      begin time.
		*/
		double get_beg(string data_type);

		/**
		* @brief return end time according to the file type.
		* @param[in]   datatype  data type(solar_flux/geomag_kp)
		* @return      end time.
		*/
		double get_end(string data_type);

	protected:
		t_mapvalue _solar_data;                     ///< solar_flux data
		double     _solar_beg_time;                 ///< begin time of solar_flux data
		double     _solar_end_time;                 ///< end time of solar_flux data
		double     _solar_itv_time;                 ///< interval of solar_flux data

		map<double, double>  _geomag_gdata;          ///< geomag_kp data
		double    _geomag_beg_time;                  ///< begin time of geomag_kp data
		double    _geomag_end_time;                  ///< end time of geomag_kp data
		double    _geomag_itv_time;                  ///< interval of geomag_kp data

	};
}  //namespace

#endif