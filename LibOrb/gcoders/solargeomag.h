/**
 * The solar_flux and geomag_kp file of LEO satellite provide the solar radiation flux information and geomagnetic kp information.\n
 * Following is an example for the header of solar_flux file.
 *
 * @verbatim
	 ##Header of solar_flux file     ! Start of header
	 %Time       58119  58238  1.0   ! Start time,end time and interval (unit:day)of solar_flux file
	 %Parnum     2                   ! par number of a parmeter-group
	 %Parpernum  4                   ! number of parmeter-group in one line
	 %Parunit    1.0  1.0            ! unit of the par in a parmeter-group
	 ##End of header                 ! End of header
	@endverbatim
 *
 * An example for the record of the solar_flux file
 * @verbatim
   =============================================================================================
   ! mjd   par11 par12 par21 par22 par31 par32 par41 par42
   =============================================================================================
	58119  69.1   0.0  69.5   0.0  70.7   0.0  69.5   0.0
	58123  69.3   0.0  69.4   0.0  69.9   0.0  70.3   0.0
	58127  70.8   0.0  70.4   0.0  70.8   0.0  70.9   0.0
	.
	.
	.
	58227  70.8   0.0  73.0   0.0  76.8   0.0  75.7   0.0
	58231  73.9   0.0  72.9   0.0  70.8   0.0  69.4   0.0
	58235  68.7   0.0  70.2   0.0  71.1   0.0  70.2   0.0
   @endverbatim

* Following is an example for the header of geomag_kp file.
 *
 * @verbatim
	 ##Header of geomag_kp file        ! Start of header
	 %Time       58119  58238  0.125   ! Start time,end time and interval (unit:day)of geomag_kp file
	 %Parnum     1                     ! par number of a parmeter-group
	 %Parpernum  8                     ! number of parmeter-group in one line
	 %Parunit    1.0                   ! unit of the par in a parmeter-group
	 ##End of header                   ! End of header
	@endverbatim
 *
 * An example for the record of the geomag_kp file
 * @verbatim
   =============================================================================================
   ! mjd   par11...                      ...par81
   =============================================================================================
	58119  3.3  3.7  2.3  2.3  2.7  1.0  1.0  1.3
	58120  1.7  0.7  0.3  0.7  1.7  1.0  2.0  1.0
	58121  0.0  0.7  1.0  0.7  0.3  1.0  0.3  0.0
	.
	.
	.
	58236  1.7  0.7  0.7  0.3  1.3  0.7  0.3  0.0
	58237  0.3  0.0  1.0  0.7  1.3  1.3  0.3  0.7
	58238  2.3  1.3  0.7  0.7  2.0  1.7  1.3  0.0
   @endverbatim


 *
 * @verbatim
 History
 -1.0 Yujie Qin	    2019-01-06 creat the file.
 -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
 -1.2 Jie Li        2019-04-02 Add abnormal throwing mechanism of functions
 -1.3 Yujie Qin     2019-07-02 Adding Doxygen Style Code Remarks
 -1.4 Jiande Huang  2020-09-12 Standard code style and Optimize code logic
 @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file			solargeomag.cpp
 * @brief		    Decode the LEO solar_flux and geomag_kp data
 *
 * @author		    Yujie Qin, Wuhan University
 * @version		    1.0.0
 * @date			2019-03-27
 *
 */


#ifndef SOLARGEOMAG_H
#define SOLARGEOMAG_H

#include <vector>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include "gcoders/gcoder.h"
#include "gutils/gtime.h"
#include "gutils/ExportOI.h"

using namespace std;
using namespace gnut;

namespace great
{
	/**
	*@brief	   Class for decoding the solar_flux and geomag_kp data
	*/
	class OI_LIBRARY_EXPORT t_solargeomag : public t_gcoder {

	public:

		/**
		 * @brief default constructor.
		 *
		 * @param[in]  s        setbase control
		 * @param[in]  version  version of the gcoder
		 * @param[in]  sz       size of the buffer
		 */
		t_solargeomag(t_gsetbase* s, string version, int sz);

		/** @brief default constructor. */
		t_solargeomag(string file_type);

		/** @brief default destructor. */
		virtual ~t_solargeomag() {};

		/**
		 * @brief decode the header of the solar_flux/geomag_kp data file.
		 *
		 * The function is used for decoding the head of the solar_flux/geomag_kp data file.\n
		 * pay attention to the buff and the size of buff which may cause some trouble when
			using the wrong value in decoding.
		 *
		 * @brief decode header of the solar_flux/geomag_kp data file
		 * @param[in]  buff        buffer of the data
		 * @param[in]  sz          buffer size of the data
		 * @param[in]  errmsg      error message of the data decoding
		 * @return
			@retval >=0 consume size of header decoding
			@retval <0  finish reading
		 */
		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg);

		/**
		 * @brief decode the data body of the solar_flux/geomag_kp data file.
		 *
		 * decode data body of solar_flux/geomag_kp file, all the data read will store in the gsolargeomag
		 *
		 * @param[in]  buff        buffer of the data
		 * @param[in]  sz          buffer size of the data
		 * @param[in]  errmsg      error message of the data decoding
		 * @return
			@retval >=0 consume size of body decoding
			@retval <0  finish reading
		 */
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

	private:
		string         _file_type;   ///< file type of the input file(solar_flux or geomag_kp)
		double         _beg_mjd;     ///< begin time(mjd) of the data for whole file 
		double         _end_mjd;     ///< end time(mjd) of the data for whole file 
		double         _interval;    ///< interval of the data 
		int            _par_type_num;  ///< number of parameter type
		int            _par_line_num;  ///< number of values in one line for one type parameter
		vector<double> _vunit;         ///< unit for each parameter type
	};
}  //namespace

#endif