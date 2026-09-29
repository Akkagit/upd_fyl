/**
 * The pannel file of LEO satellite provide the pannel information of LEO.\n
 * The pannel file does not have the header.
 *
 * An example for the record of the pannel file
 * @verbatim
  =============================================================================================
   ! name of LEO,'+'means the start of this LEO
  =============================================================================================
     +CHAMP panel information
     ##
     ## normal direction is in accelerator system, but attitude is spacecraft system   ! information about this file
     ## to inertial system. it shoulb be changed to that in SCF   
     ### -AREA- --NORMAL_DIRECTION_IN_SC_FIXED----- --RS_RD-- -DRAG- -RSi_RDi-  ------NAME_OF_PANEL---------
     ##  ACC         x            y           z
     cc  SCF        -z            x          -y
     # 1 1.2920  1.00000000  0.00000000  0.00000000 0.05 0.30  1.000 0.03 0.16  Top SiO2
     # 2 3.6239 -1.00000000  0.00000000  0.00000000 0.68 0.20  1.000 0.19 0.06  Bottom SSM/Teflon(C2F4)n
     # 3 3.1593  0.70700000  0.00000000  0.70700000 0.05 0.30  1.000 0.03 0.16  Left S.P. SiO2
     # 4 0.3020 -0.70700000  0.00000000 -0.70700000 0.40 0.26  1.000 0.23 0.15  Left S.P. back MLI Kapt/SiO2
     # 5 3.1593  0.70700000  0.00000000 -0.70700000 0.05 0.30  1.000 0.03 0.16  Right S.P. SiO2
     # 6 0.3020 -0.70700000  0.00000000  0.70700000 0.40 0.26  1.000 0.23 0.15  Right S.P. back MLI Kapt/SiO2
     # 7 0.4902 -0.34200000 -0.93970000  0.00000000 0.40 0.26  1.000 0.23 0.15  Aft Panel MLI Kapton/SiO2
     # 8 1.2199 -0.93970000  0.34200000  0.00000000 0.20 0.40  1.000 0.26 0.51  Front Panel Sandblasted Al
     # 9 0.5110  1.00000000  0.00000000  0.00000000 0.40 0.26  1.000 0.23 0.15  Top Boom MLI Kapton/SiO2
     #10 0.4190  1.00000000  0.00000000  0.00000000 0.64 0.20  1.000 0.71 0.22  Top Boom Kapton/Al
     #11 0.9300 -1.00000000  0.00000000  0.00000000 0.40 0.26  1.000 0.23 0.15  Bottom Boom MLI Kapton/SiO2
     #12 0.9300  0.00000000  0.00000000 -1.00000000 0.40 0.26  1.000 0.23 0.15  Left Boom MLI Kapton/SiO2
     #13 0.9300  0.00000000  0.00000000  1.00000000 0.40 0.26  1.000 0.23 0.15  Right Boom MLI Kapton/SiO2
     #14 0.0529  0.00000000  1.00000000  0.00000000 0.40 0.26  1.000 0.23 0.15  Front Boom MLI Kapton/SiO2
       1 1.2920  0.00000000  0.00000000 -1.00000000 0.05 0.30  1.000 0.03 0.16  Top SiO2
       2 3.6239  0.00000000  0.00000000  1.00000000 0.68 0.20  1.000 0.19 0.06  Bottom SSM/Teflon(C2F4)n
       3 3.1593  0.00000000 -0.70700000 -0.70700000 0.05 0.30  1.000 0.03 0.16  Left S.P. SiO2
       4 0.3020  0.00000000  0.70700000  0.70700000 0.40 0.26  1.000 0.23 0.15  Left S.P. back MLI Kapt/SiO2
       5 3.1593  0.00000000  0.70700000 -0.70700000 0.05 0.30  1.000 0.03 0.16  Right S.P. SiO2
       6 0.3020  0.00000000 -0.70700000  0.70700000 0.40 0.26  1.000 0.23 0.15  Right S.P. back MLI Kapt/SiO2
       7 0.4902 -0.93970000  0.00000000  0.34200000 0.40 0.26  1.000 0.23 0.15  Aft Panel MLI Kapton/SiO2
       8 1.2199  0.34200000  0.00000000  0.93970000 0.20 0.40  1.000 0.26 0.51  Front Panel Sandblasted Al
       9 0.5110  0.00000000  0.00000000 -1.00000000 0.40 0.26  1.000 0.23 0.15  Top Boom MLI Kapton/SiO2
      10 0.4190  0.00000000  0.00000000 -1.00000000 0.64 0.20  1.000 0.71 0.22  Top Boom Kapton/Al
      11 0.9300  0.00000000  0.00000000  1.00000000 0.40 0.26  1.000 0.23 0.15  Bottom Boom MLI Kapton/SiO2
      12 0.9300  0.00000000  1.00000000  0.00000000 0.40 0.26  1.000 0.23 0.15  Left Boom MLI Kapton/SiO2
      13 0.9300  0.00000000 -1.00000000  0.00000000 0.40 0.26  1.000 0.23 0.15  Right Boom MLI Kapton/SiO2
      14 0.0529  1.00000000  0.00000000  0.00000000 0.40 0.26  1.000 0.23 0.15  Front Boom MLI Kapton/SiO2
  =============================================================================================
    ! name of LEO,'-'means the end of this LEO
  =============================================================================================
	  -CHAMP panel information 
	  ## Rs  : specular reflectivity coefficient in the visible domain
      ## Rd  : diffuse reflectivity coefficient in the visible domain
      ## Rsi : specular reflectivity coefficient in the infrared domain
      ## Rdi : diffuse  reflectivity coefficient in the infrared domain
   @endverbatim
 *
 * @verbatim
     History
      -1.0 Yujie Qin     2019-01-06 creat the file.
      -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
	  -1.2 Yujie Qin     2019-07-02 Adding Doxygen Style Code Remarks
	  -1.3 Jiande Huang  2020-09-12 Standard code style and Optimize code logic
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file			pannel.h
 * @brief		    Decode the LEO pannel data
 *
 * @author		    Yujie Qin, Wuhan University
 * @version		    1.0.0
 * @date			2019-03-27
 *
 */

#ifndef PANNEL_H
#define PANNEL_H

#include <string.h>
#include <stdio.h>
#include <iostream>
#include "gcoders/gcoder.h"
#include "gdata/gpannel.h"
#include "gall/gallpannel.h"

#include "gutils/ExportOI.h"

using namespace std;

namespace great
{
	/**
	*@brief	   Class for decoding the LEO pannel data
	*/
	class OI_LIBRARY_EXPORT t_pannel : public t_gcoder
	{
	public:
		
		/**
		 * @brief default constructor.
		 *
		 * @param[in]  s        setbase control
		 * @param[in]  version  version of the gcoder
		 * @param[in]  sz       size of the buffer
		 */
		t_pannel(t_gsetbase* s, string version, int sz);

		/** @brief default destructor. */
		virtual ~t_pannel();

		/**
		 * @brief decode the header of the LEO pannel data file.
		 *
		 * The function is used for decoding the head of pannel file.\n
		 * pay attention to the buff and the size of buff which may cause some trouble when
			using the wrong value in decoding.
		 *
		 * @brief decode header of pannel file
		 * @param[in]  buff        buffer of the data
		 * @param[in]  sz          buffer size of the data
		 * @param[in]  errmsg      error message of the data decoding
		 * @return
			@retval >=0 consume size of header decoding
			@retval <0  finish reading
		 */
		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg);

		/**
		 * @brief decode the data body of the pannel data file.
		 *
		 * decode data body of pannel file, all the data read will store in the gallpannel
		 *
		 * @param[in]  buff        buffer of the data
		 * @param[in]  sz          buffer size of the data
		 * @param[in]  errmsg      error message of the data decoding
		 * @return
			@retval >=0 consume size of body decoding
			@retval <0  finish reading
		 */
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

	protected:
		string    _sat_name;  ///< name of satellte
		t_gpannel _sat_pannel;  ///< all data

	};
} // namespace

#endif