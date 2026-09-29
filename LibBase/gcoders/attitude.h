/**
 * The attitude file of LEO satellite provides the attitude information of LEO.\n
 * Following is an example for the header of attitude file.
 *
 * @verbatim
     %% Header of attitude data for LEO satellite     ! Start of header
	 % Satellite      SWARM-A                         ! Name of LEO
     % Start time     58178      0.0                  ! Start time of attitude file
     % End time       58180  86399.0                  ! End time of attitude file
     % Time interval  1                               ! The interval of record(unit:second)
     %% End of Header                                 ! End of header
	@endverbatim
 *
 * An example for the record of the attitude file
 * @verbatim
   =============================================================================================
   ! mjd   sod         q1          q2          q3          q4           flag     
   =============================================================================================
     58178     0.00000 -0.40817400 -0.44294000 -0.51730900 -0.60793800  0
     58178     1.00000 -0.40848500 -0.44326400 -0.51710000 -0.60767200  0
     58178     2.00000 -0.40879500 -0.44358700 -0.51689000 -0.60740500  0
     58178     3.00000 -0.40910500 -0.44391100 -0.51668100 -0.60713900  0
     58178     4.00000 -0.40941500 -0.44423400 -0.51647100 -0.60687200  0
     58178     5.00000 -0.40972500 -0.44455700 -0.51626100 -0.60660500  0
	 .
	 .
	 .
	 58180 86395.00000  0.68370000  0.71843900  0.06506610  0.11029800  0
     58180 86396.00000  0.68377900  0.71846100  0.06469470  0.10988200  0
     58180 86397.00000  0.68385800  0.71848300  0.06432300  0.10946600  0
     58180 86398.00000  0.68393600  0.71850500  0.06395110  0.10905000  0
     58180 86399.00000  0.68401400  0.71852700  0.06357900  0.10863400  0
   @endverbatim
 *
 * @verbatim
     History
      -1.0 Yujie Qin     2019-01-06 creat the file.
      -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
	  -1.2 Yujie Qin     2019-07-02 Adding Doxygen Style Code Remarks
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file			attitude.h
 * @brief	        Decode the LEO attitude data
 *
 * @author		    Yujie Qin, Wuhan University
 * @version		    1.0.0
 * @date			2019-03-27
 *
 */


#ifndef ATTITUDE_H
#define ATTITUDE_H

#include "gutils/ExportBase.h"

#include "gcoders/gcoder.h"
#include "gdata/gattitude.h"
#include "gutils/gtime.h"

namespace great
{
	/**
	*@brief	   Class for decoding the attitude data
	*/
	class BASE_LIBRARY_EXPORT t_attitude : public t_gcoder
	{
	public:

		/**
		 * @brief default constructor.
		 *
		 * @param[in]  s        setbase control
		 * @param[in]  version  version of the gcoder
		 * @param[in]  sz       size of the buffer
		 */
		t_attitude(t_gsetbase* s, std::string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_attitude();

		/**
		 * @brief decode the header of the attitude data file.
		 *
		 * The function is used for decoding the head of attitude file.\n
		 * pay attention to the buff and the size of buff which may cause some trouble when
			using the wrong value in decoding.
		 *
		 * @brief decode header of attitude file
		 * @param[in]  buff        buffer of the data
		 * @param[in]  sz          buffer size of the data
		 * @param[in]  errmsg      error message of the data decoding
		 * @return
			@retval >=0 consume size of header decoding
			@retval <0  finish reading
		 */
		virtual  int decode_head(char* buff, int bufLen, std::vector<std::string>& errmsg);
		
		/**
		 * @brief decode the data body of the attitude data file.
		 *
		 * decode data body of attitude file, all the data read will store in the gallattitude
		 *
		 * @param[in]  buff        buffer of the data
		 * @param[in]  sz          buffer size of the data
		 * @param[in]  errmsg      error message of the data decoding
		 * @return
			@retval >=0 consume size of body decoding
			@retval <0  finish reading
		 */
		virtual  int decode_data(char* buff, int bufLen, int& cnt, std::vector<std::string>& errmsg);

	protected:
		int         _interval_time;    ///< interval of the data
		t_gtime     _beg_time;     ///< begin time of the attitude data
		t_gtime     _end_time;     ///< end time of the attitude sata
		t_gattitude _attitude;   ///< att data
		std::string _sat_name;     ///< satellite name

	};
} // namespace

#endif