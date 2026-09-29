/**
 * The poleut1 file  provides the pole and ut1-TAI information.\n
 * Following is an example for the header of poleut1 file.
 *
 * @verbatim
     %% conveop             zhangwei            07-Mar-19                        ! the author of this file and the time of creation
     %% Bulletin A file : finals2000A.all
     %% Input EOP  file : poleut1
     %% 
     +pole&ut1
     % UT1 type = UT1R
     % Start&End%Interval =    51544   58616   1.00                              ! the start time and end time(mjd) of the file,the interval(unit: day) of the record
     % Num. of Vars&Units =     5  0.1D+01  0.1D+01  0.1D+01  0.1D+01  0.1D+01   ! the number of the parameters in a line and the unit of them
     % Format = (f9.2,1x,2f10.6,f15.7,2f10.3,6(1x,a1))                           ! the format of the record(one line)
     %% MJD        XPOLE     YPOLE      UT1-TAI        DPSI     DEPSI    PRED_ID ! the information of the record(one line)
   @endverbatim
 *
 * An example for the record of the attitude file
 * @verbatim
   =============================================================================================
   !MJD        XPOLE     YPOLE       UT1-TAI         DPSI      DEPSI PRED_ID
   =============================================================================================
	51544.00   0.043240  0.377910    -31.6444900     0.000     0.000 I I I - - -
    51545.00   0.043510  0.377750    -31.6453930     0.000     0.000 I I I - - -
    51546.00   0.043620  0.377450    -31.6461550     0.000     0.000 I I I - - -
    .
	.
	.
	58614.00   0.095760  0.418733    -37.1681963     0.000     0.000 P P P - - -
    58615.00   0.097222  0.419217    -37.1690213     0.000     0.000 P P P - - -
    58616.00   0.098694  0.419680    -37.1698288     0.000     0.000 P P P - - -
   -pole&ut1                   ! end of the file
   @endverbatim
 *
 * @verbatim
	 History
	  -1.0	Yujie Qin    2019-01-06 creat the file.
	  -1.1  lijie	     2019-03-30 adding Doxygen Style Code Remarks
	  -1.2 Yujie Qin     2019-07-02 Adding Doxygen Style Code Remarks
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		    poleut1.h
 * @brief		    The base class used to decode poleut1 file information.
*
 * @author          Yujie Qin, Wuhan University
 * @version		    1.0.0
 * @date		    2019-03-30
 */

#ifndef POLEUT1_H
#define POLEUT1_H

#include "gutils/ExportBase.h"

#include "gutils/gtime.h"
#include "gcoders/gcoder.h"
#include "gdata/gpoleut1.h"

using namespace std;
using namespace gnut;

namespace great
{
	/**
	*@brief	   Class for decoding the poleut1 data
	*/
	class BASE_LIBRARY_EXPORT t_poleut1 : public t_gcoder
	{
	public:
		
		/**
		 * @brief default constructor.
		 *
		 * @param[in]  s        setbase control
		 * @param[in]  version  version of the gcoder
		 * @param[in]  sz       size of the buffer
		 */
		t_poleut1(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_poleut1() {};

		/**
		 * @brief decode the header of the poleut1 data file.
		 *
		 * The function is used for decoding the head of poleut1 file.\n
		 * pay attention to the buff and the size of buff which may cause some trouble when
			using the wrong value in decoding.
		 *
		 * @brief decode header of poleut1 file
		 * @param[in]  buff        buffer of the data
		 * @param[in]  sz          buffer size of the data
		 * @param[in]  errmsg      error message of the data decoding
		 * @return
			@retval >=0 consume size of header decoding
			@retval <0  finish reading
		 */
		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg);

		/**
		 * @brief decode the data body of the poleut1 data file.
		 *
		 * decode data body of poleut1 file, all the data read will store in the gpoleut1
		 *
		 * @param[in]  buff        buffer of the data
		 * @param[in]  sz          buffer size of the data
		 * @param[in]  errmsg      error message of the data decoding
		 * @return
			@retval >=0 consume size of body decoding
			@retval <0  finish reading
		 */
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);
		virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

	protected:
		string          _timetype;  ///< type of time(UT1R)
		int             _begtime;   ///< begin time(modified julian day)
		int             _endtime;   ///< end time(modified julian day)
		double          _interval;  ///< time interval(unit: day)
		int             _parnum;    ///< number of parameters
		vector<double>  _parunit;   ///< unit of each parameter
		vector<string>  _parname;   ///< name of each parameter
	};
} // namespace

#endif