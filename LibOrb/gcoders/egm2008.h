/**
* The EGM file is used for calculating gravity.\n
* Following is an example for the header of EGM file.
*
* @verbatim
	 FIRST GSM-2_0110_2002214-2003242_EIGEN_G---_0001 SHM     GFZ POTSDAM  20040213
	 CMMNT GRACE-only mean gravity field
	 CMMNT Reported standard deviations are formal (not calibrated)
	 EARTH 0.3986004415E+15 0.6378136460E+07  19970101							! GM(*1E-9 km/s) radius(unit: m) and refrence time of EARTH
	 SHM     150  150 1.00 fully normalized exclusive permanent tide			! max order of the coefficient
   @endverbatim
*
* An example for the record of the EGM file
* @verbatim
  =============================================================================================
  n    m		 c						s			sigma c		sigma s
  =============================================================================================
  0    0  0.100000000000E+01  0.000000000000E+00  0.0000E+00  0.0000E+00  20020802.1200  20030830.1200  nnnn
  @endverbatim
*
* @verbatim
	History
	 -1.0 Jie Li	   2019-01-06 creat the file.
	 -1.1 Jie Li	   2019-03-30 adding doxygen style code remarks
	 -1.2 Jiande Huang 2020-09-12 Standard code style
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			egm2008.h
* @brief		decode the EGM data
*
* @author		Jie Li, Wuhan University
* @version		1.0.0
* @date			2019-03-30
*
*/

#ifndef EGM2008_H
#define EGM2008_H

#include "gcoders/gcoder.h"
#include "gutils/ExportOI.h"

using namespace gnut;

namespace great
{
	/**
	*@brief	   Class for decoding the EGM data
	*
	* The document contains the coefficiet and sigma values in different order
	*    of the gravitational field.\n
	* The gcoder t_egm2008 corresponding to the gdata egmcoef
	*/
	class  OI_LIBRARY_EXPORT t_egm2008 : public t_gcoder
	{
	public:
		/**
		* @brief constructor.
		* @param[in]  s        setbase control
		* @param[in]  version  version of the gcoder
		* @param[in]  sz       size of the buffer
		*/
		t_egm2008(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_egm2008();

		/**
		* @brief decode header of the EGM data file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return 
			@retval >=0 consume size of header decoding
			@retval <0  finish reading
		*/
		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg);

		/**
		* @brief decode the data body of the EGM data file.
		*
		* decode data body of EGM file, all the data will be stored in the gegmcoef
		*
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return 
		    @retval >=0 consume size of data decoding
			@retval <0  finish reading
		*/
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

	protected:

	};
}
#endif // !EGM2008_H
