/**
 * @verbatim
 * Name :	DE405 (May 1997)
 * Years:	1600 through 2200
 *
 * This ephemeris has been the basis for the Astronomical Almanac of the United States Naval Observatory since 2003 
 *    and powers the Jet Propulsion Laboratory's online HORIZONS service.
 *
 * Coordinates and velocities are provided for the major planets, the Sun, and Earth's Moon.
 *
 * The accuracy of this ephemeris is good; the position of the Earth or Mars, for example, should be accurate to within 2 km. 
 *
 	History
 	 -1.0 Jie Li	 2019-01-06 creat the file.
 	 -1.1 Jie Li	 2019-03-30 adding doxygen style code remarks
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file			dvpteph405.h
 * @brief			decode the binary jpleph_de405 file data 
 *
 * @author			Jie Li, Wuhan University
 * @version			1.0.0
 * @date			2019-03-30
 *
 */

#ifndef DVPTEPH405_H
#define DVPTEPH405_H

#include "gutils/ExportBase.h"

#include "gcoders/gcoder.h"
#include "gset/gsetbase.h"

namespace great
{
	/**
	 * @brief	   Class for decoding the jpleph_de405 data
	 */
	class BASE_LIBRARY_EXPORT t_dvpteph405 : public t_gcoder
	{
	public:
		/** 
		 * @brief default constructor. 
	     * 
	     * @param[in]  s        setbase control
	     * @param[in]  version  version of the gcoder
	     * @param[in]  sz       size of the buffer
		 */
		t_dvpteph405(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** 
		 * @brief default destructor. 
		 */
		virtual ~t_dvpteph405() {};

		/** 
		 * @brief decode the header of the jpleph_de405 data file. 
	     * @param[in]  buff        buffer of the data
	     * @param[in]  bufLen      buffer size of the data
	     * @param[in]  errmsg      error message of the data decoding
	     * @return consume size of header decoding
			@retval >=0 consume size of header decoding
			@retval <0  finish reading
		 */
		virtual  int decode_head(char* buff, int bufLen, vector<string>& errmsg);

		/** 
		 * @brief decode the data body of the jpleph_de405 data file. 
		 * 
	     * @param[in]  buff        buffer of the data
	     * @param[in]  bufLen      buffer size of the data
	     * @param[in]  errmsg      error message of the data decoding
	     * @return consume size for data body decoding
			@retval >=0 consume size of body decoding
		 	@retval <0  finish reading
		 */
		virtual  int decode_data(char* buff, int bufLen, int& cnt, vector<string>& errmsg);

	protected:

	private:
		int		_flag;				///< flag for decoder
		int		_version;			///< JPL file version
		t_gtime _beg;				///< start time for JPL data
		t_gtime _end;               ///< end time for JPL data
		double  _days;				///< Sampling interval of data(unit: day)
		int		_ipt[3][13];		///< index of planets for calculating position and speed
	};
}

#endif // !DVPTEPH405_H
