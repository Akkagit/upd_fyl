/**
* The REFLMOXX and EMITMOXX file provide the CERES data.\n
* Following is an example for the header of REFLMOXX/EMITMOXX file.
*
* @verbatim
##Header of emit/refl file
@endverbatim

*
* @verbatim
History
-1.0 Qian  Zhang	2019-11-9 creat the file.
-1.1 Huang Jiande   2020-09-12 Standard code style
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			albedo.h
* @brief		    Decode the CERES data named EMITMOXX/REFLMOXX
*
* @author		    Qian Zhang, Wuhan University
* @version		    1.0.0
* @date			2019-11-9
*
*/

#ifndef ALBEDO_H
#define ALBEDO_H

#include "gutils/gtime.h"
#include "gcoders/gcoder.h"
#include "gutils/ExportOI.h"

using namespace std;
using namespace gnut;

namespace great
{
	/**
	*@brief	   Class for decoding the solar_flux and geomag_kp data
	*/
	class OI_LIBRARY_EXPORT t_albedo : public t_gcoder {

	public:
		/**
		* @brief default constructor.
		*
		* @param[in]  s        setbase control
		* @param[in]  version  version of the gcoder
		* @param[in]  sz       size of the buffer
		*/
		t_albedo(t_gsetbase* s, string version, int sz);

		/** @brief default constructor. */
		t_albedo(string file_type);

		/** @brief default destructor. */
		virtual ~t_albedo() {};

		/**
		* @brief decode the header of the CERES data file.
		*
		* The function is used for decoding the head of the CERES data file.\n
		* pay attention to the buff and the size of buff which may cause some trouble when
		using the wrong value in decoding.
		*
		* @brief decode header of the CERES data file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return
		@retval >=0 consume size of header decoding
		@retval <0  finish reading
		*/
		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg);

		/**
		* @brief decode the data body of the CERES data file.
		*
		* decode data body of CERES file, all the data read will store in the galbedo
		*
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  cnt         
		* @param[in]  errmsg      error message of the data decoding
		* @return
		@retval >=0 consume size of body decoding
		@retval <0  finish reading
		*/
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

	private:
		string _file_type;    ///< file type of the input file(emit or refl)
		int    _row_index = 0;           ///< row index of the ceres data
	};
}


#endif