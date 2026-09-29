/**
*
* @verbatim
History
-1.0 LEO team      2019-01-06 creat the file.
-1.2 Huang Jiande  2020-09-12 Standard code style
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    desaiscopolecoef.h
* @brief		Storage the LEO KIN data(more than one satellite)
*
* @author       GREAT, Wuhan University
* @version		1.0.0
* @date		    2019-03-27
*
*/

#ifndef KIN_H
#define KIN_H

#include "gcoders/gcoder.h"
#include "gutils/ExportOI.h"
using namespace std;
using namespace gnut;
namespace great
{
	/**
	*@brief	   Class for decoding the LEO pannel data
	*/
	class OI_LIBRARY_EXPORT t_desaiscopolecoef : public t_gcoder
	{
	public:

		/**
		 * @brief default constructor.
		 *
		 * @param[in]  s        setbase control
		 * @param[in]  version  version of the gcoder
		 * @param[in]  sz       size of the buffer
		 */
		t_desaiscopolecoef(t_gsetbase* s, string version, int sz);

		/** @brief default destructor. */
		virtual ~t_desaiscopolecoef();

		/**
		 * @brief decode the header of the LEO pannel data file.
		 *
		 * The function is used for decoding the head of kin file.\n
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
		 * @brief decode the data body of the kin data file.
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
	
		int _tab = 0;
	};
} // namespace

#endif