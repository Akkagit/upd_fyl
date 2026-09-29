/**
*
* @verbatim
History
-1.0 BoWong  2019-02-22  creat the file.
-1.1 BoWong  2019-04-09  Adding Doxygen Style Code Remarks
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		ambupd.h
* @brief	decode and enconde ambupd file(XXXX_ambupd_yyyyddd)
*             XXXX  ---- site name
*             yyyy  ---- year
*              ddd  ---- doy
*
* @author   BoWong, Wuhan University
* @version	1.0.0
* @date		2019-04-09
*
*/

#ifndef AMBUPD_H
#define AMBUPD_H

#include <string> 
#include <vector> 

#include "gcoders/gcoder.h"
#include "gall/gallambupd.h"
using namespace std;
using namespace gnut;
namespace great
{
	/**
	*@brief	   Class for decode/encode ambupd file
	*/
	class AMB_LIBRARY_EXPORT t_ambupd : public t_gcoder
	{
	public:

		/**
		* @brief constructor.
		* @param[in]  s        setbase control
		* @param[in]  version  version of the gcoder
		* @param[in]  sz       size of the buffer
		*/
		t_ambupd(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_ambupd();

		/**
		* @brief decode header of ambupd file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return consume size of header decoding
		*/
		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg)override;

		/**
		* @brief decode body of ambupd file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return consume size of header decoding
		*/
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)override;

		/**
		* @brief encode header of ambupd file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return  size of header encoding
		*/
		virtual  int encode_head(char* buff, int sz, vector<string>& errmsg)override;

		/**
		* @brief encode header of ambupd file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return  size of data  body encoding
		*/
		virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)override;
	
	protected:
		
	private:
	};


}

#endif //AMBUPD_H
