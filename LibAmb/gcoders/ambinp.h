/**
*
* @verbatim
History
	-1.0 Jie Li  2019-10-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		ambinp.h
* @brief	decode and enconde ambinp file(ambinp_yyyyddd)
*             yyyy  ---- year
*              ddd  ---- doy
*
* @author   Jie Li, Wuhan University
* @version	1.0.0
* @date		2019-10-25
*
*/

#ifndef AMBINP_H
#define AMBINP_H

#include "gcoders/gcoder.h"
#include "gdata/gambinp.h"
#include "gutils/ExportAMB.h"
using namespace std;
using namespace gnut;
namespace great
{
	/**
	*@brief	   Class for decode/encode ambinp file
	*/
	class AMB_LIBRARY_EXPORT t_ambinp : public t_gcoder
	{
	public:

		/**
		* @brief constructor.
		* @param[in]  s        setbase control
		* @param[in]  version  version of the gcoder
		* @param[in]  sz       size of the buffer
		*/
		t_ambinp(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_ambinp();

		/**
		* @brief decode header of ambinp file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return consume size of header decoding
		*/
		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg)override;

		/**
		* @brief decode body of ambinp file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return consume size of header decoding
		*/
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)override;

		/**
		* @brief encode header of ambinp file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return  size of header encoding
		*/
		virtual  int encode_head(char* buff, int sz, vector<string>& errmsg)override;

		/**
		* @brief encode header of ambinp file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return  size of data  body encoding
		*/
		virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)override;
	
	protected:
		t_gambinp_head*  _ambinp_hd;
		int				_flag;			///< =1:read ambigulity;=2:read neq
	private:
	};


}

#endif //AMBINP_H
