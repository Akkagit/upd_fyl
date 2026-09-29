/**
*
* @verbatim
History
	-1.0 Jiande Huang  2019-11-13  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		ambinp.h
* @brief	decode and enconde ion constraint file(grtyyyydoy.ion)
*             yyyy  ---- year
*              ddd  ---- doy
*
* @author   Jiande Huang, Wuhan University
* @version	1.0.0
* @date		2019-10-25
*
*/

#ifndef ION_H
#define ION_H

#include "gutils/ExportBase.h"
#include "gcoders/gcoder.h"
#include "gdata/gion.h"

using namespace std;
using namespace gnut;

namespace great
{
	/**
	*@brief	   Class for decode/encode  amb constraint file
	*/
	class BASE_LIBRARY_EXPORT t_ion : public t_gcoder
	{
	public:

		/**
		* @brief constructor.
		* @param[in]  s        setbase control
		* @param[in]  version  version of the gcoder
		* @param[in]  sz       size of the buffer
		*/
		t_ion(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_ion();

		/**
		* @brief decode header of ambcon file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return consume size of header decoding
		*/
		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg)override;

		/**
		* @brief decode body of ambcon file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return consume size of header decoding
		*/
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)override;

		/**
		* @brief encode header of ambcon file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return  size of header encoding
		*/
		virtual  int encode_head(char* buff, int sz, vector<string>& errmsg)override;

		/**
		* @brief encode header of ambcon file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return  size of data  body encoding
		*/
		virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)override;

	protected:
		t_gtime		_beg;
		t_gtime     _end;
		double		_intv;
		int			_numb;
		set<string> _sats;
		set<string> _recs;
	private:
	};


}

#endif // !ION_H
