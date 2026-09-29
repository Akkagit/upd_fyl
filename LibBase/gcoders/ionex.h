/**
*
* @verbatim
History
	-1.0 Guolong Feng  2020-07-13  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		ionex.h
* @brief	 decode TEC grid products in The IONosphere Map EXchange Format
			 不同机构的产品：jpl/cod/wum  组合产品igsg 事后产品/预测产品 (e.g. igsg1310.18i / codg1310.18i)

* @version	1.0.0
* @date		2020-07-13
*
*/

#ifndef ION_H
#define ION_H

#include "gutils/ExportBase.h"
#include "gcoders/gcoder.h"
#include "gdata/gionex.h"

using namespace std;
using namespace gnut;

namespace great
{
	/**
	*@brief	   Class for decode/encode TEC grid products in IONEX format
	*/
	class BASE_LIBRARY_EXPORT t_ionex : public t_gcoder
	{
	public:

		/**
		* @brief constructor.
		* @param[in]  s        setbase control
		* @param[in]  version  version of the gcoder
		* @param[in]  sz       size of the buffer
		*/
		t_ionex(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_ionex();

		/**
		* @brief decode header
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return consume size of header decoding
		*/
		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg)override;

		/**
		* @brief decode body
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return consume size of header decoding
		*/
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)override;

	protected:
		t_gionex_head    _ionex_hd;
		string           _data_type;
		t_gtime          _crt_time;
		int              _ihgt, _ilat;
		int              _ilon, _nlon;
	private:
	};


}

#endif // !ION_H
