/**
*
* @verbatim
History
	-1.0 Jie Li  2019-11-13  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		ambinp.h
* @brief	decode and enconde amb constraint file(ambinp_yyyyddd)
*             yyyy  ---- year
*              ddd  ---- doy
*
* @author   Jie Li, Wuhan University
* @version	1.0.0
* @date		2019-10-25
*
*/

#ifndef AMBCON_H
#define AMBCON_H

#include "gcoders/gcoder.h"
#include "gdata/gambcon.h"
#include "gutils/ExportAMB.h"
using namespace std;
using namespace gnut;
namespace great
{
	/**
	*@brief	   Class for decode/encode  amb constraint file
	*/
	class AMB_LIBRARY_EXPORT t_ambcon : public t_gcoder
	{
	public:

		/**
		* @brief constructor.
		* @param[in]  s        setbase control
		* @param[in]  version  version of the gcoder
		* @param[in]  sz       size of the buffer
		*/
		t_ambcon(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_ambcon();

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
		double		_interval;
		int			_num_fixamb;
		set<string> _sats;
		set<string> _sites;
	private:

		int  _decoder_data_UCUD(int tmpsize, string& line, bool&  valid, t_gamb_Dd& dd);
	};


}

#endif //AMBCON_H
