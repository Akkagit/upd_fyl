/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2020-03-20	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		preproc.h
 * @brief		The base class used to decode kbr information.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2020-03-20
 */

#ifndef KBR_H
#define KBR_H

#include "gcoders/gcoder.h"
#include "gall/gallkbr.h"
#include "gdata/gkbromc.h"
#include "gutils/ExportUtl.h"

using namespace std;
using namespace gnut;
namespace great
{

	class UTL_LIBRARY_EXPORT t_kbr : public t_gcoder
	{
	public:
		/** @brief default constructor. */
		t_kbr();

		/**
		 * @brief default constructor.
		 *
		 * @param[in]  s        setbase control
		 * @param[in]  version  version of the gcoder
		 * @param[in]  sz       size of the buffer
		 */
		t_kbr(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_kbr();


		int decode_head(char* buff, int sz, vector<string>& errmsg);

		int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);
		/**
		 * @brief decode the data body of the ICS data file.
		 *
		 * decode data body of ics file, all the data read will store in the gallics_nav
		 *
		 * @param[in]  buff        buffer of the data
		 * @param[in]  sz          buffer size of the data
		 * @param[in]  errmsg      error message of the data decoding
		 * @return
			@retval >=0 consume size of body decoding
			@retval <0  finish reading
		 */
		 //virtual int decode_head(char* buff, int sz, vector<string>& errmsg);

		virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);


	protected:
		t_gtime   _start;
		t_gtime   _end;
		string    _sat_1;
		string    _sat_2;
		double    _interval;
	};
}

#endif