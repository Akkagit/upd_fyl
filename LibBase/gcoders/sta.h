/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2020-03-06	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		sta.h
 * @brief		The base class used to decode and encode range biase of slr observations.
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2020-03-06
 */

#ifndef STA_H
#define STA_H

#include "gcoders/gcoder.h"

using namespace std;
using namespace gnut;
namespace great
{

	class BASE_LIBRARY_EXPORT t_sta : public t_gcoder
	{
	public:
		/** @brief default constructor. */
		t_sta();

		/**
		 * @brief default constructor.
		 *
		 * @param[in]  s        setbase control
		 * @param[in]  version  version of the gcoder
		 * @param[in]  sz       size of the buffer
		 */
		t_sta(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_sta();


		virtual int decode_head(char* buff, int sz, vector<string>& errmsg);

		virtual int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

		virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);
	};
}

#endif