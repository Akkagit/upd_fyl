/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2019-04-25  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    giobigf.h
* @brief		header files of IO for big file
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-17
*
*/

#ifndef GIOBIGF_H
#define GIOBIGF_H


#include "gio/giof.h"

#define MAX_BUFFER_SIZE 204800000

namespace gnut
{

	class BASE_LIBRARY_EXPORT t_giobigf :public t_giof
	{
	public:

		t_giobigf(string mask="");
		virtual ~t_giobigf();

		int write(const char* buff, int size) override;
		int flush();
		
		const int buffer_size;

	private:
		char* _buffer;
		int _current;
	};

	class BASE_LIBRARY_EXPORT t_greadtemp
	{
	public:
		t_greadtemp(string tempfilename);
		~t_greadtemp();

		void read(char* dst,int size);
		void seekg_from_cur(int lensize);

	protected:
		static bool _get_filesize(FILE* file,int64_t* file_byte_size);

	private:
		FILE* _tmpfile;
		char* _buffer;
		int32_t _current;
		int32_t _endpos;
		int64_t _filesize;
	};


}


#endif // !GIOBIGF_H