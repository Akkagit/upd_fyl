/**
*
* @verbatim
History
	-1.0 Zheng Hongjie  2019-10-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gcoder_buffer.h
* @brief	main decode function for all files in base part
*
*
* @author   Jiande Huang, Wuhan University
* @version	1.0.0
* @date		2019-10-25
*
*/

#ifndef GCODER_BUFFER_H
#define GCODER_BUFFER_H

#include <string>
#include <deque>
#include <queue>

using namespace std;




namespace gnut {
	class gcoder_buffer 
	{
	public:
		gcoder_buffer();
		~gcoder_buffer();

		int size();
		int add(char* buff, int size);
		int getline(string& str, int from_pos);
		int consume(int bytes_to_eat);
		void toString(string& str);

	private:

		deque<char> _buffer;
	};

	class gcoder_char_buffer
	{
	public:
		gcoder_char_buffer();
		~gcoder_char_buffer();


		int size();
		int add(char* buff, int size);
		int getline(string& str, int from_pos);
		int getbuffer(const char *& buff);
		int consume(int bytes_to_eat);
		void toString(string& str);

	private:

		int _begpos;
		int _endpos;
		int _buffsz;
		char* _buffer;
	};
}

#endif