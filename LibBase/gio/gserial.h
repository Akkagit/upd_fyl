#ifndef GSERIAL_H
#define GSERIAL_H

/* ----------------------------------------------------------------------

* @file gserial.h
* @brief connect the serial port, receive the Data from serial port
* @verbatim
	History
	-1.0	wanghao	      2020-11-16 creat the file.
 @endverbatim

* @version 1.0
* @author wanghao
* @date 2020-11-16
* @license
-----------------------------------------------------------------------*/

#ifdef _WIN32
#pragma comment(lib,"WS2_32")
#endif // _WIN32

#define FC_DTRDSR       0x01
#define FC_RTSCTS       0x02
#define FC_XONXOFF      0x04
#define ASCII_BEL       0x07
#define ASCII_BS        0x08
#define ASCII_LF        0x0A
#define ASCII_CR        0x0D
#define ASCII_XON       0x11
#define ASCII_XOFF      0x13

#include <stdio.h>
#include <fstream>
#include <string>
#include <sys/types.h>

#ifdef _WIN32
#include <Windows.h>
#include <winsock.h>
#else
#include <unistd.h>
#include <stdlib.h>
#include <netdb.h>
#include <arpa/inet.h>
#endif // _WIN32


#include "gio/glog.h"
#include "gutils/gmutex.h"
#include "gcoders/gcoder.h"

#define BUF_SIZE 1024

using namespace std;

namespace gnut
{
	class BASE_LIBRARY_EXPORT t_gserial :public t_gio
	{
	public:
		t_gserial();
		t_gserial(int port,int baud);
		t_gserial(string path);
		~t_gserial();
		bool close(void);
		virtual void run_read();

	protected:
		virtual int _init_common();
		virtual bool glog_set(t_glog* log);

		virtual int _gio_write(const char* buff, int size) { return 1; }
		virtual int _gio_read(char* buff, int size);
	protected:
		int _port;
		int _baud;

#ifdef _WIN32
		HANDLE m_hIDComDev;
		OVERLAPPED m_OverlappedRead, m_OverlappedWrite;
#else 

#endif // _WIN32

	};

} // namespace

#endif
