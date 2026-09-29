
#ifndef GTCP_H
#define GTCP_H

/* ----------------------------------------------------------------------

* @file gtcp.h
* @brief connect the TCP Server, receive the BNC Data from port
* @verbatim
	History
	-1.0	ShenZhiheng	      2019-10-20 creat the file.
	-1.1	ShenZhiheng	      2019-11-28 add the function for getting the lcoal host.
  @endverbatim

* @version 1.0
* @author ShenZhiheng
* @date 2019-10-20
* @license 
-----------------------------------------------------------------------*/

#ifdef _WIN32
#pragma comment(lib,"WS2_32")

#endif // _WIN32

#include <stdio.h>
#include <fstream>
#include <string>
#include <sys/types.h>

#ifdef _WIN32
#include <Windows.h>
#include <winsock.h>
#else
#include <sys/socket.h>
#include <netinet/in.h>
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
	class BASE_LIBRARY_EXPORT t_gtcp :public t_gio
	{
	public:
		t_gtcp();
		t_gtcp(const int& port);
		t_gtcp(const string& address);
		t_gtcp(const string& address, bool upload); //add by xiongyun
		~t_gtcp();

		virtual void run_read();
		virtual void run_write();

		// add by xiongyun/zhshen
		void run_send(const string& Message = "");

	protected:
	
		virtual int _gio_write(const char* buff, int size);
		virtual int _gio_read(char* buff, int size);

		virtual int _localhost(string& host);

		virtual int _init_common();


		string user();
		string pwd();
		string host();
		string port();

		virtual bool glog_set(t_glog* log);


	protected:
		virtual  bool _Server_Init();
		virtual  bool _close_socket();

	protected:
		string _user;
		string _pwd;
		string _host;
		string _port;

#ifdef _WIN32  ||  defined _WIN64
		SOCKET _fd_skt;
		SOCKET _fd_acpt;
#else
		int _fd_skt;
		int _fd_acpt;
#endif
		///sockaddr_in saddr;
		sockaddr_in _caddr;

		int bandwidth;  // the bandwidth of trans
		bool _using_ntrip_2_http = false; //add by xiongyun
		bool _upload = false; //add by xiongyun
		
	};

} // namespace

#endif
