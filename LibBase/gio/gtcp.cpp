
/* ----------------------------------------------------------------------

-*/

/* paths:
 *            tcp://host:port
-*/

#include <cstring>
#include <iostream>

#include "gtcp.h"

using namespace std;

namespace gnut {

	t_gtcp::t_gtcp() :t_gio()
	{
		_user = "";
		_pwd = "";
		_host = "";
		_port = "";
		_log = nullptr;
	}

	t_gtcp::t_gtcp(const int & port)
	{
		_user = ""; _pwd = ""; _port = int2str(port);
	}

	t_gtcp::t_gtcp(const string& address) :t_gio()
	{
		string addr = address;
		for (int i = 0; i < address.size(); i++)
		{
			if (address[i] == ':' || address[i] == '/')addr[i] = ' ';
		}
		stringstream ss(addr);
		string type;
		_user = ""; _pwd = "";
		ss >> type >> _host >> _port;
		if (_host == "localhost" && _localhost(_host) < 0)
		{
			if (_log) _log->comment(0, "gtcp", " error lcoal host !");
		}
	}

	t_gtcp::t_gtcp(const string& address, bool upload) :t_gio()
	{
		_upload = upload;
		string addr = address;
		for (int i = 0; i < address.size(); i++)
		{
			if (address[i] == ':' || address[i] == '/')addr[i] = ' ';
		}
		stringstream ss(addr);
		string type;
		_user = ""; _pwd = "";
		ss >> type >> _host >> _port;
		if (_host == "localhost" && _localhost(_host) < 0)
		{
			if (_log) _log->comment(0, "gtcp", " error lcoal host !");
		}
	}

	t_gtcp::~t_gtcp()
	{

	}

	void t_gtcp::run_read()
	{
		gtrace("t_gtcp::run_read");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();

		// be sure always clean decoder
		if (_coder) { _coder->clear(); }
		// For error message
		vector<string> errmsg;

		char* loc_buff = new char[_size];

		int nbytes = 0;
		int decoded = 0;
		_stop = 0;
		_running = 1;

		while (1)
		{
			if (init_read() < 0) {
				if (_log) _log->comment(0, "gtcp", " error - connect failed");
				t_gtime::gmsleep(100); // [ms]       // [s] linux only
				continue;
			}
			while (((nbytes = _gio_read(loc_buff, _size)) > 0) && _stop != 1)
			{
				// archive the stream
				//_locf_write(loc_buff, nbytes);

				if (_coder && nbytes > 0) {
					decoded = _coder->decode_data(loc_buff, nbytes, _count, errmsg);
				}
				else {
					if (_log) _log->comment(2, "gio", "0 data decoded");
					decoded = 0;
				}

				//if (_log && _log->verb() >= 9) _log->comment(9, "gio", "read: " + int2str(nbytes)
				//	+ " decoded: " + int2str(decoded));
			}
		}
		_stop_common();
		delete[] loc_buff;

		_gmutex.unlock();
	}

	void t_gtcp::run_write()
	{
		gtrace("t_gtcp::run_write");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();

		// be sure always clean decoder
		if (_coder) { _coder->clear(); }
		// For error message
		vector<string> errmsg;

		char* loc_buff = new char[_size];

		int nbytes = 0;
		int decoded = 0;
		_stop = 0;
		_running = 1;

		while (1)
		{
			if (_Server_Init() <= 0) {
				if (_log) _log->comment(0, "gtcp", " error - connect failed");
				_opened = 0;
				t_gtime::gmsleep(100); // [ms]       // [s] linux only
				continue;
			}
			do {

				nbytes = 0;
				//// try to read from local file source
				//nbytes = _locf_read(loc_buff, _size);

				//if (_coder && nbytes < 0) {
				//	nbytes = _coder->encode_data(loc_buff, _size, _count, errmsg);
				//}
				string out = "abcdefg\n";
				//strcpy(loc_buff, out.c_str());
				loc_buff = (char*)out.c_str();
				_size = out.size();

				t_gtime::gmsleep(100);
				_stop = 0;
				nbytes = out.size();;

			} while (_stop == 0 && (_gio_write(loc_buff, nbytes) > 0 || nbytes < -1));
		}
		_stop_common();
		delete[] loc_buff;

		_gmutex.unlock();
	}

	void t_gtcp::run_send(const string& Message)
	{
		gtrace("t_gtcp::run_send");
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();

		int nbytes = 0;
		if (!_opened)//zhshen
		{
			if (_Server_Init() <= 0) {
				if (_log) _log->comment(0, "gtcp", " error - connect failed");
				_opened = 0;
				t_gtime::gmsleep(100); // [ms]       // [s] linux only
			}
		}
		nbytes = Message.length();

		if (_gio_write(Message.c_str(), nbytes) < 0)_opened = 0;

		_gmutex.unlock();
	}

	int t_gtcp::_gio_write(const char* buff, int size)
	{
		gtrace("t_gtcp::_gio_write");

		int  sz_total = size;
		int  sz_ready = 0;
		int  sz_sent = 0;
		int  sz_buffer = 1024;
		int  wait_T = 0;
		char buffer[1024]; // sz_buffer (but ISO C++ forbids variable length array)

		// wait-time in usec to provide kBps bandwidth
		if (bandwidth != 0) wait_T = (sz_buffer * 1000000 / bandwidth);

		// open failed
		if (_opened == 0) return -1;

		while (sz_total > 0) {

			// define size to be sent
			string tmp = buff;
			size_t ifirst = 0;
			if ((ifirst = tmp.find(crlf, sz_sent)) != string::npos)
				sz_ready = tmp.substr(sz_sent, ifirst + 1 - sz_sent).length();
			else
				sz_ready = sz_total;

			//strncpy( buffer, (in_buff + sz_sent), sz_ready );
			memcpy(buffer, (buff + sz_sent), sz_ready);

			sz_total -= sz_ready;

			if (wait_T > 0) t_gtime::gusleep(wait_T);  // microseconds  (e.g. 1.2kBps -> 200000 us when sz=256)
			int sz_tmp = 0;
			string new_buffer;
			if (_using_ntrip_2_http) {
				ostringstream hex_num;
				// we need the number in hexa format instead of decimal
				hex_num << std::hex << sz_ready;
				// convert the stream to string
				string hex_num_str = string(hex_num.str());
				for (unsigned int i = 0; i < hex_num_str.size(); i++) { // capitalize a-f
					// ascii 97 (decimal) is 'a', ascii 102 (decimal) is 'f'
					if (hex_num_str[i] > 96 && hex_num_str[i] < 103) {
						// 'a' is 97 (decimal) ascii, 'A' is 65 (decimal) ascii, thus 32
						hex_num_str[i] -= 32;
					}
				}
				// add header to data being sent
				new_buffer = hex_num_str + string("\r\n") + string(buffer, sz_ready) + string("\r\n");
			}
			const int new_sz_ready = _using_ntrip_2_http ? new_buffer.size() : sz_ready;
			/* buffer slightly bigger than original (original has 256), so we have space
			* for header */
			char buffer_big[1024 + 256];
			// copy the data with header to buf for call of send()
			memcpy(buffer_big, _using_ntrip_2_http ? new_buffer.c_str() : buffer, new_sz_ready);

#if  defined _WIN32 || defined _WIN64
			if (sz_ready != (sz_tmp = send(_fd_acpt, buffer_big, new_sz_ready, 0))) {            // NOT TESTED, ONLY COMPILED
#else
			if (new_sz_ready != (sz_tmp = send(_fd_skt, buffer_big, new_sz_ready, MSG_NOSIGNAL))) {  // MSG_DONTWAIT
#endif
				if (_log) _log->comment(0, "gtcp", " warning - sent to socket failed ");
				return -1;
			}
			else {
				ostringstream lg;
				lg << _fd_acpt << " sent to socket success: " << sz_ready;   // << " " << buffer << endl;
				if (_log)  _log->comment(0, "gtcp", lg.str());
				//else        cerr << lg.str() << endl;
			}
			/* so we have the size of data without header, if ntrip 2 not used then
			* new_sz_ready == sz_ready == sz_tmp and this changes nothing */
			sz_tmp = sz_ready;
			sz_sent += sz_tmp;
		}

		//if( sz_sent == 0 ) sz_sent = -1; // to set non-zero return value!
		if (_log)  _log->comment(2, "gtcp", " written bytes: " + int2str(sz_sent));
		return sz_sent;
	}

	int t_gtcp::_gio_read(char* buff, int size)
	{
		gtrace("t_gtcp::_gio_read");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		//_gmutex.lock();

		//int size = 1024;
		char* Buff = new char[size];
		int  nBytes;

#ifdef _WIN32 || defined _WIN64
		if ((nBytes = recv(_fd_skt, buff, sizeof(Buff) - 1, 0)) >= 1) {

#else
		if ((nBytes = recv(_fd_skt, buff, sizeof(Buff) - 1, MSG_WAITALL)) >= 1) {

#endif // _WIN32

			buff[nBytes] = '\0';
			//buff = Buff;
		}
		else {
			if (_log) _log->comment(0, "gtcp", " Error any response [reading data]");
		}
		delete[] Buff;
		//_gmutex.unlock();
		return nBytes;
	}

	int t_gtcp::_localhost(string & host)
	{
		try
		{
#ifdef _WIN32 || defined _WIN64
			WORD wVersionRequested;
			WSADATA wsaData;
			int err;
			wVersionRequested = MAKEWORD(1, 1);
			err = WSAStartup(wVersionRequested, &wsaData);
			if (err != 0)
			{
				if (_log) _log->comment(0, "t_gtcp", "WSAStartup error ");
			}
			PHOSTENT hostinfo;
			char hostname[255] = { 0 };
			gethostname(hostname, sizeof(hostname));
			if ((hostinfo = gethostbyname(hostname)) == NULL)
			{
				errno = GetLastError();
				fprintf(stderr, "gethostbyname Error:%d\n", errno);
				return -1;
			}
			host = inet_ntoa(*(struct in_addr*) * hostinfo->h_addr_list);
			return 1;
#endif // _WIN32
			struct hostent* hostinfolinux;
			char hostnamelinux[255] = { 0 };
			gethostname(hostnamelinux, sizeof(hostnamelinux));
			if ((hostinfolinux = gethostbyname(hostnamelinux)) == NULL)
			{
				perror("t_gtcp");
				fprintf(stderr, "gethostbyname Error:%d\n", errno);
				return -1;
			}
			host = inet_ntoa(*(struct in_addr*) * hostinfolinux->h_addr_list);
			return 1;
		}
		catch (...)
		{
			if (_log) _log->comment(0, "t_gtcp", "get local host fail! ");
			return -1;
		}
	}

	int t_gtcp::_init_common()
	{
		struct hostent* Host;
		struct sockaddr_in addr;

#ifdef _WIN32 || defined _WIN64
		WORD wVersionRequested;
		WSADATA wsaData;
		int err;
		wVersionRequested = MAKEWORD(1, 1);
		err = WSAStartup(wVersionRequested, &wsaData);
		if (err != 0)
		{
			if (_log) _log->comment(0, "t_gtcp", "WSAStartup error ");
		}
#endif // _WIN32

		_fd_skt = -1;
		_fd_skt = socket(AF_INET, SOCK_STREAM, 0);
		if (socket(AF_INET, SOCK_STREAM, 0) == -1)
		{
			if (_log)_log->comment(0, strerror(errno));
		}
		if (_fd_skt < 0) {
			if (_log)_log->comment(0, "*** Error: opening socket ***");
			return false;
		}

		// HOSTNAME
		if (!(Host = gethostbyname(_host.c_str())))
		{
			if (_log)_log->comment(0, "*** Error: no such host :", _host + "***");
			return false;
		}
#ifdef _WIN32 || defined _WIN64
		memset((char*)&addr, 0, sizeof(addr));
#else
		bzero((char*)&addr, sizeof(addr));
#endif // _WIN32

		addr.sin_family = AF_INET;

#ifdef _WIN32 || defined _WIN64
		memcpy((char*)&addr.sin_addr.s_addr, (char*)Host->h_addr, Host->h_length);
#else
		bcopy((char*)Host->h_addr, (char*)&addr.sin_addr.s_addr, Host->h_length);
#endif // _WIN32

		addr.sin_port = htons(atoi(_port.c_str()));

		// CONNECT
		try {
			if (connect(_fd_skt, (struct sockaddr*) & addr, sizeof(addr)) < 0)
			{
				//if (_log)_log->comment(0, "***Error: connecting to:", inet_ntoa(addr.sin_addr) + '***');
			}
		}
		catch (...)
		{
			if (_log)_log->comment(0, "***Error: connecting to:", inet_ntoa(addr.sin_addr) + '***');
			return false;
		}
		return true;
	}

	string t_gtcp::user()
	{
		return _user;
	}

	string t_gtcp::pwd()
	{
		return _pwd;
	}

	string t_gtcp::host()
	{
		return _host;
	}

	string t_gtcp::port()
	{
		return _port;
	}

	bool t_gtcp::glog_set(t_glog * log)
	{
		_log = log;
		return true;
	}

	bool t_gtcp::_Server_Init()
	{
		_close_socket();

		struct hostent* Host;
		struct sockaddr_in saddr;

#ifdef _WIN32 || defined _WIN64
		WORD wVersionRequested;
		WSADATA wsaData;
		int err;
		wVersionRequested = MAKEWORD(1, 1);
		err = WSAStartup(wVersionRequested, &wsaData);
		if (err != 0)
		{
			if (_log) _log->comment(0, "t_gtcp", "WSAStartup error ");
		}
#endif // _WIN32

		_fd_skt = -1;
		_fd_skt = socket(AF_INET, SOCK_STREAM, 0);
		if (socket(AF_INET, SOCK_STREAM, 0) == -1)
		{
			if (_log)_log->comment(0, strerror(errno));
		}
		if (_fd_skt < 0) {
			if (_log)_log->comment(0, "*** Error: opening socket ***");
			return false;
		}


		saddr.sin_family = AF_INET;

		saddr.sin_port = htons(atoi(_port.c_str()));
		//saddr.sin_addr.s_addr = htons(0);
		saddr.sin_addr.s_addr = htonl(INADDR_ANY);
		
		int retVal = bind(_fd_skt, (const struct sockaddr*)&saddr, int(sizeof(sockaddr_in)));
		if (retVal == -1) {
			perror("t_gtcp");
			fprintf(stderr, "_init_common Error bind failed:%d\n", errno);
			//printf("Failed bind:%d\n", WSAGetLastError());
			return false;
		}
		if (_log)_log->comment(1, "*** bind success ***");

		if (listen(_fd_skt, 10) == -1) {
			perror("t_gtcp");
			fprintf(stderr, "listen failed:%d\n", errno);
			//printf("Listen failed:%d", WSAGetLastError());
			return false;
		}
		if (_log)_log->comment(1, "*** listenning ***");

		int len = sizeof(_caddr);
		fprintf(stderr, "waiting for client\n");
#ifdef _WIN32 || defined _WIN64
		_fd_acpt = accept(_fd_skt, (sockaddr *)&_caddr, &len);
#else
		_fd_acpt = accept(_fd_skt, (sockaddr *)&_caddr, (socklen_t*)(&len));

#endif // _WIN32

		if (_fd_acpt == -1 )
		{ 
			fprintf(stderr, "accept failed:%d\n", errno);
		}

		if (_log)_log->comment(1, "*** accepted ***");
		_opened = 1;

		return true;
	}

	bool t_gtcp::_close_socket()
	{
#ifdef _WIN32 || defined _WIN64
		return closesocket(_fd_skt);
#else
		return close(_fd_skt);
#endif // _WIN32
	}

} // namespace