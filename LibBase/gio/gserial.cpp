#include <cstring>
#include <iostream>

#include "gserial.h"

using namespace std;

namespace gnut {

	t_gserial::t_gserial() : t_gio()
	{
		_port = 0;
		_baud = 0;
		_log = nullptr;
#ifdef _WIN32
		memset(&m_OverlappedRead, 0, sizeof(OVERLAPPED));
		memset(&m_OverlappedWrite, 0, sizeof(OVERLAPPED));
		m_hIDComDev = NULL;
#else 
		/// for linux
#endif
	}

	t_gserial::t_gserial(int port,int baud) : t_gio()
	{
		_port = port;
		_baud = baud;
		_log = nullptr;
#ifdef _WIN32
		memset(&m_OverlappedRead, 0, sizeof(OVERLAPPED));
		memset(&m_OverlappedWrite, 0, sizeof(OVERLAPPED));
		m_hIDComDev = NULL;
#else 
		/// for linux
#endif
	}

	t_gserial::t_gserial(string path) : t_gio()
	{
		path.erase(0, 12);  /// erase 'com'
		for (int i = 0; i < path.size(); i++)
		{
			if (path[i] == ':' || path[i] == ',' || path[i] == '/')path[i] = ' ';
		}
		stringstream ss(path);
		string  port, baud, data_bit, stop_bit, flow_contrl;
		ss >> port >> baud >> data_bit >> stop_bit >> flow_contrl;
		_port = str2int(port);
		_baud = str2int(baud);
		_log = nullptr;
		if(_port<0 || _baud<0) 
			if (_log) _log->comment(0, "gtcp", " error lcoal host !");
#ifdef _WIN32
		memset(&m_OverlappedRead, 0, sizeof(OVERLAPPED));
		memset(&m_OverlappedWrite, 0, sizeof(OVERLAPPED));
		m_hIDComDev = NULL;
#else 
		/// for linux
#endif
	}

	t_gserial::~t_gserial() 
	{
		close();
	}

	void t_gserial::run_read()
	{
		gtrace("t_gserial::run_read");

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
				if (_log) _log->comment(0, "gserial", " error - connect failed");
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


	int t_gserial::_gio_read(char* buff, int size)
	{
		gtrace("t_gserial::_gio_read");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		//_gmutex.lock();

#ifdef _WIN32
		if (!_opened || m_hIDComDev == NULL)
		{
			t_gtime::gmsleep(100); // [ms]     
			return(0);
		}
		BOOL bReadStatus;
		DWORD dwBytesRead, dwErrorFlags;
		COMSTAT ComStat;

		ClearCommError(m_hIDComDev, &dwErrorFlags, &ComStat);
		if (!ComStat.cbInQue) return(0);

		dwBytesRead = (DWORD)ComStat.cbInQue;
		if (size < (int)dwBytesRead) dwBytesRead = (DWORD)size;

		bReadStatus = ReadFile(m_hIDComDev, buff, dwBytesRead, &dwBytesRead, &m_OverlappedRead);
		if (!bReadStatus)
		{
			if (GetLastError() == ERROR_IO_PENDING)
				PurgeComm(m_hIDComDev, PURGE_TXCLEAR | PURGE_RXCLEAR | PURGE_RXABORT | PURGE_TXABORT);
			return(0);
		}

		//_gmutex.unlock();
		return((int)dwBytesRead);


#else
		///for linux
		return 0;
#endif
	}

	bool t_gserial::glog_set(t_glog* log)
	{
		_log = log;
		return true;
	}



	int t_gserial::_init_common()
	{
#ifdef _WIN32
		char szPort[15];
		char szComParams[50];
		DCB dcb;

		sprintf(szPort, "\\\\.\\COM%d", _port);
		m_hIDComDev = CreateFileA(szPort, GENERIC_READ | GENERIC_WRITE, 0, NULL,
			OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL | FILE_FLAG_OVERLAPPED, NULL);
		if (m_hIDComDev == NULL) return(0);

		memset(&m_OverlappedRead, 0, sizeof(OVERLAPPED));
		memset(&m_OverlappedWrite, 0, sizeof(OVERLAPPED));

		COMMTIMEOUTS CommTimeOuts;
		CommTimeOuts.ReadIntervalTimeout = 0xFFFFFFFF;
		CommTimeOuts.ReadTotalTimeoutMultiplier = 0;
		CommTimeOuts.ReadTotalTimeoutConstant = 0;
		CommTimeOuts.WriteTotalTimeoutMultiplier = 0;
		CommTimeOuts.WriteTotalTimeoutConstant = 5000;
		SetCommTimeouts(m_hIDComDev, &CommTimeOuts);

		sprintf(szComParams, "COM%d:%d,n,8,1", _port, _baud);

		m_OverlappedRead.hEvent = CreateEvent(NULL, TRUE, FALSE, NULL);
		m_OverlappedWrite.hEvent = CreateEvent(NULL, TRUE, FALSE, NULL);

		dcb.DCBlength = sizeof(DCB);
		GetCommState(m_hIDComDev, &dcb);
		dcb.BaudRate = _baud;
		dcb.ByteSize = 8;
		unsigned char ucSet;
		ucSet = (unsigned char)((FC_RTSCTS & FC_DTRDSR) != 0);
		ucSet = (unsigned char)((FC_RTSCTS & FC_RTSCTS) != 0);
		ucSet = (unsigned char)((FC_RTSCTS & FC_XONXOFF) != 0);
		if (!SetCommState(m_hIDComDev, &dcb) ||
			!SetupComm(m_hIDComDev, 10000, 10000) ||
			m_OverlappedRead.hEvent == NULL ||
			m_OverlappedWrite.hEvent == NULL) {
			DWORD dwError = GetLastError();
			/// debug
			cout << dwError << endl;
			if (m_OverlappedRead.hEvent != NULL) CloseHandle(m_OverlappedRead.hEvent);
			if (m_OverlappedWrite.hEvent != NULL) CloseHandle(m_OverlappedWrite.hEvent);
			CloseHandle(m_hIDComDev);
			return(0);
		}
		_opened = TRUE;
		return _opened;

#else
		///for linux 
		return 0;
#endif
	}

	bool t_gserial::close(void)
	{
#ifdef _WIN32
		if (!_opened || m_hIDComDev == NULL) return(TRUE);

		if (m_OverlappedRead.hEvent != NULL) CloseHandle(m_OverlappedRead.hEvent);
		if (m_OverlappedWrite.hEvent != NULL) CloseHandle(m_OverlappedWrite.hEvent);
		CloseHandle(m_hIDComDev);
		_opened = FALSE;
		m_hIDComDev = NULL;
		return(TRUE);
#else
		return true;
#endif
		

	}
} 