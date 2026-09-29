/**
*
* @verbatim
    History
     -1.0 Yujie Qin		2019-01-06 creat the file.
     -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
	 -1.2 Jie Li		2019-04-02 Adding abnormal throwing mechanism of functions
	 -1.3 Yujie Qin     2019-07-02 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			attitude.cpp
* @brief		Decode the LEO attitude data
*
* @author		Yujie Qin, Wuhan University
* @version		1.0.0
* @date			2019-03-27
*
*/

#include"gcoders/attitude.h"
#include "gall/gallattitude.h"

using namespace std;
namespace great {

	/** @brief default constructor. */
	t_attitude::t_attitude(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_attitude::constructor");
		_interval_time = 0;   // interval of the data

	}

	/** @brief default destructor. */
	t_attitude::~t_attitude()
	{
	}

	/**
	* @brief decode header of attitude file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size of header decoding
	*/
	int t_attitude::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

#ifdef DEBUG
		cout << " BUFFER : \n" << _buffer << "\n size = " << sz << " END OF BUFFER \n\n"; cout.flush();
#endif

		string tmp;
		int consume = 0;
		int tmpsize = 0;
		string filetype;
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp)) >= 0) 
			{
				consume += tmpsize;
				// first information
				if (tmp.substr(0, 2) == "%%")
				{
					if (tmp.substr(3, 1) == "H")
					{
						if (_log) _log->comment(2, "Begin to read attitude data.");
					}
					else
					{
						if (_log) _log->comment(2, "attitude data", "END OF HEADER");
						else       cerr << "warning: unknown header message :" << tmp << endl;
						t_gcoder::_consume(tmpsize);
						_mutex.unlock(); return -1;
					}

				}
				else if (tmp.substr(0, 2) == "% ")
				{
					istringstream istr(tmp.substr(2));
					int mjd;
					double sod;
					if (tmp.substr(2, 2) == "Sa")
					{
						string satellite;
						istr >> satellite >> _sat_name;
					}
					else if (tmp.substr(2, 2) == "St")
					{
						string start;
						istr >> start >> start >> mjd >> sod;
						_beg_time.from_mjd(mjd, sod, sod - floor(sod));
					}
					else if (tmp.substr(2, 2) == "En")
					{
						string end;
						istr >> end >> end >> mjd >> sod;
						_end_time.from_mjd(mjd, sod, sod - floor(sod));
					}
					else if (tmp.substr(2, 2) == "Ti")
					{
						string interval;
						istr >> interval >> interval >> _interval_time;
					}
					else
					{
						if (_log) _log->comment(1, "WARNING : t_attitude::decode_head unknown  head message");
					}
				}
				else {
					if (_log) _log->comment(2, "solargeomag", "END OF HEADER");
					else       cerr << "warning: unknown header message :" << tmp << endl;
					_mutex.unlock(); return -1;
				}

				t_gcoder::_consume(tmpsize);
			}

			_mutex.unlock();
			return consume;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_attitude::decode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_attitude::decode_head throw exception" << endl;
			}
			return -1;
		}
	}

	/**
	* @brief decode data body of attitude file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size for data body decoding
	*/
	int  t_attitude::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();

		if (t_gcoder::_add2buffer(buff, sz) == 0){ _mutex.unlock(); return 0; };
#ifdef DEBUG
		cout << " BUFFER : \n" << _buffer << "\n size = " << sz << " END OF BUFFER \n\n"; cout.flush();
#endif
		string tmp;
		int tmpsize = 0;
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp, 0)) >= 0)
			{

				int mjd;
				double sod;
				t_attitude_record att_data;
				istringstream istr(tmp);
				istr >> mjd >> sod;
				t_gtime T;
				T.from_mjd(mjd, sod, sod - floor(sod));

				double data1, data2, data3, data4;
				istr >> data1 >> data2 >> data3 >> data4;
				att_data.q1 = data1;
				att_data.q2 = data2;
				att_data.q3 = data3;
				att_data.q4 = data4;
				_attitude.set_attitude_record(T, att_data);
				
				if (abs(T - _end_time) <= _interval_time)
				{
					_attitude.set_time(_beg_time, _end_time);
					_attitude.set_inter(_interval_time);
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::AllATTITUDE)
						{
							((t_gallatt*)it->second)->add_attitude(_sat_name, _attitude);
						}
						it++;
					}
				}
				t_gcoder::_consume(tmpsize);
			}
			_mutex.unlock(); return tmpsize;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(2, "ERROR : t_attitude::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_attitude::decode_data throw exception" << endl;
			}
			return -1;
		}
	}
}
