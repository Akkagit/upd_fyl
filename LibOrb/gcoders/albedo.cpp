/**
*
* @verbatim
History
-1.0 Qian Zhang	   2019-11-09 creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			albedo.cpp
* @brief		Decode the CERES data
*
* @author		Qian Zhang, Wuhan University
* @version		1.0.0
* @date			2019-11-09
*
*/

#include "gcoders/albedo.h"
#include "gdata/galbedo.h"

using namespace std;
namespace great
{
	/** @brief default constructor. */
	t_albedo::t_albedo(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_dvpteph405::constructor");
	}

	/**
	* @brief default constructor.
	* @param[in] filetype  input file type(emit or refl)
	*/
	t_albedo::t_albedo(string filetype)
	{
		_file_type = filetype;
		_row_index = 0;
	}

	/**
	* @brief decode header of emit/refl file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size of header decoding
	*/
	int t_albedo::decode_head(char* buff, int sz, vector<string>& errmsg)
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
		int consume = -1;
		int tmpsize = 0;
		string filetype;
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp)) >= 0) {
				consume += tmpsize;
				if (tmp.substr(0, 1) == "#")
				{
					if (tmp.substr(2, 1) == "H")
					{
						istringstream istr(tmp);
						istr >> filetype >> filetype >> _file_type;
					}

					if (_log) _log->comment(2, "The file is a", _file_type);
				}
				t_gcoder::_consume(tmpsize);
			}
			_mutex.unlock(); return consume;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_albedo::decode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_albedo::decode_head throw exception" << endl;
			}
			return -1;
		}
	}

	/**
	* @brief decode data body of emit/refl file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size for data body decoding
	*/
	int t_albedo::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
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
		int tmpsize;
		
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp, 0)) >= 0)
			{
				string s;
				double f;
				int jj = 0;
				istringstream istr(tmp);
				_row_index++;
				while (istr >> s)
				{
					jj++;
					if (s == "NaN")
					{
						f = 0.0;
						continue;
					}						
					else
					{
						istringstream istrtmp(s);
						istrtmp >> f;
					}
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::ALLALBEDO)
						{
							((t_galbedo*)it->second)->set_data(_file_type, _row_index, jj, f);
						}
						it++;
					}
					f = 0.0;
				}
				t_gcoder::_consume(tmpsize);
			}
			_mutex.unlock(); return tmpsize;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_albedo::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_albedo::decode_data throw exception" << endl;
			}
			return -1;
		}
	}
}