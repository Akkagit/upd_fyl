/**
*
* @verbatim
History
-1.0 Yujie Qin	   2019-01-06 creat the file.
-1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
-1.2 Jie Li        2019-04-02 Aa\dd abnormal throwing mechanism of functions
-1.3 Yujie Qin     2019-07-02 Adding Doxygen Style Code Remarks
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			solargeomag.cpp
* @brief		Decode the LEO solar_flux and geomag_kp data
*
* @author		Yujie Qin, Wuhan University
* @version		1.0.0
* @date			2019-03-27
*
*/

#include"gcoders/solargeomag.h"
#include "gdata/gsolargeomag.h"

using namespace std;
namespace great
{
	/** @brief default constructor. */
	t_solargeomag::t_solargeomag(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_dvpteph405::constructor");
	}

	/**
	* @brief default constructor.
	* @param[in] filetype  input file type(solar_flux or geomag_kp)
	*/
	t_solargeomag::t_solargeomag(string filetype)
	{
		_file_type = filetype;
	}

	/**
	* @brief decode header of solar_flux/geomag_kp file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size of header decoding
	*/
	int t_solargeomag::decode_head(char* buff, int sz, vector<string>& errmsg)
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
			while ((tmpsize = t_gcoder::_getline(tmp)) >= 0) {

				consume += tmpsize;

				// first information
				if (tmp.substr(0, 1) == "#")
				{
					if (tmp.substr(2, 1) == "H")
					{
						istringstream istr(tmp);
						istr >> filetype >> filetype >> _file_type;
					}

					if (_log) _log->comment(2, "The file is a", _file_type);
				}
				else if (tmp.substr(0, 1) == "%")
				{
					istringstream istr(tmp);
					if (tmp.substr(1, 1) == "T")
					{
						string time;
						int beg, end;
						istr >> time >> beg >> end >> _interval;
						_beg_mjd = beg;
						_end_mjd = end;
					}
					else if (tmp.substr(1, 4) == "Parn")
					{
						string Parnum;
						istr >> Parnum >> _par_type_num;
					}
					else if (tmp.substr(1, 4) == "Parp")
					{
						string Parpernum;
						istr >> Parpernum >> _par_line_num;
					}
					else if (tmp.substr(1, 4) == "Paru")
					{
						string Parpernum;
						double unit;
						istr >> Parpernum;
						while (istr >> unit)
						{
							_vunit.push_back(unit);
						}
					}
				}
				else {
					if (_log) _log->comment(2, "solargeomag", "END OF HEADER");
					else       cerr << "warning: unknown header message :" << tmp << endl;
					_mutex.unlock(); return -1;
				}

				t_gcoder::_consume(tmpsize);
			}

			_mutex.unlock(); return consume;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_solargeomag::decode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_solargeomag::decode_head throw exception" << endl;
			}
			return -1;
		}
	}

	/**
	* @brief decode data body of solar_flux/geomag_kp file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size for data body decoding
	*/
	int t_solargeomag::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
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

				double t;
				double f;
				double b[2];
				//double kp;
				istringstream istr(tmp);
				istr >> t;
				/*int mjd = floor(t);
				int sec = (t - mjd) * 3600;*/
				//t_gtime test()
				/*t_gtime T;
				T.from_mjd(mjd, sec,0,false);*/
				int i = 0, j = 0;
				while (i < _par_line_num)
				{
					i++;
					j = 0;
					while (j < _par_type_num)
					{
						istr >> f;
						b[j] = f;
						j++;
					}
					map<string, t_gdata*>::iterator it = _data.begin();

					if (t > _end_mjd)
						break;

					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::ALLSOLAR)
						{
							((t_gsolargeomag*)it->second)->set_data(_file_type, t, b);
							((t_gsolargeomag*)it->second)->set_time(_file_type, _beg_mjd, _end_mjd);
							((t_gsolargeomag*)it->second)->set_intv(_file_type, _interval);
						}
						it++;
					}

					t = t + _interval;
					//T = T.from_mjd(T.mjd(false), 0, 0, false);
					memset(b, 0.0, sizeof(double));
				}

				t_gcoder::_consume(tmpsize);
			}
			_mutex.unlock(); return tmpsize;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_solargeomag::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_solargeomag::decode_data throw exception" << endl;
			}
			return -1;
		}
	}
}

