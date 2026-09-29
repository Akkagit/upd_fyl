/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2020-03-20	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		preproc.cpp
 * @brief		The base class used to decode kbr information.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2020-03-20
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>


#include "gcoders/kbr.h"
#include "gall/gallkbr.h"


using namespace std;
namespace great
{
	t_kbr::t_kbr()
	{
	}

	t_kbr::t_kbr(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_kbr::constructor");

	}

	t_kbr::~t_kbr()
	{

	}

	int t_kbr::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

#ifdef DEBUG
		cout << " BUFFER : \n" << _buffer << "\n size = " << sz << " END OF BUFFER \n\n"; cout.flush();
#endif

		string line;
		string newline;
		int pos;
		int jd_beg, jd_end;
		double sod_beg, sod_end;

		int sitsize = 0;
		int tmpsize = 0;
		try
		{
			while ((tmpsize = t_gcoder::_getline(line)) >= 0)
			{
				sitsize += tmpsize;

				if (line.find("Start and stop time") != string::npos)
				{
					pos = line.find(":");
					line.erase(0, pos + 1);

					istringstream istr(line);
					istr >> jd_beg >> sod_beg >> jd_end >> sod_end;
					_beg.from_mjd(jd_beg, (int)sod_beg, sod_beg - (int)sod_beg);
					_end.from_mjd(jd_end, (int)sod_end, sod_end - (int)sod_end);

				}
				else if (line.find("Interval") != string::npos)
				{
					pos = line.find(":");
					line.erase(0, pos + 1);

					istringstream istr(line);
					istr >> _interval;

				}
				else if (line.find("LEO SAT") != string::npos)
				{
					pos = line.find(":");
					line.erase(0, pos + 1);

					istringstream istr(line);
					istr >> _sat_1 >> _sat_2;

				}
				else if (line.find("#") != string::npos)
				{
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end()) {

						if (it->second->id_type() == t_gdata::ALLKBR) {
							((t_gallkbr*)it->second)->setHead(_beg, _end, _sat_1, _sat_2, _interval);
							
						}

						it++;
					}
					t_gcoder::_consume(tmpsize);
					_mutex.unlock();
					return -1;
				}

				t_gcoder::_consume(tmpsize);

			}
			t_gcoder::_consume(tmpsize);
			_mutex.unlock();
			return sitsize;

		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_kbr::decode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_kbr::decode_head throw exception" << endl;
			}
			return -1;
		}
	}

	int t_kbr::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

#ifdef DEBUG
		cout << " BUFFER : \n" << _buffer << "\n size = " << sz << " END OF BUFFER \n\n"; cout.flush();
#endif


		int tmpsize = 0;
		int sitsize = 0;

		string line;
		double sod;
		int jd;
		t_kbr_data kbr_one;
		t_gtime now;
		int pos;
		
		//_gicsleo = new t_gicsleo;
		try
		{
			while ((tmpsize = t_gcoder::_getline(line)) >= 0)
			{
				pos = line.find("D");
				line.replace(pos, 1, "E");
				pos = line.find("D");
				line.replace(pos, 1, "E");
				istringstream istr(line);
				
				istr >> jd >> sod >> kbr_one.bias_range >> kbr_one.range_rate >> kbr_one.iono 
					>> kbr_one.ltime_cor >> kbr_one.ltime_rate >> kbr_one.ant_cor >> kbr_one.ant_rate >> kbr_one.flag;
				now.from_mjd(jd, sod, sod - (int)sod);

				map<string, t_gdata*>::iterator it = _data.begin();
				while (it != _data.end()) {

					if (it->second->id_type() == t_gdata::ALLKBR) {
						((t_gallkbr*)it->second)->setData(now, kbr_one);

					}

					it++;
				}

				t_gcoder::_consume(tmpsize);

			}

			_mutex.unlock();
			return -1;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_kbr::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_kbr::decode_data throw exception" << endl;
			}
			_mutex.unlock();
			return -1;
		}
	}

	int t_kbr::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
		try
		{
			
			if (_ss_position == 0)
			{
				for (auto itt = _data.begin(); itt != _data.end(); itt++)
				{

					if (itt->second->id_type() == t_gdata::ALLKBR)
					{
						map<t_gtime, double> mapKBR = ((t_gallkbr*)itt->second)->getOMC();
						
						for (auto iterAll = mapKBR.begin(); iterAll != mapKBR.end(); iterAll++)
						{

							_ss << " " << setw(10) << left << iterAll->first.mjd()<< " " << setw(15) << left << iterAll->first.sod() << " " << setw(20) << right << fixed << setprecision(10) << iterAll->second;
							
							_ss << endl;

						}



					}
				}
			}


		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_kbr::encode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_kbr::encode_data throw exception" << endl;
			}
			return -1;
		}

		int size = _fill_buffer(buff, sz);

		_mutex.unlock();  return size;
	
		return 0;
	}


}