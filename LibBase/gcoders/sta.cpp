/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2020-03-06	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		sta.cpp
 * @brief		The base class used to decode range biase information.
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2020-03-06
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>


#include "gcoders/sta.h"
#include "gcoders/gcoder.h"
#include "gall/gallsta.h"


using namespace std;
namespace great
{
	t_sta::t_sta()
	{
	}

	t_sta::t_sta(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_sta::constructor");

	}

	t_sta::~t_sta()
	{

	}

	int t_sta::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		if (t_gcoder::_add2buffer(buff, sz) == 0) { return 0; };
		return -1;
	}

	int t_sta::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
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
		int tmpsize = 0;     // individual reading counter
		int consume = 0;     // total read counter
		int recsize = 0;     // single record counter
		int row, colum;
		string sta_name;
		double dx, dy, dz;

		try
		{
			while ((tmpsize = t_gcoder::_getline(line, recsize)) >= 0) {
				recsize = tmpsize; //initialize
				istringstream istr(line);
				istr >> sta_name >> dx >> dy >> dz;
				map<string, t_gdata*>::iterator it = _data.begin();
				while (it != _data.end()) {

					if (it->second->id_type() == t_gdata::ALLSTA) {
						((t_gallsta*)it->second)->addsta(sta_name, dx, dy, dz);
						
					}
					it++;
				}
				consume += t_gcoder::_consume(recsize);
				recsize = 0;
				cnt++;
			}
			_mutex.unlock(); return consume;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_sta::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_sta::decode_data throw exception" << endl;
			}
			return -1;
		}

	}

	int t_sta::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
		try
		{
			/*_ss.seekg(0, _ss.end);
			long len = (long)_ss.tellg();
			len = len - _ss_position;
			_ss.seekg(_ss_position, ios_base::beg);*/

			if (_ss_position == 0)
			{
				for (auto itt = _data.begin(); itt != _data.end(); itt++)
				{

					if (itt->second->id_type() == t_gdata::ALLSTA)
					{
						map<string, sta_cor> mapSTA = ((t_gallsta*)itt->second)->getsta();

						for (auto iterAll = mapSTA.begin(); iterAll != mapSTA.end(); iterAll++)
						{

							_ss << " " << setw(8) << left << iterAll->first << " " << setw(20) << right << fixed << setprecision(10) << iterAll->second.dx;
							_ss << " " << setw(20) << right << fixed << setprecision(10) << iterAll->second.dy;
							_ss << " " << setw(20) << right << fixed << setprecision(10) << iterAll->second.dz;
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
				_log->comment(0, "ERROR : t_result::encode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_result::encode_data throw exception" << endl;
			}
			return -1;
		}

		int size = _fill_buffer(buff, sz);

		_mutex.unlock();  return size;
	}

}