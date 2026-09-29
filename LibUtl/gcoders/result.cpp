/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-10-24	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		result.cpp
 * @brief		The base class used to decode lra file information.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-11-2
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>


#include "gcoders/result.h"
#include "gdata/gresult.h"


using namespace std;
namespace great
{
	t_result::t_result()
	{
	}

	t_result::t_result(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_result::constructor");

	}

	t_result::~t_result()
	{

	}

	void t_result::setifLEO(string ifLEO)
	{
		_ifLEO = ifLEO;
	}


	int t_result::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		

#ifdef DEBUG
		cout << " BUFFER : \n" << _buffer << "\n size = " << sz << " END OF BUFFER \n\n"; cout.flush();
#endif
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

					if (itt->second->id_type() == t_gdata::ALLRESULT)
					{
						map<string, map<int, t_gresult_one>>* mapAll = ((t_gresult*)itt->second)->getMapSat();

						for (auto iterAll = mapAll->begin(); iterAll != mapAll->end(); iterAll++)
						{
							map<int, t_gresult_one > map = iterAll->second;

							if (_ifLEO == "LEO")
							{
								for (auto iter = map.begin(); iter != map.end(); iter++)
								{
									_ss << "OMCSLR " << setw(5) << left << iter->second._now.mjd() << " " << setw(10) << right << fixed << setprecision(2) << iter->second._now.sod() + iter->second._now.dsec();
									_ss << " " << setw(5) << right << iter->second._code << " " << setw(4) << right << iter->second._name << " " << setw(7) << right << iter->second._name_long;
									_ss << " " << setw(15) << right << fixed << setprecision(5) << iter->second._omc << " " << setw(15) << right << fixed << setprecision(5) << iter->second._elev;
									_ss << " " << setw(15) << right << fixed << setprecision(5) << iter->second._dlra;
									_ss << endl;
								}
							}
							else
							{
								for (auto iter = map.begin(); iter != map.end(); iter++)
								{
									_ss << "OMCSLR " << setw(5) << left << iter->second._now.mjd() << " " << setw(10) << right << fixed << setprecision(2) << iter->second._now.sod() + iter->second._now.dsec();
									_ss << " " << setw(5) << right << iter->second._code << " " << setw(4) << right << iter->second._name << " " << setw(7) << right << iter->second._name_long;
									_ss << " " << setw(15) << right << fixed << setprecision(5) << iter->second._omc << " " << setw(15) << right << fixed << setprecision(5) << iter->second._elev;
									//_ss << " " << setw(15) << right << fixed << setprecision(5) << iter->second._dlra;
									_ss << endl;
								}
							}
							
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