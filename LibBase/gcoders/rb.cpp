/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2020-03-06	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		rb.cpp
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


#include "gcoders/rb.h"
#include "gcoders/gcoder.h"
#include "gall/gallrb.h"


using namespace std;
namespace great
{
	t_rb::t_rb()
	{
	}

	t_rb::t_rb(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_rb::constructor");

	}

	t_rb::~t_rb()
	{

	}

	int t_rb::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		if (t_gcoder::_add2buffer(buff, sz) == 0) { return 0; };
		return -1;
	}

	int t_rb::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
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
		string sta_name, satname;
		double RB;
		double std;

		try
		{
			while ((tmpsize = t_gcoder::_getline(line, recsize)) >= 0) {
				recsize = tmpsize; //initialize
				istringstream istr(line);

				istr >> satname >> sta_name >> RB >> std;
				
				map<string, t_gdata*>::iterator it = _data.begin();
				while (it != _data.end()) {

					if (it->second->id_type() == t_gdata::ALLRB) {
						((t_gallrb*)it->second)->addrb(sta_name,satname,RB);
						((t_gallrb*)it->second)->addstd(sta_name,satname, std);
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
				_log->comment(0, "ERROR : t_rb::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_rb::decode_data throw exception" << endl;
			}
			return -1;
		}

	}

	int t_rb::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
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

					if (itt->second->id_type() == t_gdata::ALLRB)
					{
						map<string, map<string,double>> mapRB = ((t_gallrb*)itt->second)->getrb();
						map<string, map<string, double>>* mapSTD = ((t_gallrb*)itt->second)->getstd();

						for (auto iterAll = mapRB.begin(); iterAll != mapRB.end(); iterAll++)
						{
							
							for (auto iter = iterAll->second.begin(); iter != iterAll->second.end(); iter++)
							{
								_ss << " " << setw(15) << left << iterAll->first << " " << setw(8) << left << iter->first << " " << setw(20) << right << fixed << setprecision(10) << iter->second;
								_ss << " " << setw(20) << right << fixed << setprecision(10) << (*mapSTD)[iterAll->first][iter->first];
								_ss << endl;
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