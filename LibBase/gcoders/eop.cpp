/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2020-05-15	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		eop.cpp
 * @brief		The base class used to decode eop results.
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2020-05-15
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>


#include "gcoders/eop.h"
#include "gcoders/gcoder.h"
#include "gall/galleop.h"


using namespace std;
namespace great
{
	t_eop::t_eop()
	{
	}

	t_eop::t_eop(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_eop::constructor");

	}

	t_eop::~t_eop()
	{

	}

	int t_eop::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		if (t_gcoder::_add2buffer(buff, sz) == 0) { return 0; };
		return -1;
	}


	int t_eop::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
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
		double time;
		t_gtime now;
		eop eop_one;

		try
		{
			while ((tmpsize = t_gcoder::_getline(line, recsize)) >= 0) {
				recsize = tmpsize; //initialize
				istringstream istr(line);
				if (line.find("TIME") != string::npos)
				{
					consume += t_gcoder::_consume(recsize);
					recsize = 0;
					cnt++;
					continue;
				}
				istr >> time >> eop_one.xpole >> eop_one.dxpole >> eop_one.ypole >> eop_one.dypole >> eop_one.ut1 >> eop_one.dut1;
				now.from_dmjd(time, true);
				map<string, t_gdata*>::iterator it = _data.begin();
				while (it != _data.end()) {

					if (it->second->id_type() == t_gdata::ALLEOP) {
						((t_galleop*)it->second)->addeop(now, eop_one);

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


	int t_eop::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
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

					if (itt->second->id_type() == t_gdata::ALLEOP)
					{
						map<t_gtime,eop>* mapEOP = ((t_galleop*)itt->second)->geteop();

						_ss << " " << setw(15) << left << "TIME" << setw(15) << right << "XPOLE  " << setw(15) << right << "DXPOLE  " << setw(15) << right << "YPOLE  " << setw(15) << right << "DYPOLE  " << setw(15) << right << "UT1  " << setw(15) << right << "dUT1  " << endl;
						for (auto iterAll = mapEOP->begin(); iterAll != mapEOP->end(); iterAll++)
						{

							_ss <<" "<<setw(15) << left << fixed<<setprecision(5)<<iterAll->first.dmjd() - 1;
							_ss << setw(15) << right << fixed << setprecision(10) << iterAll->second.xpole;
							_ss << setw(15) << right << fixed << setprecision(10) << iterAll->second.dxpole;
							_ss << setw(15) << right << fixed << setprecision(10) << iterAll->second.ypole;
							_ss << setw(15) << right << fixed << setprecision(10) << iterAll->second.dypole;
							_ss << setw(17) << right << fixed << setprecision(10) << iterAll->second.ut1;
							_ss << setw(15) << right << fixed << setprecision(10) << iterAll->second.dut1 << endl;


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