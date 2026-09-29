/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-12-24	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		preproc.cpp
 * @brief		The base class used to encode preprocessing information.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-12-24
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>


#include "gcoders/preproc.h"
#include "gall/gallpreproc.h"
#include "gdata/gpreprocslr.h"


using namespace std;
namespace great
{
	t_preproc::t_preproc()
	{
	}

	t_preproc::t_preproc(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		preproc = new t_gpreprocslr();
		gtrace("t_preproc::constructor");

	}

	t_preproc::~t_preproc()
	{

	}

	int t_preproc::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		if (t_gcoder::_add2buffer(buff, sz) == 0) { return 0; };
		return -1;
	}

	int t_preproc::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
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
		istringstream istr;
		int tmpsize = 0;     // individual reading counter
		int consume = 0;     // total read counter
		int recsize = 0;     // single record counter
		int row, colum;
		t_gtime t;
		int doy;
		int code;
		string sta, sat;
		double sod;
		double two_way;
		double atm_humidity, atm_pressure, wave, dry_tem;
		//t_gallpreobs allOBS;
		map<int, double> temp;

		try
		{
			while ((tmpsize = t_gcoder::_getline(line, recsize)) >= 0) {
				recsize += tmpsize; //initialize

				
				if (line.find("END") != string::npos)
				{
					
					consume += t_gcoder::_consume(recsize);
					recsize = 0;
					break;
				}

				if (line.substr(0, 1) == " ")
				{
					istringstream istr(line);
					istr >> doy >> sod >> sta >>sat>>code >> two_way >> atm_humidity >> atm_pressure >> dry_tem >> wave;
					t.from_mjd(doy, (int)sod, sod - (int)sod);
					preproc->setAtmHumidity(atm_humidity);
					preproc->setAtmPressure(atm_pressure);
					preproc->setCode(code);
					preproc->setSta(sta);
					preproc->setDryTem(dry_tem);
					preproc->setTwowayTime(two_way);
					preproc->setWave(wave);

					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end()) {

						if (it->second->id_type() == t_gdata::ALLPREPROC) {
							((t_gallpreproc*)it->second)->setSLROBS(sat,t, *preproc);

						}

						it++;
					}
					//allpreproc->setSLROBS();
				}
				else
				{
					continue;
				}

				
			}

			_mutex.unlock(); return consume;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_preOBS::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_preOBS::decode_data throw exception" << endl;
			}
			return -1;
		}

	}

	int t_preproc::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
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
			int code;
			double trans_t;
			if (_ss_position == 0)
			{
				for (auto itt = _data.begin(); itt != _data.end(); itt++)
				{

					if (itt->second->id_type() == t_gdata::ALLPREPROC)
					{
						map<string,map<t_gtime, t_gpreprocslr>>* mapAll = ((t_gallpreproc*)itt->second)->getallSLROBS();
						//map<int, double>* temp;
						_ss << "_TIME_" << "                         " <<"_STATION_" << "     " << "_SAT_" <<"        "<< "_CODE_" << "         " << "_TWO_WAY_TIME_" << "   " << "_ATM_HUMIDITY_" << " " << "_ATM_PRESSURE_" << " " << "_DRY_TEMPERATURE_" << " " << "_WAVE_LENGTH_" << endl;
						for (auto itersat = mapAll->begin(); itersat != mapAll->end(); itersat++)
						{
							for (auto iter = itersat->second.begin(); iter != itersat->second.end(); iter++)
							{

								_ss << " " << setw(5) << right << iter->first.mjd() << " " << setw(20) << right << fixed << setprecision(12) << iter->first.sod() + iter->first.dsec() << " " << setw(10) << right << iter->second._sta_name << setw(13) << right << itersat->first<< setw(12) << right << iter->second._code << "    " << setw(20) << right << fixed << setprecision(12) << iter->second._two_way_time << "    "
									<< setw(10) << right << fixed << setprecision(4) << iter->second._atm_humidity << "   " << setw(11) << right << fixed << setprecision(4) << iter->second._atm_pressure << "    " << setw(12) << right << fixed << setprecision(4) << iter->second._dry_tem << "    " << setw(12) << right << fixed << setprecision(4) << iter->second._wave;

								_ss << endl;

							}
						}
						_ss << " END"<<endl;

					}
				}
			}


		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_preproc::encode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_preproc::encode_data throw exception" << endl;
			}
			return -1;
		}

		int size = _fill_buffer(buff, sz);

		_mutex.unlock();  return size;

	}

}