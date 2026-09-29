/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-10-24	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		slrf.cpp
 * @brief		The base class used to decode slrf2008 file information.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-10-31
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>


#include "gcoders/slrf.h"
#include "gdata/gslrf.h"
#include "gutils/gtime.h"


using namespace std;
namespace great
{
	t_slrf::t_slrf()
	{
	}

	t_slrf::t_slrf(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_slrf::constructor");

	}

	t_slrf::~t_slrf()
	{

	}

	int t_slrf::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		if (t_gcoder::_add2buffer(buff, sz) == 0) { return 0; };
		return -1;
	}

	int t_slrf::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
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
		int inum;
		int ii;
		string tmp;
		double xv;
		int year, doy, sec;
		int month=0, day=0;
		t_gtime now;

		try
		{
			while ((tmpsize = t_gcoder::_getline(line, recsize)) >= 0) {
				recsize += tmpsize; //initialize

				if (line.find("+SITE/ID") != string::npos || _tag == 1)
				{
					if (_tag == 0)
					{
						_tag = 1;
						continue;
					}

					if (line.find("-SITE/ID") != string::npos)
					{
						_tag = 0;
						continue;
					}

					if (line.substr(0, 1).compare(" ") != 0)
					{
						continue;
					}

					_nsta = _nsta + 1;

					istringstream istr(line);
					istr >> _code;
					istr >> _domes;
					istr >> _domes;
					gslrf.setCode(_nsta - 1,_code);
					gslrf.setDomes(_nsta - 1, _domes);
					gslrf.setNameGPS(_nsta - 1, line.substr(21, 9));


					if (_domes.substr(0, 1) == "-")
					{
						_nsta = _nsta - 1;
						continue;
					}

					for (int i = 0; i < _nsta - 1; i++)
					{
						if (gslrf.getOnegslrf(i)._code == gslrf.getOnegslrf(_nsta-1)._code)
						{
							_nsta = _nsta - 1;
							break;
						}
					}
					continue;
				}
				else if (line.find("+SOLUTION/ESTIMATE") != string::npos || _tag == 2)
				{
					if (_tag == 0)
					{
						_tag = 2;
						continue;
					}

					if (line.find("-SOLUTION/ESTIMATE") != string::npos)
					{
						_tag = 3;
						continue;
					}

					if (line.substr(0, 1) != " " || (line.substr(7, 3) != "STA"&&line.substr(7, 3) != "VEL"))
					{
						continue;
					}

					inum = atoi(line.substr(14, 4).c_str());

					for ( ii = 0; ii < _nsta; ii++)
					{
						if (inum == gslrf.getOnegslrf(ii)._code)
						{
							break;
						}
					}

					if (ii >= _nsta)
					{
						continue;
					}

					if (line.substr(10, 1).compare("X") == 0)
					{
						inum = 0;
					}
					else if (line.substr(10, 1).compare("Y") == 0)
					{
						inum = 1;
					}
					else if (line.substr(10, 1).compare("Z") == 0)
					{
						inum = 2;
					}

					if (line.substr(7, 3) == "VEL")
					{
						inum = inum + 3;
					}

					istringstream istr(line);

					for (int i = 0; i < 8; i++)
					{
						istr >> tmp;
					}

					istr >> xv;

					gslrf.setXsta(ii, inum, xv);

					year = atoi(line.substr(27, 2).c_str());
					doy = atoi(line.substr(30, 3).c_str());
					sec = atoi(line.substr(34, 5).c_str());

					if (year <= 20)
					{
						year = year + 2000;
					}
					else
					{
						year = year + 1900;
					}

					now.yeardoy2monthday(year, doy, &month, &day);

					now.from_ymd(year, month, day, 0, 0);
					now.from_mjd(now.mjd(), sec, 0);
					gslrf.setREF(ii, now);
					continue;
				}

				if (_tag == 3)
				{
					int nnum = 0;
					for (int j = 0; j < _nsta; j++)
					{
						if ((gslrf.getOnegslrf(j)._xsta[0] * gslrf.getOnegslrf(j)._xsta[0] + gslrf.getOnegslrf(j)._xsta[1] * gslrf.getOnegslrf(j)._xsta[1] + gslrf.getOnegslrf(j)._xsta[2] * gslrf.getOnegslrf(j)._xsta[2])<100.0)
						{

						}
						else
						{
							nnum = nnum + 1;
						}
					}
					_nsta = nnum;
					_tag = 0;
				}
				
				if (line.find("+SITE/ECCENTRICITY") != string::npos || _tag == 4)
				{
					if (_tag == 0)
					{
						_tag = 4;
						continue;
					}

					if (line.find("-SITE/ECCENTRICITY") != string::npos)
					{
						_tag = 0;
						continue;
					}
					if (line.substr(0, 1) != " ")
					{
						continue;
					}

					inum = atoi(line.substr(1, 4).c_str());

					for (ii = 0; ii < _nsta; ii++)
					{
						if (inum == gslrf.getOnegslrf(ii)._code)
						{
							break;
						}
					}

					if (ii >= _nsta)
					{
						continue;
					}

					gslrf.getOnegslrf(ii)._offset[0] = atof(line.substr(45, 9).c_str());
					gslrf.getOnegslrf(ii)._offset[1] = atof(line.substr(54, 9).c_str());
					gslrf.getOnegslrf(ii)._offset[2] = atof(line.substr(63, 9).c_str());

					if(line.substr(43,3)=="UNE")
						_log->comment(0, "Error, SITE/ECCENTRICITY XYZ should be used !");

				}

			}

			map<string, t_gdata*>::iterator it = _data.begin();

			/*if (tot_epoch == 88)
			{
				int abc = 0;
			}*/

			while (it != _data.end()) {

				if (it->second->id_type() == t_gdata::SLRF) {
					((t_gslrf*)it->second)->setGSLRF(gslrf);
				}

				it++;
			}
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_lra::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_lra::decode_data throw exception" << endl;
			}
			_mutex.unlock();
			return -1;
		}

		t_gcoder::_consume(recsize);
		_mutex.unlock();
		return recsize;

	}

	

}