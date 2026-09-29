/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-10-24	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		lra.cpp
 * @brief		The base class used to decode lra file information.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-10-24
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>


#include "gcoders/lra.h"
#include "gall/galllra.h"


using namespace std;
namespace great
{
	t_lra::t_lra()
	{
	}

	t_lra::t_lra(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_lra::constructor");
		
	}

	t_lra::~t_lra()
	{

	}

	int t_lra::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		if (t_gcoder::_add2buffer(buff, sz) == 0) { return 0; };
		return -1;
	}

	int t_lra::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
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
		vector<string> sats;
		t_galllra::t_map_Z mapZ;
		t_galllra::t_map_A mapA;
		lra_leo lraone;

		try
		{
			while ((tmpsize = t_gcoder::_getline(line, recsize)) >= 0){
				recsize = tmpsize; //initialize

				if (line.find("START OF LRA", 60) != string::npos)
				{

					if ((tmpsize = t_gcoder::_getline(line, recsize)) <= 0){ _mutex.unlock(); return consume; }

					mapZ.clear();
					mapA.clear();
					sats.clear();
					lraone.clearall();

					while ((tmpsize = t_gcoder::_getline(line, recsize)) >= 0){
						recsize += tmpsize; // increase counter

						// END OF ANTENNA
						// end of antenna record --> process it
						if (line.find("END OF LRA", 60) != string::npos){

							// fill pcv
							map<string, t_gdata*>::iterator it = _data.begin();
							while (it != _data.end()){
								
								if (it->second->id_type() == t_gdata::ALLLRA){
									((t_galllra*)it->second)->addlra(sats, mapA);
									((t_galllra*)it->second)->azi1(sats,lraone._azi1);
									((t_galllra*)it->second)->azi2(sats,lraone._azi2);
									((t_galllra*)it->second)->dazi(sats, lraone._dazi);
									((t_galllra*)it->second)->zen1(sats, lraone._zen1);
									((t_galllra*)it->second)->zen2(sats, lraone._zen2);
									((t_galllra*)it->second)->dzen(sats, lraone._dzen);
									((t_galllra*)it->second)->wave(sats, lraone._wavelength);
								}

								it++;
							}
							consume += t_gcoder::_consume(recsize);
							recsize = 0;
							cnt++;

#ifdef DEBUG
							cout << "Pridavam antenu = " << pcv->anten() << endl << endl;
#endif

							break; // exit a single record loop

						}
						else if (line.find("AZI1 / AZI2 / DAZI", 60) != string::npos){
							lraone._azi1 = str2dbl(line.substr(2, 6));
							lraone._azi2 = str2dbl(line.substr(8, 6));
							lraone._dazi = str2dbl(line.substr(14, 6));
#ifdef DEBUG
							cout << "DAZI               [" << pcv->dazi() << "]\n";
#endif
							// ZEN1 / ZEN2 / DZEN
						}
						else if (line.find("ZEN1 / ZEN2 / DZEN", 60) != string::npos){
							lraone._zen1 = (str2dbl(line.substr(2, 6)));
							lraone._zen2 = (str2dbl(line.substr(8, 6)));
							lraone._dzen = (str2dbl(line.substr(14, 6)));
#ifdef DEBUG
							cout << "ZEN1 / ZEN2 / DZEN [" << pcv->zen1() << "] [" << pcv->zen2() << "] [" << pcv->dzen() << "]\n";
#endif

							// # OF FREQUENCIES
						}
						else if (line.find("WAVELENGTH", 60) != string::npos){
							lraone._wavelength = (str2dbl(line.substr(5, 5)));
						}
						else if (line.find("LEO", 60) != string::npos) {
							istringstream istr(line.substr(0, 60));
							string sat;
							while (istr >> sat)
							{
								sats.push_back(sat);
							}
							//continue;
						}
						else if (line.find("DESIGNER", 60) != string::npos) {
							continue;
						}
						else {
							istr.clear(); istr.str(line);
							double azi;
							istr >> azi;
							t_galllra::t_map_Z mapZ_A;
							for (double i = lraone._zen1; i <= lraone._zen2; i += lraone._dzen){
								row = (int)((i - lraone._zen1) / lraone._dzen);
								istr >> mapZ_A[row];
							}
							colum = (int)((azi - lraone._azi1) / lraone._dazi);
							mapA[colum] = mapZ_A;

							if (istr.fail()){ _mutex.unlock(); return consume; }
#ifdef DEBUG
							cout << "reading AZI-dependent fields\n";
#endif
						}

					}
				}
			}
			_mutex.unlock(); return consume;
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
			return -1;
		}

	}

}