/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-01-05	creat the file.
	  -1.1  Jiande Huang	2019-02-27	change the way to read the file.
	  -1.2  Hongjie Zheng	2019-03-22	fix the problem of limited size of buff.
	  -1.3  Jie Li			2019-04-01	add code in doxygen style.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		icsnav.cpp
 * @brief		The base class used to decode and encode navigation satellite ics file information.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>

#include "gcoders/icsnav.h"
#include "gutils/gtime.h"
#include "gall/gallics_nav.h"
#include "gforce/gfcommon.h"


using namespace std;
namespace great
{

	/** @brief constructor. */
	t_ics_nav::t_ics_nav() :
		_tab(0), _tab_data(0)
	{
	}


	t_ics_nav::t_ics_nav(t_gsetbase* s, string version, int sz) :
		t_gcoder(s, version, sz),
		_tab(0), _tab_data(0)
	{
		gtrace("t_ics::constructor");

	}

	/** @brief destructor. */
	t_ics_nav::~t_ics_nav()
	{

	}

	int t_ics_nav::decode_head(char* buff, int sz, vector<string>& errmsg)
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
		string sys_type;
		string combie_type;
		string sys_list;
		string sat_name;
		string par_name;

		int    beg_mjd = 0, end_mjd = 0;
		double beg_sod = 0.0, end_sod = 0.0;

		int sitsize = 0;
		unsigned int pos = 0;
		int tmpsize = 0;
		try
		{
			while ((tmpsize = t_gcoder::_getline(line)) >= 0)
			{
				// Line satellite
				if (line.find("%Satellite") != string::npos)
				{
					pos = line.find("=");
					line.erase(0, pos + 1);
					istringstream istr(line);
					istr >> _sys_type >> combie_type >> sys_list >> _sat_num;

					if (_sys_type != "GNS")
                    {
                        _mutex.unlock();
                        return -1;
                    }

					//single system or multi-system
					_isMixed = (combie_type == "M") ? true : false;

					for (int i = 0; i < sys_list.length(); i++)
					{
						string sys_id = sys_list.substr(i, 1);
						if (sys_id == "G" || sys_id == "E" ||
							sys_id == "R" || sys_id == "C" ||
							sys_id == "J")
						{
							_sys_list.insert(sys_id);
						}
					}
				}

				// Lines of PRN
				if (line.find("%PRN") != string::npos)
				{
					if (line.find("++") == string::npos)
					{
						if (_log) {
							_log->comment(0, "ERROR : t_ics::decode_head line %PRN format is wrong,lack ++ !");
						}
						else {
							cout << "ERROR : t_ics::decode_head line %PRN format is wrong, lack++ !" << endl;
						}
						throw exception();
					}
					pos = line.find("=");
					line.erase(0, pos + 1);
					istringstream istr(line);
					istr >> sat_name;
					//get the prn of GNSS satellites 
					while (sat_name != "++")
					{
						if (sat_name.substr(0, 1) == "G" || sat_name.substr(0, 1) == "E" ||
							sat_name.substr(0, 1) == "R" || sat_name.substr(0, 1) == "C" ||
							sat_name.substr(0, 1) == "J")
						{
							if(_sat_list.count(sat_name) == 0) _sat_list.insert(sat_name);
						}
						istr >> sat_name;

						if (_sat_list.size() > _sat_num) break;
					}
				}

				// Line of Time
				if (line.find("Time") != string::npos)
				{
					pos = line.find("=");
					line.erase(0, pos + 1);

					istringstream istr(line);
					istr >> _time_type >> beg_mjd >> beg_sod >> end_mjd >> end_sod;
					if (!istr.fail())
					{
						_beg_epo.from_mjd(beg_mjd, beg_sod, beg_sod - (int)beg_sod);
						_end_epo.from_mjd(end_mjd, end_sod, end_sod - (int)end_sod);
					}
				}

				if (line.find("## End Of Header") != string::npos)
				{
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::ALLICS)
						{
							((t_gallics_nav*)it->second)->beg_time = _beg_epo;
							((t_gallics_nav*)it->second)->end_time = _end_epo;
							((t_gallics_nav*)it->second)->sys_type = _sys_type;
							((t_gallics_nav*)it->second)->sys_list = _sys_list;
							((t_gallics_nav*)it->second)->sat_list = _sat_list;
							((t_gallics_nav*)it->second)->sat_num = _sat_num;
						}
						it++;
					}
                    _mutex.unlock();
					return -1;
				}

				t_gcoder::_consume(tmpsize);
				sitsize += tmpsize;
			}
			_mutex.unlock();
			return sitsize;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_icsnav::decode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_icsnav::decode_head throw exception" << endl;
			}
			return -1;
		}
	}


	int t_ics_nav::encode_head(char* buff, int sz, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		// init ss position
		_mutex.lock();
		int mark_par = 0;
		double sd_begin, sd_end;

		set<string>::iterator it_beg;
		set<string>::iterator it_end;
		set<string>::iterator it;
		try
		{
			if (_ss_position == 0 && _tab == 0)
			{

				// initial val from data
				auto data_iter = _data.begin();
				for (data_iter; data_iter != _data.end(); data_iter++)
				{
					if (data_iter->second->id_type() != t_gdata::ALLICS)
						continue;

					_beg_epo = ((t_gallics_nav*)data_iter->second)->beg_time;
					_end_epo = ((t_gallics_nav*)data_iter->second)->end_time;
					_sys_type = ((t_gallics_nav*)data_iter->second)->sys_type;
					_sys_list = ((t_gallics_nav*)data_iter->second)->sys_list;
					_sat_list = ((t_gallics_nav*)data_iter->second)->sat_list;
					_sat_num = ((t_gallics_nav*)data_iter->second)->sat_num;
					_time_type = t_gtime::tsys2str(_beg_epo.tsys());
				}

				_ss << "## Header Of Initial Control State Of Satellite" << endl;

				string comb_type;
				comb_type = (_sys_list.size() == 1) ? "L" : "M";

				_ss << "%Satellite        ="
					<< setw(4) << _sys_type
					<< setw(3) << comb_type << "  ";

				it_beg = _sys_list.begin();
				it_end = _sys_list.end();

				string sys = "";
				for (it = it_beg; it != it_end; it++)
				{
					//lijie _sys_list=(GPS,BDS..)
					string tmp_sys = *it;
					if (tmp_sys.length() == 3)
					{
						tmp_sys = t_gsys::str2char(tmp_sys);
					}
					sys = sys + tmp_sys;
				}

				_ss << sys << "  ";
				if (_sat_num != _sat_list.size())
				{
					_ss << setw(4) << _sat_list.size() << endl;
				}
				else
				{
					_ss << setw(4) << _sat_num << endl;
				}

				// satellite prn
				int mark = 0;
				int sat_num = _sat_list.size();


				it_beg = _sat_list.begin();
				it_end = _sat_list.end();

				while (mark < sat_num)
				{
					if (_sat_list.size() < 10)
					{
						_ss << "%PRN			  =";
						for (it = it_beg; it != it_end; it++)
						{
							_ss << setw(4) << *it;
							mark++;
						}
						_ss << setw(3) << "++" << endl;
					}
					else
					{
						for (it = it_beg; it != it_end; it++)
						{
							if(mark % 10 == 0) _ss << "%PRN			  =";
							_ss << setw(4) << *it;
							mark++;
							if (mark % 10 == 0) _ss << setw(3) << "++" << endl;
						}

						if (mark % 10 != 0) _ss << setw(3) << "++" << endl;
					}
				}

				sd_begin = _beg_epo.sod() + _beg_epo.dsec();
				sd_end = _end_epo.sod() + _end_epo.dsec();
				_ss << "%Time             ="
					<< " " << _time_type
					<< "  " << _beg_epo.mjd()
					<< "  " << fixed << setprecision(3) << sd_begin
					<< "  " << _end_epo.mjd()
					<< "  " << fixed << setprecision(3) << sd_end
					<< endl;

				_ss << "## End Of Header" << endl;

				_tab++;
			}
			int size = _fill_buffer(buff, sz);

			_mutex.unlock();  return size;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(2, "ERROR : t_icsnav::encode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_icsnav::encode_head throw exception" << endl;
			}
			return -1;
		}
	}

	int t_ics_nav::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

#ifdef DEBUG
		cout << " BUFFER : \n" << _buffer << "\n size = " << sz << " END OF BUFFER \n\n"; cout.flush();
#endif


		// tmpsize for recording the processing bytes
		int tmpsize = 0;
		// linesize fore recoding the read line size
		int linesize = 0;
		// complete for marking the record completed
		bool complete = false;

		string line;
		string sat_name;
		int    mjd_begin, mjd_end, sod_begin, sod_end = 0;
		double sd_begin, sd_end, sec_begin, sec_end = 0.0;

		t_map_sats gnss_sats = GNSS_SATS();
		try
		{
			while ((linesize = t_gcoder::_getline(line)) >= 0)
			{
				// reset tempize&&complete before one record
				tmpsize = 0;
				complete = false;
				// mark for sat system and sat type(LEO or GNS)
				string sat_sys, sat_type;

				//t_gcoder::_consume(tmpsize);
				tmpsize += linesize;

				// check the end of file
				/// @todo currently, we take the fortran format into consideration
				if (line.find("End Of File") != string::npos || line.find("END of FILE") != string::npos)
				{
					t_gcoder::_consume(tmpsize);
					break;
				}

				// fine the start of one satellite
				string PRN = line.substr(0, 3);
				if (gnss_sats[gnut::GPS].count(PRN) > 0 || gnss_sats[gnut::GLO].count(PRN) > 0 ||
					gnss_sats[gnut::GAL].count(PRN) > 0 || gnss_sats[gnut::BDS].count(PRN) > 0 ||
					gnss_sats[gnut::SBS].count(PRN) > 0 || gnss_sats[gnut::QZS].count(PRN) > 0)
				{
					// read the sat infomation
					t_gics_nav _gics;

					// first line
					istringstream istr(line);
					istr >> _gics.sat_name >> sat_sys >> sat_type;


					// read the ics data
					while ((linesize = t_gcoder::_getline(line, tmpsize)) >= 0)
					{
						//t_gcoder::_consume(tmpsize);
						tmpsize += linesize;

						// gnut && fortran
						/// @todo currently, we take the fortran format into consideration
						if (line.find("End Of Satellite") != string::npos || (line.find("END of SAT") != string::npos))
						{
							complete = true;
							break;
						}

						// check the end of file
						/// @todo currently, we take the fortran format into consideration
						if (line.find("End Of File") != string::npos || line.find("END of FILE") != string::npos) break;

						istringstream istr(line);
						istr >> _ics_record.str_name >> _ics_record.value >> _ics_record.sigma >> mjd_begin >> sd_begin >> mjd_end >> sd_end;


						sod_begin = (int)sd_begin;
						sec_begin = sd_begin - sod_begin;
						_ics_record.beg_time.from_mjd(mjd_begin, sod_begin, sec_begin);

						sod_end = (int)sd_end;
						sec_end = sd_end - sod_end;
						_ics_record.end_time.from_mjd(mjd_end, sod_end, sec_end);

						//gics.setpartype(typepar);
						// check time correction
						if (abs(_ics_record.beg_time - _beg_epo) > 86400 || abs(_ics_record.end_time - _end_epo) > 86400)
						{
							if (_log) {
								_log->comment(0, "ERROR : t_ics::decode_data", PRN + " " + _ics_record.str_name + " TIME is wrong!");
							}
							else {
								cout << "t_ics::decode_data" + PRN + " " + _ics_record.str_name + " TIME is wrong!" << endl;
							}
							throw exception();
						}

						_gics.set_par_record(_ics_record.str_name, _ics_record);
					}


					// fill ics data
					if (complete)
					{
						//save the information for one satellite
						map<string, t_gdata*>::iterator it = _data.begin();
						while (it != _data.end())
						{
							if (it->second->id_type() == t_gdata::ALLICS)
							{
								((t_gallics_nav*)it->second)->add_nav_ics(_gics.sat_name, _gics);
							}
							it++;
						}

						// consume
						t_gcoder::_consume(tmpsize);
					}
					// complete is false
					else
					{
						_mutex.unlock();
						return 0;
					}
				}
				// sat can't recognized skip!!
				else
				{
					t_gcoder::_consume(tmpsize);
				}
			}

			_mutex.unlock();
			return 0;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_icsnav::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_icsnav::decode_data throw exception" << endl;
			}
			throw exception();
			return -1;
		}
	}


	int t_ics_nav::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();

		int tab = 0;
		int size = 0;
		double sd_begin, sd_end;

		//t_gallics gallics;
		//map<ICSPAR, vector<t_ics_record>> onesat;
		//vector<t_ics_record> parrec;
		//map<string, t_gics>::iterator iter;
		//map<ICSPAR, vector<t_ics_record>>::iterator iter_onesat;

		/*_ss.seekg(0, _ss.end);
		long len = (long)_ss.tellg();
		len = len - _ss_position;
		_ss.seekg(_ss_position, ios_base::beg);*/
		try
		{
			if (_ss_position == 0)
			{
				map<string, t_gdata*>::iterator it = _data.begin();
				for (it = _data.begin(); it != _data.end(); ++it)
				{
					if (it->second->id_type() == t_gdata::ALLICS)
					{
						auto sat_iter = _sat_list.begin();
						for (sat_iter; sat_iter != _sat_list.end(); sat_iter++)
						{
							GSYS gsys = t_gsys::str2gsys((*sat_iter).substr(0, 1));
							string sys = t_gsys::gsys2str(gsys);
							_ss << setw(20) << left << *sat_iter << sys << "                GNS" << endl;
							auto onesat = ((t_gallics_nav*)it->second)->get_nav_ics(*sat_iter)->get_par_list();
							for (auto iter_satpar = onesat.begin(); iter_satpar != onesat.end(); iter_satpar++)
							{
								auto parrec = iter_satpar->second;
								sd_begin = parrec.beg_time.sod() + parrec.beg_time.dsec();
								sd_end = parrec.end_time.sod() + parrec.end_time.dsec();

								_ss << setw(9) << left << parrec.str_name
									<< setw(37) << setprecision(14) << scientific << right << parrec.value
									<< " " << setw(14) << setiosflags(ios::right) << setiosflags(ios::scientific) << setprecision(6) << parrec.sigma
									<< "  " << setw(7) << parrec.beg_time.mjd()
									<< " " << setw(11) << fixed << setprecision(5) << sd_begin
									<< "  " << setw(7) << parrec.end_time.mjd()
									<< " " << setw(11) << fixed << setprecision(5) << sd_end << endl;
							}
							_ss << "END of SAT" << endl;
						}
					}
				}
				_ss << "END of FILE" << endl;
			}

			int _size = _fill_buffer(buff, sz);

			_mutex.unlock();  return _size;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_icsnav::encode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_icsnav::encode_data throw exception" << endl;
			}
			return -1;
		}

	}


}
