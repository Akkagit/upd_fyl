/**
*
* @verbatim
	History
	 -1.0 Hongmin Zhang 2019-01-06 creat the file.
	 -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
	 -1.2 Yujie Qin     2019-04-10 Change the format of file,add encode head and data
	 -1.3 Yujie Qin     2019-07-01 Adding Doxygen Style Code Remarks
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		    icsleo.cpp
 * @brief		Decode and encode the LEO ICS data
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		    2019-03-27
 *
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>

#include "gcoders/icsleo.h"
#include "gutils/gleo.h"


using namespace std;
namespace great
{

	/** @brief default constructor. */
	t_ics_leo::t_ics_leo()
	{
	}

	/**
	* @brief default constructor.
	* @param[in]  s        setbase control
	* @param[in]  version  version of the gcoder
	* @param[in]  sz       size of the buffer
	*/
	t_ics_leo::t_ics_leo(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_ics::constructor");
		_tab = 0;
		_tab_data = 0;
		_gicsleo = new t_gics_leo();
	}

	/**
	* @brief default constructor.
	*/
	t_ics_leo::~t_ics_leo()
	{

	}

	/**
	* @brief decode header of LEO ICS file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size of header decoding
	*/
	int t_ics_leo::decode_head(char* buff, int sz, vector<string>& errmsg)
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

		int sat_num;

		string sys_type;
		string combie_type;
		string sys_list;
		string sat_name;
		string par_name;

		int beg_mjd, end_mjd = 0;
		double beg_sod, end_sod = 0.0;

		int sitsize = 0;
		int tmpsize = 0;
		try
		{
			while ((tmpsize = t_gcoder::_getline(line)) >= 0)
			{
				sitsize += tmpsize;

				if (line.substr(0, 1) == "%")
				{
					unsigned int pos = 0;
					if ((line.find("Satellite") != string::npos) && _tab == 0)
					{
						pos = line.find("=");
						line.erase(0, pos + 1);

						istringstream istr(line);
						istr >> sys_type >> sat_num;

						if (!istr.fail())
						{
							//system type
							_sys_type = sys_type;

							//the quantity of satellites
							_leo_num = sat_num;

							//in case of a newline
							if (istr >> newline)
							{
								_tab = 1;
								t_gcoder::_consume(tmpsize);
								continue;
							}
						}
					}
					else if (line.find("PRN") != string::npos && _tab == 0)
					{
						pos = line.find("=");
						line.erase(0, pos + 1);

						istringstream istr(line);
						istr >> sat_name;
						//get the prn of LEO satellites 
						while (sat_name != "++")
						{
							if (IS_SUPPORTED_LEO(sat_name))
							{
								_leo_list.insert(sat_name);
							}
							if (_leo_list.size() >= _leo_num)
							{
								break;
							}
							istr >> sat_name;
						}

						_tab = 3;
						t_gcoder::_consume(tmpsize);
						continue;
					}
					//Julian days
					else if (line.find("Time") != string::npos)
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
					else if (_tab == 2)
					{
						/*line.erase(0, 1);
						istringstream istr(line);

						while (par_name != "++")
						{
						_par_list.push_back(par_name);
						istr >> par_name;
						}*/

						_tab = 0;
					}
					else if (_tab == 3)
					{
						line.erase(0, 1);
						istringstream istr(line);

						while (sat_name != "++"&&sat_name != "")
						{
							_leo_list.insert(par_name);
							if (_leo_list.size() >= _leo_num)
							{
								break;
							}
							istr >> sat_name;
						}
						_tab = 0;
					}
				}
				else if (line.find("## End Of Header") != string::npos)
				{
					t_gcoder::_consume(tmpsize);
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::ALLICSLEO)
						{
							((t_gallics_leo*)it->second)->set_time(_beg_epo, _end_epo);
							//((t_gallics_leo*)it->second)->end = _end_epo;
							//((t_gallics_leo*)it->second)->sys_type = _sys_type;
							//((t_gallics_leo*)it->second)->leo_list = _leo_list;
							//((t_gallics_leo*)it->second)->leo_num = _leo_num;
						}

						it++;
					}
					_mutex.unlock();
					return -1;
				}


				t_gcoder::_consume(tmpsize);
			}
			//return sitsize;
			_mutex.unlock();

			// jdhuang : make sure unlock before return
			return sitsize;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_icsleo::decode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_icsleo::decode_head throw exception" << endl;
			}
			return -1;
		}

	}

	/**
	* @brief encode header of LEO ICS file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return encode_head data size of header encoding
	*/
	int t_ics_leo::encode_head(char* buff, int sz, vector<string>& errmsg)
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
					if (data_iter->second->id_type() != t_gdata::ALLICSLEO)
						continue;

					_beg_epo = ((t_gallics_leo*)data_iter->second)->get_beg_time();
					_end_epo = ((t_gallics_leo*)data_iter->second)->get_end_time();
					//_sys_type = ((t_gallics_leo*)data_iter->second)->sys_type;
					_sys_type = "LEO";
					_leo_list = ((t_gallics_leo*)data_iter->second)->get_leo_sat_list();
					_leo_num = ((t_gallics_leo*)data_iter->second)->get_leo_sat_num();
					_time_type = t_gtime::tsys2str(_beg_epo.tsys());
				}

				_ss << "## Header Of Initial Control State Of Satellite" << endl;



				_ss << "%Satellite        ="
					<< setw(4) << _sys_type
					<< setw(3) << _leo_num << endl;

				_ss << "%PRN			  =";
				int mark = 0;
				int sat_num = _leo_list.size();


				it_beg = _leo_list.begin();
				it_end = _leo_list.end();

				while (mark < sat_num)
				{
					if (_leo_list.size() < 10)
					{
						for (it = it_beg; it != it_end; it++)
						{
							_ss << setw(5) << "  "<<*it << "  ";
							mark++;
						}
						_ss << endl;
					}
					else
					{
						for (it = it_beg; it != it_end; it++)
						{
							_ss << setw(5) << *it;
							mark++;
						}

						_ss << setw(3) << "++" << endl;

						while (mark < sat_num)
						{
							_ss << setw(20) << setiosflags(ios::left) << "%";
							for (int k = 0; k < 11; k++)
							{
								if (mark >= sat_num)
								{
									break;
								}
								_ss << setw(5) << *it;
								it++;
								mark++;
							}
							if (mark >= _leo_list.size())
							{
								break;
							}
							_ss << setw(3) << "++" << endl;
						}
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


				//_ss << "%Parameter        =" << " ";

				//for (mark_par = 0; mark_par < 6; mark_par++)
				//{
				//	_ss << _par_list[mark_par] << " ";
				//}
				//_ss << "++" << endl;

				//_ss << setw(20) << setiosflags(ios::left) << "%";
				//while (mark_par < _par_list.size())
				//{
				//	_ss << _par_list[mark_par] << "  ";
				//	mark_par++;
				//}
				//_ss << endl;
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
				_log->comment(0, "ERROR : t_icsleo::encode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_icsleo::encode_head throw exception" << endl;
			}
			_mutex.unlock();
			return -1;
		}
	}

	/**
	* @brief decode data body of LEO ICS file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size for data body decoding
	*/
	int t_ics_leo::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
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
		static string sat_name, sat_type;
		int mjd_begin, mjd_end, sod_begin, sod_end;
		double sd_begin, sd_end, sec_begin, sec_end;

		//_gicsleo = new t_gicsleo;
		try
		{
			static vector<t_icsleo_record> par_all;  // for LEO PWC parameters
			while ((tmpsize = t_gcoder::_getline(line)) >= 0)
			{
				sitsize += tmpsize;
				if (line.find("END of FILE") != string::npos)
				{
					_mutex.unlock();
					return -1;
				}
				if (line.find("END of SAT") != string::npos)
				{
					if (_record.isPWC == true)
					{
						(*_gicsleo).set_par_record(par_all[0].str_name, par_all);
						par_all.clear();
						par_all.push_back(_record);
					}


					par_all.clear();
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::ALLICSLEO)
						{
							//save the information for one satellite
							((t_gallics_leo*)it->second)->add_leo_ics(sat_name, *_gicsleo);
						}
						it++;
					}
					_tab_data = 0;


					_gicsleo = new t_gics_leo();
					t_gcoder::_consume(tmpsize);
					continue;
				}

				//get the information for one satellite
				//get the first line for each satellite
				if (_tab_data == 0)
				{
					istringstream istr(line);
					istr >> sat_name >> sat_name >> sat_type;

					(*_gicsleo).sat_name = sat_name;


					//gics.setpartype(typepar);
					(*_gicsleo).sat_name = sat_name;
					(*_gicsleo).sat_type = sat_type;
					_tab_data++;
					t_gcoder::_consume(tmpsize);
					continue;
				}

				istringstream istr(line);
				string parname;
				istr >> parname >> _record.value >> _record.sig >> mjd_begin >> sd_begin >> mjd_end >> sd_end;

				sod_begin = (int)sd_begin;
				sec_begin = sd_begin - sod_begin;
				_record.beg.from_mjd(mjd_begin, sod_begin, sec_begin);

				sod_end = (int)sd_end;
				sec_end = sd_end - sod_end;
				_record.end.from_mjd(mjd_end, sod_end, sec_end);

				_record.duration = mjd_end * 86400.0 + sd_end - mjd_begin * 86400.0 - sd_begin;
				size_t ipos = parname.find(":");
				if (ipos != parname.npos)  // PWC parameter
				{
					_record.str_name = parname.substr(0, ipos);
					_record.enum_name = str2icsleo(parname.substr(0, ipos));
					_record.piecewise_time = str2int(parname.substr(ipos + 1));
					_record.isPWC = true;
					if (par_all.size() < 1)  // first one
					{
						par_all.push_back(_record);
					}
					else // not first one
					{
						if (par_all[par_all.size() - 1].str_name == _record.str_name)  // the same parameter
						{
							par_all.push_back(_record);
						}
						else  // the next different parameter
						{
							(*_gicsleo).set_par_record(par_all[0].str_name, par_all);
							par_all.clear();
							par_all.push_back(_record);
						}
					}
				}
				else
				{
					par_all.clear();
					_record.str_name = parname;
					_record.enum_name = str2icsleo(parname);
					_record.isPWC = false;
					par_all.push_back(_record);
					(*_gicsleo).set_par_record(_record.str_name, par_all);
					par_all.clear();
				}

				_tab_data++;

				t_gcoder::_consume(tmpsize);

			}

			//ofstream outFile;//创建了一个ofstream 对象
			//outFile.open("information.txt");//outFile 与一个文本文件关联

			_mutex.unlock();
			return -1;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_icsleo::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_icsleo::decode_data throw exception" << endl;
			}
			_mutex.unlock();
			return -1;
		}
	}

	/**
	* @brief encode data body of attitude file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return data size for data body encoding
	*/
	int t_ics_leo::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();

		int tab = 0;
		int size = 0;
		int sd_begin, sd_end;

		//t_gallics gallics;
		//map<ICSPAR, vector<t_ics_record>> onesat;
		//vector<t_ics_record> parrec;
		//map<string, t_gics>::iterator iter;
		//map<ICSPAR, vector<t_ics_record>>::iterator iter_onesat;
		try
		{
			_ss.seekg(0, _ss.end);
			long len = (long)_ss.tellg();
			len = len - _ss_position;
			_ss.seekg(_ss_position, ios_base::beg);
			int min;
			if (len == 0 && _tab_data < _leo_list.size())
			{
				map<string, t_gdata*>::iterator it = _data.begin();
				for (it = _data.begin(); it != _data.end(); ++it)
				{
					if (it->second->id_type() == t_gdata::ALLICSLEO)
					{
						auto sat_iter = _leo_list.begin();
						advance(sat_iter, _tab_data);
						tab = 0;
						// get onesat parlist
						
						auto onesat = ((t_gallics_leo*)it->second)->get_leo_ics(*sat_iter)->leoics();
						
						string sattype = ((t_gallics_leo*)it->second)->get_leo_ics(*sat_iter)->sattype();
						for (auto iter_satpar = onesat->begin(); iter_satpar != onesat->end(); iter_satpar++)
						{
							auto parrec = iter_satpar->second;
							
							sd_begin = parrec[0].beg.sod() + parrec[0].beg.dsec();
							sd_end = parrec[0].end.sod() + parrec[0].end.dsec();
							min = (parrec[0].end - parrec[0].beg) / 60;
							for (int rec = 0; rec < parrec.size(); rec++)
							{
								sd_begin = parrec[rec].beg.sod() + parrec[rec].beg.dsec();
								sd_end = parrec[rec].end.sod() + parrec[rec].end.dsec();

								string name = parrec[rec].str_name + ":" + int2str(min);
								if (tab == 0)
								{				
									_ss << *sat_iter << "            "<< *sat_iter <<"      LEO" << endl;
									tab++;
									if (parrec[rec].isPWC && parrec.size() != 1)
									{
										
										_ss << setw(15) << left << name;
									}
									else
									{
										_ss << setw(15) << left << parrec[rec].str_name;
									}
									_ss << setw(32) << setprecision(14) << scientific << right << parrec[rec].value
										<< " " << setw(14) << setiosflags(ios::right) << setiosflags(ios::scientific) << setprecision(6) << parrec[rec].sig
										<< "  " << setw(7) << parrec[rec].beg.mjd()
										<< " " << setw(11) << fixed << setprecision(5) << sd_begin
										<< "  " << setw(7) << parrec[rec].end.mjd()
										<< " " << setw(11) << fixed << setprecision(5) << sd_end << endl;
								}
								else
								{
									if ((*_gicsleo).is_pwc(parrec[rec]) && parrec.size() != 1)
									{
										_ss << setw(15) << left << name;
									}
									else
									{
										_ss << setw(15) << left << parrec[rec].str_name;
									}
									_ss << setw(32) << setprecision(14) << scientific << right << parrec[rec].value
										<< " " << setw(14) << setiosflags(ios::right) << setiosflags(ios::scientific) << setprecision(6) << parrec[rec].sig
										<< "  " << setw(7) << parrec[rec].beg.mjd()
										<< " " << setw(11) << fixed << setprecision(5) << sd_begin
										<< "  " << setw(7) << parrec[rec].end.mjd()
										<< " " << setw(11) << fixed << setprecision(5) << sd_end << endl;
									tab++;

								}

							}

						}
						_ss << "END of SAT" << endl;
						_tab_data++;

					}
				}
			}
			else if (len == 0 && _tab_data == _leo_list.size()) {
				_ss << "END of FILE" << endl;
				_tab_data++;
			}
			/*outFile << _ss.str();
			outFile.close();*/

			int _size = _fill_buffer(buff, sz);

			_mutex.unlock();  return _size;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_icsleo::encode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_icsleo::encode_data throw exception" << endl;
			}
			return -1;
		}
	}
}