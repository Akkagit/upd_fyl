/**
 *
 * @verbatim
	 History
	  -1.0	Jie Li			2019-01-06	creat the file.
	  -1.1  Jie Li			2019-04-01	add code in doxygen style.
	  -1.2  Jie Li			2019-04-02  add abnormal throwing mechanism of functions
	  -1.3  Hongjie Zheng	2019-04-10  add allorb decode
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		orb.cpp
 * @brief		The base class used to decode and encode orb file information.
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 */

#include <algorithm>
#include "gcoders/orb.h"
#include "gdata/gorb.h"
#include "gall/gallorb.h"

using namespace std;

namespace great
{
	/** @brief constructor. */
	t_orb_hd::t_orb_hd()
	{

	}

	/** @brief destructor. */
	t_orb_hd::~t_orb_hd()
	{

	}
	
	/**
	* @brief constructor.
	* @param[in]  s        setbase control
	* @param[in]  version  version of the gcoder
	* @param[in]  sz       size of the buffer
	*/
	t_orb::t_orb(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_orb::constructor");
	}

	/** @brief destructor. */
	t_orb:: ~t_orb(){}

	/**
	* @brief decode header of orb file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size of header decoding
	*/
	int t_orb::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		gtrace("t_orb::decode_head");
		_mutex.lock();
		if (t_gcoder::_add2buffer(buff, sz) == 0)
		{
			_mutex.unlock();
			return 0;
		}
		string line;
		int consume = 0;
		int tmpsize = 0;
		try
		{
			while ((tmpsize = t_gcoder::_getline(line)) >= 0)
			{
				consume += tmpsize;
				int pos = 0;

				if (line.find("END OF HEAD") == string::npos)
				{
					// get satellite info
					if (line.find("Sat Group:") != string::npos)
					{
						pos = line.find(":");
						line.erase(0, pos + 1);
						stringstream ss(line);
						ss >> _orb_hd.sat_group;
					}
					else if (line.find("Sat Type:") != string::npos)
					{
						pos = line.find(":");
						line.erase(0, pos + 1);
						stringstream ss(line);
						char str_tmp;
						while (ss >> str_tmp)
						{
							_orb_hd.sat_type.insert(str_tmp);
						}
					}
					else if (line.find("Sat List:") != string::npos)
					{
						pos = line.find(":");
						line.erase(0, pos + 1);
						stringstream ss(line);
						ss >> _orb_hd.sat_num;
					}
					//get sat list
					else if (line.substr(0, 1) == " ")
					{
						stringstream ss(line);
						string sat_name;
						while (ss >> sat_name)
						{
							_orb_hd.sat_prn.insert(sat_name);
						}
					}
					// get integerator control
					else if (line.find("Intege Control:") != string::npos)
					{
						pos = line.find(":");
						line.erase(0, pos + 1);
						stringstream ss(line);
						string tmp_partial;
						ss >> tmp_partial;
						if (tmp_partial == "true")
						{
							_orb_hd.partial = true;
						}
						else
						{
							_orb_hd.partial = false;
						}
						ss >> _orb_hd.admord >> _orb_hd.adminp >> _orb_hd.admout;
					}
					else if (line.find("Nequ:") != string::npos)
					{
						pos = line.find(":");
						line.erase(0, pos + 1);
						stringstream ss(line);
						ss >> _orb_hd.nequ;
					}
					else if (line.find("Nepo:") != string::npos)
					{
						pos = line.find(":");
						line.erase(0, pos + 1);
						stringstream ss(line);
						ss >> _orb_hd.nepo;
					}
					else if (line.find("Beg Time:") != string::npos)
					{
						pos = line.find(":");
						line.erase(0, pos + 1);
						stringstream ss(line);

						int mjd = 0, sod = 0;
						double tmpsod = 0.0, ds = 0.0;
						// start time
						ss >> mjd >> tmpsod;
						sod = int(tmpsod);
						ds = tmpsod - sod;
						_orb_hd.beg.from_mjd(mjd, sod, ds);
					}
					else if (line.find("End Time:") != string::npos)
					{
						pos = line.find(":");
						line.erase(0, pos + 1);
						stringstream ss(line);

						int mjd = 0, sod = 0;
						double tmpsod = 0.0, ds = 0.0;
						// end time
						ss >> mjd >> tmpsod;
						sod = int(tmpsod);
						ds = tmpsod - sod;
						_orb_hd.end.from_mjd(mjd, sod, ds);
					}
					else if (line.find("Ref Time:") != string::npos)
					{
						pos = line.find(":");
						line.erase(0, pos + 1);
						stringstream ss(line);

						int mjd = 0, sod = 0;
						double tmpsod = 0.0, ds = 0.0;
						// ref time
						ss >> mjd >> tmpsod;
						sod = int(tmpsod);
						ds = tmpsod - sod;
						_orb_hd.t_ref.from_mjd(mjd, sod, ds);
					}
					// force model number
					else
					{
						cerr << "warning: unknown header message :" << line << endl;
					}
					t_gcoder::_consume(tmpsize);
				}
				else
				{
					//end of head

					//add head data
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end()) {
						if (it->second->id_type() == t_gdata::ORB)
						{
							((t_gorb*)it->second)->add_head(_orb_hd);
						}
						it++;
					}
					t_gcoder::_consume(tmpsize);
					_mutex.unlock();
					return -1;
				}
			}
			_mutex.unlock();
			return consume;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_orb::decode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_orb::decode_head throw exception" << endl;
			}
			return -1;
		}
	}

	/**
	* @brief decode data body of orb file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size for data body decoding
	*/
	int t_orb::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
		gtrace("t_orb::decode_data");
		_mutex.lock();
		if (t_gcoder::_add2buffer(buff, sz) == 0)
		{
			_mutex.unlock();
			return 0;
		}

		string line;
		int consume = 0;
		int tmpsize = 0;
		try
		{
			while ((tmpsize = t_gcoder::_getline(line)) >= 0)
			{
				stringstream ss(line);
				consume += tmpsize;
				if (line.substr(0, 1) == "*")
				{
					stringstream ss_temp(line.substr(1));
					ss_temp >> _epo;
				}
				//else if (line.substr(0, 1) != " "&&line.substr(0, 10) != "END")
				else if (line.substr(1, 1) != " "&&line.substr(0, 3) != "END")
				{
					ss >> _sat;
				}
				else if (line.substr(0, 3) == "   ")
				{
					t_gtriple one_funct_rec;
					string mark;
					ss >> mark >> one_funct_rec[0] >> one_funct_rec[1] >> one_funct_rec[2];
					if (!ss.fail())
					{
						map<string, t_gdata*>::iterator it = _data.begin();
						while (it != _data.end()) {
							if (it->second->id_type() == t_gdata::ORB)
							{
								((t_gorb*)it->second)->add_funct(_epo, _sat, one_funct_rec);
							}
							else if (it->second->id_type() == t_gdata::ALLORB)
							{
								((t_gallorb*)it->second)->addData(_epo, _sat, _orb_hd, one_funct_rec);
							}
							it++;
						}
					}
					else
					{
						ss.clear();
						ss.str("");
						ss << line;
						ss >> one_funct_rec[0] >> one_funct_rec[1] >> one_funct_rec[2];
						if (!ss.fail())
						{
							// add funct data 
							map<string, t_gdata*>::iterator it = _data.begin();
							while (it != _data.end()) {
								if (it->second->id_type() == t_gdata::ORB)
								{
									((t_gorb*)it->second)->add_funct(_epo, _sat, one_funct_rec);
								}
								else if (it->second->id_type() == t_gdata::ALLORB)
								{
									((t_gallorb*)it->second)->addData(_epo, _sat, _orb_hd, one_funct_rec);
								}
								it++;
							}
						}
						else
						{
							cerr << "warning: unknown orb data :" << line << endl;
						}
					}
					//t_gtriple one_funct_rec;
					//ss >> one_funct_rec[0] >> one_funct_rec[1] >> one_funct_rec[2];
					//// add funct data 
					//map<string, t_gdata*>::iterator it = _data.begin();
					//while (it != _data.end()) {
					//	if (it->second->id_type() == t_gdata::ORB)
					//	{
					//		((t_gorb*)it->second)->add_funct(_epo, _sat, one_funct_rec);
					//	}
					//	else if (it->second->id_type() == t_gdata::ALLORB)
					//	{
					//		((t_gallorb*)it->second)->add_data(_epo, _sat, _orb_hd, one_funct_rec);
					//	}
					//	it++;
					//}
				}
				else if (line.substr(0, 11) == "END OF FILE")
				{
					_mutex.unlock();
					return tmpsize;
				}
				t_gcoder::_consume(tmpsize);
			}
			_mutex.unlock();
			return consume;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_orb::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_orb::decode_data throw exception" << endl;
			}
			return -1;
		}
	}

	/**
	* @brief encode header of orb file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return  size of header encoding
	*/
	int t_orb::encode_head(char* buff, int sz, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();

		try
		{
			if (_ss_position == 0)
			{
				map<string, t_gdata*>::iterator it = _data.begin();
				for (it = _data.begin(); it != _data.end(); ++it)
				{
					if (it->second->id_type() == t_gdata::ORB)
					{
						//get orb head info
						dynamic_cast<t_gorb*>(it->second)->get_head(_orb_hd);
					}
				}
				//fill head data
				_ss << "Sat Group: " << _orb_hd.sat_group << endl;

				set<char>::iterator ithd = _orb_hd.sat_type.begin();
				_ss << "Sat Type: ";
				for (ithd; ithd != _orb_hd.sat_type.end(); ++ithd)
				{
					_ss << *ithd << " ";
				}
				_ss << endl;
				//sat list
				_ss << "Sat List: " << _orb_hd.sat_num << endl;
				set<string>::iterator itsat = _orb_hd.sat_prn.begin();
				int isat = 0;
				for (itsat; itsat != _orb_hd.sat_prn.end(); ++itsat)
				{
					_ss << "  " << *itsat;
					isat++;
					if (isat == 10)
					{
						_ss << endl;
						isat = 0;
					}
				}
				_ss << endl;
				//intege contral
				string tmp_partial;
				if (_orb_hd.partial)
				{
					tmp_partial = "true";
				}
				else
				{
					tmp_partial = "false";
				}
				_ss << "Intege Control: " << tmp_partial << "  " << _orb_hd.admord << "  " << _orb_hd.adminp << "  " << _orb_hd.admout << endl;
				//Nequ/Nepo
				_ss << "Nequ: " << _orb_hd.nequ << endl;
				_ss << "Nepo: " << _orb_hd.nepo << endl;
				//start time
				int mjd = _orb_hd.beg.mjd();
				double sod = _orb_hd.beg.sod() + _orb_hd.beg.dsec();
				_ss << "Beg Time: " << mjd << " " << sod << endl;
				//end time
				mjd = _orb_hd.end.mjd();
				sod = _orb_hd.end.sod() + _orb_hd.end.dsec();
				_ss << "End Time: " << mjd << " " << sod << endl;
				//ref time
				mjd = _orb_hd.t_ref.mjd();
				sod = _orb_hd.t_ref.sod() + _orb_hd.t_ref.dsec();
				_ss << "Ref Time: " << mjd << " " << sod << endl;

				_ss << "END OF HEAD" << endl;
			}

			int size = _fill_buffer(buff, sz);
			_mutex.unlock();
			return size;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_orb::encode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_orb::encode_head throw exception" << endl;
			}
			return -1;
		}
	}

	/**
	* @brief encode data body of orb file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return  size for data body encoding
	*/
	int t_orb::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		try
		{
			if (_ss_position == 0)
			{
				//get data from _data
				//funct data
				map<string, map<int, t_sat_funct > > funct;
				map<string, t_gdata*>::iterator it = _data.begin();
				map<string, vector<string> > parlist;
				for (it = _data.begin(); it != _data.end(); ++it)
				{
					if (it->second->id_type() == t_gdata::ORB)
					{
						dynamic_cast<t_gorb*>(it->second)->get_funct(funct);
						dynamic_cast<t_gorb*>(it->second)->get_parlist(parlist);
					}
				}

				//encode
				//fill epo funct data
				map<string, map<int, t_sat_funct > >::iterator itFunct = funct.begin();
				map<int, t_sat_funct> ::iterator itEpo = (funct.begin()->second).begin();
				//circle by epoch
				int mul_step = _orb_hd.admout / _orb_hd.adminp;
				for (itEpo; itEpo != (funct.begin()->second).end(); ++itEpo)
				{
					if (((itEpo->first) % mul_step) == 0)
					{
						_ss << "*" << (itEpo->first) / mul_step << endl;
						
						//circle by sat
						for (itFunct; itFunct != funct.end(); ++itFunct)
						{
							auto sat_par = parlist[itFunct->first];
							//circle by a record
							t_sat_funct rec = (itFunct->second)[itEpo->first];
							//write sat pos
							if (rec.size()==0)
							{
								continue;
							}
							_ss << itFunct->first << endl;
							_ss << "   " << setw(20) << left << "SAT_POS";
							for (int j = 0; j < 3; j++)
							{
								_ss << setw(25) << setprecision(14) << scientific << right << rec[0][j];
							}
							_ss << endl;
							if (sat_par.size() != (rec.size() - 1))
							{
								throw;
							}
							//write partial
							for (int i = 0; i < sat_par.size(); i++)
							{
								_ss << "   " << setw(20) << left << "dPOS/d" + sat_par[i];
								for (int j = 0; j < 3; j++)
								{
									_ss << setw(25) << setprecision(14) << scientific << right << rec[i + 1][j];
								}
								_ss << endl;
							}
							_ss << "END OF SAT" << endl;
							/*t_sat_funct::iterator itRec = rec.begin();
							for (itRec; itRec != rec.end(); ++itRec)
							{
								_ss << "   ";
								for (int i = 0; i < 3; i++)
								{
									_ss << setw(25) << setprecision(14) << scientific << right << (*itRec)[i];
								}
								_ss << endl;
							}*/
						}
						itFunct = funct.begin();
					}
				}
				_ss << "END OF FILE" << endl;
			}
			int size = _fill_buffer(buff, sz);
			_mutex.unlock();
			return size;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_orb::encode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_orb::encode_data throw exception" << endl;
			}
			return -1;
		}
	}

}