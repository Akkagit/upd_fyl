/**
*
* @verbatim
History
-1.0 BoWong  2019-02-22  creat the file.
-1.1 BoWong  2019-04-09  Adding Doxygen Style Code Remarks
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		upd.cpp
* @brief	decode and enconde upd file
*
* @author   BoWong, Wuhan University
* @version	1.0.0
* @date		2019-04-09
*
*/

#include "gcoders/upd.h"

using namespace std;
namespace great {

	t_upd::t_upd(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_updlane::constructor");
		_epoch = FIRST_TIME;
		_updtype == UPDTYPE::NONE;
	}

	/** @brief destructor. */
	t_upd::~t_upd()
	{

	}

	int t_upd::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();

		if (t_gcoder::_add2buffer(buff, sz) == 0)
		{
			_mutex.unlock();
			return 0;
		};
#ifdef DEBUG
		cout << " BUFFER : \n" << _buffer << "\n size = " << sz << " END OF BUFFER \n\n"; cout.flush();
#endif
		string tmp;
		int tmpsize = 0;
		int consume = 0;
		string mode_tmp;
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp, 0)) >= 0)
			{
				consume += tmpsize;
				if (tmp.substr(0, 1) == "%")
				{
					int pos = tmp.find("_");
					tmp.erase(0, pos + 1);
					tmp = cut_crlf(tmp); // remove "\n" or "\r" at the end of line (WindowsMac/Linux)
					mode_tmp = trim(tmp);     // remove "\t" at the end of line
				    //tmp.erase(tmp.end() - 1,tmp.end());
					//mode_tmp = tmp;
					if (mode_tmp == "EWL" || mode_tmp == "ewl" || mode_tmp == "WL" || mode_tmp == "wl" ||
						mode_tmp == "NL" || mode_tmp == "nl" || mode_tmp == "EWL_epoch" || mode_tmp == "ewl_epoch" || mode_tmp == "RT_SAVE" ||
						mode_tmp == "EWL24" || mode_tmp == "ewl24" || mode_tmp == "EWL25" || mode_tmp == "ewl25")
					{
						if (mode_tmp == "EWL" || mode_tmp == "ewl")
						{
							_epoch = EWL_IDENTIFY;
							_updtype = UPDTYPE::EWL;
						}
						else if (mode_tmp == "EWL24" || mode_tmp == "ewl24")
						{
							_epoch = EWL24_IDENTIFY;
							_updtype = UPDTYPE::EWL24;
						}
						else if (mode_tmp == "EWL25" || mode_tmp == "ewl25")
						{
							_epoch = EWL25_IDENTIFY;
							_updtype = UPDTYPE::EWL25;
						}
						else if (mode_tmp == "WL" || mode_tmp == "wl")
						{
							_epoch = WL_IDENTIFY;
							_updtype = UPDTYPE::WL;
						}
						else if (mode_tmp == "NL" || mode_tmp == "nl")
						{
							_updtype = UPDTYPE::NL;
						}
						else if (mode_tmp == "EWL_epoch" || mode_tmp == "ewl_epoch")
						{
							_updtype = UPDTYPE::EWL_epoch;
						}
						else if (mode_tmp == "RT_SAVE")
						{
							_updtype = UPDTYPE::RT_SAVE;
						}
						t_gcoder::_consume(tmpsize);
						break;
					}
					else
					{
						if (_log)
						{
							_log->comment(2, "t_upd::decode_head", "ERROR: unknown upd mode");
						}
						else
						{
							cout << "ERROR: t_upd::decode_head , unknown upd mode"<<endl;
						}
						t_gcoder::_consume(tmpsize);
						_mutex.unlock();
						return -1;
					}
				}
				else
				{
					_mutex.unlock();
					return -1;
				}
			}
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(2, "t_upd::decode_head", "ERROR: unknown mistake");
			}
			else
			{
				cout << "ERROR: t_upd::decode_heada , unknown mistake" << endl;
			}
			return -1;
			throw(-1);
		}
		_mutex.unlock();
		return consume;
	}

	int t_upd::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();

		if (t_gcoder::_add2buffer(buff, sz) == 0)
		{
			_mutex.unlock();
			return 0;
		};
#ifdef DEBUG
		cout << " BUFFER : \n" << _buffer << "\n size = " << sz << " END OF BUFFER \n\n"; cout.flush();
#endif
		string tmp;
		int tmpsize = 0;
		int consume = 0;
		string str;// no use

		int mjd;
		double sod;

		string prn;
		t_updrec upd;
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp, 0)) >= 0)
			{
				istringstream istr(tmp);
				consume += tmpsize;
				if (tmp.substr(1, 10) == "EPOCH-TIME" && _updtype != UPDTYPE::RT_SAVE)
				{
					istr >> str >> mjd >> sod;
				    _epoch.from_mjd(mjd, sod, sod-floor(sod));
				}
				else if ((tmp.substr(0, 1) == "x" || tmp.substr(0, 1) == "X") && _updtype != UPDTYPE::RT_SAVE)
				{
					// glfeng
					istr >> prn >> upd.value >> upd.sigma >> upd.npoint;
					upd.sigma = 10000.0;
					//fill data loop
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end()) {
						if (it->second->id_type() == t_gdata::UPD)
							((t_gupd*)it->second)->add_sat_upd(_updtype, _epoch, prn.substr(1), upd);
						it++;
					}
				}
				else if (tmp.substr(0, 1) == " " && _updtype != UPDTYPE::RT_SAVE)
				{
					istr >> prn >> upd.value >> upd.sigma >> upd.npoint;
					//fill data loop
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end()){
						if (it->second->id_type() == t_gdata::UPD)
						    ((t_gupd*)it->second)->add_sat_upd(_updtype, _epoch, prn, upd);
						it++;
					}
				}
				else if (tmp.substr(0, 1) == "*" && _updtype == UPDTYPE::RT_SAVE)
				{
					int yyyy, mm, dd, hh, min, ss, tt; 
					string ssss;
					istr >> ssss >> yyyy >> mm >> dd >> hh >> min >> ss >> tt;
					_epoch.from_ymdhms(yyyy, mm, dd, hh, min, ss);
				}
				else if (tmp.substr(0, 3) != "EOF"&& _updtype != UPDTYPE::RT_SAVE)
				{
					stringstream ssstr(tmp);
					t_updrec wl, nl;
					int wl_all, nl_all;
					ssstr >> prn >> wl.value >> nl.value >> wl.npoint >> wl_all >> nl.npoint >> nl_all
						>> wl.sigma >> nl.sigma;
					wl.obj = prn;
					nl.obj = prn;
					if (wl.npoint * 1.0 / wl_all >= 0.9 && nl.npoint * 1.0 / nl_all >= 0.8) {
						map<string, t_gdata*>::iterator it = _data.begin();
						while (it != _data.end()) {
							if (it->second->id_type() == t_gdata::UPD) {

								((t_gupd*)it->second)->add_sat_upd(UPDTYPE::WL, _epoch, prn, wl);

								((t_gupd*)it->second)->add_sat_upd(UPDTYPE::NL, _epoch, prn, nl);
							}
							it++;
						}
					}
				}
				else if (tmp.substr(0, 3) == "EOF")
				{
					if (_log)
					{
						_log->comment(2, "t_upd::decode_data", "WARNING: End of file" + tmp);
					}
					else
					{
						cout << "WARNING: t_upd::decode_data , End of file" << tmp << endl;;
					}
					t_gcoder::_consume(tmpsize);
					_mutex.unlock();
					return -1;
				}
				else
				{
					if (_log)
					{
						_log->comment(2, "t_upd::decode_data", "WARNING: unknown upd-data message"+tmp);
					}
					else
					{
						cout << "WARNING: t_upd::decode_data , unknown upd-data message" << tmp << endl;;
					}
					t_gcoder::_consume(tmpsize);
					_mutex.unlock();
					return -1;
				}
				t_gcoder::_consume(tmpsize);
			}
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(2, "t_upd::decode_data", "ERROR: unknown mistake");
			}
			else
			{
				cout << "ERROR: t_upd::decode_data , unknown mistake"<<endl;
			}
			return -1;
			throw(-1);
		}
		_mutex.unlock();
		return consume;
	}

	int t_upd::encode_head(char* buff, int sz, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		try
		{
			if (_ss_position == 0)
			{
				_ss << "% UPD generated using upd_";
				//get upd head info
				map<string, t_gdata*>::iterator it = _data.begin();
				for (it = _data.begin(); it != _data.end(); ++it)
				{
					if (it->second->id_type() == t_gdata::UPD)
					{
						_updtype = dynamic_cast<t_gupd*>(it->second)->get_est_updtype();
					}
				}
				if (_updtype != UPDTYPE::NONE)
				{
					_ss << updmode2str(_updtype) << endl;
				}
				else
				{
					_ss << "Wrong mode of upd: " << updmode2str(_updtype) << endl;
					int size = _fill_buffer(buff, sz);
					_mutex.unlock();
					return -1;
				}
			}
			int size = _fill_buffer(buff, sz);
			_mutex.unlock();
			return size;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(2, "t_upd::encode_head", "ERROR: unknown mistake");
			}
			else
			{
				cout << "ERROR: t_upd::encode_head , unknown mistake" << endl;
			}
			return -1;
			throw(-1);
		}
	}

	int t_upd::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		try
		{
			if (_ss_position == 0)
			{
				t_gtime valid_beg;
				//get data from _data
				map<t_gtime, one_epoch_upd> upd;
				auto it = _data.begin();
				for (it = _data.begin(); it != _data.end(); ++it)
				{
					if (it->second->id_type() == t_gdata::UPD)
					{
						upd = dynamic_cast<t_gupd*>(it->second)->get_upd(_updtype);
						valid_beg = dynamic_cast<t_gupd*>(it->second)->get_valid_beg(_updtype);
					}
				}

				// encode
				// upd
				if (_updtype == UPDTYPE::NL || _updtype == UPDTYPE::EWL_epoch || _updtype == UPDTYPE::IFCB)
				{
					auto itUpd = upd.begin();
					for (itUpd; itUpd != upd.end(); ++itUpd)
					{
						if (/*itUpd->first == WL_IDENTIFY ||*/ itUpd->first < valid_beg || itUpd->second.size() == 0)
						{
							continue;
						}
						_ss << " EPOCH-TIME" << setw(8) << itUpd->first.mjd()
							<< setw(10) << fixed << setprecision(1) << itUpd->first.sod()
							+ itUpd->first.dsec() << endl;

						auto itepo = itUpd->second.begin();
						for (itepo; itepo != itUpd->second.end(); itepo++)
						{
							if (itepo->second->npoint <= 2 || itepo->second->sigma > 0.10)
							{
								_ss << "x";
							}
							else{
								_ss << " ";
							}
							_ss << setw(11) << left << itepo->first << fixed << right << setw(10) << setprecision(3)
								<< itepo->second->value << fixed << setw(10) << setprecision(3) << itepo->second->sigma
								<< setw(5) << itepo->second->npoint << endl;
						}
					}
				}
				else if (_updtype == UPDTYPE::WL || _updtype == UPDTYPE::EWL || _updtype == UPDTYPE::EWL24 || _updtype == UPDTYPE::EWL25)
				{
					t_gtime fixed_time = WL_IDENTIFY;
					if (_updtype == UPDTYPE::EWL) fixed_time = EWL_IDENTIFY;
					else if(_updtype == UPDTYPE::EWL24) fixed_time = EWL24_IDENTIFY;
					else if (_updtype == UPDTYPE::EWL25) fixed_time = EWL25_IDENTIFY;
					auto itUpd = upd.begin();
					for (itUpd; itUpd != upd.end(); ++itUpd)
					{
						if (itUpd->first != fixed_time) continue;

						auto itepo = itUpd->second.begin();
						for (itepo; itepo != itUpd->second.end(); itepo++)
						{
							if (itepo->second->npoint <= 2 || itepo->second->sigma > 0.10)
							{
								_ss << "x";
							}
							else{
								_ss << " ";
							}
							_ss << setw(9) << left << itepo->first << fixed << right << setw(10) << setprecision(3)
								<< itepo->second->value << fixed << setw(10) << setprecision(3) << itepo->second->sigma
								<< setw(5) << itepo->second->npoint << endl;
						}
					}
				}
				else
				{
					_ss << "Wrong data" << endl;
				}

				_ss << "EOF" << endl;
			}
			//else
			//{
			int size = _fill_buffer(buff, sz);
			_mutex.unlock();
			return size;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(2, "t_upd::encode_data", "ERROR: unknown mistake");
			}
			else
			{
				cout << "ERROR: t_upd::encode_data , unknown mistake" << endl;
			}
			return -1;
			throw(-1);
		}
	}


}//namespace