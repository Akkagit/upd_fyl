/**
*
* @verbatim
History
	-1.0 Jie Li  2019-11-13  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		ambinp.cpp
* @brief	decode and enconde amb constraint file(ambinp_yyyyddd)
*             yyyy  ---- year
*              ddd  ---- doy
*
* @author   Jie Li, Wuhan University
* @version	1.0.0
* @date		2019-10-25
*
*/
#include"gcoders/ambcon.h"
#include "gambfix/gambOw.h"
#include "gdata/gambcon.h"

using namespace std;
namespace great
{
	t_ambcon::t_ambcon(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_ambcon::constructor");
	}

	/** @brief destructor. */
	t_ambcon::~t_ambcon()
	{
		gtrace("t_ambcon::destructor");
	}

	int t_ambcon::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		gtrace("t_ambcon::decode_data");
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

		int consume = 0;
		int tmpsize = 0;
		string line;
		try {
			while ((tmpsize = t_gcoder::_getline(line)) >= 0)
			{
				consume += tmpsize;
				int pos = line.find("=");
				if (line.find("%Begin Time") != string::npos)
				{
					line.erase(0, pos + 1);
					stringstream ss(line);

					int beg_mjd;
					double beg_sod;
					ss >> beg_mjd >> beg_sod;
					_beg.from_mjd(beg_mjd, (int)beg_sod, beg_sod - (int)beg_sod);
				}
				else if (line.find("%Interval") != string::npos)
				{
					line.erase(0, pos + 1);
					stringstream ss(line);
					ss >> _interval;
				}
				else if (line.find("%Ambigulities") != string::npos)
				{
					line.erase(0, pos + 1);
					stringstream ss(line);

					ss >> _num_fixamb;
				}
				else if (line.find("%PRN") != string::npos)
				{
					line.erase(0, pos + 1);
					stringstream ss(line);
					string sat;
					while (ss >> sat)
					{
						_sats.insert(sat);
					}
				}
				else if (line.find("%STA") != string::npos)
				{
					line.erase(0, pos + 1);
					stringstream ss(line);
					string site;
					while (ss >> site)
					{
						_sites.insert(site);
					}
				}
				else if (line.find("%End Of Header") != string::npos)
				{
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::AMBCON)
						{
							((t_gambcon*)it->second)->set_header(_beg, _interval, _num_fixamb, _sats, _sites);
						}
						it++;
					}

					t_gcoder::_consume(tmpsize);
					_mutex.unlock();
					return -1;
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
				_log->comment(0, "ERROR : t_ambcon::decode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_ambcon::decode_head throw exception" << endl;
			}
			return -1;
		}
	}

	int t_ambcon::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
		gtrace("t_ambcon::decode_data");
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

		int   consume = 0;
		int   tmpsize = 0;
		bool  ambvalid = false;
		string line;
		try
		{
			while ((tmpsize = t_gcoder::_getline(line)) >= 0)
			{
				consume += tmpsize;

				if (line.find("LC") != string::npos ||
					line.find("L12") != string::npos ||
					line.find("L13") != string::npos ||
					line.find("L1") != string::npos ||
					line.find("L2") != string::npos ||
					line.find("L3") != string::npos)
				{
					stringstream ss(line);
					string mark, sat1, sat2, site1, site2;
					int beg_mjd, end_mjd;
					double beg_sod, end_sod, wl, sig_wl, nl, sig_nl;
					ss >> mark >> sat1 >> sat2 >> site1 >> site2
						>> beg_mjd >> beg_sod >> end_mjd >> end_sod
						>> wl >> sig_wl >> nl >> sig_nl;
					if (!ss.fail())
					{
						t_gtime beg, end;
						beg.from_mjd(beg_mjd, (int)(beg_sod), beg_sod - (int)(beg_sod));
						end.from_mjd(end_mjd, (int)(end_sod), end_sod - (int)(end_sod));
						t_gamb_Dd dd;
						dd.beg = beg;
						dd.end = end;

						if (mark.find("LC12") != string::npos)
						{
							dd.ID = t_gamb::AMB_ID::IF_L1_L2;
						}
						else if (mark.find("LC13") != string::npos)
						{
							dd.ID = t_gamb::AMB_ID::IF_L1_L3;
						}
						else if (mark.find("L12") != string::npos)
						{
							dd.ID = t_gamb::AMB_ID::ALL_L1_L2;
						}
						else if (mark.find("L13") != string::npos)
						{
							dd.ID = t_gamb::AMB_ID::ALL_L1_L3;
						}
						else if (mark.find("L1") != string::npos)
						{
							dd.ID = t_gamb::AMB_ID::ALL_L1;
						}
						else if (mark.find("L2") != string::npos)
						{
							dd.ID = t_gamb::AMB_ID::ALL_L2;
						}
						else if (mark.find("L3") != string::npos)
						{
							dd.ID = t_gamb::AMB_ID::ALL_L3;
						}
						else
						{
							dd.ID = t_gamb::AMB_ID::UNDEF;
						}

						dd.NL_value = nl; dd.NL_sigma = sig_nl;  dd.MW_value = wl; dd.MW_sigma = sig_wl;
						dd.rec_1st = site1;
						dd.rec_2nd = site2;
						dd.sat_1st = sat1;
						dd.sat_2nd = sat2;

						if (line.substr(0, 1) == " ")
						{
							dd.MW_fixed = true;
							dd.NL_fixed = true;
						}
						else if (line.substr(0, 1) == "#")
						{
							dd.MW_fixed = true;
							dd.NL_fixed = false;

						}
						else if (line.substr(0, 1) == "*")
						{
							dd.MW_fixed = false;
							dd.NL_fixed = false;
						}

						ambvalid = true;
						if (ambvalid)
						{
							//存数据
							map<string, t_gdata*>::iterator it = _data.begin();
							while (it != _data.end())
							{
								if (it->second->id_type() == t_gdata::AMBCON)
								{
									((t_gambcon*)it->second)->add_ddamb(dd);
								}
								it++;
							}
							ambvalid = false;
						}
						else
						{
							_mutex.unlock();
							consume -= tmpsize;
							return consume;
						}
					}
				}
				else if (line.find("#START ALL") != string::npos)
				{
					t_gamb_Dd dd;
					int ambsize = _decoder_data_UCUD(tmpsize, line, ambvalid, dd);
					if (ambvalid)
					{
						//存数据
						map<string, t_gdata*>::iterator it = _data.begin();
						while (it != _data.end())
						{
							if (it->second->id_type() == t_gdata::AMBCON)
							{
								((t_gambcon*)it->second)->add_ddamb(dd);
							}
							it++;
						}
						tmpsize += ambsize;
						consume += ambsize;
						ambvalid = false;
					}
					else
					{
						consume -= tmpsize;
						_mutex.unlock();
						return consume;
					}
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
				_log->comment(2, "t_ambcon::decode_data", "ERROR: unknown mistake");
			}
			else
			{
				cout << "ERROR: t_ambcon::decode_data , unknown mistake" << endl;
			}
			return -1;
			throw(-1);
		}
	}

	int t_ambcon::encode_head(char* buff, int sz, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		try
		{
			if (_ss_position == 0)
			{
				map < pair<GSYS, t_gamb::AMB_ID>, int> all_WL_num;
				map < pair<GSYS, t_gamb::AMB_ID>, int> fix_WL_num;
				map < pair<GSYS, t_gamb::AMB_ID>, int> all_NL_num;
				map < pair<GSYS, t_gamb::AMB_ID>, int> fix_NL_num;

				map<string, t_gdata*>::iterator it = _data.begin();
				for (it = _data.begin(); it != _data.end(); ++it)
				{
					if (it->second->id_type() == t_gdata::AMBCON)
					{
						//get ambinp head info
						dynamic_cast<t_gambcon*>(it->second)->get_header(_beg, _interval, _num_fixamb, _sats, _sites);
						all_WL_num = dynamic_cast<t_gambcon*>(it->second)->all_WL_num;
						fix_WL_num = dynamic_cast<t_gambcon*>(it->second)->fix_WL_num;
						all_NL_num = dynamic_cast<t_gambcon*>(it->second)->all_NL_num;
						fix_NL_num = dynamic_cast<t_gambcon*>(it->second)->fix_NL_num;
					}
				}
				//fill head data
				_ss << "## Header Of Ambcon File" << endl;
				// beg time 
				_ss << "%Begin Time 		= ";
				int beg_mjd = _beg.mjd();
				double beg_sod = _beg.sod() + _beg.dsec();
				_ss << beg_mjd << "   " << beg_sod << endl;
				//interval(seconds)
				_ss << "%Interval   		= " << _interval << endl;
				// ambigulity number
				_ss << "%Ambigulities 	 	= " << _num_fixamb << endl;

				//satellite number and list
				_ss << "%Satellites 		= " << _sats.size() << endl;
				_ss << "   %PRN  = ";
				int j = 0;
				for (auto itsat = _sats.begin(); itsat != _sats.end(); itsat++)
				{
					j++;
					_ss << *itsat << " ";
					if (j >= 20)
					{
						j = 0;
						_ss << endl;
						_ss << "   %PRN  = ";
					}
				}
				_ss << endl;

				//station number and list
				_ss << "%Stations   	  	= " << _sites.size() << endl;
				_ss << "   %STA  = ";
				int i = 0;
				for (auto itsta = _sites.begin(); itsta != _sites.end(); itsta++)
				{
					i++;
					_ss << *itsta << " ";
					if (i >= 20)
					{
						i = 0;
						_ss << endl;
						_ss << "   %STA  = ";
					}
				}
				_ss << endl;
				//end

				// write number
				for (auto& num_amb : all_WL_num)
				{
					_ss << "all_WL_num " << t_gsys::gsys2str(num_amb.first.first) << " " << t_gamb::ambId2str(num_amb.first.second) << " " << num_amb.second << endl;
				}
				for (auto& num_amb : fix_WL_num)
				{
					_ss << "fix_WL_num " << t_gsys::gsys2str(num_amb.first.first) << " " << t_gamb::ambId2str(num_amb.first.second) << " " << num_amb.second << endl;
				}
				for (auto& num_amb : all_NL_num)
				{
					_ss << "all_NL_num " << t_gsys::gsys2str(num_amb.first.first) << " " << t_gamb::ambId2str(num_amb.first.second) << " " << num_amb.second << endl;
				}
				for (auto& num_amb : fix_NL_num)
				{
					_ss << "fix_NL_num " << t_gsys::gsys2str(num_amb.first.first) << " " << t_gamb::ambId2str(num_amb.first.second) << " " << num_amb.second << endl;
				}

				_ss << "%End Of Header" << endl;
				_ss << endl;
			}
			int size = _fill_buffer(buff, sz);
			_mutex.unlock();
			return size;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_ambcon::encode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_ambcon::encode_head throw exception" << endl;
			}
			return -1;
		}

	}

	int t_ambcon::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		try
		{
			if (_ss_position == 0)
			{
				vector<t_gamb_Dd> ddambs;
				for (auto it = _data.begin(); it != _data.end(); ++it)
				{
					if (it->second->id_type() == t_gdata::AMBCON)
					{
						//get ambCON DDambs
						dynamic_cast<t_gambcon*>(it->second)->get_Dds(ddambs);
					}
				}

				// 格式说明
				_ss << "------SAT1--SAT2--SITE1-SITE2------------BEGIN TIME------------------END TIME--------------WL VALUE--------WL SIGMA-------NL VALUE---------NL SIGMA----" << endl;

				// 为了实现排序，写三次
				string mark;
				for (const auto& itamb : ddambs)
				{
					if (itamb.ID == t_gamb::AMB_ID::IF_L1_L2)
					{
						mark = " LC12";
					}
					else if (itamb.ID == t_gamb::AMB_ID::IF_L1_L3)
					{
						mark = " LC13";
					}
					else if (itamb.ID == t_gamb::AMB_ID::ALL_L1_L2)
					{
						mark = " L12";
					}
					else if (itamb.ID == t_gamb::AMB_ID::ALL_L1_L3)
					{
						mark = " L13";
					}
					else if (itamb.ID == t_gamb::AMB_ID::ALL_L1)
					{
						mark = " L1";
					}
					else if (itamb.ID == t_gamb::AMB_ID::ALL_L2)
					{
						mark = " L2";
					}
					else if (itamb.ID == t_gamb::AMB_ID::ALL_L3)
					{
						mark = " L3";
					}
					else
					{
						continue;
					}

					if (itamb.NL_fixed && itamb.MW_fixed)
					{
						_ss << setw(6) << left << mark
							<< setw(6) << itamb.sat_1st
							<< setw(6) << itamb.sat_2nd
							<< setw(6) << itamb.rec_1st
							<< setw(6) << itamb.rec_2nd;

						_ss << setw(10) << right << itamb.beg.mjd() << " "
							<< setw(16) << fixed << setprecision(8) << itamb.beg.sod() + itamb.beg.dsec();
						_ss << setw(10) << right << itamb.end.mjd() << " "
							<< setw(16) << fixed << setprecision(8) << itamb.end.sod() + itamb.end.dsec();

						_ss << setw(16) << right << fixed << setprecision(8) << itamb.MW_value
							<< setw(16) << right << fixed << setprecision(8) << itamb.MW_sigma
							<< setw(16) << right << fixed << setprecision(8) << itamb.NL_value
							<< setw(16) << right << fixed << setprecision(8) << itamb.NL_sigma << endl;

					}
				}

				for (auto& itamb : ddambs)
				{
					if (itamb.ID != t_gamb::AMB_ID::IF_L1_L2) break;
					if (!(itamb.NL_fixed && itamb.MW_fixed) && (itamb.MW_fixed))
					{
						mark = "#LC";
						_ss << setw(6) << left << mark
							<< setw(6) << itamb.sat_1st
							<< setw(6) << itamb.sat_2nd
							<< setw(6) << itamb.rec_1st
							<< setw(6) << itamb.rec_2nd;

						_ss << setw(10) << right << itamb.beg.mjd() << " "
							<< setw(16) << fixed << setprecision(8) << itamb.beg.sod() + itamb.beg.dsec();
						_ss << setw(10) << right << itamb.end.mjd() << " "
							<< setw(16) << fixed << setprecision(8) << itamb.end.sod() + itamb.end.dsec();

						_ss << setw(16) << right << fixed << setprecision(8) << itamb.MW_value
							<< setw(16) << right << fixed << setprecision(8) << itamb.MW_sigma
							<< setw(16) << right << fixed << setprecision(8) << itamb.NL_value
							<< setw(16) << right << fixed << setprecision(8) << itamb.NL_sigma << endl;
					}
				}


				for (auto& itamb : ddambs)
				{
					if (itamb.ID != t_gamb::AMB_ID::IF_L1_L2) break;
					if (!itamb.MW_fixed)
					{
						mark = "*LC";
						_ss << setw(6) << left << mark
							<< setw(6) << itamb.sat_1st
							<< setw(6) << itamb.sat_2nd
							<< setw(6) << itamb.rec_1st
							<< setw(6) << itamb.rec_2nd;

						_ss << setw(10) << right << itamb.beg.mjd() << " "
							<< setw(16) << fixed << setprecision(8) << itamb.beg.sod() + itamb.beg.dsec();
						_ss << setw(10) << right << itamb.end.mjd() << " "
							<< setw(16) << fixed << setprecision(8) << itamb.end.sod() + itamb.end.dsec();

						_ss << setw(16) << right << fixed << setprecision(8) << itamb.MW_value
							<< setw(16) << right << fixed << setprecision(8) << itamb.MW_sigma
							<< setw(16) << right << fixed << setprecision(8) << itamb.NL_value
							<< setw(16) << right << fixed << setprecision(8) << itamb.NL_sigma << endl;
					}
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
				_log->comment(0, "ERROR : t_orb::encode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_orb::encode_data throw exception" << endl;
			}
			return -1;
		}

	}

	int t_ambcon::_decoder_data_UCUD(int tmpsize, string& line, bool& valid, t_gamb_Dd& dd)
	{
		string sLC = trim(line.erase(0, line.find("#START ALL") + 10));
		int   ambtemp = 0;
		int   ambsize = 0;
		while ((ambtemp = t_gcoder::_getline(line, tmpsize + ambsize)) >= 0)
		{
			ambsize += ambtemp;

			if (line.find("#END   ALL") != string::npos)
			{
				valid = true;
				break;
			}

			string type;
			int    numb;
			string sat1;
			string sat2;
			string rec1;
			string rec2;
			int beg_mjd;
			int end_mjd;
			double beg_sod;
			double end_sod;
			double value;
			double sigma;

			stringstream ss(line);
			if (line.find("DD") != string::npos)
			{
				ss >> type >> numb
					>> sat1 >> sat2
					>> rec1 >> rec2
					>> value >> sigma
					>> beg_mjd >> beg_sod
					>> end_mjd >> end_sod;

				if (!ss.fail())
				{
					t_gtime beg;
					t_gtime end;
					beg.from_mjd(beg_mjd, (int)(beg_sod), beg_sod - (int)(beg_sod));
					end.from_mjd(end_mjd, (int)(end_sod), end_sod - (int)(end_sod));
					dd.beg = beg;
					dd.end = end;
					dd.rec_1st = rec1;
					dd.rec_2nd = rec2;
					dd.sat_1st = sat1;
					dd.sat_2nd = sat2;
					dd.value = value;
					dd.sigma = sigma;
					dd.fixed = true;

					if (sLC.find("L1") != string::npos) dd.ID = t_gamb::AMB_ID::ALL_L1;
					else if (sLC.find("L2") != string::npos) dd.ID = t_gamb::AMB_ID::ALL_L2;
					else if (sLC.find("L3") != string::npos) dd.ID = t_gamb::AMB_ID::ALL_L3;
					else if (sLC.find("L4") != string::npos) dd.ID = t_gamb::AMB_ID::ALL_L4;
					else if (sLC.find("L5") != string::npos) dd.ID = t_gamb::AMB_ID::ALL_L5;
					else dd.ID = t_gamb::AMB_ID::UNDEF;

				}
			}

			if (line.find("OW") != string::npos)
			{
				ss >> type >> numb
					>> sat1
					>> rec1
					>> value >> sigma
					>> beg_mjd >> beg_sod
					>> end_mjd >> end_sod;
				if (!ss.fail())
				{
					t_gtime beg;
					t_gtime end;
					beg.from_mjd(beg_mjd, (int)(beg_sod), beg_sod - (int)(beg_sod));
					end.from_mjd(end_mjd, (int)(end_sod), end_sod - (int)(end_sod));
					dd.Ow_time[numb] = make_pair(beg, end);
					dd.Ow_value[numb] = value;
				}
			}

		}

		return ambsize;
	}

}