/**
*
* @verbatim
History
	-1.0 Yujie Qin  2020-01-06  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		ambcon_leo.cpp
* @brief	decode and enconde amb constraint file(ambinp_yyyyddd)
*             yyyy  ---- year
*              ddd  ---- doy
*
* @author   Yujie Qin, Wuhan University
* @version	1.0.0
* @date		2020-01-06
*
*/
#include"gcoders/ambcon_leo.h"
#include "gambfix/gambOw.h"
#include "gdata/gambcon.h"

using namespace std;
namespace great
{
	t_ambcon_leo::t_ambcon_leo(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_ambcon_leo::constructor");
	}

	/** @brief destructor. */
	t_ambcon_leo::~t_ambcon_leo()
	{
		gtrace("t_ambcon_leo::destructor");
	}

	int t_ambcon_leo::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		gtrace("t_ambcon_leo::decode_data");
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
					stringstream ss(line);
					line.erase(0, pos + 1);
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

	int t_ambcon_leo::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
		gtrace("t_ambcon_leo::decode_data");
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
		bool  ambvalid = false;
		try
		{
			while ((tmpsize = t_gcoder::_getline(line)) >= 0)
			{
				consume += tmpsize;

				if (line.find("LC") != string::npos)
				{
					stringstream ss(line);
					string mark, sat1, sat2, site;
					int beg_mjd, end_mjd;
					double beg_sod, end_sod, wl, sig_wl, nl, sig_nl,rnl;
					ss >> mark >> sat1 >> sat2 >> site
						>> beg_mjd >> beg_sod >> end_mjd >> end_sod
						>> wl >> sig_wl >> nl >> sig_nl>>rnl;
					if (!ss.fail())
					{
						t_gtime beg, end;
						beg.from_mjd(beg_mjd, (int)(beg_sod), beg_sod - (int)(beg_sod));
						end.from_mjd(end_mjd, (int)(end_sod), end_sod - (int)(end_sod));
						t_gamb_Sd sd;
						sd.beg = beg;
						sd.end = end;
						sd.ID = t_gamb::AMB_ID::IF_L1_L2;

						sd.NL_value = nl; sd.NL_sigma = sig_nl;  sd.MW_value = wl; sd.MW_sigma = sig_wl;
						sd.sgl_obj = site;
						sd.dbl_obj_1st = sat1;
						sd.dbl_obj_2nd = sat2;
						sd.rnl = rnl;

						//
						t_gamb_Dd dd;
						dd.beg = beg;
						dd.end = end;
						dd.ID = t_gamb::AMB_ID::IF_L1_L2;

						dd.NL_value = nl; dd.NL_sigma = sig_nl;  dd.MW_value = wl; dd.MW_sigma = sig_wl;
						dd.rec_1st = site;
						dd.rec_2nd = "    ";
						dd.sat_1st = sat1;
						dd.sat_2nd = sat2;
						dd.rnl = rnl;

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

						//

						if (line.substr(0, 1) == " ")
						{
							sd.MW_fixed = true;
							sd.NL_fixed = true;
						}
						else if (line.substr(0, 1) == "#")
						{
							sd.MW_fixed = true;
							sd.NL_fixed = false;

						}
						else if (line.substr(0, 1) == "*")
						{
							sd.MW_fixed = false;
							sd.NL_fixed = false;
						}

						//存数据
						map<string, t_gdata*>::iterator it = _data.begin();
						while (it != _data.end())
						{
							if (it->second->id_type() == t_gdata::AMBCON)
							{
								((t_gambcon*)it->second)->add_sdamb(sd);
							}
							it++;
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
				else if (line.find("ALL") != string::npos)
				{
					stringstream ss(line);
					string mark, sLC, sat1, sat2, site;
					int beg_mjd, end_mjd;
					double beg_sod, end_sod, value, sigma, rnl;

					double Ow1 = 0.0;
					double Ow2 = 0.0;

					ss >> mark >> sLC >> sat1 >> sat2 >> site
						>> beg_mjd >> beg_sod >> end_mjd >> end_sod
						>> Ow1 >> Ow2 >> value >> sigma >> rnl;
					if (!ss.fail())
					{
						t_gtime beg, end;
						beg.from_mjd(beg_mjd, (int)(beg_sod), beg_sod - (int)(beg_sod));
						end.from_mjd(end_mjd, (int)(end_sod), end_sod - (int)(end_sod));
						t_gamb_Sd sd;
						sd.beg = beg;
						sd.end = end;
						
						sd.Ow_value[1] = Ow1;
						sd.Ow_value[2] = Ow2;
						
						
						sd.sgl_obj = site;
						sd.dbl_obj_1st = sat1;
						sd.dbl_obj_2nd = sat2;
						sd.ID = t_gamb::AMB_ID::IF_L1_L2;

						sd.value = value; sd.sigma = sigma; sd.rnl = rnl;

						if (line.substr(0, 1) == " ")
						{
							sd.fixed = true;
						}
						else
						{
							sd.fixed = false;
						}

						     if (sLC == "L1") sd.ID = t_gamb::AMB_ID::ALL_L1;
						else if (sLC == "L2") sd.ID = t_gamb::AMB_ID::ALL_L2;
						else if (sLC == "L3") sd.ID = t_gamb::AMB_ID::ALL_L3;
						else if (sLC == "L4") sd.ID = t_gamb::AMB_ID::ALL_L4;
						else if (sLC == "L5") sd.ID = t_gamb::AMB_ID::ALL_L5;
						else sd.ID = t_gamb::AMB_ID::UNDEF;

						//存数据
						map<string, t_gdata*>::iterator it = _data.begin();
						while (it != _data.end())
						{
							if (it->second->id_type() == t_gdata::AMBCON)
							{
								((t_gambcon*)it->second)->add_sdamb(sd);
							}
							it++;
						}
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
				_log->comment(2, "t_ambcon_leo::decode_data", "ERROR: unknown mistake");
			}
			else
			{
				cout << "ERROR: t_ambcon_leo::decode_data , unknown mistake" << endl;
			}
			return -1;
			throw(-1);
	}
}

	int t_ambcon_leo::encode_head(char* buff, int sz, vector<string>& errmsg)
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
				double beg_sod = _beg.sod() +_beg.dsec();
				_ss << beg_mjd << "   " << beg_sod << endl;
				//interval(seconds)
				_ss << "%Interval   		= " <<_interval << endl;
				// ambigulity number
				_ss << "%Ambigulities 	 	= " << _num_fixamb << endl;
				//satellite number and list
				_ss << "%Satellites 		= " <<_sats.size() << endl;
				_ss << "   %PRN  = ";
				for (auto itsat =_sats.begin(); itsat != _sats.end(); itsat++)
				{
					_ss << *itsat<<" ";
				}
				_ss << endl;
				//station number and list
				_ss << "%Stations   	  	= " <<_sites.size() << endl;
				_ss << "   %STA  = ";
				int i = 0;
				for (auto itsta = _sites.begin(); itsta != _sites.end(); itsta++)
				{
					i++;
					if (i > 20)
					{
						i = 0;
						_ss << endl;
						_ss<< "   %STA  = ";
					}
					_ss << *itsta << " ";
				}
				_ss << endl;
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
				//end
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
				_log->comment(0, "ERROR : t_ambcon_leo::encode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_ambcon_leo::encode_head throw exception" << endl;
			}
			return -1;
	}

	}

	int t_ambcon_leo::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		try
		{
			if (_ss_position == 0)
			{
				vector<t_gamb_Sd> sdambs;
				for (auto it = _data.begin(); it != _data.end(); ++it)
				{
					if (it->second->id_type() == t_gdata::AMBCON)
					{
						//get ambCON DDambs
						dynamic_cast<t_gambcon*>(it->second)->get_Sds(sdambs);
					}
				}

				for (const auto& itamb : sdambs)
				{
					if (itamb.ID == t_gamb::AMB_ID::IF_L1_L2)
					{
						_ss << "------SAT1--SAT2--SITE1-SITE2------------BEGIN TIME------------------END TIME----------------WL VALUE----------WL SIGMA---------NL VALUE-----------NL SIGMA----" << endl;
					}
					else
					{
						_ss << "------SAT1--SAT2--SITE1-SITE2------------BEGIN TIME------------------END TIME----------------WL VALUE----------WL SIGMA---------NL VALUE-----------NL SIGMA----" << endl;
					}
					break;
				}
				// 为了实现排序，写三次
				string mark;
				int begnum;
				int endnum;
				/*map<string, int> satt;
				satt["G01"] = 1;
				satt["G02"] = 2;
				satt["G03"] = 3;
				satt["G04"] = 4;
				satt["G05"] = 5;
				satt["G06"] = 6;
				satt["G07"] = 7;
				satt["G08"] = 8;
				satt["G09"] = 9;
				satt["G10"] = 10;
				satt["G11"] = 11;
				satt["G12"] = 12;
				satt["G13"] = 13;
				satt["G14"] = 14;
				satt["G15"] = 15;
				satt["G16"] = 16;
				satt["G17"] = 17;
				satt["G19"] = 18;
				satt["G20"] = 19;
				satt["G21"] = 20;
				satt["G22"] = 21;
				satt["G23"] = 22;
				satt["G24"] = 23;
				satt["G25"] = 24;
				satt["G26"] = 25;
				satt["G27"] = 26;
				satt["G28"] = 27;
				satt["G29"] = 28;
				satt["G30"] = 29;
				satt["G31"] = 30;
				satt["G32"] = 31;*/


				for (const auto& itamb : sdambs)
				{
					if (itamb.ID == t_gamb::AMB_ID::IF_L1_L2)
					{
						if (itamb.NL_fixed && itamb.MW_fixed)
						{
							begnum = (itamb.beg - _beg) /10 +1;
							endnum= (itamb.end - _beg) / 10 +1;
							mark = " LC";
							_ss << setw(6) << left << mark
								<< setw(6) << itamb.dbl_obj_1st
								<< setw(6) << itamb.dbl_obj_2nd
								<< setw(6) << itamb.sgl_obj
								<< setw(6) << "----";

							_ss << setw(10) << right << itamb.beg.mjd() << " "
								<< setw(16) << fixed << setprecision(8) << itamb.beg.sod() + itamb.beg.dsec();
							_ss << setw(10) << right << itamb.end.mjd() << " "
								<< setw(16) << fixed << setprecision(8) << itamb.end.sod() + itamb.end.dsec();

							_ss << setw(20) << right << fixed << setprecision(8) << itamb.MW_value
								<< setw(16) << right << fixed << setprecision(8) << itamb.MW_sigma
								<< setw(20) << right << fixed << setprecision(8) << itamb.NL_value
								<< setw(16) << right << fixed << setprecision(8) << itamb.NL_sigma 
								<< setw(16) << right << fixed << setprecision(8) << itamb.rnl << endl;
							
							/*cout<< setw(3) << left << mark
								<<"   1"
								<< setw(4) << right << satt[itamb.dbl_obj_1st]
								<<"   0"
								<< setw(4) << right << satt[itamb.dbl_obj_2nd];

							cout << setw(7) << right << begnum << " ";
							cout << setw(7) << right << endnum << " ";

							cout << setw(18) << right << fixed << setprecision(5) << itamb.WL_value
								<< setw(10) << right << fixed << setprecision(5) << itamb.WL_sigma
								<< setw(18) << right << fixed << setprecision(5) << itamb.NL_value
								<< setw(10) << right << fixed << setprecision(5) << itamb.NL_sigma 
								<< setw(10) << right << fixed << setprecision(5) <<itamb.rnl<< endl;*/

						}
					}
					else
					{
						mark = " ALL";
						string sLC;
						if (itamb.ID == t_gamb::AMB_ID::ALL_L1) sLC = "L1";
						else if (itamb.ID == t_gamb::AMB_ID::ALL_L2) sLC = "L2";
						else if (itamb.ID == t_gamb::AMB_ID::ALL_L3) sLC = "L3";
						else if (itamb.ID == t_gamb::AMB_ID::ALL_L4) sLC = "L4";
						else if (itamb.ID == t_gamb::AMB_ID::ALL_L5) sLC = "L5";
						else sLC = "UK";

						/*_ss << setw(6) << left << mark
							<< setw(6) << left << sLC
							<< setw(6) << itamb.dbl_obj_1st
							<< setw(6) << itamb.dbl_obj_2nd
							<< setw(6) << itamb.sgl_obj
							<< setw(6) << "----";

						_ss << setw(10) << right << itamb.beg.mjd() << " "
							<< setw(16) << fixed << setprecision(8) << itamb.beg.sod() + itamb.beg.dsec();
						_ss << setw(10) << right << itamb.end.mjd() << " "
							<< setw(16) << fixed << setprecision(8) << itamb.end.sod() + itamb.end.dsec();

						_ss << setw(16) << right << fixed << setprecision(8) << itamb.Ow_value.find(1)->second
						    << setw(16) << right << fixed << setprecision(8) << itamb.Ow_value.find(2)->second
							<< setw(16) << right << fixed << setprecision(8) << itamb.value
							<< setw(16) << right << fixed << setprecision(8) << itamb.sigma 
							<< setw(16) << right << fixed << setprecision(8) << itamb.rnl << endl;*/
					}
				}

				for (auto& itamb : sdambs)
				{
					if (itamb.ID != t_gamb::AMB_ID::IF_L1_L2) break;
					if (!(itamb.NL_fixed && itamb.MW_fixed) && (itamb.MW_fixed))
					{
						mark = "#LC";
						_ss << setw(6) << left << mark
							<< setw(6) << itamb.dbl_obj_1st
							<< setw(6) << itamb.dbl_obj_2nd
							<< setw(6) << itamb.sgl_obj
							<< setw(6) << "----";

						_ss << setw(10) << right << itamb.beg.mjd() << " "
							<< setw(16) << fixed << setprecision(8) << itamb.beg.sod() + itamb.beg.dsec();
						_ss << setw(10) << right << itamb.end.mjd() << " "
							<< setw(16) << fixed << setprecision(8) << itamb.end.sod() + itamb.end.dsec();

						_ss << setw(16) << right << fixed << setprecision(8) << itamb.MW_value
							<< setw(16) << right << fixed << setprecision(8) << itamb.MW_sigma
							<< setw(16) << right << fixed << setprecision(8) << itamb.NL_value
							<< setw(16) << right << fixed << setprecision(8) << itamb.NL_sigma 
							<< setw(16) << right << fixed << setprecision(8) << itamb.rnl << endl;
					}
				}
					

				for (auto& itamb : sdambs)
				{
					if (itamb.ID != t_gamb::AMB_ID::IF_L1_L2) break;
					if (!itamb.MW_fixed)
					{
						mark = "*LC";
						_ss << setw(6) << left << mark
							<< setw(6) << itamb.dbl_obj_1st
							<< setw(6) << itamb.dbl_obj_2nd
							<< setw(6) << itamb.sgl_obj
							<< setw(6) << "----";

						_ss << setw(10) << right << itamb.beg.mjd() << " "
							<< setw(16) << fixed << setprecision(8) << itamb.beg.sod() + itamb.beg.dsec();
						_ss << setw(10) << right << itamb.end.mjd() << " "
							<< setw(16) << fixed << setprecision(8) << itamb.end.sod() + itamb.end.dsec();

						_ss << setw(16) << right << fixed << setprecision(8) << itamb.MW_value
							<< setw(16) << right << fixed << setprecision(8) << itamb.MW_sigma
							<< setw(16) << right << fixed << setprecision(8) << itamb.NL_value
							<< setw(16) << right << fixed << setprecision(8) << itamb.NL_sigma 
							<< setw(16) << right << fixed << setprecision(8) << itamb.rnl << endl;
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
				_log->comment(0, "ERROR : t_ambcon_leo::encode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_ambcon_leo::encode_data throw exception" << endl;
			}
			return -1;
		}

	}
}