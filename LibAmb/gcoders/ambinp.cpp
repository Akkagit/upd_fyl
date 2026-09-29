/**
*
* @verbatim
History
	-1.0 Jie Li  2019-10-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		ambinp.cpp
* @brief	decode and enconde ambinp file(ambinp_yyyyddd)
*             yyyy  ---- year
*              ddd  ---- doy
*
* @author   Jie Li, Wuhan University
* @version	1.0.0
* @date		2019-10-25
*
*/
#include"gcoders/ambinp.h"

using namespace std;
namespace great
{
	t_ambinp::t_ambinp(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		_flag = 0;
		gtrace("t_ambinp::constructor");
	}

	/** @brief destructor. */
	t_ambinp::~t_ambinp()
	{
		gtrace("t_ambinp::destructor");
	}

	int t_ambinp::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		gtrace("t_ambinp::decode_data");
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
				if (line.find("%Time and Interval") != string::npos)
				{
					line.erase(0, pos + 1);
					stringstream ss(line);

					int beg_mjd, end_mjd;
					double beg_sod, end_sod, interval;
					ss >> beg_mjd >> beg_sod >> end_mjd >> end_sod >> interval;
					_ambinp_hd->beg.from_mjd(beg_mjd, (int)beg_sod, beg_sod - (int)beg_sod);
					_ambinp_hd->end.from_mjd(end_mjd, (int)end_sod, end_sod - (int)end_sod);
					_ambinp_hd->interval = interval;
				}
				else if (line.find("%Fix mode") != string::npos)
				{
					line.erase(0, pos + 1);
					stringstream ss(line);

					ss >> _ambinp_hd->fix_mode;
				}
				else if (line.find("%Total ambiguities") != string::npos)
				{
					line.erase(0, pos + 1);
					stringstream ss(line);

					ss >> _ambinp_hd->num_amb;
				}
				else if (line.find("%STA") != string::npos)
				{
					line.erase(0, pos + 1);
					stringstream ss(line);
					string site;
					while (ss >> site)
					{
						_ambinp_hd->site.insert(site);
					}
				}
				else if (line.find("%SAT") != string::npos)
				{
					line.erase(0, pos + 1);
					stringstream ss(line);
					string sat;
					while (ss >> sat)
					{
						_ambinp_hd->sat.insert(sat);
					}
				}
				else if (line.find("End Of Header") != string::npos)
				{
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::AMBINP)
						{
							((t_gambinp*)it->second)->add_hdinfo(_ambinp_hd);
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
				_log->comment(0, "ERROR : t_ambinp::decode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_ambinp::decode_head throw exception" << endl;
			}
			return -1;
		}
	}

	int t_ambinp::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
		gtrace("t_ambinp::decode_data");
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
		try
		{
			while ((tmpsize = t_gcoder::_getline(line)) >= 0)
			{
				consume += tmpsize;
				if (line.find("+Ambiguity") != string::npos)
				{
					_flag = 1;
				}
				else if (line.find("+NEQ") != string::npos)
				{
					_flag = 2;
				}

				if (_flag == 1)
				{
					if (line.find("-Ambiguity") != string::npos)
					{
						t_gcoder::_consume(tmpsize);
						_flag = 0;
						continue;
					}
					stringstream ss(line);
					int idx, beg_mjd, end_mjd;
					string sat, site;
					double beg_sod, end_sod, value, sigma;
					t_gtime beg, end;
					ss >> idx >> sat >> site >> beg_mjd >> beg_sod >> end_mjd >> end_sod >> value >> sigma;
					if (!ss.fail())
					{
						beg.from_mjd(beg_mjd, (int)beg_sod, beg_sod - (int)beg_sod);
						end.from_mjd(end_mjd, (int)end_sod, end_sod - (int)end_sod);
					}
					t_gpar par_amb(site, t_gpar::AMB_IF, idx, sat);
					par_amb.setTime(beg, end);
					par_amb.value(value);
					par_amb.apriori(sigma);
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::AMBINP)
						{
							((t_gambinp*)it->second)->add_ambpar(par_amb);
						}
						it++;
					}
				}
				else if (_flag == 2)
				{
					if (line.find("-NEQ") != string::npos)
					{
						_flag = 0;
						t_gcoder::_consume(tmpsize);
						continue;
					}
					stringstream ss(line);
					string type;
					double value;
					vector<double> q;
					ss >> type;
					while (ss >> value)
					{
						q.push_back(value);
					}
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::AMBINP)
						{
							((t_gambinp*)it->second)->add_neq(type, q);
						}
						it++;
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
				_log->comment(2, "t_ambupd::decode_data", "ERROR: unknown mistake");
			}
			else
			{
				cout << "ERROR: t_ambupd::decode_data , unknown mistake" << endl;
			}
			return -1;
			throw(-1);
		}
	}

	int t_ambinp::encode_head(char* buff, int sz, vector<string>& errmsg)
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
					if (it->second->id_type() == t_gdata::AMBINP)
					{
						//get ambinp head info
						dynamic_cast<t_gambinp*>(it->second)->get_hdinfo(_ambinp_hd);
					}
				}
				//fill head data
				_ss << "## Header Of Ambinp file" << endl;
				// time and interval
				_ss << "%Time and Interval = ";
				int beg_mjd = _ambinp_hd->beg.mjd();
				double beg_sod = _ambinp_hd->beg.sod() + _ambinp_hd->beg.dsec();
				int end_mjd = _ambinp_hd->end.mjd();
				double end_sod = _ambinp_hd->end.sod() + _ambinp_hd->end.dsec();
				_ss << beg_mjd << " " << beg_sod << "   " << end_mjd << " " << end_sod << "   " << _ambinp_hd->interval << endl;
				//fix mode
				_ss << "%Fix mode		   = " << _ambinp_hd->fix_mode << endl;
				//total ambigulity
				_ss << "%Total ambiguities = " << _ambinp_hd->num_amb << endl;
				//station number and list
				_ss << "%Station 		   = " << _ambinp_hd->site.size() << endl;
				_ss << "   %STA  = ";
				for (auto itsta = _ambinp_hd->site.begin(); itsta != _ambinp_hd->site.end(); itsta++)
				{
					_ss << *itsta;
				}
				_ss << endl;
				//satellite number and list
				_ss << "%Satellite = " << _ambinp_hd->sat.size() << endl;
				_ss << "   %PRN  = ";
				for (auto itsat = _ambinp_hd->sat.begin(); itsat != _ambinp_hd->sat.end(); itsat++)
				{
					_ss << *itsat;
				}
				_ss << endl;
				//end
				_ss << "## End Of Header" << endl;
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
				_log->comment(0, "ERROR : t_ambinp::encode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_ambinp::encode_head throw exception" << endl;
			}
			return -1;
		}
		
	}

	int t_ambinp::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		try
		{
			if (_ss_position == 0)
			{
				vector<t_gpar> ambpar;
				vector<vector<double> > neq;
				vector<string> parlist;
				map<string, t_gdata*>::iterator it = _data.begin();
				for (it = _data.begin(); it != _data.end(); ++it)
				{
					if (it->second->id_type() == t_gdata::AMBINP)
					{
						//get ambinp amb and neq
						dynamic_cast<t_gambinp*>(it->second)->get_allambpar(ambpar);
						dynamic_cast<t_gambinp*>(it->second)->get_neq(parlist, neq);
					}
				}
				//fill ambgulity data
				_ss << "+Ambiguity " << endl;
				_ss << setw(20) << right << "SAT" << "SITE" << "BEG_TIME" << "END_TIME" << "VALUE" << endl;
				for (int i = 0; i < ambpar.size(); i++)
				{
					_ss << setw(10) << right << i;
					_ss << setw(10) << right << ambpar[i].prn;
					_ss << setw(20) << right << ambpar[i].beg.mjd() << " " << ambpar[i].beg.sod() + ambpar[i].beg.dsec();
					_ss << setw(20) << right << ambpar[i].end.mjd() << " " << ambpar[i].end.sod() + ambpar[i].end.dsec();
					_ss << setw(20) << right << ambpar[i].value() << endl;
				}
				_ss << "-Ambiguity " << endl;
				//fill neq
				_ss << endl;
				_ss << "+NEQ" << endl;
				if (neq.size() != parlist.size())
				{
					throw exception();
				}
				for (int i = 0; i < neq.size(); i++)
				{
					_ss << setw(20) << left << parlist[i];
					for (int j = 0; j < neq[i].size(); j++)
					{
						_ss << setw(20) << right << neq[i][j];
					}
					_ss << endl;
				}
				_ss << "-NEQ" << endl;
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
}