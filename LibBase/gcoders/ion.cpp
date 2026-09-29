#include "gcoders/ion.h"
#include "gdata/gion.h"

namespace great
{
	great::t_ion::t_ion(t_gsetbase* s, string version, int sz) 
		:t_gcoder(s, version, sz)
	{
		gtrace("t_ion::constructor");
	}

	/** @brief destructor. */
	t_ion::~t_ion()
	{
		gtrace("t_ion::destructor");
	}

	int t_ion::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		gtrace("t_ion::decode_head");
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
				else if (line.find("%End   Time") != string::npos)
				{
					line.erase(0, pos + 1);
					stringstream ss(line);

					int end_mjd;
					double end_sod;
					ss >> end_mjd >> end_sod;
					_end.from_mjd(end_mjd, (int)end_sod, end_sod - (int)end_sod);

				}
				else if (line.find("%Interval"  ) != string::npos)
				{
					line.erase(0, pos + 1);
					stringstream ss(line);
					ss >> _intv;
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
					stringstream ss(line);
					line.erase(0, pos + 1);
					string site;
					while (ss >> site)
					{
						_recs.insert(site);
					}
				}
				else if (line.find("%End Of Header") != string::npos)
				{
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::ION)
						{
							((t_gion*)it->second)->set_header(_beg, _end, _intv, _numb, _sats, _recs);
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


	int t_ion::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
		gtrace("t_ion::decode_data");
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

				if (line.find("END OF FILE") != string::npos)
				{
					break;
				}

				stringstream ss(line);
				string mark;
				string rec;
				string sat;
				int beg_mjd = 0;
				int end_mjd = 0;
				double beg_sod = 0.0;
				double end_sod = 0.0;
				double value = 0.0;
				double sigma = 0.0;
				ss >> mark >> rec >> sat >> value >> sigma
					>> beg_mjd >> beg_sod >> end_mjd >> end_sod;
				if (!ss.fail())
				{
					t_gtime beg, end;
					beg.from_mjd(beg_mjd, (int)(beg_sod), beg_sod - (int)(beg_sod));
					end.from_mjd(beg_mjd, (int)(end_sod), end_sod - (int)(end_sod));

					t_gion_record record(mark, rec, sat, value, sigma, beg, end);

					//´æÊý¾Ý
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::ION)
						{
							((t_gion*)it->second)->add_record(record);
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

	int t_ion::encode_head(char* buff, int sz, vector<string>& errmsg)
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
					if (it->second->id_type() == t_gdata::ION)
					{
						//get ambinp head info
						dynamic_cast<t_gion*>(it->second)->get_header(_beg, _end, _intv, _sats, _recs);
					}
				}

				//fill head data
				_ss << "## Header Of Ion File" << endl;
				// beg time 
				_ss << "%Begin Time 		= ";
				_ss << _beg.str_mjdsod() << endl;
				_ss << "%End   Time 		= ";
				_ss << _end.str_mjdsod() << endl;
				//interval(seconds)
				_ss << "%Interval   		= " << _intv << endl;
				//satellite number and list
				_ss << "%Satellites 		= " << _sats.size() << endl;
				_ss << "   %PRN  = ";
				for (auto itsat = _sats.begin(); itsat != _sats.end(); itsat++)
				{
					_ss << *itsat << " ";
				}
				_ss << endl;
				//station number and list
				_ss << "%Stations   	  	= " << _recs.size() << endl;
				_ss << "   %STA  = ";
				int i = 0;
				for (auto itsta = _recs.begin(); itsta != _recs.end(); itsta++)
				{
					i++;
					if (i > 20)
					{
						i = 0;
						_ss << endl;
						_ss << "   %STA  = ";
					}
					_ss << *itsta << " ";
				}
				_ss << endl;
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
				_log->comment(0, "ERROR : t_ambcon::encode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_ambcon::encode_head throw exception" << endl;
			}
			return -1;
		}

	}

	int t_ion::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		try
		{
			if (_ss_position == 0)
			{
				map< pair<string, string>, map<t_gtime, t_gion_record> > tmp;
				for (auto it = _data.begin(); it != _data.end(); ++it)
				{
					if (it->second->id_type() == t_gdata::ION)
					{
						for (const auto& it_rec : _recs)
						{
							for (const auto& it_sat : _sats)
							{
								auto key = make_pair(it_rec, it_sat);
								//get ambCON DDambs
								tmp[key] = dynamic_cast<t_gion*>(it->second)->get_records(it_rec,it_sat);
							}
						}		
					}
				}

				_ss << "TYPE--SITE--SAT--------VALUE------------SIGMA------BEGIN TIME-------END TIME--------------" << endl;

				string mark;
				for (const auto& it_pair : tmp)
				{
					for (const auto& it_record : it_pair.second)
					{
						_ss << setw(6) << left << it_record.second.ion_type
							<< setw(6) << it_record.second.rec
							<< setw(6) << it_record.second.sat
							<< setw(16) << right << fixed << setprecision(8) << it_record.second.value
							<< setw(16) << right << fixed << setprecision(8) << it_record.second.sigma
							<< setw(16) << left << it_record.second.beg.str_mjdsod()
							<< setw(16) << left << it_record.second.end.str_mjdsod()
							<< endl;
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
				_log->comment(0, "ERROR : t_ion::encode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_ion::encode_head throw exception" << endl;
			}
			return -1;
		}

	}

}

