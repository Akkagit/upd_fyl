#include "gcoders/ionex.h"

namespace great
{
	great::t_ionex::t_ionex(t_gsetbase* s, string version, int sz) 
		:t_gcoder(s, version, sz)
	{
		gtrace("t_ionex::constructor");
		_data_type = "";
		_nlon = 0;
		this->clear();
	}

	/** @brief destructor. */
	t_ionex::~t_ionex()
	{
		gtrace("t_ionex::destructor");
	}

	int t_ionex::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		gtrace("t_ionex::decode_head");
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
				if (line.substr(60, 18).find("EPOCH OF FIRST MAP") != string::npos)
				{
					stringstream ss(line.substr(0, 60));

					int yyyy, mm, dd, hh, min, sec;
					ss >> yyyy >> mm >> dd >> hh >> min >> sec;
					_ionex_hd.beg = t_gtime(yyyy, mm, dd, hh, min, sec);
				}
				else if (line.substr(60, 17).find("EPOCH OF LAST MAP") != string::npos)
				{
					stringstream ss(line);

					int yyyy, mm, dd, hh, min, sec;
					ss >> yyyy >> mm >> dd >> hh >> min >> sec;
					_ionex_hd.end = t_gtime(yyyy, mm, dd, hh, min, sec);

				}
				else if (line.substr(60, 8).find("INTERVAL") != string::npos)
				{
					stringstream ss(line);
					ss >> _ionex_hd.interval; // s
				}
				else if (line.substr(60, 16).find("MAPPING FUNCTION") != string::npos)
				{
					stringstream ss(line);
					string mapfunc;
					ss >> mapfunc; 
					_ionex_hd.ion_mapfunc = dynamic_cast<t_gsetproc*>(_set)->str2ionmpfunc(mapfunc);
				}
				else if (line.substr(60, 11).find("BASE RADIUS") != string::npos)
				{
					stringstream ss(line);
					ss >> _ionex_hd.base_radius;    // kilometers
					_ionex_hd.base_radius *= 1000;  // meters
				}
				else if (line.substr(60, 13).find("MAP DIMENSION") != string::npos)
				{
					stringstream ss(line);
					ss >> _ionex_hd.map_dimension;
				}
				else if (line.substr(60, 18).find("HGT1 / HGT2 / DHGT") != string::npos)
				{
					stringstream ss(line);
					double l1, l2, dl;
					ss >> l1 >> l2 >> dl;
					l1 *= 1000; l2 *= 1000; dl *= 1000;
					_ionex_hd.hgts = make_tuple(l1, l2, dl);   // meters
					if (double_eq(l1, l2)) _ionex_hd.nhgt = 1;
					else                   _ionex_hd.nhgt = round((l2 - l1) / dl) + 1;
				}
				else if (line.substr(60, 18).find("LAT1 / LAT2 / DLAT") != string::npos)
				{
					stringstream ss(line);
					double l1, l2, dl;
					ss >> l1 >> l2 >> dl;
					_ionex_hd.lats = make_tuple(l1, l2, dl);  // degree
					_ionex_hd.nlat = round((l2 - l1) / dl) + 1;
				}
				else if (line.substr(60, 18).find("LON1 / LON2 / DLON") != string::npos)
				{
					stringstream ss(line);
					double l1, l2, dl;
					ss >> l1 >> l2 >> dl;
					_ionex_hd.lons = make_tuple(l1, l2, dl);  // degree
					_ionex_hd.nlon = round((l2 - l1) / dl) + 1;
				}
				else if (line.substr(60, 8).find("EXPONENT") != string::npos)
				{
					stringstream ss(line);
					ss >> _ionex_hd.exponent;
				}
				else if (line.substr(60, 16).find("PRN / BIAS / RMS") != string::npos)
				{
				    stringstream ss(line);
					string sat;
					double bias, rms;
					ss >> sat >> bias >> rms;
					_ionex_hd.p1p2_dcb[sat] = make_pair(bias, rms);  // ns
				}
				else if (line.substr(60, 13).find("END OF HEADER") != string::npos)
				{
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::IONEX)
						{
							((t_gionex*)it->second)->add_head(_ionex_hd);
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
				_log->comment(0, "ERROR : t_ionex::decode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_ionex::decode_head throw exception" << endl;
			}
			return -1;
		}
	}


	int t_ionex::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
		gtrace("t_ionex::decode_data");
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

				if (line.find("START OF TEC MAP") != string::npos)
				{
					_data_type = "TEC_MAP";
				}
				else if (line.find("START OF RMS MAP") != string::npos)
				{
					_data_type = "RMS_MAP";
				}
				else if (line.find("END OF TEC MAP") != string::npos)
				{
					_data_type = "";
				}
				else if (line.find("END OF RMS MAP") != string::npos)
				{
					_data_type = "";
				}
				else if (line.find("EPOCH OF CURRENT MAP") != string::npos)
				{
					stringstream ss(line.substr(0, 60));

					int yyyy, mm, dd, hh, min, sec;
					ss >> yyyy >> mm >> dd >> hh >> min >> sec;
					_crt_time = t_gtime(yyyy, mm, dd, hh, min, sec);
				}
				else if (line.find("LAT/LON1/LON2/DLON/H") != string::npos)
				{
					double lat = str2dbl(line.substr(2, 6));
					double lon1 = str2dbl(line.substr(8, 6));
					double lon2 = str2dbl(line.substr(14, 6));
					double dlon = str2dbl(line.substr(20, 6));
					double height = str2dbl(line.substr(26, 6)) * 1000;  // meters

					_ilat = round((lat - get<0>(_ionex_hd.lats)) / get<2>(_ionex_hd.lats)) + 1;
					if (_ionex_hd.map_dimension == 2)  _ihgt = 1;
					else                               _ihgt   = round((height - get<0>(_ionex_hd.hgts)) / get<2>(_ionex_hd.hgts)) + 1;

					if (_ilat <= 0 || _ihgt <= 0 || !double_eq(lon1, get<0>(_ionex_hd.lons))
						|| !double_eq(lon2, get<1>(_ionex_hd.lons)) || !double_eq(dlon, get<2>(_ionex_hd.lons))) 
					{
						_nlon = -1;
					}
					else
					{
						_nlon = _ionex_hd.nlon;
						_ilon = 1;
					}
				}
				else if (_data_type == "TEC_MAP" &&  _nlon > 0)
				{
					//存数据
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::IONEX)
						{
							int line_size = line.size() - 3;  // 换行符window/Linux -3
							for (int i_size = 0; i_size < line_size; i_size = i_size + 5)
							{
								double val = str2int(line.substr(i_size, 5)) * pow(10, _ionex_hd.exponent);
								((t_gionex*)it->second)->add_data(_crt_time, "TEC_MAP", _ilat, _ilon, _ihgt, val);
								_ilon++;
								_nlon--;
							}			
						}
						it++;
					}
				}

				else if (_data_type == "RMS_MAP" &&  _nlon > 0)
				{
					//存数据
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::IONEX)
						{
							int line_size = line.size() - 3;  // 换行符window/Linux -3
							for (int i_size = 0; i_size < line_size; i_size = i_size + 5)
							{
								double rms = str2int(line.substr(i_size, 5)) * pow(10, _ionex_hd.exponent);
								((t_gionex*)it->second)->add_data(_crt_time, "RMS_MAP", _ilat, _ilon, _ihgt, rms);
								_ilon++;
								_nlon--;
							}
						}
						it++;
					}
				}
				else if (line.substr(60, 11).find("END OF FILE") != string::npos)
				{
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
				_log->comment(2, "t_ionex::decode_data", "ERROR: unknown mistake");
			}
			else
			{
				cout << "ERROR: t_ionex::decode_data , unknown mistake" << endl;
			}
			return -1;
		}
	}

}

