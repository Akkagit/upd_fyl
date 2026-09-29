/**
*
* @verbatim
History
-1.0 BoWong  2019-02-22  creat the file.
-1.1 BoWong  2019-04-09  Adding Doxygen Style Code Remarks
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		ambupd.cpp
* @brief	decode and enconde ambupd file(XXXX_ambupd_yyyyddd)
*             XXXX  ---- site name
*             yyyy  ---- year
*              ddd  ---- doy
*
* @author   BoWong, Wuhan University
* @version	1.0.0
* @date		2019-04-09
*
*/
#include"gcoders/ambupd.h"

using namespace std;
namespace great
{
	t_ambupd::t_ambupd(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_ambupd::constructor");
	}

	/** @brief destructor. */
	t_ambupd::~t_ambupd()
	{
		gtrace("t_ambupd::destructor");
	}

	int t_ambupd::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		gtrace("t_ambupd::decode_data");
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

		_mutex.unlock();
		return -1;
	}

	int t_ambupd::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
		gtrace("t_ambupd::decode_data");
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
			  string sta;
			  string prn;
			  t_ambupd_data data_tmp;
			  t_gtime epoch;
			  int mjd;
			  double sod;
			  try
			  {
				  while ((tmpsize = t_gcoder::_getline(tmp, 0)) >= 0)
				  {
					  consume += tmpsize;
					  istringstream istr(tmp);
					  istr >> mjd >> sod >> sta >> prn >> data_tmp.amb_if >> data_tmp.amb_wl >> data_tmp.amb_l3 >> data_tmp.amb_l4 >> data_tmp.amb_l5 >> data_tmp.flag >> data_tmp.amb_wl_sigma;
					  epoch.from_mjd(mjd, sod, sod - floor(sod));

					  map<string, t_gdata*>::iterator it = _data.begin();
					  while (it != _data.end())
					  {
						  if (it->second->id_type() == t_gdata::AMBUPD)
						  {
							  ((t_gallambupd*)it->second)->addAmbUpd(sta, epoch, prn, data_tmp);
						  }
						  it++;
					  }
					  t_gcoder::_consume(tmpsize);
				  }
			  }
			  catch (...)
			  {
				  if (_log)
				  {
					  _log->comment(2, "t_ambupd::decode_data", "ERROR: unknown mistake");
				  }
				  else
				  {
					  cout << "ERROR: t_ambupd::decode_data , unknown mistake"<<endl;
				  }
				  return -1;
				  throw(-1);
			  }
			  _mutex.unlock();
			  return tmpsize;
	}

	int t_ambupd::encode_head(char* buff, int sz, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();

		int size = _fill_buffer(buff, sz);
		_mutex.unlock();
		return size;
	}

	int t_ambupd::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
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
				t_epo_ambupd    data_tmp;
				string site("    ");
				auto it = _data.begin();
				for (it = _data.begin(); it != _data.end(); ++it)
				{
					if (it->second->id_type() == t_gdata::AMBUPD)
					{
						for (const auto site_ambupd : dynamic_cast<t_gallambupd*>(it->second)->get_sites()) {
							data_tmp = dynamic_cast<t_gallambupd*>(it->second)->getOneSiteAmbUpd(site_ambupd).getAllAmbUpd();
							site = site_ambupd;
							break;
						}
					}
				}

				// encode
				auto itepo = data_tmp.begin();
				for (itepo; itepo != data_tmp.end(); ++itepo)
				{
					t_gtime epo = itepo->first;
					auto itsat = itepo->second.begin();
					for (itsat; itsat != itepo->second.end(); itsat++)
					{
						// i8,f10.1,1x,a4,1x,a3,2f19.3,i8,f10.3
						_ss << fixed << setw(8) << epo.mjd() << setw(10) << setprecision(1) << epo.sod() + 0.0
							<< " " << setw(4) << site << " " << setw(3) << itsat->first
							<< setw(19) << setprecision(3) << itsat->second->amb_if
							<< setw(19) << setprecision(3) << itsat->second->amb_wl
							<< setw(19) << setprecision(3) << itsat->second->amb_l3
							<< setw(19) << setprecision(3) << itsat->second->amb_l4
							<< setw(19) << setprecision(3) << itsat->second->amb_l5
							<< setw(8) << itsat->second->flag << setw(10) << setprecision(3)
							<< itsat->second->amb_wl_sigma << endl;
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
				_log->comment(2, "t_ambupd::encode_data", "ERROR: unknown mistake");
			}
			else
			{
				cout << "ERROR: t_ambupd::encode_data , unknown mistake" << endl;
			}
			return -1;
			throw(-1);
		}
	}
}