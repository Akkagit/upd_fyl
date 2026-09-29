/**
*
* @verbatim
	History
	 -1.0 Jie Li	 2019-01-06 creat the file.
	 -1.1 Jie Li     2019-04-01 Adding Doxygen Style Code Remarks
	 -1.2 Jie Li	 2019-04-02 Adding abnormal throwing mechanism of functions
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			egm2008.cpp
* @brief		Decode EGM data
*
* @author		lijie, Wuhan University
* @version		1.0.0
* @date			2019-01-06
*
*/

//todo move noto to header file

#include "gcoders/egm2008.h"
#include "gdata/gegmcoef.h"

using namespace std;

namespace great {

	/**
	* @brief constructor.
	* @param[in]  s        setbase control
	* @param[in]  version  version of the gcoder
	* @param[in]  sz       size of the buffer
	*/
	t_egm2008::t_egm2008(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_egm2008::constructor");
	}

	/** @brief destructor. */
	t_egm2008:: ~t_egm2008(){}

	/**
	* @brief decode header of EGM file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size of header decoding
	*/
	int t_egm2008::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		gtrace("t_egm2008::decode_head");

		_mutex.lock();

		if (t_gcoder::_add2buffer(buff, sz) == 0)
		{ 
			_mutex.unlock(); 
			return 0;
		}

		string tmp;
		int consume = 0;
		int tmpsize = 0;
		double gm = 0.0;
		double radius = 0.0;
		int i = 0, year = 0, mon = 0, day = 0;
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp)) >= 0)
			{
				stringstream ss(tmp);
				consume += tmpsize;
				if (tmp.find("+coefficient") == -1)
				{
					if (tmp.substr(0, 5) == "EARTH")
					{
						string str;
						gm = 0.0;
						radius = 0.0;
						i = 0;
						year = 0;
						mon = 0;
						day = 0;

						t_gtime mjd_ref;

						ss >> str >> gm >> radius >> i;
						gm = gm * 1e-9;
						radius = radius * 1e-3;
						year = i / 10000;
						mon = (i - year * 10000) / 100;
						day = i - year * 10000 - mon * 100;
						mjd_ref.from_ymd(year, mon, day, 0, 0.0);

						//fill time gm and radius
						map<string, t_gdata*>::iterator it = _data.begin();
						while (it != _data.end())
						{
							if (it->second->id_type() == t_gdata::EGM)
							{
								((t_gegmcoef*)it->second)->add_head(mjd_ref, gm, radius);
							}
							it++;
						}
					}

					if (tmp.substr(0, 3) == "SHM")
					{
						string str;
						int order = 0.0;
						ss >> str >> order;

						//fill max order of egm file
						map<string, t_gdata*>::iterator it = _data.begin();
						while (it != _data.end())
						{
							if (it->second->id_type() == t_gdata::EGM)
							{
								((t_gegmcoef*)it->second)->add_order(order);
							}
							it++;
						}
					}
					t_gcoder::_consume(tmpsize);
				}
				else
				{
					if (_log)
						_log->comment(2, "t_egm2008::decode_head end of head");
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
				_log->comment(0, "ERROR : t_egm2008::decode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_egm2008::decode_heada throw exception" << endl;
			}
			return -1;
		}
	}

	/**
	* @brief decode data body of EGM file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size for data body decoding
	*/
	int t_egm2008::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
		gtrace("t_egm2008::decode_data");
		_mutex.lock();

		if (t_gcoder::_add2buffer(buff, sz) == 0){ _mutex.unlock(); return 0; };
		string tmp;
		int consume = 0;
		int tmpsize = 0;
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp)) >= 0)
			{
				stringstream ss(tmp);
				consume += tmpsize;
				if (tmp.find("-coefficient") == -1)
				{
					int n = 0, m = 0;
					double c = 0.0, s = 0.0, sigma_c = 0.0, sigma_s = 0.0;
					ss >> n >> m >> c >> s >> sigma_c >> sigma_s;
					//fill coeff data
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end()) {
						if (it->second->id_type() == t_gdata::EGM)
						{
							((t_gegmcoef*)it->second)->add_egmcoef(n, m, c, s, sigma_c, sigma_s);
						}
						it++;
					}
				}
				else
				{
					//cout << "end of data" << endl;
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
				_log->comment(0, "ERROR : t_egm2008::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_egm2008::decode_data throw exception" << endl;
			}
			return -1;
		}
	}
}