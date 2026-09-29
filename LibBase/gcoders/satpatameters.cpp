/**
 *
 * @verbatim
	 History
	  -1.0	lijie	2019-01-11	creat the file.
	  -1.1  lijie	2019-04-01	add code in doxygen style.
	  -1.2  lijie   2019-04-02  add abnormal throwing mechanism of functions
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		satparameters.cpp
 * @brief		The base class used to decode sat_parameters_new file information.
 *
 * @author      lijie, Wuhan University
 * @version		1.0.0
 * @date		2019-01-11
 */


#include "gcoders/satparameters.h"
//#include "gall/gallsatparam.h" //glfeng

using namespace std;

namespace great {

	/**
	* @brief constructor.
	* @param[in]  s        setbase control
	* @param[in]  version  version of the gcoder
	* @param[in]  sz       size of the buffer
	*/
	t_satparameters::t_satparameters(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_satparameters::constructor");
		_flag = 0;
	}

	/** @brief destructor. */
	t_satparameters::~t_satparameters()
	{

	}

	/**
	* @brief decode header of sat_parameters_new file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size of header decoding
	*/
	int t_satparameters::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		gtrace("t_satparameters::decode_head");
		_mutex.lock();

		if (t_gcoder::_add2buffer(buff, sz) == 0)
		{
			_mutex.unlock();
			return 0;
		}
		string tmp;
		int consume = 0;
		int tmpsize = 0;
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp)) >= 0)
			{
				if (tmp.substr(0, 1) != "#")
				{
					_mutex.unlock();
					return -1;
				}
				t_gcoder::_consume(tmpsize);
			}
			_mutex.unlock();
			return consume;
		}
		catch(...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_satparameters::decode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_satparameters::decode_head throw exception" << endl;
			}
			return -1;
		}
	}

	/**
	* @brief decode data body of sat_parameters_new file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size for data body decoding
	*/
	int t_satparameters::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
		gtrace("t_satparameters::decode_head");
		_mutex.lock();

		if (t_gcoder::_add2buffer(buff, sz) == 0)
		{
			_mutex.unlock();
			return 0;
		}
		string tmp;
		int consume = 0;
		int tmpsize = 0;
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp)) >= 0)
			{
				//judge the data type
				if (tmp.substr(0, 12) == "+prn_indexed")
				{
					_flag = 1;
				}
				else if (tmp.substr(0, 12) == "+svn_indexed")
				{
					_flag = 2;
				}
				else if (tmp.substr(0, 12) == "-svn_indexed")
				{
					_mutex.unlock();
					return consume;
				}
				//read prn_indexed data
				if (tmp.substr(0, 1) == " ")
				{
					stringstream ss(tmp);
					consume += tmpsize;
					string prn, svn, cosparid, blocktype;
					t_gtime st, et;
					double satmass = 0.0, maxyaw = 0.0;//yaw unit is deg
					int fid = 0;
					string start, end;//tmp time
					//read prn_indexed data
					if (_flag = 1)
					{
						ss >> prn >> svn >> start >> end >> cosparid >> satmass >> maxyaw >> fid >> blocktype;
						//GPS satellites' block type has two part
						if (prn.size() == 3 && prn.substr(0, 1) == "G")
						{
							string blocktype1;
							ss >> blocktype1;
							blocktype += blocktype1;
						}
					}
					//read svn_indexed data
					if (_flag = 2)
					{
						ss >> svn >> prn >> start >> end >> cosparid >> satmass >> maxyaw >> fid >> blocktype;
						//GPS satellites' block type has two part
						if (prn.size() == 3 && prn.substr(0, 1) == "G")
						{
							string blocktype1;
							ss >> blocktype1;
							blocktype += blocktype1;
						}

						//if (ss.fail())
						//{
						//	if (_log) _log->comment(1, "sat_parameters_new", "warning: incorrect sat_parameters_new data record: " + ss.str());
						//	else       cerr << "warning: incorrect sat_parameters_new data record: " << ss.str() << endl;
						//}
					}
					//transfer the time type
					int year = 0, doy = 0, sod = 0, month = 0, day = 0;
					//start time
					year = str2int(start.substr(0, 4));
					doy = str2int(start.substr(4, 3));
					sod = str2int(start.substr(8, 7));
					// day of year to month and day
					yeardoy2monthday(year, doy, &month, &day);
					st.from_ymd(year, month, day, sod, 0.0);
					//end time
					year = str2int(end.substr(0, 4));
					doy = str2int(end.substr(4, 3));
					sod = str2int(end.substr(8, 7));
					//day of year to month and day
					yeardoy2monthday(year, doy, &month, &day);
					et.from_ymd(year, month, day, sod, 0.0);
					//if the sat haven't decommissioned
					if (year == 0 && doy == 0 && sod == 0)
					{
						et = LAST_TIME;
					}
					//fill one sat data
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end()) {
						if (it->second->id_type() == t_gdata::SATPARS)
							((t_gallsatinfo*)it->second)->addData(prn, svn, st, et, cosparid, satmass, maxyaw, fid, blocktype);
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
				_log->comment(0, "ERROR : t_satparameters::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_satparameters::decode_data throw exception" << endl;
			}
			return -1;
		}
	}

	/**
	* @brief day of year to month and day
	* @param[in]   year        year
	* @param[in]   doy         day of year
	* @param[out]  month	   month
	* @param[out]  day		   day
	*/
	void t_satparameters::yeardoy2monthday(int year, int doy, int *month, int *day)
	{
		int days_inmonth[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
		int d = doy;
		if ((year % 4) == 0 && ((year % 100) != 0 || (year % 400 == 0)))
		{
			days_inmonth[1] = 29;
		}
		for (*month = 1; *month <= 12; (*month)++)
		{
			d = d - days_inmonth[(*month) - 1];
			if (d > 0) continue;
			*day = d + days_inmonth[*month - 1];
			break;
		}
	}
}