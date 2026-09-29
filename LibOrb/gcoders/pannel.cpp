/**
*
* @verbatim
    History
     -1.0 Yujie Qin	    2019-01-06 creat the file.
     -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
	 -1.2 Jie Li 	    2019-04-02 add abnormal throwing mechanism of functions
	 -1.3 Yujie Qin     2019-07-02 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		   pannel.cpp
* @brief		   Decode the LEO pannel data
*
* @author      Yujie Qin, Wuhan University
* @version	   1.0.0
* @date		   2019-03-27
*
*/

#include "gcoders/pannel.h"

using namespace std;
namespace great {

	/** @brief default constructor. */
	t_pannel::t_pannel(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_pannel::constructor");
	}

	/** @brief default destructor. */
	t_pannel::~t_pannel()
	{

	}

	/**
	* @brief decode header of LEO pannel file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size of header decoding
	*/
	int t_pannel::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		_mutex.lock();

		if (t_gcoder::_add2buffer(buff, sz) == 0){ _mutex.unlock(); return 0; };
		
		// jdhuang : unlock when return
		_mutex.unlock();
		return -1;
	}

	/**
	* @brief decode data body of LEO pannel file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size of data body decoding
	*/
	int t_pannel::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		//_mutex.lock();

		if (t_gcoder::_add2buffer(buff, sz) == 0){ _mutex.unlock(); return 0; };
#ifdef DEBUG
		cout << " BUFFER : \n" << _buffer << "\n size = " << sz << " END OF BUFFER \n\n"; cout.flush();
#endif
		string tmp;
		int tmpsize = 0;
		int num;
		t_oneside nu;
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp, 0)) >= 0)
			{
				if (tmp.substr(0, 1) == "+")
				{
					istringstream istr(tmp.substr(1));
					istr >> _sat_name;
					_sat_pannel.add_name(_sat_name);
				}
				else if (tmp.substr(0, 1) == " ")
				{

					istringstream istr(tmp);
					istr >> num >> nu.area >> nu.scf_norm_dir[0] >> nu.scf_norm_dir[1] >> nu.scf_norm_dir[2]
						>> nu.sr_coef >> nu.dr_coef >> nu.drag_coef >> nu.isr_coef >> nu.idr_coef;
					if (tmp.size() > 76) nu.name = tmp.substr(75, tmp.size()-76);
					_sat_pannel.add_pannel(num, nu);
				}
				else if (tmp.substr(0, 1) == "-")
				{
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::ALLPANNEL)
							((t_gallpannel*)it->second)->add_sat_pannel(_sat_name, _sat_pannel);
						it++;
					}
					_sat_pannel.clear();
				}
				t_gcoder::_consume(tmpsize);
			}
			_mutex.unlock(); return tmpsize;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_pannel::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_pannel::decode_data throw exception" << endl;
			}
			return -1;
		}
	}
}