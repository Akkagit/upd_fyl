/**
*
* @verbatim
	History
	 -1.0 Yujie Qin	 2020-02-25 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			desaiscopolecoef.cpp
* @brief		Decode kin data
*
* @author		Yujie Qin, Wuhan University
* @version		1.0.0
* @date			2019-10-16
*
*/
#include "gcoders/desaiscopolecoef.h"
#include "gdata/gdeadata.h"

using namespace std;
using namespace gnut;
namespace great {

	t_desaiscopolecoef::t_desaiscopolecoef(t_gsetbase* s, string version, int sz)
		: t_gcoder(s, version, sz)
	{
		gtrace("t_desaiscopolecoef::constructor");
	}
	
	t_desaiscopolecoef::~t_desaiscopolecoef()
	{

	}

	
	

	int t_desaiscopolecoef::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		if (t_gcoder::_add2buffer(buff, sz) == 0) { return 0; };
		return -1;
	}

	int t_desaiscopolecoef::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

#ifdef DEBUG
		cout << " BUFFER : \n" << _buffer << "\n size = " << sz << " END OF BUFFER \n\n"; cout.flush();
#endif
		//t_oceantide oct[IMAXCONSTITS];
		string line;
		string _line;
		int consume = 0;
		int tmpsize = 0;
		int sitsize = 0;
		int n, m;
		double Anm, Bnm, AInm, BInm;
		int num;
		double load;
		


		tmpsize = t_gcoder::_getline(line);
		sitsize += tmpsize;
		consume += tmpsize;
		

		//find the starting place
		/*while (line.find("+coefficient") == string::npos ) {
			t_gcoder::_consume(tmpsize);
			tmpsize = t_gcoder::_getline(line, sitsize);
			sitsize += tmpsize;
			consume += tmpsize;
			complete = false;

			}*/

		
		try
		{
			//find the ending place
			while ((tmpsize = t_gcoder::_getline(line)) >= 0)
			{
				sitsize += tmpsize;
				consume += tmpsize;

				if (line.find("-coefficient") != string::npos)
				{
					_tab = 1;
				}

				if (line.substr(0, 1) == " "&&_tab == 0)
				{
	
					
					istringstream istr(line);

					istr >> n >> m >> Anm >> Bnm >> AInm >> BInm;
					
						
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::DESAISCOPOLECOEF)
						{
							((t_gdeadata*)it->second)->set_rec(n,m, Anm, Bnm, AInm, BInm);
							
						}

						it++;
					}
						
				}

				

				if (line.substr(0, 1) == " "&&_tab == 1)
				{
					if (line.substr(0, 1) == " ")
					{
						string line3;
						istringstream istr(line);
						istr >> num >> load;
						
						
						if (!istr.fail())
						{
							//
							//t_oceantd.setload(_load);
							map<string, t_gdata*>::iterator it = _data.begin();
							while (it != _data.end())
							{
								if (it->second->id_type() == t_gdata::DESAISCOPOLECOEF)
								{
									((t_gdeadata*)it->second)->set_load(load);
								}

								it++;
							}
						}
					}

				}


				t_gcoder::_consume(tmpsize);
			}


			
			_mutex.unlock(); return sitsize;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_desaiscopolecoef::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_desaiscopolecoef::decode_data throw exception" << endl;
			}
			return -1;
		}
	}
	
	


}
