#include"gcoders/aug.h"

using namespace std;
namespace great
{
	t_aug::t_aug(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
	}


	t_aug::~t_aug()
	{
	}

	int t_aug::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		_mutex.lock();
		if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };
		string tmp;
		int consume = 0;
		int tmpsize = 0;
		while ((tmpsize = t_gcoder::_getline(tmp)) >= 0){
			consume += tmpsize;
			istringstream istr(tmp);
			if (tmp.find("SYS / # / AUG TYPES") != string::npos)
			{
				string s_tmp, str_sys, str_num;
				istr >> str_sys >> str_num;
				GSYS sys = t_gsys::str2gsys(str_sys);
				vector<t_augtype> augtype;
				while (istr >> s_tmp)
				{
					AUGTYPE type(AUGTYPE::TYPE_NON);
					string s_type = s_tmp.substr(0, 1);
					if (s_type == "P") {
						type = AUGTYPE::TYPE_P;
					}
					else if (s_type == "L") {
						type = AUGTYPE::TYPE_L;
					}
					else if (s_type == "I") {
						type = AUGTYPE::TYPE_ION;
					}
					else if (s_type == "T") {
						type = AUGTYPE::TYPE_TRP;
					}
					else break;
					string str_band = s_tmp.substr(s_tmp.size() - 1, s_tmp.size() - 0);//lvhb test
					GOBSBAND band = (GOBSBAND)std::atoi(str_band.c_str());
					augtype.push_back(make_pair(type, band));
				}
				_augtype[sys] = augtype;
			}
			else if (tmp.find("MARKER NAME") != string::npos) {
				istr >> _site;
			}
			else if (tmp.find("END OF HEADER") != string::npos)
			{
				if (_log) _log->comment(2, "aug", "END OF HEADER");
				t_gcoder::_consume(tmpsize);//lvhb added in 202007
				_mutex.unlock(); return -1;
			}
			t_gcoder::_consume(tmpsize);
		}
		_mutex.unlock(); return consume;
	}

	int t_aug::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
		_mutex.lock();

		if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

		string tmp;
		int consume = 0;
		int tmpsize = 0;
		
		while ((tmpsize = t_gcoder::_getline(tmp, 0)) >= 0) {
			consume += tmpsize;
			istringstream istr(tmp);
			if (tmp.substr(0, 1) == ">") {
				string str;
				int year, month, day, hour, min;
				double second;

				istr >> str >> year >> month >> day >> hour >> min >> second;
				_epoch.from_ymdhms(year, month, day, hour, min, second);
			}
			else if (t_gsys::str2gsys(tmp.substr(0, 1)) != GSYS::GNS) {
				string sat;
				istr >> sat;
				GSYS sys = t_gsys::sat2gsys(sat);
				vector<t_augtype> augtype = _augtype[sys];
				for (auto iter = augtype.begin(); iter != augtype.end(); iter++)
				{
					int pos = distance(augtype.begin(), iter);
					//string str_value = tmp.substr(3 + pos * 12, 15 + pos * 12);
					string str_value; istr >> str_value;//lvhb changed in npp
					if (str_value.empty())continue;
					double value = atof(str_value.c_str());

					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::AUG)
						{
							((t_gaug*)it->second)->add_data(_site, _epoch, sat, *iter, value);
							break;
						}
						it++;
					}
				}
			}
			else
			{
				if (_log) _log->comment(1, "aug", "warning: incorrect AUG data record: " + _ss.str());
				else       cerr << "warning: incorrect AUG data record: " << _ss.str() << endl;
				t_gcoder::_consume(tmpsize);
				_mutex.unlock(); return -1;
			}
			t_gcoder::_consume(tmpsize);
		}
		_mutex.unlock(); return consume;
	}
}