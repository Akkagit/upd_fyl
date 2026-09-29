/**
*
* @verbatim
History
-1.0 xjhan  2019-11-22  creat the file.

@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		rinexc_leo.cpp
* @brief	read rinexc containing leo data.
*
* @author   xjhan , Wuhan University
* @version	1.0.0
* @date		2019-11-22
*
*/

#include "rinexc_leo.h"
#include "gall/gallprec.h"

namespace gnut {

	t_rinexc_leo::t_rinexc_leo(t_gsetbase* s, string version, int sz)
		:t_rinexc(s, version, sz)
	{
	}

	t_rinexc_leo::~t_rinexc_leo()
	{
	}

	int t_rinexc_leo::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
		gtrace(_class + "::decode_data");

		t_gtime epoch(t_gtime::GPS);       // epoch
		double data[6];                    // clk/rms values (off,vel,acc)
		string  flg, id, line;
		int tmpsize = 0;                   // individual reading counter
		int consume = 0;                   // total read counter

		for (int i = 0; i < 6; i++) data[i] = UNDEFVAL_CLK; // initialize

		_mutex.lock();

		// complete main buffer
		if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

		while ((tmpsize = t_gcoder::_getline(line)) >= 0) {

			istringstream istr(line);
			istr.clear();

			int yr, mn, dd, hr, mi;
			double sc;
			int ncol = 1;

			// new record
			istr >> flg;

			if (flg.find("AS") != string::npos) {
				char sat[3 + 1]; sat[3] = '\0';
				char dummy;
				istr >> noskipws >> dummy >> sat[0] >> sat[1] >> sat[2]
					>> skipws >> yr >> mn >> dd >> hr >> mi >> sc >> ncol;
				id = t_gsys::eval_sat_addleo(string(sat));
			}
			else { istr >> id >> yr >> mn >> dd >> hr >> mi >> sc >> ncol; }

			if (istr.fail()) break;

			epoch.from_ymd(yr, mn, dd, (hr * 3600 + mi * 60), sc);

			// loop over n-column
			for (int i = 0; i < ncol; i++) {

				// first line
				istr >> data[i];
				if (istr.fail()) break;

				// additional (second) line
				if (i > 2) {
					if ((tmpsize += t_gcoder::_getline(line)) < 0) break;

					istr >> data[i];
					if (istr.fail()) break;
				}
			}

			// ZHJ:add for recclk
			if ((!_data.empty()) && (flg.find("AS") != string::npos || flg.find("AR") != string::npos)) {

				double clk[3] = { data[0], data[2], data[4] };
				double var[3] = { data[0], data[2], data[4] };

				// fill clk
				map<string, t_gdata*>::iterator it = _data.begin();
				while (it != _data.end()) {
					if (_log && _log->verb() >= 3) {
						ostringstream ltmp;
						ltmp << "add sat clk " << id
							<< epoch.str_ymdhms(" ")
							<< scientific
							<< setw(14) << clk[0] << setw(14) << var[0]
							<< setw(14) << clk[1] << setw(14) << var[1]
							<< setw(14) << clk[2] << setw(14) << var[2];
						_log->comment(3, "rinexc", ltmp.str());
					}
					if (it->second->id_type() == t_gdata::ALLPREC) {
						((t_gallprec*)it->second)->addclk(id, epoch, clk, var);
					}

					it++;
				}
				cnt++;
			}
			consume += t_gcoder::_consume(tmpsize);
		}

		_mutex.unlock(); return consume;
	}


}