
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
	  Research Institute of Geodesy, Topography and Cartography
	  Ondrejov 244, 251 65, Czech Republic

  This file is part of the G-Nut C++ library.

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 3 of
  the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, see <http://www.gnu.org/licenses>.

-*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <memory>

#include "gcoders/bncobs.h"
#include "gio/gtcp.h"
#include "gutils/gconst.h"
#include "gutils/gnss.h"
#include "gutils/gsys.h"
#include "gall/gallobj.h"

using namespace std;

namespace gnut {

	// constructor
	t_bncobs::t_bncobs(t_gsetbase* s, string version, int sz)
		: t_gcoder(s, version, sz)
	{
		_begepoch = false;
		//  _endepoch = false;

		_tt = FIRST_TIME;
	}

	// ---------
	// synchronized real-time observation in BNC structure (BNC:bnchelp.html)
	// ---------
	int t_bncobs::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();

		// no header expected, but fill the buffer
		t_gcoder::_add2buffer(buff, sz);
		_mutex.unlock(); return -1;
	}


	// ---------
	// synchronized real-time observation in BNC structure (BNC:bnchelp.html)
	// ---------
	int t_bncobs::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();

		//   if( _log ) _log->comment(2,"bncobs","Start decode_data");

		  // complete gcoder buffer
		if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; }

		//  cout << " BUFFER : \n" << _buffer << "\n size = " << sz  << " END OF BUFFER \n\n"; cout.flush();

		string line;
		int tmpsize = 0;     // individual reading counter
		int consume = 0;     // total read counter

		while ((tmpsize = t_gcoder::_getline(line, 0)) >= 0) {
			//cout << line << endl;
			if (line.find(crlf) == string::npos) return 0; // end of buffer

			if (line.compare(crlf) == 0) {  // empty line
				consume += t_gcoder::_consume(tmpsize);
				continue;
			}

			istringstream ss(line);

			string c, site;
			int gw, satnum;
			double sow;
			char satsys;

			if (line.compare(0, 1, ">") == 0) { // read epoch
				ss.clear();
				ss >> c >> gw >> sow;

				this->available_true();
				//cout << line << endl;       
				//cout << "gw: " << gw << " sow: " << sow << endl;
				if (!ss.fail()) {
					t_gtime t(t_gtime::GPS);
					t.from_gws(gw, sow);
					_tt = t;
					//		cout << _tt.str_ymdhms() << " " << ss.str() << endl;
					consume += t_gcoder::_consume(tmpsize);
					continue;
				}
			}

			if (!_validepo(_tt)) {
				consume += t_gcoder::_consume(tmpsize);
				continue;
			}

			ss.clear();
			ss >> site >> satsys >> satnum;
			if (!ss.fail() && (satsys == 'G' ||
				satsys == 'R' ||
				satsys == 'E' ||
				satsys == 'C' ||
				satsys == 'S' ||
				satsys == 'J')) {
				string prn = t_gsys::eval_sat(int2str(satnum), t_gsys::char2gsys(satsys));
				GSYS gsys = t_gsys::char2gsys(satsys);

				// filtering REC
				string site_4ch = site.substr(0, 4);
				site = site.substr(0, 4);
				if (_rec.size() == 0 || (_rec.find(site) == _rec.end() && _rec.find(site_4ch) == _rec.end())) {
					//        cout << "Filtering out Site " << site << endl;
					consume += t_gcoder::_consume(tmpsize);
					continue;
				}

				shared_ptr<t_gobsgnss> obs = make_shared<t_gobsgnss>(site, prn, _tt);

				map<GOBS, double> gobs;
				map<GOBS, int> glli;

				// checking if receiver object exists
				_check_recobj(site);

				// filtering GNSS systems and satellites
				if ((_sys.size() > 0 && _sys.find(t_gsys::char2str(satsys)) == _sys.end()) ||
					(_sat[gsys].size() > 0 && _sat[gsys].find(prn) == _sat[gsys].end())) {
					//        cout << "Filtering out sat/sys " << prn << endl;
					consume += t_gcoder::_consume(tmpsize);
					continue;
				}

				/*
				if (satsys == 'R'){
				  int slotnum = 999;
				  ss.clear();
				  ss >> slotnum;
				  obs->channel(slotnum);
				}
				*/

				string obstyp;
				double observ;
				int lli;
				ss.clear();

				while (ss >> obstyp >> observ && !ss.fail()) {
					if (gsys == BDS && obstyp[1] == '1')  obstyp[1] = '2'; // change B1 -> B2 !!!
					GOBS typ = str2gobs(obstyp);
					gobs[typ] = observ;

					if (obstyp.compare(0, 1, "L") == 0) {
						ss >> lli;
						glli[typ] = lli;
					}
				}

				for (auto it = gobs.begin(); it != gobs.end(); ++it)
				{
					if (_obs[gsys].size() == 0 ||                 // filtering signals
						_obs[gsys].find(trim(gobs2str(it->first))) != _obs[gsys].end())
					{
						// cout << "Signal passed filter " << gobs2str(it->first) << endl;
						if (it->second != 0.0) { obs->addobs(it->first, it->second); }
					}
				}
				/*
					  // if BNC provides non-zero always values (cycles slips), PPP has problem
					  // ----------------------------------------------------------------------
					  for( auto it = glli.begin(); it != glli.end(); ++it )
					  {
						if( _obs[gsys].size() == 0 ||                // filtering signals
							_obs[gsys].find(trim(gobs2str(it->first))) != _obs[gsys].end() )
						{
						  obs->addlli(it->first, it->second );
						}
					  }
				*/
#ifdef DEBUG       
				cout << fixed << setprecision(3)
					<< " DECODED: " << site
					<< " " << prn
					<< " C1C: " << setw(16) << obs->getobs(C1C)
					<< " L1C: " << setw(16) << obs->getobs(L1C)
					<< " S1C: " << setw(16) << obs->getobs(S1C)
					<< obs->epoch().str(" %Y-%m-%d %H:%M:%S[%T] %W %w %v %J %s") << "\n";
				//    int ooo; cin >> ooo;
#endif

				map<string, t_gdata*>::iterator it;

				// adding GLO channel number to gobsgnss
				it = _data.begin();
				while (it != _data.end()) {
					if (it->second->id_type() != t_gdata::ALLOBJ) {
						it++;
						continue;
					}
					t_gallobj* all_obj = (t_gallobj*)it->second;
					shared_ptr<t_gtrn> one_obj = dynamic_pointer_cast<t_gtrn>(all_obj->obj(prn));
					if (one_obj != 0 && one_obj->id_type() == t_gdata::TRN) {
						int ch = one_obj->channel();
						obs->channel(ch);
					}
					it++;
				}

				// fill data
				it = _data.begin();
				while (it != _data.end()) {
					if (it->second->id_type() == t_gdata::ALLOBS && _validepo(_tt))
					{
						((t_gallobs*)it->second)->addobs(obs);

						if (_log) {
							ostringstream ltmp;
							double diff = t_gtime::current_time(t_gtime::GPS) - obs->epoch();
							ltmp << "add site observations " << site << " " << obs->sat() << " "
								<< " " << obs->epoch().str("%Y-%m-%d %H:%M:%S[%T]") << " Diff " << setprecision(2) << diff;
							_log->comment(3, "bncobs", ltmp.str());
						}
					}
					it++;
				}
			}
			cnt++;

			consume += t_gcoder::_consume(tmpsize);

		} // while
		_mutex.unlock(); return consume;

	} //end method


	/* -------- *
	 * ENCODER
	 * -------- */

	int t_bncobs::encode_head(char* buff, int sz, vector<string>& errmsg)
	{
		return 0;
	}

	// encode data
	// ---------
	int t_bncobs::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();

		// fill data
		if (_ss_position == 0)
		{
			for (auto it = _data.begin(); it != _data.end(); ++it)
			{
				// data types
				if (it->second->id_type() != t_gdata::ALLOBS &&
					it->second->id_group() != t_gdata::GRP_OBSERV) {
					continue;
				}

				t_gallobs* gallobs = dynamic_cast<t_gallobs*>(it->second);

				set<string> sites = gallobs->stations();
				vector<t_gtime> epochs = gallobs->epochs(*sites.begin());

				if (sites.size() == 0 || epochs.size() == 0) {
					if (_log) _log->comment(2, "bncobs", "no data to be saved "
						+ dbl2str(epochs.size()) + dbl2str(sites.size()));
					continue;
				}

				// epochs
				for (auto itEPO = epochs.begin(); itEPO != epochs.end(); ++itEPO)
				{
					t_gtime epo = *itEPO;
					if (epo < _beg || epo > _end) { continue; } // filter BEG/END
					if (_out_smp > 0 && int(round(epo.sod() + epo.dsec())) % int(_out_smp) != 0) { continue; } // output SAMPLING

					_ss << ">" << fixed << setw(5) << epo.gwk()
						<< setprecision(7) << setw(15) << epo.sow() + epo.dsec()
						<< setw(22) << epo.str_ymdhms()
						<< endl;

					// stations
					for (auto itSIT = sites.begin(); itSIT != sites.end(); ++itSIT)
					{
						string site = *itSIT;

						// satellites/observations
						vector<t_spt_gobs> vec_obs = gallobs->obs_pt(site, epo);
						for (auto itOBS = vec_obs.begin(); itOBS != vec_obs.end(); ++itOBS)
						{
							t_spt_gobs pt_gobs = *itOBS;
							string sat = pt_gobs->sat();

							_ss << site << setw(4) << sat;

							// signals
							vector<GOBS> vec_sig = pt_gobs->obs();
							for (auto itSIG = vec_sig.begin(); itSIG != vec_sig.end(); ++itSIG)
							{
								GOBS  gobs = *itSIG;
								t_gobs go(gobs);
								string obs = gobs2str(gobs);
								if (obs[0] == 'S') {
									_ss << setw(4) << obs << setprecision(3) << setw(9) << pt_gobs->getobs(gobs);
								} // SNR
								else { _ss << setw(4) << obs << setprecision(3) << setw(15) << pt_gobs->getobs(gobs); } // OBS
								if (go.is_phase()) { _ss << setprecision(0) << setw(5) << pt_gobs->getlli(gobs); } // LLI
							} // signals
							_ss << endl;
						} // satellites/observations
					} // stations
				} // epochs
				_ss << endl;
			} // data type
		} // condition

		int size = _fill_buffer(buff, sz);

		_mutex.unlock(); return size;
	}

	bool t_bncobs::available(const t_gtime& now)
	{


		if (now <= _tt)return true;
		else
		{
			//return _available;
			return false;
		}

	}

	void t_bncobs::available_false()
	{
		gtrace("t_gtcp::available_false");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		//_mutex.lock();

		_available = false;

		//_mutex.unlock();
	}

	void t_bncobs::available_true()
	{
		gtrace("t_gtcp::available_true");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		//_mutex.lock();

		_available = true;

		//_mutex.unlock();
	}

	bool t_bncobs::_validepo(const t_gtime& t)
	{
		if (t > _end || t < _beg) return false;

		if (_filter_epoch(t)) return true;
		else return false;
	}

	// Checking if a receiver object exists
	// ----------------------------------
	void t_bncobs::_check_recobj(string& site)
	{
		map<string, t_gdata*>::iterator it = _data.begin();
		for (it = _data.begin(); it != _data.end(); ++it) {
			if (it->second->id_type() == t_gdata::ALLOBJ) {
				t_gallobj* all_obj = (t_gallobj*)it->second;
				shared_ptr<t_gobj> rec = all_obj->obj(site);
				if (rec == nullptr) {
					rec = make_shared<t_grec>();
					rec->id(site);
					rec->name(site);
					all_obj->add(rec);
				}
			}
		}

	}

} // namespace
