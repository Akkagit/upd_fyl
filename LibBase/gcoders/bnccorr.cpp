
/* ----------------------------------------------------------------------
* @author ShenZhiheng
* @version 1.0
* @date 2019-11-10
-*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <memory>

#include "bnccorr.h"
#include "gio/gtcp.h"
#include "gutils/gconst.h"
#include "gutils/gnss.h"
#include "gutils/gsys.h"
#include "gall/gallobj.h"
#include "gdata/gupd.h"
#include "gdata/gaug.h"

using namespace std;
using namespace great;

namespace gnut {

	// constructor
	t_bnccorr::t_bnccorr(t_gsetbase* s, string version, int sz)
		: t_gcoder(s, version, sz)
	{
		_begepoch = false;
		//  _endepoch = false;

		_tt = FIRST_TIME;

		_band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(s)->band_index(gnut::GPS);
		_band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(s)->band_index(gnut::GAL);
		_band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(s)->band_index(gnut::GLO);
		_band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(s)->band_index(gnut::BDS);
	}

	// ---------
	// synchronized real-time observation in BNC structure (BNC:bnchelp.html)
	// ---------
	int t_bnccorr::decode_head(char* buff, int sz, vector<string>& errmsg)
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
	int t_bnccorr::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();

		//   if( _log ) _log->comment(2,"bncobs","Start decode_data");

		  // complete gcoder buffer
		if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; }

		//  cout << " BUFFER : \n" << _buffer << "\n size = " << sz  << " END OF BUFFER \n\n"; cout.flush();
		string flag;
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

			string c, site, stream;
			int yr, mn, dd, hr, mi, indicator, size;
			int gpsw, nsat;

			float sec, gpssec;

			if (line.compare(0, 1, ">") == 0) { // read epoch
				ss.clear();
				//ss >> c >> flag >> yr >> mn >> dd >> hr >> mi >> sec >> indicator >> size >> stream;
				ss >> c >> flag;
				_type = str2corr(flag);
				this->available_true();

				if     (_type == UPD_COR) ss >> yr >> mn >> dd >> hr >> mi >> sec >> nsat;
				else if(_type == AUG_COR) ss >> gpsw >> gpssec >> site;
				else                      ss >> yr >> mn >> dd >> hr >> mi >> sec >> indicator >> size >> stream;

				//cout << line << endl;       
				//cout << "gw: " << gw << " sow: " << sow << endl;
				if (!ss.fail()) {
					t_gtime t(t_gtime::GPS);
					if (_type != AUG_COR) {
						t.from_ymdhms(yr, mn, dd, hr, mi, sec);
					}
					else {
						t.from_gws(gpsw, gpssec);
						_decode_aug(site, t, ss);
						cnt++;
					}
					_tt = t;
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end()) {
						if (it->second->id_type() == t_gdata::ALLPREC)
						{
							if (_type == ORBIT) { ((t_gallprec*)it->second)->add_orb_interval(intv[indicator]); _udorbInt = intv[indicator]; }
							if (_type == CLOCK) { ((t_gallprec*)it->second)->add_clk_interval(intv[indicator]); _udclkInt = intv[indicator]; }
						}
						it++;
					}
					consume += t_gcoder::_consume(tmpsize);
					continue;
				}

			}

			if (!_validepo(_tt)) {
				consume += t_gcoder::_consume(tmpsize);
				continue;
			}

			ss.clear();
			int tmp;

			if (_type == ORBIT)tmp = _decode_orbit(_tt,ss);
			else if (_type == CLOCK)tmp = _decode_clock(_tt,ss);
			else if (_type == CODE_BIAS)tmp = _decode_cb(_tt,ss);
			else if (_type == PHASE_BIAS)tmp = _decode_pb(_tt,ss);
			else if (_type == VTEC)tmp = _decode_vtec(_tt,ss);
			else if (_type == UPD_COR) tmp = _decode_upd(_tt, ss);
			else {
				if (_log) _log->comment(4, "bnc", "SSR NO " + corr2str(_type));
			}

			cnt++;
			consume += t_gcoder::_consume(tmpsize);

		} // while
		_mutex.unlock(); return consume;

	} //end method


	int t_bnccorr::_decode_orbit(const t_gtime& epo, istringstream & ss)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		//_mutex.lock();
		int satnum;
		char satsys;
		int IOD;
		t_gtriple dxyz(0.0, 0.0, 0.0);
		t_gtriple dvxyz(0.0, 0.0, 0.0);
		ss >> satsys >> satnum;
		if (!ss.fail() && (satsys == 'G' ||
			satsys == 'R' ||
			satsys == 'E' ||
			satsys == 'C' ||
			satsys == 'S' ||
			satsys == 'J')) {
			string prn = t_gsys::eval_sat(int2str(satnum), t_gsys::char2gsys(satsys));
			GSYS gsys = t_gsys::char2gsys(satsys);

			ss >> IOD >> dxyz[0] >> dxyz[1] >> dxyz[2] >> dvxyz[0] >> dvxyz[1] >> dvxyz[2];

			// fill single data record
			if (!_filter_gnss(prn)) {
				if (_log) _log->comment(4, "SSR ORBIT", "skip " + prn);
			}
			else {
				map<string, t_gdata*>::iterator it = _data.begin();
				while (it != _data.end()) {
					if (it->second->id_type() == t_gdata::ALLPREC)
					{
						((t_gallprec*)it->second)->add_delta_pos_vel(prn, epo, IOD, dxyz, dvxyz);
						((t_gallprec*)it->second)->set_end(epo);
					}

					it++;
				}
			}

		}
		//_mutex.unlock();
		return 0;
	}

	int t_bnccorr::_decode_clock(const t_gtime& epo, istringstream & ss)
	{
		int satnum;
		char satsys;
		int IOD;
		double C0 = UNDEFVAL_CLK, C1 = UNDEFVAL_CLK, C2 = UNDEFVAL_CLK;
		ss >> satsys >> satnum;
		if (!ss.fail() && (satsys == 'G' ||
			satsys == 'R' ||
			satsys == 'E' ||
			satsys == 'C' ||
			satsys == 'S' ||
			satsys == 'J')) {
			string prn = t_gsys::eval_sat(int2str(satnum), t_gsys::char2gsys(satsys));
			GSYS gsys = t_gsys::char2gsys(satsys);

			ss >> IOD >> C0 >> C1 >> C2;

			// fill single data record
			if (!_filter_gnss(prn)) {
				if (_log) _log->comment(4, "SSR CLK", "skip " + prn);
			}
			else {
				map<string, t_gdata*>::iterator it = _data.begin();
				while (it != _data.end()) {
					if (it->second->id_type() == t_gdata::ALLPREC)
					{
						((t_gallprec*)it->second)->add_delta_clk(prn, epo, IOD, C0, C1, C2);
						((t_gallprec*)it->second)->set_end(epo);
					}

					it++;
				}
			}

		}

		return 0;
	}

	int t_bnccorr::_decode_cb(const t_gtime& epo, istringstream & ss)
	{
		int satnum;
		char satsys;
		int number;
		string channel;
		double val;
		ss >> satsys >> satnum;
		if (!ss.fail() && (satsys == 'G' ||
			satsys == 'R' ||
			satsys == 'E' ||
			satsys == 'C' ||
			satsys == 'S' ||
			satsys == 'J')) {
			string prn = t_gsys::eval_sat(int2str(satnum), t_gsys::char2gsys(satsys));
			GSYS gsys = t_gsys::char2gsys(satsys);

			// fill single data record
			if (!_filter_gnss(prn)) {
				if (_log) _log->comment(4, "SSR Code Bias", "skip " + prn);
			}
			else {

				ss >> number;
				for (int i = 0; i < number; i++)
				{
					shared_ptr<t_gbias> code_bias = make_shared<t_gbias>();                            // create first reference bias 
					ss >> channel >> val;
					code_bias->set(val, str2gobs('C' + channel));

					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end()) {
						if (it->second->id_type() == t_gdata::ALLBIAS)
						{
							((t_gallbias*)it->second)->add(prn, epo, channel, code_bias);
						}

						it++;
					}
				}
			}

		}

		return 0;
	}

	int t_bnccorr::_decode_pb(const t_gtime& epo, istringstream & ss)
	{
		int satnum;
		char satsys;
		double yaw, dyaw;
		int number;
		string channel;
		double val;
		ss >> satsys >> satnum;
		if (!ss.fail() && (satsys == 'G' ||
			satsys == 'R' ||
			satsys == 'E' ||
			satsys == 'C' ||
			satsys == 'S' ||
			satsys == 'J')) {
			string prn = t_gsys::eval_sat(int2str(satnum), t_gsys::char2gsys(satsys));
			GSYS gsys = t_gsys::char2gsys(satsys);

			// fill single data record
			if (!_filter_gnss(prn)) {
				if (_log) _log->comment(4, "SSR Phase Bias", "skip " + prn);
			}
			else {

				ss >> yaw >> dyaw >> number;
				for (int i = 0; i < number; i++)
				{
					shared_ptr<t_gbias> phase_bias = make_shared<t_gbias>();                            // create first reference bias 
					ss >> channel >> val;
					phase_bias->set(val, str2gobs('L' + channel));

					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end()) {
						if (it->second->id_type() == t_gdata::ALLBIAS)
						{
							((t_gallbias*)it->second)->add(prn, epo, channel, phase_bias);
						}

						it++;
					}
				}
			}

		}

		return 0;
	}

	int t_bnccorr::_decode_vtec(const t_gtime& epo, istringstream & ss)
	{
		return 0;
	}

	int t_bnccorr::_decode_upd(const t_gtime& epo, istringstream& ss)
	{
		int satnum;
		char satsys;
		int valid;
		double wl_val, wl_sig, wl_nrec, nl_val, nl_sig, nl_nrec;
		t_updrec upd;
		ss >> satsys >> satnum;
		if (!ss.fail() && (satsys == 'G' ||
			satsys == 'R' ||
			satsys == 'E' ||
			satsys == 'C' ||
			satsys == 'S' ||
			satsys == 'J')) {
			string prn = t_gsys::eval_sat(int2str(satnum), t_gsys::char2gsys(satsys));
			GSYS gsys = t_gsys::char2gsys(satsys);

			// fill single data record
			if (!_filter_gnss(prn)) {
				if (_log) _log->comment(4, "SSR UPD Cors", "skip " + prn);
			}
			else {

				ss >> valid >> wl_val >> wl_sig >> wl_nrec >> nl_val >> nl_sig >> nl_nrec;

				if (valid == 1)
				{
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end()) {
						if (it->second->id_type() == t_gdata::UPD) {
							upd.obj = prn;
							// WL
							upd.sigma = wl_sig;
							upd.value = wl_val;
							upd.npoint = wl_nrec;
							((t_gupd*)it->second)->add_sat_upd(UPDTYPE::WL, epo, prn, upd);

							// NL
							upd.sigma = nl_sig;
							upd.value = nl_val;
							upd.npoint = nl_nrec;
							((t_gupd*)it->second)->add_sat_upd(UPDTYPE::NL, epo, prn, upd);
						}
						it++;
					}

				}
			}

		}

		return 0;
	}

	int t_bnccorr::_decode_aug(const string& site, const t_gtime& epo, istringstream& ss)
	{
		/*
		          sgpsw  sow          site  prn      p1       l1        p2        l2
		   > AUG  2103 345689.000000  HKST  G03    4.4449  -2.9724    6.8491   -5.3766
		   > AUG  2103 345689.000000  HKST  G16    8.5141  -7.3809   13.6521  -12.5189
		   > AUG  2103 345689.000000  HKST  E21    1.4067  -0.5427    2.1592   -1.2953
		*/
		int satnum;
		char satsys;
		double p1, l1, p2, l2;

		ss >> satsys >> satnum;
		if (!ss.fail() && (satsys == 'G' ||
			satsys == 'R' ||
			satsys == 'E' ||
			satsys == 'C' ||
			satsys == 'S' ||
			satsys == 'J')) {
			string prn = t_gsys::eval_sat(int2str(satnum), t_gsys::char2gsys(satsys));
			GSYS gsys = t_gsys::char2gsys(satsys);

			// fill single data record
			if (!_filter_gnss(prn)) {
				if (_log) _log->comment(4, "SSR UPD Cors", "skip " + prn);
			}
			else {

				ss >> p1 >> l1 >> p2 >> l2;

				t_augtype augtype;
				map<string, t_gdata*>::iterator it = _data.begin();
				while (it != _data.end()) {
					if (it->second->id_type() == t_gdata::AUG) {
						augtype = make_pair(AUGTYPE::TYPE_P, _band_index[gsys][FREQ_1]);
						((t_gaug*)it->second)->add_data(site, epo, prn, augtype, p1);

						augtype = make_pair(AUGTYPE::TYPE_P, _band_index[gsys][FREQ_2]);
						((t_gaug*)it->second)->add_data(site, epo, prn, augtype, p2);

						augtype = make_pair(AUGTYPE::TYPE_L, _band_index[gsys][FREQ_1]);
						((t_gaug*)it->second)->add_data(site, epo, prn, augtype, l1);

						augtype = make_pair(AUGTYPE::TYPE_L, _band_index[gsys][FREQ_2]);
						((t_gaug*)it->second)->add_data(site, epo, prn, augtype, l2);			
					}
					it++;
				}
			}

		}
		return 1;
	}

	bool t_bnccorr::available(const t_gtime& now)
	{

		if (now < _tt || now.diff(_tt) < _udorbInt)return true;
		else
		{
			//return _available;
			return false;
		}

	}

	void t_bnccorr::available_false()
	{
		gtrace("t_gtcp::available_false");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		//_mutex.lock();

		_available = false;

		//_mutex.unlock();
	}

	void t_bnccorr::available_true()
	{
		gtrace("t_gtcp::available_true");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		//_mutex.lock();

		_available = true;

		//_mutex.unlock();
	}

	string t_bnccorr::corr2str(const CORR_TYPE & type)
	{
		switch (type)
		{
		case ORBIT:return "ORBIT";
		case CLOCK:return "CLOCK";
		case CODE_BIAS:return "CODE_BIAS";
		case PHASE_BIAS:return"PHASE_BIAS";
		case VTEC:return "VTEC";
		case UPD_COR:return "UPD";
		case AUG_COR:return "AUG";
		default:
			break;
		}

		return "DEFAULT";
	}

	t_bnccorr::CORR_TYPE t_bnccorr::str2corr(const string & str)
	{
		string tmp = str;
		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
		if (tmp == "ORBIT") return ORBIT;
		else if (tmp == "CLOCK") return CLOCK;
		else if (tmp == "CODE_BIAS") return CODE_BIAS;
		else if (tmp == "PHASE_BIAS") return PHASE_BIAS;
		else if (tmp == "VTEC") return VTEC;
		else if (tmp == "UPD")  return UPD_COR;
		else if (tmp == "AUG")  return AUG_COR;
		else return DEFAULT;
	}

	bool t_bnccorr::_validepo(const t_gtime& t)
	{
		if (t > _end || t < _beg) return false;

		if (_filter_epoch(t)) return true;
		else return false;
	}

	// Checking if a receiver object exists
	// ----------------------------------
	void t_bnccorr::_check_recobj(string& site)
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
