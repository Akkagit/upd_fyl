/**
*
* @verbatim
	History
	 -1.0 liyuhao		2020-07-27  creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gclkdif.cpp
* @brief		clock differ processing for GNSS satellite
* @author       liyuhao, Wuhan University
* @version		1.0.0
* @date		    2020-07-27
*
*/

#include <algorithm>
#include <set>
#include "gclkdif.h"

#include "gset/gsetout.h" 
#include "gset/gsetorbdif.h"
using namespace std;

namespace great
{

	t_gclkdif::t_gclkdif()
	{
	}

	t_gclkdif::t_gclkdif(t_gsetbase* set, t_glog* log)
	{
		_set = set;
		_glog = log;
	}

	t_gclkdif::~t_gclkdif()
	{
	}

	bool t_gclkdif::processBatch()
	{
		string function_name = "t_gclkdif::processBatch";

		if (!_gclkprd || !_gclkref)
		{
			_glog->comment(0, "clkdif", "ERROR : clkdif ==> t_gclkdif::no comeplete file !");
			return false;
		}
		// get begin/end time,sampling, sat list from xml
		if (!_get_setting())
		{
			_glog->comment(0, "clkdif", "ERROR : clkdif ==> t_gclkdif::there is something wrong witn xml !");
			return false;
		}
		_set_out();

		//check time
		if (!_check_time(""))
		{
			return false;
		}

		////get and skip satellites
		if (!_get_sat())
		{
			return false;
		}

		_init_par();

		t_gtime epo_now = _beg;
		int nepoch = 0;
		double clk_rms = 0.0;
		double dclk = 0.0;

		// loop by epoch
		while (epo_now < _end)
		{
			int ncount = 0;
			// loop by sat
			for (auto iter = _sat.begin(); iter != _sat.end();iter++)
			{
				//get clock in rinexc reference file and rinexc product file
				double clk_tmp[2], clkdif;
				if (dynamic_cast<t_gallprec*>(_gclkprd)->clk(*iter, epo_now, &clk_tmp[0], &clk_rms, &dclk) == -1 || dynamic_cast<t_gallprec*>(_gclkref)->clk(*iter, epo_now, &clk_tmp[1], &clk_rms, &dclk) == -1)
				{
					clk_tmp[0] = 0.0;
					clk_tmp[1] = 0.0;
					_flag[nepoch][ncount] = 0;
				}
				clkdif = clk_tmp[0] - clk_tmp[1];
				clkdif = clkdif * 1e9;
				_res[nepoch][ncount]=clkdif;
				ncount++;
			}
			nepoch++;
			epo_now = epo_now + _sampling;
		}
		
		//clock difference between satellites
		for (int ncount = 0; ncount < _numsat; ncount++)
		{
			for (int nepoch = 0; nepoch < _nepoch; nepoch++)
			{
				if (ncount ==_ncount)
				{
					break;
				}
				_res[nepoch][ncount] = _res[nepoch][ncount] - _res[nepoch][_ncount];
			}
		}
		//rms statistics
		_rms_statistics();

		//generate result
		_generate_result();
	}

	void t_gclkdif::setCLKPRD(t_gallnav* gclkprd)
	{
		_gclkprd = gclkprd;
	}

	void t_gclkdif::setCLKREF(t_gallnav* gclkref)
	{
		_gclkref = gclkref;
	}

	bool t_gclkdif::_get_setting()
	{
		_beg = dynamic_cast<t_gsetgen*>(_set)->beg();
		_end = dynamic_cast<t_gsetgen*>(_set)->end();
		_sampling = dynamic_cast<t_gsetgen*>(_set)->sampling();

		_end = _beg + ((int)((_end - _beg) / _sampling)) * _sampling;

		_refsat = dynamic_cast<t_gsetgen*>(_set)->refsat();

		// get satlist from xml
		set<string> sys = dynamic_cast<t_gsetgen*>(_set)->sys();
		if (sys.empty())
		{
			_glog->comment(0, "clkdif", "ERROR : clkdif ==> No satellite systerm, check your xml file(sys) !");
			return false;
		}
		// sats list from gset
		_sat = dynamic_cast<t_gsetgnss*>(_set)->sat();
		if (_sat.empty())
		{
			t_map_sats gnss_sats = GNSS_SATS();
			for (auto itGNS = sys.begin(); itGNS != sys.end(); ++itGNS) {
				GSYS gsys = t_gsys::str2gsys(*itGNS);
				_sat.insert(gnss_sats[gsys].begin(), gnss_sats[gsys].end());
			}
		}
		_numsat = _sat.size();
		return true;
	}

	bool t_gclkdif::_set_out()
	{
		string tmp = dynamic_cast<t_gsetout*>(_set)->outputs("clkdif");

		if (tmp.empty())
		{
			string tmp = dynamic_cast<t_gsetout*>(_set)->outputs("dif");
		}

		if (tmp.empty())
		{
			tmp = "clkdif" + int2str(_beg.year()) + int2str(_beg.doy(), 3);
		}
		else
		{
			substitute(tmp, "$(date)", int2str(_beg.year()) + int2str(_beg.doy(), 3), false);
		}
		_clkdif = new t_giof;
		_clkdif->mask(tmp);
		_clkdif->append(false);

		return true;
	}

	bool t_gclkdif::_check_time(string sat)
	{
		t_gtime beg_prd = dynamic_cast<t_gallprec*>(_gclkprd)->beg_clk(sat);
		t_gtime end_prd = dynamic_cast<t_gallprec*>(_gclkprd)->end_clk(sat);

		t_gtime beg_ref = dynamic_cast<t_gallprec*>(_gclkref)->beg_clk(sat);
		t_gtime end_ref = dynamic_cast<t_gallprec*>(_gclkref)->end_clk(sat);

		if (_beg < beg_prd) _beg = beg_prd;
		if (_beg < beg_ref) _beg = beg_ref;
		if (_end > end_prd) _end = end_prd;
		if (_end > end_ref) _end = end_ref;
		if (_beg >= _end)
		{
			_glog->comment(0, "clkdif", "ERROR : clkdif ==> time is otrange, please chelk time in xml file or rinexc reference file or rinexc product file  !");

			return false;
		}
		return true;
	}

	bool t_gclkdif::_get_sat()
	{
		string function_name = "t_gclkdif::_get_sat";
		set<string> sat_prd = dynamic_cast<t_gallprec*>(_gclkprd)->clk_objs();
		set<string> sat_ref = dynamic_cast<t_gallprec*>(_gclkref)->clk_objs();
		for (auto iter = _sat.begin(); iter != _sat.end();)
		{
			//skip sat
			if (sat_prd.find(*iter) == sat_prd.end())
			{
				_glog->comment(1, function_name, "WARNING : " + *iter + "can not be found in rinexc product file!");
				_sat.erase(iter++);
				continue;
			}
			if (sat_ref.find(*iter) == sat_ref.end())
			{
				_glog->comment(1, function_name, "WARNING : " + *iter + "can not be found in rinexc reference file!");
				_sat.erase(iter++);
				continue;
			}
			iter++;
		}
		return true;
	}

	bool t_gclkdif::_init_par()
	{
		string function_name = "t_gclkdif::_init_par";
		_nepoch = (int)((_end - _beg) / _sampling);
		_res.resize(_nepoch);
		_numsat = _sat.size();
		for (int i = 0; i < _nepoch; i++)
		{
			_res[i].resize(_numsat);
		}

		_flag.resize(_nepoch);
		for (int tmp = 0; tmp < _nepoch; tmp++)
		{
			_flag[tmp].resize(_numsat);
			for (int ncount = 0; ncount < _numsat; ncount++)
			{
				_flag[tmp][ncount] = 1;
			}
		}

		set<string>::iterator refsat;
		if ((refsat = _sat.find(_refsat)) == _sat.end())
		{
			_glog->comment(1, function_name, "ERROR : clkdif ==> t_gclkdif::there is no reference satellite in satellite list!");
			return false;
		}
		int ncount = 0;
		for (auto iter = _sat.begin(); iter != _sat.end(); iter++)
		{
			if (*iter == _refsat)
			{
				_ncount = ncount;
				break;
			}
			ncount++;
		}
		return true;
	}

	bool t_gclkdif::_rms_statistics()
	{
		//3 for MEAN,RMS,STD
		_res_rms.resize(3);
		for (int tmp = 0; tmp < 3; tmp++)
		{
			_res_rms[tmp].resize(_numsat);
		}
		double mean_rms = 0.0;
		double mean_std = 0.0;
		//1 for DEL
		_res_del.resize(_numsat);

		for (int ncount = 0; ncount < _numsat;)
		{
			double mean = 0.0;
			double rms = 0.0;
			double std = 0.0;
			int m = 0;
			int ndel = 0;

			for (int nepoch = 0; nepoch < _nepoch; nepoch++)
			{
				if (_flag[nepoch][ncount])
				{
					if (nepoch >= 1)
					{
						if (_res[nepoch][ncount] - _res[nepoch - 1][ncount] > 100.0)
						{
							_flag[nepoch][ncount] = 0;
							continue;
						}
					}
					mean = mean + _res[nepoch][ncount];
					rms = rms + _res[nepoch][ncount] * _res[nepoch][ncount];
					m++;
				}
			}
			mean = mean / m;

			for (int nepoch = 0; nepoch < _nepoch; nepoch++)
			{
				if (_flag[nepoch][ncount])
				{
					std = std + (_res[nepoch][ncount] - mean) * (_res[nepoch][ncount] - mean);
				}
			}

			if (m == 1)
			{
				rms = 999;
				std = 999;
			}
			else
			{
				rms = sqrt(rms / m);
				std = sqrt(std / (m - 1));
			}

			for (int nepoch = 0; nepoch < _nepoch; nepoch++)
			{
				if (!_flag[nepoch][ncount] || fabs(_res[nepoch][ncount] - mean) <= 3 * std)
					continue;
				_flag[nepoch][ncount] = 0;
				ndel++;
			}
			_res_del[ncount] = _res_del[ncount] + ndel;

			if (!ndel)
			{
				_res_rms[0][ncount] = mean;
				_res_rms[1][ncount] = rms;
				_res_rms[2][ncount] = std;
				if (ncount != _ncount)
				{
					mean_rms = mean_rms + _res_rms[1][ncount];
					mean_std = mean_std + _res_rms[2][ncount];
				}
				ncount++;
			}
		}
		mean_rms = mean_rms / (_numsat - 1);
		mean_std = mean_std / (_numsat - 1);
		_mean_rms = mean_rms;
		_mean_std = mean_std;

		return true;
	}

	bool t_gclkdif::_generate_result()
	{
		ostringstream os;
		os << setw(4) << right << "NAME";
		for (auto iter = _sat.begin(); iter != _sat.end();iter++)
		{
			os << setw(9) << right << *iter;
		}
		os << endl;
		for (int nepoch = 0; nepoch < _nepoch; nepoch++)
		{
			os << setw(4) << right << nepoch + 1;
			for (int i = 0; i < _sat.size(); i++)
			{
				os << setw(9) << fixed << setprecision(2) << right << _res[nepoch][i];
			}
			os << endl;
		}
		os << "-------------------------------------------" << endl;
		os << setw(4) << right << "MEAN";
		for (int ncount = 0; ncount < _numsat; ncount++)
		{
			os << setw(9) << fixed << setprecision(2) << right << _res_rms[0][ncount];
		}
		os << endl;
		os << setw(4) << right << "RMS:" ;
		for (int ncount = 0; ncount < _numsat; ncount++)
		{
			os << setw(9) << fixed << setprecision(2) << right << _res_rms[1][ncount];
		}
		os << endl;
		os << setw(4) << right << "STD:";
		for (int ncount = 0; ncount < _numsat; ncount++)
		{
			os << setw(9) << fixed << setprecision(2) << right << _res_rms[2][ncount];
		}
		os << endl;
		os << setw(4) << right << "DEL:";
		for (int ncount = 0; ncount < _numsat; ncount++)
		{
			os << setw(9) << right << _res_del[ncount];
		}
		os << endl;
		os << setw(4) << right << "NAME";
		for (auto iter = _sat.begin(); iter != _sat.end(); iter++)
		{
			os << setw(9) << right << *iter;
		}
		os << endl;
		os << setw(3) << right << "EOF" << endl;
		os << " MEAN RMS =";
		os << setw(7) << fixed << setprecision(2) << right << _mean_rms;
		os << "  MEAN STD";
		os << setw(7) << fixed << setprecision(2) << right << _mean_std;
		os << endl;

		_clkdif->write(os.str().c_str(), os.str().size());
		_clkdif->flush();
		if (_clkdif)
		{
			delete _clkdif;
			_clkdif = NULL;
		}

		return true;
	}
}