

#include "gambfix/gambfixchk.h"
#include "gambfix/gallamb.h"
#include "gcoders/ambcon_leo.h"
#include "gset/gsetout.h"
#include "gio/gfile.h"

#include "math.h"
#include <string>
#include <stdio.h>
using namespace std;

namespace great
{
	t_gambchk::t_gambchk()
	{

	}

	t_gambchk::t_gambchk(t_glog* log, t_gsetbase* set)
	{
		cerr << "================> constructor of t_gambchk" << endl;

		_gset = set;
		_glog = log;
	}


	t_gambchk::~t_gambchk()
	{

	}

	bool t_gambchk::ProcessBatch()
	{
		string func_id = "ProcessBatch";
		if (!get_setting())
		{
			return false;
		}
		// 判断是否有数据
		if (!this->_glog || !this->_gset || !this->_grecover1 || !this->_grecover2 || !this->_gambcon)
		{
			if (this->_glog) this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, this->_class_id, func_id, "check the data");
			return false;
		}
		
		int iepo = 0, mepo = 0, beg_epo = 0, nepo = 0;
		double res1=0.0, res2=0.0, wgt=0.0;
		t_gtime epo = _beg;
		int epo_beg, epo_end;
		set<string>::iterator site_beg = _sites.begin();
		set<string>::iterator site_end = _sites.end();
		set<string>::iterator sat_beg = _sats.begin();
		set<string>::iterator sat_end = _sats.end();
		t_map_site_equ allrecover1 = _grecover1->get_map_site_equ();
		t_map_site_equ allrecover2 = _grecover2->get_map_site_equ();
		int isnew = 0;
		for (set<string>::iterator itsite = site_beg; itsite != site_end; itsite++)
		{
			for (set<string>::iterator itsat = sat_beg; itsat != sat_end; itsat++)
			{
				epo = _beg;
				mepo = 0;
				res1 = 0;
				res2 = 0;
				wgt = 0;
				beg_epo = 0;
				while (epo <= _end)
				{
					auto tmp1 = allrecover1[*itsite][*itsat][epo];
					if (tmp1.size() != 0) isnew = tmp1[0]->is_newamb;
					else isnew = 0;
					
					if (isnew == 1 || epo == _end)
					{
						if (mepo > 0)
						{
							res1 = sqrt(res1 / wgt);
							res2 = sqrt(res2 / wgt);
							if (res2 > 10 && res2 / res1 > 1.7)
							{
								epo_end = epo_beg + mepo - 1;
								//cout << *itsite << "   " << *itsat << "   " << epo_beg << "   " << epo_end << "   " << res1 << "   " << res2 << endl;
								for(int i=0;i<mepo;i++)
									_flag[*itsite][*itsat][epo_beg+i] = 1;
								

							}
							mepo = 0;
							res1 = 0;
							res2 = 0;
							wgt = 0;
						}
						epo_beg = (epo - _beg) / _intv + 1;
	
					}

					

					tmp1 = allrecover1[*itsite][*itsat][epo];
					auto tmp2 = allrecover2[*itsite][*itsat][epo];
					if (tmp1.size() == 0 || tmp2.size() == 0)
					{
						epo = epo + _intv;
						continue;
					}
					/*if (mepo == 0)
					{
						epo_beg = epo;
					}*/
					mepo++;
					auto iter1 = tmp1[0];
					res1 += pow(iter1->resuidal * 1000, 2)*iter1->weight;
					wgt += iter1->weight;
					auto iter2 = tmp2[0];
					res2 += pow(iter2->resuidal * 1000, 2)*iter2->weight;
					epo = epo + _intv;
				}
				
			}

		}

		
		write_ambcon();

		
	}

	bool t_gambchk::get_setting()
	{
		_beg = dynamic_cast<t_gsetgen*>(this->_gset)->beg();
		_end = dynamic_cast<t_gsetgen*>(this->_gset)->end();
		_bias = 100;
		_sites = dynamic_cast<t_gsetgen*>(this->_gset)->recs();
		_sats = dynamic_cast<t_gsetgnss*>(this->_gset)->sat();
		t_gtime beg_ambcon;
		set<string> sat_ambcon, site_ambcon, sats, sites;
		_gambcon->get_header(beg_ambcon, _intv, _ambnum, sat_ambcon, site_ambcon);


		set<string>::const_iterator site_beg = _sites.begin();
		set<string>::const_iterator site_end = _sites.end();
		set<string>::const_iterator site2_beg = site_ambcon.begin();
		set<string>::const_iterator site2_end = site_ambcon.end();
		set_intersection(site_beg, site_end, site2_beg, site2_end,
			inserter(sites, sites.begin()));
		_sites = sites;

		set<string>::const_iterator sat_beg = _sats.begin();
		set<string>::const_iterator sat_end = _sats.end();
		set<string>::const_iterator sat2_beg = sat_ambcon.begin();
		set<string>::const_iterator sat2_end = sat_ambcon.end();
		set_intersection(sat_beg, sat_end, sat2_beg, sat2_end,
			inserter(sats, sats.begin()));
		_sats = sats;

		if (_sites.size() == 0 || _sats.size() == 0) {
			if (_glog) {
				_glog->comment(0, "need site or sat!");
			}
			return false;
		}
		/*_res1 = _grecover1->get_map_time_equ;
		_res2 = _grecover2->get_map_time_equ;*/
	}

	bool t_gambchk::add_data(t_gallrecover* res1, t_gallrecover* res2, t_gambcon* ambcon)
	{
		if (!res1 && !res2 && !ambcon) return false;
		// 获取交集
		try
		{
			_grecover1 = res1;
			_grecover2 = res2;
			_gambcon = ambcon;
		}
		catch (...)
		{
			return false;
		}

		return true;

	}

	void t_gambchk::write_ambcon()
	{
		

		t_gambcon     _gambcon_save;
		t_gtime epo = _beg;
		


		while (epo <= _end)
		{
			auto sdambs = _gambcon->get_Sds(epo);
			if (sdambs.empty())
			{
				epo = epo + _intv;
				continue;
			}
			for (const auto& iter : sdambs)
			{
				t_gamb_Sd Sd_amb = iter;
				const t_gtime& beg_time = Sd_amb.beg;
				const t_gtime& end_time = Sd_amb.end;
				string rec = Sd_amb.sgl_obj;
				string sat_1st = Sd_amb.dbl_obj_1st;
				string sat_2nd = Sd_amb.dbl_obj_2nd;
				int beg = (beg_time - _beg) / _intv + 1;
				int end = (end_time - _beg) / _intv + 1;
				if (_flag[rec][sat_1st][beg] == 1 || _flag[rec][sat_2nd][beg] == 1)
				{
					//cout << sat_1st << "  " << sat_2nd <<"  "<<beg<< endl;
					_ambnum--;
					continue;				
					
				}
					
				_gambcon_save.add_sdamb(Sd_amb);

			}
			epo = epo + _intv;
		}
		_gambcon_save.set_header(_beg, _intv, _ambnum, _sats, _sites);
		t_gfile gout;
		t_ambcon_leo gamb_coder(_gset, "", 4096);
		string path = (dynamic_cast<t_gsetout*>(this->_gset)->outputs("ambcon_leo"));
		if (path.empty())
		{
			path = "grt$(date).con";
		}
		substitute(path, "$(date)", _beg.str_yyyydoy(), false);
		gout.path(path); 
		gout.glog(this->_glog);
		gamb_coder.clear(); 
		gamb_coder.glog(this->_glog); 
		gamb_coder.path(path);
		gamb_coder.add_data("ID1", &_gambcon_save);
		gout.coder(&gamb_coder);
		gout.run_write();
	}
}



