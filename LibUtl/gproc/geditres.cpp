/**
*
* @verbatim
	History
	 -1.0 Hongjie Zheng		2019-09-06  creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    geditres.cpp
* @brief		Edit residuals function
* @author       Hongjie Zheng, Wuhan University
* @version		1.0.0
* @date		    2019-09-06
*
*/

#include <gproc/geditres.h>
#include <gcoders/ambflag.h>
#include <gcoders/recover.h>
#include <gset/gsetinp.h>
#include <gset/gcfg_editres.h>
#include <gset/gsetinp.h>
#include <gset/gsetout.h>
#include <gio/gfile.h>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include "gutils/gstring.h"
#include "gall/gallsum.h"
#include "gutils/gleo.h"


#define FREQ_L12 "L12"
#define FREQ_L13 "L13"
#define FREQ_L14 "L14"
#define FREQ_L15 "L15"

#ifdef DEBUG
void print_flag_debuginfo(string site, string sat, vector<int>& flag_data, int number = 0) {
	//cout << TEST_INFO << " initial flag value " << endl;
	cout << TEST << setw(5) << site << setw(5) << sat;
	if (number == 0) {
		for (int value : flag_data) {
			cout << setw(10) << value;
		}
	}
	else {
		for (int i = 0; i < number; i++) {
			cout << setw(10) << flag_data[i];
		}
	}
	cout << endl;
}
#endif // DEBUG


namespace great
{

	const  int t_geditres::GOOD = 0;
	const  int t_geditres::OLDAMB = 1;
	const  int t_geditres::NEWAMB = 2;
	const  int t_geditres::LOWELE = 3;
	const  int t_geditres::NOTALL = 4;
	const  int t_geditres::DELSHRT = 5;
	const  int t_geditres::DELBAD = 6;
	const  int t_geditres::DELRMS = 7;
	const  int t_geditres::NEWBAD = 8;
	const  int t_geditres::NODATA = 9;



	t_geditres::t_geditres(t_glog * log) :
		_log(log),
		beg_time(t_gtime()),
		end_time(t_gtime()),
		interval(0),
		bad_elisp(50),
		jump_elisp(50),
		short_elisp(600),
		res_max(0.2),
		rec_num(10),
		_amb_ptrflag(false),
		_recover_ptrflag(false)
	{
	}

	t_geditres::t_geditres(t_gsetbase* set, t_glog * log) :
		_log(log),
		beg_time(t_gtime()),
		end_time(t_gtime()),
		interval(0),
		freq("LC"),
		bad_elisp(50),
		jump_elisp(50),
		short_elisp(600),
		res_max(0.2),
		rec_num(10),
		_amb_ptrflag(false),
		_recover_ptrflag(false)
	{

		t_gsetres* temp_setres = dynamic_cast<t_gsetres*>(set);
		if (temp_setres) 
		{
			bad_elisp = temp_setres->get_bad_elisp();
			jump_elisp = temp_setres->get_jump_elisp();
			short_elisp = temp_setres->get_short_elisp();
			freq = temp_setres->get_freq();
			mode = temp_setres->get_freq_mode();
			res_max = temp_setres->get_res_max();
			rec_num = temp_setres->get_sta_num();
		}

		_gset = set;
		_all_sum = new t_gallsum();
		t_gsetinp* temp_input = dynamic_cast<t_gsetinp*>(set);
		if (temp_input)
		{
			multimap<IFMT, string> inp = temp_input->inputs_all();
			multimap<IFMT, string>::const_iterator itINP = inp.begin();
			for (size_t i = 0; i < inp.size() && itINP != inp.end(); ++i, ++itINP)
			{
				if( 
					(itINP->first == AMBFLAG12_INP && mode == FREQ_L12) ||
					(itINP->first == AMBFLAG13_INP && mode == FREQ_L13) ||
					(itINP->first == AMBFLAG14_INP && mode == FREQ_L14) ||
					(itINP->first == AMBFLAG15_INP && mode == FREQ_L15)	
				   )
				{
					_logfiles.push_back(itINP->second);
				}
				else if (itINP->first == RECOVER_INP)
				{
					this->set_recover_data_by_great_resfile(itINP->second);
				}
				else if (itINP->first == PREOBS_INP) 
				{
					set_preOBS(itINP->second);
				}
			}
			this->set_ambflag_data(_logfiles);
		}
	}

	t_geditres::~t_geditres()
	{
		if (_amb_ptrflag && _all_ambflag) {
			delete _all_ambflag;
			_all_ambflag = nullptr;
		}

		if (_recover_ptrflag && _all_recover) {
			delete _all_recover;
			_all_recover = nullptr;
		}


	}

	void t_geditres::set_recover_data(t_gallrecover * recover_data)
	{
		if (_recover_ptrflag && _all_recover) {
			delete _all_recover;
			_all_recover = nullptr;
		}
		_recover_ptrflag = false;
		_all_recover = recover_data;
	}

	void t_geditres::set_preOBS(string filename)
	{
		_all_preproc = new t_gallpreproc();

		_log->comment(t_glog::LOG_LV::LOG_INFO, "t_geditres", "set_preOBS", "Reading slrfile: " + filename);
	
		t_gio* gin = new t_gfile;
		gin->path(filename); gin->glog(_log);
		t_gcoder* read_coder = new t_preproc(nullptr, "", 4096);
		read_coder->clear(); read_coder->glog(_log); read_coder->path(filename);
		read_coder->add_data("ID1", _all_preproc);
		gin->coder(read_coder);
		gin->run_read();
		delete gin;
		delete read_coder;

		_log->comment(t_glog::LOG_LV::LOG_INFO, "t_geditres", "set_preOBS", "Reading slrfile: Finish!");
	}
	void t_geditres::set_recover_data_by_great_resfile(string resfile_name)
	{
		if (_recover_ptrflag && _all_recover) {
			delete _all_recover;
			_all_recover = nullptr;
		}
		_recover_ptrflag = true;
		_all_recover = new t_gallrecover();

		_log->comment(t_glog::LOG_LV::LOG_INFO, "t_geditres", "set_recover_data_by_resfile", "Reading resfile: " + resfile_name);

		t_gio* gin = new t_gfile;
		gin->path(resfile_name); gin->glog(_log);
		t_gcoder* read_coder = new t_resfile(nullptr, "", 4096);
		read_coder->clear(); read_coder->glog(_log); read_coder->path(resfile_name);
		read_coder->add_data("ID1", _all_recover);
		gin->coder(read_coder);
		gin->run_read();
		delete gin;
		delete read_coder;

		_log->comment(t_glog::LOG_LV::LOG_INFO, "t_geditres", "set_recover_data_by_resfile", "Reading resfile: Finish!");
	}

	void t_geditres::set_recover_data_by_panda_resfile(string resfile_name)
	{
		if (_recover_ptrflag && _all_recover) {
			delete _all_recover;
			_all_recover = nullptr;
		}
		_recover_ptrflag = true;
		_all_recover = new t_gallrecover();

		fstream resfile;

		resfile.open(resfile_name, ios::in);

		string line;
		vector<string> sat_list;
		vector<string> site_list;
		t_gtime beg_time;
		double interval = 30.0;
		while (getline(resfile, line)) {
			int idx = 0;
			if (idx = line.find("Time&Interval :") != -1) {
				stringstream temp_line(line.substr((int64_t)idx + 15));
				int mjd;
				double sod;
				temp_line >> mjd >> sod >> interval;
				_all_recover->set_interval(interval);
				beg_time = t_gtime(mjd, int(sod), sod - int(sod));
			}
			else if (idx = line.find("SAT:") != -1) {
				stringstream temp_line(line.substr((int64_t)idx + 4));
				string sat;
				while (true) {
					temp_line >> sat;
					if (sat == "") {
						break;
					}
					sat_list.push_back(sat);
					sat = "";
				}
			}
			else if (idx = line.find("STA:") != -1) {
				stringstream templine(line.substr((int64_t)idx + 4));
				string site;
				while (true) {
					templine >> site;
					if (site == "") {
						break;
					}
					site_list.push_back(site);
					site = "";
				}
			}
			else if (idx = line.find("RES") != -1) {
				int find_d;
				while ((find_d = line.find("D")) != string::npos) {
					line.replace(find_d, 1, "E");
				}
				stringstream templine(line.substr((int64_t)idx + 3));
				int epo, site_id, sat_id, is_newamb;
				double azim, elev;
				double Pweight, Pres, Lweight, Lres;
				templine >> epo >> site_id >> sat_id >> is_newamb >> azim >> elev
					>> Lres >> Lweight >> Pres >> Pweight;

				t_grecover_equation temp_data(beg_time + ((int64_t)epo - 1) * interval, site_list[(int64_t)site_id - 1], sat_list[(int64_t)sat_id - 1]);
				temp_data.set_recover_equation(t_gobscombtype("LC"), make_pair(Lweight, Lres), is_newamb);
				_all_recover->add_recover_equation(temp_data);
				//add for test
				temp_data.set_recover_equation(t_gobscombtype("PC"), make_pair(Pweight, Pres), 9);
				_all_recover->add_recover_equation(temp_data);
			}
		}

		resfile.close();
	}


	void t_geditres::set_ambflag_data(t_gallambflag * ambflag_data)
	{
		if (_amb_ptrflag && _all_ambflag) {
			delete _all_ambflag;
			_all_ambflag = nullptr;
		}
		_amb_ptrflag = false;
		_all_ambflag = ambflag_data;
	}
	void t_geditres::set_ambflag_data(const vector<string>& logfiles)
	{
		if (_amb_ptrflag && _all_ambflag) 
		{
			delete _all_ambflag;
			_all_ambflag = nullptr;
		}

		_recover_ptrflag = true;
		_all_ambflag = new t_gallambflag(t_gdata::AMBFLAG);
		_logfiles = logfiles;

		for (string logfile : _logfiles)
		{
			_log->comment(t_glog::LOG_LV::LOG_INFO, "t_geditres", "set_ambflag_data_by_logfiles", "Reading logfile: " + logfile);
			t_gio* gin = new t_gfile;
			gin->path(logfile);
			gin->glog(_log);
			t_gcoder* gcoder = new t_ambflag(nullptr, "");
			gcoder->clear(); gcoder->path(logfile);
			gcoder->glog(_log);
			gcoder->add_data("ID1", _all_ambflag);
			gin->coder(gcoder);
			gin->run_read();

			delete gin;
			delete gcoder;

			_log->comment(t_glog::LOG_LV::LOG_INFO, "t_geditres", "set_ambflag_data_by_logfiles", "Reading logfile: Finish!");
		}
	}

	void t_geditres::set_logfiles(const vector<string>& logfiles)
	{
		_logfiles = logfiles;
	}

	t_gallambflag t_geditres::get_gambflag_data() const
	{
		return *_all_ambflag;
	}

	void t_geditres::processBatch()
	{

		_log->comment(t_glog::LOG_LV::LOG_INFO, "t_geditres", "processBatch", "Starting editres process....");

		if (!_all_ambflag) {
			_log->comment(t_glog::LOG_LV::LOG_ERROR, "t_geditres", "processBatch", "No Ambflag Data!!!!");
			return;
		}
		if (!_all_recover) {
			_log->comment(t_glog::LOG_LV::LOG_ERROR, "t_geditres", "processBatch", "No Recover Data!!!!");
			return;
		}

		if ((freq == "L2" && mode == FREQ_L13) || (freq == "L3" && mode == FREQ_L12))
		{
			_log->logError("t_geditres", "processBatch", "freq and mode is conflict !!!!");
			return;
		}

		//!! change from seconds to epochs
		short_elisp = round(short_elisp / _all_recover->get_interval());
		if (short_elisp < 1) 
		{
			short_elisp = 2;
		}

		t_gsetres* temp_setres = dynamic_cast<t_gsetres*>(_gset);
		string mode = temp_setres->get_mode();
		if (mode == "SLR")
		{
			// edit residuals for slr
			_log->comment(t_glog::LOG_LV::LOG_INFO, "t_geditres", "processBatch", " edit residuals of slr ....");
			_slr_res_edit();
		}
		if (mode == "SLR" || mode == "LEO")
		{
			// write sum file ... 
			_log->comment(t_glog::LOG_LV::LOG_INFO, "t_geditres", "processBatch", " Write sum file ....");
			_rms_table();
		}
		_log->comment(t_glog::LOG_LV::LOG_INFO, "t_geditres", "processBatch", " Checking residuals ....");
		_check_resid();

		_log->comment(t_glog::LOG_LV::LOG_INFO, "t_geditres", "processBatch", " Update ambflag data....");
		_update_ambflagdata();

		_log->comment(t_glog::LOG_LV::LOG_INFO, "t_geditres", "processBatch", "Editres: Normal End!");
	}

	void t_geditres::_slr_res_edit()
	{
		set<string> sat_list_t;
		set<string> site_list_t;
		vector<t_gtime> time_list;
		set<t_gtime> time_list_t;
		vector<t_gtime> time_list_slr;
		t_map_site_equ site_equ;
		map<int, int> _num_slr;
		//map<int, int> _nsatepo_tmp;

		//_nsatepo = 0;

		sat_list_t = _all_recover->get_sat_list();
		site_list_t = _all_recover->get_site_list();
		time_list = _all_recover->get_all_obs_time_list();
		site_equ = _all_recover->get_map_site_equ();
		time_list_t = _all_recover->get_time_list();

		vector<string> sat_list;
		vector<string> site_list;
		vector<string> sat_list_v;
		vector<string> site_list_v;
		vector<string> sat_list_slr;
		vector<string> site_list_slr;

		set<string>::iterator it_sat;
		int i = 0;
		for (it_sat = sat_list_t.begin(); it_sat != sat_list_t.end(); it_sat++)
		{
			if (IS_SUPPORTED_LEO(strs2strl(*it_sat)))
			{
				sat_list_slr.push_back(*it_sat);
			}
			sat_list_v.push_back(*it_sat);
		}

		set<string>::iterator it_site;
		i = 0;
		for (it_site = site_list_t.begin(); it_site != site_list_t.end(); it_site++)
		{
			if (!IS_SUPPORTED_LEO(strs2strl(*it_site)))
			{
				site_list_slr.push_back(*it_site);
			}
			
			site_list_v.push_back(*it_site);
		}

		set<t_gtime>::iterator it_time;
		for (it_time = time_list_t.begin(); it_time != time_list_t.end(); it_time++)
		{
			time_list_slr.push_back(*it_time);
		}

		
		for (int ista = 0; ista <= site_list_slr.size(); ista++)
		{
			//_rmstbl_site_sat[ista][isat] = 0.0;
			_num_slr[ista] = 0;
		}
		
		map<string,map<t_gtime, t_gpreprocslr>>* slrobs;
		slrobs = _all_preproc->getallSLROBS();

		if (sat_list_slr.size() != 0)
		{
			for (size_t isat = 1; isat <= sat_list_slr.size(); isat++)
			{

				for (size_t ista_ = 1; ista_ <= site_list_slr.size(); ista_++)
				{

					for (size_t iepo = 0; iepo < time_list_slr.size(); iepo++)
					{
						vector<t_grecover_equation*> recover_equ = site_equ[site_list_slr[ista_ - 1]][sat_list_slr[isat - 1]][time_list_slr[iepo]];
						if (recover_equ.size() != 0)
						{
							//逐步收敛，后面改0.15
							if (abs(recover_equ[0]->resuidal) > res_max)
							{
								for (auto iter = (*slrobs)[strs2strl(sat_list_slr[isat - 1])].begin(); iter != (*slrobs)[strs2strl(sat_list_slr[isat - 1])].end(); iter++)
								{
									if (recover_equ[0]->time.mjd() == iter->first.mjd()&& recover_equ[0]->time.sod() == iter->first.sod())
									{
										(*slrobs)[strs2strl(sat_list_slr[isat - 1])].erase(iter);
										break;
									}
								}
							}
						}

					}
					int abc = 1;
				}
			}


		}


		//delete the slr station which only has less than two observations
		for (auto iterall = slrobs->begin(); iterall != slrobs->end(); iterall++)
		{
			for (auto iter = iterall->second.begin(); iter != iterall->second.end(); iter++)
			{
				for (int inum = 0; inum < site_list_slr.size(); inum++)
				{
					if (iter->second._sta_name == site_list_slr[inum])
					{
						_num_slr[inum] += 1;
					}
				}
			}
		}
		
		for (int ista = 0; ista < site_list_slr.size(); ista++)
		{
			//10 wait for testing
			if (_num_slr[ista] <= rec_num)
			{
				for (auto iterall = slrobs->begin(); iterall != slrobs->end(); iterall++)
				{
					for (auto iter = iterall->second.begin(); iter != iterall->second.end();)
					{
						if (iter->second._sta_name == site_list_slr[ista])
						{
							(*slrobs)[iterall->first].erase(iter++);
							continue;
						}
						++iter;
					}
				}
			}
		}

		for (auto iterall = slrobs->begin(); iterall != slrobs->end(); iterall++)
		{
			for (auto iter = iterall->second.begin(); iter != iterall->second.end(); iter++)
			{
				t_gtime time = iter->first;
				_all_preproc->setSLROBS(iterall->first, time, iter->second);
			}
		}
		_generate_pre_obs();

	}
	void t_geditres::_generate_pre_obs()
	{
		string tmp = dynamic_cast<t_gsetout*>(_gset)->outputs("preobs");
		if (tmp.empty())
		{
			tmp = "preOBS";
		}

		t_gio* gout_pre = new t_gfile;  t_gcoder* gpre_coder = new t_preproc(_gset, "", 4096);
		string pre_path = tmp;
		//string ics_path = (dynamic_cast<t_gsetinp*>(_set)->inputs("ics"))[0];
		gout_pre->path(pre_path); gout_pre->glog(_log);
		gpre_coder->clear(); gpre_coder->glog(_log); gpre_coder->path(pre_path);
		gpre_coder->add_data("ID1", _all_preproc);
		gout_pre->coder(gpre_coder);
		gout_pre->run_write();
		delete gout_pre; gout_pre = nullptr;
		delete gpre_coder; gpre_coder = nullptr;
	}

	void t_geditres::_rms_table()
	{
		
		vector<int> flag;
		set<string> sat_list_t;
		set<string> site_list_t;
		vector<t_gtime> time_list;
		set<t_gtime> time_list_t;
		vector<t_gtime> time_list_slr;
		t_map_site_equ site_equ;
		//map<int, int> _nsatepo_tmp;

		//_nsatepo = 0;

		sat_list_t=_all_recover->get_sat_list();
		site_list_t = _all_recover->get_site_list();
		time_list = _all_recover->get_all_obs_time_list();
		site_equ = _all_recover->get_map_site_equ();
		time_list_t = _all_recover->get_time_list();

		vector<string> sat_list;
		vector<string> site_list;
		vector<string> sat_list_v;
		vector<string> site_list_v;
		vector<string> sat_list_slr;
		vector<string> site_list_slr;
		
		set<string>::iterator it_sat;
		int i = 0;
		for (it_sat = sat_list_t.begin(); it_sat != sat_list_t.end(); it_sat++)
		{
			if (!IS_SUPPORTED_LEO(strs2strl(*it_sat)))
			{
				sat_list.push_back(*it_sat);
			}
			else
			{
				sat_list_slr.push_back(*it_sat);
			}
			sat_list_v.push_back(*it_sat);
		}

		set<string>::iterator it_site;
		i = 0;
		for (it_site = site_list_t.begin(); it_site != site_list_t.end(); it_site++)
		{
			if (IS_SUPPORTED_LEO(strs2strl(*it_site)))
			{
				site_list.push_back(*it_site);
			}
			else
			{
				site_list_slr.push_back(*it_site);
			}
			site_list_v.push_back(*it_site);
		}

		set<t_gtime>::iterator it_time;
		for (it_time = time_list_t.begin(); it_time != time_list_t.end(); it_time++)
		{
			time_list_slr.push_back(*it_time);
		}

		for (int isat = 0; isat <= sat_list.size(); isat++)
		{
			for (int ista = 0; ista <= site_list.size(); ista++)
			{
				//_rmstbl_site_sat[ista][isat] = 0.0;
				_nepotbl[ista][isat] = 0;
				_rmstable[t_gobscombtype("PC")][ista][isat] = 0.0;
				_rmstable[t_gobscombtype("LC")][ista][isat] = 0.0;
			}
		}

		for (int isat = 0; isat <= sat_list_slr.size(); isat++)
		{
			for (int ista = 0; ista <= site_list_slr.size(); ista++)
			{
				//_rmstbl_site_sat[ista][isat] = 0.0;
				_nepotbl_slr[ista][isat] = 0;
				_rmstable_slr["SLR"][ista][isat] = 0.0;
			}
		}
		
		int bcd = 0;
		int mn = 0;
		//add for slr
		if (sat_list_slr.size() != 0)
		{
			for (size_t isat = 1; isat <= sat_list_slr.size(); isat++)
			{
				
				for (size_t ista = 1; ista <= site_list_slr.size(); ista++)
				{
					bcd = 0;
					mn = 0;
					//flag = vector<int>(time_list_slr.size(), 9);

					//for (const t_grecover_equation& recover_equ_t : _all_recover->get_first_phase_recover_equation(site_list_slr[ista - 1], sat_list_slr[isat - 1], freq))
					//{
					//	//int idx = (recover_equ_t.time - _all_recover->get_equ_beg_time()) / _all_recover->get_interval();
					//	flag[] = 10;
					//	bcd++;
					//}

					for (size_t iepo = 0; iepo < time_list_slr.size(); iepo++)
					{
						vector<t_grecover_equation*> recover_equ = site_equ[site_list_slr[ista - 1]][sat_list_slr[isat - 1]][time_list_slr[iepo]];
						if (recover_equ.size() != 0)
						{
							_nepotbl_slr[ista][isat] += 1;
							_rmstable_slr["SLR"][ista][isat] += (recover_equ[0]->resuidal*recover_equ[0]->resuidal / time_list_slr.size());
							mn++;
						}
						
					}

					_nepotbl_slr[ista][0] = _nepotbl_slr[ista][0] + _nepotbl_slr[ista][isat];
					_nepotbl_slr[0][isat] = _nepotbl_slr[0][isat] + _nepotbl_slr[ista][isat];
					_rmstable_slr["SLR"][ista][0] = _rmstable_slr["SLR"][ista][0] + _rmstable_slr["SLR"][ista][isat];
					_rmstable_slr["SLR"][0][isat] = _rmstable_slr["SLR"][0][isat] + _rmstable_slr["SLR"][ista][isat];
				}

				_nepotbl_slr[0][0] = _nepotbl_slr[0][0] + _nepotbl_slr[0][isat];
				_rmstable_slr["SLR"][0][0] = _rmstable_slr["SLR"][0][0] + _rmstable_slr["SLR"][0][isat];

			}

			for (int ista = 0; ista <= site_list_slr.size(); ista++)
			{
				for (int isat = 0; isat <= sat_list_slr.size(); isat++)
				{
					if (_nepotbl_slr[ista][isat] > 0)
					{
						_rmstable_slr["SLR"][ista][isat] = sqrt(_rmstable_slr["SLR"][ista][isat] / (_nepotbl_slr[ista][isat] * 1.0 / time_list_slr.size()));
					}
					else
					{
						_rmstable_slr["SLR"][ista][isat] = 0;
					}
				}
			}

		}

		for (size_t isat = 1; isat <= sat_list.size(); isat++)
		{
			for (size_t ista = 1; ista <= site_list.size(); ista++)
			{
				flag = vector<int>(time_list.size(), 9);

				for (const t_grecover_equation& recover_equ_t : _all_recover->get_first_phase_recover_equation(site_list[ista - 1], sat_list[isat - 1], freq))
				{
					int idx = (recover_equ_t.time - _all_recover->get_equ_beg_time()) / _all_recover->get_interval();
					flag[idx] = recover_equ_t.is_newamb;
				}

				for (int iepo = 0; iepo < time_list.size(); iepo++)
				{
					if (t_geditres::istrue(flag[iepo], "OK"))
					{
						_nepotbl[ista][isat] += 1;
						vector<t_grecover_equation*> recover_equ = site_equ[site_list[ista - 1]][sat_list[isat - 1]][time_list[iepo]];
						for (int iobs = 0; iobs < recover_equ.size(); iobs++)
						{
							//pseudorange code
							if (recover_equ[iobs]->obstype.is_code())
							{
								_rmstable[t_gobscombtype("PC")][ista][isat] += (recover_equ[iobs]->resuidal*recover_equ[iobs]->resuidal / time_list.size());
							}
							else if(recover_equ[iobs]->obstype.is_phase())
							{
								_rmstable[t_gobscombtype("LC")][ista][isat] += (recover_equ[iobs]->resuidal*recover_equ[iobs]->resuidal / time_list.size());
							}
							
						}
						
					}
					
				}

				_nepotbl[ista][0] = _nepotbl[ista][0] + _nepotbl[ista][isat];
				_nepotbl[0][isat] = _nepotbl[0][isat] + _nepotbl[ista][isat];
				_rmstable[t_gobscombtype("LC")][ista][0] = _rmstable[t_gobscombtype("LC")][ista][0] + _rmstable[t_gobscombtype("LC")][ista][isat];
				_rmstable[t_gobscombtype("LC")][0][isat] = _rmstable[t_gobscombtype("LC")][0][isat] + _rmstable[t_gobscombtype("LC")][ista][isat];
				_rmstable[t_gobscombtype("PC")][ista][0] = _rmstable[t_gobscombtype("PC")][ista][0] + _rmstable[t_gobscombtype("PC")][ista][isat];
				_rmstable[t_gobscombtype("PC")][0][isat] = _rmstable[t_gobscombtype("PC")][0][isat] + _rmstable[t_gobscombtype("PC")][ista][isat];
				
			}

			_nepotbl[0][0] = _nepotbl[0][0] + _nepotbl[0][isat];
			_rmstable[t_gobscombtype("LC")][0][0] = _rmstable[t_gobscombtype("LC")][0][0] + _rmstable[t_gobscombtype("LC")][0][isat];
			_rmstable[t_gobscombtype("PC")][0][0] = _rmstable[t_gobscombtype("PC")][0][0] + _rmstable[t_gobscombtype("PC")][0][isat];

			
		}

		for (int ista = 0; ista <= site_list.size(); ista++)
		{
			for (int isat = 0; isat <= sat_list.size(); isat++)
			{
				if (_nepotbl[ista][isat] > 0)
				{
					_rmstable[t_gobscombtype("LC")][ista][isat] = sqrt(_rmstable[t_gobscombtype("LC")][ista][isat] / (_nepotbl[ista][isat] * 1.0 / time_list.size()));
					_rmstable[t_gobscombtype("PC")][ista][isat] = sqrt(_rmstable[t_gobscombtype("PC")][ista][isat] / (_nepotbl[ista][isat] * 1.0 / time_list.size()));
					
				}
				else
				{
					_rmstable[t_gobscombtype("LC")][ista][isat] = 0;
					_rmstable[t_gobscombtype("PC")][ista][isat] = 0;
				}

				
			}
		}

		_all_sum->setSite(site_list_v);
		_all_sum->setSat(sat_list_v);
		_all_sum->setNEPO(_nepotbl);
		_all_sum->setRMS(_rmstable);
		_all_sum->setNEPO_SLR(_nepotbl_slr);
		_all_sum->setRMS_SLR(_rmstable_slr);

		t_gfile gout;
		string sum_path = "sum";

		gout.path(sum_path);
		gout.glog(_log);

		t_sum sum_coder(_gset);
		sum_coder.add_data("ID1", _all_sum);
		gout.coder(&sum_coder);

		gout.run_write();
	}

	void t_geditres::update_logfiles()
	{
		for (string logfile : _logfiles)
		{
			string site_name = logfile.substr(logfile.rfind(".log") - 12, 4);
			if (_all_ambflag->getAllAmbFlag().find(site_name) == _all_ambflag->getAllAmbFlag().end()) {
				continue;
			}
			t_gio* gout = new t_gfile;
			//set I/O
			gout->path(logfile); gout->glog(_log);
			//set coder
			t_gcoder* gcoder = new t_ambflag(nullptr, "");
			gcoder->clear(); gcoder->glog(_log); gcoder->path(logfile);
			gcoder->add_data("ID1", _all_ambflag);
			gout->coder(gcoder);
			// write
			gout->run_write();
			delete gout;
			delete gcoder;
		}
	}

	void t_geditres::_check_resid()
	{
		stringstream loginfo;

		loginfo << " Check jump and remove short piece (" << jump_elisp << bad_elisp << short_elisp;
		_log->comment(t_glog::LOG_LV::LOG_INFO, "t_geditres", "_check_resid", loginfo.str());

		vector<t_gtime> time_list = _all_recover->get_all_obs_time_list();

		for (string site : _all_recover->get_site_list())
		{
			_log->comment(t_glog::LOG_LV::LOG_INFO, "t_geditres", "_check_resid", "Check resuid of Site: " + site);
			_res_ambflag[site] = map<string, vector<int> >();
			for (string sat : _all_recover->get_sat_list())
			{
				_res_ambflag[site][sat] = vector<int>(time_list.size(), 9);
				vector<int> flag = vector<int>(time_list.size(), 9);
				vector<double> value = vector<double>(time_list.size(), 0.0);
				for (const t_grecover_equation& recover_equ : _all_recover->get_first_phase_recover_equation(site, sat, freq))
				{
					int idx = (recover_equ.time - _all_recover->get_equ_beg_time()) / _all_recover->get_interval();
					flag[idx] = recover_equ.is_newamb;
					_res_ambflag[site][sat][idx] = (recover_equ.is_newamb);
					value[idx] = recover_equ.resuidal * 1E3;
				}


#ifdef DEBUG
				cout << TEST_INFO << " initial flag &&  value " << endl;
				print_flag_debuginfo(site, sat, flag);
				cout << TEST << setw(5) << site << setw(5) << sat;
				for (double v : value) {
					cout << setw(15) << setprecision(6) << setiosflags(ios::fixed) << v;
				}
				cout << endl;
#endif // DEBUG

				//remove short piece
				_remove_short(flag);

#ifdef DEBUG
				cout << TEST_INFO << "After first remove short" << endl;
				print_flag_debuginfo(site, sat, flag);
#endif // DEBUG


				//check jump and bad (0/1 degree poly. fitting)
				for (int iter = 0; iter < 10; iter++) {
					if (!_check_jump(flag, value)) {
						break;
					}
				}

#ifdef DEBUG
				cout << TEST_INFO << "After check jump " << endl;
				print_flag_debuginfo(site, sat, flag);
#endif // DEBUG



				//remove new short piece because of newly found ambiguities
				_remove_short(flag);


#ifdef DEBUG
				cout << TEST_INFO << "After second check jump" << endl;
				print_flag_debuginfo(site, sat, flag);

				cout << TEST_INFO << "inital logfile value" << endl;
				print_flag_debuginfo(site, sat, _res_ambflag[site][sat]);

#endif // DEBUG


				//TODO: is need to give the value to resls aging??????
				//restore flg, and  add amb. flag at first obs.
				bool is_found = false;
				for (int idx = 0; idx < time_list.size(); idx++) {

					_res_ambflag[site][sat][idx] = flag[idx] - int(flag[idx] / 10) * 10;

					if (!is_found && t_geditres::istrue(_res_ambflag[site][sat][idx], "ok")) {
						is_found = true;
						_res_ambflag[site][sat][idx] = t_geditres::OLDAMB;
					}
				}

#ifdef DEBUG

				cout << TEST_INFO << "update log value:" << endl;
				print_flag_debuginfo(site, sat, _res_ambflag[site][sat]);

#endif // DEBUG

			}
		}
		return;

	}

	int t_geditres::_find_flag(int istart, int istop, const vector<int>& flg, string flag_description, string direction)
	{
		int  i1, i2, istep;

		if (direction == "+") {
			i1 = istart;
			i2 = istop;
			istep = 1;
		}
		else {
			i1 = istop;
			i2 = istart;
			istep = -1;
		}

		// return value;
		int iepo = flg.size();
		for (int i = i1; i != i2 + istep; i += istep) {
			if (!t_geditres::istrue(flg[i], flag_description)) {
				continue;
			}
			iepo = i;
			break;
		}

		return iepo;

	}

	void t_geditres::_find_flag_upd(int istart, int istop, const vector<int>& flg, string flag_desription, string direction, int & iepo, int i0, int RLscale, int & j1, int & j2)
	{
		bool is_first = true;

		int  i1, i2, istep;


		if (direction == "+") {
			i1 = istart;
			i2 = istop;
			istep = 1;
		}
		else {
			i1 = istop;
			i2 = istart;
			istep = -1;
		}

		iepo = flg.size();

		for (int i = i1; i != i2 + istep; i += istep) {
			if (!t_geditres::istrue(flg[i], flag_desription)) {
				continue;
			}

			if (is_first) {
				j1 = i;
				is_first = false;
			}

			if (direction == "-") {
				if ((i - i0) % RLscale != 0) {
					continue;
				}
			}

			iepo = i;
			j2 = i;
			break;
		}

		return;
	}

	bool t_geditres::_check_jump(vector<int>& flag_data, const vector<double>& val_data)
	{

		size_t nepoch = flag_data.size();
		// return value
		bool is_jump = false;
		int imax, ndif, jump1, jbad1, first_epo;
		imax = ndif = jump1 = jbad1 = first_epo = 0;
		double sig, mean, difmax, signew, meannew, dummy, res;
		sig = mean = difmax = signew = meannew = dummy = res = 0.0;
		vector<double> dif(flag_data.size(), 0.0);
		vector<int> ipt(nepoch + 1, 0);
		bool  found_bad;

		for (int i = nepoch - 1; i >= 1; i--) {
			if (!t_geditres::istrue(flag_data[i], "OK")) {
				continue;
			}
			for (int j = i - 1; j >= 0; j--) {
				if (!t_geditres::istrue(flag_data[j], "OK")) {
					continue;
				}
				first_epo = i;
				dif[i] = val_data[i] - val_data[j];
				ipt[i] = j + 1;

				//sigma, and the largest difference

				if (t_geditres::istrue(flag_data[i], "GOOD")) {
					if (abs(dif[i]) > difmax) {
						difmax = abs(dif[i]);
						imax = i;
					}
					ndif++;
					mean += dif[i] / nepoch;
					sig += dif[i] * dif[i] / nepoch;
				}
				break;
			}
		}

		if (ndif <= 1) {
			return is_jump;
		}

		// mean and sigma of all 
		mean *= nepoch * 1.0 / ndif;
		sig = sig * (nepoch*1.0 / ndif) - mean * mean;

		if (sig <= (int64_t)bad_elisp * bad_elisp / 9.0) {
			sig = bad_elisp / 3.0;
		}
		else {
			sig = sqrt(sig);
		}

		// should make it work for extremely bad case(e.g. with more than
		//  one large bad or jump)
		if (sig > bad_elisp / 3) {
			meannew = mean - dif[imax] / nepoch;
			signew = sig - dif[imax] * dif[imax] / nepoch;

			meannew *= nepoch * 1.0 / ((int64_t)ndif - 1);
			signew = signew * (nepoch*1.0 / ((int64_t)ndif - 1)) - meannew * meannew;

			if (signew <= (int64_t)bad_elisp * bad_elisp / 9.0) {
				signew = bad_elisp / 3.0;
			}
			else {
				signew = sqrt(signew);
			}
			if (sig > 1.6*signew) {
				if (difmax > jump_elisp) {
					writeLogInfo(_log, 1, " new large jump at pos.:" + to_string(imax) + to_string(difmax) + to_string(jump_elisp));
					flag_data[imax] = t_geditres::NEWAMB + (ipt[imax]) * 10;
					is_jump = true;
					_check_jump(flag_data, val_data);
					return is_jump;
				}
			}
			sig = signew;
		}

		//!! adapt jump0 and jbad0 according to the sigma
		if (jump_elisp / 3.0 < sig) {
			jump1 = int(sig*3.0);
		}
		else {
			jump1 = jump_elisp;
		}
		if (bad_elisp / 3.0 < sig) {
			jbad1 = int(sig*3.0);
		}
		else {
			jbad1 = bad_elisp;
		}

		// check jump
		for (int i = first_epo; i < nepoch; i++) {
			if (t_geditres::istrue(flag_data[i], "OK")) {
				res = abs(dif[i]);

				if (res > jump1) {
					if (t_geditres::istrue(flag_data[i], "GOOD")) {
						writeLogInfo(_log, 1, "possible new jump.at pos.:" + to_string(i) + to_string(res) + to_string(jump1));
						flag_data[i] = t_geditres::NEWAMB + 10 * (ipt[i]);
						is_jump = true;
					}
				}
				else if (res > jbad1) {
					if (t_geditres::istrue(flag_data[i], "GOOD")) {
						writeLogInfo(_log, 1, "possible new bad.at pos.:" + to_string(i) + to_string(res) + to_string(jbad1));
						flag_data[i] = t_geditres::NEWBAD + 10 * (ipt[i]);
						found_bad = true;
					}
				}
				else {
					if (t_geditres::istrue(flag_data[i], "NEWAMB")) {
						flag_data[i] = t_geditres::GOOD;
					}
				}
			}
		}

		//!!handling for newbad and amb.flags one after another or flags at
		//	!!begin or end
		//	!!1. the 2th obs.is flagged as amb.newbad, the 1th obs.is removed
		//	!!2. the last  is flagged as amb. new bad, then remove it
		//	!!3. a single new bad, then set it to GOOD
		//	!!4. two flags one after another
		//	!!amb    amb, remove first one as bad.
		//	!!newbad, remove first one as bad, set the second as amb
		//	!!newbad amb, remove first one as bad.
		//	!!newbad, remove first one and set next one to GOOD
		for (int i = first_epo; i < nepoch; i++) {
			if (!t_geditres::istrue(flag_data[i], "AMB") &&
				!t_geditres::istrue(flag_data[i], "NEWBAD")) {
				continue;
			}
			int j = i + 1;
			while (ipt[j] == 0 && j < nepoch) {
				j++;
			}
			if (ipt[i] - 1 < first_epo) {
				flag_data[ipt[i] - 1] = t_geditres::DELBAD + 10 * (i + 1);
				flag_data[i] = t_geditres::GOOD;
			}
			else if (j >= nepoch) {
				flag_data[i] = t_geditres::DELBAD + ipt[i] * 10;
			}
			else if (t_geditres::istrue(flag_data[j], "AMB") || t_geditres::istrue(flag_data[j], "NEWBAD")) {
				if (t_geditres::istrue(flag_data[i], "OLDAMB")) {
					flag_data[i] = t_geditres::DELBAD + (j + 1) * 10;
					flag_data[j] = t_geditres::OLDAMB;
				}
				else if (t_geditres::istrue(flag_data[i], "NEWAMB")) {
					if (t_geditres::istrue(flag_data[j], "NEWAMB")) {
						flag_data[i] = t_geditres::DELBAD + ipt[i] * 10;
						flag_data[j] = t_geditres::NEWAMB + ipt[i] * 10;
					}
					else if (t_geditres::istrue(flag_data[j], "OLDAMB")) {
						flag_data[i] = t_geditres::DELBAD + ipt[i] * 10;
					}
					else {
						flag_data[j] = t_geditres::NEWAMB + ipt[i] * 10;
						flag_data[i] = t_geditres::DELBAD + ipt[i] * 10;
					}
				}
				else if (t_geditres::istrue(flag_data[i], "NEWBAD")) {
					if (t_geditres::istrue(flag_data[j], "AMB")) {
						flag_data[i] = t_geditres::DELBAD + ipt[i] * 10;
					}
					else {
						flag_data[i] = t_geditres::DELBAD + ipt[i] * 10;
						flag_data[j] = t_geditres::GOOD;
					}
				}
			}
			else {
				if (t_geditres::istrue(flag_data[i], "NEWBAD")) {
					flag_data[i] = t_geditres::GOOD;
				}
			}
		}

		return is_jump;
	}


	bool t_geditres::_remove_short(vector<int>& flag_data)
	{

		int nepoch = flag_data.size();
		bool is_removed = false;
		bool found = true;
		int j = 0;
		int i, k;
		while (found) {
			if (j == 0) {
				i = _find_flag(j, nepoch - 1, flag_data, "OK", "+");
				if (i >= nepoch) {
					return false;
				}
				if (!t_geditres::istrue(flag_data[i], "GOOD")) {
					flag_data[i] = t_geditres::GOOD;
				}
			}
			else {
				i = _find_flag(j, nepoch - 1, flag_data, "OLDAMB", "+");
			}

			if (i != nepoch) {
				j = _find_flag(i + 1, nepoch - 1, flag_data, "OLDAMB", "+");
				k = _find_flag(i, j - 1, flag_data, "OK", "-");
				int ngood = 1;
				for (int iepo = i + 1; iepo <= k; iepo++) {
					if (flag_data[iepo] == t_geditres::GOOD) {
						ngood++;
					}
				}
				if (ngood <= short_elisp) {
					is_removed = true;

					for (int iepo = i; iepo <= k; iepo++) {
						if (t_geditres::istrue(flag_data[iepo], "OK")) {
							flag_data[iepo] = t_geditres::DELSHRT + ngood * 10;
							writeLogInfo(_log, 1, "epoch " + to_string(iepo) + " lremoved as short");
						}
					}
				}
			}
			else {
				found = false;
			}
		}
		return is_removed;
	}

	void t_geditres::_update_ambflagdata()
	{

		double window = _all_recover->get_interval()*0.1;
		// cycle all the site 
		size_t nepoch = 0;
		set<string> sat_list;
		auto& all_site_ambflag = _all_ambflag->getAllAmbFlag();
		for (string site : _all_recover->get_site_list()) {
			string strlow_site = str2lower(site);
			//transform(site.begin(), site.end(), site.begin(), tolower);
			if (all_site_ambflag.find(strlow_site) == all_site_ambflag.end()) {
				_log->comment(t_glog::LOG_LV::LOG_ERROR, "t_geditres", "_update_ambflagdata", site + " logfile don't eixst! Skip!!!!");
				continue;
			}
			_log->comment(t_glog::LOG_LV::LOG_INFO, "t_geditres", "_update_ambflagdata", "Updating the ambflagdata of the site: " + site);

			t_gambflag& site_ambflag = all_site_ambflag[strlow_site];

			auto ambflag_header = site_ambflag.getAmbFlagHead();

			// get site_list
			set<string> sat_list_now;
			sat_list_now = site_ambflag.getAllSatSet();
			set<string> sat_temp;
			set_union(sat_list.begin(), sat_list.end(),
				sat_list_now.begin(), sat_list_now.end(),
				inserter(sat_temp, sat_temp.begin()));
			sat_list = sat_temp;

			// get number of epoch [note:begin from beg epoch in log header]
			//t_gtime ambflag_beg_time = t_gtime(ambflag_header->beg_mjd, int(ambflag_header->beg_sod), ambflag_header->beg_sod - int(ambflag_header->beg_sod));
			t_gtime ambflag_beg_time;
			ambflag_beg_time.from_mjd(ambflag_header->beg_mjd, int(ambflag_header->beg_sod), ambflag_header->beg_sod - int(ambflag_header->beg_sod));
			int nepoch_now = (site_ambflag.getMaxEndEpoch() - ambflag_beg_time) / ambflag_header->intv + 1;
			if (nepoch_now > nepoch) {
				nepoch = nepoch_now;
			}

			// rest ambflag_header
			ambflag_header->old_rm_obs += ambflag_header->new_rm_obs;
			ambflag_header->exist_amb += ambflag_header->new_amb;

			//Difference of then beginning times of log-file and res-file (log-res)
			double difs = ambflag_beg_time - _all_recover->get_equ_beg_time();
			int istart = int(difs / ambflag_header->intv);
			int RLScale = round(_all_recover->get_interval() / ambflag_header->intv);
			int nshort = short_elisp * RLScale;

			// intersection sat list in log and recover

			set<string> sat_final;
			set<string> recover_sat_list = _all_recover->get_sat_list();
			set_union(sat_list.begin(), sat_list.end(),
				recover_sat_list.begin(), recover_sat_list.end(),
				inserter(sat_final, sat_final.begin()));
			sat_list = sat_final;

#ifdef DEBUG
			cout << TEST
				<< setw(8) << site
				<< setw(8) << "HDINFO"
				<< setw(5) << nepoch
				<< setw(5) << sat_list.size()
				<< setw(5) << ambflag_header->old_rm_obs
				<< setw(5) << ambflag_header->exist_amb
				<< setw(5) << istart + 1
				<< setw(5) << RLScale
				<< setw(5) << nshort;
			cout << endl;
#endif // DEBUG

			// initial data struct for record ambflag 
			map<string, vector<int> > logamb, flag;
			for (string sat : sat_list) {
				logamb[sat] = vector<int>((int64_t)nepoch + 1, 9); // why +1 ?
				flag[sat] = vector<int>((int64_t)nepoch + 1, 0);
			}

			// set logamb by logfile ambflag data
			for (auto sat_records : site_ambflag.getAmbFlagData()) {
				string sat = sat_records.first;
				for (auto record : sat_records.second) {
					logamb[sat].at((int64_t)record->beg_epo - 1) = atoi(record->iflag.c_str());

					if (record->identify == "DEL") {
						for (int i = record->beg_epo; i <= record->end_epo; i++)
							logamb[sat].at((int64_t)i - 1) = atoi(record->iflag.c_str());
					}
					else if (record->identify == "AMB") {
						for (int i = record->beg_epo + 1; i <= record->end_epo; i++)
							logamb[sat].at((int64_t)i - 1) = 0;
					}
					else if (record->identify == "BAD") {
						// nothing to do
					}
					else {
						writeLogInfo(_log, 1, " unknown flag type in log-file :" + record->identify);
					}
				}
			}

#ifdef DEBUG
			for (string sat : sat_list) {
				print_flag_debuginfo(site, sat, logamb[sat], nepoch);
			}

#endif // DEBUG

			// merge flags from logfile and that from the newly editing(read from resfile)
			// initial the new ambflag header
			ambflag_header->avaiable_obs = 0;
			ambflag_header->new_rm_obs = 0;
			ambflag_header->new_amb = 0;

			// merge flags ,storage in logamb data structor
			size_t nepo_res = _all_recover->get_all_obs_time_list().size();
			for (size_t iepo = 0; iepo < nepo_res; ++iepo) {
				// aligned to firts epoch of logfile
				double dt = iepo * _all_recover->get_interval() - difs;
				// epoch of resfile is before the beg time of logfile,skipped
				if (dt < -window) continue;
				// epoch time does not match, skipped
				size_t jepo = round(dt / ambflag_header->intv);
				// epoch time bigger than nepoch skipped
				if (jepo > nepoch - 1) break;
				// last epoch of residuals earlier than that of ambflag file, which may lead to unexpected errors
				// e.g. interval of residuals is 300s, while that of ambflag file is 30s, 287 => 2871 288 => 2881, 
				// data from 2871-2880 now will be flaged as delet, Aug. 2020, by jqwu
				if (iepo == nepo_res - 1) {
					if (jepo < nepoch - RLScale) {
						for (int i = jepo + RLScale; i < nepoch; ++i) {
							for (string sat : sat_list)
								logamb[sat].at(i) = 6;
						}
					}
				}
				if (fabs(jepo*ambflag_header->intv - dt) > window) continue;

				for (string sat : _all_recover->get_sat_list()) {
					//cout << setw(6) << site << sat << iepo << _res_ambflag[site][sat].size() << endl;
					int k = _res_ambflag[site][sat][iepo];

					if (t_geditres::istrue(k, "nodata")) {
						if (t_geditres::istrue(logamb[sat].at(jepo), "ok")) {
							// move the new ambiguity start flag to next epoch, Aug 2020, jqwu
							if (t_geditres::istrue(logamb[sat].at(jepo), "amb") && iepo < nepo_res-1 && jepo < nepoch-1) {
								if (t_geditres::istrue(logamb[sat].at(jepo+1), "amb"))
									logamb[sat].at(jepo+1) = 1;
								if (_res_ambflag[site][sat][iepo + 1])
									_res_ambflag[site][sat][iepo + 1] = 1;
							}
							//!'DELBAD'   !! delete all those obs not used in this solution at the right epoch
							logamb[sat].at(jepo) = 6;
						}
						continue;
					}

					if (t_geditres::istrue(k, "ok")) {
						ambflag_header->avaiable_obs++;
					}

					if (t_geditres::istrue(k, "newdel")) {
						ambflag_header->new_rm_obs++;
					}

					if (t_geditres::istrue(k, "newamb")) {
						ambflag_header->new_amb++;
					}

					if (!t_geditres::istrue(k, "ok") && t_geditres::istrue(logamb[sat].at(jepo), "amb")) {
						for (int i = jepo + 1; i < nepoch; i++) {
							if (t_geditres::istrue(logamb[sat].at(i), "ok")) {
								// 'SHTAMB'
								logamb[sat].at(i) = 3;
								break;
							}
						}
					}

					logamb[sat].at(jepo) = k;
				}
			}

#ifdef DEBUG
			cout << TEST
				<< setw(8) << site
				<< setw(8) << "HDINFO"
				<< setw(5) << ambflag_header->avaiable_obs
				<< setw(5) << ambflag_header->new_rm_obs
				<< setw(5) << ambflag_header->new_amb
				<< endl;

			for (string sat : sat_list) {
				print_flag_debuginfo(site, sat, logamb[sat], nepoch);
			}

#endif // DEBUG


			//!! prepare for output log-file
			//!!compute the largest # of ambiguities per epoch and the span	
			// TODO for What??????
			bool lend = false;
			int j1, j2;
			int namb = 0;
			vector<int> namb_epo((int64_t)nepoch + 1, 0);
			for (string sat : sat_list) {
				// init flg storage logamb old data
				vector<int> flg = vector<int>(logamb[sat].begin(), logamb[sat].end());
				// delete the good epoches in ambflag file before the begin epoch of residuals, Aug 2020, jqwu
				for (size_t iepo = 0; iepo < nepoch; ++iepo) {
					if (t_geditres::istrue(flg[iepo], "good"))
						flg[iepo] = 6;
					if (t_geditres::istrue(flg[iepo], "amb"))
						break;
				}
				// cycle the all epoch
				size_t iepo = 0;
				while (iepo < nepoch) {
					if (!t_geditres::istrue(flg[iepo], "nodata")) {

						//removed observation, find out the last epoch of the removed piece
						if (!t_geditres::istrue(flg[iepo], "ok")) {

#ifdef DEBUG
							cout << "ok before:" << to_string(iepo) << endl;
#endif
							int jepo = iepo;
							while (jepo < nepoch && !t_geditres::istrue(flg[jepo], "ok")) {
								if (!t_geditres::istrue(flg[jepo], "nodata")) {
									flag[sat][iepo] = jepo + 1;
								}
								jepo++;
							}
							iepo = jepo;
#ifdef DEBUG
							cout << "ok after: iepo " << to_string(iepo) << " jepo " << to_string(jepo) << endl;

#endif
						}

						//ambiguity, find out the next ambiguity, then the last good obs. before
						//this ambiguity.That is all the continous data following the ambiguity
						else {
							if (t_geditres::istrue(flg[iepo], "amb")) {
#ifdef DEBUG
								cout << "amb : " << to_string(iepo) << endl;
#endif // DEBUG

								int jepo = _find_flag(iepo + 1, nepoch - 1, flg, "amb", "+");

#ifdef DEBUG
								cout << jepo << endl;
#endif // DEBUG

								if (jepo >= nepoch) {
									lend = true;
								}

								int jnstart = jepo;
								_find_flag_upd(iepo + 1, jepo - 1, flg, "good", "-", jepo, istart, RLScale, j1, j2);

#ifdef DEBUG
								cout << jepo << endl;
#endif // DEBUG


								if (jepo >= nepoch) {
									jepo = iepo;
								}

#ifdef DEBUG
								cout << "jepo -ipeo before : " << to_string(iepo) << " jepo " << to_string(jepo) << endl;
#endif // DEBUG


								if (jepo - iepo < nshort) {
									stringstream loginfo("DELETING short ambiguity ");
									loginfo << " " << iepo << " " << jepo << " " << sat;
									writeLogInfo(_log, 1, loginfo.str());
									if (ambflag_header->new_amb != 0) {
										ambflag_header->new_amb--;
									}
									else if (ambflag_header->exist_amb != 0) {
										ambflag_header->exist_amb--;
									}

									flag[sat][iepo] = jepo + 1;
									for (int i = iepo; i <= jepo; i++) {
										logamb[sat][i] = 6;
									}

									for (int jj = jepo; jj <= jnstart - 1; jj++) {
										if (t_geditres::istrue(logamb[sat][jj], "ok")) {
											logamb[sat][jj] = 6;
											flg[jj] = 6;
										}
									}
									iepo = jepo;

#ifdef DEBUG
									cout << "after deleting short : " << to_string(iepo) << " jepo " << to_string(jepo) << endl;
#endif // DEBUG

								}
								else {
									namb++;
									for (int i = iepo; i <= jepo; i++) {
										namb_epo[i]++;
									}
									flag[sat][iepo] = -(jepo + 1);

									if (lend && nepoch - jepo <= RLScale + 1) {
										flag[sat][iepo] = -nepoch;
										lend = false;
										iepo++;
										continue;
									}

									for (int jj = j2 + 1; jj <= j1; jj++) {
										if (t_geditres::istrue(flg[jj], "amb")) {
											writeLogInfo(_log, 1, "ERROR find_flag_upd missed amb flag");
											throw exception(logic_error(""));
										}
										if (t_geditres::istrue(flg[jj], "good")) {
											flg[jj] = 6;
											logamb[sat][jj] = 6;
										}
									}

#ifdef DEBUG
									cout << "no  deleting  : " << to_string(iepo) << " jepo " << to_string(jepo) << endl;
#endif // DEBUG


								}

								// 100 GOTO continue
						}
							iepo++;
					}
				}
					else {

#ifdef DEBUG
						cout << "nodata:" << to_string(iepo) << endl;
#endif // DEBUG

						iepo++;
					}

			} // site cycle

#ifdef DEBUG

				cout << TEST
					<< setw(8) << site
					<< setw(8) << sat
					<< setw(5) << ambflag_header->new_rm_obs
					<< setw(5) << ambflag_header->new_amb
					<< endl;


#endif // DEBUG

		}// sat cycle


#ifdef DEBUG
			for (string sat : sat_list) {
				print_flag_debuginfo(site, sat, logamb[sat], nepoch);
				print_flag_debuginfo(site, sat, flag[sat], nepoch);
			}

#endif // DEBUG

			ambflag_header->exist_amb = namb;
			ambflag_header->new_amb = 0;
			ambflag_header->max_amb_1epo = namb_epo[0];
			for (int number_amb : namb_epo) {
				if (ambflag_header->max_amb_1epo < number_amb) {
					ambflag_header->max_amb_1epo = number_amb;
				}
			}
#ifdef DEBUG
			cout << TEST
				<< setw(8) << site
				<< setw(8) << "HDINFO"
				<< setw(10) << ambflag_header->exist_amb
				<< setw(10) << ambflag_header->max_amb_1epo
				<< endl;
#endif // DEBUG


			// write_log header 
			// write log file
			// Don't need  chang into the logfile decoder
			// But need to change the record int all ambflagdata
			// generate ambflagdata

			// update ambflag data
			//auto site_ambflag_data = site_ambflag.getAmbFlagData();
			site_ambflag.clearAllAmbFlagData();
			for (int iepo = 0; iepo < nepoch; iepo++) {
				for (string sat : sat_list) {

					ambflag_data temp_ambflag_data;
					string cflag = "";
					string line = "";

					// get cflag and line value
					if (flag[sat][iepo] == iepo + 1) {
						cflag = "BAD";
						line = "delbad";
					}
					else if (flag[sat][iepo] > iepo + 1) {
						cflag = "DEL";
						line = "delmix";
					}
					else if (flag[sat][iepo] < 0) {
						cflag = "AMB";
						flag[sat][iepo] = -flag[sat][iepo];
						if (t_geditres::istrue(logamb[sat][iepo], "oldamb")) {
							line = "oldamb";
						}
						else if (t_geditres::istrue(logamb[sat][iepo], "newamb")) {
							line = "newamb";
						}
						else {
							writeLogInfo(_log, 0, "***ERROR: unknown data flag in res - file");
							throw exception(logic_error(""));
						}
					}

					// add new ambflag data
					if (cflag != "") {
						int iflag = 0;
						if (logamb[sat][iepo] >= 5 || logamb[sat][iepo] == 3) {
							iflag = 5;
						}
						else if (logamb[sat][iepo] < 5) {
							iflag = 1;
						}
						temp_ambflag_data.identify = cflag;
						temp_ambflag_data.beg_epo = iepo + 1;
						temp_ambflag_data.end_epo = flag[sat][iepo];
						temp_ambflag_data.iflag = to_string(iflag);
						temp_ambflag_data.C1 = 0.0;
						temp_ambflag_data.C2 = 0.0;
						temp_ambflag_data.reason = "RN_" + line;

						site_ambflag.addAmbFlagData(sat, temp_ambflag_data);
					}

				}
			}

			// 鎶婁慨鏀硅祴鍊煎洖鍘?
			all_site_ambflag[strlow_site] = site_ambflag;
		}
		return;
	}

	bool t_geditres::istrue(int flag_code, string flag_type)
	{
		bool istrue = false;

		int i = flag_code;
		string type = flag_type;

		if (flag_code >= 9) {
			i = flag_code - int(flag_code / 10) * 10;
		}
		transform(flag_type.begin(), flag_type.end(), type.begin(), ::toupper);

		type = trim(type);

		if (type == "OK") {
			if (i == t_geditres::GOOD || i == t_geditres::NEWAMB || i == t_geditres::OLDAMB) {
				istrue = true;
			}
		}
		else if (type == "GOOD") {
			if (i == t_geditres::GOOD) {
				istrue = true;
			}
		}
		else if (type == "NODATA") {
			if (i == t_geditres::NODATA) {
				istrue = true;
			}
		}
		else if (type == "AMB") {
			if (i == t_geditres::NEWAMB || i == t_geditres::OLDAMB) {
				istrue = true;
			}
		}
		else if (type == "NEWFLG") {
			if (i == t_geditres::NEWAMB || i == t_geditres::DELBAD
				|| i == t_geditres::DELSHRT || i == t_geditres::DELRMS)
			{
				istrue = true;
			}
		}
		else if (type == "NEWDEL") {
			if (i == t_geditres::DELBAD || i == t_geditres::DELSHRT
				|| i == t_geditres::DELRMS || i == t_geditres::NEWBAD)
			{
				istrue = true;
			}
		}
		else if (type == "NEWAMB") {
			if (i == t_geditres::NEWAMB) {
				istrue = true;
			}
		}
		else if (type == "OLDAMB") {
			if (i == t_geditres::OLDAMB) {
				istrue = true;
			}
		}
		else if (type == "DELRMS") {
			if (i == t_geditres::DELRMS) {
				istrue = true;
			}
		}
		else if (type == "DELBAD") {
			if (i == t_geditres::DELBAD) {
				istrue = true;
			}
		}
		else if (type == "DELSHRT") {
			if (i == t_geditres::DELSHRT) {
				istrue = true;
			}
		}
		else if (type == "NEWBAD") {
			if (i == t_geditres::NEWBAD) {
				istrue = true;
			}
		}
		else {
			string error_info = "***ERROR: unknow c input," + type;
			throw exception(logic_error(error_info.c_str()));
		}

		return istrue;
	}

}