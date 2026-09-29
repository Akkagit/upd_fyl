/**
* @verbatim
	History
	 -1.0 Guolong Feng	 2020-07-20  creat the file.
  @endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file			gclockrepair.cpp
* @author       Guolong Feng, Wuhan University
* @version		1.0.0
* @date			2020-07-20
*
*/

#include <iomanip>
#include <memory>
#include "gproc/gclockrepair.h"


namespace great {
	t_gclock_repair::t_gclock_repair(t_gsetbase* settings, t_gallobs* allobs)
	{
		_set = settings;
		_allobs = allobs;
		set<string> sys_list = dynamic_cast<t_gsetgen*>(_set)->sys();
		for (auto sys : sys_list)
		{
			GSYS gsys = t_gsys::str2gsys(sys);
			_band_index[gsys] = dynamic_cast<t_gsetgnss*>(_set)->band_index(gsys);
		}

		_site_list = dynamic_cast<t_gsetgen*>(_set)->recs();

		_jumpCor = 0;
	}

	void t_gclock_repair::repairClockJump(string site, const t_gtime& begT, const t_gtime& endT, double sampling, bool isRealtime)
	{
		_jumpCor = 0;
		t_gtime crt_time = begT;
		vector<shared_ptr<t_gobsgnss>>       epoData;
		vector<t_gsatdata>  epoDataOrig;
		while (crt_time <= endT)
		{
			epoData = _allobs->obs_pt(site, crt_time);
			epoDataOrig= _allobs->obs(site, crt_time);
			if (epoDataOrig.empty())
			{
				crt_time = crt_time + sampling;
				continue;
			}
			if (crt_time == begT)
			{
				_epoDataPre = epoData;
				_epoOrigPre = epoDataOrig;
				crt_time = crt_time + sampling;
				continue;
			}
			if(isRealtime)
			{ 
				this->_repairClockJump(_epoDataPre, epoData);
				_epoDataPre = epoData;
			}
			else
			{
				this->_repairClockJump(_epoOrigPre, epoData);
				_epoOrigPre = epoDataOrig;
			}

			// Next
			crt_time = crt_time + sampling;
		}

	}

	int t_gclock_repair::ProcessBatch(string site, const t_gtime& begT, const t_gtime& endT, double sampling)
	{
		// repairClockJump
		this->repairClockJump(site, begT, endT, sampling);

		// Write O

		string obspath = dynamic_cast<t_gsetout*>(_set)->outputs("obs_dir");
		if (obspath.empty())   obspath = "file://obs_repair";
		if (ACCESS(obspath.substr(7).c_str(), 0) != 0)
		{
			MKDIR(obspath.substr(7).c_str());
		}

		string site_tmp = site;
		transform(site_tmp.begin(), site_tmp.end(), site_tmp.begin(), ::tolower);

		stringstream obs_file;
		obs_file << obspath << PATH_SEPARATOR << site_tmp << setw(3) << setfill('0') << begT.doy() << "0."
			<< setw(2) << setfill('0') << begT.yr() << "o";

		if (ACCESS(obs_file.str().c_str(), 0) == 0)
		{
			string obs_file_orig = obs_file.str() + ".orig";
			int isok = rename(obs_file.str().c_str(), obs_file_orig.c_str());
			cout << site << "  " << isok << endl;
		}

		t_gio* gout = new t_gfile;
		//set I/O
		gout->path(obs_file.str()); gout->glog(_glog);
		//set coder
		t_gcoder* gcoder = new t_rinexo(nullptr, "", 4096);
		dynamic_cast<t_rinexo*>(gcoder)->setSite(site);
		gcoder->clear(); gcoder->glog(_glog); gcoder->path(obs_file.str());
		gcoder->add_data("IDX", _allobs);
		gcoder->add_data("OBJ", _allobj);
		gout->coder(gcoder);
		// write
		gout->run_write();
		delete gout; gout = nullptr;
		delete gcoder; gcoder = nullptr;
		return 1;
		
	}

	void t_gclock_repair::_repairClockJump(vector<t_spt_gobs> AllObsPre, vector<t_spt_gobs> AllObsCrt)
	{
		int Jump_flag_1 = this->_clock_jump_detection(AllObsPre, AllObsCrt);
		if (Jump_flag_1 > 0)
		{
			int Jump_flag_2 = this->_clock_jump_repair(Jump_flag_1, AllObsCrt);
		}
	}

	void t_gclock_repair::_repairClockJump(const vector<t_gsatdata>& AllDataOrig, vector<t_spt_gobs> AllObsCrt)
	{
		int Jump_P = this->_check_msec_jump(AllDataOrig, AllObsCrt, false);//check jump for range obs
		int Jump_L = this->_check_msec_jump(AllDataOrig, AllObsCrt, true); //check jump for phase ons

		if (Jump_P != Jump_L && (Jump_P != 0 || Jump_L != 0)) _jumpCor = _jumpCor + (Jump_P - Jump_L);

		if (_jumpCor == 0) return;
		//correct obs
		for (auto it = AllObsCrt.begin(); it != AllObsCrt.end(); it++)
		{
			(*it)->mod_L(1e-3 * _jumpCor * CLIGHT, X);  // modify all phases by dL [m]
		}

	}

	int t_gclock_repair::_check_msec_jump(const vector<t_gsatdata>& AllObsPre, vector<t_spt_gobs> AllObsCrt, bool isPhase)
	{
		int numAllSat = 0;    // all sats
		int numJumpSat = 0;   // Jump sats
		vector<int> jump;     // jump data for every sat
		for (auto crtObs : AllObsCrt)
		{
			GSYS gsys = crtObs->gsys();
			string sat = crtObs->sat();
			GOBSBAND b1 = _band_index[gsys][FREQ_1];

			double pre_value, now_value;
			pre_value = 0.0; now_value = 0.0;
			if (isPhase)
			{
				GOBS L_now = crtObs->select_phase(b1, true);
				now_value = crtObs->obs_L(L_now);
				for (auto preObs : AllObsPre)
				{
					if (preObs.sat() != sat) continue;
					GOBS L_pre = preObs.select_phase(b1, true);
					pre_value = preObs.obs_L(L_pre);
					break;
				}
			}
			else
			{
				GOBS P_now = crtObs->select_range(b1, true);
				now_value = crtObs->obs_C(P_now);
				for (auto preObs : AllObsPre)
				{
					if (preObs.sat() != sat) continue;
					GOBS P_pre = preObs.select_range(b1, true);
					pre_value = preObs.obs_C(P_pre);
					break;
				}
			}

			if (double_eq(now_value, 0.0) || double_eq(pre_value, 0.0)) continue;

			numAllSat++;

			double delt_jump = (now_value - pre_value)*1e3 / CLIGHT;
			int jump_tmp = round(delt_jump);
			if (abs(delt_jump - jump_tmp) > 0.15) jump_tmp = 0;
			else numJumpSat++;
			jump.push_back(jump_tmp);
			
		}
		if (numAllSat == 0 || numJumpSat == 0) return 0;

		int max_numSame = 0;
		int jump_final = 0;

		for (auto it = jump.begin(); it != jump.end(); it++)
		{
			if (*it == 0) continue;
			int tmp_same = 0;
			for (auto it_tmp = jump.begin(); it_tmp != jump.end(); it_tmp++)
			{
				if (*it_tmp == 0) continue;
				if (*it == *it_tmp) tmp_same++;
			}
			if (tmp_same > max_numSame && tmp_same > 1)
			{
				max_numSame = tmp_same;
				jump_final = *it;
				if (max_numSame == numJumpSat) return jump_final;
			}
		}

		if (max_numSame > (numJumpSat - max_numSame)) return jump_final; //found but not the same for all

		return 0;
				
	}


	int t_gclock_repair::_clock_jump_detection(vector<t_spt_gobs> AllObsPre, vector<t_spt_gobs> AllObsCrt)
	{
		_sumS = 0.0;
		int numAllSat = 0;    // all sats
		int numJumpSat = 0;   // Jump sats
		for (auto crtObs : AllObsCrt)
		{
			GOBS L_now, P_now;
			double Lv_now, Pv_now, Lv_pre, Pv_pre;
			Lv_now = Pv_now = Lv_pre = Pv_pre = 0.0;

			GSYS gsys = crtObs->gsys();
			GOBSBAND b1 = _band_index[gsys][FREQ_1];
			
			L_now = crtObs->select_phase(b1, true);
			P_now = crtObs->select_range(b1, true);

			Lv_now = crtObs->obs_L(L_now);
			Pv_now = crtObs->obs_C(P_now);

			if (double_eq(Lv_now, 0.0) || double_eq(Pv_now, 0.0)) continue;

			string sat = crtObs->sat();
			for (auto preObs : AllObsPre)
			{
				if (preObs->sat() != sat) continue;
				GOBS L_pre = preObs->select_phase(b1, true);
				GOBS P_pre = preObs->select_range(b1, true);
				Lv_pre = preObs->obs_L(L_pre);
				Pv_pre = preObs->obs_C(P_pre);
				break;
			}

			if (double_eq(Lv_pre, 0.0) || double_eq(Pv_pre, 0.0)) continue;

			numAllSat++;

			double delta_PL = (Pv_now - Pv_pre) - (Lv_now - Lv_pre);

			double k1 = 1e-3 * CLIGHT;

			// 对于某一历元，当且仅当所有可用卫星满足时，认为该历元时刻可能存在钟跳或所有卫星同时发生大周跳(进一步判断)
			if (abs(delta_PL) > k1) 
			{
#ifdef DEBUG
				cout << crtObs->epoch().str_hms() << " " << crtObs->sat() << " delta_P: " << fixed << setw(16) << setprecision(6) << (Pv_now - Pv_pre)
					<< " delta_L: " << setw(16) << setprecision(6) << (Lv_now - Lv_pre) << " delta_PL: " << setw(16) << setprecision(6) << delta_PL 
					<< " k1: " << setw(16) << setprecision(6) << k1 << endl;
#endif
				_sumS += delta_PL;
				numJumpSat++;
			}
		}

		if (numAllSat == numJumpSat) return numJumpSat;   // Maybe Jump
		else                         return 0;   // Not Jump
	}

	int t_gclock_repair::_clock_jump_repair(int numJump, vector<t_spt_gobs> AllObs)
	{
		int ClockJump = 0;   // integer ms clk jump  

		// 钟跳候选值M_candidate
		double M_candidate = (_sumS * 1e3) / (numJump * CLIGHT);
		// 阈值判断    1e-5 ~ 1e-7
		double thres_CJ = abs(M_candidate - int(M_candidate));
		// Not Jump
		if (thres_CJ > 1e-5)
		{
			return 0;
		}
		else
		{
			ClockJump = int(M_candidate);
			for (vector<t_spt_gobs>::iterator it = AllObs.begin(); it != AllObs.end(); it++)
			{
				(*it)->mod_L(1e-3 * ClockJump * CLIGHT, X);  // modify all phases by dL [m]
			}
#ifdef DEBUG
			cout  << fixed << " ClockJump: " << setw(18) << setprecision(7) << ClockJump
				<< " ms " << setw(18) << setprecision(7) << 1e-3 * ClockJump * CLIGHT << " meters " << endl;
#endif
			return 1;
		}
	}

} // namespace
