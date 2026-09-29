#include "gproc/gupdateparALL.h"
#include "gutils/gturboedit.h"

namespace great
{
	t_gupdateparALL::t_gupdateparALL(shared_ptr<t_gcycleslip> cycleslip, const map< GSYS, map<FREQ_SEQ, GOBSBAND> > & band_list, const IONO_ORDER& iono_order) :
		_band_list(band_list)
	{
		_cycleslip = cycleslip;
		_iono_order = iono_order;
		_freq = 2;
	}

	t_gupdateparALL::t_gupdateparALL(shared_ptr<t_gcycleslip> cycleslip_12, shared_ptr<t_gcycleslip> cycleslip_13, map< GSYS, map<FREQ_SEQ, GOBSBAND> > band_list, const IONO_ORDER& iono_order)
		: _band_list(band_list)
	{
		_cycleslip = cycleslip_12;
		_cycleslip_13 = cycleslip_13;
		_iono_order = iono_order;
		_freq = 3;
		//bool mjd = dynamic_cast<t_gturboedit*>(_cycleslip.get())->beg_mjd() == dynamic_cast<t_gturboedit*>(_cycleslip_13.get())->beg_mjd();
		//bool sod = dynamic_cast<t_gturboedit*>(_cycleslip.get())->beg_sod() == dynamic_cast<t_gturboedit*>(_cycleslip_13.get())->beg_sod();
		//bool itv = dynamic_cast<t_gturboedit*>(_cycleslip.get())->intv() == dynamic_cast<t_gturboedit*>(_cycleslip_13.get())->intv();
	
		//if (!mjd || !sod || !itv)
		//{
		//	throw runtime_error("Can't judge cycle slip unless mjd sod and intv is the same");
		//}
	}
	
	t_gupdateparALL::t_gupdateparALL(shared_ptr<t_gcycleslip> cycleslip_12, shared_ptr<t_gcycleslip> cycleslip_13, shared_ptr<t_gcycleslip> cycleslip_14, shared_ptr<t_gcycleslip> cycleslip_15, map< GSYS, map<FREQ_SEQ, GOBSBAND> > band_list, const IONO_ORDER& iono_order)
		: _band_list(band_list)
	{
		_cycleslip = cycleslip_12;
		_cycleslip_13 = cycleslip_13;
		_cycleslip_14 = cycleslip_14;
		_cycleslip_15 = cycleslip_15;
		_iono_order = iono_order;

		if (_cycleslip) _freq = 2;
		if (_cycleslip && _cycleslip_13) _freq = 3;
		if (_cycleslip && _cycleslip_13 && _cycleslip_14 && _cycleslip_15) _freq = 5;

	}

	t_gupdateparALL::~t_gupdateparALL()
	{
	}

	
	t_gupdateparinfo t_gupdateparALL::get_all_update_parameters(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata)
	{
		t_gupdateparinfo update_par_info;
		// update_ion_info
		_udpate_ion_pars(epoch, allpars, obsdata, update_par_info);

		// update gps rec ifb pars
		if(_cycleslip_13 != nullptr && allpars.orbParNumber() > 0) _udpate_gps_rec_ifb_pars(epoch, allpars, obsdata, update_par_info);

		// update_amb_info
		if (_lite_update_amb)
		{
			_lite_update_amb_pars(epoch, allpars, obsdata, update_par_info);
		}
		else
		{
			_update_amb_pars(epoch, allpars, obsdata, update_par_info);
		}
		

		// update_process_info
		_update_process_pars(epoch, allpars, update_par_info);

		return update_par_info;
	}

	//// jdhuang =======================================================================================================================================
	//void t_gupdateparALL::_update_amb_pars_triple_freq(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo & update_info)
	//{
	//	map<string, set<string> > mapPrn;
	//	// add new amb
	//	for (auto data : obsdata)
	//	{
	//		string rec = data.site();
	//		string sat = data.sat();
	//		GSYS   sys = data.gsys();

	//		mapPrn[rec].insert(sat);

	//		auto tb_slip = dynamic_cast<t_gturboedit*>(_cycleslip.get());
	//		auto tb_slip3 = dynamic_cast<t_gturboedit*>(_cycleslip_13.get());
	//		int idx_amb_L1 = allpars.getParam(rec, t_gpar::AMB_L1, sat);
	//		if (
	//			(idx_amb_L1 < 0 && 
	//			                   (_cycleslip->num_of_amb_arc(rec, sat, epoch) || 
	//							    _cycleslip_13->num_of_amb_arc(rec, sat, epoch))
	//			)
	//			||
	//			(idx_amb_L1 >= 0 && 
	//			                   (_cycleslip->cycle_slip123(data, epoch)) || 
	//				                _cycleslip_13->cycle_slip123(data, epoch))
	//			)
	//		{
	//			t_gtime end_time;

	//			// log tb for 1 and 2 or 1 and 3

	//			int  isFlag12 = _cycleslip->num_of_amb_arc(rec, sat, epoch);
	//			int  isFlag13 = _cycleslip_13->num_of_amb_arc(rec, sat, epoch);

	//			if(isFlag12 > 0) _cycleslip->set_amb_flag(rec, sat, isFlag12);
	//			if(isFlag13 >0) _cycleslip_13->set_amb_flag(rec, sat, isFlag13);

	//			if (isFlag12  > 0) end_time = tb_slip ? tb_slip->get_crt_amb_end(rec, sat) : LAST_TIME;
	//			if (isFlag13 > 0) end_time = tb_slip3 ? tb_slip3->get_crt_amb_end(rec, sat) : LAST_TIME;

	//			if (isFlag12 > 0 && isFlag13 > 0)
	//			{
	//				t_gtime end_time1 = tb_slip ? tb_slip->get_crt_amb_end(rec, sat) : LAST_TIME;
	//				t_gtime end_time2 = tb_slip3 ? tb_slip3->get_crt_amb_end(rec, sat) : LAST_TIME;
	//				end_time = end_time2.diff(end_time1) > 0 ? end_time1 : end_time2;
	//			}

	//			// 两个都没有观测弧段，直接跳过
	//			if (isFlag12 == 0 && isFlag13 == 0) continue;

	//			// get the band used
	//			for (const auto& it : _band_list[sys])
	//			{
	//				t_gpar::t_type amb_type = t_gpar::AMB_IF;

	//				if (t_gsys::band2freq(sys, it.second) == FREQ_2 && isFlag12 == 0) continue;
	//				if (t_gsys::band2freq(sys, it.second) == FREQ_3 && isFlag13 == 0) continue;

	//				     if (t_gsys::band2freq(sys, it.second) == FREQ_1) amb_type = t_gpar::AMB_L1;
	//				else if (t_gsys::band2freq(sys, it.second) == FREQ_2) amb_type = t_gpar::AMB_L2;
	//				else if (t_gsys::band2freq(sys, it.second) == FREQ_3) amb_type = t_gpar::AMB_L3;
	//				else if (t_gsys::band2freq(sys, it.second) == FREQ_4) amb_type = t_gpar::AMB_L4;
	//				else if (t_gsys::band2freq(sys, it.second) == FREQ_5) amb_type = t_gpar::AMB_L5;
	//				else {
	//					throw std::logic_error("Unknown band!!");
	//				};
	//				int idx_max = allpars[allpars.parNumber() - 1].index;
	//				t_gpar par_newamb = t_gpar(rec, amb_type, idx_max + 1, sat);

	//				par_newamb.value(0.0);
	//				par_newamb.setTime(epoch, end_time);
	//				par_newamb.apriori(9000.0);

	//				// add new amb
	//				int idx_amb = allpars.getParam(rec, amb_type, sat);
	//				if (idx_amb < 0)
	//				{
	//					//add_parameter(par_newamb);
	//					update_info.add(par_newamb);
	//				}
	//				// reset the new amb
	//				else
	//				{
	//					allpars[idx_amb].end = epoch - _intv;
	//					update_info.add(idx_amb + 1, par_newamb);
	//				}
	//			}
	//		}
	//	}

	//	// remove old amb
	//	for (unsigned int i = 0; i < allpars.parNumber(); i++)
	//	{
	//		if (allpars[i].parType == t_gpar::AMB_L1 ||
	//			allpars[i].parType == t_gpar::AMB_L2 ||
	//			allpars[i].parType == t_gpar::AMB_L3 ||
	//			allpars[i].parType == t_gpar::AMB_L4 ||
	//			allpars[i].parType == t_gpar::AMB_L5)
	//		{
	//			if (allpars[i].end < epoch)
	//			{
	//				if (!update_info.exist(i + 1))
	//				{
	//					update_info.add(i + 1);
	//				}
	//			}
	//			else
	//			{
	//				continue;
	//			}
	//		}
	//	}
	//}

	//void t_gupdateparALL::_update_amb_pars_triple_freq_new(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo & update_info)
	//{
	//	// add new amb
	//	auto tb_slip_12 = dynamic_cast<t_gturboedit*>(_cycleslip.get());
	//	auto tb_slip_13 = dynamic_cast<t_gturboedit*>(_cycleslip_13.get());

	//	if (!tb_slip_12 || !tb_slip_13) throw runtime_error("tb is empty!");
	//	for (auto data : obsdata)
	//	{
	//		string rec = data.site();
	//		string sat = data.sat();
	//		GSYS   sys = data.gsys();
	//		auto   addIdx = make_pair(rec, sat);

	//		bool isUpL1 = false; //update L1
	//		bool isUpL2 = false; //update L2
	//		bool isUpL3 = false; //update L3

	//		bool isFinalUpL1 = false;
	//		bool isFinalUpL2 = false;
	//		bool isFinalUpL3 = false;

	//		//case 1.1 idx_L2 < 0 && idx_L3 < 0
	//		//case 1.2 idx_L2 < 0 && idx_L3 > 0
	//		//case 1.2 idx_L2 > 0 && idx_L3 < 0
	//		//case 1.2 idx_L2 > 0 && idx_L3 > 0
	//		int idx_amb_L2 = allpars.getParam(rec, t_gpar::AMB_L2, sat);
	//		int idx_amb_L3 = allpars.getParam(rec, t_gpar::AMB_L3, sat);

	//		if (idx_amb_L2 >= 0 && idx_amb_L3 < 0) throw runtime_error("L1 L2 L3 should exist or remove at the same time.");;
	//		if (idx_amb_L2 < 0 && idx_amb_L3 >= 0) throw runtime_error("L1 L2 L3 should exist or remove at the same time.");;

	//		int Flag12 = tb_slip_12->num_of_amb_arc(rec, sat, epoch);
	//		int Flag13 = tb_slip_13->num_of_amb_arc(rec, sat, epoch);

	//		t_gtime end_time    = FIRST_TIME;
	//		t_gtime end_time_12 = FIRST_TIME;
	//		t_gtime end_time_13 = FIRST_TIME;

	//		if (idx_amb_L2 < 0 && idx_amb_L3 < 0)
	//		{
	//			if (Flag12 > 0)
	//			{
	//				tb_slip_12->set_amb_flag(rec, sat, Flag12);
	//				end_time_12 = tb_slip_12 ? tb_slip_12->get_crt_amb_end(rec, sat) : LAST_TIME;
	//				isUpL1 = true;
	//				isUpL2 = true;
	//				isAddL2ForL3[addIdx] = false;
	//			}

	//			if (Flag13 > 0)
	//			{
	//				tb_slip_13->set_amb_flag(rec, sat, Flag13);
	//				end_time_13 = tb_slip_13 ? tb_slip_13->get_crt_amb_end(rec, sat) : LAST_TIME;
	//				isUpL1 = true;
	//				isUpL3 = true;
	//				isAddL3ForL2[addIdx] = false;
	//			}
	//		}

	//		if (idx_amb_L2 >= 0 && idx_amb_L3 >= 0)
	//		{
	//			if (isAddL2ForL3[addIdx])
	//			{
	//				if (Flag12 > 0)
	//				{
	//					tb_slip_12->set_amb_flag(rec, sat, Flag12);
	//					end_time_12 = tb_slip_12 ? tb_slip_12->get_crt_amb_end(rec, sat) : LAST_TIME;
	//					isUpL1 = true;
	//					isUpL2 = true;
	//					isAddL2ForL3[addIdx] = false;
	//				}
	//			}
	//			else
	//			{
	//				bool Slip_12 = tb_slip_12->cycle_slip(data, epoch);
	//				if (Slip_12)
	//				{
	//					tb_slip_12->set_amb_flag(rec, sat, Flag12);
	//					end_time_12 = tb_slip_12 ? tb_slip_12->get_crt_amb_end(rec, sat) : LAST_TIME;
	//					isUpL1 = true;
	//					isUpL2 = true;
	//					if (Flag12 == 0) throw runtime_error("Flag12 is zero while idx_L2 >= 0");
	//				}
	//			}

	//			if (isAddL3ForL2[addIdx])
	//			{
	//				if (Flag13 > 0)
	//				{
	//					tb_slip_13->set_amb_flag(rec, sat, Flag13);
	//					end_time_13 = tb_slip_13 ? tb_slip_13->get_crt_amb_end(rec, sat) : LAST_TIME;
	//					isUpL1 = true;
	//					isUpL3 = true;
	//					isAddL3ForL2[addIdx] = false;
	//				}
	//			}
	//			else
	//			{
	//				bool Slip_13 = tb_slip_13->cycle_slip(data, epoch);
	//				if (Slip_13)
	//				{
	//					tb_slip_13->set_amb_flag(rec, sat, Flag13);
	//					end_time_13 = tb_slip_13 ? tb_slip_13->get_crt_amb_end(rec, sat) : LAST_TIME;
	//					isUpL1 = true;
	//					isUpL3 = true;
	//					if (Flag13 == 0) throw runtime_error("Flag13 is zero while idx_L3 >= 0");
	//				}
	//			}
	//		}

	//		     if (!isUpL2 && !isUpL3) continue;      // 无需更新
	//		else if (isUpL2  && !isUpL3) end_time = end_time_12;
	//		else if (!isUpL2 &&  isUpL3) end_time = end_time_13;
	//		else end_time = end_time_13.diff(end_time_12) < 0 ? end_time_12 : end_time_13;

	//		isFinalUpL1 = isUpL1;
	//		isFinalUpL2 = isUpL2;
	//		isFinalUpL3 = isUpL3;
	//		// when L2 and L3 not exist, if only one add, we both add another
	//		if (idx_amb_L2 < 0 && !isUpL2 &&  isUpL3)
	//		{
	//			isFinalUpL1 = true;
	//			isFinalUpL2 = true;
	//			isAddL2ForL3[make_pair(rec, sat)] = true;
	//		}

	//		if (idx_amb_L3 < 0 && isUpL2 && !isUpL3)
	//		{
	//			isFinalUpL1 = true;
	//			isFinalUpL3 = true;
	//			isAddL3ForL2[make_pair(rec, sat)] = true;
	//		}

	//		// when L2 and L3 exist, if only one update, we both update another
	//		if (idx_amb_L2 >= 0 && !isUpL2 &&  isUpL3)
	//		{
	//			isFinalUpL1 = true;
	//			isFinalUpL2 = true; // 没有更新，直接移除
	//			t_gtime end_time_L2 = allpars[idx_amb_L2].end;
	//			end_time = end_time.diff(end_time_L2) < 0 ? end_time_L2 : end_time;
	//		}
	//		if (idx_amb_L3 >= 0 &&  isUpL2 && !isUpL3)
	//		{
	//			isFinalUpL1 = true;
	//			isFinalUpL3 = true;
	//			t_gtime end_time_L3 = allpars[idx_amb_L3].end;
	//			end_time = end_time.diff(end_time_L3) < 0 ? end_time_L3 : end_time;
	//		}

	//		if (end_time <= epoch) end_time = epoch;

	//		if(isFinalUpL1 != isFinalUpL2 || isFinalUpL1 != isFinalUpL3)  throw runtime_error("L1 L2 L3 should update at the same time.");
	//		// get the band used
	//		for (const auto& it : _band_list[sys])
	//		{
	//			t_gpar::t_type amb_type = t_gpar::AMB_IF;

	//			     if (it.first == FREQ_1) amb_type = t_gpar::AMB_L1;
	//			else if (it.first == FREQ_2) amb_type = t_gpar::AMB_L2;
	//			else if (it.first == FREQ_3) amb_type = t_gpar::AMB_L3;
	//			else {
	//				throw std::logic_error("Unknown band!!");
	//			};

	//			// add new amb
	//			int idx_amb = allpars.getParam(rec, amb_type, sat);

	//			if (amb_type == t_gpar::AMB_L1 && !isFinalUpL1) continue;
	//			if (amb_type == t_gpar::AMB_L2 && !isFinalUpL2) continue;
	//			if (amb_type == t_gpar::AMB_L3 && !isFinalUpL3) continue;

	//			int idx_max = allpars[allpars.parNumber() - 1].index;
	//			t_gpar par_newamb = t_gpar(rec, amb_type, idx_max + 1, sat);

	//			par_newamb.value(0.0);
	//			par_newamb.setTime(epoch, end_time);
	//			par_newamb.apriori(9000.0);

	//			if (idx_amb < 0)
	//			{
	//				//add_parameter(par_newamb);
	//				update_info.add(par_newamb);
	//			}
	//			// reset the new amb
	//			else
	//			{
	//				allpars[idx_amb].end = epoch - _intv;
	//				update_info.add(idx_amb + 1, par_newamb);
	//			}
	//		}
	//	}

	//	// remove old amb
	//	for (unsigned int i = 0; i < allpars.parNumber(); i++)
	//	{
	//		if (allpars[i].parType == t_gpar::AMB_L1 ||
	//			allpars[i].parType == t_gpar::AMB_L2 ||
	//			allpars[i].parType == t_gpar::AMB_L3 ||
	//			allpars[i].parType == t_gpar::AMB_L4 ||
	//			allpars[i].parType == t_gpar::AMB_L5)
	//		{
	//			if (allpars[i].end < epoch)
	//			{
	//				if (!update_info.exist(i + 1))
	//				{
	//					update_info.add(i + 1);
	//				}
	//			}
	//			else
	//			{
	//				continue;
	//			}
	//		}
	//	}
	//}

	//void t_gupdateparALL::_update_amb_pars_double_freq(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo & update_info)
	//{
	//	map<string, set<string> > mapPrn;
	//	// add new amb
	//	for (auto data : obsdata)
	//	{
	//		string rec = data.site();
	//		string sat = data.sat();
	//		GSYS   sys = data.gsys();
	//		mapPrn[rec].insert(sat);

	//		auto tb_slip = dynamic_cast<t_gturboedit*>(_cycleslip.get());
	//		int idx_amb_L1 = allpars.getParam(rec, t_gpar::AMB_L1, sat);
	//		if (idx_amb_L1 < 0 && _cycleslip->num_of_amb_arc(rec, sat, epoch) ||
	//			idx_amb_L1 >= 0 && _cycleslip->cycle_slip(data, epoch))
	//		{
	//			int  isFlag = _cycleslip->num_of_amb_arc(rec, sat, epoch);

	//			// set flag should put outside loop
	//			_cycleslip->set_amb_flag(rec, sat, isFlag);
	//			t_gtime end_time = tb_slip ? tb_slip->get_crt_amb_end(rec, sat) : LAST_TIME;

	//			// get the band used
	//			for (const auto& it : _band_list[sys])
	//			{
	//				t_gpar::t_type amb_type = t_gpar::AMB_IF;
	//				     if (it.first == FREQ_1) amb_type = t_gpar::AMB_L1;
	//				else if (it.first == FREQ_2) amb_type = t_gpar::AMB_L2;
	//				else if (it.first == FREQ_3) amb_type = t_gpar::AMB_L3;
	//				else if (it.first == FREQ_4) amb_type = t_gpar::AMB_L4;
	//				else if (it.first == FREQ_5) amb_type = t_gpar::AMB_L5;
	//				else {
	//					throw std::logic_error("Unknown band!!");
	//				};
	//				//if (amb_type == t_gpar::AMB_IF) return false;

	//				int idx_max = allpars[allpars.parNumber() - 1].index;
	//				t_gpar par_newamb = t_gpar(rec, amb_type, idx_max + 1, sat);

	//				par_newamb.value(0.0);
	//				par_newamb.setTime(epoch, end_time);
	//				par_newamb.apriori(9000.0);

	//				// add new amb
	//				int idx_amb = allpars.getParam(rec, amb_type, sat);
	//				if (idx_amb < 0)
	//				{
	//					//add_parameter(par_newamb);
	//					update_info.add(par_newamb);
	//				}
	//				// reset the new amb
	//				else
	//				{
	//					allpars[idx_amb].end = epoch - _intv;
	//					update_info.add(idx_amb + 1, par_newamb);
	//				}
	//			}
	//		}
	//	}

	//	// remove old amb
	//	for (unsigned int i = 0; i < allpars.parNumber(); i++)
	//	{
	//		if (allpars[i].parType == t_gpar::AMB_L1 ||
	//			allpars[i].parType == t_gpar::AMB_L2 ||
	//			allpars[i].parType == t_gpar::AMB_L3 ||
	//			allpars[i].parType == t_gpar::AMB_L4 ||
	//			allpars[i].parType == t_gpar::AMB_L5)
	//		{
	//			if (allpars[i].end < epoch)
	//			{
	//				if (!update_info.exist(i + 1))
	//				{
	//					update_info.add(i + 1);
	//				}
	//			}
	//			else
	//			{
	//				continue;
	//			}
	//		}
	//	}
	//}
 //   // ===============================================================================================================================================
	
    void t_gupdateparALL::_update_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo & update_info)
	{
		// add new amb
		auto tb_slip_12 = dynamic_cast<t_gturboedit*>(_cycleslip.get());
		auto tb_slip_13 = dynamic_cast<t_gturboedit*>(_cycleslip_13.get());
		auto tb_slip_14 = dynamic_cast<t_gturboedit*>(_cycleslip_14.get());
		auto tb_slip_15 = dynamic_cast<t_gturboedit*>(_cycleslip_15.get());

		map<string, set<string> > mapPrn;

		//	if (!tb_slip_12 || !tb_slip_13) throw runtime_error("tb is empty!");
		if (!tb_slip_12) throw runtime_error("tb is empty!");
		for (auto data : obsdata)
		{
			string rec = data.site();
			string sat = data.sat();
			GSYS   sys = data.gsys();

			mapPrn[rec].insert(sat);
			bool isUpL1 = false; //update L1
			bool isUpL2 = false; //update L2
			bool isUpL3 = false; //update L3
			bool isUpL4 = false;
			bool isUpL5 = false;

			int idx_amb_L2 = allpars.getParam(rec, t_gpar::AMB_L2, sat);
			int idx_amb_L3 = allpars.getParam(rec, t_gpar::AMB_L3, sat);
			int idx_amb_L4 = allpars.getParam(rec, t_gpar::AMB_L4, sat);
			int idx_amb_L5 = allpars.getParam(rec, t_gpar::AMB_L5, sat);


			int Flag12 = 0, Flag13 = 0, Flag14 = 0, Flag15 = 0;
			                Flag12 = tb_slip_12->num_of_amb_arc(rec, sat, epoch);
			if (tb_slip_13) Flag13 = tb_slip_13->num_of_amb_arc(rec, sat, epoch);
			if (tb_slip_14) Flag14 = tb_slip_14->num_of_amb_arc(rec, sat, epoch);
			if (tb_slip_15) Flag15 = tb_slip_15->num_of_amb_arc(rec, sat, epoch);

			t_gtime end_time    = LAST_TIME;
			t_gtime end_time_12 = LAST_TIME;
			t_gtime end_time_13 = LAST_TIME;
			t_gtime end_time_14 = LAST_TIME;
			t_gtime end_time_15 = LAST_TIME;

			if (idx_amb_L2 < 0 || allpars[idx_amb_L2].end < epoch)
			{
				if (Flag12 > 0)
				{
					tb_slip_12->set_amb_flag(rec, sat, Flag12);
					end_time_12 = tb_slip_12 ? tb_slip_12->get_crt_amb_end(rec, sat) : LAST_TIME;
					isUpL1 = true;
					isUpL2 = true;
				}
			}
			else
			{
				bool Slip_12 = tb_slip_12->cycle_slip(data, epoch);
				if (Slip_12)
				{
					tb_slip_12->set_amb_flag(rec, sat, Flag12);
					end_time_12 = tb_slip_12 ? tb_slip_12->get_crt_amb_end(rec, sat) : LAST_TIME;
					isUpL1 = true;
					isUpL2 = true;
					if (Flag12 == 0) throw runtime_error("Flag12 is zero while idx_L2 >= 0");
				}
			}

			if (tb_slip_13 || tb_slip_14 || tb_slip_15)
			{
				if (idx_amb_L3 < 0)
				{
					if (Flag13 > 0)
					{
						tb_slip_13->set_amb_flag(rec, sat, Flag13);
						end_time_13 = tb_slip_13 ? tb_slip_13->get_crt_amb_end(rec, sat) : LAST_TIME;
						isUpL1 = true;
						isUpL3 = true;
					}
				}
				else
				{
					bool Slip_13 = tb_slip_13->cycle_slip(data, epoch);
					if (Slip_13)
					{
						tb_slip_13->set_amb_flag(rec, sat, Flag13);
						end_time_13 = tb_slip_13 ? tb_slip_13->get_crt_amb_end(rec, sat) : LAST_TIME;
						isUpL1 = true;
						isUpL3 = true;
						if (Flag13 == 0) throw runtime_error("Flag13 is zero while idx_L3 >= 0");
					}
				}

				if (idx_amb_L4 < 0)
				{
					if (Flag14 > 0)
					{
						tb_slip_14->set_amb_flag(rec, sat, Flag14);
						end_time_14 = tb_slip_14 ? tb_slip_14->get_crt_amb_end(rec, sat) : LAST_TIME;
						isUpL1 = true;
						isUpL4 = true;
					}
				}
				else
				{
					bool Slip_14 = tb_slip_14->cycle_slip(data, epoch);
					if (Slip_14)
					{

						tb_slip_14->set_amb_flag(rec, sat, Flag14);
						end_time_14 = tb_slip_14 ? tb_slip_14->get_crt_amb_end(rec, sat) : LAST_TIME;
						isUpL1 = true;
						isUpL4 = true;
						if (Flag14 == 0) throw runtime_error("Flag14 is zero while idx_L3 >= 0");
					}

				}				
				
				if (idx_amb_L5 < 0)
				{
					if (Flag15 > 0)
					{
						tb_slip_15->set_amb_flag(rec, sat, Flag15);
						end_time_15 = tb_slip_15 ? tb_slip_15->get_crt_amb_end(rec, sat) : LAST_TIME;
						isUpL1 = true;
						isUpL5 = true;
					}
				}
				else
				{
					bool Slip_15 = tb_slip_15->cycle_slip(data, epoch);
					if (Slip_15)
					{
						tb_slip_15->set_amb_flag(rec, sat, Flag15);
						end_time_15 = tb_slip_15 ? tb_slip_15->get_crt_amb_end(rec, sat) : LAST_TIME;
						isUpL1 = true;
						isUpL5 = true;
						if (Flag15 == 0) throw runtime_error("Flag12 is zero while idx_L5 >= 0");
					}

				}
				if (!isUpL2 && !isUpL3 && !isUpL4 && !isUpL5) continue;      // 无需更新
				end_time = end_time.diff(end_time_12) < 0 ? end_time : end_time_12;
				end_time = end_time.diff(end_time_13) < 0 ? end_time : end_time_13;
				end_time = end_time.diff(end_time_14) < 0 ? end_time : end_time_14;
				end_time = end_time.diff(end_time_15) < 0 ? end_time : end_time_15;
			}
			else
			{
				if (!isUpL2) continue;      // 无需更新
				end_time = end_time_12;
			}


			// get the band used
			for (const auto& it : _band_list[sys])
			{

				t_gpar::t_type amb_type = t_gpar::AMB_IF;


				if (it.first == FREQ_1) amb_type = t_gpar::AMB_L1;
				else if (it.first == FREQ_2) amb_type = t_gpar::AMB_L2;
				else if (it.first == FREQ_3) amb_type = t_gpar::AMB_L3;
				else if (it.first == FREQ_4) amb_type = t_gpar::AMB_L4;
				else if (it.first == FREQ_5) amb_type = t_gpar::AMB_L5;
				else {
					throw std::logic_error("Unknown band!!");
				};


				// add new amb
				int idx_amb = allpars.getParam(rec, amb_type, sat);


				if (amb_type == t_gpar::AMB_L2 && Flag12 == 0) continue;
				if (amb_type == t_gpar::AMB_L3 && Flag13 == 0) continue;
				if (amb_type == t_gpar::AMB_L4 && Flag14 == 0) continue;
				if (amb_type == t_gpar::AMB_L5 && Flag15 == 0) continue;


				int idx_max = allpars[allpars.parNumber() - 1].index;
				t_gpar par_newamb = t_gpar(rec, amb_type, idx_max + 1, sat);


				par_newamb.value(0.0);
				par_newamb.setTime(epoch, end_time);
				par_newamb.apriori(9000.0);


				if (idx_amb < 0)
				{
					update_info.add(par_newamb);
				}
				// reset the new amb
				else
				{
					allpars[idx_amb].end = epoch - _intv;
					update_info.add(idx_amb + 1, par_newamb);
				}
			}
		}

		// remove old amb
		for (unsigned int i = 0; i < allpars.parNumber(); i++)
		{
			if (allpars[i].parType == t_gpar::AMB_L1 ||
				allpars[i].parType == t_gpar::AMB_L2 ||
				allpars[i].parType == t_gpar::AMB_L3 ||
				allpars[i].parType == t_gpar::AMB_L4 ||
				allpars[i].parType == t_gpar::AMB_L5)
			{
				if (allpars[i].end < epoch)
				{
					if (!update_info.exist(i + 1))
					{
						update_info.add(i + 1);
					}
				}
				else if (mapPrn[allpars[i].site].count(allpars[i].prn) == 0)
				{   
					if (!update_info.exist(i + 1))
					{
						allpars[i].end = epoch - _intv; // glfeng 2018.5.1
						update_info.add(i + 1);
					}
				}
				else
				{
					//*******************************
			//   if (allpars[i].beg + -1*_intv < epoch)
			//   {
			//   	t_gpar newamb = allpars[i];
			//   	newamb.beg = epoch;
			//   	if (!update_info.exist(i + 1))
			//   	{
			//   		update_info.add(i + 1, newamb);
			//   	}
			//   }

					if (allpars[i].amb_ini == true && allpars[i].beg + -1 * _intv < epoch)
					{
						t_gpar newamb = allpars[i];
						newamb.beg = epoch;
						if (!update_info.exist(i + 1))
						{
							update_info.add(i + 1, newamb);
						}
					}
					//*****************************
					continue;
				}
			}
		}
	}

	void t_gupdateparALL::_lite_update_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info)
	{
		// add/reset amb
		for (auto data : obsdata)
		{
			int num_band = 0;
			string rec = data.site();
			string sat = data.sat();
			GSYS   sys = data.gsys();

			for (const auto& it : _band_list[sys])
			{
				t_gpar::t_type amb_type = t_gpar::AMB_IF;
				if (it.first == FREQ_1)      amb_type = t_gpar::AMB_L1;
				else if (it.first == FREQ_2) amb_type = t_gpar::AMB_L2;
				else if (it.first == FREQ_3) amb_type = t_gpar::AMB_L3;
				else if (it.first == FREQ_4) amb_type = t_gpar::AMB_L4;
				else if (it.first == FREQ_5) amb_type = t_gpar::AMB_L5;
				else {
					throw std::logic_error("Unknown band!!");
				};
				num_band++;
				if (num_band > _freq) continue;

				GOBS obsL = data.select_phase(it.second, true);

				int idx_amb = allpars.getParam(rec, amb_type, sat);
				if (idx_amb < 0)
				{
					if (obsL != GOBS::X)
					{
						int idx_max = allpars[allpars.parNumber() - 1].index;
						t_gpar par_newamb = t_gpar(rec, amb_type, idx_max + 1, sat);

						par_newamb.value(0.0);
						par_newamb.setTime(epoch, epoch);
						par_newamb.apriori(9000.0);
						update_info.add(par_newamb);
					}
				}
				else
				{
					if (obsL != GOBS::X)
					{
						if (data.getlli(obsL) < 1)
						{
							allpars[idx_amb].end = epoch;
						}
						else
						{
							allpars[idx_amb].end = epoch - _intv;

							int idx_max = allpars[allpars.parNumber() - 1].index;
							t_gpar par_newamb = t_gpar(rec, amb_type, idx_max + 1, sat);

							par_newamb.value(0.0);
							par_newamb.setTime(epoch, epoch);
							par_newamb.apriori(9000.0);
							update_info.add(idx_amb + 1, par_newamb);
						}
					}
					else  // obsL==X
					{
						update_info.add(idx_amb + 1);
					}
				}
			}
		}

		// remove old amb
		for (unsigned int i = 0; i < allpars.parNumber(); i++)
		{
			if (allpars[i].parType == t_gpar::AMB_L1 || allpars[i].parType == t_gpar::AMB_L2 ||
				allpars[i].parType == t_gpar::AMB_L3 || allpars[i].parType == t_gpar::AMB_L4 || allpars[i].parType == t_gpar::AMB_L5)
			{
				if (allpars[i].end < epoch)
				{
					if (!update_info.exist(i + 1))
					{
						update_info.add(i + 1);
					}
				}
				else if (allpars[i].amb_ini == true && allpars[i].beg + -1 * _intv < epoch)
				{
					t_gpar newamb = allpars[i];
					newamb.beg = epoch;
					if (!update_info.exist(i + 1))
					{
						update_info.add(i + 1, newamb);
					}
				}
			}
		}
	}

	void t_gupdateparALL::_udpate_ion_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo & update_info)
	{
		map<string, set<string> > mapPrn;
		for (auto& data : obsdata)
		{
			mapPrn[data.site()].insert(data.sat());

			// add new sion par
			int idx_ion = allpars.getParam(data.site(), t_gpar::SION, data.sat());
			if (idx_ion < 0)
			{
				int idx_max = allpars[allpars.parNumber() - 1].index;
				t_gpar par_ion = t_gpar(data.site(), t_gpar::SION, idx_max + 1, data.sat());
				par_ion.value(0.0);
				par_ion.setTime(epoch, epoch);
				par_ion.apriori(_sig_ion);
				update_info.add(par_ion);
			}

			if (_iono_order == IONO_ORDER::THIRD_ORDER)
			{
				int idx_ion_3 = allpars.getParam(data.site(), t_gpar::SION_3, data.sat());
				if (idx_ion_3 < 0)
				{
					int idx_max = allpars[allpars.parNumber() - 1].index;
					t_gpar par_ion = t_gpar(data.site(), t_gpar::SION_3, idx_max + 1, data.sat());
					par_ion.value(0.0);
					par_ion.setTime(epoch, epoch);
					par_ion.apriori(9000);
					update_info.add(par_ion);
				}
			}
		}

		// remove old sion par
		for (unsigned int i = 0; i < allpars.parNumber(); i++)
		{
			if (allpars[i].parType == t_gpar::SION && mapPrn[allpars[i].site].count(allpars[i].prn) == 0)
			{
				update_info.add(i + 1);
			}

			if (_iono_order == IONO_ORDER::THIRD_ORDER)
			{
				if (allpars[i].parType == t_gpar::SION_3 && mapPrn[allpars[i].site].count(allpars[i].prn) == 0)
				{
					update_info.add(i + 1);
				}
			}
		}

	}

	void t_gupdateparALL::_check_log(shared_ptr<t_gcycleslip> cycleslip_12, shared_ptr<t_gcycleslip> cycleslip_13)
	{
	}
}
