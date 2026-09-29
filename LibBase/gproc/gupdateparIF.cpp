/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2019-12-23  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gupdatepar.cpp
* @brief		define update parameter IF model
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-12-23
*
*/

#include "gproc/gupdateparIF.h"
#include "gutils/gturboedit.h"

namespace great
{
	t_gupdateparIF::t_gupdateparIF(shared_ptr<t_gcycleslip> cycleslip, const map< GSYS, map<FREQ_SEQ, GOBSBAND> >& band_list)
	{
		_cycleslip = cycleslip;
		_band_list = band_list;
	}

	//t_gupdateparIF::t_gupdateparIF(shared_ptr<t_gcycleslip> cycleslip, t_gsetbase* set)
	//{
	//	_cycleslip = cycleslip;
	//	_set = set;
	//}

	t_gupdateparIF::~t_gupdateparIF()
	{
	}

	void t_gupdateparIF::_update_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata,t_gupdateparinfo& update_info)
	{
		//bool realtime = dynamic_cast<t_gsetproc*>(_set)->realtime();
		if (_lite_update_amb)
		{
			_update_one_type_amb_pars(epoch, allpars, obsdata, t_gpar::AMB_IF, update_info);
		}
		else
		{
			_update_one_type_amb_pars(epoch, allpars, obsdata, t_gpar::AMB_IF, _cycleslip.get(), update_info);
		}

	}

	void t_gupdateparIF::_update_one_type_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata,t_gpar::t_type ambtype, t_gcycleslip* _slip, t_gupdateparinfo& update_info)
	{

		map<string, set<string> > mapPrn;
		for (auto obs_data : obsdata)
		{
			string rec = obs_data.site();
			string sat = obs_data.sat();
			mapPrn[rec].insert(sat);


			int idx_amb = allpars.getParam(rec, ambtype, sat);

			// judge cycle  and else
			t_gtime crt_epoch = epoch;
			// first add new amb
			// second update old amb if cycleslip
			if (idx_amb < 0 && _slip->num_of_amb_arc(rec, sat, crt_epoch) ||
				idx_amb >= 0 && _slip->cycle_slip(obs_data, crt_epoch))
			{
				int idx_max = allpars[allpars.parNumber() - 1].index;
				t_gpar par_newamb = t_gpar(rec, ambtype, idx_max + 1, sat);

				int amb_flag = _slip->num_of_amb_arc(rec, sat, epoch);
				_slip->set_amb_flag(rec, sat, amb_flag);
				auto tb_slip = dynamic_cast<t_gturboedit*>(_slip);

				t_gtime end_time = LAST_TIME;
				if (tb_slip) {
					end_time = tb_slip->get_crt_amb_end(rec, sat);
				}
				par_newamb.value(0.0);
				par_newamb.setTime(epoch, end_time);
				par_newamb.apriori(10000.0);

				// add new amb
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

		// delete old amb
		for (unsigned int i = 0; i < allpars.parNumber(); i++)
		{
			if (allpars[i].parType == ambtype)
			{
				// jdhuang : 都以log文件为准
				if (allpars[i].end < epoch)
				{
					if (!update_info.exist(i + 1))
					{
						//allpars[i].end = epoch - _intv; // glfeng 2018.5.1
						update_info.add(i + 1);
					}
				}
				//// jdhuang
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
					//****每个一段时间初始化模糊度
					//if (allpars[i].beg + -1*_intv < epoch)
					//{
					//	t_gpar newamb = allpars[i];
					//	newamb.beg = epoch;
					//	if (!update_info.exist(i + 1))
					//	{
					//		update_info.add(i + 1, newamb);
					//	}
					//}
					//*****************************
					continue;
				}
			}
		}

	}

	void t_gupdateparIF::_update_one_type_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gpar::t_type ambtype, t_gupdateparinfo& update_info)
	{
        //// origin ======================================================================
		//map<string, set<string> > mapPrn;
		//for (auto obs_data : obsdata)
		//{
		//	string rec = obs_data.site();
		//	string sat = obs_data.sat();
		//	mapPrn[rec].insert(sat);


		//	int idx_amb = allpars.getParam(rec, ambtype, sat);

		//	// judge cycle  and else
		//	t_gtime crt_epoch = epoch;
		//	map<GOBS, int> slip = obs_data.slip();
		//	bool cycleslip = false;
		//	for (auto it = slip.begin(); it != slip.end(); it++)
		//	{
		//		if (it->second != 0)
		//		{
		//			cycleslip = true;
		//			break;
		//		}
		//	}
		//	map<GOBS, int> lli = obs_data.lli();
		//	for (auto it = lli.begin(); it != lli.end(); it++)
		//	{
		//		if (it->second != 0)
		//		{
		//			cycleslip = true;
		//			break;
		//		}
		//	}

		//	// first add new amb
		//	// second update old amb if cycleslip
		//	if (idx_amb < 0 || idx_amb >= 0 && cycleslip)
		//	{
		//		int idx_max = allpars[allpars.parNumber() - 1].index;
		//		t_gpar par_newamb = t_gpar(rec, ambtype, idx_max + 1, sat);

		//		//int amb_flag = _slip->num_of_amb_arc(rec, sat, epoch);
		//		//_slip->set_amb_flag(rec, sat, amb_flag);
		//		//auto tb_slip = dynamic_cast<t_gturboedit*>(_slip);

		//		t_gtime end_time = LAST_TIME;
		//		//if (tb_slip) {
		//		//	end_time = tb_slip->get_crt_amb_end(rec, sat);
		//		//}
		//		par_newamb.value(0.0);
		//		par_newamb.setTime(epoch, end_time);
		//		par_newamb.apriori(10000.0);

		//		// add new amb
		//		if (idx_amb < 0)
		//		{
		//			update_info.add(par_newamb);
		//		}
		//		// reset the new amb
		//		else
		//		{
		//			allpars[idx_amb].end = epoch - _intv;
		//			update_info.add(idx_amb + 1, par_newamb);
		//		}
		//	}
		//}

		//// delete old amb
		//for (unsigned int i = 0; i < allpars.parNumber(); i++)
		//{
		//	if (allpars[i].parType == ambtype)
		//	{
		//		// jdhuang : 都以log文件为准
		//		if (allpars[i].end < epoch)
		//		{
		//			if (!update_info.exist(i + 1))
		//			{
		//				//allpars[i].end = epoch - _intv; // glfeng 2018.5.1
		//				update_info.add(i + 1);
		//			}
		//		}
		//		else
		//		{
		//			continue;
		//		}
		//	}
		//}
		//// origin ======================================================================

		for (auto obs_data : obsdata)
		{
			string rec = obs_data.site();
			string sat = obs_data.sat();
			GSYS   sys = obs_data.gsys();

			int idx_amb = allpars.getParam(rec, ambtype, sat);
			GOBS obsL1, obsL2;
			if (ambtype == t_gpar::AMB_IF)
			{
				obsL1 = obs_data.select_phase(_band_list[sys][FREQ_1]);
				obsL2 = obs_data.select_phase(_band_list[sys][FREQ_2]);
			}
			else if (ambtype == t_gpar::AMB13_IF)
			{
				obsL1 = obs_data.select_phase(_band_list[sys][FREQ_1]);
				obsL2 = obs_data.select_phase(_band_list[sys][FREQ_3]);
			}

			if (idx_amb < 0)
			{
				if (obsL1 != GOBS::X && obsL2 != GOBS::X)
				{
					int idx_max = allpars[allpars.parNumber() - 1].index;
					t_gpar par_newamb = t_gpar(rec, ambtype, idx_max + 1, sat);

					par_newamb.value(0.0);
					par_newamb.setTime(epoch, epoch);
					par_newamb.apriori(9000.0);
					update_info.add(par_newamb);
				}
			}
			else
			{
				if (obsL1 != GOBS::X && obsL2 != GOBS::X)
				{
					if (obs_data.getlli(obsL1) < 1 && obs_data.getlli(obsL2) < 1)
					{
						allpars[idx_amb].end = epoch;
					}
					else
					{
						allpars[idx_amb].end = epoch - _intv;

						int idx_max = allpars[allpars.parNumber() - 1].index;
						t_gpar par_newamb = t_gpar(rec, ambtype, idx_max + 1, sat);

						par_newamb.value(0.0);
						par_newamb.setTime(epoch, epoch);
						par_newamb.apriori(9000.0);

						update_info.add(idx_amb + 1, par_newamb);
					}
				}
				else  // obsL1==X || obsL2==X
				{
					update_info.add(idx_amb + 1);
				}
			}
		}

		// delete old amb
		for (unsigned int i = 0; i < allpars.parNumber(); i++)
		{
			if (allpars[i].parType == ambtype)
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

	t_gupdateparIF13::t_gupdateparIF13(shared_ptr<t_gcycleslip> cycleslip12, shared_ptr<t_gcycleslip> cycleslip13, const map< GSYS, map<FREQ_SEQ, GOBSBAND> >& band_list):
		t_gupdateparIF(cycleslip12, band_list)
	{
		_cycleslip13 = cycleslip13;
	}
	t_gupdateparIF13::~t_gupdateparIF13()
	{
	}

	void t_gupdateparIF13::set_cycleslip13(shared_ptr<t_gcycleslip> cycleslip13)
	{
		_cycleslip13 = cycleslip13;
	}

	t_gupdateparinfo t_gupdateparIF13::get_all_update_parameters(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata)
	{
		t_gupdateparinfo update_info;
		// update AMB
		t_gupdateparIF13::_update_amb_pars(epoch, allpars, obsdata, update_info);

		// update gps rec ifb pars
		if (allpars.orbParNumber() > 0)
		{
			_udpate_gps_rec_ifb_pars(epoch, allpars, obsdata, update_info);
		}

		// update Other pars
		_update_process_pars(epoch, allpars, update_info);

		return update_info;
	}

	void t_gupdateparIF13::_update_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata,t_gupdateparinfo& update_info)
	{
		t_gupdateparIF::_update_amb_pars(epoch, allpars, obsdata, update_info);
		if (_lite_update_amb)
		{
			_update_one_type_amb_pars(epoch, allpars, obsdata, t_gpar::AMB13_IF, update_info);
		}
		else
		{
			_update_one_type_amb_pars(epoch, allpars, obsdata, t_gpar::AMB13_IF, _cycleslip13.get(), update_info);

		}
	}

	t_updateparWL::t_updateparWL(shared_ptr<t_gcycleslip> cycleslip, const map< GSYS, map<FREQ_SEQ, GOBSBAND> >& band_list):
		t_gupdateparIF(cycleslip, band_list)
	{
	}
	t_updateparWL::~t_updateparWL()
	{
	}
	void t_updateparWL::_update_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo & update_info)
	{
		_update_one_type_amb_pars(epoch, allpars, obsdata, t_gpar::AMB_WL, _cycleslip.get(), update_info);
	}
}