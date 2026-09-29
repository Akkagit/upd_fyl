/**
*
* @verbatim
History
-1.0 BoWong  2019-02-25  creat the file.
-1.1 BoWong  2019-04-08  Adding Doxygen Style Code Remarks
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambflag.cpp
* @brief	Storage the XXXXddd0.yyo.log/log13 files' data(only one site)
*				XXXX  ---- SITE name
*				 ddd  ---- Doy of the file
*				  yy  ---- year
*
* @author   BoWong, Wuhan University
* @version	1.0.0
* @date		2019-04-08
*
*/

#include"gdata/gambflag.h"

using namespace std;
namespace great {

	/** @brief default constructor. */
	t_gambflag::t_gambflag() //:_ambflag_data()
	{
		gtrace("t_gambflag::constructor");
		id_type(t_gdata::AMBFLAG);
	}

	/** @brief default destructor. */
	t_gambflag::~t_gambflag()
	{
		gtrace("t_gambflag::destructor");
	}

	void t_gambflag::setAmbFlagHead(ambflag_hd head)
	{
		gtrace("t_gambflag::set_ambflag_head");
		_ambflag_head = make_shared<ambflag_hd>(head);
	}

	shared_ptr<ambflag_hd>& t_gambflag::getAmbFlagHead()
	{
		return _ambflag_head;
	}

	void t_gambflag::addAmbFlagData(string prn, const ambflag_data& data)
	{
		gtrace("t_gambflag::add_ambflag_data");
		try	
		{
			
			_gmutex.lock();
			_all_sats_ambflag[prn].push_back(make_shared<ambflag_data>(data));
			_gmutex.unlock();
		}
		catch (exception ex)
		{
			std::cout << ex.what();
		}
	}

	t_map_sat_ambflag& t_gambflag::getAmbFlagData()
	{
		return _all_sats_ambflag;
	}

	t_vec_ambflag& t_gambflag::getSatAmbFlag(const string& sat)
	{
		return _all_sats_ambflag[sat];
	}

	bool t_gambflag::isValid(const string& prn, const t_gtime& time, int& pos)
	{
		try
		{
			pos = 0;
			if (_all_sats_ambflag.find(prn) == _all_sats_ambflag.end())
				return false;
			
			bool is_found = false;
			for (auto iter = _all_sats_ambflag[prn].begin(); iter != _all_sats_ambflag[prn].end(); ++iter) {
				if (time >= epoch2time((*iter)->beg_epo) && time <= epoch2time((*iter)->end_epo)) {
					if ((*iter)->identify == "BAD" || (*iter)->identify == "DEL") { 
						return false; 
					}
					else { 
						is_found = true; 
						pos = distance(_all_sats_ambflag[prn].begin(), iter); break; 
					}
				}
			}
			return is_found;
			//t_gtime beg;
			//t_gtime end;

			// find in unused
			//for (auto iter = _all_sats_ambflag[prn].begin(); iter != _all_sats_ambflag[prn].end(); iter++)
			//{
			//	if ((*iter)->identify == "BAD" || (*iter)->identify == "DEL") {
			//		//beg.from_mjd(_ambflag_head->beg_mjd, _ambflag_head->beg_sod + _ambflag_head->intv*((*iter)->beg_epo - 1), 0.0);
			//		//end.from_mjd(_ambflag_head->beg_mjd, _ambflag_head->beg_sod + _ambflag_head->intv*((*iter)->end_epo - 1), 0.0);
			//		if (time >= epoch2time((*iter)->beg_epo) && time <= epoch2time((*iter)->end_epo)) {
			//			return false;
			//		}
			//	}
			//}

			//// find in used
			//for (auto iter = _all_sats_ambflag[prn].begin(); iter != _all_sats_ambflag[prn].end(); iter++)
			//{
			//	if ((*iter)->identify == "AMB") {
			//		//beg.from_mjd(_ambflag_head->beg_mjd, _ambflag_head->beg_sod + _ambflag_head->intv*((*iter)->beg_epo - 1), 0.0);
			//		//end.from_mjd(_ambflag_head->beg_mjd, _ambflag_head->beg_sod + _ambflag_head->intv*((*iter)->end_epo - 1), 0.0);
			//		if (time >= epoch2time((*iter)->beg_epo) && time <= epoch2time((*iter)->end_epo)) {
			//			pos = distance(_all_sats_ambflag[prn].begin(), iter);
			//			return true;
			//		}
			//	}
			//}
			//return false;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(2, "t_gambflag::is_available", "ERROR: unknown mistake");
			}
			else
			{
				cout << "ERROR: t_gambflag::is_available , unknown mistake" << endl;
			}
			return false;
			throw(-1);
		}
	}

	bool t_gambflag::is_carrier_range(const string& prn, const t_gtime& time)
	{
		if (_all_sats_ambflag.find(prn) == _all_sats_ambflag.end())
			return false;

		for (auto iter = _all_sats_ambflag[prn].begin(); iter != _all_sats_ambflag[prn].end(); ++iter) {
			if (time >= epoch2time((*iter)->beg_epo) && time <= epoch2time((*iter)->end_epo)) {
				return (*iter)->identify == "IAM" ? true: false;
			}
		}
		return false;
	}

	bool t_gambflag::is_carrier_range(const string& prn, const t_gtime& time, double &c1, double &c2)
	{
		c1 = 0;
		c2 = 0;
		if (_all_sats_ambflag.find(prn) == _all_sats_ambflag.end())
			return false;

		for (auto iter = _all_sats_ambflag[prn].begin(); iter != _all_sats_ambflag[prn].end(); ++iter) {
			if (time >= epoch2time((*iter)->beg_epo) && time <= epoch2time((*iter)->end_epo)) {
				if ((*iter)->identify != "IAM") {
					return false;
				}
				else {
					c1 = (*iter)->C1;
					c2 = (*iter)->C2;
					return true;
				}
			}
		}
		return false;
	}

	int t_gambflag::get_amb_pos(const string& prn, const t_gtime& time, bool apply_carrier_range)
	{
		if (_all_sats_ambflag.find(prn) == _all_sats_ambflag.end())
			return -1;

		int pos = -1;
		for (auto iter = _all_sats_ambflag[prn].begin(); iter != _all_sats_ambflag[prn].end(); ++iter) {
			if (time >= epoch2time((*iter)->beg_epo) && time <= epoch2time((*iter)->end_epo)) {
				if ((*iter)->identify == "BAD" || (*iter)->identify == "DEL") {
					return -1;
				}
				else if ((*iter)->identify == "IAM") {
					if (!apply_carrier_range) {
						pos = distance(_all_sats_ambflag[prn].begin(), iter); break;
					}
					else {
						return -1;
					}
				}
				else {
					pos = distance(_all_sats_ambflag[prn].begin(), iter); break;
				}
			}
		}
		return pos;
	}

	t_gtime t_gambflag::epoch2time(const int &epo)
	{
		t_gtime time;
		time.from_mjd(_ambflag_head->beg_mjd, _ambflag_head->beg_sod + _ambflag_head->intv * (epo - 1), 0.0);
		return time;
	}

	void t_gambflag::reset_iflag(const string& prn, const string& flag, const int& pos)
	{ 
		_all_sats_ambflag[prn][pos]->iflag = flag;
	}

	set<string> t_gambflag::getAllSatSet()
	{
		set<string> satlist;
		auto it_sat = _all_sats_ambflag.begin();
		for (; it_sat != _all_sats_ambflag.end(); it_sat++)
		{
			satlist.insert(it_sat->first);
		}
		return satlist;
	}

	t_gtime t_gambflag::getMaxEndEpoch()
	{
		int max_end_epo = 0;
		for(auto sat_record : _all_sats_ambflag) {
			for (auto record : sat_record.second) {
				if (record->end_epo > max_end_epo) {
					max_end_epo = record->end_epo;
				}
			}
		}
		t_gtime beg_time(this->_ambflag_head->beg_mjd, int(this->_ambflag_head->beg_sod), this->_ambflag_head->beg_sod - int(this->_ambflag_head->beg_sod));
		t_gtime end_time = beg_time + (max_end_epo - 1) * this->_ambflag_head->intv;
		return end_time;
	}

	void t_gambflag::clearAllAmbFlagData()
	{
		_all_sats_ambflag.clear();
	}

}//namespace