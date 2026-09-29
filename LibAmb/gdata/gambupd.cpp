/**
*
* @verbatim
History
-1.0 BoWong  2019-02-20  creat the file.
-1.1 BoWong  2019-04-08  Adding Doxygen Style Code Remarks
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambupd.cpp
* @brief	Storage the ambiguity upd files' data(only one site) for narrow-lane upd compute
*
* @author   BoWong, Wuhan University
* @version	1.0.0
* @date		2019-04-08
*
*/

#include"gdata/gambupd.h"

using namespace std;
namespace great
{
	/** @brief default constructor. */
	t_gambupd::t_gambupd() :_ambupd()
	{
		gtrace("t_gambupd::constructor");
		id_type(t_gdata::AMBUPD);
		_sat_list.clear();
	}

	/** @brief default destructor. */
	t_gambupd::~t_gambupd()
	{
		gtrace("t_gambupd::destructor");
	}

	/**
	* @brief set the station of the class data.
	* @param[in]  sta		station name.
	* @return	  void
	*/
	//void  t_gambupd::set_station(string sta)
	//{
	//	gtrace("t_gambupd::set_station");
	//	_station = sta;
	//}

	void t_gambupd::addAmbUpd(t_gtime epoch, string sat, const t_ambupd_data& data)
	{
		gtrace("t_gambupd::add_sat_ambupd");
		if (_ambupd.find(epoch) == _ambupd.end())
		{
			_ambupd[epoch][sat] = make_shared<t_ambupd_data>(data);
		}
		else if (_ambupd[epoch].find(sat) == _ambupd[epoch].end())
		{
			_ambupd[epoch][sat] = make_shared<t_ambupd_data>(data);
		}
		else
		{//LX modify the bug from <= to ==  
			if (_ambupd[epoch][sat]->amb_if == 0.0&&data.amb_if !=0) _ambupd[epoch][sat]->amb_if = data.amb_if;
			if (_ambupd[epoch][sat]->amb_wl == 0.0&&data.amb_wl != 0) _ambupd[epoch][sat]->amb_wl = data.amb_wl;
			if (_ambupd[epoch][sat]->amb_l3 == 0.0&&data.amb_l3 != 0) _ambupd[epoch][sat]->amb_l3 = data.amb_l3;
			if (_ambupd[epoch][sat]->amb_l4 == 0.0&&data.amb_l4 != 0) _ambupd[epoch][sat]->amb_l4 = data.amb_l4;
			if (_ambupd[epoch][sat]->amb_l5 == 0.0&&data.amb_l5 != 0) _ambupd[epoch][sat]->amb_l5 = data.amb_l5;
		}

		_sat_list.insert(sat);
	}

	shared_ptr<t_ambupd_data>& t_gambupd::getSatAmbUpd(t_gtime epoch, string sat)
	{
		return _ambupd[epoch][sat];
	}

	t_sat_ambupd& t_gambupd::getEpoAmbUpd(t_gtime epoch)
	{
		return _ambupd[epoch];
	}

	t_satepo_ambupd t_gambupd::getSatepoAmbUpd(t_gtime& beg, t_gtime& end, string sat)
	{
		t_satepo_ambupd satepo_ambupd;
		for (auto iter = _ambupd.upper_bound(beg); iter != _ambupd.lower_bound(end); ++iter) {
			auto iter_sat = iter->second.find(sat);
			if (iter_sat != iter->second.end())
				satepo_ambupd[iter->first] = iter_sat->second;
		}
		return satepo_ambupd;
	}

	t_epo_ambupd& t_gambupd::getAllAmbUpd()
	{
		return _ambupd;
	}
    
	t_gtime t_gambupd::getSatEndTime()
	{
		t_gtime end_time;
		end_time = _ambupd.rbegin()->first;		
		return end_time;
	}

}