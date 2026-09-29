/**
*
* @verbatim
	History
	 -1.0 Jie Li		2019-01-11 creat the file.
	 -1.1 Jie Li     2019-04-01 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gallsatparam.cpp
* @brief		Storage the sat parameters data
*
* @author		Jie Li, Wuhan University
* @version		1.0.0
* @date			2019-01-11
*
*/



#include "gall/gallsatparam.h"

using namespace std;
using namespace gnut;

namespace great
{
	/** @brief constructor. */
	t_gallsatinfo::t_gallsatinfo()
	{
		gtrace("t_gallsatparam::constructor");
		id_type(t_gdata::SATPARS);
	}

	/** @brief  destructor. */
	t_gallsatinfo:: ~t_gallsatinfo()
	{

	}

	/**
	* @brief push one satellite parameters data into the data map.
	* @param[in]  prn		satellite name.
	* @param[in]  svn		satellite svn.
	* @param[in]  start		launched time of the satellite.
	* @param[in]  end		decommissioned time of the satellite.
	* @param[in]  cosparid	cospar_id of the satellite.
	* @param[in]  mass		mass of the satellite.
	* @param[in]  maxyaw	max_yaw of the satellite.
	* @param[in]  fid		frequence id of the satellite.
	* @param[in]  blocktype	BLOCK-TYPE of the satellite.
	* @return	whether exit normally
	*/
	int t_gallsatinfo::addData(string prn, string svn, t_gtime start, t_gtime end, string cosparid, double mass, double maxyaw, int fid, string blocktype)
	{
		// jdhuang, 有内存泄漏危险，导致读文件崩溃后被发现，并修正
		t_gsatparam tmp;
		tmp.set_prn(prn);
		tmp.set_svn(svn);
		tmp.set_start(start);
		tmp.set_end(end);
		tmp.set_cosparid(cosparid);
		tmp.set_mass(mass);
		tmp.set_maxyaw(maxyaw);
		tmp.set_fid(fid);
		tmp.set_blocktype(blocktype);

		_allsatinfo[prn][svn] = tmp;
		//_allsatinfo[prn][svn].set_prn(prn);
		//_allsatinfo[prn][svn].set_svn(svn);
		//_allsatinfo[prn][svn].set_start(start);
		//_allsatinfo[prn][svn].set_end(end);
		//_allsatinfo[prn][svn].set_cosparid(cosparid);
		//_allsatinfo[prn][svn].set_mass(mass);
		//_allsatinfo[prn][svn].set_maxyaw(maxyaw);
		//_allsatinfo[prn][svn].set_fid(fid);
		//_allsatinfo[prn][svn].set_blocktype(blocktype);
		return 0;
	}

	/** @brief whether the satellite parameters data is empty.
	* @return  bool
	*	@retval   true    satellite parameters data is empty
	*   @retval   false   satellite parameters data is existent
	*/
	bool t_gallsatinfo::isEmpty()
	{
		if (_allsatinfo.size() == 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	/**
	* @brief get valid satellite parameters data.
	* @param[in]  sat			satellite name.
	* @param[in]  ics_start		begin time from ics data.
	* @param[in]  ics_end		end time from ics data.
	* @param[out] sat_param		corresponding satellite parameters data.
	* @return	if find the correct satellite parameters data.
	*	 @retval true	get data successfully
	*	 @retval false	get data failed
	*/
	bool t_gallsatinfo::getSatinfoUsed(string sat, t_gtime ics_start, t_gtime ics_end, t_gsatparam& sat_param)
	{
		if (_allsatinfo.count(sat) == 0)
		{
			cout << "ERROR: " << sat << "is not existed in sat_paramters file" << endl;
			return false;
		}

		map<string, t_gsatparam>::iterator itSvn = _allsatinfo[sat].begin();
		for (itSvn; itSvn != _allsatinfo[sat].end(); ++itSvn)
		{
			t_gtime satpars_st = (itSvn->second).start();
			t_gtime satpars_et = (itSvn->second).end();
			double delt_st = satpars_st.dmjd() - ics_start.dmjd();
			double delt_et = satpars_et.dmjd() - ics_end.dmjd();
			double eps = 5.0 / 86400.0;
			if (delt_st > eps || delt_et < -eps)
			{
				continue;
			}
			itSvn->second.valid = true;
			sat_param = itSvn->second;
			return true;
		}
		if (_log)
		{
			_log->comment(0, "ERROR : ics time is out of range,can't find corresponding sat information ");
		}
		else
		{
			cerr << "ERROR : ics time is out of range,can't find corresponding sat information " << endl;
		}
		return false;
	}
	int t_gallsatinfo::getGloFid(string sat,t_gtime epoch)
	{
		t_gsatparam tmp_satinfo;
		if (getSatinfoUsed(sat, epoch, epoch, tmp_satinfo))
		{
			return tmp_satinfo.fid();
		}
		return 0;
	}
}