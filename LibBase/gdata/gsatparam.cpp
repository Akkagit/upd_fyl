/**
 *
 * @verbatim
	 History
	  -1.0	Jie Li	2019-01-11	creat the file.
	  -1.1  Jie Li	2019-04-01	add code in doxygen style.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gsatparam.cpp
 * @brief		The base class used to storage one satellite information.
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 */

#include "gdata/gsatparam.h"

using namespace std;

namespace great
{
	/** @brief constructor. */
	t_gsatparam::t_gsatparam()
	{
		gtrace("t_gsatparam::t_gsatparam");
		id_type(t_gdata::SATPARS);
		_mass = 0.0;
	}

	/** @brief destructor. */
	t_gsatparam:: ~t_gsatparam()
	{

	}

	/**
	* @brief set satellite prn.
	* @param[in]  prn	satellite prn
	*/
	void t_gsatparam::set_prn(string prn)
	{
		_prn = prn;
	}

	/**
	* @brief set satellite svn.
	* @param[in]  svn	satellite svn
	*/
	void t_gsatparam::set_svn(string svn)
	{
		_svn = svn;
	}

	/**
	* @brief set satellite launch time.
	* @param[in]  start 	satellite launch time
	*/
	void t_gsatparam::set_start(t_gtime start)
	{
		_start = start;
	}

	/**
	* @brief set satellite decommissioned time.
	* @param[in]  end	  satellite decommissioned time
	*/
	void t_gsatparam::set_end(t_gtime end)
	{
		_end = end;
	}

	/**
	* @brief set satellite cospar_id.
	* @param[in]  cosparid	  satellite cospar_id
	*/
	void t_gsatparam::set_cosparid(string cosparid)
	{
		_cosparid = cosparid;
	}

	/**
	* @brief set satellite mass.
	* @param[in]  mass	  satellite mass
	*/
	void t_gsatparam::set_mass(double mass)
	{
		_mass = mass;
	}

	/**
	* @brief set satellite max_yaw.
	* @param[in]  maxyaw	  satellite max_yaw
	*/
	void t_gsatparam::set_maxyaw(double maxyaw)
	{
		_maxyaw = maxyaw;
	}

	/**
	* @brief set satellite frequency id.
	* @param[in]  fid	  satellite frequency id
	*/
	void t_gsatparam::set_fid(int fid)
	{
		_fid = fid;
	}

	/**
	* @brief set satellite block type.
	* @param[in]  blocktype	  satellite block type
	*/
	void t_gsatparam::set_blocktype(string blocktype)
	{
		_blocktype = blocktype;
	}

	/**
	* @brief get satellite prn.
	* @return	 satellite prn.
	*/
	string t_gsatparam::prn()
	{
		return _prn;
	}

	/**
	* @brief get satellite svn.
	* @return	 satellite svn.
	*/
	string t_gsatparam::svn()
	{
		return _svn;
	}

	/**
	* @brief get satellite launch time.
	* @return	 satellite launch time.
	*/
	t_gtime t_gsatparam::start()
	{
		return _start;
	}

	/**
	* @brief get satellite decommissioned time.
	* @return	 satellite decommissioned time.
	*/
	t_gtime t_gsatparam::end()
	{
		return _end;
	}

	/**
	* @brief get satellite cospar_id.
	* @return	 satellite cospar_id.
	*/
	string t_gsatparam::cosparid()
	{
		return _cosparid;
	}

	/**
	* @brief get satellite mass.
	* @return	 satellite mass.
	*/
	double t_gsatparam::mass()
	{
		return _mass;
	}

	/**
	* @brief get satellite max_yaw.
	* @return	 satellite max_yaw.
	*/
	double t_gsatparam::maxyaw()
	{
		return _maxyaw;
	}

	/**
	* @brief get satellite frequency id.
	* @return	 satellite frequency id.
	*/
	int t_gsatparam::fid()
	{
		return _fid;
	}

	/**
	* @brief get satellite block type.
	* @return	 satellite block type.
	*/
	string t_gsatparam::blocktype()
	{
		return _blocktype;
	}
}