/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-01-05	creat the file.
	  -1.1  Jie Li			2019-04-01	add code in doxygen style.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gics_nav.cpp
 * @brief		The base class used to storage the GNS ICS data for one satellite.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-01-05
 */


#include "gdata/gics_nav.h"
#include "gdata/gdata.h"

#include "gforce/gfsolar.h"

using namespace std;

namespace great
{
	t_gics_nav::t_gics_nav(string sat)
		: t_gdata()
	{
		sat_name = sat;
	}

	t_gics_nav::t_gics_nav()
		: t_gdata()
	{

	}

	t_gics_nav::~t_gics_nav()
	{

	}

	/**
	* @brief judge whether the sat record is empty or not.
	* @return  bool
		@retval   true    record is empty
		@retval   false   record is existent
	*/
	bool t_gics_nav::is_empty()
	{
		return _par_list.empty();
	}

	/**
	* @brief set parameter data for one parameter.
	* @param[in]   par_name   par name of string
	* @param[in]   record     par values(maybe more than one record, PWC param)
	*/
	void t_gics_nav::set_par_record(string par_name, t_ics_record & record)
	{		
		_par_list[str2ics(par_name)] = record;
	}

	/**
	* @brief get par value according to the par name(string).
	* @param[in]   par_name      par name of enum
	* @return  value of the par
	*/
	double t_gics_nav::get_par_value(string par_name)
	{
		ICS_PAR par = str2ics(par_name);
		return _par_list.count(par) > 0 ? _par_list[par].value : 0.0;
	}

	/**
	* @brief get par value according to the par name(ICS_PAR).
	* @param[in]   par_name      par name of enum
	* @return  value of the par
	*/
	double t_gics_nav::get_par_value(ICS_PAR par_name)
	{
		if (_par_list.count(par_name) > 0)
		{
			return _par_list[par_name].value;
		}
		else
		{
			if (_log)
			{
				_log->comment(1, "WARNNING : can not get the ics value of " + ics2str(par_name));
			}
			else
			{
				cout << "WARNNING : can not get the ics value of " + ics2str(par_name) << endl;
			}
			return 0.0;
		}	
	}

	/**
	* @brief get all par value of the satellite.
	* @return  par name and its'value of the satellite
	*/
	t_gics_nav::t_map_pars& t_gics_nav::get_par_list()
	{
		return _par_list;
	}

	map<ICS_PAR, double> t_gics_nav::get_ecom3_solar_pars_values()
	{
		map<ICS_PAR, double> solar_pars;
		solar_pars[ICS_PAR::ECOM_D0] = (_par_list.count(ICS_PAR::ECOM_D0) > 0) ? _par_list[ICS_PAR::ECOM_D0].value : 0;
		solar_pars[ICS_PAR::ECOM_Y0] = (_par_list.count(ICS_PAR::ECOM_Y0) > 0) ? _par_list[ICS_PAR::ECOM_Y0].value : 0;
		solar_pars[ICS_PAR::ECOM_X0] = (_par_list.count(ICS_PAR::ECOM_X0) > 0) ? _par_list[ICS_PAR::ECOM_X0].value : 0;
		return solar_pars;
	}

	/**
	* @brief get 5 solar radiation parameters value of ECOM model.
	* @return  5 solar radiation parameters value of ECOM model
	*/
	map<ICS_PAR, double> t_gics_nav::get_ecom5_solar_pars_values()
	{
		map<ICS_PAR, double> solar_pars;
		solar_pars[ICS_PAR::ECOM_D0] = (_par_list.count(ICS_PAR::ECOM_D0) > 0) ? _par_list[ICS_PAR::ECOM_D0].value : 0;
		solar_pars[ICS_PAR::ECOM_Y0] = (_par_list.count(ICS_PAR::ECOM_Y0) > 0) ? _par_list[ICS_PAR::ECOM_Y0].value : 0;
		solar_pars[ICS_PAR::ECOM_X0] = (_par_list.count(ICS_PAR::ECOM_X0) > 0) ? _par_list[ICS_PAR::ECOM_X0].value : 0;
		solar_pars[ICS_PAR::ECOM_Xc] = (_par_list.count(ICS_PAR::ECOM_Xc) > 0) ? _par_list[ICS_PAR::ECOM_Xc].value : 0;
		solar_pars[ICS_PAR::ECOM_Xs] = (_par_list.count(ICS_PAR::ECOM_Xs) > 0) ? _par_list[ICS_PAR::ECOM_Xs].value : 0;
		return solar_pars;
	}

	/**
	* @brief get 9 solar radiation parameters value of ECOM model.
	* @return  9 solar radiation parameters value of ECOM model
	*/
	map<ICS_PAR, double> t_gics_nav::get_ecom9_solar_pars_values()
	{
		map<ICS_PAR, double> solar_pars;
		solar_pars[ICS_PAR::ECOM_D0] = (_par_list.count(ICS_PAR::ECOM_D0) > 0) ? _par_list[ICS_PAR::ECOM_D0].value : 0;
		solar_pars[ICS_PAR::ECOM_Dc] = (_par_list.count(ICS_PAR::ECOM_Dc) > 0) ? _par_list[ICS_PAR::ECOM_Dc].value : 0;
		solar_pars[ICS_PAR::ECOM_Ds] = (_par_list.count(ICS_PAR::ECOM_Ds) > 0) ? _par_list[ICS_PAR::ECOM_Ds].value : 0;
		solar_pars[ICS_PAR::ECOM_Y0] = (_par_list.count(ICS_PAR::ECOM_Y0) > 0) ? _par_list[ICS_PAR::ECOM_Y0].value : 0;
		solar_pars[ICS_PAR::ECOM_Yc] = (_par_list.count(ICS_PAR::ECOM_Yc) > 0) ? _par_list[ICS_PAR::ECOM_Yc].value : 0;
		solar_pars[ICS_PAR::ECOM_Ys] = (_par_list.count(ICS_PAR::ECOM_Ys) > 0) ? _par_list[ICS_PAR::ECOM_Ys].value : 0;
		solar_pars[ICS_PAR::ECOM_X0] = (_par_list.count(ICS_PAR::ECOM_X0) > 0) ? _par_list[ICS_PAR::ECOM_X0].value : 0;
		solar_pars[ICS_PAR::ECOM_Xc] = (_par_list.count(ICS_PAR::ECOM_Xc) > 0) ? _par_list[ICS_PAR::ECOM_Xc].value : 0;
		solar_pars[ICS_PAR::ECOM_Xs] = (_par_list.count(ICS_PAR::ECOM_Xs) > 0) ? _par_list[ICS_PAR::ECOM_Xs].value : 0;
		return solar_pars;
	}

	/**
	* @brief get 9 solar radiation parameters value of ECOM2 model.
	* @return  9 solar radiation parameters value of ECOM2 model
	*/
	map<ICS_PAR, double> t_gics_nav::get_ecom2_solar_pars_values()
	{
		map<ICS_PAR, double> solar_pars;
		solar_pars[ICS_PAR::ECOM2_D0]  = (_par_list.count(ICS_PAR::ECOM2_D0)  > 0) ? _par_list[ICS_PAR::ECOM2_D0].value  : 0;
		solar_pars[ICS_PAR::ECOM2_D2c] = (_par_list.count(ICS_PAR::ECOM2_D2c) > 0) ? _par_list[ICS_PAR::ECOM2_D2c].value : 0;
		solar_pars[ICS_PAR::ECOM2_D2s] = (_par_list.count(ICS_PAR::ECOM2_D2s) > 0) ? _par_list[ICS_PAR::ECOM2_D2s].value : 0;
		solar_pars[ICS_PAR::ECOM2_D4c] = (_par_list.count(ICS_PAR::ECOM2_D4c) > 0) ? _par_list[ICS_PAR::ECOM2_D4c].value : 0;
		solar_pars[ICS_PAR::ECOM2_D4s] = (_par_list.count(ICS_PAR::ECOM2_D4s) > 0) ? _par_list[ICS_PAR::ECOM2_D4s].value : 0;
		solar_pars[ICS_PAR::ECOM2_Y0]  = (_par_list.count(ICS_PAR::ECOM2_Y0)  > 0) ? _par_list[ICS_PAR::ECOM2_Y0].value  : 0;
		solar_pars[ICS_PAR::ECOM2_X0]  = (_par_list.count(ICS_PAR::ECOM2_X0)  > 0) ? _par_list[ICS_PAR::ECOM2_X0].value  : 0;
		solar_pars[ICS_PAR::ECOM2_X1c] = (_par_list.count(ICS_PAR::ECOM2_X1c) > 0) ? _par_list[ICS_PAR::ECOM2_X1c].value : 0;
		solar_pars[ICS_PAR::ECOM2_X1s] = (_par_list.count(ICS_PAR::ECOM2_X1s) > 0) ? _par_list[ICS_PAR::ECOM2_X1s].value : 0;
		return solar_pars;
	}

	/**
	* @brief get 5 solar radiation parameters value of the a priori box-wing model.
	* @return  5 solar radiation parameters value of the a priori box-wing model
	*/
	map<ICS_PAR, double> t_gics_nav::get_pbw_solar_pars_values()
	{
		map<ICS_PAR, double> solar_pars;
		solar_pars[ICS_PAR::ECOM_D0] = (_par_list.count(ICS_PAR::ECOM_D0) > 0) ? _par_list[ICS_PAR::ECOM_D0].value : 0;
		solar_pars[ICS_PAR::ECOM_Dc] = (_par_list.count(ICS_PAR::ECOM_Dc) > 0) ? _par_list[ICS_PAR::ECOM_Dc].value : 0;
		solar_pars[ICS_PAR::ECOM_Ds] = (_par_list.count(ICS_PAR::ECOM_Ds) > 0) ? _par_list[ICS_PAR::ECOM_Ds].value : 0;
		solar_pars[ICS_PAR::ECOM_Y0] = (_par_list.count(ICS_PAR::ECOM_Y0) > 0) ? _par_list[ICS_PAR::ECOM_Y0].value : 0;
		solar_pars[ICS_PAR::ECOM_Yc] = (_par_list.count(ICS_PAR::ECOM_Yc) > 0) ? _par_list[ICS_PAR::ECOM_Yc].value : 0;
		solar_pars[ICS_PAR::ECOM_Ys] = (_par_list.count(ICS_PAR::ECOM_Ys) > 0) ? _par_list[ICS_PAR::ECOM_Ys].value : 0;
		solar_pars[ICS_PAR::ECOM_X0] = (_par_list.count(ICS_PAR::ECOM_X0) > 0) ? _par_list[ICS_PAR::ECOM_X0].value : 0;
		solar_pars[ICS_PAR::ECOM_Xc] = (_par_list.count(ICS_PAR::ECOM_Xc) > 0) ? _par_list[ICS_PAR::ECOM_Xc].value : 0;
		solar_pars[ICS_PAR::ECOM_Xs] = (_par_list.count(ICS_PAR::ECOM_Xs) > 0) ? _par_list[ICS_PAR::ECOM_Xs].value : 0;
		return solar_pars;
	}

	/**
	* @brief get parameters' size of the satellite.
	* @return  parameters' size of the satellite
	*/
	int t_gics_nav::get_par_num()
	{
		return _par_list.size();
	}

	/**
	* @brief ICS_PAR to t_gpar.
	* @return  t_gpar
	*/
	t_gpar::t_type t_gics_nav::_icspar2gpar(ICS_PAR par_name)
	{
		switch (par_name)
		{
		case great::ICS_PAR::PXSAT:
			return t_gpar::PXSAT;
			break;
		case great::ICS_PAR::PYSAT:
			return t_gpar::PYSAT;
			break;
		case great::ICS_PAR::PZSAT:
			return t_gpar::PZSAT;
			break;
		case great::ICS_PAR::VXSAT:
			return t_gpar::VXSAT;
			break;
		case great::ICS_PAR::VYSAT:
			return t_gpar::VYSAT;
			break;
		case great::ICS_PAR::VZSAT:
			return t_gpar::VZSAT;
			break;
		case great::ICS_PAR::ECOM_D0:
			return t_gpar::ECOM_D0;
			break;
		case great::ICS_PAR::ECOM_Dc:
			return t_gpar::ECOM_Dc;
			break;
		case great::ICS_PAR::ECOM_Ds:
			return t_gpar::ECOM_Ds;
			break;
		case great::ICS_PAR::ECOM_Y0:
			return t_gpar::ECOM_Y0;
			break;
		case great::ICS_PAR::ECOM_Yc:
			return t_gpar::ECOM_Yc;
			break;
		case great::ICS_PAR::ECOM_Ys:
			return t_gpar::ECOM_Ys;
			break;
		case great::ICS_PAR::ECOM_X0:
			return t_gpar::ECOM_B0;
			break;
		case great::ICS_PAR::ECOM_Xc:
			return t_gpar::ECOM_Bc;
			break;
		case great::ICS_PAR::ECOM_Xs:
			return t_gpar::ECOM_Bs;
			break;
		case great::ICS_PAR::ECOM2_D0:
			return t_gpar::ECOM2_D0;
			break;
		case great::ICS_PAR::ECOM2_D2c:
			return t_gpar::ECOM2_D2c;
			break;
		case great::ICS_PAR::ECOM2_D2s:
			return t_gpar::ECOM2_D2s;
			break;
		case great::ICS_PAR::ECOM2_D4c:
			return t_gpar::ECOM2_D4c;
			break;
		case great::ICS_PAR::ECOM2_D4s:
			return t_gpar::ECOM2_D4s;
			break;
		case great::ICS_PAR::ECOM2_Y0:
			return t_gpar::ECOM2_Y0;
			break;
		case great::ICS_PAR::ECOM2_X0:
			return t_gpar::ECOM2_X0;
			break;
		case great::ICS_PAR::ECOM2_X1c:
			return t_gpar::ECOM2_X1c;
			break;
		case great::ICS_PAR::ECOM2_X1s:
			return t_gpar::ECOM2_X1s;
			break;
		case great::ICS_PAR::NO_DEF:
			throw runtime_error("ICS Par NO DEFINE!!! CAN'T CHANGE!");
			break;
		}
		return t_gpar::NO_DEF;
	}

	t_ics_record::t_ics_record()
	{
	}
}