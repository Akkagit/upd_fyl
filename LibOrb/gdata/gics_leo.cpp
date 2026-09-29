/**
*
* @verbatim
    History
     -1.0 Hongmin Zhang	 2019-01-06 creat the file.
     -1.1 Wei Zhang      2019-03-27 Adding Doxygen Style Code Remarks
	 -1.2 Yujie Qin		 2019-04-10 change the code adjusting to the new iceleo.h
	 -1.3 Yujie Qin		 2019-07-02 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gics_leo.cpp
* @brief		Storage the LEO ICS data for one satellite
*
* @author       Hongmin Zhang, Wuhan University
* @version		1.0.0
* @date		    2019-03-27
*
*/

#include "gdata/gics_leo.h"

using namespace std;

namespace great
{
	/** @brief default constructor. */
	t_icsleo_record::t_icsleo_record()
	{
		enum_name = ICS_PAR::NO_DEF;
		str_name = "";
		isPWC = false;
		value = 0.0;
		sig = 0.0;
		duration = 0.0;
		piecewise_time = 0;
	}

	/** @brief default destructor. */
	t_icsleo_record::~t_icsleo_record()
	{

	}


	/** @brief default constructor. */
	t_gics_leo::t_gics_leo()
	{
		gtrace("t_gicsleo::constructor");
		id_type(t_gdata::ALLICSLEO);
	}

	/** @brief default destructor. */
	t_gics_leo::~t_gics_leo()
	{

	}


	string t_gics_leo::sattype()
	{
		return sat_type;
	}
	/**
	* @brief set parameter data for one parameter.
	* @param[in]   par_name   par name of string
	* @param[in]   record     par values(maybe more than one record, PWC param)
	*/
	void t_gics_leo::set_par_record(string par_name, vector<t_icsleo_record>& record)
	{
		ICS_PAR par = str2icspar(par_name);
		_satleo[par] = record;
		return;
	}

	/**
	* @brief get par value according to the par name(string) and time.
	* @param[in]   par_name      par name of enum
	* @param[in]   time          requested time
	* @return  value of the par in the requested time
	*/
	double t_gics_leo::get_par_value(ICS_PAR par_name, t_gtime time)
	{
		//ICS_PAR par = _str2icspar(par_name);
		vector<t_icsleo_record> record = _satleo[par_name];
		int n = record.size();
		if (n == 0) return 0.0;
		if (time < record[0].beg)
		{
			return record[0].value;
		}
		if (time > record[n - 1].end)
		{
			return record[n - 1].value;
		}
		for (int i = 0; i < n; i++)
		{
			if (time >= record[i].beg && time <= record[i].end)
			{
				return record[i].value;
			}
		}
		return 0.0;
	}


	/**
	* @brief get pars value according to the par name(string) and time(mainly for EMP params).
	* @param[in]   par_first      first three characters of par name(for EMP param)
	* @param[in]   time           requested time
	* @return  value of pars in the requested time(all EMP params)
	*/
	map<ICS_PAR, double> t_gics_leo::get_pars_value(string par_first, t_gtime time)
	{
		map<ICS_PAR, double> parvalue;
		if (par_first == "EMP")
		{
			ICS_PAR cuspar[18] = 
			{
				ICS_PAR::EMP_Ac,  ICS_PAR::EMP_Bc,  ICS_PAR::EMP_Aa, 
				ICS_PAR::EMP_Ba,  ICS_PAR::EMP_Ar,  ICS_PAR::EMP_Br, 
				ICS_PAR::EMP_Sc1, ICS_PAR::EMP_Cc1, ICS_PAR::EMP_Sa1, 
				ICS_PAR::EMP_Ca1, ICS_PAR::EMP_Sr1, ICS_PAR::EMP_Cr1, 
				ICS_PAR::EMP_Sc2, ICS_PAR::EMP_Cc2, ICS_PAR::EMP_Sa2, 
				ICS_PAR::EMP_Ca2, ICS_PAR::EMP_Sr2, ICS_PAR::EMP_Cr2 };
			for (int i = 0; i < 18; i++)
			{
				ICS_PAR EMP = cuspar[i];
				if (_satleo.find(EMP) == _satleo.end())
				{
					continue;
				}
				vector<t_icsleo_record> record = _satleo[EMP];
				int n = _satleo[EMP].size();
				if (time < record[0].beg)
				{
					parvalue[EMP] = record[0].value;
					continue;
				}
				if (time > record[n - 1].end)
				{
					parvalue[EMP] = record[n - 1].value;
					continue;
				}
				for (int i = 0; i < n; i++)
				{
					if (time >= record[i].beg && time <= record[i].end)
					{
						parvalue[EMP] = record[i].value;
						break;
					}
				}
			}
		}
		return parvalue;
	}

	/**
	* @brief return number of parameters to be estimated.
	* @return   number of the par
	*/
	int t_gics_leo::get_par_num()
	{
		return _satleo.size();
	}

	/**
	* @brief return par number of one PWC param.
	* @param[in]  parname      par string name
	* @return  number of PWC param
	*/
	int t_gics_leo::get_pwc_num(string parname)
	{
		ICS_PAR par = str2icspar(parname);
		return _satleo[par].size();
	}

	/**
	* @brief get all param for one LEO satellite
	* @return  all param for one LEO satellite
	*/
	typedef map<ICS_PAR, vector<t_icsleo_record>> t_map_pars_leo;
	t_map_pars_leo* t_gics_leo::leoics()
	{
		return &_satleo;
	}

	/**
	* @brief convert par string name to par enum name.
	* @param[in]  par_name     par string name
	* @return  par enum name
	*/
	ICS_PAR t_gics_leo::str2icspar(string par_name)
	{
		if (par_name == "PXSAT") return ICS_PAR::PXSAT;
		if (par_name == "PYSAT") return ICS_PAR::PYSAT;
		if (par_name == "PZSAT") return ICS_PAR::PZSAT;
		if (par_name == "VXSAT") return ICS_PAR::VXSAT;
		if (par_name == "VYSAT") return ICS_PAR::VYSAT;
		if (par_name == "VZSAT") return ICS_PAR::VZSAT;

		if (par_name == "EMP_Ac") return ICS_PAR::EMP_Ac;
		if (par_name == "EMP_Aa") return ICS_PAR::EMP_Aa;
		if (par_name == "EMP_Ar") return ICS_PAR::EMP_Ar;
		if (par_name == "EMP_Bc") return ICS_PAR::EMP_Bc;
		if (par_name == "EMP_Ba") return ICS_PAR::EMP_Ba;
		if (par_name == "EMP_Br") return ICS_PAR::EMP_Br;

		if (par_name == "EMP_Sc1") return ICS_PAR::EMP_Sc1;
		if (par_name == "EMP_Cc1") return ICS_PAR::EMP_Cc1;
		if (par_name == "EMP_Sa1") return ICS_PAR::EMP_Sa1;
		if (par_name == "EMP_Ca1") return ICS_PAR::EMP_Ca1;
		if (par_name == "EMP_Sr1") return ICS_PAR::EMP_Sr1;
		if (par_name == "EMP_Cr1") return ICS_PAR::EMP_Cr1;

		if (par_name == "EMP_Sc2") return ICS_PAR::EMP_Sc2;
		if (par_name == "EMP_Cc2") return ICS_PAR::EMP_Cc2;
		if (par_name == "EMP_Sa2") return ICS_PAR::EMP_Sa2;
		if (par_name == "EMP_Ca2") return ICS_PAR::EMP_Ca2;
		if (par_name == "EMP_Sr2") return ICS_PAR::EMP_Sr2;
		if (par_name == "EMP_Cr2") return ICS_PAR::EMP_Cr2;

		if (par_name == "SR_scale") return ICS_PAR::SR_scale;
		if (par_name == "DRAG_c")   return ICS_PAR::DRAG_c;
		
		return ICS_PAR::NO_DEF;
	}

	/**
	* @brief convert par enum name to par string name.
	* @param[in]  par_name     par enum name
	* @return  par string name
	*/
	string t_gics_leo::icspar2str(ICS_PAR par_name)
	{
		if (par_name == ICS_PAR::PXSAT) return "PXSAT";
		if (par_name == ICS_PAR::PYSAT) return "PYSAT";
		if (par_name == ICS_PAR::PZSAT) return "PZSAT";
		if (par_name == ICS_PAR::VXSAT) return "VXSAT";
		if (par_name == ICS_PAR::VYSAT) return "VYSAT";
		if (par_name == ICS_PAR::VZSAT) return "VZSAT";
		if (par_name == ICS_PAR::EMP_Sc1) return "EMP_Sc1";
		if (par_name == ICS_PAR::EMP_Cc1) return "EMP_Cc1";
		if (par_name == ICS_PAR::EMP_Sa1) return "EMP_Sa1";
		if (par_name == ICS_PAR::EMP_Ca1) return "EMP_Ca1";
		if (par_name == ICS_PAR::EMP_Sr1) return "EMP_Sr1";
		if (par_name == ICS_PAR::EMP_Cr1) return "EMP_Cr1";

		if (par_name == ICS_PAR::EMP_Sc2) return "EMP_Sc2";
		if (par_name == ICS_PAR::EMP_Cc2) return "EMP_Cc2";
		if (par_name == ICS_PAR::EMP_Sa2) return "EMP_Sa2";
		if (par_name == ICS_PAR::EMP_Ca2) return "EMP_Ca2";
		if (par_name == ICS_PAR::EMP_Sr2) return "EMP_Sr2";
		if (par_name == ICS_PAR::EMP_Cr2) return "EMP_Cr2";

		if (par_name == ICS_PAR::EMP_Ac) return "EMP_Ac";
		if (par_name == ICS_PAR::EMP_Bc) return "EMP_Bc";
		if (par_name == ICS_PAR::EMP_Aa) return "EMP_Aa";
		if (par_name == ICS_PAR::EMP_Ba) return "EMP_Ba";
		if (par_name == ICS_PAR::EMP_Ar) return "EMP_Ar";
		if (par_name == ICS_PAR::EMP_Br) return "EMP_Br";

		if (par_name == ICS_PAR::SR_scale) return "SR_scale";
		if (par_name == ICS_PAR::DRAG_c) return "DRAG_c";

		return "";
	}

	t_gpar::t_type t_gics_leo::icspar2gpar(ICS_PAR par)
	{
		switch (par)
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
		case great::ICS_PAR::SR_scale:
			return t_gpar::SR_scale;
			break;
		case great::ICS_PAR::DRAG_c:
			return t_gpar::DRAG_c;
			break;
		case great::ICS_PAR::EMP_Sa1:
			return t_gpar::EMP_Sa1;
			break;
		case great::ICS_PAR::EMP_Ca1:
			return t_gpar::EMP_Ca1;
			break;
		case great::ICS_PAR::EMP_Sc1:
			return t_gpar::EMP_Sc1;
			break;
		case great::ICS_PAR::EMP_Cc1:
			return t_gpar::EMP_Cc1;
			break;
		case great::ICS_PAR::EMP_Sr1:
			return t_gpar::EMP_Sr1;
			break;
		case great::ICS_PAR::EMP_Cr1:
			return t_gpar::EMP_Cr1;
			break;
		case great::ICS_PAR::EMP_Sa2:
			return t_gpar::EMP_Sa2;
			break;
		case great::ICS_PAR::EMP_Ca2:
			return t_gpar::EMP_Ca2;
			break;
		case great::ICS_PAR::EMP_Sc2:
			return t_gpar::EMP_Sc2;
			break;
		case great::ICS_PAR::EMP_Cc2:
			return t_gpar::EMP_Cc2;
			break;
		case great::ICS_PAR::EMP_Sr2:
			return t_gpar::EMP_Sr2;
			break;
		case great::ICS_PAR::EMP_Cr2:
			return t_gpar::EMP_Cr2;
			break;
		case great::ICS_PAR::EMP_Aa:
			return t_gpar::EMP_Aa;
			break;
		case great::ICS_PAR::EMP_Ba:
			return t_gpar::EMP_Ba;
			break;
		case great::ICS_PAR::EMP_Ac:
			return t_gpar::EMP_Ac;
			break;
		case great::ICS_PAR::EMP_Bc:
			return t_gpar::EMP_Bc;
			break;
		case great::ICS_PAR::EMP_Ar:
			return t_gpar::EMP_Ar;
			break;
		case great::ICS_PAR::EMP_Br:
			return t_gpar::EMP_Br;
			break;
		case great::ICS_PAR::NO_DEF:
			throw runtime_error("ICS Par NO DEFINE!!! CAN'T CHANGE!");
			break;
		}
		return t_gpar::NO_DEF;
	}

	void t_gics_leo::setoneparam(t_icsleo_record& onerec)
	{
		_par_rec.push_back(onerec);
	}

	bool t_gics_leo::is_pwc(t_icsleo_record record)
	{
		
		if (record.str_name.substr(2, 3) != "SAT") {
			record.isPWC = true;
			return true;
		}
		return false;
	}

}