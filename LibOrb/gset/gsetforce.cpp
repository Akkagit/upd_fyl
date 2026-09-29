/**
*
* @verbatim
	History
	 -1.0 Jie Li		 2019-01-10 creat the file.
	 -1.1 Jie Li      2019-04-01 add code in doxygen style.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			getforce_nav.cpp
* @brief		force model set from XML
*
* @author       Jie Li, Wuhan University
* @version		1.0.0
* @date			2019-01-10
*
*/

#include <iomanip>
#include <sstream>
#include <algorithm>

#include "gset/gsetforce.h"
#include "gutils/gsys.h"
using namespace std;
using namespace pugi;


namespace great
{
	/** @brief constructor for navigation satellite. */
	t_gsetforce::t_gsetforce()
		: t_gsetbase()
	{
		_gmutex.lock();

		_set.insert(XMLKEY_FORCE);


		//init GNS setting
		t_map_sats gnss_sats = GNSS_SATS();

		for (auto itGNSS = gnss_sats.begin(); itGNSS != gnss_sats.end(); ++itGNSS)
		{
			GSYS gsys = itGNSS->first;
			set<string> sats = itGNSS->second;

			for (auto itSAT = sats.begin(); itSAT != sats.end(); ++itSAT)
			{
				string sat = *itSAT;
				_map_force_valid[sat][GRAVITY]    = true;
				_map_force_valid[sat][SOLAR_RADI] = true;
				_map_force_valid[sat][RELATIVITY] = true;
				_map_force_valid[sat][POINT_MASS] = true;
				_map_force_valid[sat][EARTH_RADI] = false;

				_map_force_model[sat][GRAVITY]    = default_model_GRAVITY   ;
				_map_force_model[sat][SOLAR_RADI] = default_model_SOLAR_RADI;
				_map_force_model[sat][RELATIVITY] = default_model_RELATIVITY;
				_map_force_model[sat][POINT_MASS] = default_model_POINT_MASS;
				_map_force_model[sat][EARTH_RADI] = default_model_EARTH_RADI;

				_inp_step[sat] = default_inp_setp;
				_out_step[sat] = default_out_step;
				_integra_degree[sat] = default_integra_degree;

				_ot_degree[sat] = 30;
				_opt_valid[sat] = false;
				_opt_degree[sat] = 30;
				_flag = 0;
			}
		}
		//init LEO setting
		set<GLEOLONG> all_leo = LEO_SUPPORTED();
		set<GLEOLONG>::iterator it_leo;
		for (it_leo = all_leo.begin(); it_leo != all_leo.end(); it_leo++)
		{
			string sat = long2str(*it_leo);

			_map_force_valid[sat][GRAVITY] = true;
			_map_force_valid[sat][RELATIVITY] = true;
			_map_force_valid[sat][CUSTOMER] = true;
			_map_force_valid[sat][POINT_MASS] = true;
			_map_force_valid[sat][SOLAR_RADI] = true;
			_map_force_valid[sat][ATMOSPHERE] = true;
			_map_force_valid[sat][EARTH_RADI] = false;

			_map_force_model[sat][EARTH_RADI] = "NONE";
			_map_force_model[sat][GRAVITY] = "EGM08_120";
			_map_force_model[sat][RELATIVITY] = "Default";
			_map_force_model[sat][POINT_MASS] = "Default";
			_map_force_model[sat][SOLAR_RADI] = "BOX-WING";
			_map_force_model[sat][ATMOSPHERE] = "DTM94";

			map<string, int> customer_list;
			map<string, int> drag_list;
			customer_list["Sc1"] = 90;
			customer_list["Cc1"] = 90;
			customer_list["Sa1"] = 90;
			customer_list["Ca1"] = 90;
			customer_list["Sr1"] = 90;
			customer_list["Cr1"] = 90;
			drag_list["DRAG_c"] = 360;
			_leo_dg_list[sat] = drag_list;
			_leo_cs_list[sat] = customer_list;

			_inp_step[sat] = 5;
			_out_step[sat] = 30;
			_integra_degree[sat] = 7;

			_ot_degree[sat] = 30;
			_opt_valid[sat] = true;
			_opt_degree[sat] = 30;

			_flag = 0;
		}

		_gmutex.unlock();
	}

	/**
	* @brief get used force for GNS satellite(valid)
	* @param[in]  sat   satellite name
	* @return  the forces used
	*/
	set<FORCE> t_gsetforce::getForceUsedSet(string sat)
	{
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}

		//不会存在不是卫星的string输入，因为一开始从xml读卫星列表的时候就已经判断是不是一颗卫星了
		set<FORCE> force_used;
		map< FORCE, bool > force_tmp = _map_force_valid[sat];

		map< FORCE, bool >::iterator itUsed = force_tmp.begin();
		for (itUsed; itUsed != force_tmp.end(); ++itUsed)
		{
			if (itUsed->second)
			{
				force_used.insert(itUsed->first);
			}
		}
		return force_used;

		/*GSYS gsys = str2gsys(sat.substr(0, 1));
		set<FORCE> force_used;
		if (gsys == GNS)
		{
			force_used.insert(GRAVITY);
			force_used.insert(RELATIVITY);
			force_used.insert(POINT_MASS);
			force_used.insert(SOLAR_RADI);
			return force_used;
		}
		else
		{
			map< FORCE, bool > force_tmp = _map_force_valid[sat];

			map< FORCE, bool >::iterator itUsed = force_tmp.begin();
			for (itUsed; itUsed != force_tmp.end(); ++itUsed)
			{
				if (itUsed->second)
				{
					force_used.insert(itUsed->first);
				}
			}
			return force_used;
		}*/
	}

	/**
	* @brief get used force model for GNS satellite
	* @param[in]  sat  satellite name
	* @return  force model used
	*/
	map<FORCE, string> t_gsetforce::getForceModelMap(string sat)
	{
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		//两种卫星设置不一样
		return _map_force_model[sat];
		/*GSYS gsys = str2gsys(sat.substr(0, 1));
		if (gsys == GNS)
		{
			map<FORCE, string> tmp_force_model;
			tmp_force_model[GRAVITY]    = "EGM08_12";
			tmp_force_model[SOLAR_RADI] = "ECOM_5";
			return tmp_force_model;
		}
		else
		{
			return _map_force_model[sat];
		}*/
	}

	/**
	* @brief get all force models
	* @param[in]  sat              satellite name
	* @return    map of force model(true or false)
	*/
	map<FORCE, bool> t_gsetforce::getForcesValidMap(string sat)
	{
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		if (_map_force_valid.count(sat) == 0)
		{
			map<FORCE, bool> tmp_vaild;
			tmp_vaild[GRAVITY] = true;
			tmp_vaild[RELATIVITY] = true;
			tmp_vaild[CUSTOMER] = true;
			tmp_vaild[POINT_MASS] = true;
			tmp_vaild[SOLAR_RADI] = true;
			tmp_vaild[ATMOSPHERE] = true;
			tmp_vaild[EARTH_RADI] = false;
			return tmp_vaild;
		}
		else
		{
			return _map_force_valid[sat];
		}
		
		/*GSYS gsys = str2gsys(sat.substr(0, 1));
		set<string> st_pars;
		if (gsys == GNS)
		{
			map<FORCE, bool> tmp_vaild;
			tmp_vaild[GRAVITY] = true;
			tmp_vaild[RELATIVITY] = true;
			tmp_vaild[POINT_MASS] = true;
			tmp_vaild[SOLAR_RADI] = true;
			tmp_vaild[EARTH_RADI] = false;
			return tmp_vaild;
		}
		else
		{
			return _map_force_valid[sat];
		}*/
	}

	/**
	* @brief update all force models for GNS satellite
	*/
	void t_gsetforce::_updataForceWithXml()
	{
		//get force model node
		xml_node force_model = _doc.child(XMLKEY_ROOT).child(XMLKEY_FORCE);

		//reset forcemodel info by xml
		if (!force_model)
		{
			return;
		}

		map<int, set<string> > id_all;
		for (xml_node spacecraft = force_model.first_child(); spacecraft; spacecraft = spacecraft.next_sibling())
		{
			//get force model ID and classificate
			string ID = spacecraft.first_attribute().value();
			if (ID == "GNS"||ID=="LEO")
			{
				id_all[1].insert(ID);
			}
			else if (ID == "GPS" || ID == "GLO" || ID == "GAL" || ID == "BDS" || ID == "QZS")
			{
				id_all[2].insert(ID);
			}
			else
			{
				id_all[3].insert(ID);
			}
		}

		//from low-level sat group to sat sys to sat
		for (auto itID = id_all.begin(); itID != id_all.end(); ++itID)
		{
			set<string> tmp = itID->second;

			if (tmp.empty())
			{
				continue;
			}

			for (auto itVALUE = tmp.begin(); itVALUE != tmp.end(); ++itVALUE)
			{
				xml_node tmp_sc = force_model.find_child_by_attribute("spacecraft", "ID", (*itVALUE).c_str());
				int switch_case = itID->first;
				switch (switch_case)
				{
				case 1:
					//reset GNS forcemodel
					_resetAll(tmp_sc);
					break;
				case 2:
					//reset single sys forcemodel
					_resetSys(t_gsys::str2gsys(*itVALUE), tmp_sc);
					break;
				case 3:
					//reset single sat forcemodel
					_resetSat(*itVALUE, tmp_sc);
					break;
				default:
					break;
				}
			}
		}
	}

	/**
	* @brief get step-in of adams orbit integration(unit: s)
	* @param[in]  sat      satellite name
	* @return    OI input step
	*/
	int t_gsetforce::getInpStep(string sat)
	{
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _inp_step[sat];
		/*GSYS gsys = str2gsys(sat.substr(0, 1));

		if (gsys == GNS)
		{
			return 60;
		}
		else
		{
			return _inp_step[sat];
		}*/
	}

	/**
	* @brief get step-out of orbit integration(unit: s)
	* @param[in]  sat      satellite name
	* @return    OI output step
	*/
	int t_gsetforce::getOutStep(string sat)
	{
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _out_step[sat];
		/*GSYS gsys = str2gsys(sat.substr(0, 1));

		if (gsys == GNS)
		{
			return 300;
		}
		else
		{
			return _out_step[sat];
		}*/
	}

	/**
	* @brief get order of adams orbit integration
	* @param[in]  sat      satellite name
	* @return    OI adams order
	*/
	int t_gsetforce::getIntegraDegree(string sat)
	{
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _integra_degree[sat];
		/*GSYS gsys = str2gsys(sat.substr(0, 1));

		if (gsys == GNS)
		{
			return 11;
		}
		else
		{
			return _integra_degree[sat];
		}*/
	}
	/**
	* @brief get gravity(ocean_tide and ocean_pole_tide) set 
	* @param[in]  sat                  satellite name
	* @param[out] ot_degree            ocean_tide degree
	* @param[out] _opt_valid           ocean_pole_tide valid
	* @param[out] _opt_degree          ocean_pole_tide degree
	* @retrun
	*/
	void t_gsetforce::getGravitySet(string sat, int& ot_degree, bool& opt_valid, int& opt_degree)
	{
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		ot_degree = _ot_degree[sat];
		opt_valid = _opt_valid[sat];
		opt_degree = _opt_degree[sat];
	}

	/**
	* @brief get atmosphere model list
	* @param[in]  sat                  satellite name
	* @retrun    LEO atmosphere params
	*/
	map<string, int> t_gsetforce::getAtmospherePars(string sat)
	{
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _leo_dg_list[sat];
	}

	/**
    * @brief get solar radiation pressure model list
    * @param[in]  sat                  satellite name
    * @retrun    LEO solar radiation params
    */
	map<string, int> t_gsetforce::getSolarRadiationPars(string sat) {
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _leo_sr_list[sat];
	}

	/**
	* @brief get customer model list
	* @param[in]  sat               satellite name
	* @return    LEO customer params
	*/
	map<string, int> t_gsetforce::getCustomerPars(string sat)
	{
		if (_flag == 0)
		{
			_updataForceWithXml();
			_flag = 1;
		}
		return _leo_cs_list[sat];
	}

	/**
	* @brief reset force models for all GNSS satellite according to the GNS commom settings in config file
	* @param[in]  spacecraft  set node in XML file
	*/
	void t_gsetforce::_resetAll(xml_node spacecraft)
	{
		t_map_sats gnss_data = GNSS_SATS();
		for (auto itGNSS = gnss_data.begin(); itGNSS != gnss_data.end(); ++itGNSS)
		{
			if (itGNSS->first == LEO) continue;
			_resetSys(itGNSS->first, spacecraft);
		}

		set<GLEOLONG> leos = LEO_SUPPORTED();
		set<GLEOLONG>::iterator it_leo;
		for (it_leo = leos.begin(); it_leo != leos.end(); it_leo++)
		{
			_resetSat(long2str(*it_leo), spacecraft);
		}
	}

	/**
	* @brief reset force models for a single system satellite according to the system commom settings in config file
	* @param[in]  spacecraft  set node in XML file
	*/
	void t_gsetforce::_resetSys(GSYS sys, xml_node spacecraft)
	{
		t_map_sats gnss_data = GNSS_SATS();
		for (auto itSAT = gnss_data[sys].begin(); itSAT != gnss_data[sys].end(); ++itSAT)
		{
			_resetSat(*itSAT, spacecraft);
		}
	}

	/**
	* @brief reset force models for one satellite according to the sat commom settings in config file
	* @param[in]  spacecraft  set node in XML file
	*/
	void t_gsetforce::_resetSat(string sat, xml_node spacecraft)
	{
		//GSYS sys = str2gsys(sat.substr(0, 1));
		for (xml_node tmp_set = spacecraft.first_child(); tmp_set; tmp_set = tmp_set.next_sibling())
		{

			if ((strcmp("integrator", tmp_set.name())) == 0)
			{
				_inp_step[sat] = str2int(tmp_set.attribute("step_in").value());
				_out_step[sat] = str2int(tmp_set.attribute("step_out").value());
				_integra_degree[sat] = str2int(tmp_set.attribute("order").value());
			}

			if ((strcmp("gravity", tmp_set.name())) == 0)
			{
				_map_force_valid[sat][GRAVITY] = str2bool(tmp_set.attribute("valid").value());
				_map_force_model[sat][GRAVITY] = tmp_set.attribute("model").value();
				//lijie add ocean_tide and ocean_pole_tide select
				for (xml_node child_set = tmp_set.first_child(); child_set; child_set = child_set.next_sibling())
				{
					if ((strcmp("ocean_tide", child_set.name())) == 0)
					{
						_ot_degree[sat] = str2int(child_set.attribute("degree").value());
					}
					if ((strcmp("ocean_pole_tide", child_set.name())) == 0)
					{
						_opt_valid[sat] = str2bool(child_set.attribute("valid").value());
						_opt_degree[sat] = str2int(child_set.attribute("degree").value());
					}
				}


			}

			if ((strcmp("relativity", tmp_set.name())) == 0)
			{
				_map_force_valid[sat][RELATIVITY] = str2bool(tmp_set.attribute("valid").value());
			}
			if ((strcmp("point_mass", tmp_set.name())) == 0)
			{
				_map_force_valid[sat][POINT_MASS] = str2bool(tmp_set.attribute("valid").value());

				string word;
				//_map_point_mass_list[sys][sat].clear();
				istringstream is(tmp_set.first_child().value());
				while (is >> word)
				{
					//_map_point_mass_list[sys][sat].insert(word);
				}
			}
			if ((strcmp("solar_radiation", tmp_set.name())) == 0)
			{
				_map_force_valid[sat][SOLAR_RADI] = str2bool(tmp_set.attribute("valid").value());
				_map_force_model[sat][SOLAR_RADI] = tmp_set.attribute("model").value();

				string word;
				map<string, int> srp;
				set<string> srppar;
				_leo_sr_list[sat].clear();
				istringstream is(tmp_set.first_child().value());
				while (is >> word) {
					srppar.insert(word);
				}
				set<string>::iterator ipara = srppar.begin();
				for (ipara; ipara != srppar.end(); ipara++)
				{
					size_t m = (*ipara).find_first_of(":");
					size_t n = (*ipara).size();
					string name = (*ipara).substr(0, m);
					int min = str2int((*ipara).substr(m + 1, n));
					srp[name] = min;
				}
				_leo_sr_list[sat] = srp;
			}
			if ((strcmp("earth_radiation", tmp_set.name())) == 0)
			{
				_map_force_valid[sat][EARTH_RADI] = str2bool(tmp_set.attribute("valid").value());
				_map_force_model[sat][EARTH_RADI] = tmp_set.attribute("model").value();
			}

			
			if ((strcmp("customer", tmp_set.name())) == 0)
			{
				_map_force_valid[sat][CUSTOMER] = str2bool(tmp_set.attribute("valid").value());
				string word, parname;
				int gap;
				set<string> cuspar;
				map<string, int> cust;
				_leo_cs_list[sat].clear();
				istringstream is(tmp_set.first_child().value());
				while (is >> word) {
					cuspar.insert(word);
				}
				set<string>::iterator ipara = cuspar.begin();
				for (ipara; ipara != cuspar.end(); ipara++)
				{
					size_t m = (*ipara).find_first_of(":");
					size_t n = (*ipara).size();
					string name = (*ipara).substr(0, m);
					int min = str2int((*ipara).substr(m + 1, n));
					cust[name] = min;
				}
				_leo_cs_list[sat] = cust;
			}
			if ((strcmp("atmosphere", tmp_set.name())) == 0)
			{
				_map_force_valid[sat][ATMOSPHERE] = str2bool(tmp_set.attribute("valid").value());
				_map_force_model[sat][ATMOSPHERE] = tmp_set.attribute("model").value();
				string word, parname;
				int gap;
				map<string, int> drag;
				set<string> dragpar;
				_leo_dg_list[sat].clear();
				istringstream is(tmp_set.first_child().value());
				while (is >> word) {
					dragpar.insert(word);
				}
				set<string>::iterator ipara = dragpar.begin();
				for (ipara; ipara != dragpar.end(); ipara++)
				{
					size_t m = (*ipara).find_first_of(":");
					size_t n = (*ipara).size();
					string name = (*ipara).substr(0, m);
					int min = str2int((*ipara).substr(m + 1, n));
					drag[name] = min;
				}
				_leo_dg_list[sat] = drag;
			}

		}
	}

	/**
	* @brief settings check
	*/
	void t_gsetforce::check()
	{
		_gmutex.lock();
		_gmutex.unlock();
		return;
	}

	/**
	* @brief settings help
	*/
	void t_gsetforce::help()
	{
		_gmutex.lock();
		
		// jdhuang
		cerr << "<force_model_nave>                                                       " << endl
			 << "	<spacecraft ID = \"GNS\">                                             " << endl
			 << "		<integrator step_in = \"60\" step_out = \"300\"  order = \"11\" />" << endl
			 << "		<gravity          valid = \"True\"  model = \"EGM_12\" />         " << endl
			 << "		<relativity       valid = \"True\"  model = \"Default\"/>         " << endl
			 << "		<point_mass       valid = \"True\"  model = \"Default\"/>         " << endl
			 << "		<solar_radiation  valid = \"True\"  model = \"ECOM_5\" />         " << endl
			 << "		<earth_radiation  valid = \"False\" model = \"NONE\"   />         " << endl
			 << "	</spacecraft>                                                         " << endl
			 << endl
			 << "<!--> if your have special setting for one system, you can set as following : <!-->" << endl
			 << "	<spacecraft ID = \"GPS\">                                             " << endl
			 << "		<integrator step_in = \"60\" step_out = \"300\"  order = \"11\" />" << endl
			 << "		<gravity    valid = \"True\" model = \"EGM_12\" />                " << endl
			 << "		<relativity valid = \"True\" model = \"Default\"  />              " << endl
			 << "		<point_mass valid = \"True\" model = \"Default\"  />              " << endl
			 << "		<solar_radiation  valid = \"True\"  model = \"ECOM_5\" />         " << endl
			 << "		<earth_radiation  valid = \"False\" model = \"NONE\"   />         " << endl
			 << "	</spacecraft>                                                         " << endl 
			 << endl
			 << "<!--> if your have special setting for one satellite, you can set as following : <!-->" << endl
			 << "	<spacecraft ID = \"G01\">                                             " << endl
			 << "		<integrator step_in = \"60\" step_out = \"300\"  order = \"11\" />" << endl
			 << "		<gravity    valid = \"True\" model = \"EGM_12\" />                " << endl
			 << "		<relativity valid = \"True\" model = \"Default\"  />              " << endl
			 << "		<point_mass valid = \"True\" model = \"Default\"  />              " << endl
			 << "		<solar_radiation  valid = \"True\"  model = \"ECOM_5\" />         " << endl
			 << "		<earth_radiation  valid = \"False\" model = \"NONE\"   />         " << endl
			 << "	</spacecraft>  		                                                  " << endl	
			 << endl;
		_gmutex.unlock();
		return;
	}

}