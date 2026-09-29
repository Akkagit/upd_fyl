/**
*
* The format :
* @verbatim
	<force_model_nav>
		<spacecraft ID="GNS">
			<integrator order="11" step_in="60" step_out="300" />
			<gravity model="EGM_12" valid="True" />
			<relativity model="Default" valid="True" />
			<point_mass model="Default" valid="True" />
			<solar_radiation model="ECOM_5" valid="True" />
			<earth_radiation model="NONE" valid="False" />
		</spacecraft>
	</force_model_nav>
  @endverbatim
*
* @verbatim
	History
	 -1.0 Jiande Huang	 2019-01-04 creat the file.
	 -1.1 Jie Li		 2019-03-31 add doxygen style cade
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			getforce_nav.h
* @brief		force model set for navigation satellite orbit integral from XML
*
* @author       Jiande Huang, Wuhan University
* @version		1.0.0
* @date			2019-01-04
*
*/

#ifndef GSETFORCE_H
#define GSETFORCE_H

#define XMLKEY_FORCE             "force_model"
#define XMLKEY_FORCE_ID          "ID"
#define XMLKEY_FORCE_VALID       "valid"
#define XMLKEY_FORCE_MODEL       "model"
#define XMLKEY_FORCE_SPACE_CRAFT "space_craft"

#include <map>
#include <string>
#include <iostream>

#include "gutils/gnss.h"
#include "gutils/gleo.h"
#include "gset/gsetbase.h"

#include "gutils/ExportOI.h"
#include "gutils/gconst_oi.h"
#include "gforce/gfcommon.h"

using namespace std;

namespace great
{
	// string : sat name
	typedef map< string, map< FORCE, bool   > >	t_map_sys_force_model_valid;
	// first string : sat name, second string : force model name
	typedef map< string, map< FORCE, string > >	t_map_sys_force_model_name;

	/**
	 *@brief Class for GNS OI force model setting
	 */
	class OI_LIBRARY_EXPORT t_gsetforce : public virtual t_gsetbase
	{
	public:
		/** @brief default constructor. */
		t_gsetforce();

		/** @brief default destructor. */
		virtual ~t_gsetforce() {};

		/** @brief settings check. */
		void check();

		/** @brief settings help. */
		void help();

		/**
		* @brief get step-in of adams orbit integration(unit: s)
		* @param[in]  sat      satellite name
		* @return    OI input step
		*/
		int		getInpStep(string sat);

		/**
		* @brief get step-out of orbit integration(unit: s)
		* @param[in]  sat      satellite name
		* @return    OI output step
		*/
		int		getOutStep(string sat);

		/**
		* @brief get the integra degree of adams orbit integration
		* @param[in]  sat      satellite name
		* @return    OI adams order
		*/
		int		getIntegraDegree(string sat);

		/**
		* @brief get gravity(ocean_tide and ocean_pole_tide) set
		* @param[in]  sat                  satellite name
		* @param[out] ot_degree            ocean_tide degree
		* @param[out] _opt_valid           ocean_pole_tide valid
		* @param[out] _opt_degree          ocean_pole_tide degree
		* @retrun
		*/
		void    getGravitySet(string sat, int& ot_degree, bool& opt_valid, int& opt_degree);

		/**
		* @brief get solar radiation pressure model list
		* @param[in]  sat                  satellite name 
        * @retrun    LEO solar radiation params
        */
		map<string, int> getSolarRadiationPars(string sat);
		/**
		* @brief get atmosphere model list
		* @param[in]  sat                  satellite name
		* @retrun    LEO atmosphere params
		*/
		map<string, int> getAtmospherePars(string sat);
		/**
		* @brief get customer model list
		* @param[in]  sat               satellite name
		* @return    LEO customer params
		*/
		map<string, int> getCustomerPars(string sat);

		/** @brief get the gravity  degree for one sat.*/
		// int     gravity_degree(string sat);	

		/**
		* @brief get used force for GNS satellite(valid)
		* @param[in]  sat   satellite name
		* @return  the forces used
		*/
		set<FORCE>			getForceUsedSet(string sat);	

		/**
		*@brief get used force model for GNS satellite
		* @param[in]  sat  satellite name
		* @return  force model used
		*/
		map<FORCE, string>	getForceModelMap(string sat);

		/**
		* @brief get all force models
		* @param[in]  sat              satellite name
		* @return    map of force model(true or false)
		*/
		map<FORCE, bool>	getForcesValidMap(string sat);


	protected:
		t_map_sys_force_model_valid _map_force_valid;		///< if use the model
		t_map_sys_force_model_name  _map_force_model;		///< the model that used
		map<string, int> _inp_step;							///< adams inp step
		map<string, int> _out_step;							///< adams out step
		map<string, int> _integra_degree;					///< adams order
		map<string, int> _ot_degree;   ///< ocean_tide degree
		map<string, int> _opt_degree;  ///< ocean_pole_tide_degree
		map<string, bool>_opt_valid;   ///< if use ocean_pole_tide
	private:
		// defalut settings, list here for the convenience of changing
		const int default_inp_setp = 60;
		const int default_out_step = 300;
		const int default_integra_degree = 11;
		const string default_model_GRAVITY    = "EGM08_12";
		const string default_model_SOLAR_RADI = "ECOM_5";
		const string default_model_RELATIVITY = "Default";
		const string default_model_POINT_MASS = "Default";
		const string default_model_EARTH_RADI = "NONE";
		//leo defult
		//map<string, int> _drag_list;
		//map<string, int> _customer_list;
		map<string, map<string, int>> _leo_dg_list;
		map<string, map<string, int>> _leo_cs_list;
		map<string, map<string, int>> _leo_sr_list;
		
		int _flag;								///< flag for updata hole force model info

		/** @brief updata info if xml defined. */
		void _updataForceWithXml();

		/**
		* @brief updata all GNSS sats settings
		*
		* reset force models for all GNSS satellite according to the GNS commom settings in config file
		*
		* @param[in]  spacecraft  set node in XML file
		*/
		void _resetAll(xml_node spcecraft);

		/**
		* @brief updata one system sats settings
		*
		* reset force models for a single system satellite according to the system commom settings in config file
		*
		* @param[in]  spacecraft  set node in XML file
		*/
		void _resetSys(GSYS sys, xml_node spacecraft);	

		/** 
		* @brief updata one sat settings
		*
		* reset force models for one satellite according to the sat commom settings in config file
		*
		* @param[in]  spacecraft  set node in XML file
		*/
		void _resetSat(string sat, xml_node spacecraft);
	};
}
#endif // !GSETFORCE_H
