/**
 * @verbatim
	 History
	  -1.0	Jie Li        2019-01-06 creat the file.
	  -1.2  Jiande Huang  2020-09-30 Standard code style and Optimize code logic
   @endverbatim
 *
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gics_nav.h
 * @brief		The base class for save ics file information.
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 */

#ifndef GICS_NAV_H
#define GICS_NAV_H

#include <vector>
#include <map>

#include "gdata/gdata.h"
#include "gforce/gfsolar.h"
#include "gforce/gfcommon.h"
#include "gutils/ExportOI.h"
#include "gutils/gconst_oi.h"
#include "gutils/gtime.h"
#include "gmodels/gpar.h"

using namespace std;

namespace great
{
	/**
	*@brief	   Class for save One satellite record of the ICS
	*
	* This class save one line record of a satellite,include par type 
	*	and corresponding value.
	*/
	class OI_LIBRARY_EXPORT t_ics_record
	{
	public:
		/** @brief default constructor. */
		t_ics_record();

		/** @brief default destructor. */
		~t_ics_record() {};

		ICS_PAR  enum_name;			///< par_type
		string   str_name;			///< par_type

		double  value = 0.0;		///< the value of parameter
		double  sigma = 0.0;		///< sigma of parameter

		t_gtime beg_time;				///< begin time of parameter
		t_gtime end_time;				///< end time of parameter

	};

	/**
	*@brief	   Class for save One satellite record of the ICS
	*
	* This class save one line record of a satellite,include par type
	*	and corresponding value and provide the interface to get data.
	*/
	class OI_LIBRARY_EXPORT t_gics_nav : public t_gdata
	{

	typedef map<ICS_PAR, t_ics_record> t_map_pars;

	public:
		t_gics_nav(string sat);
		/** @brief default constructor. */
		t_gics_nav();
		/** @brief default destructor. */
		virtual ~t_gics_nav();

		string	sat_name;	///< satellite name(GNS: G01, G02, C01)
		/**
		* @brief judge whether the sat record is empty or not.
		* @return  bool
			@retval   true    record is empty
		    @retval   false   record is existent
		*/
		bool is_empty();

		/**
		* @brief set parameter data for one parameter.
		* @param[in]   par_name   par name of string
		* @param[in]   record     par values(maybe more than one record, PWC param)
		*/
		void set_par_record(string par_name, t_ics_record& record);   

		/**
		* @brief get par value according to the par name(string).
		* @param[in]   par_name      par name of enum
		* @return  value of the par 
		*/
		double get_par_value(string par_name);

		/**
		* @brief get par value according to the par name(ICS_PAR).
		* @param[in]   par_name      par name of enum
		* @return  value of the par
		*/
		double get_par_value(ICS_PAR par_name);

		/**
		* @brief get all par value of the satellite.
		* @return  par name and its'value of the satellite
		*/
		t_map_pars& get_par_list();

		/**
		* @brief get 3 solar radiation parameters value of ECOM model.
		* @return  3 solar radiation parameters value of ECOM model
		*/
		map<ICS_PAR, double> get_ecom3_solar_pars_values();

		/**
		* @brief get 5 solar radiation parameters value of ECOM model.
		* @return  5 solar radiation parameters value of ECOM model
		*/
		map<ICS_PAR, double> get_ecom5_solar_pars_values();

		/**
		* @brief get 9 solar radiation parameters value of ECOM model.
		* @return  9 solar radiation parameters value of ECOM model
		*/
		map<ICS_PAR, double> get_ecom9_solar_pars_values();

		/**
		* @brief get 9 solar radiation parameters value of ECOM2 model.
		* @return  9 solar radiation parameters value of ECOM2 model
		*/
		map<ICS_PAR, double> get_ecom2_solar_pars_values();

		/**
		* @brief get 5 solar radiation parameters value of the a priori box-wing model.
		* @return  5 solar radiation parameters value of the a priori box-wing model
		*/
		map<ICS_PAR, double> get_pbw_solar_pars_values();

		/**
		* @brief get parameters' size of the satellite.
		* @return  parameters' size of the satellite
		*/
		int get_par_num();

		/**
		* @brief ICS_PAR to t_gpar.
		* @return  t_gpar
		*/
		static t_gpar::t_type _icspar2gpar(ICS_PAR par_name);

	protected:

		t_map_pars	_par_list;		///< list of all parameters for one satellite
	};

} //namespace

#endif