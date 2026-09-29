/**
*
* @verbatim
    History
     -1.0 Hongmin Zhang 2019-01-06 creat the file.
     -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
	 -1.2 Yujie Qin     2019-04-10 change the code adjusting to the new iceleo.h
	 -1.3 Yujie Qin     2019-07-02 Adding Doxygen Style Code Remarks
	 -1.4 Jiande Huang  2020-09-12 Standard code style and Optimize code logic
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gics_leo.h
* @brief		Storage the LEO ICS data for one satellite
*
* @author       Hongmin Zhang, Wuhan University
* @version		1.0.0
* @date		    2019-01-06
*
*/




#ifndef GICS_LEO_H
#define GICS_LEO_H

#include "gdata/gdata.h"
#include "gforce/gfsolar.h"
#include "gforce/gfcommon.h"
#include "gutils/ExportOI.h"
#include "gutils/gconst_oi.h"
#include "gmodels/gpar.h"
#include "gutils/gtime.h"

using namespace std;

namespace great
{
	/**
	*@brief	   Class for one record of LEO ICS data
	*/
	class OI_LIBRARY_EXPORT t_icsleo_record
	{
	public:
		/** @brief default constructor. */
		t_icsleo_record();

		/** @brief default destructor. */
		~t_icsleo_record();

		ICS_PAR enum_name;              ///< par name(enum)
		string  str_name;               ///< par name(string)
									    
		bool    isPWC;				    ///< parameter is the PWC(piecewise constant) parameter or not
		double  value;				    ///< the value of parameter
		double  sig;				    ///< sigma of parameter

		t_gtime beg;				    ///< begin time of parameter
		t_gtime end;				    ///< end time of parameter

		double  duration;			    ///< last time of parameter
		int     piecewise_time;		    ///< time of piecewise linearity

	};

	/**
	*@brief	   Class for one LEO ICS data storaging
	*/
	class OI_LIBRARY_EXPORT t_gics_leo : public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gics_leo();

		/** @brief default destructor. */
		virtual ~t_gics_leo();

		string  sat_name;    ///< satellite name(GNS: G01, G02, C01,.../LEO: GRACE-A, SWARM-A,...)
		string  sat_type;    ///< satellite type(SWARM-A,SWARM-B,SWARM-C,...)
		
		typedef map<ICS_PAR, vector<t_icsleo_record>> t_map_pars_leo;

		/**
	     * @brief set parameter data for one parameter.
	     * @param[in]   par_name   par name of string
	     * @param[in]   record     par values(maybe more than one record, PWC param)
	     */
		void set_par_record(string par_name, vector<t_icsleo_record>& record);
		
		//void set_par_record_leo(string par_name);   // set parameter data for one parameter (for LEO)
		/** @brief get par value according to the par name(string) and time. */
		//double get_par_value(string par_name, t_gtime time);

		/**
	     * @brief get par value according to the par name(string) and time.
	     * @param[in]   par_name      par name of enum
	     * @param[in]   time          requested time
	     * @return  value of the par in the requested time
	     */
		double get_par_value(ICS_PAR par_name, t_gtime time);

		/**
	     * @brief get pars value according to the par name(string) and time(mainly for EMP params).
	     * @param[in]   par_first      first three characters of par name(for EMP param)
	     * @param[in]   time           requested time
	     * @return  value of pars in the requested time(all EMP params)
	     */
		map<ICS_PAR, double> get_pars_value(string par_first, t_gtime time);

		/**
		 * @brief set one record for one param.
		 * @param[in]   onerec      par values(maybe more than one record, PWC param)
		 */
		void setoneparam(t_icsleo_record& onerec);

		/**
	     * @brief return number of parameters to be estimated.
	     * @return   number of the par
	     */
		int get_par_num();

		/**
	     * @brief convert par string name to par enum name.
	     * @param[in]  par_name     par string name
	     * @return  par enum name
	     */
		static ICS_PAR str2icspar(string par_name);

		/**
	     * @brief convert par enum name to par string name.
	     * @param[in]  par_name     par enum name
	     * @return     par string name
	     */
		string icspar2str(ICS_PAR par_name);

		/**
	     * @brief get all param for one LEO satellite
	     * @return  all param for one LEO satellite
	     */
		t_map_pars_leo* leoics();

		/**
	     * @brief return par number of one PWC param.
	     * @param[in]  parname      par string name
	     * @return     number of PWC param
	     */
		int get_pwc_num(string parname);

		/**
		 * @brief convert ICSPAR to great partype.
	     * @param[in]  par parameter to convert
	     * @return     the parameter in t_type
	     */
		static t_gpar::t_type icspar2gpar(ICS_PAR par);

		/**
		 * @brief get the sat type. 
		 * @return    the type of sat(SWARM-A,SWARM-B,...)
		 */
		string  sattype();

		bool is_pwc(t_icsleo_record record);
	protected:

		t_map_pars_leo _satleo;            ///< map of all params for one LEO satellite
		vector<t_icsleo_record> _par_rec;  ///< data of one parameter(maybe more than one record, PWC param)

	};
} // namespace

#endif