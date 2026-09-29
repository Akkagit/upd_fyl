/**
 *
 * @verbatim
	 History
	  -1.0  Jie Li       2019-01-06 Creat the file.
	  -1.1  Jie Li       2019-03-30 Code comments are added in Doxygen format.
	  -1.2  Huang Jiande 2020-09-12 Add function get sat list to fix the format of the orb file 
   @endverbatim
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		gorb.h
 * @brief		The class for storaging orb data.
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-03-30
 *
 */


#ifndef GORB_H
#define GORB_H

#include "gutils/ExportBase.h"
#include "gutils/gleo.h"
#include "gdata/gdata.h"
#include "gcoders/orb.h"
#include "gall/gallpar.h"

using namespace gnut;
using namespace great;
namespace great
{
	/**
	*@brief	   Class for storaging one orb file data
	*
	* The class contains orb header information and posation and partial matrix of satellites,
	*		and provide the interface to get orb data.
	*/
	class BASE_LIBRARY_EXPORT t_gorb :public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gorb();

		/** @brief default denstructor. */
		virtual ~t_gorb();

		typedef vector<t_gtriple>  t_sat_funct;

		/**
		* @brief add orb head data.
		* @param[in]   orb_hd   orb head data
		*/
		void add_head(const t_orb_hd& orb_hd);
		
		/**
		* @brief add orb sat data(include sat list,sat num).
		* @param[in]   sys   system
		* @param[in]   sats  sat list
		*/
		void add_sat_info(set<string>& sys, set<string>& sats);

		/**
		* @brief add one record funct data for one epoch one sat to _funct map.
		* @param[in]   epo			  epoch(=(epoch-ref)/out_step)
		* @param[in]   sat			  satellite prn
		* @param[in]   one_funct_rec  one record(one line) of funct of a sat in epo
		*/
		void add_funct(int epo, string sat, t_gtriple& funct);
	
		/**
		* @brief add all record funct data for one epoch from calculate result(stste) to map.
		* @param[in]   epo		epoch(=(epoch-ref)/out_step)
		* @param[in]   sat		satellite prn
		* @param[in]   funct	all funct of one sat in epo
		*/
		void add_funct(int epo, string sat, const vector<double>& funct);
		
		/**
		* @brief add all record funct data for one sat to map.
		* @param[in]   sat		satellite prn
		* @param[in]   funct	all funct of one sat
		*/
		void add_funct(string sat, map<int, t_sat_funct >& funct);

		/**
		* @brief add parlist for one sat to map.
		* @param[in]   sat		satellite prn
		* @param[in]   parlist	parlistt of one sat
		*/
		void add_parlist(string sat, vector<string> parlist);

		/** @brief whether the orb data is empty.
		* @return  bool
		*	@retval   true   orb data is empty
		*   @retval   false  orb data is existent
		*/
		bool is_empty();

		/**
		* @brief get orb head data.
		* @param[out]   orb_hd	orb head data
		*/
		void get_head(t_orb_hd& orb_hd);
		
		/**
		* @brief get whole orb funct data.
		* @param[out]   funct	whole orb funct data
		*/
		void get_funct(map<string, map<int, t_sat_funct > >& funct);

		/**
		* @brief get orb funct data of a satellite.
		* @param[in]   sat		satellite prn
		* @param[out]  funct	orb funct data of the sat
		*/
		void get_funct(string sat, map<int, t_sat_funct >& funct);

		/**
		* @brief interpolate a satellite's position(unit: m) in specific epoch.
		* @param[in]   epoch	epoch time
		* @param[in]   sat		satellite prn
		* @param[out]  px_sat	x coordinates
		* @param[out]  py_sat	y coordinates
		* @param[out]  pz_sat	z coordinates
		* @return   if interpolate successfully
		*	@retval   true   interpolate successfully
		*   @retval   false  interpolate unsuccessfully
		*/
		bool get_xyz(const t_gtime& epoch, string sat, double& px_sat, double& py_sat, double& pz_sat);

		/**
		* @brief interpolate a satellite's velocity(unit: m/s) in specific epoch.
		* @param[in]   epoch	epoch time
		* @param[in]   sat		satellite prn
		* @param[out]  vx_sat	x velocity
		* @param[out]  vy_sat	y velocity
		* @param[out]  vz_sat	z velocity
		* @return   if interpolate successfully
		*	@retval   true   interpolate successfully
		*   @retval   false  interpolate unsuccessfully
		*/
		bool get_vel(const t_gtime& epoch, string sat, double& vx_sat, double& vy_sat, double& vz_sat);
		
		/**
		* @brief interpolate a satellite's partial matrix(only dxdx0..dxdvx0...dxdnp) in specific epoch.
		* @param[in]   epoch	epoch time
		* @param[in]   sat		satellite prn
		* @param[out]  partial	partial matrix
		* @return   if interpolate successfully
		*	@retval   true   interpolate successfully
		*   @retval   false  interpolate unsuccessfully
		*/
		bool get_partial(const t_gtime& epo, string sat, Matrix& partial); // interpolating

		/**
		* @brief interpolate a LEO satellite's partial matrix(especially for PWC parameters) in specific epoch.
		* @param[in]   epoch	epoch time
		* @param[in]   sat		satellite prn
		* @param[in]   orbpar	par list
		* @param[out]  partial	partial matrix
		* @return   if interpolate successfully
		*	@retval   true   interpolate successfully
		*   @retval   false  interpolate unsuccessfully
		*/
		bool get_partial_leo(const t_gtime& epo, string sat, Matrix& partial, vector<t_gpar>& orbpar);
		
		/**
		* @brief interpolate a satellite's partial matrix(both current pos and vel) in specific epoch.
		* @param[in]   epo_left		first interpolate point
		* @param[in]   epo_right	last interpolate point
		* @param[in]   epo_now		current epoch time
		* @param[in]   sat			satellite prn
		* @param[out]  partial		partial matrix(pos and vel)
		* @return   if interpolate successfully
		*	@retval   true   interpolate successfully
		*   @retval   false  interpolate unsuccessfully
		*/
		bool get_partial_interp(int epo_left, int epo_right, double epo_now, string sat, Matrix& partial);

		/**
		* @brief get the interpolate point
		* @param[in]   epo			epoch input(int)
		* @param[in]   norder		lagrange interpolate order
		* @param[in]   sat			satellite prn
		* @param[out]  epo_left		first interpolate point
		* @param[out]  epo_right	last interpolate point
		* @param[out]  epo_now		current epoch time(mjd)
		* @return   if get interpolate point successfully
		*	@retval   true   get interpolate poin successfully
		*   @retval   false  get interpolate poin unsuccessfully
		*/
		bool get_interp_points(const t_gtime& epo, int norder, string sat, int& epo_left, int& epo_right, double& epo_now);

		/**
		* @brief get orb funct data of a satellite.
		* @param[in]   sat		satellite prn
		* @param[out]  parlist	parlist data of the sat
		*/
		void get_parlist(map<string,vector<string> >& parlist);

	protected:
		t_orb_hd	_orb_hd;	///< orb head data	
		map<string, vector<string> > _par_list; ///< partial order of sat(sat,parlist)
		map<string, map<int, t_sat_funct > > _funct;  ///< orb funct (sat , epo  ,data)
	};

}

#endif