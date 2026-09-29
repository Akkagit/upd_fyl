/**
 *
 * @verbatim
	 History
	  -1.0	Jie Li   2019-01-11 Creat the file.
	  -1.1  Jie Li	 2019-03-30 Code comments are added in Doxygen format.
   @endverbatim
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		gsatparam.h
 * @brief		The class for storaging one sat parameters data.
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-03-30
 *
 */

#ifndef GSATPARAM_H
#define GSATPARAM_H

#include "gdata/gdata.h"
#include "gutils/ExportBase.h"

using namespace std;
using namespace gnut;

namespace great
{
	/**
	*@brief	   Class for storaging one satellite's parameters data
	*
	* The class contains one satellite's parameters for file
	*		and provide the interface to get parameters' value.
	*/
	class  BASE_LIBRARY_EXPORT t_gsatparam :public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gsatparam();

		/** @brief default destructor. */
		virtual ~t_gsatparam();

		bool valid = false;		///< effectiveness of the sat data

		/// set satelllite information(include prn,svn,start time and so on)
		/**
		* @brief save satellite prn from file.
		* @param[in]  prn	satellite prn
		*/
		void set_prn(string prn);

		/**
		* @brief save satellite svn from file.
		* @param[in]  svn	satellite svn
		*/
		void set_svn(string svn);

		/**
		* @brief save satellite launched time from file.
		* @param[in]  start		satellite launched time
		*/
		void set_start(t_gtime start);

		/**
		* @brief save satellite decommissioned time from file.
		*
		* When the end time is equen to 0 or LAST_TIME, that means 
		*	the satellite is still in service.
		*
		* @param[in]  end		decommissioned time
		*/
		void set_end(t_gtime end);

		/**
		* @brief save satellite's cosparid from file.
		* @param[in]  cosparid		satellite's cosparid
		*/
		void set_cosparid(string cosparid);

		/**
		* @brief save satellite's weight from file.
		* @param[in]  mass		satellite's weight(kg)
		*/
		void set_mass(double mass);

		/**
		* @brief save satellite's max_yaw from file.
		* @param[in]  maxyaw		satellite's max_yaw
		*/
		void set_maxyaw(double maxyaw);

		/**
		* @brief save satellite's frequency id from file.
		*
		* Only the GLONASS satellites need to be noticed
		*
		* @param[in]  fid		satellite's frequency id
		*/
		void set_fid(int fid);

		/**
		* @brief save satellite's block type from file.
		* @param[in]  blocktype		satellite's block type
		*/
		void set_blocktype(string blocktype);

		/// return satelllite information
		/**
		* @brief get satellite prn.
		* @return	 satellite prn.
		*/
		string prn();

		/**
		* @brief get satellite svn.
		* @return	 satellite svn.
		*/
		string svn();

		/**
		* @brief get satellite launched time.
		* @return	 satellite launched time.
		*/
		t_gtime start();

		/**
		* @brief get satellite decommissioned time.
		* @return	 satellite decommissioned time.
		*/
		t_gtime end();

		/**
		* @brief get satellite's cosparid.
		* @return	 satellite's cosparid.
		*/
		string cosparid();

		/**
		* @brief get satellite's weight.
		* @return	 satellite's weight.
		*/
		double mass();

		/**
		* @brief get satellite's max yaw.
		* @return	 satellite's max yaw.
		*/
		double maxyaw();

		/**
		* @brief get satellite's frequency id.
		* @return	 satellite's frequency id.
		*/
		int fid();

		/**
		* @brief get satellite's block type.
		* @return	 satellite's block type.
		*/
		string blocktype();

	protected:
		string  _prn;		///< satellite prn
		string  _svn;		///< satelllite svn
		t_gtime _start;		///< launch time of satellite 
		t_gtime _end;		///< retired time of satellite
		string _cosparid;	///< cospar_id of satellite
		double _mass;		///< weight of satellite(unit: kg)
		double _maxyaw;		///< max_yaw of satellite(unit: degree)
		int _fid;			///< frequence id of satellite
		string _blocktype;	///< block type
	};
}
#endif // !GSATPARAM_H
