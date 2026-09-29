/**
*
* @verbatim
	History
	 -1.0 Jie Li	 2019-01-11 creat the file.
	 -1.1 Jie Li     2019-03-30 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gallsatparam.h
* @brief		Storage the satellites parameters data
*
* @author       Jie Li, Wuhan University
* @version		1.0.0
* @date		    2019-03-30
*
*/

#ifndef GALLSATPARAM_H
#define GALLSATPARAM_H

#include "gutils/ExportBase.h"

#include "gdata/gdata.h"
#include "gdata/gsatparam.h"

namespace great
{
	/**
	*@brief	   Class for all satellites parameters data storaging
	*/
	class BASE_LIBRARY_EXPORT t_gallsatinfo : public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gallsatinfo();

		/** @brief default destructor. */
		virtual ~t_gallsatinfo();

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
		int addData(std::string sat_prn, std::string sat_svn, 
			            t_gtime beg_time, t_gtime end_end, 
						std::string cosparid, double sat_mass, double sat_max_yaw, 
						int sat_fid, std::string block_type);

		/** @brief whether the satellite parameters data is empty.
		* @return  bool
		*	@retval   true    satellite parameters data is empty
		*   @retval   false   satellite parameters data is existent
		*/
		bool isEmpty();

		/**
		* @brief get valid satellite parameters data.
		*
		* Get paramter data in corresponding time of one satellite from the data map
		*
		* @param[in]  sat			satellite name.
		* @param[in]  ics_start		begin time from ics data.
		* @param[in]  ics_end		end time from ics data.
		* @param[out] sat_param		corresponding satellite parameters data.
		* @return	if find the correct satellite parameters data.
		*	 @retval true	get data successfully
		*	 @retval false	get data failed
		*/
		bool getSatinfoUsed(std::string sat, t_gtime beg_time, t_gtime end_time, t_gsatparam& satinfo);

		/**
		* @brief get valid glonass frequency id.
		* @param[in]  sat			satellite name.
		* @param[in]  epoch			time from sat.
		* @param[out] fid			frequency id.
		*/
		int getGloFid(string sat,t_gtime epoch);

	protected:
		std::map<std::string, std::map<std::string, t_gsatparam> > _allsatinfo;		///< map for all satellites parameters data(prn,svn,record)	
	};

}
#endif // !GALLSATPARAM_H