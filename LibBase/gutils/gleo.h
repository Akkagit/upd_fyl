/**
*
* @verbatim
    History
     -1.0 Wei Zhang     2019-01-18 creat the file.
     -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
	 -1.2 Yujie Qin     2019-07-03 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gleo.h
* @brief		all supported LEO
*
* @author		Wei Zhang, Wuhan University
* @version		1.0.0
* @date			2019-03-27
*
*/

#ifndef GLEO_H
#define GLEO_H

#include <map>
#include <vector>
#include <string>
#include <set>
#include "gutils/ExportBase.h"

using namespace std;

namespace great
{
	/** @brief full name for all supported LEO satellites. */
	enum GLEOLONG
	{
		CHAMP,
		//GNOS,
		HY_2A,
		TIANGONG_2,
		ZIYUAN_3,
		KOMPSAT5,
		LINGQIAO,
		FY_3C, FY_3D,
		JASON_2, JASON_3,
		GRACE_A, GRACE_B,
		GRACE_C, GRACE_D,
		METOP_A, METOP_B, METOP_C,
		SWARM_A, SWARM_B, SWARM_C,
		COSMIC_1, COSMIC_2, COSMIC_3, COSMIC_4, COSMIC_5, COSMIC_6,
		SENTINEL_1A, SENTINEL_1B, SENTINEL_2A, SENTINEL_2B, SENTINEL_3A, SENTINEL_3B,
		TERRASAR_X, TANDEM_X,
		PAZSAT,
		UNDEF
	};

	/** @brief abbreviation for all supported LEO satellites(four character, just like a station name) */
	enum GLEOSHORT
	{
		CHAM,
		//GNOS,
		HY2A,
		TG02,
		ZY03,
		KOM5,
		LING,
		TERX,
		FY3C, FY3D,
		JAS2, JAS3,
		GRAA, GRAB,
		GRAC, GRAD,
		META, METB, METC,
		SWAA, SWAB, SWAC,
		COS1, COS2, COS3, COS4, COS5, COS6,
		SE1A, SE1B, SE2A, SE2B, SE3A, SE3B,
		TESX, TADX,
		PAZS,
		UNDF
	};

	// ---------------------------------
	// global funcitons
	// ---------------------------------
	
	/**
	* @brief convert lower string to upper string of LEO satellite name
	* @param[in]  s  lower string name of satellite
	* @return     the upper string name of satellite
	*/
	BASE_LIBRARY_EXPORT std::string low2upp(std::string s);
	
	/**
	* @brief convert lower string to upper string of LEO satellite name
	* @param[in]  s  upper string name of satellite
	* @return     the lower string name of satellite
	*/
	BASE_LIBRARY_EXPORT string upp2low(string s);

	/**
	* @brief convert long string to GLEOLONG enum
	* @parma[in]  s   long string name of satellite
	* @return    the GLEOLONG name of satellite
	*/
	BASE_LIBRARY_EXPORT GLEOLONG   str2long(string s);
	
	/**
	* @brief convert GLEOSHORT enum to GLEOLONG
	* @param[in]  leo   GLEOSHORT name of satellite
	* @return      GLEOLONG name of satellite
	*/
	BASE_LIBRARY_EXPORT GLEOLONG   short2long(GLEOSHORT leo);
	
	/**
	* @brief convert short string to GLEOSHORT
	* @param[in]  s  short string name of satellite
	* @return     GLEOSHORT name of satellite
	*/
	BASE_LIBRARY_EXPORT GLEOSHORT  str2short(string s);
	
	/**
	* @brief convert GLEOLONG enum to GLEOSHORT
	* @param[in]  leo  GLEOLONG name of satellite
	* @return    GLEOSHORT name of satellite
	*/
	BASE_LIBRARY_EXPORT GLEOSHORT  long2short(GLEOLONG leo);

	/**
	* @brief convert GLEOLONG enum to long string
	* @param[in]  leo  GLEOLONG name of satellite
	* @return   long string name of satellite
	*/
	BASE_LIBRARY_EXPORT string long2str(GLEOLONG leo);
	
	/**
	* @brief convert GLEOSHORT enum to short string
	* @param[in]  leo  GLEOSHORT name of satellite
	* @return      short string name of satellite
	*/
	BASE_LIBRARY_EXPORT string short2str(GLEOSHORT leo);
	
	/**
	* @brief convert short string to long string
	* @param[in]  s  short string name of satellite
	* @return   long string name of satellite
	*/
	BASE_LIBRARY_EXPORT string strs2strl(string s);
	
	/**
	* @brief convert long string to short string
	* param[in]  s  long string name of satellite
	* @return  short string name of satellite
	*/
	BASE_LIBRARY_EXPORT string strl2strs(string s);

	/**
	* @brief get all supported LEO satellites
	* @return  all supported LEO satellites(GLEOLONG name)
	*/
	BASE_LIBRARY_EXPORT set<GLEOLONG> LEO_SUPPORTED();
	
	/**
	* @brief judge the LEO satellite is supported LEO or not
	* @param[in]  s  long string name of satellite
	* @return   bool
	*	@retval   true   the satellite is supported LEO
	*   @retval   false   the satellite isnot supported LEO
	*/
	BASE_LIBRARY_EXPORT bool IS_SUPPORTED_LEO(std::string s);

} // namespace

#endif