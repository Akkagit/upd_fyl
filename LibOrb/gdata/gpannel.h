/**
*
* @verbatim
    History
     -1.0 Yujie Qin	    2019-01-06 creat the file.
     -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
	 -1.2 Yujie Qin     2019-07-03 Adding Doxygen Style Code Remarks
     -1.3 Jiande Huang  2020-09-30 Standard code style and Optimize code logic
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gpannel.h
* @brief		Storage the LEO pannel data for one satellite
*
* @author       Yujie Qin, Wuhan University
* @version		1.0.0
* @date		    2019-03-27
*
*/

#ifndef GPANNEL_H
#define GPANNEL_H

#include <map>
#include <vector>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include "gdata/gdata.h"
#include "gutils/gtriple.h"

#include "gutils/ExportOI.h"

using namespace std;
using namespace gnut;

namespace great
{
	/**
	*@brief	   Class for one pannel of LEO pannel data
	*/
	class OI_LIBRARY_EXPORT t_oneside
	{
	public:
		/** @brief default constructor. */
		t_oneside();

		/** @brief default destructor. */
		virtual ~t_oneside();

		double       area;             ///< the area of one pannel side
		t_gtriple    scf_norm_dir;     ///< the norminal direction of the pannel side(x, y, z) in the SCF 面板侧的法线方向（x，y，z）
		ColumnVector crs_norm_dir;    ///< the norminal direction of the pannel side(x, y, z) in the CRS
		double       sr_coef   = 0.0;         ///< specular reflection coefficient of visible light 可见光的镜面反射系数
		double       dr_coef   = 0.0;         ///< diffuse reflection coefficient of visible light 可见光的漫反射系数
		double       isr_coef  = 0.0;         ///< specular reflection coefficient of infrared ray 红外线的镜面反射系数
		double       idr_coef  = 0.0;         ///< diffuse reflection coefficient of infrared ray 红外线的漫反射系数
		double       drag_coef = 0.0;         ///< drag coefficient  阻力系数
		string       name;                    ///< name of pannel
	};

	/**
	*@brief	   Class for all pannel data of one LEO
	*/
	class OI_LIBRARY_EXPORT t_gpannel : public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gpannel();

		/** @brief default destructor. */
		virtual ~t_gpannel();

		/**
	    * @brief push the LEO satellite name into the class
	    * @param[in]  name   LEO satellite name
	    */
		void add_name(string name);              

		/**
	    * @brief push the pannel index and the data of one side into the map.
	    * @param[in]  i       the index of pannel
	    * @param[in]  oneside the data of one pannel
	    */
		void add_pannel(int i, t_oneside oneside); 
		
		/**
		* @brief whether the pannel data is empty.
		* @return  bool
		@retval   true    pannel data is empty
		@retval   false   pannel data is existent
		*/
		bool isEmpty();

		/**
	    * @brief clear pannel data for one LEO satellite.
	    */
		void clear();
		
		/**
	    * @brief calculate the nordir in the CRS
	    * @param[in]  Rotmat  rotation matrix from TRS to CRS
		* @param[in]  unit_sat2sun unit vector from sat to Sun in CRS
	    */
		void get_rotmat(Matrix Rotmat, ColumnVector unit_sat2sun);
		
		/**
	    * @brief return all pannels information for one satellite
	    * @retrun   all pannels information for one satellite
	    */
		map<int, t_oneside>* get_pannels();           
		
		/**
	    * @brief return satellite name
	    * @return   LEO satellite name
	    */
		string get_sat_name();
	protected:
		string              _sat_name;  ///< satellite name
		map<int, t_oneside> _map_pannels;  ///< all pannel information for one satellite
	};
} //namespace

#endif