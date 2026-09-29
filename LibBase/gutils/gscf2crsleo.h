/**
*
* @verbatim
    History
     -1.0 Wei Zhang     2019-01-20 creat the file.
     -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
	 -1.2 Yujie Qin     2019-07-04 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gscf2crsleo.h
* @brief		calculate the rotation matrix from SCF to CRS
*
* @author		Wei Zhang, Wuhan University
* @version		1.0.0
* @date			2019-03-27
*
*/


#ifndef GSCF2CRSLEO_H
#define GSCF2CRSLEO_H

#include <map>

#include "newmat/newmat.h"
#include "newmat/newmatio.h"
#include "gall/gallattitude.h"
#include "gdata/gattitude.h"
#include "gutils/gtime.h"
#include "gutils/gscf2crs.h"
using namespace std;

namespace great
{
	/**
	*@brief	   Class for rotation matrix from SCF to CRS
	*/
	class BASE_LIBRARY_EXPORT t_gscf2crsleo : public t_gscf2crs
	{
	public:
		/** @brief default constructor. */
		t_gscf2crsleo(string satname);

		/** @brief default destructor. */
		virtual ~t_gscf2crsleo();

		/** 
		* @brief get the time of matrix 
		* @return  epoch of the rot Matrix
		*/
		t_gtime getCurtEpoch();

		/**
		* @brief get the name of satname 
		* @return the sat name
		*/
		string getSatName();

		/**
	    * @brief main function for rotation calculation
	    * @param[in]  time    the current epoch
	    * @param[in]  attdata the attitude data for the LEO satellite
	    * @param[in]  pos     the position of satellite in CRS
	    * @param[in]  vel     the velocity of satellite in CRS
	    * @param[in]  trs2crs the rotation matrix from TRS to CRS
	    * @param[out] scf2crs the rotaiton matrix from SCF to CRS
	    * @return   
		*	@retval>0 calculate the rotation matrix correctly
		*	@retval<0 calculate the rotation matrix not correctly
	    */
		int calcRotMat(t_gtime time, t_gattitude* attdata,const ColumnVector& pos,const ColumnVector& vel, Matrix& trs2crs, Matrix& scf2crs);  // main function

	protected:
		
		/**
	    * @brief calculate the nominal rotation matrix from SCF to CRS(without attitude data)
	    * @param[in]  pos  the position of satellite in CRS
	    * @param[in]  vel  the velocity of satellite in CRS
	    * @return  
		*	@retval>0 calculate the rotation matrix correctly
		*	@retval<0 calculate the rotation matrix not correctly 
	    */
		int  _calcNomMat(ColumnVector pos, ColumnVector vel); 

		/**
	    * @brief calculate the true rotation matrix from SCF to CRS(with attitude data)
	    * @param[in]  time     the current epoch
	    * @param[in]  attdata  the attitude data for the LEO satellite
	    * @return   
		*	@retval>0 calculate the rotation matrix correctly
		*	@retval<0 calculate the rotation matrix not correctly 
	    */
		int  _calcAttMat(t_gtime time, t_gattitude* attdata);    

		// variables
		string     _sat_name;       ///< LEO satellite name
		bool       _isAttValid;     ///< the satellite has attitude data or not
		t_gtime	   _epoch;	        ///< epoch of the rot Matrix

		Matrix     _rotMatScf2crs;   ///< rotation matrix from SCF to CRS for LEO satellite(Nominal attitude)
		Matrix     _rotMat;          ///< the final rotation matrix
	};
} //namespace

#endif