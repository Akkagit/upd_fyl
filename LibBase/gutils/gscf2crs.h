/**
*
* @verbatim
	History
	 -1.0 Bo Wong	    2019-01-14 creat the file.
	 -1.1 Jie Li	    2019-03-31 add doxygen style code
	 -1.2 Yujie Qin     2019-07-04 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			    gscf2crs.h
* @brief		    GNSS satellite-fixed coordinate transform to CRS coordinate
*
* @author           Bo Wong, Wuhan University
* @version		    1.0.0
* @date			    2019-01-14
*
*/


#ifndef GROT_SCF2CRS_H
#define GROT_SCF2CRS_H

#include "newmat/newmat.h"
#include "gutils/ExportBase.h"
#include "gmodels/gpppmodel.h"

using namespace std;

namespace great
{
	/**
	*@brief	   Class for computing rotation from GNSS satellite-fixed system to CRS
	*/
	class BASE_LIBRARY_EXPORT t_gscf2crs
	{
	public:
		/** @brief default constructor. */
		t_gscf2crs();

		/** @brief default destructor. */
		virtual ~t_gscf2crs();

		/**
	    * @brief calculate rotation matrix of scf to crs
	    * @param[in]  xsat		position of satellite
	    * @param[in]  vsat		velocity of satellite
	    * @param[in]  xsun		position of satellite
	    * @param[out] rotmat	rotation matrix of scf to crs
	    */
		void calcRotMat(ColumnVector& xsat, ColumnVector& vsat, ColumnVector& xsun, Matrix& rotmat);                                       // main funtion , calculate rotation matrix scf to crs
		// extended by yqyuan
		void calcRotMat(string blocktype, string prn, ColumnVector& xsat, ColumnVector& vsat, ColumnVector& xsun, Matrix& rotmat);           // main funtion , calculate rotation matrix scf to crs
			   

	protected:
		/**
	    * @brief get number of block type(not used for now)
	    * @param[in]  type		block type
	    * @return   number
	    */
		int _blocktype_nam2num(string type);       
		
		/**
	    * @brief Remove string trailing spaces
	    * @param[in]  str		string need to be removed trailing space
	    * @return   string after removing trailing space
	    */
		string _trimR(string str);               
	};
}; // namespace

#endif