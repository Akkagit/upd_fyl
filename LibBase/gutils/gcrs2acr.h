/**
*
* @verbatim
	History
	 -1.0 Bo Wong	    2019-01-04 creat the file.
	 -1.1 Jie Li		2019-03-31 add doxygen style code
	 -1.2 Yujie Qin     2019-07-04 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			  gcrs2acr.h
* @brief		  CRS coordinate transform to ACR coordinate
*
* @author         Bo Wong, Wuhan University
* @version		  1.0.0
* @date			  2019-01-04
*
*/


#ifndef GCRS2ARC_H
#define GCRS2ARC_H

#include "newmat/newmat.h"
#include "gutils/ExportBase.h"
using namespace std;

namespace great
{
	/**
	*@brief	   Class for calculating CRS to ACR rotation matrix
	*/
	class BASE_LIBRARY_EXPORT t_gcrs2acr
	{
	public:
		/** @brief default constructor. */
		t_gcrs2acr();		

		/** @brief default destructor. */
		virtual ~t_gcrs2acr();

		//todo º¯Êý·µ»ØÖµ£¿
		/**
	    * @brief calculate rotation matrix of CRS to ACR
	    * @param[in]  xsat		position of satellite
	    * @param[in]  vsat		velocity of satellite
	    * @param[out] rotmat	rotation matrix of CRS to ACR
	    */
		void calcRotMat(ColumnVector& xsat, ColumnVector& vsat, Matrix& rotmat);

	protected:


	};
}
#endif