/**
*
* @verbatim
	History
	 -1.0 Bo Wong	2019-01-04 creat the file.
	 -1.1 Jie Li	2019-04-02 add doxygen style code
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gcrs2acr.cpp
* @brief		CRS coordinate transform to ACR coordinate
*
* @author       Bo Wong, Wuhan University
* @version		1.0.0
* @date			2019-01-04
*
*/

#include "gutils/gcrs2acr.h"

using namespace std;
namespace great
{
	/** @brief default constructor. */
	t_gcrs2acr::t_gcrs2acr()
	{

	}

	/** @brief default destructor. */
	t_gcrs2acr::~t_gcrs2acr()
	{

	}

	/**
	* @brief calculate rotation matrix of CRS to ACR
	* @param[in]  xsat		position of satellite
	* @param[in]  vsat		velocity of satellite
	* @param[out] rotmat	rotation matrix of CRS to ACR
	*/
	void t_gcrs2acr::calcRotMat(ColumnVector& xsat, ColumnVector& vsat, Matrix& rotmat)
	{
		ColumnVector radial(0.0, 3);
		ColumnVector along(0.0, 3);
		ColumnVector cross(0.0, 3);

		radial = xsat / xsat.NormFrobenius();
		along  = vsat / vsat.norm_Frobenius();

		cross = crossproduct(radial, along);
		cross = cross / cross.norm_Frobenius();

		along = crossproduct(cross, radial);
		along = along / along.norm_Frobenius();

		for (int i = 1; i <= 3; i++)
		{
			rotmat(1, i) = along(i);
			rotmat(2, i) = cross(i);
			rotmat(3, i) = radial(i);
		}
	}

}

