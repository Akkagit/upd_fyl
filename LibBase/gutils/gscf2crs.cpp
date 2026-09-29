/**
*
* @verbatim
	History
	 -1.0 Bo Wong	2019-01-14 creat the file.
	 -1.1 Jie Li	2019-03-31 add doxygen style code
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gscf2crs.cpp
* @brief		GNSS satellite-fixed coordinate transform to CRS coordinate
*
* @author       Bo Wong, Wuhan University
* @version		1.0.0
* @date			2019-01-14
*
*/

#include "gutils/gscf2crs.h"
#include <math.h>
using namespace std;
using namespace gnut;
namespace great
{
	/** @brief default constructor. */
	t_gscf2crs::t_gscf2crs()
	{

	}

	/** @brief default destructor. */
	t_gscf2crs::~t_gscf2crs()
	{

	}

	/**
	* @brief calculate rotation matrix of scf to crs
	* @param[in]  xsat		position of satellite
	* @param[in]  vsat		velocity of satellite
	* @param[in]  xsun		position of satellite
	* @param[out] rotmat	rotation matrix of scf to crs
	*/
	void t_gscf2crs::calcRotMat(ColumnVector& xsat, ColumnVector& vsat, ColumnVector& xsun, Matrix& rotmat)
	{
		double fx, fz, u;

		ColumnVector scf_X(0.0, 3);
		ColumnVector scf_Y(0.0, 3);
		ColumnVector scf_Z(0.0, 3);

		ColumnVector unit_vsat(0.0, 3);
		ColumnVector unit_sc2sun(0.0, 3);
		ColumnVector unit_orbnorm(0.0, 3);
		ColumnVector unit_xsun(0.0, 3);
		ColumnVector unit_sunnorm(0.0, 3);
		ColumnVector unit_node(0.0, 3);

		// sc-fixed z-axis, from sc to earth
		scf_Z = -xsat * 1e3;
		scf_Z = scf_Z / scf_Z.norm_Frobenius();

		// unint vector from sc to sun
		unit_sc2sun = xsun * 1e3 - xsat * 1e3;
		unit_sc2sun = unit_sc2sun / unit_sc2sun.norm_Frobenius();

		// unit vector of orbit plane
		unit_vsat = vsat / vsat.NormFrobenius();

		// dx1
		unit_orbnorm = crossproduct(-scf_Z, unit_vsat);
		unit_orbnorm = unit_orbnorm / unit_orbnorm.norm_Frobenius();

		// beta angel between orbit plane and sun (not used)
		unit_xsun = xsun / xsun.norm_Frobenius();

		// u angle from midnight
		unit_sunnorm = crossproduct(unit_xsun, unit_orbnorm);
		unit_sunnorm = unit_sunnorm / unit_sunnorm.norm_Frobenius();
		unit_node = crossproduct(unit_sunnorm, unit_orbnorm);
		unit_node = unit_node / unit_node.norm_Frobenius();
		fx = dotproduct(unit_node, -scf_Z);
		fz = dotproduct(crossproduct(unit_node, -scf_Z), unit_orbnorm);
		u = atan2(fz, fx);

		// in shadow
		scf_Y = crossproduct(scf_Z, unit_sc2sun);

		// for IGS antex-file not necessary. 20160607
		scf_Y = scf_Y / scf_Y.norm_Frobenius();

		// the sc-fixed z-axis
		scf_X = crossproduct(scf_Y, scf_Z);

		//rotation-matrix
		for (int i = 1; i < 4; i++)
		{
			rotmat(i,1 ) = scf_X(i);
			rotmat(i,2 ) = scf_Y(i);
			rotmat(i,3 ) = scf_Z(i);
		}

		/*if (sat[0] == 'G' || sat[0] == 'R')
		{
			
		}*/
		return;
	}
	
	void t_gscf2crs::calcRotMat(string blocktype, string prn, ColumnVector& xsat, ColumnVector& vsat, ColumnVector& xsun, Matrix& rotmat)
	{
		ColumnVector scf_X(0.0, 3);
		ColumnVector scf_Y(0.0, 3);
		ColumnVector scf_Z(0.0, 3);

		gnut::t_gpppmodel* nav_attitude = new gnut::t_gpppmodel();

		nav_attitude->attitude(blocktype, prn, xsat, vsat, xsun, scf_X, scf_Y, scf_Z);

		//rotation-matrix
		for (int i = 1; i < 4; i++)
		{
			rotmat(i, 1) = scf_X(i);
			rotmat(i, 2) = scf_Y(i);
			rotmat(i, 3) = scf_Z(i);
		}

		delete nav_attitude; 
		nav_attitude = NULL;

		return;
	}



	/**
	* @brief get number of block type(not used for now)
	* @param[in]  type		block type
	* @return   number
	*/
	int t_gscf2crs::_blocktype_nam2num(string type)
	{
		string case_tmp = _trimR(type);
		int iblock;
		//GPS
		if (case_tmp == "BLOCK I")
		{
			iblock = 1;
		}
		else if (case_tmp == "BLOCK II")
		{
			iblock = 2;
		}
		else if (case_tmp == "BLOCK IIA")
		{
			iblock = 3;
		}
		else if (case_tmp == "BLOCK IIR")
		{
			iblock = 4;
		}
		else if (case_tmp == "BLOCK IIR-A")
		{
			iblock = 5;
		}
		else if (case_tmp == "BLOCK IIR-B")
		{
			iblock = 6;
		}
		else if (case_tmp == "BLOCK IIR-M")
		{
			iblock = 7;
		}
		else if (case_tmp == "BLOCK IIF")
		{
			iblock = 8;
		}

		//GLONSS
		else if (case_tmp == "GLONASS")
		{
			iblock = 101;
		}
		else if (case_tmp == "GLONASS-M")
		{
			iblock = 102;
		}
		else if (case_tmp == "GLONASS-K" || case_tmp == "GLONASS-K1")
		{
			iblock = 103;
		}

		//Galileo
		else if (case_tmp == "GALILEO GIOVEA" || case_tmp == "GALILEO-0A")
		{
			iblock = 201;
		}
		else if (case_tmp == "GALILEO GIOVEB" || case_tmp == "GALILEO-0B")
		{
			iblock = 201;
		}
		else if (case_tmp == "GALILEO IOV" || case_tmp == "GALILEO-1")
		{
			iblock = 202;
		}
		else if (case_tmp == "GALILEO FOC" || case_tmp == "GALILEO-2")
		{
			iblock = 203;
		}
		//BDS
		else if (case_tmp == "BEIDOU-2G")
		{
			iblock = 301;
		}
		else if (case_tmp == "BEIDOU-2I")
		{
			iblock = 302;
		}
		else if (case_tmp == "BEIDOU-2M")
		{
			iblock = 303;
		}
		//QZSS
		else if (case_tmp == "QZSS")
		{
			iblock = 401;
		}
		else
		{
			return -1;
		}
		return iblock;
	}
	
	// jdhuang : remove the warning, make some change
	/**
	* @brief Remove string trailing spaces
	* @param[in]  str		string need to be removed trailing space
	* @return   string after removing trailing space
	*/
	string t_gscf2crs::_trimR(string str)
	{
		int nlen = str.length();
		for (int i = nlen - 1; i >= 0; i--)
		{
			if (!isspace(str[i]))
			{
				return str.substr(0, i + 1);
			}
		}
		return str;
	}

}// namespace
