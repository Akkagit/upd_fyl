/**
*
* @verbatim
    History
     -1.0 Yujie Qin     2019-01-06 creat the file.
     -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gpannel.cpp
* @brief		    Storage the LEO pannel data for one satellite
*
* @author       Yujie Qin, Wuhan University
* @version		1.0.0
* @date		    2019-03-27
*
*/

#include "gdata/gpannel.h"
#include "gall/gallpannel.h"

using namespace std;
namespace great {

	/** @brief default constructor. */
	t_oneside::t_oneside()
	{
		area = 0.0;
		crs_norm_dir.resize(3);
		crs_norm_dir = 0;
		sr_coef = 0.0;
		dr_coef = 0.0;
		drag_coef = 0.0;
		isr_coef = 0.0;
		idr_coef = 0.0;
	}

	/** @brief default destructor. */
	t_oneside::~t_oneside()
	{

	}

	/** @brief default constructor. */
	t_gpannel::t_gpannel()
	{

	}

	/** @brief default destructor. */
	t_gpannel::~t_gpannel()
	{}

	/**
	* @brief push the LEO satellite name into the class
	* @param[in]  name   LEO satellite name
	*/
	void t_gpannel::add_name(string name)
	{
		_sat_name = name;
	}

	/**
	* @brief push the pannel index and the data of one side into the map.
	* @param[in]  i       the index of pannel
	* @param[in]  oneside the data of one pannel
	*/
	void t_gpannel::add_pannel(int i, t_oneside oneside)
	{
		_map_pannels[i] = oneside;
	}

	bool t_gpannel::isEmpty()
	{
		if (_map_pannels.size() == 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	/**
	* @brief return all pannels information for one satellite
	* @retrun   all pannels information for one satellite
	*/
	map<int, t_oneside> *t_gpannel::get_pannels()
	{
		return &_map_pannels;
	}

	/**
	* @brief return satellite name
	* @return   LEO satellite name
	*/
	string t_gpannel::get_sat_name()
	{
		return _sat_name;
	}

	/**
	* @brief clear pannel data for one LEO satellite.
	*/
	void t_gpannel::clear()
	{
		_sat_name = "";
		_map_pannels.clear();
	}

	/**
	* @brief calculate the nordir in the CRS
	* @param[in]  Rotmat  rotation matrix from TRS to CRS
	* @param[in]  unit_sat2sun unit vector from sat to Sun in CRS
	*/
	void t_gpannel::get_rotmat(Matrix Rotmat, ColumnVector unit_sat2sun)
	{
		map<int, t_oneside>::iterator it;
		try
		{
			for (it = _map_pannels.begin(); it != _map_pannels.end(); it++)
			{
				Matrix norm;
				norm = it->second.scf_norm_dir.crd_cvect();
				norm.resize_keep(3, 1);
				//it->second.normdirj = Rotmat * norm;
				double cosa = dotproduct(norm, norm);
				// the front of solar pan is pointing to the sun
				if (it->second.name.substr(0, 10) == "Solar_PAN+" && cosa < 1e-10) {
					it->second.crs_norm_dir = unit_sat2sun;
				}
				else if (it->second.name.substr(0, 10) == "Solar_PAN-" && cosa < 1e-10) {
					it->second.crs_norm_dir = -1.0 * unit_sat2sun;
				}
				else {
					it->second.crs_norm_dir = Rotmat * norm;
				}
			}
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_gpannel::rotmat throw exception");
			}
			else
			{
				cout << "ERROR : t_gpannel::rotmat throw exception" << endl;
			}
			return;
		}
	}
}