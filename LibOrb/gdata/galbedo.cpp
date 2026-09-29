/**
*
* @verbatim
History
-1.0 Qian Zhang    2019-11-09 creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    galbedo.cpp
* @brief		    Storage the CERES data for one month
*
* @author       Qian Zhang, Wuhan University
* @version		1.0.0
* @date		    2019-11-09
*
*/

#include "galbedo.h"
using namespace std;

namespace great {
	/** @brief default constructor. */
	t_galbedo::t_galbedo()
		:t_gdata()
	{
		id_type(t_gdata::ALLALBEDO);
		_ceresr_data.resize(72, 144);
		_cerese_data.resize(72, 144);
		for (int i = 1; i <= 72; i++)
		{
			for (int j = 1; j <= 144; j++)
			{
				_ceresr_data(i, j) = 0.0;
				_cerese_data(i,j) = 0.0;
			}
		}
	}

	/**
	* @brief set one record with time and corresponding value(solar_flux/geomag_kp).
	* @param[in]    datatype      data type(solar/geomag)
	* @param[in]    ii            row index
	* @param[in]    jj            column index
	* @param[in]    value         value with the specified index
	*/
	void t_galbedo::set_data(string filetype, int ii, int jj, double value)
	{
		if (filetype == "refl")
		{
			_ceresr_data(ii, jj) = value;
		}
		else if (filetype == "emit")
		{
			_cerese_data(ii, jj) = value;
		}
	}

	/** @brief whether the CERES data is empty.
	* @return  bool
	*	@retval   true    CERES data is empty
	*   @retval   false   CERES data is existent
	*/
	bool t_galbedo::is_empty()
	{
		if (_ceresr_data.is_zero() || _cerese_data.is_zero())
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}