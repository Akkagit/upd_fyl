/**
*
* @verbatim
History
-1.0 Qian Zhang		 2019-11-09 creat the file.
-1.1 Jiande Huang    2020-09-12 Standard code style and Optimize code logic
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    galbedo.h
* @brief		Storage the CERES data for one month
*
* @author       Qian Zhang, Wuhan University
* @version		1.0.0
* @date		    2019-11-09
*
*/

#ifndef GALBEDO_H
#define GALBEDO_H

#include <vector>
#include <map>
#include <string.h>
#include "gdata/gdata.h"
#include "gutils/gtime.h"
#include "gutils/ExportOI.h"
#include "newmat/newmat.h"
#include "newmat/newmatio.h"

using namespace std;
using namespace gnut;
namespace great{
	/**
	*@brief	   Class for CERES data
	*/
	class OI_LIBRARY_EXPORT t_galbedo : public t_gdata {

	public:
		/** @brief default constructor. */
		t_galbedo();

		/** @brief default destructor. */
		virtual ~t_galbedo(){};

		/**
		* @brief set one record with time and corresponding value(solar_flux/geomag_kp).
		* @param[in]    filetype      file type(solar/geomag)
		* @param[in]    ii            row index
		* @param[in]    jj            column index
		* @param[in]    value         value in the specified time
		*/
		void set_data(string file_type, int ii, int jj, double value);

		/**
		* @brief whether the ceres data is empty.
		* @return  bool
		*	@retval   true    ceres data is empty
		*   @retval   false   ceres data is existent
		*/
		bool is_empty();

		/**
		* @brief return refl data.
		* @return  refl data
		*/
		Matrix& get_ceresr_data() { return _ceresr_data; };

		/**
		* @brief return emit data.
		* @return  emit data
		*/
		Matrix& get_cerese_data() { return _cerese_data; };

	protected:
		Matrix _ceresr_data;      ///< refl data
		Matrix _cerese_data;      ///< emit data
	};
}

#endif