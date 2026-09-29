#ifndef GLIDARDATA_H
#define GLIDARDATA_H

/**
* @file glidardata.h
* @brief lidar data structure for storing lidar data
* @details Store data by vector

* @verbatim
	History
	-1.0	WangHuidan	      2020-10-08 create the file.
  @endverbatim

* @author WangHuidan
* @version 1.0
* @date 2020-10-08
* @license
*/

#include "gdata/gdata.h"
#include "glidar/gutility.h"
#include "gutils/ExportLidar.h"

using namespace gnut;

namespace glidar
{
	/**
	*@brief t_glidardata Class for storing the LIDAR data,just the timestamp and lidarpath
	*
	* lidarfile is used for Front end icp register.
	* The data contains information about  timestamp and lidarpath.
	* The gcoder t_glidardata corresponding to the gcoder t_lidarfile.
	*/
	class LIDAR_LIBRARY_EXPORT t_glidardata : public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_glidardata();
		/** @brief default destructor. */
		virtual ~t_glidardata() {};
	
		/**
		* @brief add one data into vector
		* @note data includes time and img0_path and img1_path
		* @param[in] 't' is time
		* @param[in] 'lidar_path' is lidar data path
		* @return value reprents success or failure of add
		*/
		int add_LIDAR(const double& t,const string& lidar_path);

		/**
		* @brief load data into lidar_path according to imu_t
		* @param[in] imu_t
		* @param[out] lidar_t
		* @param[out] lidar_path
		* @return void
		*/
		void load(double& imu_t, double& lidar_t, lidarPath& lidar_path);
		
		/** 
		* @brief return lidarpath data vector size
		* @note not use
		* @return int
		*/
		int size();
	
		private:
			vector<lidarPath> _veclidar;	//lidarpath data
	};

}
#endif