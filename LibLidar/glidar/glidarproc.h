#ifndef GLIDARPROC_H
#define GLIDARPROC_H


/**
* @file glidarproc.h
* @details 
		define some common functions and all parameters that may be used in subclasses
* @endverbatim
* @verbatim
	History
	-1.0	WangHuidan	      2020-10-07 create the file.
  @endverbatim

* @author WangHuidan
* @version 1.0
* @date 2020-10-07
* @license
*/

#include <Eigen/Eigen>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include "glidar/gutility.h"
#include "gutils/gstring.h"
#include "gutils/ExportLIDAR.h"
#include "gset/gsetbase.h"


using namespace std;
using namespace gnut;

namespace glidar {
	/**
	* @class t_lidarproc
	* @brief t_lidarproc Class for reading parameters and virtual functions
	*/
	class LIDAR_LIBRARY_EXPORT t_glidarproc
	{
	public:

		/**
		* @brief constructor.
		*
		* @param[in]  _set     setbase control
		*/
		t_glidarproc(t_gsetbase* _set);

		/** @brief destructor. */
		~t_glidarproc() {}

		/**
		* @brief save imu observation in _vecimu
		* all imu observations' timestamps are between cur frame and pre frame
		*
		* @param[in] t		 imu observation time
		* @param[in] gyro	 corresponding angular velocity observation
		* @param[in] acc     corresponding acceleration observation
		*/
		void load_imuobs(const double &t, const vector<Eigen::Vector3d> &gyro, const vector<Eigen::Vector3d> &acc);

		/**
		* @brief get cur lidar observation
		* @param[in] t			 lidar observation time
		* @param[in] lidar_path	 corresponding lidar path
		*/
		void load_lidarobs(const double &t,const lidarPath& lidar_path);

		/**
		* @brief lidar processing main function
		*
		* @param[in] t		 cur lidar time
		*/
		virtual LidarFrame ProcessBatch(const double& t);


		/**
		* @struct imumsg
		* @brief  imudata vector only used in glidarproc class
		*/
		struct imumsg
		{
			double t;
			Eigen::Vector3d linear_acceleration;
			Eigen::Vector3d angular_velocity;
		};

	public:
		//new add;
		lidarPath cur_lidar_path;

		vector<imumsg> _vecimu;			// vector store imu data between pre frame and cur frame
		dataLIDAR _pointCloud;			// vector store pointcloude

		double ts;						// camera sample interval. //
		int freq;						// camera sample freq. //
		double imu_ts;					// imu sample interval. //
		int imu_freq;					// imu sample freq //

		bool use_scan, use_map;
		bool proj_lastmap;
		int map_gap;

		/* Basic parameters */
		Eigen::Matrix3d R_lidar_imu;
		Eigen::Vector3d t_lidar_imu;
		Eigen::Isometry3d T_lidar_imu;
		/* Basic parameters */
		

	};
}

#endif