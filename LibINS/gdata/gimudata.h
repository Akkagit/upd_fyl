#ifndef GIMUDATA_H
#define GIMUDATA_H


/**
* @file gimudata.h
* @brief imu data structure for storing imu data
* @details Store data by vector

* @verbatim
	History
	-1.0	ShenZhiheng	      2018-12-15 creat the file.
  @endverbatim

* @author ShenZhiheng
* @version 1.0
* @date 2018-12-15
* @license
*/

#include <queue>
#include <stack>
#include "gdata/gdata.h"
#include "gins/gutility.h"
#include "gutils/ExportINS.h"
using namespace gnut;

namespace gsins
{

	/**
	*@brief t_gimudata Class for storing the IMU data
	*
	* t_gimudata is used for sins mechanical arrangement.
	* The data contains information about gyro and acce observation.
	* The gcoder t_gimudata corresponding to the gcoder t_imufile.
	*/
	class INS_LIBRARY_EXPORT t_gimudata : public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gimudata();

		/** @brief default destructor. */
		virtual ~t_gimudata() {}

		/**
		* @brief add one data into vector
		* @note data includes time and wm and vm
		* @param[in] 't' is time
		* @param[in] 'wm' is Angular increment
		* @param[in] 'vm' is Velocity increment
		* @return value reprents success or failure of add
		*/
		int add_IMU(const double& t, const Eigen::Vector3d& wm, const Eigen::Vector3d& vm);

		int add_IMU(const double& t, const Eigen::Vector3d& wm, const Eigen::Vector3d& vm, const Eigen::Vector3d& mm);

		/**
		* @brief set imu data interval
		* @param[in] 'ts' is interval
		* @return void
		*/
		void set_ts(double ts);

		void set_smoother(bool b) { _FBS = b; }

		/**
		* @brief load data into wm and vm and t accroding to subsamples 
		* @param[in] 'nSamples' is subsamples
		* @param[in] 'direction' is the direction
		* @param[out] 'wm' is Angular increment
		* @param[out] 'vm' is Velocity increment
		* @param[out] 'ts' is computing intervals
		* @return bool
		*/

		bool load(vector<Eigen::Vector3d>& wm, vector<Eigen::Vector3d>& vm, double& t, double& ts, int nSamples,bool _beg_end = true);
		bool load(vector<Eigen::Vector3d>& wm, vector<Eigen::Vector3d>& vm, vector<Eigen::Vector3d>& mm, double& t, double& ts, int nSamples, bool _beg_end = true);
		

		/**
		* @brief erase imu data before t
		*/
		t_gtime erase_bef(t_gtime t, bool _beg_end = true);

		/**
		* @brief return imu data vector size
		* @return int
		*/
		int size(bool _beg_end = true);


		virtual bool available(const t_gtime& now, bool _beg_end = true);

		/**
		* @struct dataIMU
		* @brief describe imu data information
		*/
		struct dataIMU 
		{ 
			double t; /**< time */
			Eigen::Vector3d wm, vm, mm;
		};

		virtual double beg_obs(bool _beg_end = true);  // get first epoch
		virtual double end_obs(bool _beg_end = true);  // get end epoch

	private:
		double _ts;					 /// imu interval 
		deque<dataIMU> _imu_back;    /// imu date stored for backward processing
		deque<dataIMU> _imu_forward;		 /// imu date stored for backward processing
		int first, end;				 /// imu starting and ending time 
		double _pre_time;			 /// previous epoch time 
		bool _FBS = false;
	};

}



#endif