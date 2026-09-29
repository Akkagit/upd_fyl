
#ifndef GPUBLISH_H
#define GPUBLISH_H

/**
* @file gpublish.h
* @brief visualization algorithm
* @details
		VIS define the t_gpublish class for trajectory display
* @verbatim
		History
		-1.0	LiShengyu	       2020-05-20 create the file.
  @endverbatim

* @author LiShengyu
* @email
* @version 1.0
* @date 2020-05-20
* @license
*/

#include "gviewer.h"
#include "gutility.h"
#include "gutils/ExportVIS.h"
#include <Eigen/Eigen>
#include <fstream>

using namespace std;

namespace gvis
{

	/**
	* @class t_gpublish
	* @brief t_gpublish Class is used for trajectory display
	*/
	class VIS_LIBRARY_EXPORT t_gpublish
	{
	public:
		/** @brief default constructor. */
		t_gpublish();

		/** @brief default destructor. */
		~t_gpublish();

		/**
		* @brief initialize the viewer for display
		*/
		void Initialize();

		/**
		* @brief update the parameters for trajectory display
		*
		* @param[in] imu_state		store the imu information to show on the viewer
		* @param[in] cam_states		store the cam information to show on the viewer
		* @param[in] init_campos	the initial value of the camera position
		*/
		void UpdateNewState(const IMUState &imu_state,const CamStateServer &cam_states,const Eigen::Vector3d & init_campos);
	private:
		
		gviewer *viewer;				// the core class of trajectory display
		bool firstFlag_imu = true;		// A indicator to determine the value of _init_imupos
		bool firstFlag_cam = true;		// A indicator to determine the value of _init_campos,
										// maybe imgdata is later than imudata
		Eigen::Vector3d _init_campos;	// transfer the zero of the IMU coordinate
		Eigen::Vector3d _init_imupos;	// transfer the zero of the IMU coordinate
		Eigen::Matrix3d R_e_n;          // rotation from ECEF frame to ENU frame for first imudata

		
	};
}


#endif
