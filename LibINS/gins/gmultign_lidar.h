#ifdef USE_LIDAR
#ifndef GMULTIGN_LIDAR_H
#define GMULTIGN_LIDAR_H

/**
* @file gmultign_lidar.h
* @brief ppp/ins/lidar tightly coupled integration key algorithm
* @note
* @verbatim
		History
		-1.0	WangHuidan	       2020-06-18 creat the file.
  @endverbatim

* @author WangHuidan
* @date 2020-10-12
* @license
*/

#include "gins/gins.h"
#include "gins/gintegration.h"
#include "gutils/gmutex.h"
#include "gset/gsetbase.h"
#include "gdata/gimudata.h"
#include "gdata/glidardata.h"
#include "gio/glog.h"
#include "gio/giof.h"
#include "gset/gsetign.h"
#include "gset/gcfg_multign.h"
#include "gproc/gpvtflt.h"
#include "gdata/gposdata.h"
#include "gutils/ExportLIDAR.h"
#include "glidar/glidar.h"
using namespace gnut;
using namespace glidar;


namespace gsins
{

	/**
	* @class t_gmulting
	*
	* t_gmulting Class for ppp/ins/vis tightly integration
	*
	* loosely-coupled integration based on Adaptive and Sequential kalman filter
	* Inherited from t_gintegration
	*
	*/
	class INS_LIBRARY_EXPORT t_gmultign_lidar :virtual public t_gintegration,virtual public t_glidar
	{
	public:
		/**
		* @brief Constructor
		* @note set parameter value by gset
		*/
		t_gmultign_lidar(string site, string site_base, t_gsetbase* gset, t_gallproc *allproc);

		/**
		* @brief set whether t_gallpar is used or not
		*/
		void setgallparuse(bool status);

		/**
		* @brief struct store lidar data
		*/
		virtual void Add_LIDAR(t_glidardata* lidar);
	    
		/**
		* @brief determines whether the interval between two timestamps is less than the delay
		*/
		virtual bool _lidar_time_valid(t_gtime lidart, t_gtime inst);

		/** @brief processBatch. */
		virtual int processBatch(const t_gtime& beg, const t_gtime& end);

		/** @brief process one epoch.	*/
		virtual int _processGNSSEpoch(const t_gtime& runEpoch);
	
		void removeLidarState(vector<StateIDType>& rm_lidar_state_ids);

		/**
		* @brief process lidar observation
		*
		* @param[in] runEpoch			lidar time
		*/
		int _processLidarEpoch(const t_gtime& runEpoch);

		/**
		* @brief process cam observation
		*
		* @param[in] dx		correction of params	
		*/
		virtual int _lidar_feedback(const ColumnVector & dx);

		/**
		* @brief process cam observation
		*
		* @param[in] dx		correction of params
		*/
		virtual int _lidar_feedback(const Eigen::VectorXd & dx);

		/**
		* @brief StateIDType to string
		*
		* @param[in] id			lidar state id
		*/
		string state_id2str(const StateIDType& id);
	

		/**
		* @brief perform cam measurement update
		*
		* @param[in] H			measurement matrix vector
		* @param[in] r			residual vector
		*/
		void lidar_meas_update(const Eigen::MatrixXd& H, const Eigen::VectorXd& r);



		
		/**
		* @brief state augment for a new cam state and complete time update of cov 
		*  between cam and imu
		*
		* @param[in] runEpoch		cam time
		*/
		void stateAugmentation(const t_gtime& runEpoch);


		/**  !!!functions below are not used in PPP/INS/VISION tightly coupled£¡£¡£¡        */

		/**
		* @brief put covariance into Pk PIC and PCC
		*/
		void assign_statecov(const Eigen::MatrixXd &P);

		/**
		* @brief get covariance from Pk PIC and PCC
		*/
		Eigen::MatrixXd get_statecov();
		
		/**
		* @brief mostly used for gnss update in STC and LC mode
		*/
		int _mult_meas_update(double fading = 1.0);

		// this functions just for test
		//int test_processBatch(const t_gtime& beg, const t_gtime& end, bool prtOut);
		//virtual int test_processGNSSEpoch(const t_gtime& runEpoch);
		//int test_processCamEpoch(const t_gtime& runEpoch);
		//int test_cam_feedback();
		//void test_stateAugmentation(const t_gtime& runEpoch);
		
		//void test_assign_statecov(const Eigen::MatrixXd &P);
		//void test_removeLostFeatures();
		//bool test_gatingTest(const Eigen::MatrixXd& H, const Eigen::VectorXd& r, const int& dof);
		//void test_cam_meas_update(const Eigen::MatrixXd& H, const Eigen::VectorXd& r);
		//void test_pruneCamStateBuffer();

	protected:
		t_gtime _lidar_crt;					/// record cur image time
		t_glidardata* _lidardata;

		int _framecount;
		//int _scanPeriod;

		/* Make the structure clear and reduce the temporary definition of intermediate variables in the structure  */
		LidarFrame _lidarframe;             /// point cloud result,from font-end
		LidarFrame _lastframe;				/// point cloud of last lidar frame
		bool mIsFirstLidar;			/// A indicator to daetermine if it is the first image 

		//whd+ test
		bool pos_jump=false;

									        /// Maintain numerical and corresponding derivative stability
		double _lidar_t;						/// get cur image time
		lidarPath _lidar_path;					/// get cur image path
		t_glidarMapping _lidarmap;
		t_glidarOdometry _lidarOdo;

	private:
		bool is_allTC = false;			// whether t_gallpar is used or not
	};

}




#endif

#endif // USE_LIDAR