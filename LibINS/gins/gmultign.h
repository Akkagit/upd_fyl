#ifdef USE_VIS
#ifndef GMULTIGN_H
#define GMULTIGN_H

/**
* @file gmultign.h
* @brief ppp/ins/vis tightly coupled integration key algorithm
* @note
* @verbatim
		History
		-1.0	LiShengyu	       2020-06-18 creat the file.
  @endverbatim

* @author LiShengyu
* @date 2020-06-18
* @license
*/

#include "gins/gins.h"
#include "gins/gintegration.h"
#include "gutils/gmutex.h"
#include "gset/gsetbase.h"
#include "gdata/gimudata.h"
#include "gdata/gimgdata.h"
#include "gio/glog.h"
#include "gio/giof.h"
#include "gset/gsetign.h"
#include "gset/gcfg_multign.h"
#include "gproc/gpvtflt.h"
#include "gdata/gposdata.h"
#include "gutils/ExportINS.h"
#include "gvis/gvis.h"
using namespace gnut;
using namespace gvis;

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
	class INS_LIBRARY_EXPORT t_gmultign :virtual public t_gintegration,virtual public t_gvis
	{
	public:
		/**
		* @brief Constructor
		* @note set parameter value by gset
		*/
		t_gmultign(string site, string site_base, t_gsetbase* gset, t_gallproc *allproc);

		/**
		* @brief set whether t_gallpar is used or not
		*/
		void setgallparuse(bool status);

		/**
		* @brief struct store image data
		*/
		virtual void Add_IMG(t_gimgdata* img);
	    
		/**
		* @brief determines whether the interval between two timestamps is less than the delay
		*/
		virtual bool _img_time_valid(t_gtime imgt, t_gtime inst);

		/** @brief processBatch. */
		virtual int processBatch(const t_gtime& beg, const t_gtime& end, bool prtOut)override;
		virtual int processBatch(const t_gtime& beg, const t_gtime& end);

		

		/** @brief process one epoch.	*/
		virtual int _processGNSSEpoch(const t_gtime& runEpoch);

		/** @brief obtain all observation types 	*/
		virtual int _getMeas();

		/** @brief GNSS_update for integrated stystem 	*/
		virtual int _GNSS_Update();

		/** @brief Cam_update for integrated stystem 	*/
		virtual int _Cam_Update();
		/**
		* @brief process cam observation
		*
		* @param[in] runEpoch			cam time
		*/
		int _processCamEpoch(const t_gtime& runEpoch);

		bool _extrinsic_init();

		/**
		* @brief process cam observation
		*
		* @param[in] dx		correction of params	
		*/
		virtual int _cam_feedback(const ColumnVector & dx);

		/**
		* @brief process cam observation
		*
		* @param[in] dx		correction of params
		*/
		virtual int _cam_feedback(const Eigen::VectorXd & dx);

		/**
		* @brief StateIDType to string
		*
		* @param[in] id			cam state id
		*/
		string state_id2str(const StateIDType& id);
	


		/**
		* @brief perform cam measurement update
		*
		* @param[in] H			measurement matrix vector
		* @param[in] r			residual vector
		*/
		void cam_meas_update(const Eigen::MatrixXd& H, const Eigen::VectorXd& r);

		/**
		* @brief check whether is static or not
		*
		* @return
			retval = true.   static
			retval = false.  not static
		*/
		bool checkStaticMotion();

		/**
		* @brief state augment for a new cam state and complete time update of cov 
		*  between cam and imu
		*
		* @param[in] runEpoch		cam time
		*/
		void stateAugmentation(const t_gtime& runEpoch);


		/**
		* @brief perform chi-square test
		*
		* @param[in] H		measurement matrix vector
		* @param[in] r		residual vector
		* @param[in] dof	number of observation
		*/
		bool gatingTest(const Eigen::MatrixXd& H, const Eigen::VectorXd& r, const int& dof);

		/**
		* @brief process lost feature update for keyframe
		*/
		void removeLostFeatures();

		/**
		* @brief process sliding window update for keyframe
		*/
		void pruneCamStateBuffer();


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
		t_gtime _img_crt;					/// record cur image time
		t_gimgdata* _imgdata;

		/* Make the structure clear and reduce the temporary definition of intermediate variables in the structure  */
		PointCloud _pointcloud;             /// point cloud result,from font-end
		Eigen::Vector3d initCamPos = Eigen::Vector3d::Zero(); /// initial camera pos used in display
		bool mIsFirstImg = true;			/// A indicator to daetermine if it is the first image 
									        /// Maintain numerical and corresponding derivative stability
		double _img_t;						/// get cur image time
		imgPath _img_path;					/// get cur image path
		bool is_GNSS = false;
		bool is_VISION = false;
	private:
		bool is_allTC = false;			// whether t_gallpar is used or not
	};

}




#endif

#endif // USE_VIS