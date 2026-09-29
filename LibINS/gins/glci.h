#ifndef GLCI_H
#define GLCI_H

/**
* @file glci.h
* @brief loosely coupled integration key algorithm
* @note
		lci is loosely-coupled integration;

* @verbatim
		History
		-0.5	ShenZhiheng	       2018-10-15 creat the file.
		-0.7	ShenZhiheng	       2019-12-26 restruct the file
  @endverbatim

* @author ShenZhiheng
* @date 2018-10-15
* @license
*/

#include "gins/gins.h"
#include "gutils/gmutex.h"
#include "gset/gsetbase.h"
#include "gset/gsetout.h"
#include "gdata/gimudata.h"
#include "gio/glog.h"
#include "gio/giof.h"
#include "gset/gsetign.h"
#include "gset/gcfg_ign.h"
#include "gdata/gposdata.h"
#include "gall/gallpcv.h"
#include "gall/gallproc.h"
#include "gall/gallotl.h"
#include "gall/gallbias.h"
#include "gutils/ExportINS.h"
using namespace gnut;

namespace gsins
{

	/**
	* @class t_glci
	*
	* t_glci Class for loosely-coupled integration
	*
	* loosely-coupled integration based on Adaptive and Sequential kalman filter
	* Inherited from t_gsinskf
	*
	*/
	class INS_LIBRARY_EXPORT t_glci : virtual public t_gsinskf
	{
	public:

		/** @brief default constructor. */
		t_glci();


		/**
		* @brief Constructor
		* @note set parameter value by dimension
		*/
		t_glci(int m0, int n0);


		/**
		* @brief Constructor
		* @note set parameter value by gset
		*/
		t_glci(t_gsetbase* gset);


		/**
		* @brief Initial Function
		*
		* virual function,inherited from t_gsinskf
		*/
		virtual void init(void) {}


		/** @brief init filter. */
		void init15();


		/**
		* @brief set measurement matrix Hk
		*
		* virual function,and inherited from t_gsinskf
		*
		*/
		virtual void set_Hk();


		/**
		* @brief set measurement Zk
		*
		* virual function,and inherited from t_gsinskf
		*
		*/
		virtual void set_Meas();


		/**
		* @brief set measurement Rk
		*
		* @param[in] Qvn			velocity variance
		* @param[in] Qpos			position variance
		*
		*/
		virtual void set_Rk(const Eigen::Vector3d& Qvn, const Eigen::Vector3d& Qpos);

		/**
		* @brief set external measurement
		*
		* @param[in] vn				velocity
		* @param[in] pos			position
		* @param[in] t				time
		*
		*/
		void set_Meas(const Eigen::Vector3d& vn, const Eigen::Vector3d& pos, double t);

		/**
		* @brief detect outlier by post-residual
		*
		* @param[in] _v_post			post-residual
		* @param[in] threshold			threshold
		* @param[out] index				index of outlier
		* @return
			outlier's index
		*/
		virtual int outlier_detect(double threshold);


		/**
		* @brief Kalman Update
		*
		* @param[in] wm		angle increment
		* @param[in] vm		velocity increment
		* @param[in] scm	scheme of processing
		* @param[out] fkf	output file
		*
		*/
		void Update(const vector<Eigen::Vector3d>& wm, const vector<Eigen::Vector3d>& vm, const t_scheme& scm);	 // kalman update

		/** @brief Kalman Feedback. */
		virtual void feedback(double kfts);
		virtual int _prt_ins_kml() { return 0; }

	};

}




#endif