#ifndef GIMU_H
#define GIMU_H

/**
* @file gimu.h
* @brief imu core algorithm
* @details
		IMU update algorithm,including cone,poly error compensation

* @verbatim
		History
		-0.1	ShenZhiheng	       2018-10-15 creat the file.
		-0.5	ShenZhiheng	       2019-03-20 add transition class t_gsinskf to inheriting for lci and tci.
		-0.7	ShenZhiheng	       2019-12-26 restruction.
  @endverbatim

* @author ShenZhiheng
* @email
* @version 1.0
* @date 2018-10-15
* @license
*/


#include "gbase.h"
#include "gquat.h"
#include "gearth.h"
#include "gavar.h"
#include "gkalman.h"

namespace gsins
{

	/**
	* @class t_gimu
	* @brief t_gimu Class is IMU core class for imu updating.
	*
	* t_gimu is used for simu core algorithms,including cone,poly error compensation
	* and the main function is imu info updating.
	*/
	class INS_LIBRARY_EXPORT t_gimu
	{
	public:

		/** @brief default constructor. */
		t_gimu();

		/**
		* @brief Phi cone compensation
		* @param[in] wm				angle increments
		* @return 
			result of cone error.
		*/
		Eigen::Vector3d phi_cone_crt(const vector<Eigen::Vector3d>& wm);


		/**
		* @brief Phi poly compensation
		* @param[in] wm				angle increments
		* @return 
			result of poly error.
		*/
		Eigen::Vector3d phi_poly_crt(const vector<Eigen::Vector3d>& wm);

		/**
		* @brief Vm cone compensation
		* @param[in] wm				angle increments
		* @param[in] vm				velocity increments
		* @return 
			result of cone error.
		*/
		Eigen::Vector3d vm_cone_crt(const vector<Eigen::Vector3d>& wm, const vector<Eigen::Vector3d>& vm);


		/**
		* @brief Vm poly compensation
		* @param[in] wm				angle increments
		* @param[in] vm				velocity increments
		* @return
			result of poly error.
		*/
		Eigen::Vector3d vm_poly_crt(const vector<Eigen::Vector3d>& wm, const vector<Eigen::Vector3d>& vm);


		/**
		* @brief IMU error compensation
		*
		* @param[in] wm				angle increments
		* @param[in] vm				velocity increments
		* @param[in] scm			scheme of processing
		*/
		void Update(const vector<Eigen::Vector3d>& wm, const vector<Eigen::Vector3d>& vm,const t_scheme& scm);

		
		Eigen::Vector3d phim, dvbm;		/// increment because of poly or cone.

	private:
		Eigen::Vector3d wm_1, vm_1;		/// last increment.
		double pf1, pf2;				/// flag of first time.

	};

}



#endif
