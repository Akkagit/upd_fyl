#ifndef GEARTH_H
#define GEARTH_H

/**
* @file gearth.h
* @details 
* @note t_gearth Class for Earth parameter updating.

* @verbatim
		History
		-0.1	ShenZhiheng	       2018-10-15 creat the file.
		-0.3	ShenZhiheng	       2019-03-20 add some base functions
		-0.5	ShenZhiheng	       2019-06-10 rename class name
		-0.7	ShenZhiheng	       2019-12-02 restruct the file
  @endverbatim

* @author ShenZhiheng
* @version 0.7
* @date 2018-10-15
* @license
*/



#include "gutility.h"
#include "gutils/ExportINS.h"

using namespace std;

namespace gsins
{

	/**
	* @class t_gearth
	* @brief t_gearth Class for Earth parameter updating.
	*
	* The Earth Parameter,such as e,a,b et al.
	* The main function is Earth Updating.
	*/
	class INS_LIBRARY_EXPORT t_gearth
	{
	public:

		/**
		* @brief constructor.
		*
		* @param[in]  a0   major semi axis  
		* @param[in]  f0   flattening
		* @param[in]  g0   gravity
		*/
		t_gearth(double a0 = glv.Re, double f0 = glv.f, double g0 = glv.g0);


		/**
		* @brief The Earth Parameter Update
		* The Earth Parameter mainly includes sb, sb2, sb4, cb, tb, sl, cl.
		*
		* @param[in]  pos        position of body
		* @param[in]  vn         velocity of body
		*/
		void Update(const Eigen::Vector3d &pos, const Eigen::Vector3d &vn);


		/**
		* @brief Compute dpos according to vn and ts
		*
		* @param[in]  vn         velocity of body
		* @param[in]  ts         interval
		* @return 
			Position Increments
		*/
		Eigen::Vector3d v2dp(const Eigen::Vector3d& vn, double ts);

	public:
		double a, b;							             /// major & Minor Axissemi axis 
		double f, e, e2, wie;					             /// flattening,eccentricity,..,angular of rotation
		double sb, sb2, sb4, cb, tb, sl, cl;				 /// sinB,..,..,cosB,tanB,sinL,cosL
		double RMh, RNh, cbRNh, f_RMh, f_RNh, f_cbRNh;		 /// RM+H,RN+H,RMH*cosB
		Eigen::Vector3d pos, vn;							 /// position,velocity
		Eigen::Vector3d wnie, wnen, wnin, weie, gn, gcc;
		Eigen::Matrix3d Cen, Cne;							 /// DCM from n-fram to e-frame

	};


}



#endif