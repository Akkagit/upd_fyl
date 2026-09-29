/**
*
* @verbatim
History
       -1.0 Qian Zhang   2019-11-15 Creat the file.
	   -1.2 Jiande Huang  2020-09-30 Standard code style and Optimize code logic
* Copyright (c) 2018, Wuhan University.All rights reserved.
*
* @file		    gdtm2013.h
* @brief		    The class for dtm2013 model, serving for atmospheric drag model by providing the density.
*
* @author          Qian Zhang, Wuhan University
* @version		    1.0.0
* @date		    2019-11-15
*
*/
#ifndef GDTM2013_H_
#define GDTM2013_H_

#include "gforce/gfatmosphere.h"
#include "gutils/gconst_oi.h"
#include "gdata/gnavde.h"



using namespace std;

namespace great
{
	/**
	* @brief dtm2013 model, provide the density for the calculation of acc and cmat of atmosphere.
	* @details an old model which serving for atmospheric drag model, need updating.
	*/
	class OI_LIBRARY_EXPORT t_gfatmos_dtm2013 : public t_gforce
	{
	public:
		/** @brief constructor. */
		t_gfatmos_dtm2013(t_glog * log);

		/** @brief destructor. */
		virtual ~t_gfatmos_dtm2013();

		/**
		* @brief coordinate transforming, get the BLH according to the position in trs.
		* @param[in] crd the position in trs
		* @param[out] ell the BLH
		* @return
		@retval >= 0 calculate normally
		@retval < 0  have some error
		*/
		bool xyz2blh(ColumnVector &crd, ColumnVector &ell);
		
		/**
		* @brief calculate the dendity.
		* @param[in] sat_pos the sat pos in trs
		* @param[in] sun_pos the sun pos in trs
		* @param[in] epoch the dynamic time
		* @param[in] kp the kp data interpolated from geomag_kp file
		* @param[in] solar the solar data interpolated from solar_flux file
		* @return    density
		*/
		double calDensity(float day, float f, float fbar, double * akp, float alt, float stl, float lat, float lon);  // int mon, int year, int hr, int min, double sec,
	};
}
#endif // ! GDTM2013_H
