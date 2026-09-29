/**
 *
 * @verbatim
     History
      -1.0 Yujie Qin and Bo Wong   2019-01-06 Creat the file.
	  -1.1 Hongmin Zhang           2019-03-04 update.
	  -1.2 Yujie Qin               2019-03-05 edit for final version.
      -1.3 Yujie Qin               2019-03-15 Code comments are added in Doxygen format.
	  -1.4 Yujie Qin               2019-07-03 Adding Doxygen Style Code Remarks
	  -1.5 Jiande Huang            2020-09-30 Standard code style and Optimize code logic
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		    gdtm94.h
 * @brief		    The class for dtm94 model, serving for atmospheric drag model by providing the density.
*
 * @author          Yujie Qin, Wuhan University
 * @version		    1.0.0
 * @date		    2019-01-06
 *
 */
#ifndef GDTM94_H
#define GDTM94_H

#include "gforce/gfatmosphere.h"
#include "gutils/gconst_oi.h"
#include "gdata/gnavde.h"

using namespace std;

namespace great
{
	/**
	* @brief dtm94 model, provide the density for the calculation of acc and cmat of atmosphere.
	* @details an old model which serving for atmospheric drag model, need updating.
	*/
	class OI_LIBRARY_EXPORT t_gfatmos_dtm94 : public t_gforce
    {
	 public:
		 /** @brief constructor. */
	    t_gfatmos_dtm94(t_glog * log);

		/** @brief destructor. */
	    virtual ~t_gfatmos_dtm94();

		/** 
		* @brief calculate the dendity.
	    * @param[in] sat_pos the sat pos in trs
	    * @param[in] sun_pos the sun pos in trs
	    * @param[in] epoch the dynamic time
	    * @param[in] kp the kp data interpolated from geomag_kp file
	    * @param[in] solar the solar data interpolated from solar_flux file
	    * @return    density
	    */
	    double cal_density(ColumnVector &sat_pos, ColumnVector &sun_pos, double &epoch, double &kp, vector<double> &solar);

		/**
		* @brief check the force class
		*/
		bool check();

	  protected:
		double         _kp;        ///< the kp data interpolated from geomag_kp file according the _epoch.
		double         _epoch;     ///< the current gps time.
		double         _pnm[12];   ///< the non-normalized Legendre associated functions
		double         _GL;        ///< the expansion function G(L) for model quantities
		vector<double> _solar;     ///< the solar data interpolated from solar_flux file according the _epoch.
	    ColumnVector   _sat_geo_crd;    ///< geodetic coordinate of satellite(BLH)
	    ColumnVector   _sun_geo_crd;    ///< geodetic coordinate of sun(BLH)

		/**
		* @brief calculate Legendre polynomial.
		* @return
			@retval >= 0 calculate normally
			@retval < 0  have some error
		*/
	    bool cal_Pnm();

		/** 
		* @brief calculate the Gl, Gl:Empirical Function.
	    * @param[in] gas the gas to calculate
	    * @param[in] temp use index counting or not
	    * @param[in] horhe use solar_activity or not
		* @return
			@retval >= 0 calculate normally
			@retval < 0  have some error
	    */
	    bool cal_GL(const double* gas, bool temp = true, bool horhe = false ); 

		/**
		* @brief coordinate transforming, get the BLH according to the position in trs.
		* @param[in] crd the position in trs
		* @param[out] ell the BLH
		* @return
		@retval >= 0 calculate normally
		@retval < 0  have some error
		*/
		bool xyz2blh(ColumnVector &crd, ColumnVector &ell);
  };
} // namespace

#endif // !GDTM94_H