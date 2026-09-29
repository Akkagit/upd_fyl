/**
*
* @verbatim
History
	-1.0 Guolong Feng  2020-07-13  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @version	1.0.0
* @date		2020-07-13
*
*/

#ifndef GIONEX_H
#define GIONEX_H

#include "gutils/ExportBase.h"
#include "gutils/gsysconv.h"
#include "gmodels/ginterp.h"

using namespace gnut;

namespace great
{
	class BASE_LIBRARY_EXPORT t_gionex_head
	{
	public:
		t_gionex_head();
		virtual ~t_gionex_head();
	
		/** @brief assignment operator . */
		void operator=(const t_gionex_head& Other);

		t_gtime                             beg;              // EPOCH OF FIRST MAP
		t_gtime                             end;              // EPOCH OF LAST MAP
		int                                 interval;         // interval between the TEC maps
		int                                 map_dimension;    // Dimension of TEC/RMS maps -- 2/3
		tuple<double, double, double>       hgts;             // definition of an equidistant grid in height [HGT1 / HGT2 / DHGT]
			                                                  // HGT1 -> HGT2 with increment DHGT (in km), e.g. 50.0 450.0  50.0
			                                                  // for 2-dimension, HGT1 = HGT2. here, stored in meters
		tuple<double, double, double>       lons;             // definition of an equidistant grid in longitude [LON1 / LON2 / DLON]
													          // LON1 -> LON2 with increment DLON (in degree), e.g. -180.0 180.0   5.0
		tuple<double, double, double>       lats;             // definition of an equidistant grid in latitude [LAT1 / LAT2 / DLAT]
															  // LAT1 -> LAT2 with increment DLAT (in degree), e.g. 87.5 -87.5  -2.5
		int                                 nhgt;
		int                                 nlon;
		int                                 nlat;
		int                                 exponent;         // default -1
		double                              base_radius;      // mean earth radius or bottom of height, here, stored in meters
		IONMPFUNC                           ion_mapfunc;      // COSZ: 1/cos(z) / QFAC: Q-factor
		map<string, pair<double, double> >  p1p2_dcb;         // sat  bias[ns]  rms
	};


	class BASE_LIBRARY_EXPORT t_gionex_data
	{
	public:
		t_gionex_data() {};
		virtual ~t_gionex_data() {};

		map< tuple<int, int, int>, double >  tec_val;   // lat_index/lon_index/height_index , unit TECU
		map< tuple<int, int, int>, double >  tec_rms;   // lat_index/lon_index/height_index , unit TECU

	};

	/**
	*@brief	 Class for storaging ionex file data
	*/
	class BASE_LIBRARY_EXPORT t_gionex : public t_gdata
	{
	public:

		/** @brief default constructor. */
		t_gionex();

		/** @brief default destructor. */
		virtual ~t_gionex();

		/**
		* @brief add head of GIM grid files.
		* @param[in] ionex_hd   head structure.
		*/
		void add_head(t_gionex_head ionex_hd);

		/**
		* @brief add TEC grid data.
		* @param[in] time             curruent epoch.
		* @param[in] type             TEC_MAP/RMS_MAP
		* @param[in] ilat/ilon/ihgt   position of t_gionex_data structure.
		* @param[in] v                corresponding TEC.
		* @return	                  void
		*/
		void add_data(const t_gtime& time, string type, int ilat, int ilon, int ihgt, double v);
		
		/**
		* @brief get STEC by GIM grid data (position and time interpolation)
		* @param[in] time        curruent epoch.
		* @param[in] satdata     satellite data, provide elevation and azimuth.
		* @param[in] site_ell    receiver position {lat,lon,h}.
		* @param[out] value      STEC value for the satellite, unit : TECU.
		* @param[out] rms        rms value, unit : TECU.
		* @return	             true success, false fail
		*/
		bool getSTEC(const t_gtime& time, t_gsatdata& satdata, t_gtriple& site_ell, double& value, double& rms);

	protected:
		t_gionex_head                  _ionex_hd;
		map<t_gtime, t_gionex_data>    _map_ionex_data;

		/**
		* @brief interpolate tec grid data according ipp coordinate (position interpolation)
		* @param[in] time        curruent epoch.
		* @param[in] ipp_ell     ionospheric pierce point position(B/L).
		* @param[in] ihgt        which ionosphere layer
		* @param[out] value      VTEC at ipp, unit : TECU.
		* @param[out] rms        rms of VTEC at ipp, unit : TECU.
		* @return	             true success, false fail
		*/
		bool _interpolate_ion(const t_gtime& time, const t_gtriple& ipp_ell, const int& ihgt, double& value, double& rms);
		
	};
}

#endif // !GION_H
