
/**
*
* @verbatim
History
	-1.0 Guolong Feng  2020-07-17  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		giono.h
* @brief	implements ionosphere model class

* @version	1.0.0
* @date		2020-07-17
*
*/


#ifndef GIONO_H
#define GIONO_H 

#include "gall/gallnav.h"
#include "gdata/gionex.h"
#include "gdata/gsatdata.h"
#include "gutils/ExportBase.h"

using namespace std;
using namespace gnut;

namespace great
{

	class BASE_LIBRARY_EXPORT t_giono
	{
	public:
		t_giono();
		virtual ~t_giono();

		// mapping function
		double ionMapFunc(IONMPFUNC mapfunc, double ele);

	protected:
	};

	class BASE_LIBRARY_EXPORT t_giono_brdc : public t_giono {
	public:
		t_giono_brdc() {}
		~t_giono_brdc() {}

		/**
		* @brief compute ionospheric delay by broadcast ionosphere model (klobuchar model)
		* @param[in] nav         store iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3}
		* @param[in] epo         curruent epoch.
		* @param[in] satdata     provide azimuth/elevation angle {az,el} (rad)
		* @param[in] site_ell    receiver position {lat,lon,h} (rad,m)
		
		* @return	             ionospheric delay (L1) (m)
        **/
		double getIonoDelay(t_gallnav* nav, t_gsatdata& satdata, const t_gtime& epo, const t_gtriple& site_ell);

	};


	class BASE_LIBRARY_EXPORT t_giono_sbas : public t_giono {
	public:
		t_giono_sbas() {}
		~t_giono_sbas() {}

	};
	
	
	
	// get ionospheric delay by tec grid data 
	class BASE_LIBRARY_EXPORT t_giono_tecgrid : public t_giono {
	public:
		t_giono_tecgrid() {}
		~t_giono_tecgrid() {}

		bool getIonoDelay(t_gionex* ionexdata, t_gsatdata& satdata, const t_gtime& epo, t_gtriple& site_pos, double& value, double& rms);

	};
	
	
	

} // namespace

#endif
