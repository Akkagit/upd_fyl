/**
*
* @verbatim
	History
	 -1.0 Huang Jiande  2019-07-02  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gtideIERS.h
* @brief		precise model for computing  correction
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-18
*
*/

#ifndef GTIDEIERS_H
#define GTIDEIERS_H

#include "gmodels/gtide.h"
#include "gdata/gnavde.h"

namespace great
{

	/// The class for gtide of IERS model
	class BASE_LIBRARY_EXPORT t_gtideIERS : public t_gtide
	{
	public:
		/// constructor
		t_gtideIERS(t_glog* l);
		/// desconstructor
		virtual ~t_gtideIERS() {};

		/// Mark sure the pos in J2000,
		/// solid earth tides
		t_gtriple tide_solid(const t_gtime& epo, t_gtriple& xyz, Matrix& rot_trs2crs, t_gnavde* nav_planet);

		/// pole tides
		t_gtriple tide_pole();
		/// pole tides
		t_gtriple tide_pole_pod(const t_gtime& epo, double xpole, double ypole, t_gtriple& xyz);
		/// ocean tide loading
		t_gtriple load_ocean(const t_gtime& epoch, const string& site, const t_gtriple& xRec) override;
		/// atmospheric tide loading
		t_gtriple load_atmosph() override;
		t_gtriple tide_freq(const string& site, const t_gtriple & xRec, double gast);
		void getMeanPole(double mjd, double& xpm, double& ypm);

	protected:

	}; 
}// namespace

#endif // !GTIDEPOD_H
