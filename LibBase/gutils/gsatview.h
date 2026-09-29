
#ifndef GSATVIEW_H
#define GSATVIEW_H 

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: satellite sky view
  Version: $ Rev: $

  2015-01-03 /JD: created

-*/

#include <map>
#include <string>
#include <iomanip>
#include <iostream>

#include "gutils/gnss.h"
#include "gutils/gconst.h"
#include "gutils/ExportBase.h"
#include "gall/gallnav.h"

using namespace std;

namespace gnut {

	typedef map<t_gtime, t_gtime>              t_map_sky_epo; // sat visibility of a single sat/site
	typedef map<string, t_map_sky_epo>        t_map_sky_sat; // sat visibility of a single satellite

	int sat_view(t_map_sky_sat& sat_view,
		t_gtriple      xyz_rec,
		t_gallnav* all_nav,
		const t_gtime& beg,                        // should be setup reasonably !
		const t_gtime& end,                        // should be setup reasonably !
		double mask = 0.0);

	int sat_elev(t_map_sky_epo& sat_elev,
		t_gtriple      xyz_rec,
		t_gallnav* all_nav,
		const t_gtime& beg,                        // should be setup reasonably !
		const t_gtime& end,                        // should be setup reasonably !
		const  string& prn,
		double mask);

} // namespace

#endif // # GCOMMON_H
