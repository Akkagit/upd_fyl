#ifndef GTIDE_H
#define GTIDE_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: implements tides
  Version: $ Rev: $

  2012-11-05 /JD: created

-*/

#include <vector>

#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

#include "gdata/gdata.h"
#include "gutils/gconst.h"
#include "gutils/gtime.h"
#include "gutils/gmutex.h"
#include "gutils/gtriple.h"
#include "gmodels/gephplan.h"
#include "gall/gallotl.h"
#include "gutils/ExportBase.h"

using namespace std;

namespace gnut {

	// ----------
	class BASE_LIBRARY_EXPORT t_gtide
	{

	public:
		t_gtide(t_glog* l);
		virtual ~t_gtide();

		void   glog(t_glog* l) { _log = l; }

		virtual void setOTL(t_gallotl* gallotl);

		virtual t_gtriple tide_searth(const t_gtime& epo, t_gtriple& xyz);       // solid earth tides
		virtual t_gtriple tide_pole();                                           // pole tides
		virtual t_gtriple load_ocean(const t_gtime& epoch, const string& site, const t_gtriple& xRec);  // ocean tide loading
		virtual t_gtriple load_atmosph();                                        // atmospheric tide loading

	protected:
		t_gephplan   _planEph;
		t_gallotl*   _gotl;
		t_gmutex     _mutex;

		t_glog*      _log;

#ifdef BMUTEX
		boost::mutex _mutex;
#endif   
	};

} // namespace

#endif
