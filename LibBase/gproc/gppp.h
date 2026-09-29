
#ifndef GPPP_H
#define GPPP_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
	  Research Institute of Geodesy, Topography and Cartography
	  Ondrejov 244, 251 65, Czech Republic

  Purpose: implements PPP client abstract class
  Version: $ Rev: $

  2012-06-06 /JD: created
  2013-06-18 /JD: joined gpppfilter & glsq

-*/
#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

#include "gproc/gspp.h"
#include "gdata/gsatdata.h"
#include "gall/gallotl.h"
#include "gall/gallrslt.h"
#include "gmodels/gtropo.h"
#include "gmodels/gpar.h"
#include "gmodels/gpar.h"
#include "gmodels/gtropo.h"
#include "gmodels/gtide.h"
#include "gmodels/gtide96.h"
#include "gmodels/gtide2010.h"
#include "gmodels/gephplan.h"

#ifdef _WIN32
#pragma warning(disable:4250)                   // suppress Visual Studio WARNINGS about inheritance via dominance
#endif

namespace gnut
{

	class BASE_LIBRARY_EXPORT t_gppp : public virtual t_gspp
	{
	public:
		t_gppp(string mark, t_gsetbase* set);
		virtual ~t_gppp();

		virtual void setOTL(t_gallotl*   gotl);
		virtual void setOBJ(t_gallobj*   gobj);

		virtual bool isRunning();

	protected:

		/** @brief Get settings from XML file and set local variables. */
		virtual int   _get_settings();

		/** @brief Print results. */
		virtual void  _prt_results(t_giof* giof, const t_gallrslt& rslt);

		t_gallotl*    _gotl;             /**< Ocean loading date. */
		t_gtide*      _tides;            /**< Tides. */

		bool          _running;          /**< Running status. */
		double        _sigAmbig;         /**< Ambiguity sigma (constraints). */

#ifdef BMUTEX
		boost::mutex  _mutex;
#endif
		t_gmutex      _gmutex;

	private:
	};

} // namespace

#endif //GPPP_H
