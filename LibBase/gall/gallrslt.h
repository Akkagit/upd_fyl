
#ifndef GALLRSLT_H
#define GALLRSLT_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
	  Research Institute of Geodesy, Topography and Cartography
	  Ondrejov 244, 251 65, Czech Republic

  Purpose: container for all biases
  Version: $ Rev: $

  2017-11-05 /PV: created

-*/

#include <vector>
#include <set>
#include <iostream>

#include "gutils/gtime.h"

namespace gnut
{

	class BASE_LIBRARY_EXPORT t_gallrslt
	{
	public:
		t_gallrslt();

		struct result
		{
			string type;
			string prn;

			t_gtime beg;
			t_gtime end;

			int index;

			double adj;
			double rms;
			double val;

			bool  operator <(const t_gallrslt::result&) const;
		};

		vector<t_gallrslt::result> v_rslt;

		set<t_gallrslt::result> set_crd;
		set<t_gallrslt::result> set_trp;
		set<t_gallrslt::result> set_clk;
		set<t_gallrslt::result> set_amb;

		void append(const t_gallrslt::result &rslt);
		void append(string type, t_gtime beg, t_gtime end, int idx,
			string prn, double adj, double rms, double val);

	};

} // namespace

#endif // GALLRSLT_H
