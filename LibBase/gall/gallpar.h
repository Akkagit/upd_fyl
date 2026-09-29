
#ifndef GALLPAR_H
#define GALLPAR_H 

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: parametres container
  Version: $ Rev: $

  2011-04-18 /PV: created

-*/

#include <string>
#include <set>

#include "gdata/gsatdata.h"
#include "gutils/gtriple.h"
#include "gmodels/gpar.h"
#include <unordered_map>
#include <algorithm>

using namespace std;

namespace gnut
{

	class BASE_LIBRARY_EXPORT  t_gallpar
	{
	public:
		void addParam(t_gpar);
		void delParam(int i);
		void delAllParam();
		void resetAllParam();

		set<int> delEmpty();
		vector<int> delAmb();

		int getParam(string site, t_gpar::t_type, string prn,
			t_gtime beg = FIRST_TIME, t_gtime end = LAST_TIME) const;
		int getParam(int index);
		int getParIndex(int idx);
		double getParValue(int idx);
		void setParValue(int idx, double value);
		void setParRemove(int idx, bool judge);

		unsigned int parNumber() const;
		unsigned int ambNumber() const;
		unsigned int orbParNumber() const;

		int getAmbParam(string site, string prn, t_gpar::t_type  type, 
			const t_gtime& beg = FIRST_TIME, const t_gtime& end = LAST_TIME) const;
		t_gpar& getAmbParam(int idx);

		int getCrdParam(string site, t_gtriple& crd,
			t_gtime beg = FIRST_TIME, t_gtime end = LAST_TIME) const;
		int getVelParam(string site, t_gtriple& crd,
			t_gtime beg = FIRST_TIME, t_gtime end = LAST_TIME) const;

		vector<int> getPartialIndex(string site, string sat);

		ColumnVector get_cvect(t_gallpar& par);
		t_gpar&   getPar(int idx);
		t_gpar&   operator[](const size_t idx);
		t_gallpar operator-(t_gallpar& par);
		t_gallpar operator+(t_gallpar& par);
		static int mult(const Matrix& K, t_gallpar& par, t_gallpar& res);
		void   reIndex();
		void   decIndex(int i);
		void   incIndex(int i);
		void   setSite(string site);
		set<string> amb_prns();
		set<string> amb_recs();
		map<string, string> rec_sys_map();
		int sum(t_gallpar& X1, t_gallpar& X2);

		//add for LEO
		void setOrbPar(string sat, vector<t_gpar> pars);
		vector<t_gpar>  getOrbPar(string sat);
		vector<t_gpar>  getAllPar();
		map<string, vector<t_gpar>> getAllOrbPar();
		friend ostream& operator<<(ostream& os, t_gallpar& x);
		// µ¥Æµ/Ë«Æµ/ÈýÆµ/... ÎÀÐÇÊý
		map<string, int> freq_sats_num(int freq);

	private:
		vector<t_gpar> _vParam;
		map<string, vector<t_gpar>> _vOrbParam;  // add for orb param

		// Fast get param index
		map<t_gparhead, map<t_gtimearc, long > > _index_par;  // index -> point
		long _max_point=0;
		vector<long> _point_par;  // point -> vParam

		// Fast get parital  
		map<pair<string, string>, vector<int> > _index_for_parital; // index->param
		pair<long, int> _last_point;
		void _update_partial_index();

		t_gmutex _allpar_mtx;
	};

} // namespace

#endif
