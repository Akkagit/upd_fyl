#ifndef GALLBIAS_H
#define GALLBIAS_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: container for all biases
  Version: $ Rev: $

  2012-11-05 /JD: created
  
  History GREAT:
   -1.0 Jiande Huang : Optimizing the efficiency of "get" function

  -*/

#include "gdata/gdata.h"
#include "gmodels/gbias.h"
#include "gutils/gtime.h"

using namespace std;

namespace gnut
{

	typedef shared_ptr<t_gbias> t_spt_bias;

	// ----------
	class BASE_LIBRARY_EXPORT t_gallbias : public t_gdata
	{

		typedef map<GOBS, t_spt_bias>   t_map_gobs;
		typedef map<string, t_map_gobs> t_map_sat;
		typedef map<t_gtime, t_map_sat> t_map_epo;
		typedef map<string, t_map_epo>  t_map_ac;

	public:
		t_gallbias();
		virtual ~t_gallbias();

		void   add(const string& ac, const t_gtime& epo, const string& obj, t_spt_bias pt_bias);                          // set single bias element value
		double get(const t_gtime& epo, const string& obj, const GOBS& gobs1, const GOBS& gobs2, string ac = "");          // get single bias element
		double get(const string prd, const t_gtime& epo, const string& obj, const GOBS& gobs1, const bool meter = true);  // get single bias element
		vector<string> get_ac();
		string ac_priority();
		void   set_used_ac(string ac);
		string get_used_ac();
		bool is_osb(); // if observation specific bias

		void overwrite(bool b) { _overwrite = b; } // set/get overwrite mode
		bool overwrite() { return _overwrite; }
		void clean_outer(const t_gtime& beg, const t_gtime& end);
	protected:

		t_spt_bias _find(const string& ac, const t_gtime& epo, const string& obj, const GOBS& gobs);   // get single bias element pointer
		vector<t_spt_bias> _find_ref(const string& ac, const t_gtime& epo, const string& obj, const GOBS& ref);    // get all biased pointers with particular reference signal

		void	_convert_obstype(const string& ac, const string& obj, GOBS& obstype);

		void       _connect_first(const t_spt_bias pt_cb1, t_spt_bias pt_cb2);
		void       _connect_second(const t_spt_bias pt_cb1, t_spt_bias pt_cb2);
		void       _consolidate(const string& ac, const string& obj, const t_spt_bias pt_cb1, t_spt_bias pt_cb2);

		t_map_ac   _mapbias;                  // map of all satellite biases (all ACs & all period & all objects)   
		string     _used_ac;
		bool       _overwrite = false;

		map<string, int> _ac_order;
		bool             _isOrdered = false;
		string           _pri_ac = "DLR_R";
	};

} // namespace

#endif
