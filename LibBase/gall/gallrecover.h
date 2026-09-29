/**
*
* @verbatim
History
-1.0 ZhengHJ  2019-09-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		  gallrecover.h
* @brief	  Storage all recover info
* @author     ZhengHJ, Wuhan University
* @version	  1.0.0
* @date		  2019-09-25
*
*/

#ifndef GALLRECOVER_H
#define GALLRECOVER_H

#include "gutils/ExportBase.h"

#include "gdata/grecoverdata.h"
#include "gall/gallrb.h"
#include "gall/gallsta.h"
#include "gall/galleop.h"
#include "gdata/gpoleut1.h"
#include "gall/gallpar.h"
#include <gall/gallprec.h>

namespace great
{

	typedef map<t_gtime, vector<t_grecover_equation*> > t_map_time_equ;
	typedef map<string, t_map_time_equ>	t_map_sat_equ;
	typedef map<string, t_map_sat_equ> t_map_site_equ;

	typedef map<t_gtime, vector<t_grecover_par*> > t_map_time_par;
	typedef map<t_gpar::t_type, vector<t_grecover_par*> > t_map_type_par;

	BASE_LIBRARY_EXPORT typedef tuple<string, string, string, double, double, t_gtime, t_gtime> t_tuple_ion; // type, rec, sat, value, sigma,beg,end

	class BASE_LIBRARY_EXPORT t_gallrecover: public t_gdata
	{
	public:
		t_gallrecover();
		~t_gallrecover();

		void add_recover_equation(const t_grecover_equation& recover_equ);
		void add_recover_par(const t_grecover_par& recover_par);


		//// get recover equation
		//vector<t_grecover_equation> get_recover_equation(string site);
		//vector<t_grecover_equation> get_recover_equation(string site, string sat);
		void get_clkdata(t_gallprec& clkdata,t_gallprec::clk_type type = t_gallprec::UNDEF);
		void get_clk13data(t_gallprec& clkdata, t_gallprec::clk_type type = t_gallprec::UNDEF);//add by xiongyun
		void get_iondata(vector< t_tuple_ion>& ion_data);
		void get_rbdata(t_gallrb& rbdata);
		void get_stadata(t_gallsta& stadata);
		void get_eopdata(t_galleop& eopdata);
		void get_poledata(t_gpoleut1& poledata);
		//void get_crddata(t_gallpar& pars);

		// jdhuang : for edit res of L1-L5
		vector<t_grecover_equation> get_first_phase_recover_equation(string site, string sat, string freq = "LC");


		//// get recover par
		//vector<t_grecover_par> get_recover_par(t_gtime time);
		vector<t_grecover_par> get_recover_par(t_gpar::t_type parType);

		t_gtime get_beg_time() const;
		t_gtime get_end_time() const;
		t_gtime get_equ_beg_time() const;
		double get_interval() const;
		void set_interval(double intv);

		// get all record
		const vector<t_grecover_data*>& get_all_recover_data() const;
		double get_recover_data_value(t_gpar::t_type parType);
		// get map_par
		const t_map_time_par& get_map_time_par() const;
		const t_map_time_equ& get_map_time_equ() const;
		const t_map_site_equ& get_map_site_equ() const;

		set<string> get_sat_list() const;
		set<string> get_site_list() const;
		set<t_gtime> get_time_list() const;
		vector<t_gtime> get_all_obs_time_list() const;
		
		// ½«recover par×ª»»Îª gallpar
		t_gallpar get_all_pars();

		// get nsat per epoch
		map<t_gtime, int> get_sat_number_per_epo() const;

	protected:
		t_grcover_head _recover_head;
		vector<t_grecover_data*> _recover_data;

	private:

		void _add_common_data(t_grecover_data*  data);

		// Mainly for index for time for and sat but no for storaging 
		t_map_time_equ _time_equmap;
		t_map_site_equ _site_sat_time_equmap;
		t_map_time_par _time_parmap;
		map<t_gtime, vector<t_grecover_par> > _time_parmap_pce;

		t_map_type_par _type_parmap;
	};


}


#endif // !GALLRECOVER_H