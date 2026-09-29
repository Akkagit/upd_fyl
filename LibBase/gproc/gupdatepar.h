/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2019-12-23  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gupdatepar.h
* @brief		define update parameter base model
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-12-23
*
*/

#ifndef GUPDATEPAR_H
#define GUPDATEPAR_H

#include "gall/gallpar.h"
#include "gutils/gtime.h"
#include "gset/gsetbase.h"
#include "gutils/gcycleslip.h"
#include "gproc/glsqmatrix.h"
#include "gmodels/gstochasticmodel.h"

using namespace gnut;

namespace great
{
	class BASE_LIBRARY_EXPORT t_gupdateparinfo
	{
	public:

		t_gupdateparinfo();
		~t_gupdateparinfo();

		bool exist(int id);

		// add remove info 
		void add(int id);

		void add(t_gpar newpar);

		// add update info
		void add(int id,t_gpar newpar);

		// add update info with state equ
		void add(vector<pair<int,t_gpar> > update_par, t_glsqEquationMatrix state_equ);

		// get remove info
		void get(vector<int>& remove_id);
		// get new parlist info
		void get(vector<t_gpar>& newparlist);
		// get state equ info
		void get(vector<t_gpar>& update_newpar,t_glsqEquationMatrix& equ);//change by xiongyun
		// remove all info
		void get(vector<int>& remove_id,vector<t_gpar>& newparlist,vector<t_gpar>& equ_parlist, t_glsqEquationMatrix& equ);
	


	private:
		set<int> _remove_id;
		vector<t_gpar> _new_parlist;
		vector<t_gpar> _equ_parlist;
		t_glsqEquationMatrix _state_equ;
	};



	class BASE_LIBRARY_EXPORT t_gupdatepar
	{
	public:
		t_gupdatepar();
		virtual ~t_gupdatepar();

		t_gupdatepar(const t_gupdatepar& Other);

		void set_interval(double interval);
		void set_sig_ion(double sigion);
		void set_cycleslip(shared_ptr<t_gcycleslip> cycleslip);
		void set_amb_update_way(bool way);
		
		void set_par_state_mode(t_gpar::t_type type, int order, double dt, double noise);

		void update_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info);

		/**
		* @brief  all lsq par with now obs data
		* @note according to obsdata update amb par , and time update outsate par
		* @param[in] update epoch
		* @param[in] all lsq updateparameters
		* @param[in] obsdata observ info data
		* @return all update idx and newpar( if newpar is none represent only remove but not add) 
		*/
		virtual t_gupdateparinfo get_all_update_parameters(const t_gtime& epoch,t_gallpar& allpars,const vector<t_gsatdata>& obsdata);

		

	protected:

		virtual void _update_amb_pars(const t_gtime& epoch,t_gallpar& allpars,const vector<t_gsatdata>& obsdata,t_gupdateparinfo& update_info) = 0;

		// must last use
		void _update_process_pars(const t_gtime& epoch, t_gallpar& allpars,t_gupdateparinfo& update_info);

		// update one process par
		void _update_state_par(const t_gtime& epoch,int update_id, t_gallpar& allpars, t_gupdateparinfo& update_info);
		void _update_process_par(const t_gtime& epoch, int update_id, t_gallpar& allpars, t_gupdateparinfo& update_info);
		void _udpate_gps_rec_ifb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info);

	protected:
		
		double _intv = 0.0;
		double _sig_ion = 9000;

		shared_ptr<t_gcycleslip> _cycleslip;
		map<t_gpar::t_type, t_statemode> _state_mode;
		//t_gsetbase* _set;//add by xiongyun
		bool _lite_update_amb;
	};



}


#endif /* GUPDATEPAR_H */