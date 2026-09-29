/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2019-12-23  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gupdatepar.h
* @brief		define update parameter IF model
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-12-23
*
*/

#ifndef GUPDATEPARIF_H
#define GUPDATEPARIF_H

#include <map>
#include "gproc/gupdatepar.h"

namespace great {
	class BASE_LIBRARY_EXPORT t_gupdateparIF :public t_gupdatepar
	{
	public:
		t_gupdateparIF(shared_ptr<t_gcycleslip> cycleslip, const map< GSYS, map<FREQ_SEQ, GOBSBAND> >& band_list);
		//t_gupdateparIF(shared_ptr<t_gcycleslip> cycleslip, t_gsetbase* set);//add by xiongyun
		virtual ~t_gupdateparIF();

	protected:

		virtual void _update_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata,t_gupdateparinfo& update_info) override;

		void _update_one_type_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata,t_gpar::t_type amb_type,t_gcycleslip* _cycleslip,t_gupdateparinfo& update_info);
		// add by xiongyun for realtime
		void _update_one_type_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gpar::t_type amb_type, t_gupdateparinfo& update_info);
	
	protected:
		map<pair<string, string>, int> _amb_id;
		map<GSYS, map<FREQ_SEQ, GOBSBAND> >  _band_list;

	};

	class BASE_LIBRARY_EXPORT t_gupdateparIF13 :public t_gupdateparIF
	{
	public:
		t_gupdateparIF13(shared_ptr<t_gcycleslip> cycleslip12,shared_ptr<t_gcycleslip> cycleslip13, const map< GSYS, map<FREQ_SEQ, GOBSBAND> >& band_list);
		virtual ~t_gupdateparIF13();

		void set_cycleslip13(shared_ptr<t_gcycleslip> cycleslip13);
		t_gupdateparinfo get_all_update_parameters(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata) override;

	protected:

		void _update_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info) override;


		shared_ptr<t_gcycleslip> _cycleslip13;
	};

	class BASE_LIBRARY_EXPORT t_updateparWL : public t_gupdateparIF
	{
	public:
		t_updateparWL(shared_ptr<t_gcycleslip> cycleslip, const map< GSYS, map<FREQ_SEQ, GOBSBAND> >& band_list);
		virtual ~t_updateparWL();

	protected:
		virtual void _update_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info) override;

	};

}
#endif // !GUPDATEPARIF_H
