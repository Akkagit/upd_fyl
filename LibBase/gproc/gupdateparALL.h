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

#ifndef GUPDATEPARALL_H
#define GUPDATEPARALL_H

#include <map>
#include "gproc/gupdatepar.h"

using namespace gnut;

namespace great
{


	class BASE_LIBRARY_EXPORT t_gupdateparALL :public t_gupdatepar
	{
	public:
		//t_gupdateparALL(shared_ptr<t_gcycleslip> cycleslip, t_gsetbase* _gset);
		t_gupdateparALL(shared_ptr<t_gcycleslip> cycleslip, const map< GSYS, map<FREQ_SEQ, GOBSBAND> >& band_list, const IONO_ORDER& iono_order);
		t_gupdateparALL(shared_ptr<t_gcycleslip> cycleslip_12, shared_ptr<t_gcycleslip> cycleslip_13, map< GSYS, map<FREQ_SEQ, GOBSBAND> > band_list, const IONO_ORDER& iono_order);
		t_gupdateparALL(shared_ptr<t_gcycleslip> cycleslip_12, shared_ptr<t_gcycleslip> cycleslip_13, shared_ptr<t_gcycleslip> cycleslip_14, shared_ptr<t_gcycleslip> cycleslip_15, map< GSYS, map<FREQ_SEQ, GOBSBAND> > band_list, const IONO_ORDER& iono_order);
		virtual ~t_gupdateparALL();

		t_gupdateparinfo get_all_update_parameters(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata) override;

	protected:
		//void _update_amb_pars_triple_freq(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info);
		//void _update_amb_pars_triple_freq_new(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info);
		//void _update_amb_pars_double_freq(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info);
		void _update_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info)override;
		void _lite_update_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info);
		void _udpate_ion_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info);

		void _check_log(shared_ptr<t_gcycleslip> cycleslip_12, shared_ptr<t_gcycleslip> cycleslip_13);
		// bands list;
		map< GSYS, map<FREQ_SEQ, GOBSBAND> >  _band_list;
		shared_ptr<t_gcycleslip>              _cycleslip_13 = nullptr;
		shared_ptr<t_gcycleslip>              _cycleslip_14 = nullptr;
		shared_ptr<t_gcycleslip>              _cycleslip_15 = nullptr;

		IONO_ORDER                            _iono_order = IONO_ORDER::NONE;
		int                                   _freq;
		 
	//	map< pair<string, string>, bool> isAddL2ForL3; // rec_sat_bool
	//	map< pair<string, string>, bool> isAddL3ForL2; // rec_sat_bool
	};


}
#endif // !GUPDATEPARALL_H
