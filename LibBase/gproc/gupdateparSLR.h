/**
*
* @verbatim
	History
	 -1.0 ZhangHongmin  2020-05-14  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gupdateparSLR.h
* @brief		update parameter of SLR
*
* @author       ZhangHongmin, Wuhan University
* @version		1.0.0
* @date		    2020-05-14
*
*/

#ifndef GUPDATESLR_H
#define GUPDATESLR_H

#include <map>
#include "gproc/gupdatepar.h"
#include "gall/gallproc.h"
#include "gdata/gpoleut1.h"

namespace great {
	class BASE_LIBRARY_EXPORT t_gupdateparSLR :public t_gupdatepar
	{
	public:
		t_gupdateparSLR(t_gallproc* data);
		virtual ~t_gupdateparSLR();

	protected:

		bool _update_erp_par(const t_gtime& epoch, int update_id, t_gallpar& allpars, t_gupdateparinfo& update_info);
		t_gupdateparinfo get_all_update_parameters(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata);
		void _update_process_pars(const t_gtime& epoch, t_gallpar& allpars, t_gupdateparinfo & update_info);
		void _update_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo & update_info);

	protected:
		map<pair<string, string>, int> _amb_id;
		t_gallproc*      _alldata = nullptr;
		t_gpoleut1*      _allerp = nullptr;
		map<string, double> data1, data2;

	};
}

#endif