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

#ifndef GUPDATEPARIRC_H
#define GUPDATEPARIRC_H

#include <map>
#include "gproc/gupdatepar.h"

namespace great {
	class BASE_LIBRARY_EXPORT t_gupdateparIRC :public t_gupdatepar
	{
	public:
		t_gupdateparIRC(shared_ptr<t_gcycleslip> cycleslip);
		virtual ~t_gupdateparIRC();

	protected:

		virtual void _update_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata,t_gupdateparinfo& update_info) override;
	


	};


}
#endif // !GUPDATEPARIRC_H
