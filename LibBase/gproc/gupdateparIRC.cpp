#include "gproc/gupdateparIRC.h"
#include "gutils/gturboedit.h"

namespace great
{
	t_gupdateparIRC::t_gupdateparIRC(shared_ptr<t_gcycleslip> cycleslip)
	{
	}
	t_gupdateparIRC::~t_gupdateparIRC()
	{
	}
	void t_gupdateparIRC::_update_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo & update_info)
	{
		// Don't need remove amb parameters
		return;
	}
}