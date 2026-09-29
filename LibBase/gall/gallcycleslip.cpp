#include "gall/gallcycleslip.h"

namespace great
{
	great::t_gallcycleslip::t_gallcycleslip()
	{

	}

	great::t_gallcycleslip::t_gallcycleslip(t_gsetbase* set, t_glog* log)
	{
		//_slip_model = dynamic_cast<t_gsetproc*>(set)->slip_model();
		_set = set;
		_log = log;
	}

	great::t_gallcycleslip::~t_gallcycleslip()
	{
	}

	void great::t_gallcycleslip::read_amb_info_files(const FREQ_SEQ& freq_1, const FREQ_SEQ& freq_2, const vector<string>& _amb_info_files, t_gallambflag* _all_ambflag)
	{

	}
}

