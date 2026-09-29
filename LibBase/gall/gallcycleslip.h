#ifndef GALLCYCLESLIP_H
#define GALLCYCLESLIP_H

#include "gutils/gcycleslip.h"
#include "gall/gallambflag.h"

using namespace gnut;
namespace great
{
	class t_gallcycleslip
	{
	public:
		t_gallcycleslip();
		t_gallcycleslip(t_gsetbase* set, t_glog* log);
		virtual ~t_gallcycleslip();

	public:
	protected:
		void read_amb_info_files(const FREQ_SEQ& freq_1, const FREQ_SEQ& freq_2, const vector<string>& _amb_info_files, t_gallambflag* _all_ambflag);

	protected:
		t_glog*      _log;
		t_gsetbase*  _set;
		map<pair<FREQ_SEQ, FREQ_SEQ> , SLIPMODEL>               _slip_model;
		map< pair<FREQ_SEQ, FREQ_SEQ>, shared_ptr<t_gcycleslip>> _amb_info;
		set<string>                                              _rec_list;
		map< pair<FREQ_SEQ, FREQ_SEQ>, set<string>>              _amb_info_files; 
	};

}
#endif // !GALLCYCLESLIP_H
