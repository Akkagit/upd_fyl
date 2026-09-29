

#ifndef GAMBFIXD_H
#define GAMBFIXD_H 

#include "gutils/ExportAMB.h"
#include "gambfix/gambfix.h"
#include "gset/gsetamb.h"
#include "gdata/gambcon.h"

using namespace std;

namespace great
{

	/**
	*@ brief class for LAMADA search finding fixing ambiguity.
	*/
	class AMB_LIBRARY_EXPORT t_gambfixD : public t_gambfix
	{
	public:
		t_gambfixD();
		t_gambfixD(t_glog* log, t_gsetbase* set);
		
		virtual ~t_gambfixD();

		
		// 确定beg到end时间内的单差模糊度
		bool ProcessBatch();
		
		// 确定测站基线
		vector<pair<string, string>> baselines();

	protected:
		t_gamb_baseline_recs  _gRecBl;
		t_gamb_baseline_sats  _gSatBl;

		// 通过测站/卫星列表确定地面测站的基线
		bool _screen_rec_and_sat(t_gallpar* pars, t_gsetbase* set, t_gallobj* obj);
		
		void _write_ambcon();
	private:
		string _class_id = "t_gambfixD";
	};
}

#endif