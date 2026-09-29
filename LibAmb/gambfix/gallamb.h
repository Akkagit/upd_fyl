/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gallamb.h
* @brief	function for ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#ifndef GALLAMB_H
#define GALLAMB_H 

#include <tuple>
#include "gutils/ExportAMB.h"
#include "gio/glog.h"
#include "gutils/gnss.h"
#include "gmodels/gpar.h"
#include "gambfix/gambOw.h"
#include "gall/gallpar.h"
#include "gall/gallobs.h"
#include "gall/gallbias.h"
#include "gdata/gupd.h"
#include "gutils/gmutex.h"
#include "gset/gsetamb.h"
#include "gmodels/gprecisemodel.h"

using namespace std;
using namespace gnut;

namespace great
{
	class AMB_LIBRARY_EXPORT t_gallambs
	{
		
	public:
		typedef map<string, set<string> >  t_rec_sats;
		typedef pair<GSYS, t_gamb::AMB_ID> t_dd_sys_index;
		typedef pair<string, string>       t_dd_bl_index;
		typedef map<t_dd_sys_index, map<string, vector<t_gamb_Sd>> > t_d;
		typedef map<t_dd_sys_index, map<t_dd_bl_index, vector<t_gamb_Dd>> > t_dd;
		enum class AMB_COMB_OPTION 
		{ 
			CMB_REC_SD, // ²âÕ¾µ¥²î
			CMB_SAT_SD, // ÎÀÐÇµ¥²î
			CMB_REC_DD, // ÎÀÐÇµ¥²î»ù´¡ÉÏ×ö²âÕ¾Ë«²î
			CMB_SAT_DD, // ²âÕ¾µ¥²î»ù´¡ÉÏ×öÎÀÐÇË«²î
			UNDEF
		};

		t_gallambs();
		t_gallambs(t_gsetbase* gset, t_glog* log);
		virtual ~t_gallambs();


		void add_set(t_gsetbase* gset);
		void add_log(t_glog* log);
		void add_data(t_gallpar* pars, t_gallobs* obs, t_gallbias* allbias,t_gallobj* obj=nullptr);
		void add_data(t_gallpar* pars, t_gallobs* obs, t_gallbias* allbias, t_gupd* upd);
		void add_id(const t_gamb::AMB_ID& ID);
		void add_WL_bdeci(t_gbdeci* WL);
		void add_NL_bdeci(t_gbdeci* NL);
		t_d get_Sds();

		t_dd get_Dds();
		set<t_gamb::AMB_ID> IDs();
		bool _add_upd(string type,t_gamb_Sd& sd);
		bool cmb_amb(const AMB_COMB_OPTION& opt, const t_gamb::AMB_ID& ID, const vector<pair<string, string>>& baselines, double shortest_arc);
		bool cmb_amb(const AMB_COMB_OPTION& opt, const t_gamb::AMB_ID& ID, const vector<pair<string, string>>& baselines, 
			double shortest_arc, t_gsetamb::UPD_MODE upd_mode);

	public:
		map < pair<GSYS, t_gamb::AMB_ID>, int> all_WL_num;
		map < pair<GSYS, t_gamb::AMB_ID>, int> fix_WL_num;
		map < pair<GSYS, t_gamb::AMB_ID>, int> all_NL_num;
		map < pair<GSYS, t_gamb::AMB_ID>, int> fix_NL_num;
	protected:
		// for new 
		bool _cmb_intersat_Sd_amb(const pair<string, string>& sat_baseline, const  t_gamb::AMB_ID& ID,  double shortest_arc);
		bool _cmb_interrec_Dd_amb(const pair<string, string>& rec_baseline, const  t_gamb::AMB_ID& ID,  double shortest_arc);
		
		bool _cmb_intersat_Sd_amb(const pair<string, string>& sat_baseline, const  t_gamb::AMB_ID& ID1, const  t_gamb::AMB_ID& ID2, double shortest_arc);
		bool _cmb_interrec_Dd_amb(const pair<string, string>& rec_baseline, const  t_gamb::AMB_ID& ID1, const  t_gamb::AMB_ID& ID2, double shortest_arc);

		void _get_ows_table_from_pars();
		void _get_obs_table_from_gobs();

		void _get_dependent_IF_Sd(const t_gamb::AMB_ID& ID);
		void _get_dependent_IF_Sd2(const t_gamb::AMB_ID& ID);
		bool _checkAmbDepend(bool isFirst, int iNamb, int *iNdef, int iN_oneway, int *arriIpt2ow, int iMaxamb_ow, int iMaxamb_for_check);
		void _get_dependent_IF_Dd(const t_gamb::AMB_ID& ID);
		void _get_dependent_1X_Dd(const t_gamb::AMB_ID& ID);

		void _get_fix_rate(const t_gamb::AMB_ID& ID);
		void _get_SDfix_rate(const t_gamb::AMB_ID& ID);
		void math_dependent_IF_Dd(vector<t_gamb_Dd>& input, int num_Ow_1st, int num_Ow_2nd,                                vector<t_gamb_Dd>& result, Matrix& vector_space);
		void math_dependent_1X_Dd(vector<t_gamb_Dd>& input, int num_Ow_1st, int num_Ow_2nd, int num_Ow_3rd, int num_Ow_4th, vector<t_gamb_Dd>& result, Matrix& vector_space);
		void _select_Dds_value(const t_gamb::AMB_ID& ID);
		void _select_Sds_value();
		void _get_Dds_value(const t_gamb::AMB_ID& ID);
		void _get_Sds_value(t_gsetamb::UPD_MODE upd_mode);
		//void _get_Dds_value(const t_gamb::AMB_ID& ID1, const t_gamb::AMB_ID& ID2);

		bool _get_NL_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& Dd);
		bool _get_IF_NL_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& Dd);
		bool _get_1X_NL_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& Dd);
		
		bool _get_MW_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& Dd);

		//void _get_WL_value(const pair<GSYS, t_gamb::AMB_ID>& ID1, const pair<GSYS, t_gamb::AMB_ID>& ID2);
		bool _get_WL_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& Dd);
		//bool _get_WL_value(const pair<GSYS, t_gamb::AMB_ID>& index_key1, const pair<GSYS, t_gamb::AMB_ID>& index_key2, const pair<string, string>& bl, t_gamb_Dd& amb1, t_gamb_Dd& amb2);

		bool _get_IF_Dd_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& Dd);
		bool _get_1X_Dd_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& Dd);
		bool _get_LX_Dd_value(const pair<GSYS, t_gamb::AMB_ID>& index_key, const pair<string, string>& bl, t_gamb_Dd& Dd);


		void _get_MW(const string& rec, const string& sat,const int& freq, const set<t_gtime>& time_table, set<t_gtime>::iterator beg, set<t_gtime>::iterator end,
			double& value, double& sigma);
		bool _getSingleUpd(string mode, t_gtime t, string sat, double& value, double& sigma);
		bool _getSingleUpd(string mode, t_gtime t, string sat, double& value, double& sigma, double& rnl);
		void _setWaveLen(map<string, map<string, double>>& sys_wavelen);
		set< pair<GSYS, t_gamb::AMB_ID>> get_Dds_index_key(const t_gamb::AMB_ID& ID);
	private:
		t_glog*         _glog  = nullptr;
		t_gsetbase*		_gset  = nullptr;
		t_gallobs*      _gobs  = nullptr;
		t_gallbias*     _gbias = nullptr;
		t_gallpar*      _gpars = nullptr;
		t_gallobj*      _gobj  = nullptr;
		t_gupd*         _gupd  = nullptr;
		double*         _pdE   = nullptr;
		double*         _pdC   = nullptr;

		t_gbdeci*       _evaluator_MW = nullptr;
		t_gbdeci*       _evaluator_NL = nullptr;
		t_gprecisemodel* _gprecise = nullptr;


		set<t_gamb::AMB_ID> _IDs;
		DD_MODEL    _dd_mode;
		DCB_MODEL  _dcb_model;
		FREQ_MODEL _freq_model;
		map< GSYS, map<FREQ_SEQ, GOBSBAND> > _band_index;
		map< GSYS, map<GOBSBAND, FREQ_SEQ> > _freq_index;
		int _frequency;

		set<string>     _rec_list;
		set<string>     _sat_list;
		string          _class_id = "t_gallambs";	
		t_gmutex        _gmutex;

		set<t_gtime>     time_table;
		map<GSYS, map<string, map<pair<string, t_gtime>, map<int, double>>>>                            _MW_used;  // sys, rec, sat, time
		map<pair<GSYS, t_gamb::AMB_ID>, map<string, vector<t_gamb_Ow>> >                                _nOw_ambs; // map< pair<sys,obs_type>, map< rec, vector<ow>> >
		map<pair<GSYS, t_gamb::AMB_ID>, map<string, vector<t_gamb_Sd>> >                                _nSd_ambs; // map< pair<sys,obs_type>, map< rec, vector<sd>> >
		map<pair<GSYS, t_gamb::AMB_ID>, map<string, vector<t_gamb_Sd>> >                                _nSd_ambs_fix_WL; // map< pair<sys,obs_type>, map< rec, vector<sd>> >
		map<pair<GSYS, t_gamb::AMB_ID>, map<string, vector<t_gamb_Sd>> >                                _nSd_ambs_fix_NO; // map< pair<sys,obs_type>, map< rec, vector<sd>> >
		map<pair<GSYS, t_gamb::AMB_ID>, map<pair<string, string>, vector<t_gamb_Dd>> >                  _nDd_ambs; // map< pair<sys,obs_type>, map< rec1,rec2,sat1,sat2, vector<Dd>> >
		map<pair<GSYS, t_gamb::AMB_ID>, map<pair<string, string>, vector<t_gamb_Dd>> >                  _nDd_ambs_fix_WL; // map< pair<sys,obs_type>, map< rec1,rec2,sat1,sat2, vector<Dd>> >
		map<pair<GSYS, t_gamb::AMB_ID>, map<pair<string, string>, vector<t_gamb_Dd>> >                  _nDd_ambs_fix_NO; // map< pair<sys,obs_type>, map< rec1,rec2,sat1,sat2, vector<Dd>> >

	};

	bool sort_Dd(const t_gamb_Dd& A, const t_gamb_Dd& B);
	bool sort_Sd(const t_gamb_Sd& A, const t_gamb_Sd& B);
}

#endif