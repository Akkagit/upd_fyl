/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambOw.h
* @brief	one way ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#ifndef GAMBOW_H
#define GAMBOW_H 


#include "gutils/ExportAMB.h"
#include "gio/glog.h"
#include "gutils/gnss.h"
#include "gmodels/gpar.h"
#include "gambfix/gbdeci.h"
#include <map>
using namespace std;
using namespace gnut;

namespace great
{
	class AMB_LIBRARY_EXPORT t_gamb
	{
	public:
		enum class AMB_ID 
		{ 
			IF_L1_L2,    // 无电离层12   频率
			IF_L1_L3,	 // 无电离层13   频率，只有在12频存在的时候才使用来标定13频率无电离层组合
			IF_L1_L2_L3, // 无电离层12+13频率，只有在三频联合固定的时候才使用
			ALL_L1,      // 非组合 L1
			ALL_L2,      // 非组合 L2
			ALL_L3,      // 非组合 L3
			ALL_L4,      // 非组合 L4
			ALL_L5,      // 非组合 L5
			ALL_L1_L2,   // 非组合 L1 + L2     , 求双差时采用
			ALL_L1_L3,   // 非组合 L1 + L3     , 求双差时采用
			ALL_L1_L2_L3,// 非组合 L1 + L2 + L3, 求联合三频双差时采用 
			WL_L1_L2, 
			WL_L1_L3, 
			UNDEF 
		};

		static GOBS_LC    PAR2LC(const t_gpar::t_type& ambtype);
		static GOBSBAND  LC2BAND(const GSYS&   gsys, const GOBS_LC& LC);
		static double    LC2FREQ(const string& sat,  const GOBS_LC& LC);
		static string    ambId2str(AMB_ID ID);

		t_gamb();
		t_gamb(AMB_ID ID, string sat, string rec, t_gtime beg, t_gtime end);
		t_gamb(const t_gamb& amb);
		virtual ~t_gamb();

		AMB_ID        ID   = AMB_ID::UNDEF;       // IF_L1L2, ALL_L1, ALL_L2, ALL_L3, ALL_L4, ALL_L5, UNDEF
		GOBS_LC       freq = GOBS_LC::LC_UNDEF;   // LC_UNDEF, LC_L1, LC_L2, LC_L3, LC_L4, LC_L5, LC_IF, LC_MW, LC_NL, LC_WL, LC_GF
		t_gtime       beg;
		t_gtime       end;

		string		  sat;
		string		  rec;
		GSYS		  gsys;

		double        value = 0.0;
		double        sigma = 0.0;
		bool	      fixed = false;
	};

	class AMB_LIBRARY_EXPORT t_gamb_Ow : public t_gamb
	{
	public:
		t_gamb_Ow(const t_gamb::AMB_ID& ID, const string& sat, const string& rec, const t_gtime& beg, const t_gtime& end);
		t_gamb_Ow(const t_gamb_Ow& ow);
		virtual ~t_gamb_Ow();
	};


	// 单差模糊度
	class AMB_LIBRARY_EXPORT t_gamb_Sd : public t_gamb
	{
	public:
		t_gamb_Sd();
		t_gamb_Sd(const int& idx1, const int& idx2, const t_gamb::AMB_ID& ID);
		t_gamb_Sd(const map<int, int>& idx1, const map<int, int>& idx2, const t_gamb::AMB_ID& ID);
		t_gamb_Sd(const t_gamb_Sd& Sd);
		virtual ~t_gamb_Sd();

		bool       NL_fixed = false;
		bool       MW_fixed = false;
		double     NL_value = 0.0;
		double     NL_sigma = 0.0;
		double     MW_value = 0.0;
		double     MW_sigma = 0.0;
		double     rnl = 0.0;

		
		map<int, int>            idx2ow_1st;
		map<int, int>            idx2ow_2nd;

		string                   sgl_obj;
		string                   dbl_obj_1st;
		string                   dbl_obj_2nd;
		map<int, double>                 Ow_value;
		map<int, pair<t_gtime, t_gtime>> Ow_time;
	};

	// 双差模糊度
	class AMB_LIBRARY_EXPORT t_gamb_Dd : public t_gamb
	{
	public:
		t_gamb_Dd() {};
		t_gamb_Dd(const t_gamb::AMB_ID& ID, const int& idx_1st, const int& idx_2nd, const int& idx_3rd, const int& idx_4th);
		t_gamb_Dd(const t_gamb::AMB_ID& ID, const map<int, int>& idx_1st, const map<int, int>& idx_2nd);
		t_gamb_Dd(const t_gamb_Dd& Dd);
		virtual ~t_gamb_Dd();


		bool      isNormal  = false;
		bool       NL_fixed = false;
		bool       MW_fixed = false;
		bool      EMW_fixed = false;
		bool       WL_fixed = false;
		bool      EWL_fixed = false;

		double     NL_value = 0.0;
		double     NL_sigma = 0.0;

		double     WL_value = 0.0;
		double     WL_sigma = 0.0;

		double    EWL_value = 0.0;
		double    EWL_sigma = 0.0;

		double     MW_value = 0.0;
		double     MW_sigma = 0.0;

		double    EMW_value = 0.0;
		double    EMW_sigma = 0.0;

		double     rnl = 0.0;
		map<t_gamb::AMB_ID, double>  DD_LX_value;
		map<t_gamb::AMB_ID, bool>    DD_LX_fixed;
		map<t_gamb::AMB_ID, double>  DD_LX_sigma;

		map<int, int>       idx2ow_1st;
		map<int, int>       idx2ow_2nd;

		string sat_1st;
		string sat_2nd;
		string rec_1st;
		string rec_2nd;

		map<int, double>                 Ow_value;
		map<int, double>                 Mw_value;
		map<int, double>                EMw_value;
		map<int, double>                 Mw_sigma;
		map<int, double>                EMw_sigma;
		map<int, pair<t_gtime, t_gtime>> Ow_time;
	};

}

#endif