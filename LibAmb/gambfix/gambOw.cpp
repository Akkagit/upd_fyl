/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambOw.cpp
* @brief	function for ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#include "gambfix/gambOw.h"
#include "gutils/gtime.h"

namespace great
{
	t_gamb::t_gamb()
	{

	}

	GOBS_LC t_gamb::PAR2LC(const t_gpar::t_type& ambtype)
	{
		if (ambtype == t_gpar::AMB_IF) return GOBS_LC::LC_IF;
		else if (ambtype == t_gpar::AMB_L1) return GOBS_LC::LC_L1;
		else if (ambtype == t_gpar::AMB_L2) return GOBS_LC::LC_L2;
		else if (ambtype == t_gpar::AMB_L3) return GOBS_LC::LC_L3;
		else if (ambtype == t_gpar::AMB_L4) return GOBS_LC::LC_L4;
		else if (ambtype == t_gpar::AMB_L5) return GOBS_LC::LC_L5;
		else return GOBS_LC::LC_UNDEF;
	}

	GOBSBAND t_gamb::LC2BAND(const GSYS& gsys, const GOBS_LC& LC)
	{
		t_map_band map_band = GNSS_BAND_PRIORITY;
		return map_band[gsys][LC];
	}


	double t_gamb::LC2FREQ(const string& sat, const GOBS_LC& LC)
	{
		t_gobsgnss gnss;
		gnss.sat(sat);
		auto sys  = t_gsys::sat2gsys(sat);
		auto band = t_gamb::LC2BAND(sys, LC);
		return gnss.frequency(band);
	}

	string t_gamb::ambId2str(AMB_ID ID)
	{
		switch (ID)
		{
			    case AMB_ID::IF_L1_L2    : return "IF_L1_L2    ";
				case AMB_ID::IF_L1_L3	 : return "IF_L1_L3	   ";
				case AMB_ID::IF_L1_L2_L3 : return "IF_L1_L2_L3 ";
				case AMB_ID::ALL_L1      : return "ALL_L1      ";
				case AMB_ID::ALL_L2      : return "ALL_L2      ";
				case AMB_ID::ALL_L3      : return "ALL_L3      ";
				case AMB_ID::ALL_L4      : return "ALL_L4      ";
				case AMB_ID::ALL_L5      : return "ALL_L5      ";
				case AMB_ID::ALL_L1_L2   : return "ALL_L1_L2   ";
				case AMB_ID::ALL_L1_L3   : return "ALL_L1_L3   ";
				case AMB_ID::ALL_L1_L2_L3: return "ALL_L1_L2_L3";
				case AMB_ID::WL_L1_L2	 : return "WL_L1_L2	   ";
				case AMB_ID::WL_L1_L3	 : return "WL_L1_L3	   ";
				default                  : return "UNDEF	   ";
		}
	}

	t_gamb::t_gamb(AMB_ID ID, string sat, string rec, t_gtime beg, t_gtime end)
	{
		this->ID = ID;
		this->beg = beg;
		this->end = end;
		this->sat = sat;
		this->rec = rec;
		this->gsys = gsys;

		if (!this->sat.empty() && this->gsys != t_gsys::str2gsys(this->sat.substr(0, 1)))
		{
			this->gsys = t_gsys::str2gsys(this->sat.substr(0, 1));
		}
	}

	t_gamb::t_gamb(const t_gamb & amb)
	{
		this->ID = ID;
		this->freq = freq;
		this->beg = beg;
		this->end = end;
		this->sat = sat;
		this->rec = rec;
		this->gsys = gsys;
		this->value = value;
		this->sigma = sigma;
		this->fixed = fixed;
	}

	t_gamb::~t_gamb()
	{
	}


	t_gamb_Ow::t_gamb_Ow(const t_gamb::AMB_ID& ID, const string & sat, const string & rec, const t_gtime & beg, const t_gtime & end) :
		t_gamb(ID,sat,rec,beg,end)
	{		

	}

	t_gamb_Ow::t_gamb_Ow(const t_gamb_Ow & ow)
	{
		this->ID   = ow.ID;
		this->beg  = ow.beg;
		this->end  = ow.end;
		this->sat  = ow.sat;
		this->rec  = ow.rec;
		this->gsys = ow.gsys;
		this->value = ow.value;
		this->fixed = ow.fixed;
		this->freq = ow.freq;
		this->sigma = ow.sigma;
	}

	great::t_gamb_Ow::~t_gamb_Ow()
	{
	}

	t_gamb_Sd::t_gamb_Sd()
	{
		/*bool       NL_fixed = false;
		bool       MW_fixed = false;
		double     NL_value = 0.0;
		double     NL_sigma = 0.0;
		double     MW_value = 0.0;
		double     MW_sigma = 0.0;
		double     rnl = 0.0;*/
	}

	t_gamb_Sd::t_gamb_Sd(const int& idx1, const int& idx2, const t_gamb::AMB_ID& ID)
	{
		idx2ow_1st[1] = idx1;
		idx2ow_1st[2] = idx2;
		this->ID =  ID;
	}

	t_gamb_Sd::t_gamb_Sd(const map<int, int>& idx1, const map<int, int>& idx2, const t_gamb::AMB_ID& ID)
	{
		idx2ow_1st = idx1;
		idx2ow_2nd = idx2;
		this->ID = ID;
	}

	t_gamb_Sd::t_gamb_Sd(const t_gamb_Sd & Sd)
	{
		this->beg = Sd.beg;
		this->end = Sd.end;
		this->fixed = Sd.fixed;
		this->freq = Sd.freq;
		this->gsys = Sd.gsys;
		this->ID = Sd.ID;
		this->rec = Sd.rec;
		this->sat = Sd.sat;
		this->value = Sd.value;
		this->sigma = Sd.sigma;
		this->MW_fixed = Sd.MW_fixed;
		this->MW_value = Sd.MW_value;
		this->MW_sigma = Sd.MW_sigma;
		this->NL_fixed = Sd.NL_fixed;
		this->NL_value = Sd.NL_value;
		this->NL_sigma = Sd.NL_sigma;
		this->idx2ow_1st  = Sd.idx2ow_1st;;
		this->idx2ow_2nd  = Sd.idx2ow_2nd;;
		this->sgl_obj     = Sd.sgl_obj;
		this->dbl_obj_1st = Sd.dbl_obj_1st;
		this->dbl_obj_2nd = Sd.dbl_obj_2nd;
	}						

    t_gamb_Sd::~t_gamb_Sd()
	{

	}

	t_gamb_Dd::t_gamb_Dd(const t_gamb::AMB_ID& ID, const int& idx_1st, const int& idx_2nd, const int& idx_3rd, const int& idx_4th)
	{
		this->ID = ID;
		idx2ow_1st[1] = idx_1st;
		idx2ow_1st[2] = idx_2nd;
		idx2ow_1st[3] = idx_3rd;
		idx2ow_1st[4] = idx_4th;
	}

	t_gamb_Dd::t_gamb_Dd(const t_gamb::AMB_ID& ID, const map<int, int>& idx_1st, const map<int, int>& idx_2nd)
	{
		idx2ow_1st = idx_1st;
		idx2ow_2nd = idx_2nd;
		this->ID = ID;
	}

	t_gamb_Dd::t_gamb_Dd(const t_gamb_Dd & Dd)
	{
		this->NL_fixed = Dd.NL_fixed;
		this->MW_fixed = Dd.MW_fixed;
		this->WL_fixed = Dd.WL_fixed;
		this->NL_value = Dd.NL_value;
		this->NL_sigma = Dd.NL_sigma;
		this->WL_value = Dd.WL_value;
		this->WL_sigma = Dd.WL_sigma;
		this->MW_value = Dd.MW_value;
		this->MW_sigma = Dd.MW_sigma;
		this->DD_LX_value = Dd.DD_LX_value;
		this->DD_LX_fixed = Dd.DD_LX_fixed;
		this->DD_LX_sigma = Dd.DD_LX_sigma;
		this->idx2ow_1st = Dd.idx2ow_1st;
		this->idx2ow_2nd = Dd.idx2ow_2nd;
		this->sat_1st = Dd.sat_1st;
		this->sat_2nd = Dd.sat_2nd;
		this->rec_1st = Dd.rec_1st;
		this->rec_2nd = Dd.rec_2nd;
		this->Ow_value = Dd.Ow_value;
		this->Mw_value = Dd.Mw_value;
		this->Mw_sigma = Dd.Mw_sigma;
		this->Ow_time = Dd.Ow_time;
		this->ID = Dd.ID;
		this->freq = Dd.freq;
		this->beg = Dd.beg;
		this->end = Dd.end;
		this->sat = Dd.sat;
		this->rec = Dd.rec;
		this->gsys = Dd.gsys;
		this->value = Dd.value;
		this->sigma = Dd.sigma;
		this->fixed = Dd.fixed;
	}


	great::t_gamb_Dd::~t_gamb_Dd()
	{

	}

}





