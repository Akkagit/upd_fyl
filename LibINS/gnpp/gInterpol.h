/**
*
* @verbatim
	History
	 -1.0 Lvhb  2020-06-15  created (Network Precise Position,NPP)
  @endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file		    gposproc.h
* @brief		virtual class for process
*
* @author       Lvhb, Wuhan University
* @version		1.0.0
* @date		    2020-06-15
*				2020-10-24 lvhb added NRTK/URT interpolation for testing
*				2020-11-06 lvhb modified interpolation for NPP+INS
*		
*/


#ifndef GINTERPOL_H
#define GINTERPOL_H

#include "gutils/ExportINS.h"
#include "gio/gfile.h"

#include "gall/gallproc.h"
#include "gio/glog.h"
#include "gset/gsetbase.h"
#include "gall/gallobj.h"
#include "gdata/gaug.h"
#include "gall/gallobs.h"

namespace great
{

	/**
	* @brief virtual class for interpolation method
	*/
	class  INS_LIBRARY_EXPORT t_gInterpol; typedef shared_ptr<t_gInterpol> t_spt_gInterpol;
	class INS_LIBRARY_EXPORT t_gInterpol
	{
	public:
		t_gInterpol(t_gsetbase* set, t_gallproc* data, t_glog* log);
		virtual ~t_gInterpol();
	
		/** @brief  Get Aug-data for NRTK/URTK/PPP-RTK */
		bool interpolAug(const t_gtime& nowT, vector<t_gsatdata>& data_rover, vector<t_gsatdata>& data_base);
		void printAugFile(t_giof* augfile);

	private:
		void _setCrd(ColumnVector& dx, ColumnVector& dy, ColumnVector& dz);

		double _MLCM(ColumnVector& AugData);			///< for NRTK/URTK/PPP-RTK,lvhb
		double _URTKInterpol(ColumnVector& dAugData);	///< for URTK: refer zou xuan paper, only support 3 sites,lvhb in 202010
		double _NRTKInterpol(ColumnVector& dAugData);	///< for NRTK:twm, only support 3 sites,lvhb test in 202010
	
		ColumnVector _crd_dx;
		ColumnVector _crd_dy;
		ColumnVector _crd_dz;
		ColumnVector _coef_a;

	private:	///< Lvhb: mainly for NRTK solution about following functions

		/* @brief Select major site according baseline length,Just for NRTK,lvhb */
		bool _selMajorSite(const t_gtime& nowT);

		/* @brief Check whether current VRS data exist,Just for NRTK,lvhb */
		t_spt_gobs _check_vrs_prn_pt(const t_gtime& nowT, string& sat);

		/** @brief Get SPP approximately coordinates */
		virtual bool _SppProcOneEpoch(const t_gtime& nowT, const string& site, t_gtriple* crd, double* reclk, map<string, t_gsatdata>* satdata = NULL);

		/** @brief Calculate satellite position/clk/rho/elevation/azel,lvhb added for NRTK */
		void _cal_sat_inf(t_gsatdata& obs_sat, const t_gtriple& xyz_r);

		/** @brief Get Tropospheric delay according SAAS or other model,lvhb added for NRTK(URTK) */
		double _getTropModValue(const t_gtime& nowT, t_gsatdata& obsdata);

		virtual bool _avaliable(t_gtime nowT, string* site = NULL, string* site_base = NULL);

		/** @brief Judge interpol data range (>= 3 sites)
		  _Rectangle: 1:x_max 2: x_min; 3:y_max; 4:y_min
		* @Note not using now
		* @param[in] base site coordinate
		* @return	 True  - successfully
		*            False - failure
		*/
		bool _isInterpolRange(const t_gtriple* crdT = NULL);

		/*  @brief Not using now by as a result of using comprehensive augmentation in current NRTK*/
		double _corrTropHeight(const t_gtime& nowT, int sign, const t_gtriple& bcrd, t_gtriple& rcrd, const t_gtriple& crd,
			t_gsatdata bdata_sat, t_gsatdata bdata_bsat, t_gsatdata rdata_sat, t_gsatdata rdata_bsat);
		
		/** @brief Calculate rho/ele according to positions of Sat/Site,lvhb added to for NRTK */
		void cal_rho_azel(const string& site_name, t_gtriple& xyz_s, const t_gtriple& xyz_r, t_gsatdata& obs_sat);

	private:
		t_gsetbase* _set; t_gallproc* _alldata; t_glog* _log;
		t_gallobj* _allobj;
		t_gaug* _allaug;							///< Augmentation
		bool _isBase;								///< For DD/RTK model
		bool _isCompAug;							///< 目前强制设为true by lvhb, 202010
		string _site;
		string _strMajorSite, _site_vrs;			///< For NRTK

		int _nipSite;								///<  site number
		vector<string> _list_ipsite;				///<  sites name

		//for one epoch augmentation information
		vector<t_gsatdata> _data_rover, _data_base;
		t_gtriple _crd_rover;
		ColumnVector _Rectangle;					///< 1:x_max 2: x_min; 3:y_max; 4:y_min

		/* For one epoch processing,lvhb created for URTK/PPP-RTK SSR */
		//map<int, map<string, map<GOBSBAND, pair<double, double>>>> _AllCorr; ///< for one epoch: first-range,second-phase //int site,string sat
		//map<int, map<string, map<GOBSBAND, double>>> _AllIonAug;			 ///< int site,string sat
		//map<int, map<string, double>> _AlltropAug;							 ///< int site,string sat

		t_giof* _augfile = nullptr;
	};

}
#endif