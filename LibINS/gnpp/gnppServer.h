/**
*
* @verbatim
History
-1.0  2020-06-11 lvhb created (Network Precise Position,NPP)
@endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file		    gnppServer.cpp
* @brief		main control class for ppp/rtk estimate in network Server
*
* @author       lvhb,BoWang Wuhan University
* @version		1.0.0
				1.1.0 lvhb finished NRTK testing (Comprehensive augmentation) in 2020-10-15
				1.1.1 lvhb added and finished URTK testing (Comprehensive augmentation) in 2020-10-24
				1.1.2 lvhb tested in real-time TCP/IP solution about NRTK/URTK in 2020-11-10
* @date		    2020-06-11
*
*/


#ifndef GNPPSERVER_H
#define GNPPSERVER_H

#include "gutils/ExportINS.h"
#include "gnpp/gposproc.h"

namespace great
{
	class INS_LIBRARY_EXPORT t_gnppServer; typedef shared_ptr<t_gnppServer> t_spt_gnppServer;
	class INS_LIBRARY_EXPORT t_gnppServer
	{
	public:
		t_gnppServer(t_gsetbase* set, t_gallproc* data, t_glog* log);

		/** @brief default destructor. */
		virtual ~t_gnppServer();

		/**
		* @brief Batch processing loop over epoch.
		* @param[in] data  All data PPP need such as obs sp3 etc.
		* @param[in] beg   Loop begin time.
		* @param[in] end   Loop end time.
		* @return	 True  - compute successfully
		*            False - failure
		*/
		bool ProcessBatch(const t_gtime& begT, const t_gtime& endT);

	private:

		//1、Initialization for NPP solution
		bool _InitProcAll(const t_gtime& begT, const t_gtime& endT);

		//2、Get UD information according to DD ambiguity, which just for NRTK/URTK,created by lvhb
		bool _selMajorSite(const t_gtime& nowT, vector<int>& ind_state);//for NRTK/URTK
		bool _getOneBLuducParam(const t_gtime& nowT, const int& i, bool bmajor_site);//for NRTK/URTK;
		bool _getAllBLuducParam(const t_gtime& nowT, vector<int>& ind_state);//for NRTK/URTK;

		//3、Get augmentation data for NRTK/URTK/PPP-RTK,created by lvhb
		bool _getAugData(const t_gtime& nowT, vector<int>& ind_state);

		//3.1、Not using now, lvhb added for NRTK/URTK/PPP-RTK in 202010
		bool _getUDIonAug(const int& i, t_gsatdata& data, /*t_gsatdata& data_base,*/ map<FREQ_SEQ, GOBSBAND>& crt_bands, t_gallpar& param);
		bool _getUDTropAug(const t_gtime& now, const int& i, t_gsatdata& data, t_gallpar& param);

		//4、Output augmentation data for post processing,（wb）
		void _prtAugFileHeader(const t_gtime& begT, const t_gtime& endT, string site);

	private:
		t_gsetbase* _set; t_gallproc* _alldata; t_glog* _log;
		bool _isBase;    ///< NRTK/URTK 是否利用基站数据
		bool _isCompAug; ///< 综合误差, 目前强制使用综合误差，lvhb
		NPP_MODEL _nppmodel;
		shared_ptr<t_gbasemodel> _base_model = nullptr;

		string _strMajorSite;												 ///< For NRTK/URTK,lvhb
		string _strMajorSite_pre;											 ///< For NRTK/URTK,added by lvhb, 20201025

		int _nSite;
		vector<string> _list_site;											 ///< rover site name
		vector<string> _list_site_base;										 ///< base site name for NRTK/URTK
		vector<string> _list_rec_bl;										 ///< all fix site name for NRTK/URTK
		map<int, t_spt_gposproc> _mapProc;									 ///< map rtkflt/ppplsq
		map<string, shared_ptr<t_giof>> _augfile;
		shared_ptr<t_gtcp> _augtcp;  // zhshen

		/* For one epoch augmentation information */
		map<int, map<string, map<GOBSBAND, pair<double, double>>>> _AllCorr; ///< for one epoch: first-range,second-phase //int site,string sat
		map<int, map<string, map<GOBSBAND, double>>> _AllIonAug;			 ///< int site,string sat,band
		map<int, map<string, double>> _AlltropAug;							 ///< int site,string sat
		
		/* For NRTK/URTK in every epoch,created by lvhb */
		map<string, map<string, map<GOBSBAND, double>>> _Alluducamb;		///< string site, sat, for NRTK/URTK [m]
		map<string, t_gallpar>	_Alluducpar;								///< string site, for NRTK/URTK
		map<int, t_gallpar> _Allsducpar;									///< int base-rover, for NRTK/URTK
	};
}

#endif