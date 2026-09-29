/**
*
* @verbatim
History
-1.0  2020-06-11  Lvhb created (Network Precise Position,NPP)
@endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file		    gppprtkServer.cpp
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


#ifndef GNPPCLIENT_H
#define GNPPCLIENT_H

#include "gnpp/gposproc.h"
#include "gdata/gaug.h"
#include "gnpp/gInterpol.h"

namespace great
{
	class INS_LIBRARY_EXPORT t_gnppClient; typedef shared_ptr<t_gnppClient> t_spt_gnppClient;
	class INS_LIBRARY_EXPORT t_gnppClient
	{
	public:
		t_gnppClient(t_gsetbase* set, t_gallproc* data, t_glog* log);

		/** @brief default destructor. */
		virtual ~t_gnppClient();

		/**
		* @brief Batch processing loop over epoch.
		* @param[in] data  All data PPP need such as obs sp3 etc.
		* @param[in] beg   Loop begin time.
		* @param[in] end   Loop end time.
		* @return	 True  - compute successfully
		*            False - failure
		*/
		bool ProcessBatch(const t_gtime& begT, const t_gtime& endT) ;

	private:

		/* For post-processing testing,wb */
		void _prtAugFileHeader(const t_gtime& begT, const t_gtime& endT, string site);

	private:
		t_gsetbase* _set; t_gallproc* _alldata; t_glog* _log;
		t_gallobj* _allobj;
		t_gaug* _allaug;							///< Augmentation
		bool _isBase;								///< For DD/RTK model
		bool _isCompAug;							///< 目前强制设为true by lvhb, 202010
		bool _isWLProc;
		string _site;
		string _site_vrs;			///< For NRTK
		NPP_MODEL _nppmodel;///< not using now

		int _nipSite;								///<  site number
		vector<string> _list_ipsite;				///<  sites name

		t_spt_gInterpol _gInterpol;

		t_spt_gposproc _Proc;						///<  rtkflt/ppplsq
		t_spt_gposproc _Proc_WL;					///<  ppplsq

		//for one epoch augmentation information
		vector<t_gsatdata> _data_rover,_data_base;

		t_giof* _augfile = nullptr;


	};
}

#endif