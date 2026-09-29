/**
*
* @verbatim
	History
	 -1.0 Lvhb  2020-06-15  created
  @endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file		    gposproc.h
* @brief		virtual class for process
*
* @author       Lvhb, Wuhan University
* @version		1.0.0
* @date		    2020-06-15
*
*/


#ifndef GPOSPROC_H
#define GPOSPROC_H

#include "gproc/gppplsq.h"
#include "gproc/gpvtflt.h"

//#define DEBUG_BPL_LOG

namespace great
{

	/**
	* @brief virtual class for npp process
	*/
	class INS_LIBRARY_EXPORT t_gposproc; typedef shared_ptr<t_gposproc> t_spt_gposproc;
	class INS_LIBRARY_EXPORT t_gposproc
	{
	public:
		t_gposproc(string site, string site_base, t_gsetbase* set, t_gallproc* data, t_glog* log);
		virtual ~t_gposproc();

		/** @brief Initializing NPP setting,lvhb created in 20200615 */
		virtual bool InitProc(const t_gtime& begT, const t_gtime& endT);
		
		/** @brief Preprocessing for FLT/LSQ,mainly about cycle slip detection */
		virtual bool Preprocess(t_gtime& begT, t_gtime& endT, double& sampling, vector<t_gsatdata>* data_rover = NULL, vector<t_gsatdata>* data_base = NULL, bool doppler_check = false);
		
		/** @brief One epoch processing of RTK/PPP */
		virtual int  ProcessOneEpoch(const t_gtime& nowT, vector<t_gsatdata>* data_rover=NULL, vector<t_gsatdata>* data_base=NULL);
		
		/** @brief main process batch */
		virtual bool ProcessBatch( const t_gtime& beg, const t_gtime& end);//no use for now
	
		//////////////////////////////////////////////////////////////////////////
	public:
		/** @brief For real-time processing */
		virtual bool avaliable(t_gtime nowT,string* site = NULL,string* site_base =NULL);

		/** @brief Get UC RTK/PPP parameters */
		virtual t_gallpar GetParam(const t_gtime& nowT, set<string>* satref = NULL);
		
		/** @brief Get SPP approximately coordinates */
		virtual bool SppProcOneEpoch(const t_gtime& nowT, const string& site, t_gtriple* crd, double* reclk, map<string, t_gsatdata>* satdata = NULL);

		//using for ppplsq wb
		void setWLParams(t_gallpar& params);
		void setClient(bool isClient);
		void setObsMode(OBSCOMBIN mode);
		void resetLsqUpdate(t_gsetbase* set, t_gallproc* data, t_glog* log);
		virtual bool _ExtractResfile();

	private:
		bool _isBase;			//For DD/RTK model
		bool _isFLT;			//For filter estimation
		string _strEst;
		string _site, _site_base; 
		t_gsetbase* _set; 
		t_gallproc* _alldata;
		t_glog* _log;

		//process
		shared_ptr<t_gppplsq> _ppplsq;
		shared_ptr<t_gpvtflt> _pvtflt;
		//preprocess
		shared_ptr<t_gexeturboedit> _gturboedit;
		shared_ptr<t_gpreproc> _gpre;
		//model
		shared_ptr<t_gbiasmodel> _precise_bias;
	};


}
#endif