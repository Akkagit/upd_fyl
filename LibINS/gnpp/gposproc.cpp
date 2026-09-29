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
#include "gposproc.h"

//#define DEBUG_LOG
//#define DEBUG_BPL_LOG
// Judge whether file exsit
#if  defined _WIN32 || defined _WIN64
#include <io.h>
#include <direct.h> 
#define ACCESS(fileName,accessMode) _access(fileName,accessMode)
#define MKDIR(path) _mkdir(path)
#else
#include <unistd.h>
#include <sys/stat.h>
#define ACCESS(fileName,accessMode) access(fileName,accessMode)
#define MKDIR(path) mkdir(path,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif


namespace great
{

	t_gposproc::t_gposproc(string site, string site_base, t_gsetbase* set, t_gallproc* data, t_glog* log)
		:_site(site),
		_site_base(site_base),
		_set(set),
		_alldata(data),
		_log(log),
		_isBase(false)
	{
		if (!site_base.empty()) _isBase = true;
		_strEst = dynamic_cast<t_gsetgen*>(set)->estimator();
		_isFLT = (_isBase || _strEst == "FLT");

		auto gobs = dynamic_cast<t_gallobs*>((*data)[t_gdata::ALLOBS]);
		auto gorb = dynamic_cast<t_gallnav*>((*data)[t_gdata::GRP_EPHEM]);

		_precise_bias = make_shared<t_gprecisebias>(_alldata, _log, _set);

		if (!_isFLT){
			_ppplsq = make_shared<t_gppplsq>(site, set, data, log);
		}
		else {
			_pvtflt = make_shared<t_gpvtflt>(site, site_base,set,data);
		}

		// set and prepare 
		if (dynamic_cast<t_gsetproc*>(set)->slip_model() == SLIPMODEL::TURBO_EDIT){
			_gturboedit = make_shared<t_gexeturboedit>(set, log, data);
		}
		else  {// DEFAULT
			_gpre = make_shared <t_gpreproc>(gobs, set);
			_gpre->glog(log);
			_gpre->setNav(gorb);
		}

	}

	t_gposproc::~t_gposproc()
	{
		//if (_log) { if (_log->is_open()) { _log->close(); }; delete _log; _log = nullptr; }
	}

	bool t_gposproc::InitProc(const t_gtime& begT, const t_gtime& endT)
	{
		if (_isFLT){
			_pvtflt->InitProc(begT,endT);
		}
		else {//lsq
			_ppplsq->InitProc(_alldata, begT, endT);
		}

		return true;
	}

	bool t_gposproc::Preprocess(t_gtime& begT,t_gtime& endT,double& sampling, vector<t_gsatdata>* data_rover, vector<t_gsatdata>* data_base, bool doppler_check)
	{
		try
		{
			//auto nppmodel = dynamic_cast<t_gsetproc*>(_set)->npp_model();
			bool isClient = dynamic_cast<t_gsetgen*>(_set)->isClient();
			// set and prepare 
			if (dynamic_cast<t_gsetproc*>(_set)->slip_model() == SLIPMODEL::TURBO_EDIT)
			{
				if (doppler_check) {
					_gturboedit->checkSlipWithDoppler(_site, begT, endT, sampling);
				}
				else {
					if(isClient) _gturboedit->setSingleEpoData(data_rover);
					_gturboedit->ProcessBatch(_site, begT, endT, sampling);
					if (_isBase) {
						if (!isClient) {
							_gturboedit->setSingleEpoData(data_base);
							_gturboedit->ProcessBatch(_site_base, begT, endT, sampling);
						}
					}
				}
			}
			else  // DEFAULT
			{
				if (isClient) _gpre->setSingleEpoData(data_rover);
				_gpre->preprocess(_site, begT, endT, sampling, false);
				if (_isBase) {
					if (!isClient){
						_gpre->setSingleEpoData(data_base);
						_gpre->preprocess(_site_base, begT, endT, sampling, false);
					}
				}
			}
		}
		catch (exception e)
		{
			if (_log) _log->comment(1, "t_gposproc", "preprocess failed!");
			return false;
		}

		return true;
	}

	int t_gposproc::ProcessOneEpoch(const t_gtime& nowT,vector<t_gsatdata>* data_rover, vector<t_gsatdata>* data_base)
	{
		int epoch_valid = -1;

		if (_isFLT) {//_pvtflt
			epoch_valid = _pvtflt->ProcessOneEpoch(nowT, data_rover, data_base);
		}
		else {//_ppplsq
			epoch_valid = _ppplsq->ProcessOneEpoch(nowT, data_rover);
		}
		return epoch_valid;
	}

	bool t_gposproc::_ExtractResfile()
	{
		bool res_valid = false;
		if (_isFLT) {//_pvtflt
		}
		else {//_ppplsq
			res_valid = _ppplsq->ExtractResfile();
		}
		return res_valid;
	}

	bool t_gposproc::avaliable(t_gtime nowT, string* site, string* site_base)
	{
		auto gobs = dynamic_cast<t_gallobs*>((*_alldata)[t_gdata::ALLOBS]);
		auto gorb = dynamic_cast<t_gallnav*>((*_alldata)[t_gdata::GRP_EPHEM]);
		auto gupd = dynamic_cast<t_gupd*>((*_alldata)[t_gdata::UPD]);
		auto nppmodel = dynamic_cast<t_gsetproc*>(_set)->npp_model();
		bool isClient = dynamic_cast<t_gsetgen*>(_set)->isClient();

		bool b = nowT < gobs->end_obs(site? (*site):_site);
		if (_isBase && b) {
			if (!isClient) {
				b = nowT < gobs->end_obs(site_base ? (*site_base) : _site_base);
			}
			else if (site_base) {
				b = nowT < gobs->end_obs(site_base ? (*site_base) : _site_base);
			}
		}

		if (b && !isClient && nppmodel == NPP_MODEL::PPP_RTK)
		{
			if(!gupd->upd_avail(nowT)) b = false;
		}

		return b;
	}

	t_gallpar t_gposproc::GetParam(const t_gtime& nowT,set<string>* satref)
	{
		t_gallpar param;
		if (_isFLT){
			param = _pvtflt->getUCparam(satref);//return UC param
		}
		else {
			param = _ppplsq->getAllPar();
		}
		return param;
	}

	bool t_gposproc::SppProcOneEpoch(const t_gtime& nowT, const string& site, t_gtriple* crd, double* reclk, map<string, t_gsatdata>* satdata)
	{
		auto gobs = dynamic_cast<t_gallobs*>((*_alldata)[t_gdata::ALLOBS]);
		auto gorb = dynamic_cast<t_gallnav*>((*_alldata)[t_gdata::GRP_EPHEM]);
		auto gobj = dynamic_cast<t_gallobj*>((*_alldata)[t_gdata::ALLOBJ]);
		auto gbias = dynamic_cast<t_gallbias*>((*_alldata)[t_gdata::ALLBIAS]);
		shared_ptr<t_gspplsq> spplsq = make_shared<t_gspplsq>(site, _set, gobj);
		spplsq->setDAT(gobs, gorb);
		spplsq->setDCB(gbias);
		if (spplsq->processBatch(nowT, nowT) <= 0) {
			if (_log) {
				_log->comment(0, site + "spplsq failed");
			}
			return false;
		}
		if(crd) *crd = spplsq->getCrd(nowT);
		if(reclk)*reclk = spplsq->getRecClk(nowT);
		if (satdata) spplsq->outcomsat(*satdata);
		return true;
	}

	bool t_gposproc::ProcessBatch(const t_gtime& beg, const t_gtime& end)
	{
		if (!_isFLT) _ppplsq->ProcessBatch(_alldata, beg, end);
		else {
			_pvtflt->processBatch(beg,end,true);
		}
		return true;
	}

	void t_gposproc::setWLParams(t_gallpar& params)
	{
		if (!_isFLT)
		{
			_ppplsq->setWLParams(params);
		}
	}

	void t_gposproc::setClient(bool isClient)
	{
		if (!_isFLT)
		{
			_ppplsq->setClient(isClient);
		}
	}

	void t_gposproc::setObsMode(OBSCOMBIN mode)
	{
		if (!_isFLT)
		{
			_ppplsq->setObsMode(mode);
		}
	}

	void t_gposproc::resetLsqUpdate(t_gsetbase* set, t_gallproc* data, t_glog* log)
	{
		if (!_isFLT)
		{
			_ppplsq->resetLsqUpdate(set, data, log);
		}
	}

}