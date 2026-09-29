
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  This file is part of the G-Nut C++ library.
 
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 3 of
  the License, or (at your option) any later version.
 
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, see <http://www.gnu.org/licenses>.

-*/

#include <iostream>
#include <sstream>
#include <chrono>
#include "gio/gtcp.h"
#include "gproc/gpppflt.h"
#include "gproc/gpreproc.h"
#include "gmodels/gpppmodel.h"
#include "gutils/gmatrixconv.h"
#include "gutils/gtimesync.h"
#include "gcoders/bnccorr.h"
#include "gcoders/bncobs.h"

using namespace std;

namespace gnut {

// constructor
// ----------
t_gpppflt::t_gpppflt(string mark, t_gsetbase* gset)
: t_gspp(mark, gset),
  t_gppp(mark, gset),
  t_gsppflt(mark, gset),
  t_gxml("kml"),
   _read(false),
   _flt(0),
   _enufile(0),
	_smt(0),
   _kml(false),
   _restarted(false),
   _beg_end(true)
{
  // Nothing but cout a string here
  gtrace("t_gpppflt::t_gpppflt");   

  t_gppp::_get_settings();

  this->_setOut();       // set output && identify if smoothing requested

  _epoch.tsys(t_gtime::GPS);

  _grdStoModel = new t_randomwalk();
  _ambStoModel = new t_randomwalk();
  _grdStoModel->setq(dynamic_cast<t_gsetflt*>(_set)->rndwk_grd());
  _ambStoModel->setq(dynamic_cast<t_gsetflt*>(_set)->rndwk_amb());

  int npar = _param.parNumber();
  
  // Add tropospheric gradient parameters
  if( _tropo_grad ){
    t_gpar par_grdN(_site, t_gpar::GRD_N, ++npar, "");
    t_gpar par_grdE(_site, t_gpar::GRD_E, ++npar, "");     
    par_grdN.setMF( _grd_mf );
    par_grdN.setMF( _ztd_mf );
    par_grdE.setMF( _grd_mf );
    par_grdE.setMF( _ztd_mf );     
    _param.addParam(par_grdN);
    _param.addParam(par_grdE);     
  }
  
  npar = _param.parNumber();
  
  //// Receiver P1P2 code biases parameters
  //if(_observ == RAW_DOUBLE || _observ == RAW_ALL){
  //  set<string> systems = dynamic_cast<t_gsetgen*>(_set)->sys();
  //  for(auto it = systems.begin(); it != systems.end(); it++){
  //    GSYS gsys = t_gsys::str2gsys(*it);
  //    // Add GPS P1P2 REC DCB
  //    if(gsys == GPS){
  //      t_gpar par_dcb_G(_site, t_gpar::P1P2G_REC, ++npar, "");
  //      _param.addParam(par_dcb_G);
  //    }
  //    // Add GAL P1P2 REC DCB
  //    if(gsys == GAL){
  //      t_gpar par_dcb_E(_site, t_gpar::P1P2E_REC, ++npar, "");
  //      _param.addParam(par_dcb_E);
  //    }      
  //  }
  //}          
  
  // Filling init parameter covariance matrix
  SymmetricMatrix tmp = _Qx; 
  _Qx.ReSize( _param.parNumber() ); _Qx = 0.0;
  
  for(int i = 1; i <= tmp.Nrows(); i++) _Qx(i,i) = tmp(i,i);
  
  for(size_t i=1; i<=_param.parNumber(); i++) {
    if (_param[i-1].parType == t_gpar::GRD_N)
      _Qx(i,i) = _sig_init_grd * _sig_init_grd;
    else if (_param[i-1].parType == t_gpar::GRD_E)
      _Qx(i,i) = _sig_init_grd * _sig_init_grd;
    else if (_param[i-1].parType == t_gpar::P1P2G_REC)
      _Qx(i,i) = 100 * 100;
    else if (_param[i-1].parType == t_gpar::P1P2E_REC)
      _Qx(i,i) = 100 * 100;    
  }
  
  _reset_amb = dynamic_cast<t_gsetflt*>(_set)->reset_amb();
  _reset_par = dynamic_cast<t_gsetflt*>(_set)->reset_par();   
  
} // end constructor


// destructor
// ----------
t_gpppflt::~t_gpppflt()
{ 
  gtrace("t_gpppflt::~t_gpppflt");   
   
  if( _grdStoModel ) delete _grdStoModel;
  
  if( _flt ){ if( _flt->is_open() ){ _flt->close(); }; delete _flt; }
  //LX added the enufile
  if (_enufile) { if (_enufile->is_open()) { _enufile->close(); }; delete _enufile; }

  map<OFMT, t_gtcp*>::const_iterator it = _maptcp.begin();
  while (it != _maptcp.end())
  {
	  delete it->second; it++;
  }

  this->write(_kml_name);
}

void t_gpppflt::add_coder(const vector<t_gcoder*>& coder)
{
	_gcoder = coder;
}

// Process observation batch
// ----------
int t_gpppflt::processBatch(const t_gtime& beg_r, const t_gtime& end_r, bool prtOut)
{

  gtrace("t_gpppflt::processBatch");
        
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif 
  _gmutex.lock();

  // check for site, no object means no station to be processed
  if(_grec == nullptr){
    ostringstream os;
    os << "ERROR: No object found (" << _site << "). Processing terminated!!! " << beg_r.str_ymdhms() << " -> " << end_r.str_ymdhms() << endl;
    if(_log) _log->comment(0, "t_gpppflt", os.str());
    return -1;
  }
  
   int sign = 1;
   t_gtime begT;
   t_gtime endT;
   
   // check the processing direction
   if(!_beg_end){	//from begin to end
      begT = end_r;
      endT = beg_r;
      sign = -1;      
      if(_log) _log->comment(2, "t_gpppflt", "Filtering in end -> begin direction!");
   }else{			// from end to begin
      begT = beg_r;
      endT = end_r;
      if(_log) _log->comment(2, "t_gpppflt", "Filtering in begin -> end direction!");
   }

   // set the observation, navigation and station data
  _dop.set_data(_gnav, _gobs, _site);   
   // check whether running or not
  _running = true;
   // Confirm the parameters to be estimated and set site name for all pars in gallpar
  _param.setSite(_site);
  // intval scale factor and sampling
  double subint = 0.1;
  if (_scale    > 0) subint = 1.0 / _scale;
  if (_sampling > 1) subint = pow(10, floor(log10(_sampling)));
  // Get the settings for filter
  string smtModStr(dynamic_cast<t_gsetflt*>(_set)->method_smt() );
  string fltModStr(dynamic_cast<t_gsetflt*>(_set)->method_flt() );
  int smt_delay = dynamic_cast<t_gsetflt*>(_set)->smt_delay();      
  
  // Init current time
  t_gtime now( begT );

#ifdef DEBUG
   if( _log ) _log->comment(1,"t_gpppflt", begT.str_ymdhms( "Filtering: " ) + endT.str_ymdhms(" "));
#endif

   // Cout the information for beginning
   cout << _site << ": Start filtering: " << _site << " " << now.str_ymdhms() << " " << endT.str_ymdhms() << endl;
   // Peugeot of Continuous Cycle
   bool time_loop = true;
//   high_resolution_clock::time_point chrono_t1 = high_resolution_clock::now();
   while(time_loop)     
   {
	  // Check whether the cycle is valid or not
	  // case : begin to end
      if(_beg_end && (now < endT || now == endT)) time_loop = true;
      else if(_beg_end && now > endT) {time_loop = false; break; }
      // case : end to begin     
      if(!_beg_end && (now > endT || now == endT)) time_loop = true;
      else if(!_beg_end && now < endT) {time_loop = false; break; }
      
      // synchronization
     if(now != endT){     
        if (!time_sync(now, _sampling, _scale, _log)){ // now.sod()%(int)_sampling != 0 ){      
           now.add_dsec(sign * subint / 100); // add_dsec used for synchronization!

           continue;
        }
        if(_sampling > 1) now.reset_dsec();
        
     }      
	 // added by zhShen
	 if (dynamic_cast<t_gsetproc*>(_set)->realtime())
	 {
		 while (!((t_bncobs*)_gcoder[0])->available(now) || !((t_bnccorr*)_gcoder[1])->available(now))
			 t_gtime::gmsleep(1);

		 t_gpreproc prep(_gobs, _set);
		 prep.glog(_log);
		 t_gtime epo_prev = now - 10.0;
		 if (_sampling > 1) epo_prev.add_secs(-int(sign * 2 * _sampling));  // =<1Hz data
		 else               epo_prev.add_dsec(-sign * 2 * _sampling);       //  >1Hz data     
		 prep.preprocess(_site, epo_prev, now, _sampling * 2, false);
	 }
     // clean/collect/filter epoch data
     _data.erase( _data.begin(), _data.end() );
	 // Get obs data for new epoch 
     _data = _gobs->obs( _site, now );                
	 // Check whether the data for current epoch is empty or not 
     if( _data.size() == 0 ){
       if( _log && _log->verb() >= 1 ){
           _log->comment(1,"gpppflt", _site + now.str_ymdhms(" no observation found at epoch: "));
       }
       if (_sampling > 1) now.add_secs(int(sign*_sampling));  // =<1Hz data
       else               now.add_dsec(sign*_sampling);       //  >1Hz data
       continue;
     }
	 // Update time
      t_gtime obsEpo = _data.begin()->epoch();
     _timeUpdate(obsEpo); 

     // reseting params according to conf
     if (_reset_amb > 0) {
       if (now.sod() % _reset_amb == 0) {
         _reset_ambig();
       }
     }
     if (_reset_par > 0) {
       if (now.sod() % _reset_par == 0) {
         _reset_param();
       }
     }

     // save apriory coordinates     
     _saveApr(obsEpo, _param, _Qx);

#ifdef _WIN32
	 system("cls");
#else
	 system("clear");
#endif
	 cerr << now.str() << endl;
	// t_gtime::gmsleep(1000*1.5);
     int irc_epo = _processEpoch(obsEpo);      

     if( irc_epo < 0 ) {
       _success = false;
       _removeApr(obsEpo);
       if( _log && _log->verb() >= 1 ) 
       _log->comment(1,"gpppflt", _site + now.str_ymdhms(" epoch ") + " was not calculated" );
       
       if (_sampling > 1) now.add_secs(int(sign*_sampling));  // =<1Hz data
       else               now.add_dsec(sign*_sampling);       //  >1Hz data
       
       // repeat preprocessing in case of skipped epoch -> sampling changed
       t_gpreproc prep(_gobs, _set);
       prep.glog(_log);
       t_gtime epo_prev = now;
       if (_sampling > 1) epo_prev.add_secs(-int(sign*2*_sampling));  // =<1Hz data
       else               epo_prev.add_dsec(-sign*2*_sampling);       //  >1Hz data     
       prep.preprocess(_site, epo_prev, now, _sampling*2, false);
       //       cout << "New preprocessing due to skipping epoch: " << epo_prev.str_hms() << " " << now.str_hms() << " " << _sampling*2 << endl;

       continue;
     }else _success = true;

     if( _log && _log->verb() >= 1 )
     _log->comment(1,"gpppflt", _site + now.str_ymdhms(" processing epoch: "));
     
     ostringstream os;
     bool savePrd = true;
     if(_smooth) savePrd = false;
     _prtOut(obsEpo, _param, _Qx, _data, os, line, savePrd);

     // Print flt results
     if( _flt && prtOut){
       _flt->write(os.str().c_str(), os.str().size());
       _flt->flush();
     }
     
     _initialized = true;
     
     if (_sampling > 1) now.add_secs(int(sign*_sampling));  // =<1Hz data
     else               now.add_dsec(sign*_sampling);       //  >1Hz data
     
     // RTS Smooth
     if ( _smooth ){
       if (smt_delay != 0 && now.sod() % smt_delay == 0 ){
         _backflt(prtOut);
         _vfltdata.clear();
       }      
     }      
     
   }//end now

  _running = false;  
//   high_resolution_clock::time_point chrono_t2 = high_resolution_clock::now();
//   auto duration = duration_cast<microseconds>( chrono_t2 - chrono_t1 ).count();
//   cout << "Duration FLT: " << duration << endl;
   
//   chrono_t1 = high_resolution_clock::now();
  if (_smooth && smt_delay == 0){
     _backflt(prtOut);
  }
//   chrono_t2 = high_resolution_clock::now();
//   duration = duration_cast<microseconds>( chrono_t2 - chrono_t1 ).count();
//   cout << "Duration SMT: " << duration << endl;

   if(beg_r != end_r){   // do not print in real-time
      double npd_perc = 0;
      npd_perc = (double(_n_NPD_flt) / double(_n_ALL_flt)) * 100;
      if( _log ) _log->comment(1,"gpppflt", _site + ": Forward filter -     Not positive definite covariance matrices: " + dbl2str(npd_perc,0) + " %   ("
                               + int2str(_n_NPD_flt) + ", "  + int2str(_n_ALL_flt) + ")");
     
      if( _smooth ){
        npd_perc = (double(_n_NPD_smt) / double(_n_ALL_smt)) * 100;
        if( _log ) _log->comment(1,"gpppflt", _site + ": Backward smoothing - Not positive definite covariance matrices: " + dbl2str(npd_perc,0) + " %   ("
                               + int2str(_n_NPD_smt) + ", "  + int2str(_n_ALL_smt) + ")");
      }
   }

  _gmutex.unlock();

  return 1;
}

// Process observation batch in forward <-> backward direction
// ----------
int t_gpppflt::processBatchFB(const t_gtime& beg_r, const t_gtime& end_r)
{
   gtrace("t_gpppflt::processBatchFB");
        
   this->processBatch(beg_r, end_r, false);
     
   _beg_end = false;
     
   this->processBatch(beg_r, end_r, true);
     
   _beg_end = true;     
   
  return 1;
}
   
   
// ----------------------------------------------------------------------------------
// PROTECTED
// ----------------------------------------------------------------------------------
// 

// Update time for RDW processes
// ---------------------------
void t_gpppflt::_timeUpdate(const t_gtime& epo)
{
   t_gsppflt::_timeUpdate(epo);
   
   //tropo gradient
   _grdStoModel->updateTime(epo);
   _grdStoModel->setTcurr(epo);

   // amb rndwk
   _ambStoModel->updateTime(epo);
   _ambStoModel->setTcurr(epo);

}                   


// Process one epoch
// -----------------------------
int t_gpppflt::_processEpoch(const t_gtime& runEpoch)
{
   gtrace("t_gpppflt::_processEpoch");   
   
   int irc = t_gsppflt::_processEpoch(runEpoch);
/*   
   // restart ambiguities in case of many cycle slips 
   int newamb = 0;
   for(map<string, int>::iterator it = _newAMB.begin(); it != _newAMB.end(); it++) {      
        if(it->second == 2 ) newamb++;
   }

   if(_data.size() - newamb <= 4 && _initialized && !_restarted && irc > 0){
      cout << "Restart " << runEpoch.str_hms() << " " << _data.size() << " " << _param.ambNumber() << " " << newamb << " " << irc << " " << _initialized << endl;
      _reset_ambig(); 
      _restarted = true;
      irc = -1;
   }
   
   _restarted = false;
*/ 
   
   return irc;
}


// add one satellite data to A, l, P - carrier phase
// ---------------------------
int t_gpppflt::_addObsL(t_gsatdata& satdata, unsigned int& iobs, t_gtriple& ell, Matrix& A, ColumnVector& l, DiagonalMatrix& P)
{
  t_gsys gsys(satdata.gsys());
  GSYS gs = gsys.gsys();

  GOBSBAND b1, b2, b3, b4, b5;
  b1 = b2 = b3 = b4 = b5 = BAND;
  
  if(_auto_band){    // automatic dual band selection -> for Anubis purpose
    set<GOBSBAND> bands = satdata.band_avail();
    auto itBAND = bands.begin();
    if(bands.size() < 2) return -1;
    b1 = *itBAND;
    itBAND++;
    b2 = *itBAND;
  }else {           // select fix defined band according the table
    b1 = t_gsys::band_priority(gs, FREQ_1);
    b2 = t_gsys::band_priority(gs, FREQ_2);
    b3 = t_gsys::band_priority(gs, FREQ_3);
    b4 = t_gsys::band_priority(gs, FREQ_4);
    b5 = t_gsys::band_priority(gs, FREQ_5);
  }
  
  t_gobs gobs1, gobs2, gobs3, gobs4, gobs5;                       
  
  _getgobs(satdata.sat(), TYPE_L, b1, gobs1);
  _getgobs(satdata.sat(), TYPE_L, b2, gobs2);
  _getgobs(satdata.sat(), TYPE_L, b3, gobs3);
  _getgobs(satdata.sat(), TYPE_L, b4, gobs4);
  _getgobs(satdata.sat(), TYPE_L, b5, gobs5);
  
  if( _phase ) {
         
    bool com = true;
    if( _gnav ) com = _gnav->com();
       
    double LIF, L1, L2, L3, L4, L5;
    LIF = L1 = L2 = L3 = L4 = L5 = 0.0;
    
    LIF = satdata.L3(gobs1, gobs2);
    L1  = satdata.obs_L(gobs1);
    L2  = satdata.obs_L(gobs2);  
    if(b3 != BAND) L3 = satdata.obs_L(gobs3);
    if(b4 != BAND) L4 = satdata.obs_L(gobs4);
    if(b5 != BAND) L5 = satdata.obs_L(gobs5);
    
    if(_observ == IONO_FREE  &&  double_eq(LIF, 0.0)) return -1;
    if(_observ == RAW_SINGLE &&  double_eq(L1,  0.0)) return -1;
    if((_observ == RAW_DOUBLE || _observ == RAW_ALL) && (double_eq(L1,  0.0) || double_eq(L2, 0.0))) return -1;
       
    double weight_coef = _weightObs(satdata, gobs1);
          
    int addIobs = 1;
    if(L3 > 0) addIobs++;
    if(L4 > 0) addIobs++;
    if(L5 > 0) addIobs++;         
    
    double sigPhase = 0.0;
    switch( gs ){
           case GPS : sigPhase = _sigPhaseGPS; break;
           case GLO : sigPhase = _sigPhaseGLO; break;
           case GAL : sigPhase = _sigPhaseGAL; break;
           case BDS : sigPhase = _sigPhaseBDS; break;
           case QZS : sigPhase = _sigPhaseQZS; break;
           default  : sigPhase = 0.0;
    }     
    
          // Create weight matrix
    P(iobs,  iobs  ) = weight_coef * (1/(sigPhase * sigPhase));
    if(_observ == RAW_DOUBLE ||
                 _observ == RAW_ALL)      P(iobs+1,iobs+1) = weight_coef * (1/(sigPhase * sigPhase));
    if(_observ == RAW_ALL)
    for(int i = 2; i <= addIobs; i++) P(iobs+i,iobs+i) = weight_coef * (1/(sigPhase * sigPhase));
        
    GOBS type1, type2, type3, type4, type5;
    
    if(gobs1.gattr().attr() == ATTR) type1 = satdata.id_phase(gobs1.gband().band());
    else type1 = gobs1.gobs();
    if(gobs2.gattr().attr() == ATTR) type2 = satdata.id_phase(gobs2.gband().band());
    else type2 = gobs2.gobs();
    if(gobs3.gattr().attr() == ATTR) type3 = satdata.id_phase(gobs3.gband().band());
    else type3 = gobs3.gobs();
    if(gobs4.gattr().attr() == ATTR) type4 = satdata.id_phase(gobs4.gband().band());
    else type4 = gobs4.gobs();
    if(gobs5.gattr().attr() == ATTR) type5 = satdata.id_phase(gobs5.gband().band());
    else type5 = gobs5.gobs();
    
    for( size_t iPar=0; iPar<_param.parNumber(); iPar++ ){
      if(_param[iPar].prn.compare( satdata.sat() ) == 0 ) {       
                 
        double modObsL = 0.0;                
        
        double obs_L = 0.0;
        if(     _param[iPar].parType == t_gpar::AMB_IF) {obs_L = LIF; modObsL = _gModel->cmpObs(_epoch, _param, satdata, gobs1, com);}
        else if(_param[iPar].parType == t_gpar::AMB_L1) {obs_L = L1;  modObsL = _gModel->cmpObs(_epoch, _param, satdata, gobs1, com);}
        else if(_param[iPar].parType == t_gpar::AMB_L2) {obs_L = L2;  modObsL = _gModel->cmpObs(_epoch, _param, satdata, gobs2, com);}
        else if(_param[iPar].parType == t_gpar::AMB_L3) {obs_L = L3;  modObsL = _gModel->cmpObs(_epoch, _param, satdata, gobs3, com);}
        else if(_param[iPar].parType == t_gpar::AMB_L4) {obs_L = L4;  modObsL = _gModel->cmpObs(_epoch, _param, satdata, gobs4, com);}
        else if(_param[iPar].parType == t_gpar::AMB_L5) {obs_L = L5;  modObsL = _gModel->cmpObs(_epoch, _param, satdata, gobs5, com);}
   
        if (modObsL < 0) return -1;

        // Reset ambiguity in case of cycle slip
        satdata.addslip(false);
        if( (_param[iPar].parType == t_gpar::AMB_L1 && satdata.getlli(type1) >= 1) ||
            (_param[iPar].parType == t_gpar::AMB_L2 && satdata.getlli(type2) >= 1) || 
            (_param[iPar].parType == t_gpar::AMB_L3 && satdata.getlli(type3) >= 1) ||
            (_param[iPar].parType == t_gpar::AMB_L4 && satdata.getlli(type4) >= 1) ||
            (_param[iPar].parType == t_gpar::AMB_L5 && satdata.getlli(type5) >= 1) ||
            (_param[iPar].parType == t_gpar::AMB_IF && (satdata.getlli(type1) >= 1 || satdata.getlli(type2) >= 1)) ){ 
              satdata.addslip(true);
              string prn = satdata.sat();
              _param[iPar].value( obs_L - modObsL );
			  //cout << "cycle slip !!!!" << endl;
              _Qx(iPar+1, iPar+1) = _sigAmbig * _sigAmbig;
              _newAMB[prn] = 1;
              if(_smooth) _slips.insert(prn);
//              cout << "Cycle slip: " << prn << " " << gobs2str(type1) << " " << gobs2str(type2) << " " << _epoch.str_hms() << fixed << setprecision(3) << " " << _param[iPar].value() << " " << LIF << " " << modObsL << endl;
              if (_smooth){
                SymmetricMatrix Qp = _fltdata.Qp();
                Qp(iPar+1, iPar+1) = _sigAmbig * _sigAmbig;
                _fltdata.Qp(Qp);
              }           
            }

        // Create reduced measurement (prefit residuals)                        
        if(_param[iPar].parType == t_gpar::AMB_L1 || _param[iPar].parType == t_gpar::AMB_IF){
          l(iobs) = obs_L - _param[iPar].value() - modObsL;
		/*  cout << "obs_L modObsL amb" << setprecision(5)
			  << setw(20) << obs_L
			  << setw(20) << modObsL
			  << setw(20) << _param[iPar].value() << endl;*/


        }
        if((_observ == RAW_DOUBLE || _observ == RAW_ALL) && _param[iPar].parType == t_gpar::AMB_L2){
          l(iobs+1) = obs_L - _param[iPar].value() - modObsL;
        }
        if(_observ == RAW_ALL){
          if(_param[iPar].parType == t_gpar::AMB_L3) l(iobs+2) = obs_L - _param[iPar].value() - modObsL;
          if(_param[iPar].parType == t_gpar::AMB_L4) l(iobs+3) = obs_L - _param[iPar].value() - modObsL;
          if(_param[iPar].parType == t_gpar::AMB_L5) l(iobs+4) = obs_L - _param[iPar].value() - modObsL;                           
        }

        // check satellite post shadow period
        if(satdata.ecl()) return -1;
      }
                                 
    }      
          
    // Create first design matrix 
    for( unsigned int ipar=1; ipar<=_param.parNumber(); ipar++ ){
      A(iobs,  ipar) = _param[ipar-1].partial( satdata, _epoch, ell, gobs1);
      if(_observ == RAW_DOUBLE ||
         _observ == RAW_ALL)      A(iobs+1,ipar) = _param[ipar-1].partial( satdata, _epoch, ell, gobs2);
      if(_observ == RAW_ALL){
        for(int i = 2; i <= addIobs; i++) {                       
          if(i == 2) A(iobs+i,ipar) = _param[ipar-1].partial( satdata, _epoch, ell, gobs3);
          if(i == 3) A(iobs+i,ipar) = _param[ipar-1].partial( satdata, _epoch, ell, gobs4);
          if(i == 4) A(iobs+i,ipar) = _param[ipar-1].partial( satdata, _epoch, ell, gobs5);
        }
      }
    }   
    
    if(_observ == RAW_DOUBLE) iobs += 2;
    else if(_observ == RAW_ALL) iobs = iobs + addIobs + 1;
    else iobs++;          
    
  } 
  
  return 1;
}
   
// Applying tides
// -------------------------------------------------------
int t_gpppflt::_apply_tides(t_gtime& _epoch, t_gtriple& xRec)
{
   gtrace("t_gpppflt::_apply_tides");

    if( _tides ){
      xRec = xRec + _tides->tide_searth(_epoch, xRec)
                  + _tides->load_ocean(_epoch, _site, xRec)
                  + _tides->tide_pole()
                  + _tides->load_atmosph();
    } else return -1;
//cout << fixed << setprecision(10) << "xRec = " << xRec[0] << " " << xRec[1] << " " << xRec[2] << endl;   
   return 1;
}

// Predict state vector and covariance matrix
// =======================================================
void t_gpppflt::_predict()
{
   gtrace("t_gpppflt::_predict");

  // Predict coordinates, clock and troposphere
  // state vector is the same, covariance matrix is predicted with noise   
   int i;      
   
   t_gsppflt::_predict();   
   
   i = _param.getParam(_site, t_gpar::GRD_N, "");
   if ( i >= 0 )  {
     if (_cntrep == 1 && _initialized) _Qx(i+1,i+1) += _grdStoModel->getQ();
     if(_smooth) _Noise(i+1,i+1)  = _grdStoModel->getQ();
   }

   i = _param.getParam(_site, t_gpar::GRD_E, "");
   if ( i >= 0 )  {
     if (_cntrep == 1 && _initialized) _Qx(i+1,i+1) += _grdStoModel->getQ();
     if(_smooth) _Noise(i+1,i+1)  = _grdStoModel->getQ();
   }  
   
  // Predict ambiguity 
  // Add or remove ambiguity parameter and appropriate rows/columns covar. matrix
  if( _phase )  _syncAmb();

//   cout << "Predicted Q: \n" << _Qx << endl;
//   cout << "Predicted X: \n" << _param << endl;
   
  // Create _Noise matrix (DiagonalMatrix)
  if( _smooth ){
    Matrix tmp = _Noise;
    _Noise.ReSize(_Qx.Nrows()); _Noise = 0;     
    for (int i = 1; i <= _Noise.Ncols(); i++){
      if (i > tmp.Ncols()) break;
      _Noise(i,i) = tmp(i,i);
    }
  }
  
} // end prediction


// Satelite position
// ----------
int t_gpppflt::_satPos(t_gtime& epo, t_gsatdata& gsatdata)
{
  gtrace("t_gpppflt::_satPos");      

  if( t_gsppflt::_satPos(epo, gsatdata) < 0 ) return -1;

  /* PCO correction moved to t_gpppmodel

  string satname  =  gsatdata.sat();
  
  // phase center offset correction
  // phase center variations are applied in gpppmodel.cpp
  if( _gallobj != 0 ){
    t_gtriple crd           =  gsatdata.satcrd();

    shared_ptr<t_gobj> obj  =  _gallobj->obj( satname );
    
    shared_ptr<t_gpcv> apc;
    if ( obj != 0 ) apc =  obj->pcv( epo );
    
    t_gtriple pco(0,0,0);
    string antenna = "";
    
    if ( apc != 0 ) {                
      apc->pcoS(gsatdata, pco, LC_IF);
      antenna = apc->anten();
    }else{
      t_gpcv apc_def;
      apc_def.pcoS(gsatdata, pco, LC_IF);
    }
    
    // transformation to ECEF sytem
    t_gtriple dx(0,0,0);
    if( _gnav && _gnav->com() ){
      ColumnVector i(3);
      ColumnVector j(3);
      ColumnVector k(3);
      
      dynamic_cast<t_gpppmodel*>(_gModel)->attitude(gsatdata, antenna, i, j, k);
      dx[0] = pco[0]*i(1) + pco[1]*j(1) + pco[2]*k(1);
      dx[1] = pco[0]*i(2) + pco[1]*j(2) + pco[2]*k(2);
      dx[2] = pco[0]*i(3) + pco[1]*j(3) + pco[2]*k(3);
      
      for (int i=0; i<3; i++) crd[i] += dx[i];
    }                
    gsatdata.addcrd(crd);
    gsatdata.addpco(dx);
    
  }
  */
  
  return 1;
}


// Synchronize ambiguity with respect to observation data.
// First new ambiguity parameters are added.
// Second old ambiguity parameters are deleted;
// -----------------------------------------------------
void t_gpppflt::_syncAmb()
{
  gtrace("t_gpppflt::_syncAmb");

  for(map<string, int>::iterator it = _newAMB.begin(); it != _newAMB.end(); ){
     if(it->second > 2) _newAMB.erase(it++);
     else ++it;
  }

  // Add ambiguity parameter and appropriate rows/columns covar. matrix   
  set<string> mapPRN;

  vector<t_gsatdata>::iterator it;
  for( it = _data.begin(); it != _data.end(); it++){       // loop over all satellites
    mapPRN.insert( it->sat() );

    GSYS gs = it->gsys();

    GOBSBAND b1, b2, b3, b4, b5;
    b1 = b2 = b3 = b4 = b5 = BAND;
    
    if(_auto_band){    // automatic dual band selection -> for Anubis purpose
      set<GOBSBAND> bands = it->band_avail();
      auto itBAND = bands.begin();
      if(bands.size() < 2) continue;
      b1 = *itBAND;
      itBAND++;
      b2 = *itBAND;
    }else{            // select fix defined band according the table    
      b1 = t_gsys::band_priority(gs, FREQ_1);
      b2 = t_gsys::band_priority(gs, FREQ_2);
      b3 = t_gsys::band_priority(gs, FREQ_3);
      b4 = t_gsys::band_priority(gs, FREQ_4);
      b5 = t_gsys::band_priority(gs, FREQ_5);    
    }
     
    t_gobs gobs1, gobs2, gobs3, gobs4, gobs5;

     _getgobs(it->sat(), TYPE_L, b1, gobs1);
     _getgobs(it->sat(), TYPE_L, b2, gobs2);
     _getgobs(it->sat(), TYPE_L, b3, gobs3);
     _getgobs(it->sat(), TYPE_L, b4, gobs4);
     _getgobs(it->sat(), TYPE_L, b5, gobs5);     

     double LIF, L1, L2, L3, L4, L5;
     LIF = L1 = L2 = L3 = L4 = L5 = 0.0;

     LIF = it->L3(gobs1, gobs2);
     L1 = it->obs_L(gobs1);
     L2 = it->obs_L(gobs2);
     if(b3 != BAND) L3 = it->obs_L(gobs3);
     if(b4 != BAND) L4 = it->obs_L(gobs4);
     if(b5 != BAND) L5 = it->obs_L(gobs5);

     if(_observ == IONO_FREE  &&  double_eq(LIF, 0.0)) continue;
     if(_observ == RAW_SINGLE &&  double_eq(L1, 0.0)) continue;
     if(_observ == RAW_DOUBLE && (double_eq(L1, 0.0) || double_eq(L2, 0.0))) continue;
     
     bool com = true;
     if( _gnav ) com = _gnav->com();
         
     if(_observ == IONO_FREE && _param.getParam(_site, t_gpar::AMB_IF, it->sat()) < 0 ) {

       t_gpar newPar(it->site(), t_gpar::AMB_IF, _param.parNumber()+1, it->sat());
       double cmpObs = _gModel->cmpObs(_epoch, _param, *it, gobs1, com);
       if (cmpObs < 0) continue;
       
       newPar.value( LIF - cmpObs );      // first ambiguity value
       
       _param.addParam(newPar);
       _newAMB[it->sat()] = 1;
       
       Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
       _Qx(_param.parNumber(), _param.parNumber()) = _sigAmbig * _sigAmbig;
//cout << "ADD AMB IF " << it->sat() << " " << it->epoch().str_ymdhms() << endl;
       if( _log ) _log->comment(3, "gppp", "AMB_IF was added! For Sat PRN " + it->sat() + " Epoch: " + _epoch.str_ymdhms());           

     }else if(_observ == RAW_SINGLE || _observ == RAW_DOUBLE || _observ == RAW_ALL){                                
       int newAmb = 0;
       // add L1 amb - everytime when RAW observations
       if( _param.getParam(_site, t_gpar::AMB_L1, it->sat()) < 0 ){
         double cmpObs = _gModel->cmpObs(_epoch, _param, *it, gobs1, com);
         if (cmpObs < 0) continue;
         t_gpar newPar1(it->site(), t_gpar::AMB_L1, _param.parNumber()+1, it->sat());         
         newPar1.value( L1 - cmpObs );
         _param.addParam(newPar1);
         Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
         _Qx(_param.parNumber(), _param.parNumber()) = _sigAmbig * _sigAmbig;
         //         cout << "ADD AMB L1 " << it->sat() << " " << it->epoch().str_ymdhms()  << " " << L1 - cmpObs << endl;
         if( _log ) _log->comment(3, "gppp", "RAW AMB_L1 was added! For Sat PRN " + it->sat() + " Epoch: " + _epoch.str_ymdhms());
         newAmb = 1;
       }
                
       // add L2 amb - when RAW DOUBLE or ALL frequency
       if( _observ != RAW_SINGLE && _param.getParam(_site, t_gpar::AMB_L2, it->sat()) < 0 ){
         double cmpObs = _gModel->cmpObs(_epoch, _param, *it, gobs2, com);
         if (cmpObs < 0) continue;
         t_gpar newPar2(it->site(), t_gpar::AMB_L2, _param.parNumber()+1, it->sat());
         newPar2.value( L2 - cmpObs );
         _param.addParam(newPar2);
         Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
         _Qx(_param.parNumber(), _param.parNumber()) = _sigAmbig * _sigAmbig; 
         //         cout << "ADD AMB L2 " << it->sat() << " " << it->epoch().str_ymdhms() << " " << L2 - cmpObs << endl;        
         if( _log ) _log->comment(3, "gppp", "RAW AMB_L2 was added! For Sat PRN " + it->sat() + " Epoch: " + _epoch.str_ymdhms());
         newAmb = 1;
       }
       
       // add L3 amb - when RAW ALL frequency
       if( _observ != RAW_SINGLE && _observ != RAW_DOUBLE && _param.getParam(_site, t_gpar::AMB_L3, it->sat()) < 0 && L3 > 0.0){
         double cmpObs = _gModel->cmpObs(_epoch, _param, *it, gobs3, com);
         if (cmpObs < 0) continue;
         t_gpar newPar3(it->site(), t_gpar::AMB_L3, _param.parNumber()+1, it->sat());
         newPar3.value( L3 - cmpObs );
         _param.addParam(newPar3);
         Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
         _Qx(_param.parNumber(), _param.parNumber()) = _sigAmbig * _sigAmbig;         
         if( _log ) _log->comment(3, "gppp", "RAW AMB_L3 was added! For Sat PRN " + it->sat() + " Epoch: " + _epoch.str_ymdhms());               
         newAmb = 1;
       }
       
       // add L4 amb - when RAW ALL frequency
       if( _observ != RAW_SINGLE && _observ != RAW_DOUBLE && _param.getParam(_site, t_gpar::AMB_L4, it->sat()) < 0 && L4 > 0.0){
         double cmpObs = _gModel->cmpObs(_epoch, _param, *it, gobs4, com);
         if (cmpObs < 0) continue;
         t_gpar newPar4(it->site(), t_gpar::AMB_L4, _param.parNumber()+1, it->sat());
         newPar4.value( L4 - cmpObs );
         _param.addParam(newPar4);
         Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
         _Qx(_param.parNumber(), _param.parNumber()) = _sigAmbig * _sigAmbig;         
         if( _log ) _log->comment(3, "gppp", "RAW AMB_L4 was added! For Sat PRN " + it->sat() + " Epoch: " + _epoch.str_ymdhms());               
         newAmb = 1;
       }
      
       // add L5 amb - when RAW ALL frequency
       if( _observ != RAW_SINGLE && _observ != RAW_DOUBLE && _param.getParam(_site, t_gpar::AMB_L5, it->sat()) < 0 && L5 > 0.0){
         double cmpObs = _gModel->cmpObs(_epoch, _param, *it, gobs5, com);
         if (cmpObs < 0) continue;
         t_gpar newPar5(it->site(), t_gpar::AMB_L5, _param.parNumber()+1, it->sat());
         newPar5.value( L5 - cmpObs );
         _param.addParam(newPar5);
         Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
         _Qx(_param.parNumber(), _param.parNumber()) = _sigAmbig * _sigAmbig;         
         if( _log ) _log->comment(3, "gppp", "RAW AMB_L5 was added! For Sat PRN " + it->sat() + " Epoch: " + _epoch.str_ymdhms());               
         newAmb = 1;
       }      
       
       if(newAmb) _newAMB[it->sat()] = 1;
     }
         
           
  }  // end loop over all observations

// Remove ambiguity parameter and appropriate rows/columns covar. matrix

  for (unsigned int i=0; i<_param.parNumber(); i++) {
    if (_param[i].parType == t_gpar::AMB_IF ||
        _param[i].parType == t_gpar::AMB_L1 ||
        _param[i].parType == t_gpar::AMB_L2 ||
        _param[i].parType == t_gpar::AMB_L3 ||
        _param[i].parType == t_gpar::AMB_L4 ||        
        _param[i].parType == t_gpar::AMB_L5  ) {      
       
      set<string>::iterator prnITER = mapPRN.find(_param[i].prn);
      if (prnITER == mapPRN.end()) {
        
        if( _log )
          _log->comment(3, "gppp", "AMB will be removed! For Sat PRN " + _param[i].prn + " Epoch: " + _epoch.str_ymdhms());

#ifdef DEBUG
   cout << "AMB will be removed! For Sat PRN " << _param[i].prn
        << " Epoch: " << _epoch.str_ymdhms() << endl;
#endif        

   _newAMB.erase(_param[i].prn);
        
   Matrix_remRC(_Qx, _param[i].index, _param[i].index);
   _param.delParam(i);
   _param.reIndex();
   i--;
      }   
    }   
  }   
    
  return;
}


// Save apriory coordinates to products
// -----------------------------
void t_gpppflt::_saveApr(t_gtime& epoch, t_gallpar& X, const SymmetricMatrix& Q)
{
  gtrace("t_gpppflt::_saveApr");
   
  // get CRD params
  t_gtriple xyz;
  X.getCrdParam(_site, xyz);

  // CRD using eccentricities
  t_gtriple xyz_ecc = xyz - _grec->eccxyz(epoch); // MARKER + ECC = ARP
   
  // get CRD rms  (XYZ)
  double Xrms = 0.0, Yrms = 0.0, Zrms = 0.0;
  int icrdx = _param.getParam(_site, t_gpar::CRD_X, "");
  int icrdy = _param.getParam(_site, t_gpar::CRD_Y, "");
  int icrdz = _param.getParam(_site, t_gpar::CRD_Z, "");
  if (Q(icrdx + 1, icrdx + 1) < 0) Xrms = -1; else Xrms = sqrt(Q(icrdx + 1, icrdx + 1));
  if (Q(icrdy + 1, icrdy + 1) < 0) Yrms = -1; else Yrms = sqrt(Q(icrdy + 1, icrdy + 1));
  if (Q(icrdz + 1, icrdz + 1) < 0) Zrms = -1; else Zrms = sqrt(Q(icrdz + 1, icrdz + 1));
  t_gtriple crd_rms(Xrms, Yrms, Zrms);
  
  if( _allprod != 0 ){
    shared_ptr<t_gprod> prdcrd = _allprod->get(_site, t_gdata::POS, epoch);                      
    if( prdcrd ) {       
      dynamic_pointer_cast<t_gprodcrd>(prdcrd)->apr(xyz_ecc);
      dynamic_pointer_cast<t_gprodcrd>(prdcrd)->apr_rms(crd_rms);
    }else{
      prdcrd = make_shared<t_gprodcrd>(epoch);
      dynamic_pointer_cast<t_gprodcrd>(prdcrd)->apr(xyz_ecc);
      dynamic_pointer_cast<t_gprodcrd>(prdcrd)->apr_rms(crd_rms);
      _allprod->add(prdcrd, _site);        
    }            
  }
}
   
// Remove apriory coordinates from products
// -----------------------------
void t_gpppflt::_removeApr(t_gtime& epoch)
{
  gtrace("t_gpppflt::_removeApr");   
   
  if( _allprod != 0 ){
    _allprod->rem(_site, t_gdata::POS, epoch);                      
  }
}
   
// Set output & return if smooth
// -----------------------------
void t_gpppflt::_prtOut(t_gtime& epoch, t_gallpar& X, const SymmetricMatrix& Q, vector<t_gsatdata>& data, ostringstream& os, xml_node& node, bool saveProd)
{
  gtrace("t_gpppflt::_prtOut");
    
  // get CRD params
  t_gtriple xyz, ell;
  X.getCrdParam(_site, xyz);
  xyz2ell(xyz, ell, false);

  // CRD using eccentricities
  t_gtriple xyz_ecc = xyz - _grec->eccxyz(epoch); // MARKER + ECC = ARP

  // get CRD rms  (XYZ)
  double Xrms = 0.0, Yrms = 0.0, Zrms = 0.0;
  double cov_xy = 0.0, cov_xz = 0.0, cov_yz = 0.0;
  if( Q(1,1) < 0 ) Xrms = -1; else Xrms = sqrt(Q(1,1));
  if( Q(2,2) < 0 ) Yrms = -1; else Yrms = sqrt(Q(2,2));      
  if( Q(3,3) < 0 ) Zrms = -1; else Zrms = sqrt(Q(3,3));
  t_gtriple crd_rms(Xrms, Yrms, Zrms);
  cov_xy = Q(2,1);
  cov_xz = Q(3,1);
  cov_yz = Q(3,2);   
   
  // get CLK param
  double clk = 0.0;
  int iclk = X.getParam(_site, t_gpar::CLK, "");
  if( iclk >= 0 ) clk = X[iclk].value();
  
  // get ZTD params
  double ZHD = 0.0, ZWD = 0.0, ZTD = 0.0, ZWDrms = 0.0;
  int iztd = X.getParam(_site, t_gpar::TRP, "");
  if( iztd >= 0 ){
    ZHD = X[iztd].apriori();//_gModel->tropoModel()->getZHD(ell, epoch);
    ZWD = X[iztd].value();
    ZTD = ZHD + ZWD;
    if( Q(iztd+1,iztd+1) < 0 ) ZWDrms = -1;
    else ZWDrms = sqrt(Q(iztd+1,iztd+1));          // _param is indexed from 0, but Q from 1 !
  }else{
     // test purpose only (if ZTD not estimated)
    if(_gModel->tropoModel() != 0){
      ZWD = _gModel->tropoModel()->getZWD(ell, epoch);
      ZHD = _gModel->tropoModel()->getZHD(ell, epoch);
      ZTD = ZHD + ZWD;
    }          
  }

  // get Tropo gradient
  int ign = X.getParam(_site, t_gpar::GRD_N, "");
  int ige = X.getParam(_site, t_gpar::GRD_E, "");
  double grdN = 0.0, grdE = 0.0;
//  double rmsN = 0.0, rmsE = 0.0;
  if( ign >= 0 && ige >= 0 ){     
    grdN = X[ign].value();
    grdE = X[ige].value();
//    rmsN = sqrt(Q(ign+1,ign+1));                  // _param is indexed from 0, but Q from 1 !
//    rmsE = sqrt(Q(ige+1,ige+1));                  // _param is indexed from 0, but Q from 1 !
  }

  set<string> sat_list;
  for(auto it = _data.begin(); it != _data.end(); it++) sat_list.insert(it->sat());
  _dop.set_sats(sat_list);
  double gdop = -1;
  if (_dop.calculate(epoch, xyz) >=0 ) gdop = _dop.gdop();  
      
  set<string> ambs = X.amb_prns();
  int nsat = ambs.size();   
   
  if( _allprod != 0 && saveProd ){
    shared_ptr<t_gprod> prdcrd = _allprod->get(_site, t_gdata::POS, epoch);                      
    if( prdcrd ) {
      dynamic_pointer_cast<t_gprodcrd>(prdcrd)->xyz(xyz_ecc);
      dynamic_pointer_cast<t_gprodcrd>(prdcrd)->xyz_rms(crd_rms);
      dynamic_pointer_cast<t_gprodcrd>(prdcrd)->cov(COV_XY, cov_xy);
      dynamic_pointer_cast<t_gprodcrd>(prdcrd)->cov(COV_XZ, cov_xz);
      dynamic_pointer_cast<t_gprodcrd>(prdcrd)->cov(COV_YZ, cov_yz);
    }else{
      prdcrd = make_shared<t_gprodcrd>(epoch);
      dynamic_pointer_cast<t_gprodcrd>(prdcrd)->xyz(xyz_ecc);
      dynamic_pointer_cast<t_gprodcrd>(prdcrd)->xyz_rms(crd_rms);
      dynamic_pointer_cast<t_gprodcrd>(prdcrd)->cov(COV_XY, cov_xy);
      dynamic_pointer_cast<t_gprodcrd>(prdcrd)->cov(COV_XZ, cov_xz);
      dynamic_pointer_cast<t_gprodcrd>(prdcrd)->cov(COV_YZ, cov_yz);   
      _allprod->add(prdcrd, _site);              
    }

    t_gtriple Ell, XYZ;   
    if (_param.getCrdParam(_site, XYZ) > 0){
    }else if ( _valid_crd_xml ) { XYZ = _grec->crd_arp(_epoch); }
    xyz2ell( XYZ, Ell, false);
    int itrp = _param.getParam(_site, t_gpar::TRP, "");
    if(itrp >= 0) _param[itrp].apriori(_gModel->tropoModel()->getZHD(Ell, _epoch) );     
  }


  // write params
  if(_kml){
    string val = dbl2str(ell[1]*R2D, 11) + "," + dbl2str(ell[0]*R2D, 11);
    this->_default_node(node, "coordinates", val.c_str(), false);
  }

  string str_dsec = dbl2str(epoch.dsec());
    
  os <<
	  // "=FLT " << 
	  fixed << setprecision(4)
	  //<< epoch.str_ymdhms() << str_dsec.substr(2) << setprecision(4)
	  << epoch.str() << setprecision(4)
	 // << " " << epoch.sow() + epoch.dsec()
	  //<< " " << setw(13) << ell[0] * R2D          // [m]
	  //<< " " << setw(13) << ell[1] * R2D         // [m]
	  //<< " " << setw(13) << ell[2]       // [m]
	  << " " << setw(13) << xyz_ecc[0]          // [m]
	  << " " << setw(13) << xyz_ecc[1]          // [m]
	  << " " << setw(13) << xyz_ecc[2]          // [m]
     << " "   << setw(7)  << ZTD             // [m]
     << " "   << setw(7)  << grdN*1000       // [mm]
     << " "   << setw(7)  << grdE*1000       // [mm]
     << " "   << setw(15) << clk //<< (clk/CLIGHT)*1e6             // [m]
     << fixed << setprecision(5)
     << " "   << setw(9)  << Xrms            // [m]
     << " "   << setw(9)  << Yrms            // [m]
     << " "   << setw(9)  << Zrms            // [m]
     << " "   << setw(7)  << ZWDrms          // [m]
     << fixed << setprecision(0)
     << " "   << setw(2)  << nsat            // nsat
     << fixed << setprecision(1)     
     << " "   << setw(3)  << gdop            // gdop
     << fixed << setprecision(1)     
     << " "   << setw(5)  << _sig_unit       // m0
     << endl;

 return;
}

// Set output & return if smooth
// -----------------------------
void t_gpppflt::_setOut()
{
  string tmp;
  tmp = dynamic_cast<t_gsetout*>(_set)->outputs("flt");
  if( ! tmp.empty() && !_read){
   substitute(tmp,"$(rec)",_site, false);
    _flt = new t_giof;
    _flt->tsys(t_gtime::GPS);
    _flt->mask(tmp);
    _flt->append( dynamic_cast<t_gsetout*>(_set)->append() );
  }
   //LX added 
   // set enu file
  tmp = dynamic_cast<t_gsetout*>(_set)->outputs("enu");
  if (!tmp.empty() && !_read) {
	  substitute(tmp, "$(rec)", _site, false);
	  _enufile = new t_giof;
	  _enufile->tsys(t_gtime::GPS);
	  _enufile->mask(tmp);
	  _enufile->append(dynamic_cast<t_gsetout*>(_set)->append());
  }

  string port = dynamic_cast<t_gsetout*>(_set)->output_port("flt");
  if (port.size() > 0)
  {
	  t_gtcp* flt = new t_gtcp(str2int(port));
	  _maptcp.insert(make_pair(FLT_OUT, flt));
  }

  tmp = dynamic_cast<t_gsetout*>(_set)->outputs("kml");
  if( !tmp.empty() && !_read){
    substitute(tmp,"$(rec)",_site, false);
    _kml_name = tmp;
    _kml = true;
     
    xml_node node, document, placemark, multigeometry;
    xml_node root = _doc;
    node = this->_default_node(root, _root.c_str() );
     
    document = this->_default_node(node, "Document");
     
    string mode;
    string smooth = dynamic_cast<t_gsetout*>(_set)->outputs("smt");
    if(!smooth.empty()) mode = "_SMT";
    else mode = "_FLT";
    string desc = _site + mode;
    this->_default_node(document, "name", desc.c_str(), false);

	// for style in document
	xml_node doc_style1 = this->_default_node(document, "Style");
	string id1 = "P1";
	this->_default_attr(doc_style1, "id", id1);
	xml_node IconStyle = this->_default_node(doc_style1, "IconStyle");
	this->_default_node(IconStyle, "color", "ff008800");
	this->_default_node(IconStyle, "scale", "0.3");
	xml_node Icon = this->_default_node(IconStyle, "Icon");
	this->_default_node(Icon, "href", "http://maps.google.com/mapfiles/kml/pal2/icon18.png");
	
	xml_node doc_style2 = document.append_copy(doc_style1);
	string id2 = "P2";
	this->_default_attr(doc_style2, "id", id2, true);
	IconStyle = this->_default_node(doc_style2, "IconStyle");
	this->_default_node(IconStyle, "color", "fff97306", true);

	xml_node doc_style3 = document.append_copy(doc_style1);
	string id3 = "P3";
	this->_default_attr(doc_style3, "id", id3, true);
	IconStyle = this->_default_node(doc_style3, "IconStyle");
	this->_default_node(IconStyle, "color", "ff0000ff", true);

	// for style in placemark
	placemark = this->_default_node(document, "Placemark");
    this->_default_node(placemark, "name", "Trajection", false);
	xml_node style = this->_default_node(placemark, "Style");
	xml_node LineStyle = this->_default_node(style, "LineStyle");
	this->_default_node(LineStyle, "color", "ff00ffff", false);
	this->_default_node(LineStyle, "width", "2", false);
	// multigeometry = this->_default_node(placemark, "MultiGeometry");
	line = this->_default_node(placemark, "LineString");
	this->_default_node(line, "coordinates");
  }   
     
  _read = true;
  return;
}

// Run smoothing
// ---------------------
void t_gpppflt::_backflt(bool prtOut)
{
  gtrace("t_gpppflt::_backflt");
  return;
}

} // namespace


// reset ambiguity
// ------------------
void t_gpppflt::_reset_ambig()
{
   vector<int> ind = _param.delAmb();
   Matrix_rem(_Qx, ind);
}


// Check and update  ambiguity  struct when the "_prepareData" called failure
// Called when the observation is not enough or else
// ------------------
void t_gpppflt::_check_ambig()
{
  // Find number of ambiguities in the observation data
  set<string> mapPRN;
  vector<t_gsatdata>::iterator it;
  for (it = _data.begin(); it != _data.end(); it++){       // loop over all observations
    mapPRN.insert(it->sat());
  }
  
  //Delete the ambiguity parameter is the satellite is missed in the epoch
  for (size_t i = 0; i < _param.parNumber(); i++)
    {
      if (_param[i].parType == t_gpar::AMB_IF) {
        
        set<string>::iterator prnITER = mapPRN.find(_param[i].prn);
        if (prnITER == mapPRN.end()) {
          
          if (_log)
            _log->comment(2, "gppp", "AMB will be removed for missed Sat PRN " + _param[i].prn + " Epoch: " + _epoch.str_ymdhms());
          
          _newAMB.erase(_param[i].prn);
          Matrix_remRC(_Qx, _param[i].index, _param[i].index);
          _param.delParam(i);
          _param.reIndex();
          i--;    //One ambiguity is deleted, i --
        }
      }
    }
}

// reset all params
// ------------------
void t_gpppflt::_reset_param()
{
   vector<int> ind = _param.delAmb();
   Matrix_rem(_Qx, ind);           
   _initialized = false;  
}
