
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "gutils/gtypeconv.h"
#include "gall/gallprec.h"

using namespace std;

namespace gnut {  

// constructor
// ----------
t_gallprec::t_gallprec() : t_gallnav(),
   _degree_sp3(9),
   _sec(3600.0*6),
   _ref(t_gtime::GPS),
   _clkref(t_gtime::GPS),
   _clkrnx(true),
   _clksp3(false),
   _clknav(false),
   _posnav(false),
	_realtime(false)
{
  gtrace("t_gallprec::t_gallprec");

  id_type(  t_gdata::ALLPREC );
  id_group( t_gdata::GRP_EPHEM );

  _com = 1;
}


// destructor
// ----------
t_gallprec::~t_gallprec()
{
   gtrace("t_gallprec::~t_gallprec");
   
  _mapprec.clear();
}


// return satellite health
// ----------
bool t_gallprec::health( string sat, const t_gtime& t )
{
  gtrace("t_gallprec::health");
  
  if( _mapsat.size() > 0 ) return t_gallnav::health( sat, t );
  
  return true; // default healthy if no presence of NAV
}


// t .. GPS time of transmission,
// i.e. GPS time corrected for transit time (range/speed of light)
// ----------
int t_gallprec::pos( string sat, const t_gtime& t, double xyz[], double var[], double vel[], bool chk_mask )
{       
   gtrace("t_gallprec::pos");   

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();


  shared_ptr<t_geph> tmp = t_gallprec::_find( sat, t );
  if( tmp == _null ){
    for(int i = 0; i<3; i++){
                xyz[i] = 0.0;
      if( var ) var[i] = 0.0;
      if( vel ) vel[i] = 0.0;
    }     
    _gmutex.unlock();
    if( _posnav && t_gallnav::pos( sat, t, xyz, var, vel, chk_mask ) >= 0 ){ return 1; } // cout << "USING POS NAV:\n"; return 1; }
    return -1;
  }

  int irc = tmp->pos( t, xyz, var, vel, _chk_health && chk_mask );
  _gmutex.unlock(); return irc;
}


bool t_gallprec::posnew(string sat, const t_gtime& t, double xyz[])
{
	gtrace("t_gallprec::posnew");

#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	map<t_gtime, t_map_dat>::iterator itReq = _mapsp3[sat].lower_bound(t); // 1st equal|greater [than t]
	if (itReq == _mapsp3[sat].end())
	{
		_gmutex.unlock();
		return false;
	}
	t_map_dat  data_t;
	if (abs(itReq->first - t)<0.01)
	{
		data_t = itReq->second;
		if (data_t.find("X") != data_t.end())
		{
			xyz[0] = data_t["X"];
			xyz[1] = data_t["Y"];
			xyz[2] = data_t["Z"];
		}
		_gmutex.unlock();
		return true;
    }

	double var[3];
	double vel[3];
	bool chk_mask;
	shared_ptr<t_geph> tmp = t_gallprec::_find(sat, t);
	if (tmp == _null) {
		for (int i = 0; i < 3; i++) {
			xyz[i] = 0.0;
			if (var) var[i] = 0.0;
			if (vel) vel[i] = 0.0;
		}
		_gmutex.unlock();
		if (_posnav && t_gallnav::pos(sat, t, xyz, var, vel, chk_mask) >= 0) { return 1; } // cout << "USING POS NAV:\n"; return 1; }
		return false;
	}

	int irc = tmp->pos(t, xyz, var, vel, _chk_health && chk_mask);
	_gmutex.unlock(); return irc==0;
}
// t .. GPS time of transmission,
// i.e. GPS time corrected for transit time (range/speed of light)
// ----------
int t_gallprec::pos_int( string sat, const t_gtime& t, double xyz[], double var[], double vel[] )
{   
   gtrace("t_gallprec::pos_int");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  shared_ptr<t_geph> tmp = t_gallprec::_find( sat, t );

  if( tmp == _null ){
    for(int i = 0; i<3; i++){
                xyz[i] = 0.0;
      if( var ) var[i] = 0.0;
      if( vel ) vel[i] = 0.0;
    }     
    _gmutex.unlock(); return -1;
  }

  int irc = dynamic_pointer_cast<t_gephprec>(tmp)->pos_int( t, xyz, var, vel );
  _gmutex.unlock(); return irc;
}

// t .. GPS time of transmission,
// i.e. GPS time corrected for transit time (range/speed of light)
// POSITION directly interpolated from array of SP3 positions
// ----------
int t_gallprec::pos_alt( string sat, const t_gtime& t, double xyz[], double var[], double vel[] )
{      
   gtrace("t_gallprec::pos_alt");
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  if( _get_crddata( sat, t ) < 0 ){
    for(int i = 0; i<3; i++){
                xyz[i] = 0.0;
      if( var ) var[i] = 0.0;
      if( vel ) vel[i] = 0.0;
    }     
    _gmutex.unlock(); return -1;
  }

  t_gpoly poly;
  poly.interpolate( _PT, _X, t.diff(_ref), xyz[0], var[0] );
  poly.interpolate( _PT, _Y, t.diff(_ref), xyz[1], var[1] );
  poly.interpolate( _PT, _Z, t.diff(_ref), xyz[2], var[2] );

  _gmutex.unlock(); return 1;
}


// t .. GPS time of transmission,
// i.e. GPS time corrected for transit time (range/speed of light)
// NOT CORRECTED FOR 2nd relativistic effect !!!
// ----------
int t_gallprec::clk( string sat, const t_gtime& t, double* clk, double* var, double* dclk, bool chk_mask )
{   gtrace("t_gallprec::clk");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  if( !_clkrnx || _get_clkdata( sat, t ) < 0 ){
              *clk  = 0.0;
    if( var ) *var  = 0.0;
    if(dclk ) *dclk = 0.0;

    _gmutex.unlock();

    // REMOVED, used instead gallnav.health() function explicitly!
    // ==========================================================
//     //Before calling clk_int, t_gallnav::clk is called to ensure the satellite is healthy
//    int irc = t_gallnav::clk(sat, t, clk, var, dclk);
//    if (irc >= 0){
//      *clk = 0.0;
//      if( var  ) *var  = 0.0;
//      if( dclk ) *dclk = 0.0;
//    }
//    else{
//      if (_log) _log->comment(0, "Navigation information is not valid for satellites: " + sat + t.str(" %Y-%m-%d %H:%M:%S"));
//      return -1;
//    }

//  don't call clk_sp3() -> a recurrent call due to clk() will happen !
    if( _clksp3 &&  this->clk_int( sat, t, clk, var, dclk           ) >= 0 ){ return 1; }
    if( _clknav && t_gallnav::clk( sat, t, clk, var, dclk, chk_mask ) >= 0 ){ return 1; }
    return -1;
  }
   
  t_gpoly poly;
  poly.interpolate( _CT, _C, t.diff(_clkref), *clk, *dclk );
  // added by zhshen
  *dclk = *dclk / (_CT.back() - _CT.front());

  _gmutex.unlock(); return 1;
}

int t_gallprec::clk_ifcb(string sat, const t_gtime& t, double* clk, double* var, double* dclk, bool chk_mask)
{
    gtrace("t_gallprec::clk");

#ifdef BMUTEX   
    boost::mutex::scoped_lock lock(_mutex);
#endif
    _gmutex.lock();

    if (!_clkrnx || _get_clkdata(sat, t) < 0) {
        *clk = 0.0;
        if (var) *var = 0.0;
        if (dclk) *dclk = 0.0;

        _gmutex.unlock();

        return -1;
    }

    t_gpoly poly;
    poly.interpolate(_CT, _IFCB_F3, t.diff(_clkref), *clk, *dclk);

    _gmutex.unlock(); return 1;
}

int t_gallprec::clk_cdr(string sat, const t_gtime& t, double* clk, double* var, double* dclk, double* ifcb, bool chk_mask)
{
    if (_mapclk.find(sat) == _mapclk.end()) return -1;

    map<t_gtime, t_map_dat>::iterator itReq = _mapclk[sat].lower_bound(t); // 1st equal|greater [than t]   

    if (itReq == _mapclk[sat].end()) return -1;

    if (itReq->first > t) return -1;
    else
    {
        *clk  = itReq->second["C0"];
        *var  = 0.0;
        *dclk = 0.0;
        *ifcb = itReq->second["IFCB_F3"];
    }
    return 0;
}


// t .. GPS time of transmission,
// i.e. GPS time corrected for transit time (range/speed of light)
// NOT CORRECTED FOR 2nd relativistic effect !!!
// ----------
int t_gallprec::clk_sp3( string sat, const t_gtime& t, double* clk, double* var, double* dclk )
{   
   gtrace("t_gallprec::clk_sp3");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  shared_ptr<t_geph> tmp = t_gallprec::_find( sat, t );

  if( tmp == _null ){
              *clk  = 0.0;
    if( var ) *var  = 0.0;
    if(dclk ) *dclk = 0.0;
    _gmutex.unlock(); return -1; 
  }
   
  int irc = tmp->clk( t, clk, var, dclk );
  _gmutex.unlock(); return irc;
}


// t .. GPS time of transmission,
// i.e. GPS time corrected for transit time (range/speed of light)
// NOT CORRECTED FOR 2nd relativistic effect !!!
// ----------
int t_gallprec::clk_int( string sat, const t_gtime& t, double* clk, double* var, double* dclk )
{   
   gtrace("t_gallprec::clk_int");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  shared_ptr<t_geph> tmp = t_gallprec::_find( sat, t );

  if( tmp == _null ){
              *clk  = 0.0;
    if( var ) *var  = 0.0;
    if(dclk ) *dclk = 0.0;
    _gmutex.unlock(); return -1;
  }
  int irc = dynamic_pointer_cast<t_gephprec>(tmp)->clk_int( t, clk, var, dclk );
   
  _gmutex.unlock(); return irc;

}



// t .. GPS time of transmission,
// i.e. GPS time corrected for transit time (range/speed of light)
// CLOCK CORRECTIONS directly interpolated from array of CLOCK-RINEX
// NOT CORRECTED FOR 2nd relativistic effect !!!
// ----------
int t_gallprec::clk_alt( string sat, const t_gtime& t, double* clk, double* var, double* dclk )
{   
   gtrace("t_gallprec::clk_alt");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  if( _get_clkdata( sat, t ) < 0 ){
              *clk  = 0.0;
    if( var ) *var  = 0.0;
    if(dclk ) *dclk = 0.0;
    _gmutex.unlock(); return -1;
  }
   
  t_gpoly poly;
  poly.interpolate( _CT, _C, t.diff(_clkref), *clk, *dclk );
   
  // add relativistic correction ! -2rs/c/c
//  *clk -= 2.0*;
  _gmutex.unlock(); return 1;
}

//yjqin
void t_gallprec::add_interval(int intv)
{
	_intv = intv;
}


void t_gallprec::add_interval(string sat,int intv)
{
	_intvm[sat] = intv;
}


void t_gallprec::add_agency( string agency)
{
	_agency = agency;
}

void t_gallprec::add_clk_interval(double intv)
{
	_udclkInt = intv;
}

int t_gallprec::intv()
{
	return _intv;

}

int t_gallprec::intv(string sat)
{
	return _intvm[sat];

}

string t_gallprec::get_agency()
{
	return _agency;

}

t_gtime t_gallprec::get_beg()
{
	return _tbeg;
}

t_gtime t_gallprec::get_end()
{
	return _tend;
}

vector<string> t_gallprec::get_sats()
{
	return _sats;
}

vector<string> t_gallprec::get_sat3()
{
	return _sat3;
}

string t_gallprec::get_data_type()
{
	return _datatype;
}

string t_gallprec::get_sat_type()
{
	return _sattype;
}

int t_gallprec::get_sp3_size()
{
	return _mapsp3.size();
}

bool t_gallprec::get_pos_vel(string sat, t_gtime epoch, double pos[3], double vel[3],int &obs_num)
{
	if (_mapsp3[sat].size() == 0) return false; // make_shared<t_geph>();
	// if not exists satellite not in cache
	t_map_dat data = _mapsp3[sat][epoch];
	if (data.size() != 0) {
		pos[0] = data["X"];
		pos[1] = data["Y"];
		pos[2] = data["Z"];
		obs_num = data["OBS"];
		if (vel)
		{
			vel[0] = data["VX"];
			vel[1] = data["VY"];
			vel[2] = data["VZ"];
		}
		
		return true;
	}

    return false;
}

void t_gallprec::set_beg(t_gtime beg)
{
	_tbeg = beg;
}

void t_gallprec::set_end(t_gtime end)
{
	_tend = end;
}

void t_gallprec::set_sat(vector<string> sat)
{
	_sats = sat;
}

void t_gallprec::set_sat3(vector<string> sat)
{
	_sat3 = sat;
}

void t_gallprec::set_data_type(string type)
{
	_datatype = type;
}

void t_gallprec::set_sat_type(string type)
{
	_sattype = type;
}

void t_gallprec::set_agency(string agency)
{
	_agency = agency;
}

void t_gallprec::add_orb_interval(double intv)
{
	_udorbInt = intv;
}

int t_gallprec::add_delta_pos_vel(string sat, const t_gtime & ep, const int& iod, const t_gtriple & dxyz, const t_gtriple & dvxyz)
{
	gtrace("t_gallprec::add_delta_pos_vel");

#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	if (_overwrite || _mapsp3[sat].find(ep) == _mapsp3[sat].end()) {

		_mapsp3[sat][ep]["IOD"] = iod;
		_mapsp3[sat][ep]["DX"] = dxyz[0];
		_mapsp3[sat][ep]["DY"] = dxyz[1];
		_mapsp3[sat][ep]["DZ"] = dxyz[2];

		_mapsp3[sat][ep]["DVX"] = dvxyz[0];
		_mapsp3[sat][ep]["DVY"] = dvxyz[1];
		_mapsp3[sat][ep]["DVZ"] = dvxyz[2];
		_realtime = true;
		if (_log && _log->verb() >= 4) {
			ostringstream lg;
			lg << "add sat dxyz,dvxyz " << sat << fixed << setprecision(6)
				<< " " << ep.str("%Y-%m-%d %H:%M:%S[%T] ")
				<< " " << setw(16) << _mapsp3[sat][ep]["DX"]
				<< " " << setw(16) << _mapsp3[sat][ep]["DY"]
				<< " " << setw(16) << _mapsp3[sat][ep]["DZ"]
				<< " " << setw(16) << _mapsp3[sat][ep]["DVX"]
				<< " " << setw(16) << _mapsp3[sat][ep]["DVY"]
				<< " " << setw(16) << _mapsp3[sat][ep]["DVZ"];
			_log->comment(4, "gallprec", lg.str());
		}

	}
	else {
		if (_log) _log->comment(3, "gallprec", ep.str_ymdhms(sat + " not overwriting position for "));
		_gmutex.unlock(); return -1;
	}

	_gmutex.unlock();
	return 0;
}
//yjqin



// add position
// ----------
//int t_gallprec::addpos( string sat, const t_gtime& ep, double xyzt[], double dxyz[] )
int t_gallprec::addpos( string sat, const t_gtime& ep, t_gtriple  xyz, double t,
			                               t_gtriple dxyz, double dt )
{
   gtrace("t_gallprec::addpos");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  
  if( _overwrite || _mapsp3[sat].find(ep) == _mapsp3[sat].end() ){

    _mapsp3[sat][ep]["X"]  = xyz[0];
    _mapsp3[sat][ep]["Y"]  = xyz[1];
    _mapsp3[sat][ep]["Z"]  = xyz[2];
    _mapsp3[sat][ep]["C"]  = t;
   
    _mapsp3[sat][ep]["SX"] = dxyz[0];
    _mapsp3[sat][ep]["SY"] = dxyz[1];
    _mapsp3[sat][ep]["SZ"] = dxyz[2];
    _mapsp3[sat][ep]["SC"] = dt;
     
    if( _log && _log->verb() >= 4 ){
      ostringstream lg;
      lg << "add sat xyz, t " << sat << fixed << setprecision(6)
         << " " << ep.str("%Y-%m-%d %H:%M:%S[%T] ")
         << " " << setw(16) << _mapsp3[sat][ep]["X"]
	 << " " << setw(16) << _mapsp3[sat][ep]["Y"]
	 << " " << setw(16) << _mapsp3[sat][ep]["Z"] << setprecision(9)
	 << " " << setw(16) << _mapsp3[sat][ep]["C"];
      _log->comment(4,"gallprec",lg.str());
    }
     
  }else{
    if(_log) _log->comment(3,"gallprec",ep.str_ymdhms(sat + " not overwriting position for "));
    _gmutex.unlock(); return -1;
  }  
  _gmutex.unlock(); return 0;
}


// add velocity
// ----------
int t_gallprec::addvel( string sat, const t_gtime& ep, double xyzt[], double dxyz[] )
{
   gtrace("t_gallprec::addvel");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  
  if( _overwrite || _mapsp3[sat].find(ep) == _mapsp3[sat].end() ){

    _mapsp3[sat][ep]["VX"] = xyzt[0];
    _mapsp3[sat][ep]["VY"] = xyzt[1];
    _mapsp3[sat][ep]["VZ"] = xyzt[2];
    _mapsp3[sat][ep]["VC"] = xyzt[3];
   
    _mapsp3[sat][ep]["SVX"] = dxyz[0];
    _mapsp3[sat][ep]["SVY"] = dxyz[1];
    _mapsp3[sat][ep]["SVZ"] = dxyz[2];
    _mapsp3[sat][ep]["SVC"] = dxyz[3];

  }else{
    if( _log && _log->verb() >= 4 ){
      ostringstream lg;
      lg << "add sat dxyz,dt " << sat << fixed << setprecision(6)
         << " " << ep.str("%Y-%m-%d %H:%M:%S[%T] ")
         << " " << setw(16) << _mapsp3[sat][ep]["VX"]
	     << " " << setw(16) << _mapsp3[sat][ep]["VY"]
	     << " " << setw(16) << _mapsp3[sat][ep]["VZ"] << setprecision(9)
	     << " " << setw(16) << _mapsp3[sat][ep]["VT"];
      _log->comment(4,"gallprec",lg.str());
    }
    if( _log && _log->verb() >= 3 ) _log->comment(3,"gallprec",ep.str_ymdhms(sat + " not overwriting velocity for "));
    _gmutex.unlock(); return -1;
  }
  _gmutex.unlock(); return 0;
}

int t_gallprec::add_delta_clk(string sat, const t_gtime & ep, const int & iod, const double & dt, const double & dot_dt, const double & dot_dot_dt)
{
	gtrace("t_gallprec::add_delta_clk");

#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	if (_overwrite || _mapclk[sat].find(ep) == _mapclk[sat].end()) {

		_mapclk[sat][ep]["IOD"] = iod;
		_mapclk[sat][ep]["DCLK"] = dt;
		_mapclk[sat][ep]["DOTCLK"] = dot_dt;
		_mapclk[sat][ep]["DOTDOTCLK"] = dot_dot_dt;
		_realtime = true;

		if (_log && _log->verb() >= 4) {
			ostringstream lg;
			lg << "add sat dclk " << sat << fixed << setprecision(6)
				<< " " << ep.str("%Y-%m-%d %H:%M:%S[%T] ")
				<< " " << setw(16) << _mapclk[sat][ep]["DCLK"]
				<< " " << setw(16) << _mapclk[sat][ep]["DOTCLK"]
				<< " " << setw(16) << _mapclk[sat][ep]["DOTDOTCLK"];
			_log->comment(4, "gallprec", lg.str());
		}

	}
	else {
		if (_log) _log->comment(3, "gallprec", ep.str_ymdhms(sat + " not overwriting clk for "));
		_gmutex.unlock(); return -1;
	}

	_gmutex.unlock();
	return 0;
}

int t_gallprec::get_pos_clk_correction(const string & sat, const t_gtime & t, double *xyz, double *vxyz, double &clk, double &dclk)
{
	gtrace("t_gallprec::get_pos_clk_correction");

#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	if (_get_delta_pos_vel(sat, t) < 0 ) {
		xyz[0] = xyz[1] = xyz[2] = 0.0;
		if (_log) _log->comment(3, "gallprec", t.str_ymdhms(sat + " not ssr position correction for "));
	}	
	if (_get_delta_clk(sat, t) < 0) {
		clk = 0.0;
		if (_log) _log->comment(3, "gallprec", t.str_ymdhms(sat + " not ssr clk correction for "));
	}

	t_gpoly poly;
	poly.interpolate(_CTCorr, _CCorr, t.diff(_clkref), clk, dclk);

	// the vxyz&dclk may be changed to another way, reserved temporary
	poly.interpolate(_PTCorr, _XCorr, t.diff(_ref), xyz[0], vxyz[0]);
	poly.interpolate(_PTCorr, _YCorr, t.diff(_ref), xyz[1], vxyz[1]);
	poly.interpolate(_PTCorr, _ZCorr, t.diff(_ref), xyz[2], vxyz[2]);
	_gmutex.unlock();
	return 1;
}


// added by zhShen
int t_gallprec::get_pos_clk_correction(const string & sat, const t_gtime & t, int iod, double * xyz, double * vxyz, double & clk, double & dclk)
{
	gtrace("t_gallprec::get_pos_clk_correction");
	//return -1;
#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	t_gtime tRefOrb, tRefClk;
	t_map_dat orbcorr, clkcorr;
	int f1 = 0, f2 = 0;

	if (f1 = _get_delta_pos_vel(sat, t, iod, tRefOrb, orbcorr) < 0) {
		xyz[0] = xyz[1] = xyz[2] = 0.0;
		if (_log) _log->comment(3, "gallprec", t.str_ymdhms(sat + " not ssr position correction for "));
	}
	if (f2 = _get_delta_clk(sat, t, iod, tRefClk, clkcorr) < 0) {
		clk = 0.0; dclk = 0.0;
		if (_log) _log->comment(3, "gallprec", t.str_ymdhms(sat + " not ssr clk correction for "));
	}
	if (f1 == -1 && f2 == -1) {
		_gmutex.unlock(); return -1;
	}

	double tdifforb = t.diff(tRefOrb);
	if (_udorbInt)tdifforb -= 0.5*_udorbInt;
	if (fabs(tdifforb) > 5 * 60) { _gmutex.unlock(); return -1; }
	// get the position correction, no ratation
	xyz[0] = orbcorr["DX"] + orbcorr["DVX"] * tdifforb;
	xyz[1] = orbcorr["DY"] + orbcorr["DVY"] * tdifforb;
	xyz[2] = orbcorr["DZ"] + orbcorr["DVZ"] * tdifforb;

	// get the velocity correction, no ratation
	vxyz[0] = orbcorr["DVX"];
	vxyz[1] = orbcorr["DVY"];
	vxyz[2] = orbcorr["DVZ"];

	double tdiffclk = t.diff(tRefClk);
	if (_udclkInt)tdiffclk -= 0.5*_udclkInt;
	if (fabs(tdiffclk) > 10) { _gmutex.unlock(); return -1; }
	// get the clock correction
	clk = clkcorr["DCLK"] + clkcorr["DOTCLK"] * tdiffclk + clkcorr["DOTDOCLK"] * tdiffclk*tdiffclk;
	dclk = clkcorr["DOTCLK"] + clkcorr["DOTDOCLK"] * tdiffclk;
	
#ifdef DEBUG

	cout << t.str() << "  " << sat << "  " << tdifforb << "  " << tdiffclk << "  ";

	cout << fixed << setprecision(3) << xyz[0] << "  "
		<< xyz[1] << "  "
		<< xyz[2] << "  "
		<< vxyz[0] << "  "
		<< vxyz[1] << "  "
		<< vxyz[2] << "  "
		<< clk << "  "
		<< dclk << "  " << endl;

#endif // DEBUG


	_gmutex.unlock();
	return 1;
}

bool t_gallprec::corr_avali()
{
	return (_mapsp3.size() > 0 && _mapclk.size() > 0);
}

bool t_gallprec::corr_avali(t_gtime now)
{
	if (now < _tend || now.diff(_tend) < _udorbInt)return true;
	else
	{
		return false;
	}
}

//yjqin
int t_gallprec::add_pos_vel(string sat, const t_gtime& ep, t_gtriple xyz, double t, t_gtriple dxyz, double dt, double xyzt[4], double var[4], int obs_num)
{
	gtrace("t_gallprec::add_pos_vel");

#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();
	if (_overwrite || _mapsp3[sat].find(ep) == _mapsp3[sat].end()) {

		_mapsp3[sat][ep]["X"] = xyz[0];
		_mapsp3[sat][ep]["Y"] = xyz[1];
		_mapsp3[sat][ep]["Z"] = xyz[2];
		_mapsp3[sat][ep]["C"] = t;

		_mapsp3[sat][ep]["SX"] = dxyz[0];
		_mapsp3[sat][ep]["SY"] = dxyz[1];
		_mapsp3[sat][ep]["SZ"] = dxyz[2];
		_mapsp3[sat][ep]["SC"] = dt;

		_mapsp3[sat][ep]["VX"] = xyzt[0];
		_mapsp3[sat][ep]["VY"] = xyzt[1];
		_mapsp3[sat][ep]["VZ"] = xyzt[2];
		_mapsp3[sat][ep]["VC"] = xyzt[3];

		_mapsp3[sat][ep]["SVX"] = var[0];
		_mapsp3[sat][ep]["SVY"] = var[1];
		_mapsp3[sat][ep]["SVZ"] = var[2];
		_mapsp3[sat][ep]["SVC"] = var[3];
		_mapsp3[sat][ep]["OBS"] = obs_num;
	}
	else {
		if (_log) _log->comment(3, "gallprec", ep.str_ymdhms(sat + " not overwriting position and velocity for "));
		_gmutex.unlock(); return -1;
	}
	_gmutex.unlock(); return 0;
}
// add clocks
// ----------
int t_gallprec::addclk( string sat, const t_gtime& ep, double clk[], double dxyz[] )
{
   gtrace("t_gallprec::addclk");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  
  if( _overwrite || _mapclk[sat].find(ep) == _mapclk[sat].end() ){

	// change by ZhengHJ for slow
	  t_map_dat data = {
		  {"C0",clk[0]},
		  {"C1",clk[1]},
		  {"C2",clk[2]},
		  {"SC0",dxyz[0]},
		  {"SC1",dxyz[1]},
		  {"SC2",dxyz[2]},
		};
	  _mapclk[sat][ep] = data;
    //_mapclk[sat][ep]["C0"] = clk[0];
    //_mapclk[sat][ep]["C1"] = clk[1];
    //_mapclk[sat][ep]["C2"] = clk[2];
   
    //_mapclk[sat][ep]["SC0"] = dxyz[0];
    //_mapclk[sat][ep]["SC1"] = dxyz[1];
    //_mapclk[sat][ep]["SC2"] = dxyz[2];


	if (sat.length() == 3) {
		_clk_type_list.insert(AS);
	}
	else {
		_clk_type_list.insert(AR);
	}

  }else{
    if( _log && _log->verb() >= 3 ) _log->comment(3,"gallprec",ep.str_ymdhms(sat + " not overwriting clocks for "));
    _gmutex.unlock(); return 1;
  }

#ifdef DEBUG
  if( _log ) _log->comment(0,"gallprec",ep.str_ymdhms(sat + " add sat clk " + int2str(_mapclk[sat].size()) + " for " ));
#endif
  _gmutex.unlock(); return 0;
}

int t_gallprec::addclk_tri(string sat, const t_gtime& ep, double clk[], double dxyz[])
{
    gtrace("t_gallprec::addclk_tri");

#ifdef BMUTEX   
    boost::mutex::scoped_lock lock(_mutex);
#endif
    _gmutex.lock();

    if (_overwrite || _mapclk[sat].find(ep) == _mapclk[sat].end()) {

        // change by ZhengHJ for slow
        // add ifcb by xiongyun
        t_map_dat data = {
            {"C0",clk[0]},
            {"C1",clk[1]},
            {"C2",clk[2]},
            {"IFCB_F3",clk[3]},
            {"SC0",dxyz[0]},
            {"SC1",dxyz[1]},
            {"SC2",dxyz[2]},
        };
        _mapclk[sat][ep] = data;
        //_mapclk[sat][ep]["C0"] = clk[0];
        //_mapclk[sat][ep]["C1"] = clk[1];
        //_mapclk[sat][ep]["C2"] = clk[2];

        //_mapclk[sat][ep]["SC0"] = dxyz[0];
        //_mapclk[sat][ep]["SC1"] = dxyz[1];
        //_mapclk[sat][ep]["SC2"] = dxyz[2];


        if (sat.length() == 3) {
            _clk_type_list.insert(AS);
        }
        else {
            _clk_type_list.insert(AR);
        }

    }
    else {
        if (_log && _log->verb() >= 3) _log->comment(3, "gallprec", ep.str_ymdhms(sat + " not overwriting clocks for "));
        _gmutex.unlock(); return 1;
    }

#ifdef DEBUG
    if (_log) _log->comment(0, "gallprec", ep.str_ymdhms(sat + " add sat clk " + int2str(_mapclk[sat].size()) + " for "));
#endif
    _gmutex.unlock(); return 0;
}


// return number of epochs
// ----------
unsigned int t_gallprec::nepochs( const string& prn )
{
   gtrace("t_gallprec::nepochs");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  unsigned int tmp = 0;
  if( _mapsp3.find(prn) != _mapsp3.end() ) tmp = _mapsp3[prn].size();

  _gmutex.unlock(); return tmp;
}


// return list of available satellites
// ----------
set<string> t_gallprec::satellites()
{
  gtrace("t_gallprec::satellites");

  set<string> all_sat = t_gallnav::satellites();

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  t_map_prn::const_iterator itSP3 = _mapsp3.begin();
  while( itSP3 != _mapsp3.end() ){
    if( all_sat.find(itSP3->first) == all_sat.end() ) all_sat.insert(itSP3->first);
    itSP3++;
  }

  _gmutex.unlock(); return all_sat;
}

// clean clk function
// clean sp3 derived from gnav
// ----------
void t_gallprec::clean_outer( const t_gtime& beg, const t_gtime& end )
{
   gtrace("t_gallprec::clean_outer");   
   
  if( end < beg ) return;

  // first clean navigation messages
  t_gallnav::clean_outer( beg, end );

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  // prec ephemeris - loop over all satellites
  // -----------------------------------------
  map<string,t_map_epo>::const_iterator itPRN = _mapsp3.begin();
  while( itPRN != _mapsp3.end() ){
    string prn = itPRN->first;

    // find and CLEAN all data (epochs) out of the specified period !
    map<t_gtime,t_map_dat>::iterator it;     
    map<t_gtime,t_map_dat>::iterator itFirst = _mapsp3[prn].begin();
    map<t_gtime,t_map_dat>::iterator itLast  = _mapsp3[prn].end();
    map<t_gtime,t_map_dat>::iterator itBeg   = _mapsp3[prn].lower_bound(beg);  // greater|equal
    map<t_gtime,t_map_dat>::iterator itEnd   = _mapsp3[prn].upper_bound(end);  // greater only!


    // remove before BEGIN request
    if( itBeg != itFirst ){
       
      // begin is last
      if( itBeg == itLast ){
        itBeg--;
        if( _log && _log->verb() >= 2 ) _log->comment(2,"gallprec",itBeg->first.str_ymdhms(prn + " sp3 removed before "));
        _mapsp3[prn].erase(itFirst,itLast);
       
      // begin is not last
      }else{
        if( _log && _log->verb() >= 2 ) _log->comment(2,"gallprec",itBeg->first.str_ymdhms(prn + " sp3 removed before "));
         _mapsp3[prn].erase(itFirst,itBeg);
      }
    }

    // remove after END request
    if( itEnd != itLast ){
        if( _log && _log->verb() >= 2 ) _log->comment(2,"gallprec",itEnd->first.str_ymdhms(prn + " sp3 removed after  "));
        _mapsp3[prn].erase(itEnd,itLast);
    }
    itPRN++;
  }     

  // prec clocks - loop over all satellites
  // --------------------------------------
  itPRN = _mapclk.begin();
  while( itPRN != _mapclk.end() ){
    string prn = itPRN->first;

    // find and CLEAN all data (epochs) out of the specified period !
    map<t_gtime,t_map_dat>::iterator it;     
    map<t_gtime,t_map_dat>::iterator itFirst = _mapclk[prn].begin();
    map<t_gtime,t_map_dat>::iterator itLast  = _mapclk[prn].end();
    map<t_gtime,t_map_dat>::iterator itBeg   = _mapclk[prn].lower_bound(beg);  // greater|equal
    map<t_gtime,t_map_dat>::iterator itEnd   = _mapclk[prn].upper_bound(end);  // greater only!

    // remove before BEGIN request
    if( itBeg != itFirst ){
       
      // begin is last
      if( itBeg == itLast ){
        itBeg--;
        if( _log && _log->verb() >= 2 ) _log->comment(3,"gallprec",itFirst->first.str_ymdhms(prn + " clk removed from ")
                                            + itBeg->first.str_ymdhms(" to "));

         _mapclk[prn].erase(itFirst,itLast);
       
      // begin is not last
      }else{
        if( _log && _log->verb() >= 2 ) _log->comment(3,"gallprec",itFirst->first.str_ymdhms(prn + " clk removed from ")
                                            + itBeg->first.str_ymdhms(" to "));
         _mapclk[prn].erase(itFirst,itBeg);     
      }
    }

    // remove after END request
    if( itEnd != itLast ){ // && ++itEnd != itLast ){
      if( _log && _log->verb() >= 2 ) _log->comment(2,"gallprec",itEnd->first.str_ymdhms(prn + " clk removed after "));

      _mapclk[prn].erase(itEnd,itLast);
    }
    itPRN++;
  }     
  _gmutex.unlock(); return;
}


// return first epoch of sp3 position/clocks
// ----------
t_gtime t_gallprec::beg_data( string prn )
{
   gtrace("t_gallprec::beg_data");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtime tmp = LAST_TIME;
  if( ! prn.empty() ){
    if( _mapsp3.find(prn) != _mapsp3.end() ) tmp = _mapsp3[prn].begin()->first;
  }else{
    for( auto itSAT = _mapsp3.begin(); itSAT != _mapsp3.end(); ++itSAT ){
      for( auto it = itSAT->second.begin(); it != itSAT->second.end(); ++it ){
	if( _mapsp3[itSAT->first].begin()->first < tmp ) tmp = _mapsp3[itSAT->first].begin()->first;
      }
    }
  }

  _gmutex.unlock(); return tmp;
}


// return last epoch of sp3 position/clocks
// ----------
t_gtime t_gallprec::end_data( string prn )
{
   gtrace("t_gallprec::end_data");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtime tmp = FIRST_TIME;
  if( ! prn.empty() ){
    if( _mapsp3.find(prn) != _mapsp3.end() ) tmp = _mapsp3[prn].rbegin()->first;
  }else{
    for( auto itSAT = _mapsp3.begin(); itSAT != _mapsp3.end(); ++itSAT ){
      for( auto it = itSAT->second.begin(); it != itSAT->second.end(); ++it ){
	if( _mapsp3[itSAT->first].rbegin()->first > tmp ) tmp = _mapsp3[itSAT->first].rbegin()->first;
      }
    }
  }

  _gmutex.unlock(); return tmp;
}


// return first epoch of rinex clocks
// ----------
t_gtime t_gallprec::beg_clk( string prn )
{
   gtrace("t_gallprec::beg_clk");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtime tmp = LAST_TIME;
  if( ! prn.empty() ){
    if( _mapclk.find(prn) != _mapclk.end() ) tmp = _mapclk[prn].begin()->first;
  }else{
    for( auto itSAT = _mapclk.begin(); itSAT != _mapclk.end(); ++itSAT ){
      for( auto it = itSAT->second.begin(); it != itSAT->second.end(); ++it ){
	if( _mapclk[itSAT->first].begin()->first < tmp ) tmp = _mapclk[itSAT->first].begin()->first;
      }
    }
  }

  _gmutex.unlock(); return tmp;
}


// return last epoch of rinex clocks
// ----------
t_gtime t_gallprec::end_clk( string prn )
{
   gtrace("t_gallprec::end_clk");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtime tmp = FIRST_TIME;
  if( ! prn.empty() ){
    if( _mapclk.find(prn) != _mapclk.end() ) tmp = _mapclk[prn].rbegin()->first;
  }else{
    for( auto itSAT = _mapclk.begin(); itSAT != _mapclk.end(); ++itSAT ){
      for( auto it = itSAT->second.begin(); it != itSAT->second.end(); ++it ){
	if( _mapclk[itSAT->first].rbegin()->first > tmp ) tmp = _mapclk[itSAT->first].rbegin()->first;
      }
    }
  }

  _gmutex.unlock(); return tmp;
}


// return first epoch of polynomials
// ----------
t_gtime t_gallprec::beg_prec( string prn )
{
   gtrace("t_gallprec::beg_prec");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtime tmp = LAST_TIME;
  if( ! prn.empty() ){
    if( _mapprec.find(prn) != _mapprec.end() ) tmp = _mapprec[prn].begin()->first;
  }else{
    for( auto itSAT = _mapprec.begin(); itSAT != _mapprec.end(); ++itSAT ){
      for( auto it = itSAT->second.begin(); it != itSAT->second.end(); ++it ){
	if( _mapprec[itSAT->first].begin()->first < tmp ) tmp = _mapprec[itSAT->first].begin()->first;
      }
    }
  }

  _gmutex.unlock(); return tmp;
}


// return last epoch of polynomials
// ----------
t_gtime t_gallprec::end_prec( string prn )
{
   gtrace("t_gallprec::end_prec");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  
  t_gtime tmp = FIRST_TIME;
  if( ! prn.empty() ){
    if( _mapprec.find(prn) != _mapprec.end() ) tmp = _mapprec[prn].rbegin()->first;
  }else{
    for( auto itSAT = _mapprec.begin(); itSAT != _mapprec.end(); ++itSAT ){
      for( auto it = itSAT->second.begin(); it != itSAT->second.end(); ++it ){
	if( _mapprec[itSAT->first].rbegin()->first > tmp ) tmp = _mapprec[itSAT->first].rbegin()->first;
      }
    }
  }

  _gmutex.unlock(); return tmp;
}

   
// Approximative position
// ----------   
int t_gallprec::nav(string sat, const t_gtime& t, double  xyz[3], double  var[3], double  vel[3], bool chk_mask )
{
  gtrace("t_gallprec::nav");

  int fitdat = 24;           // fitting samples
  int fitdeg = 12;           // fitting degree

  for(int i = 0; i<3; i++){
              xyz[i] = 0.0;
    if( var ) var[i] = 0.0;
    if( vel ) vel[i] = 0.0;
  }
   
  // alternative use of gnav
  if( _mapsp3[sat].size() == 0 )
    return (( _clknav && t_gallnav::nav( sat, t, xyz, var, vel, chk_mask ) >= 0 ) ? 1 : -1 );

  _gmutex.lock();

  t_gtime beg(_mapsp3[sat].begin()->first);
  t_gtime end(_mapsp3[sat].rbegin()->first);

  if( t < beg-900 || t > end+900 ){
    if(_log) _log->comment(1,"gallprec","no position available prior/after "+beg.str_ymdhms()+"/"+end.str_ymdhms() );
     _gmutex.unlock(); return (( _clknav && t_gallnav::nav( sat, t, xyz, var, vel ) >= 0 ) ? 1 : -1 );
  }

  if( _poly_beg.find(sat) == _poly_beg.end() ) _poly_beg[sat] = LAST_TIME;
  if( _poly_end.find(sat) == _poly_end.end() ) _poly_end[sat] = FIRST_TIME;

  // use existing approximative estimates from cached polynomials
  if( !( t > _poly_end[sat] || t < _poly_beg[sat] ) )
  { _sec = _poly_end[sat] - _poly_beg[sat];
    _ref = _poly_beg[sat] + _sec/2;
    
    _poly_x[sat].evaluate( t.diff(_ref)/_sec, 0, xyz[0] );
    _poly_y[sat].evaluate( t.diff(_ref)/_sec, 0, xyz[1] );
    _poly_z[sat].evaluate( t.diff(_ref)/_sec, 0, xyz[2] );

#ifdef DEBUG     
  cout << " PRN: " << sat
       << " req: " << t.str_ymdhms()
       << " ref: " << _ref.str_ymdhms()
       << " sec: " << _sec
       << " beg: " << _poly_beg[sat].str_ymdhms()
       << " end: " << _poly_end[sat].str_ymdhms() << endl;
#endif
    _gmutex.unlock(); return 1;
  }

  // prepare approximative estimates
  _PT.clear(); _X.clear(); _Y.clear(); _Z.clear(); 

  map<t_gtime,t_map_dat>::iterator itBeg = _mapsp3[sat].begin();
  map<t_gtime,t_map_dat>::iterator itEnd = _mapsp3[sat].end();
  map<t_gtime,t_map_dat>::iterator itReq = _mapsp3[sat].upper_bound(t);

  int dst = distance(itBeg, itEnd);                   // tab values [#]
  if( dst < fitdeg ){ 
     if(_log) _log->comment(1,"gallprec","no position available, few samples: "+int2str(dst) );
     _gmutex.unlock(); return (( _clknav && t_gallnav::nav( sat, t, xyz, var, vel, chk_mask ) >= 0 ) ? 1 : -1 ); 
  }
  if( dst < fitdat ){ fitdat = dst; }                 // shorten window

  int sign = 1;  // towards future
  double diffEnd = t-_poly_end[sat];
  double diffBeg = t-_poly_beg[sat];
  if( _poly_beg[sat] == LAST_TIME ){ itReq = _mapsp3[sat].upper_bound(t); --itReq; } // towards future, initialize
  else if( diffEnd > 0 && diffEnd < +900*fitdat ){ sign =  1; itReq = _mapsp3[sat].lower_bound(_poly_end[sat]); } // towards future
  else if( diffBeg < 0 && diffBeg > -900*fitdat ){ sign = -1; itReq = _mapsp3[sat].lower_bound(_poly_beg[sat]); } // towards past

  if(      sign > 0 && distance(itReq, itEnd) <= fitdat ){ itReq=itEnd; advance(itReq, -fitdat-1); } // towards future, but shift
  else if( sign < 0 && distance(itBeg, itReq) <= fitdat ){ itReq=itBeg; advance(itReq, +fitdat  ); } // towards past,   but shift

  _poly_beg[sat] = itReq->first; advance(itReq, +fitdat);
  _poly_end[sat] = itReq->first; advance(itReq, -fitdat);
   
  _sec = _poly_end[sat] - _poly_beg[sat];
  _ref = _poly_beg[sat] + _sec/2;

#ifdef DEBUG
  cout << " PRN: " << sat
       << " epo: " << t.str_ymdhms()
       << " req: " << itReq->first.str_ymdhms()
       << " ref: " << _ref.str_ymdhms()
       << " sec: " << _sec
       << " beg: " << _poly_beg[sat].str_ymdhms()
       << " end: " << _poly_end[sat].str_ymdhms() << setprecision(0)
       << " dat: " << fitdat << endl;
#endif

  while( _PT.size() < static_cast<unsigned int>(fitdat) )
  {
    ++itReq;
    t_gtime tt = itReq->first;
    _PT.push_back( tt.diff(_ref)/_sec );
     _X.push_back( _mapsp3[sat][tt]["X"] );
     _Y.push_back( _mapsp3[sat][tt]["Y"] );
     _Z.push_back( _mapsp3[sat][tt]["Z"] );

#ifdef DEBUG
     cout << fixed << setprecision(3)
          << " PRN " << sat << tt.str_ymdhms(" epo:") << " " <<  tt.diff(_ref)/_sec
          << " " <<  _mapsp3[sat][tt]["X"]
          << " " <<  _mapsp3[sat][tt]["Y"]
          << " " <<  _mapsp3[sat][tt]["Z"] 
          << " " <<   _X.size() << " " <<  _PT.size() << endl;
#endif
  }
   
  _poly_x[sat].fitpolynom( _PT, _X, fitdeg, _sec, _ref ); _poly_x[sat].evaluate( t.diff(_ref)/_sec, 0, xyz[0] );
  _poly_y[sat].fitpolynom( _PT, _Y, fitdeg, _sec, _ref ); _poly_y[sat].evaluate( t.diff(_ref)/_sec, 0, xyz[1] );
  _poly_z[sat].fitpolynom( _PT, _Z, fitdeg, _sec, _ref ); _poly_z[sat].evaluate( t.diff(_ref)/_sec, 0, xyz[2] );

//  _poly_x[sat].polynomials( _PT, _X ); _poly_x[sat].evaluate( t.diff(_ref), 0, xyz[0] );
//  _poly_y[sat].polynomials( _PT, _Y ); _poly_y[sat].evaluate( t.diff(_ref), 0, xyz[1] );
//  _poly_z[sat].polynomials( _PT, _Z ); _poly_z[sat].evaluate( t.diff(_ref), 0, xyz[2] );

#ifdef DEBUG   
  vector<double> v_xpoly = _poly_x[sat].polynomials();
  cout << " PRN: " << sat << fixed << setprecision(3) 
       << " req: " <<    t.str_ymdhms()
       << " ref: " << _ref.str_ymdhms()
       << " beg: " << _poly_beg[sat].str_ymdhms()
       << " end: " << _poly_end[sat].str_ymdhms() << setprecision(0)
       << " " << t.diff(_ref) << " " << xyz[0] << "  " << xyz[1] << "  " << xyz[2] << "  " << v_xpoly.size() << endl;
#endif

  _gmutex.unlock(); return 1;
}

set<string> t_gallprec::clk_objs()
{
	set<string> allobj;
	for (auto iter : _mapclk) {
		allobj.insert(iter.first);
	}
	return allobj;
}

set<t_gtime> t_gallprec::clk_epochs()
{
	set<t_gtime> alltime;
	for (auto iter : _mapclk) {
		for (auto time : iter.second) {
			alltime.insert(time.first);
		}
	}
	return alltime;
}

set<t_gallprec::clk_type> t_gallprec::get_clk_type() const
{
	return _clk_type_list;
}

// lagrange_interpolate
// -------------------
int t_gallprec::lagrange_pos(string sat, const t_gtime& t, t_gtriple& crd)
{
	gtrace("t_gallprec::_get_crddata");

	_T.clear();
	_PT.clear(); _X.clear(); _Y.clear(); _Z.clear();
	_CT.clear(); _C.clear();

	if (_mapsp3.find(sat) == _mapsp3.end()) return -1;

	map<t_gtime, t_map_dat>::iterator itReq = _mapsp3[sat].lower_bound(t); // 1st equal|greater [than t]

	if (itReq == _mapsp3[sat].end()) return -1;

	// add for nearest epoch
	map<t_gtime, t_map_dat>::iterator itReq_tmp = --(_mapsp3[sat].lower_bound(t));
	double sample = _intvm[sat];//yjqin
	if (abs(t.diff(itReq_tmp->first)) < abs(t.diff(itReq->first)))  itReq = itReq_tmp;

	//_ref = itReq->first; // get the nearest epoch after t as reference

	map<t_gtime, t_map_dat>::iterator itBeg = _mapsp3[sat].begin();
	map<t_gtime, t_map_dat>::iterator itEnd = _mapsp3[sat].end();
	map<t_gtime, t_map_dat>::iterator it = itReq;

	if (itReq == itEnd) {
		//    cerr << "itReq = " << itReq->first.str("%Y-%m-%d %H:%M:%S " <<
		//         << "itEnd = " << itEnd->first.str("%Y-%m-%d %H:%M:%S " << endl;
		if (_log && _log->verb() >= 1) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - no ephemeris found for "));
		return -1;
	}

#ifdef DEBUG
	cerr << "EPH FOUND: " << sat << " " << itReq->first.str("%Y-%m-%d %H:%M:%S[%T]")
		<< fixed
		<< setprecision(3)
		<< "  x= " << setw(13) << itReq->second["X"]
		<< "  y= " << setw(13) << itReq->second["Y"]
		<< "  z= " << setw(13) << itReq->second["Z"]
		<< setprecision(6)
		<< "  c= " << setw(13) << itReq->second["C"] * 1000000
		<< endl;
#endif

	// DISTANCE() NEED TO BE A POSITIVE DIFFERENCE !!!
	int limit = static_cast<int>(_degree_sp3 / 2); // round (floor)

	// too few data
	if (distance(itBeg, itEnd) < static_cast<int>(_degree_sp3)) {
		if (_log && _log->verb() >= 1) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - not enough eph data found for "));
		return -1;

		// start from the first item
	}
	else if (distance(itBeg, itReq) <= limit) {
		it = itBeg;

		// start from the last item
	}
	else if (distance(itReq, itEnd) <= static_cast<int>(_degree_sp3 - limit)) {
		it = itEnd;
		for (int i = 0; i <= static_cast<int>(_degree_sp3); i++) it--;

		// around requested item (standard case)
	}
	else {
		for (int i = 0; i < limit; i++) it--;
	}

	//Added by lewen
	if (it == itEnd) {
		if (_log && _log->verb() >= 1) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - no ephemeris found for "));
		return -1;
	}
	t_gtime t_beg = it->first;
	double epo_now = (t - t_beg) / sample + 1;
	if (epo_now > 9)//yjqin
	{
		if (_log) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - no ephemeris found for "));
		return -1;
	}
	// vector for polynomial
	for (unsigned int i = 0; i < _degree_sp3; it++, i++) {
		
		if (it->second["X"] != UNDEFVAL_POS) {

			_T.push_back(it->first);
			_X.push_back(it->second["X"]);
			_Y.push_back(it->second["Y"]);
			_Z.push_back(it->second["Z"]);
			_C.push_back(it->second["C"]);

#ifdef DEBUG
			cout << "EPH:" << i << " " << sat
				<< fixed
				<< setprecision(0)
				<< "  r=" << setw(6) << tdiff
				<< "  t=" << it->first.str("%Y-%m-%d %H:%M:%S ")
				<< setprecision(3)
				<< "  x=" << setw(13) << _X[i] // it->second["X"] 
				<< "  y=" << setw(13) << _Y[i] // it->second["Y"] 
				<< "  z=" << setw(13) << _Z[i] // it->second["Z"]
				<< setprecision(6)
				<< "  c=" << setw(13) << _C[i] * 1000000 // it->second["C"]*10000000
				<< " " << _X.size() << " " << _degree_sp3 * MAXDIFF_EPH
				<< endl;
#endif

		}
	}

	if (_X.size() != _degree_sp3) {
		if (_log && _log->verb() >= 1) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - not enough eph data for "));
		return -1;
	}


	vector<double> X;
	for (unsigned int i = 1; i <= _degree_sp3; i++) {
		X.push_back(i);
	}

	crd[0] = lagrange_interpolate(X, _X, epo_now);
	crd[1] = lagrange_interpolate(X, _Y, epo_now);
	crd[2] = lagrange_interpolate(X, _Z, epo_now);

	return 1;
}

int t_gallprec::lagrange_vel(string sat, const t_gtime& t, t_gtriple& crd)
{
	gtrace("t_gallprec::_get_crddata");

	_T.clear();
	_PT.clear(); _X.clear(); _Y.clear(); _Z.clear();
	_CT.clear(); _C.clear();

	if (_mapsp3.find(sat) == _mapsp3.end()) return -1;

	map<t_gtime, t_map_dat>::iterator itReq = _mapsp3[sat].lower_bound(t); // 1st equal|greater [than t]

	if (itReq == _mapsp3[sat].end()) return -1;

	// add for nearest epoch
	map<t_gtime, t_map_dat>::iterator itReq_tmp = --(_mapsp3[sat].lower_bound(t));
	double sample = itReq->first.diff(itReq_tmp->first);
	if (abs(t.diff(itReq_tmp->first)) < abs(t.diff(itReq->first)))  itReq = itReq_tmp;

	//_ref = itReq->first; // get the nearest epoch after t as reference

	map<t_gtime, t_map_dat>::iterator itBeg = _mapsp3[sat].begin();
	map<t_gtime, t_map_dat>::iterator itEnd = _mapsp3[sat].end();
	map<t_gtime, t_map_dat>::iterator it = itReq;

	if (itReq == itEnd) {
		//    cerr << "itReq = " << itReq->first.str("%Y-%m-%d %H:%M:%S " <<
		//         << "itEnd = " << itEnd->first.str("%Y-%m-%d %H:%M:%S " << endl;
		if (_log && _log->verb() >= 1) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - no ephemeris found for "));
		return -1;
	}

#ifdef DEBUG
	cerr << "EPH FOUND: " << sat << " " << itReq->first.str("%Y-%m-%d %H:%M:%S[%T]")
		<< fixed
		<< setprecision(3)
		<< "  Vx= " << setw(13) << itReq->second["VX"]
		<< "  Vy= " << setw(13) << itReq->second["VY"]
		<< "  Vz= " << setw(13) << itReq->second["VZ"]
		<< setprecision(6)
		<< "  c= " << setw(13) << itReq->second["C"] * 1000000
		<< endl;
#endif

	// DISTANCE() NEED TO BE A POSITIVE DIFFERENCE !!!
	int limit = static_cast<int>(_degree_sp3 / 2); // round (floor)

	// too few data
	if (distance(itBeg, itEnd) < static_cast<int>(_degree_sp3)) {
		if (_log && _log->verb() >= 1) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - not enough eph data found for "));
		return -1;

		// start from the first item
	}
	else if (distance(itBeg, itReq) <= limit) {
		it = itBeg;

		// start from the last item
	}
	else if (distance(itReq, itEnd) <= static_cast<int>(_degree_sp3 - limit)) {
		it = itEnd;
		for (int i = 0; i <= static_cast<int>(_degree_sp3); i++) it--;

		// around requested item (standard case)
	}
	else {
		for (int i = 0; i < limit; i++) it--;
	}

	//Added by lewen
	if (it == itEnd) {
		if (_log && _log->verb() >= 1) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - no ephemeris found for "));
		return -1;
	}
	t_gtime t_beg = it->first;
	double epo_now = (t - t_beg) / sample + 1;
	// vector for polynomial
	for (unsigned int i = 0; i < _degree_sp3; it++, i++) {

		if (it->second["X"] != UNDEFVAL_POS) {

			_T.push_back(it->first);
			_X.push_back(it->second["X"]);
			_Y.push_back(it->second["Y"]);
			_Z.push_back(it->second["Z"]);
			_C.push_back(it->second["C"]);

#ifdef DEBUG
			cout << "EPH:" << i << " " << sat
				<< fixed
				<< setprecision(0)
				<< "  r=" << setw(6) << tdiff
				<< "  t=" << it->first.str("%Y-%m-%d %H:%M:%S ")
				<< setprecision(3)
				<< "  x=" << setw(13) << _X[i] // it->second["X"] 
				<< "  y=" << setw(13) << _Y[i] // it->second["Y"] 
				<< "  z=" << setw(13) << _Z[i] // it->second["Z"]
				<< setprecision(6)
				<< "  c=" << setw(13) << _C[i] * 1000000 // it->second["C"]*10000000
				<< " " << _X.size() << " " << _degree_sp3 * MAXDIFF_EPH
				<< endl;
#endif

		}
	}

	if (_X.size() != _degree_sp3) {
		if (_log && _log->verb() >= 1) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - not enough eph data for "));
		return -1;
	}


	vector<double> X;
	for (unsigned int i = 1; i <= _degree_sp3; i++) {
		X.push_back(i);
	}
	crd[0] = lagrange_interpolate(X, _X, epo_now,true);
	crd[1] = lagrange_interpolate(X, _Y, epo_now,true);
	crd[2] = lagrange_interpolate(X, _Z, epo_now,true);
	return 1;
}
// find t_geph element
// ---------



shared_ptr<t_geph> t_gallprec::_find( string sat, const t_gtime& t )
{
  gtrace("t_gallprec::_find");

  if( _mapsp3[sat].size() == 0 ) return _null; // make_shared<t_geph>();

  // if not exists satellite not in cache
  t_map_sp3::iterator it = _prec.find(sat);
  if( it == _prec.end() ){
    if( _get_crddata( sat, t ) < 0 ) return _null; // make_shared<t_geph>();
//    cout << " CACHE INIT ! : " << sat << " " << t.str("%Y-%m-%d %H:%M:%S") << endl;
  }

  // could not find the data at all --- SHOULDN'T OCCURE, SINCE _get_crddata already return !
  it = _prec.find(sat);
  if( it == _prec.end() ){
    if( _log && _log->verb() >= 1 ) _log->comment(1,"gallprec",t.str_ymdhms(sat + " warning - gephprec element not found for "));
    return _null; //make_shared<t_geph>();
  }   

  double t_minus_ref = t - (it->second)->epoch();

//  cout   << " CACHE        : " << sat << " " << t.str("%Y-%m-%d %H:%M:%S") << endl;

  // standard case: cache - satellite found and cache still valid!
  if( fabs((float)t_minus_ref) < (it->second)->interval()/_degree_sp3 &&  // more frequently updated cache
      (it->second)->valid(t)                                              // internal gephprec validity
  ){
//    cout << " CACHE USED ! : " << sat << " " << t.str_ymdhms()
//         << " " << t_minus_ref << " < " << (it->second)->interval()/_degree_sp3 << endl;
    return it->second;

  // specific case: cache - satellite is not within its standard validity (try to update cache)
  }else{

    t_gtime beg(_mapsp3[sat].begin()->first);
    t_gtime end(_mapsp3[sat].rbegin()->first);
     
    // update cache only if not close to the prec data boundaries
    if( ( fabs(t.diff(beg)) > (it->second)->interval()/2 &&
          fabs(t.diff(end)) > (it->second)->interval()/2 ) ||
       ! (it->second)->valid(t)
    ){
/*       
      cout << " CACHE UPDATE : " << sat << " " << t.str_ymdhms()
//	   << " beg: " << beg.str_ymdhms()
//	   << " end: " << end.str_ymdhms()
	   << "   " << fabs(t.diff(beg))
	   << "   " << fabs(t.diff(end))
	   << " int1: " << (it->second)->interval()/2
	   << " int2: " << (it->second)->interval()/_degree_sp3
	   << " val:" << (it->second)->valid(t)
	   << " ref:" << t_minus_ref
	   << endl;
*/
      if( _get_crddata( sat, t ) < 0 ) return _null; //make_shared<t_geph>();
      it = _prec.find(sat);
      if( _log && _log->verb() >= 3 && _log->verb() >= 3 ) _log->comment(4,"gallprec",t.str_ymdhms(sat + " updated cache for "));
       
//    }else{
//      cout << " CACHE !! UPD : " << sat << " " << t.str("%Y-%m-%d %H:%M:%S") << endl;
    }
  }

  return it->second;

}


// fill PT,X,Y,Z vectors
// ----------
int t_gallprec::_get_crddata( string sat, const t_gtime& t )
{
   gtrace("t_gallprec::_get_crddata");   

  _T.clear();
  _PT.clear(); _X.clear(); _Y.clear(); _Z.clear(); 
  _CT.clear(); _C.clear();

  if( _mapsp3.find(sat) == _mapsp3.end() ) return -1;
   
  map<t_gtime,t_map_dat>::iterator itReq = _mapsp3[sat].lower_bound(t); // 1st equal|greater [than t]
   
  if( itReq == _mapsp3[sat].end() ) return -1;

  // add for nearest epoch add by glfeng
  map<t_gtime, t_map_dat>::iterator itReq_tmp = --(_mapsp3[sat].lower_bound(t));

  // jdhuang : fix the problem caused by glfeng add above
  // 别问，问就是不知道他为什么要这样改
  if (itReq_tmp != std::end(_mapsp3[sat]))
  {
	  if (abs(t.diff(itReq_tmp->first)) < abs(t.diff(itReq->first)))  itReq = itReq_tmp;
  }
   
  _ref = itReq->first; // get the nearest epoch after t as reference

  map<t_gtime,t_map_dat>::iterator itBeg = _mapsp3[sat].begin();
  map<t_gtime,t_map_dat>::iterator itEnd = _mapsp3[sat].end();
  map<t_gtime,t_map_dat>::iterator it    = itReq;

  if( itReq == itEnd ){
//    cerr << "itReq = " << itReq->first.str("%Y-%m-%d %H:%M:%S " <<
//         << "itEnd = " << itEnd->first.str("%Y-%m-%d %H:%M:%S " << endl;
    if( _log && _log->verb() >= 1 ) _log->comment(1,"gallprec",t.str_ymdhms(sat + " warning - no ephemeris found for "));
    return -1;
  }

#ifdef DEBUG
    cerr << "EPH FOUND: " << sat << " " << itReq->first.str("%Y-%m-%d %H:%M:%S[%T]")
         << fixed 
         << setprecision(3)
         << "  x= "  << setw(13) << itReq->second["X"]
         << "  y= "  << setw(13) << itReq->second["Y"] 
         << "  z= "  << setw(13) << itReq->second["Z"]
         << setprecision(6)
         << "  c= "  << setw(13) << itReq->second["C"]*1000000
         << endl;
#endif

  // DISTANCE() NEED TO BE A POSITIVE DIFFERENCE !!!
  int limit = static_cast<int>(_degree_sp3/2); // round (floor)
   
  // too few data
  if( distance(itBeg,itEnd) < static_cast<int>(_degree_sp3) ){
    if( _log && _log->verb() >= 1 ) _log->comment(1,"gallprec",t.str_ymdhms(sat + " warning - not enough eph data found for "));
    return -1;

 // =============================  modified by glfeng  <= change to < /// add return -1
  // start from the first item
  }else if( distance(itBeg,itReq) < limit ){
    it = itBeg;
	//return -1; // add by glfeng // jdhuang : fix
// =============================  modified by glfeng  <= change to < /// add return -1
// start from the last item
  }else if( distance(itReq,itEnd) < static_cast<int>(_degree_sp3 - limit) ){
    it = itEnd;
	//return -1; // add by glfeng // jdhuang : fix
    for(int i=0; i <= static_cast<int>(_degree_sp3); i++) it--;
// =============================
  // around requested item (standard case)
  }else{ 
     for(int i = 0; i < limit; i++ ) it--;
  } 

  //Added by lewen
  if (it == itEnd){
	  if (_log && _log->verb() >= 1) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - no ephemeris found for "));
	  return -1;
  }

  // vector for polynomial
  for(unsigned int i = 0; i <= _degree_sp3; it++, i++ ){
    double tdiff = it->first - _ref;
    
    // check maximum interval allowed between reference and sta/end epochs
    if( fabs(tdiff) > static_cast<double>(_degree_sp3*MAXDIFF_EPH) ) continue;
       
    if( it->second["X"] != UNDEFVAL_POS ){
       
      _PT.push_back( tdiff );
       _T.push_back( it->first );
       _X.push_back( it->second["X"] );
       _Y.push_back( it->second["Y"] );
       _Z.push_back( it->second["Z"] );
      _CT.push_back( tdiff );
       _C.push_back( it->second["C"] );

#ifdef DEBUG
       cout << "EPH:" << i << " " << sat
            << fixed 
	    << setprecision(0)
            << "  r=" << setw(6) << tdiff
            << "  t=" << it->first.str("%Y-%m-%d %H:%M:%S ")
            << setprecision(3)
            << "  x=" << setw(13) << _X[i] // it->second["X"] 
            << "  y=" << setw(13) << _Y[i] // it->second["Y"] 
            << "  z=" << setw(13) << _Z[i] // it->second["Z"]
            << setprecision(6)
            << "  c=" << setw(13) << _C[i]*1000000 // it->second["C"]*10000000
	    << " " << _X.size() << " " << _degree_sp3*MAXDIFF_EPH
            << endl;
#endif

    }
  }

  if( _X.size() != _degree_sp3 + 1 ){
    if( _log && _log->verb() >= 1 ) _log->comment(1,"gallprec",t.str_ymdhms(sat + " warning - not enough eph data for "));
    return -1;
  }

  if( _prec.find(sat) != _prec.end() ){
    _prec[sat]->degree(_degree_sp3);
    _prec[sat]->add( sat, _T, _X, _Y, _Z, _C );
    if( _log && _log->verb() >= 3 ) _log->comment(3,"gallprec",_prec[sat]->epoch().str_ymdhms(sat + " updating cache for "));
  }else{
    shared_ptr<t_gephprec> tmp(new t_gephprec); if( _log ) tmp->glog(_log);
    tmp->degree(_degree_sp3);
    tmp->add( sat, _T, _X, _Y, _Z, _C );
     _prec[sat] = tmp;
    if( _log && _log->verb() >= 3 ) _log->comment(3,"gallprec",_prec[sat]->epoch().str_ymdhms(sat + " creating cache for "));
  }

  return 1;
}


// fill CT,C vectors
// ----------
int t_gallprec::_get_clkdata( string sat, const t_gtime& t )
{
   gtrace("t_gallprec::_get_clkdata");   
 

   _CT.clear(); _C.clear(); _IFCB_F3.clear();
   
  if( _mapclk.find(sat) == _mapclk.end() ) return -1;
  map<t_gtime,t_map_dat>::iterator itBeg  = _mapclk[sat].begin();
  map<t_gtime,t_map_dat>::iterator itEnd  = _mapclk[sat].end();
  map<t_gtime,t_map_dat>::iterator itReq  = _mapclk[sat].lower_bound(t); // 1st equal|greater [than t]   

  if( itReq == _mapclk[sat].end() ) return -1;    // too old products
  if( t < itBeg->first ) return -1;               // too new products

  _clkref = itReq->first;                         // get the nearest epoch after t as reference

  map<t_gtime,t_map_dat>::iterator it     = itReq;
//  map<t_gtime,t_map_dat>::iterator itPlus = itReq; itPlus++;
//  map<t_gtime,t_map_dat>::iterator itMin  = itReq; if( itReq != itBeg ) itMin--;

  if( itReq == itEnd ){
    if( _log && _log->verb() >= 1 ) _log->comment(1,"gallprec",t.str_ymdhms(sat + " warning - no clock found for "));
    return -1;
  }

#ifdef DEBUG
    cout << "CLK FOUND: " << itReq->first.str_ymdhms(" ")
         << " c1= "   << itReq->second["C0"]
         << " c2= "   << itReq->second["C1"] 
         << " c3= "   << itReq->second["C2"]
         << endl;     
#endif

/*
  double dist = 0.0, dist1 = 0.0, dist2 = 0.0; // seconds
  if( itPlus != itEnd ) dist1 = fabs(itReq->first - itPlus->first );
  if( itMin  != itBeg ) dist2 = fabs(itReq->first -  itMin->first );

  if( dist1 != dist2 ){
     cout << sat <<  t.str(" %Y-%m-%d %H:%M:%S ")
          << " dist1 = " << dist1
          << " dist2 = " << dist2 
          << " itReq: "  << itReq->first.str(" %Y-%m-%d %H:%M:%S ")
          << " itBeg: "  << itBeg->first.str(" %Y-%m-%d %H:%M:%S ")
          << " itMin: "  << itMin->first.str(" %Y-%m-%d %H:%M:%S ")
          << endl;
    if( _log ) _log->comment(1,"gallprec",t.str_ymdhms(sat + " warning - no clocks distance clear for "));
    return -1;
  }else{
    dist = dist1;
  }
*/  
//  cerr << sat <<  t.str(" %Y-%m-%d %H:%M:%S  distance = ") << dist << endl;
   
  unsigned int degree_clk = 1; // MOZNA DAT 3 degree polinom kvuli 1Hz datum (ale cekovat ktery 3)
   
//  if( dist > 30 ) degree_clk = 7; // use higher polynomial degree
   
  // DISTANCE() NEED TO BE A POSITIVE DIFFERENCE !!!
  int limit = static_cast<int>(degree_clk/2); // round (floor)

  // too few data
  // disatance(beg,end) == size () distance is too low! ZHJ
  //if( distance(itBeg,itEnd) < static_cast<int>(degree_clk) ){

  //int distance_left = distance(itBeg, itReq);
  //int distance_right = _mapclk[sat].size() - distance_left;
  bool flag_left = false;
  auto itleft = itReq;
  for (int i = 0; i <= limit; i++) {
	  itleft--;
	  if (itleft == _mapclk[sat].end()) {
		  flag_left = true;
		  break;
	  }
  }
  bool flag_right = false;
  auto itright = itReq;
  for (int i = 0; i < static_cast<int>(degree_clk - limit); i++) {
	  itright++;
	  if (itright == _mapclk[sat].end()) {
		  flag_right = true;
		  break;
	  }
  }

  if (_mapclk[sat].size() < static_cast<int>(degree_clk)) {
    if( _log && _log->verb() >= 1 ) _log->comment(1,"gallprec",t.str_ymdhms(sat + " warning - not enough clk data found for "));
    return -1;

  // start from the first item
  //}else if( distance(itBeg,itReq) <= limit ){  
  // change by ZHJ
  }else if (flag_left) {
    it = itBeg;
     
  // start from the last item
  //}else if( distance(itReq,itEnd) <= static_cast<int>(degree_clk - limit) ){
  // change by ZHJ
  }else if (flag_right) {
    it = itEnd;
    for(int i=0; i <= static_cast<int>(degree_clk); i++) it--;

  // around requested item (standard case)
  }else{ 
     for(int i = 0; i <= limit; i++ ) it--;
  } 

  // calculate
  if( it->second["C1"] < UNDEFVAL_CLK ){
    if( it->second["C2"] < UNDEFVAL_CLK ){
        // xiongyun:add for ifcb
        if (it->second.find("IFCB_F3") == it->second.end())
        {
            //      double tdiff = it->first - _clkref;  // seconds !!!!!!!!!!!!!!!!!!!! instead of _clkref should be ReqT = t !!
            double tdiff = it->first - t;  // seconds !!!!!!!!!!!!!!!!!!!! instead of _clkref should be ReqT = t !!
            _CT.push_back(tdiff);
            _C.push_back(it->second["C0"] + it->second["C1"] * tdiff + it->second["C2"] * tdiff * tdiff);
        }
        else
        {
            it++;
            double tdiff = it->first - t;  // seconds !!!!!!!!!!!!!!!!!!!! instead of _clkref should be ReqT = t !!
            _CT.push_back(tdiff);
            _C.push_back(it->second["C0"]);
            _IFCB_F3.push_back(it->second["IFCB_F3"]);//add by xiongyun
        }
//      cout << "POCITAM C0+C1+C2 " << it->second["C0"] << " " <<  it->second["C1"] << " " << tdiff << "\n";
      return 1;
    }else{
//      double tdiff = it->first - _clkref;  // seconds !!!!!!!!!!!!!!!!!!!! instead of _clkref should be ReqT = t !!
      double tdiff = it->first - t;  // seconds !!!!!!!!!!!!!!!!!!!! instead of _clkref should be ReqT = t !!
      _CT.push_back( tdiff );
      _C.push_back( it->second["C0"] + it->second["C1"]*tdiff);
//      cout << "POCITAM C0+C1 " << it->second["C0"] << " " <<  it->second["C1"] << " " << tdiff << "\n";
      return 1;
    }

  // interpolate
  }else{	   
   
    // vector for polynomial
    for(unsigned int i = 0; i <= degree_clk; it++, i++ ){
      double tdiff = it->first - _clkref;
    
      // check maximum interval allowed between reference and sta/end epochs
      if( fabs(tdiff) > static_cast<double>(degree_clk*MAXDIFF_CLK) ){
	continue;
      }
       
      if( it->second["C0"] != UNDEFVAL_CLK ){
        _CT.push_back( tdiff );
         _C.push_back( it->second["C0"] );
	 
#ifdef DEBUG
         cout << "CLK:" << i 
	      << " "    << _CT.size()
	      << " "    << _C.size()
	      << fixed  << setprecision(0)
              << "  r=" << setw(6) << tdiff
              << "  t=" << it->first.str("%Y-%m-%d %H:%M:%S ")
	      << scientific << setprecision(8)
              << " c0=" << setw(14) << it->second["C0"]
              << " c1=" << setw(14) << it->second["C1"] 
              << " c2=" << setw(14) << it->second["C2"]
	      << endl;
#endif    

      }
    }

    if( _C.size() != degree_clk + 1 ){
		
		//Use Liner Poly by ZHJ
		if (itReq != itBeg) {
			_C.clear();
			_CT.clear();
			// Add First
			auto itFirst = itReq;
			itFirst--;
			_C.push_back(itFirst->second["C0"]);
			_CT.push_back(itFirst->first - _clkref);

			// Add Second
			_C.push_back(itReq->second["C0"]);
			_CT.push_back(0.0);

			return 1;
		}
		
      if( _log && _log->verb() >= 1 ) _log->comment(1,"gallprec",t.str_ymdhms(sat + " warning - not enough clk data found for "));
      return -1;
    }
  }

  return 1;
}


int t_gallprec::_get_delta_pos_vel(const string& sat, const t_gtime& t)
{
	gtrace("t_gallprec::_get_delta_pos_vel");

	_TCorr.clear();
	_PTCorr.clear(); _XCorr.clear(); _YCorr.clear(); _ZCorr.clear();
	_CTCorr.clear(); _CCorr.clear();

	if (_mapsp3.find(sat) == _mapsp3.end()) return -1;

	map<t_gtime, t_map_dat>::iterator itReq = _mapsp3[sat].lower_bound(t); // 1st equal|greater [than t]


	// if itReq=end, so the itReq is end--.
	if (itReq == _mapsp3[sat].end()) { itReq--; }

	_ref = itReq->first; // get the nearest epoch after t as reference

	map<t_gtime, t_map_dat>::iterator itBeg = _mapsp3[sat].begin();
	map<t_gtime, t_map_dat>::iterator itEnd = _mapsp3[sat].end();
	map<t_gtime, t_map_dat>::iterator it = itReq;

	if (itReq == itEnd) {
		//    cerr << "itReq = " << itReq->first.str("%Y-%m-%d %H:%M:%S " <<
		//         << "itEnd = " << itEnd->first.str("%Y-%m-%d %H:%M:%S " << endl;
		if (_log && _log->verb() >= 1) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - no SRR Correction found for "));
		return -1;
	}

#ifdef DEBUG
	cerr << "EPH FOUND: " << sat << " " << itReq->first.str("%Y-%m-%d %H:%M:%S[%T]")
		<< fixed
		<< setprecision(3)
		<< "  x= " << setw(13) << itReq->second["DX"]
		<< "  y= " << setw(13) << itReq->second["DY"]
		<< "  z= " << setw(13) << itReq->second["DZ"]
		<< setprecision(6)
		<< "  c= " << setw(13) << itReq->second["DCLK"] * 1000000
		<< endl;
#endif

	// DISTANCE() NEED TO BE A POSITIVE DIFFERENCE !!!
	int limit = static_cast<int>(_degree_sp3 / 2); // round (floor)

												   // too few data
	if (distance(itBeg, itEnd) < static_cast<int>(_degree_sp3)) {
		if (_log && _log->verb() >= 1) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - not enough SSR data found for "));
		return -1;

		// =============================  modified by glfeng  <= change to < /// add return -1
		// start from the first item
	}
	else if (distance(itBeg, itReq) < limit) {
		it = itBeg;
		return -1; // add by glfeng
				   // =============================  modified by glfeng  <= change to < /// add return -1
				   // start from the last item
	}
	//else if (distance(itReq, itEnd) < static_cast<int>(_degree_sp3 - limit)) {
	//	it = itEnd;
	//	return -1; // add by glfeng
	//	for (int i = 0; i <= static_cast<int>(_degree_sp3); i++) it--;
	//	// =============================
	//	// around requested item (standard case)
	//}
	else {
		for (int i = 0; i < limit; i++) it--;
	}

	// vector for polynomial
	for (unsigned int i = 0; i <= _degree_sp3; it++, i++) {
		double tdiff = it->first - _ref;

		// check maximum interval allowed between reference and sta/end epochs
		if (fabs(tdiff) > static_cast<double>(_degree_sp3*MAXDIFF_EPH)) continue;

		if (it->second["DX"] != UNDEFVAL_POS) {

			_PTCorr.push_back(tdiff);
			_TCorr.push_back(it->first);
			_XCorr.push_back(it->second["DX"]);
			_YCorr.push_back(it->second["DY"]);
			_ZCorr.push_back(it->second["DZ"]);
			_CTCorr.push_back(tdiff);
			_CCorr.push_back(it->second["DCLK"]);

#ifdef DEBUG
			cout << "EPH:" << i << " " << sat
				<< fixed
				<< setprecision(0)
				<< "  r=" << setw(6) << tdiff
				<< "  t=" << it->first.str("%Y-%m-%d %H:%M:%S ")
				<< setprecision(3)
				<< "  x=" << setw(13) << _XCorr[i] // it->second["DX"] 
				<< "  y=" << setw(13) << _YCorr[i] // it->second["DY"] 
				<< "  z=" << setw(13) << _ZCorr[i] // it->second["DZ"]
				<< setprecision(6)
				<< "  c=" << setw(13) << _CCorr[i] * 1000000 // it->second["C"]*10000000
				<< " " << _X.size() << " " << _degree_sp3*MAXDIFF_EPH
				<< endl;
#endif

		}
	}

	if (_XCorr.size() != _degree_sp3 + 1) {
		if (_log && _log->verb() >= 1) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - not enough SSR data for "));
		return -1;
	}

	if (_prec.find(sat) != _prec.end()) {
		_prec[sat]->degree(_degree_sp3);
		_prec[sat]->add(sat, _TCorr, _XCorr, _YCorr, _ZCorr, _CCorr);
		if (_log && _log->verb() >= 3) _log->comment(3, "gallprec", _prec[sat]->epoch().str_ymdhms(sat + " updating cache for "));
	}
	else {
		shared_ptr<t_gephprec> tmp(new t_gephprec); if (_log) tmp->glog(_log);
		tmp->degree(_degree_sp3);
		tmp->add(sat, _TCorr, _XCorr, _YCorr, _ZCorr, _CCorr);
		_prec[sat] = tmp;
		if (_log && _log->verb() >= 3) _log->comment(3, "gallprec", _prec[sat]->epoch().str_ymdhms(sat + " creating cache for "));
	}

	return 1;
}

int t_gallprec::_get_delta_pos_vel(const string & sat, const t_gtime & t, int iod, t_gtime& tRef, t_map_dat& orbcorr)
{
	gtrace("t_gallprec::_get_delta_pos_vel");
#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif
	/*_gmutex.lock();*/
	if (_mapsp3.find(sat) == _mapsp3.end() || _mapsp3[sat].size() == 0) {
		/*_gmutex.unlock();*/ return -1;
	}

	map<t_gtime, t_map_dat>::iterator itLast = _mapsp3[sat].lower_bound(t); // 1st equal|greater [than t]

	map<t_gtime, t_map_dat>::iterator itPrev = --(_mapsp3[sat].lower_bound(t)); // 1st equal|greater [than t]

	if (itLast != _mapsp3[sat].end() && int(itLast->second["IOD"]) == iod)
	{
		tRef = itLast->first;
		orbcorr = itLast->second;
		/*_gmutex.unlock();*/ return 1;
	}
	else if (itPrev != _mapsp3[sat].end() && int(itPrev->second["IOD"]) == iod) {
		tRef = itPrev->first;
		orbcorr = itPrev->second;
		/*_gmutex.unlock();*/ return 1;
	}
	else
	{
		/*_gmutex.unlock(); */return -1;
	}
}

int t_gallprec::_get_delta_clk(const string& sat, const t_gtime& t)
{
	gtrace("t_gallprec::_get_delta_clk");

	_CTCorr.clear(); _CCorr.clear();

	if (_mapclk.find(sat) == _mapclk.end()) return -1;
	map<t_gtime, t_map_dat>::iterator itBeg = _mapclk[sat].begin();
	map<t_gtime, t_map_dat>::iterator itEnd = _mapclk[sat].end();
	map<t_gtime, t_map_dat>::iterator itReq = _mapclk[sat].lower_bound(t); // 1st equal|greater [than t]   

	if (itReq == _mapclk[sat].end()) return -1;    // too old products
	if (t < itBeg->first) return -1;               // too new products

	_clkref = itReq->first; // get the nearest epoch after t as reference

	map<t_gtime, t_map_dat>::iterator it = itReq;
	//  map<t_gtime,t_map_dat>::iterator itPlus = itReq; itPlus++;
	//  map<t_gtime,t_map_dat>::iterator itMin  = itReq; if( itReq != itBeg ) itMin--;

	if (itReq == itEnd) {
		if (_log && _log->verb() >= 1) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - no clock correction found for "));
		return -1;
	}

#ifdef DEBUG
	cout << "CLK FOUND: " << itReq->first.str_ymdhms(" ")
		<< " c1= " << itReq->second["DCLK"]
		<< " c2= " << itReq->second["DOTCLK"]
		<< " c3= " << itReq->second["DOTDOTCLK"]
		<< endl;
#endif

	unsigned int degree_clk = 1; // MOZNA DAT 3 degree polinom kvuli 1Hz datum (ale cekovat ktery 3)

								 //  if( dist > 30 ) degree_clk = 7; // use higher polynomial degree

								 // DISTANCE() NEED TO BE A POSITIVE DIFFERENCE !!!
	int limit = static_cast<int>(degree_clk / 2); // round (floor)

												  // too few data
												  // disatance(beg,end) == size () distance is too low! ZHJ
												  //if( distance(itBeg,itEnd) < static_cast<int>(degree_clk) ){

												  //int distance_left = distance(itBeg, itReq);
												  //int distance_right = _mapclk[sat].size() - distance_left;
	bool flag_left = false;
	auto itleft = itReq;
	for (int i = 0; i <= limit; i++) {
		itleft--;
		if (itleft == _mapclk[sat].end()) {
			flag_left = true;
			break;
		}
	}
	bool flag_right = false;
	auto itright = itReq;
	for (int i = 0; i < static_cast<int>(degree_clk - limit); i++) {
		itright++;
		if (itright == _mapclk[sat].end()) {
			flag_right = true;
			break;
		}
	}

	if (_mapclk[sat].size() < static_cast<int>(degree_clk)) {
		if (_log && _log->verb() >= 1) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - not enough clk data found for "));
		return -1;

		// start from the first item
		//}else if( distance(itBeg,itReq) <= limit ){  
		// change by ZHJ
	}
	else if (flag_left) {
		it = itBeg;

		// start from the last item
		//}else if( distance(itReq,itEnd) <= static_cast<int>(degree_clk - limit) ){
		// change by ZHJ
	}
	else if (flag_right) {
		it = itEnd;
		for (int i = 0; i <= static_cast<int>(degree_clk); i++) it--;

		// around requested item (standard case)
	}
	else {
		for (int i = 0; i <= limit; i++) it--;
	}

	// calculate
	if (it->second["DOTCLK"] < UNDEFVAL_CLK) {
		if (it->second["DOTDOTCLK"] < UNDEFVAL_CLK) {
			//      double tdiff = it->first - _clkref;  // seconds !!!!!!!!!!!!!!!!!!!! instead of _clkref should be ReqT = t !!
			double tdiff = it->first - t;  // seconds !!!!!!!!!!!!!!!!!!!! instead of _clkref should be ReqT = t !!
			_CTCorr.push_back(tdiff);
			_CCorr.push_back(it->second["DCLK"] + it->second["DOTCLK"] * tdiff + it->second["DOTDOTCLK"] * tdiff*tdiff);
			//      cout << "POCITAM C0+C1+C2 " << it->second["C0"] << " " <<  it->second["C1"] << " " << tdiff << "\n";
			return 1;
		}
		else {
			//      double tdiff = it->first - _clkref;  // seconds !!!!!!!!!!!!!!!!!!!! instead of _clkref should be ReqT = t !!
			double tdiff = it->first - t;  // seconds !!!!!!!!!!!!!!!!!!!! instead of _clkref should be ReqT = t !!
			_CTCorr.push_back(tdiff);
			_CCorr.push_back(it->second["DCLK"] + it->second["DOTCLK"] * tdiff);
			//      cout << "POCITAM C0+C1 " << it->second["C0"] << " " <<  it->second["C1"] << " " << tdiff << "\n";
			return 1;
		}

		// interpolate
	}
	else {

		// vector for polynomial
		for (unsigned int i = 0; i <= degree_clk; it++, i++) {
			double tdiff = it->first - _clkref;

			// check maximum interval allowed between reference and sta/end epochs
			if (fabs(tdiff) > static_cast<double>(degree_clk*MAXDIFF_CLK)) {
				continue;
			}

			if (it->second["DCLK"] != UNDEFVAL_CLK) {
				_CTCorr.push_back(tdiff);
				_CCorr.push_back(it->second["DCLK"]);

#ifdef DEBUG
				cout << "CLK:" << i
					<< " " << _CTCorr.size()
					<< " " << _CCorr.size()
					<< fixed << setprecision(0)
					<< "  r=" << setw(6) << tdiff
					<< "  t=" << it->first.str("%Y-%m-%d %H:%M:%S ")
					<< scientific << setprecision(8)
					<< " c0=" << setw(14) << it->second["DCLK"]
					<< " c1=" << setw(14) << it->second["DOTCLK"]
					<< " c2=" << setw(14) << it->second["DOTDOTCLK"]
					<< endl;
#endif    

			}
		}

		if (_CCorr.size() != degree_clk + 1) {
			if (_log && _log->verb() >= 1) _log->comment(1, "gallprec", t.str_ymdhms(sat + " warning - not enough clk correction found for "));
			return -1;
		}
	}

	return 1;
}

int t_gallprec::_get_delta_clk(const string & sat, const t_gtime & t, int iod, t_gtime& tRef, t_map_dat& clkcorr)
{
	gtrace("t_gallprec::_get_delta_pos_vel");

#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif
	//_gmutex.lock();

	if (_mapclk.find(sat) == _mapclk.end() || _mapclk[sat].size() == 0) {
		/*_gmutex.unlock();*/ return -1;
	}
	map<t_gtime, t_map_dat>::iterator itLast = _mapclk[sat].lower_bound(t); // 1st equal|greater [than t]

	map<t_gtime, t_map_dat>::iterator itPrev = --_mapclk[sat].lower_bound(t); // 1st equal|greater [than t]


	if (itLast != _mapclk[sat].end() && int(itLast->second["IOD"]) == iod)
	{
		tRef = itLast->first;
		clkcorr = itLast->second;
		/*_gmutex.unlock();*/ return 1;
	}
	else if (itPrev != _mapclk[sat].end() && int(itPrev->second["IOD"]) == iod) {
		tRef = itPrev->first;
		clkcorr=itPrev->second;
		/*_gmutex.unlock();*/ return 1;
	}
	else
	{
		/*_gmutex.unlock();*/ return -1;
	}
}

} // namespace
