/*
 * rtcm_base.cpp
 *
 *  Created on: Jan 5, 2015
 *      Author: jakub
 */

#include <cmath>
#include "gcoders/rtcm_base.h"
#include "gutils/gtime.h"
using namespace gnut;

t_rtcm_base::t_rtcm_base()
{
  _buffer = NULL;
  _buffsize = NULL;
  _data = NULL;
  _log = NULL;
  _maxdteph = 2000.0;
  _maxdtobs = 1000.0;
  _tref = t_gtime();
  _trefptr = NULL;
}

void t_rtcm_base::rtcminit(char** targetbuffer, int* buffsize, map<string,t_gdata*>* datapointer, t_glog** logpointer)
{
  _buffer = targetbuffer;
  _buffsize = buffsize;
  _data = datapointer;
  _log = logpointer;

}

t_rtcm_base::~t_rtcm_base() {
  // TODO Auto-generated destructor stub
}

void t_rtcm_base::setmaxdteph(double maxdteph){
  _maxdteph = maxdteph;
}
double t_rtcm_base::getmaxdteph(){
  return _maxdteph;
}
void t_rtcm_base::setmaxdtobs(double maxdtobs){
  _maxdtobs = maxdtobs;
}
double t_rtcm_base::getmaxdtobs(){
  return _maxdtobs;
}

void t_rtcm_base::settref(t_gtime tref){
  _tref = tref;
  _trefptr = &_tref;
}
void t_rtcm_base::settrefrealtime(){
  _trefptr = NULL;
}
const t_gtime* t_rtcm_base::gettref() const{
  return _trefptr;
}

t_gtime t_rtcm_base::estimateabstime(t_gtime tref, double tepoch, double maxdelta, double epochsz, t_gtime tstart, t_gtime::t_tsys ts){
  if (!(2*maxdelta < epochsz)){
    //TODO: put warning into log that maxdelta has exceeded allowed limits and was set to epochsz/2 => no over-limit warnings shall be produced
	maxdelta = epochsz/2;
  }
  int32_t nepochs;
  t_gtime tmph; // temporary high (for epoch after tref)
  t_gtime tmpm; // temporary middle (for epoch with tref)
  t_gtime tmpl; // temporary low (for epoch before tref)

  nepochs = (int32_t)((tref + maxdelta - tstart)/epochsz);
  tmpm = tstart + tepoch + ((double)nepochs * epochsz);
  if (tmpm > (tref + maxdelta)){
	tmpl = tstart + tepoch + ((double)(nepochs - 1) * epochsz);
	if (tmpl < (tref - maxdelta)){ // there is no match within given maxdelta interval
	  //TODO: Raise over-limit warning in log
	  // finish looking for closest match outside the interval
	  tmph = tstart + tepoch + ((double)(nepochs + 1) * epochsz);
	  double deltah = abs(tref - tmph);
	  double deltam = abs(tref - tmpm);
	  double deltal = abs(tref - tmpl);
	  //pick the closest match and return it
	  if(deltal <= deltam && deltal <= deltah){
		return tmpl;
	  }else if(deltam <= deltal && deltam <= deltah){
	    return tmpm;
	  }else{
	    return tmph;
	  }
	}
	return tmpl;
  }
  return tmpm;
}

//complement time of unknown GPS epoch to most recent possible estimate (weekepoch: true=>epoch=1week; false=>epoch=1024weeks)
//maxdelta, tepoch... [sec]
t_gtime t_rtcm_base::estimateepochGPS(double tepoch, t_gtime tref, double maxdelta, bool weekepoch){
  t_gtime tstart;
  tstart = t_gtime(1980, 1, 6, 0, 0, 0, 0.0, t_gtime::GPS);
  double epochsz = 604800; // seconds in week
  if(!weekepoch){ epochsz = (int64_t)1024 * 604800; }
  return estimateabstime(tref, tepoch, maxdelta, epochsz, tstart, t_gtime::GPS);
}
t_gtime t_rtcm_base::estimateepochGPS(double tepoch, t_gtime tref, bool weekepoch){
  double maxdt = (weekepoch)? _maxdtobs : _maxdteph;
  return estimateepochGPS(tepoch, tref, maxdt, weekepoch);
}
t_gtime t_rtcm_base::estimateepochGPS(double tepoch, double maxdelta, bool weekepoch){
  t_gtime tref;
  if(_trefptr != NULL){
    tref = t_gtime(_tref);
    tref.tsys(t_gtime::GPS);
  }else{
    tref = t_gtime(t_gtime::GPS); //sets the current time
  }
  return estimateepochGPS(tepoch, tref, maxdelta, weekepoch);
}
t_gtime t_rtcm_base::estimateepochGPS(double tepoch, bool weekepoch){
  double maxdt = (weekepoch)? _maxdtobs : _maxdteph;
  return estimateepochGPS(tepoch, maxdt, weekepoch);
}
t_gtime t_rtcm_base::estimateepochGPShour(double tepoch){
  double maxdt = _maxdtobs;
  t_gtime tref;
  if(_trefptr != NULL){
    tref = t_gtime(_tref);
    tref.tsys(t_gtime::GPS);
  }else{
    tref = t_gtime(t_gtime::GPS); //sets the current time
  }
  t_gtime tstart;
  tstart = t_gtime(1980, 1, 6, 0, 0, 0, 0.0, t_gtime::GPS);
  double epochsz = 3600.0; // hour
  return estimateabstime(tref, tepoch, maxdt, epochsz, tstart, t_gtime::GPS);
}

//complement time of unknown GLONASS epoch to most recent possible estimate (weekepoch: true=>epoch=1week; false=>epoch=4years)
//maxdelta, tepoch... [sec]
t_gtime t_rtcm_base::estimateepochGLO(double tepoch, t_gtime tref, double maxdelta, bool weekepoch){
  t_gtime tstart;//TODO: when GLO time is implemented replace all UTC values for GLO and remove the + 3h offsets
  tstart = t_gtime(1996, 1, 1, 0, 0, 0, 0.0, t_gtime::UTC);
  double epochsz = 604800; // seconds in week
  if(!weekepoch){
    epochsz = (int64_t)1461 * 86400;
  }else{
    //NOTE: GLONASS day of week in RTCM 3.2 standard starts with Sunday however GLONASS start date is a Monday, see DF416
    tstart = tstart - 86400.0;//move GLONASS start date to Sunday closest to start date to match weekday numbering
  }
  return estimateabstime(tref, tepoch, maxdelta, epochsz, tstart, t_gtime::UTC) + 10800.0;// + 3h offset
}
t_gtime t_rtcm_base::estimateepochGLO(double tepoch, t_gtime tref, bool weekepoch){
  double maxdt = (weekepoch)? _maxdtobs : _maxdteph;
  return estimateepochGLO(tepoch, tref, maxdt, weekepoch);
}
t_gtime t_rtcm_base::estimateepochGLO(double tepoch, double maxdelta, bool weekepoch){
  t_gtime tref;
  if(_trefptr != NULL){
    tref = t_gtime(_tref);
    tref.tsys(t_gtime::UTC);
  }else{
    tref = t_gtime(t_gtime::UTC); //sets the current time
  }
  return estimateepochGLO(tepoch, tref, maxdelta, weekepoch);
}
t_gtime t_rtcm_base::estimateepochGLO(double tepoch, bool weekepoch){
  double maxdt = (weekepoch)? _maxdtobs : _maxdteph;
  return estimateepochGLO(tepoch, maxdt, weekepoch);
}
t_gtime t_rtcm_base::estimateepochGLOday(double tepoch){
  double maxdt = _maxdtobs;
  t_gtime tref;
  if(_trefptr != NULL){
    tref = t_gtime(_tref);
    tref.tsys(t_gtime::UTC);
  }else{
    tref = t_gtime(t_gtime::UTC); //sets the current time
  }
  t_gtime tstart;//TODO: when GLO time is implemented replace all UTC values for GLO and remove the + 3h offsets
  tstart = t_gtime(1996, 1, 1, 0, 0, 0, 0.0, t_gtime::UTC);
  double epochsz = 86400.0; // one day
  return estimateabstime(tref, tepoch, maxdt, epochsz, tstart, t_gtime::UTC) + 10800.0;// + 3h offset
}
t_gtime t_rtcm_base::estimateepochGLOhour(double tepoch){
  double maxdt = _maxdtobs;
  t_gtime tref;
  if(_trefptr != NULL){
    tref = t_gtime(_tref);
    tref.tsys(t_gtime::UTC);
  }else{
    tref = t_gtime(t_gtime::UTC); //sets the current time
  }
  t_gtime tstart;//TODO: when GLO time is implemented replace all UTC values for GLO and remove the + 3h offsets
  tstart = t_gtime(1996, 1, 1, 0, 0, 0, 0.0, t_gtime::UTC);
  double epochsz = 3600.0; // hour
  return estimateabstime(tref, tepoch, maxdt, epochsz, tstart, t_gtime::UTC) + 10800.0;// + 3h offset
}

//complement time of unknown GAL epoch to most recent possible estimate (weekepoch: true=>epoch=1week; false=>epoch=4096weeks)
//maxdelta, tepoch... [sec]
t_gtime t_rtcm_base::estimateepochGAL(double tepoch, t_gtime tref, double maxdelta, bool weekepoch){
  t_gtime tstart; //TODO: replace GPS times with GAL when implemented
  tstart = t_gtime(1999, 8, 22, 0, 0, 0, 0.0, t_gtime::GPS);
  double epochsz = 604800; // seconds in week
  if(!weekepoch){ epochsz = 4096.0 * 604800; }
  return estimateabstime(tref, tepoch, maxdelta, epochsz, tstart, t_gtime::GPS);
}
t_gtime t_rtcm_base::estimateepochGAL(double tepoch, t_gtime tref, bool weekepoch){
  double maxdt = (weekepoch)? _maxdtobs : _maxdteph;
  return estimateepochGAL(tepoch, tref, maxdt, weekepoch);
}
t_gtime t_rtcm_base::estimateepochGAL(double tepoch, double maxdelta, bool weekepoch){
  t_gtime tref;
  if(_trefptr != NULL){
    tref = t_gtime(_tref);
    tref.tsys(t_gtime::GPS);
  }else{
    tref = t_gtime(t_gtime::GPS); //sets the current time
  }
  return estimateepochGAL(tepoch, tref, maxdelta, weekepoch);
}
t_gtime t_rtcm_base::estimateepochGAL(double tepoch, bool weekepoch){
  double maxdt = (weekepoch)? _maxdtobs : _maxdteph;
  return estimateepochGAL(tepoch, maxdt, weekepoch);
}

//complement time of unknown QZS epoch to most recent possible estimate (weekepoch: true=>epoch=1week; false=>epoch=1024weeks)
//maxdelta, tepoch... [sec]
t_gtime t_rtcm_base::estimateepochQZS(double tepoch, t_gtime tref, double maxdelta, bool weekepoch){
  t_gtime tstart;
  tstart = t_gtime(1980, 1, 6, 0, 0, 0, 0.0, t_gtime::GPS);
  double epochsz = 604800; // seconds in week
  if(!weekepoch){ epochsz = (int64_t)1024 * 604800; }
  return estimateabstime(tref, tepoch, maxdelta, epochsz, tstart, t_gtime::GPS);
}
t_gtime t_rtcm_base::estimateepochQZS(double tepoch, t_gtime tref, bool weekepoch){
  double maxdt = (weekepoch)? _maxdtobs : _maxdteph;
  return estimateepochQZS(tepoch, tref, maxdt, weekepoch);
}
t_gtime t_rtcm_base::estimateepochQZS(double tepoch, double maxdelta, bool weekepoch){
  t_gtime tref;
  if(_trefptr != NULL){
    tref = t_gtime(_tref);
    tref.tsys(t_gtime::GPS);
  }else{
    tref = t_gtime(t_gtime::GPS); //sets the current time
  }
  return estimateepochQZS(tepoch, tref, maxdelta, weekepoch);
}
t_gtime t_rtcm_base::estimateepochQZS(double tepoch, bool weekepoch){
  double maxdt = (weekepoch)? _maxdtobs : _maxdteph;
  return estimateepochQZS(tepoch, maxdt, weekepoch);
}

//complement time of unknown BDS epoch to most recent possible estimate (weekepoch: true=>epoch=1week; false=>epoch=1024weeks)
//maxdelta, tepoch... [sec] // TODO: remove the 14 leap second offset on returned value, when BDS is available in t_gtime
t_gtime t_rtcm_base::estimateepochBDS(double tepoch, t_gtime tref, double maxdelta, bool weekepoch){
  t_gtime tstart; //NOTE: BDS offset of 14 leap seconds between GPS a BDS start, update this when BDS is available in t_gtime
  tstart = t_gtime(2006, 1, 1, 0, 0, 0, 0.0, t_gtime::GPS);
  double epochsz = 604800; // seconds in week
  if(!weekepoch){ epochsz = (int64_t)1024 * 604800; }
  return estimateabstime(tref, tepoch, maxdelta, epochsz, tstart, t_gtime::GPS) - 14; //TODO: see note above...
}
t_gtime t_rtcm_base::estimateepochBDS(double tepoch, t_gtime tref, bool weekepoch){
  double maxdt = (weekepoch)? _maxdtobs : _maxdteph;
  return estimateepochBDS(tepoch, tref, maxdt, weekepoch);
}
t_gtime t_rtcm_base::estimateepochBDS(double tepoch, double maxdelta, bool weekepoch){
  t_gtime tref;
  if(_trefptr != NULL){
    tref = t_gtime(_tref);
    tref.tsys(t_gtime::GPS);
  }else{
    tref = t_gtime(t_gtime::GPS); //sets the current time
  }
  return estimateepochBDS(tepoch, tref, maxdelta, weekepoch);
}
t_gtime t_rtcm_base::estimateepochBDS(double tepoch, bool weekepoch){
  double maxdt = (weekepoch)? _maxdtobs : _maxdteph;
  return estimateepochBDS(tepoch, maxdt, weekepoch);
}

//binary magnitude bmag(number, exp) = number * 2^exp
double t_rtcm_base::bmag(double number, int exp){
  if(exp < -63 || exp > 63)return 0;
  if(exp < 0){
    return number / (double)((uint64_t)1 << (-exp));
  }else{
    return number * (double)((uint64_t)1 << (exp));
  }
}

//true if specified bit in the buffer is logical '1'
bool t_rtcm_base::getonebit(int bitnumber) {
  int b = bitnumber/8;//number of byte within buffer
  int n = bitnumber%8;//number of bit within byte (MSB => 0)
  return ((*_buffer)[b] & BITSELECTOR[n]) != 0;
}

//read number of bits from specified location in buffer as unsigned integer number
uint64_t t_rtcm_base::getuint(int startbit, int numberofbits) {
  uint64_t result;
  result = 0;
  for (int i = startbit; i<(startbit + numberofbits); i++) {
    int b = i/8;//number of byte within buffer
    int n = i%8;//number of bit within byte (MSB => 0)
    result = result << 1;
    if (((*_buffer)[b] & BITSELECTOR[n]) != 0) {
      result |= 1;
    }
  }
  return result;
}

//read number of bits from specified location in buffer as signed magnitude integer number
int64_t t_rtcm_base::getsint(int startbit, int numberofbits) {
  int64_t result = (int64_t)t_rtcm_base::getuint(startbit+1, numberofbits-1);
  if (((*_buffer)[startbit/8] & BITSELECTOR[startbit%8]) != 0) {
    result = -result;
  }
  return result;
}

//read number of bits from specified location in buffer as 2's complement integer number
int64_t t_rtcm_base::get2int(int startbit, int numberofbits) {
  int64_t result = (int64_t)t_rtcm_base::getuint(startbit+1, numberofbits-1);
  if (((*_buffer)[startbit/8] & BITSELECTOR[startbit%8]) != 0) {
    result = ((!0) << ((int64_t)numberofbits -1)) | result;
  }
  return result;
}

//interpret numberofbits LSB as 2's complement
int64_t t_rtcm_base::decode2int(uint64_t unumber, int numberofbits) {
  uint64_t filter = (1 << ((int64_t)numberofbits-1));
  int64_t result = (int64_t)(unumber ^ filter);
  if ((filter & unumber) != 0) {
    result =   ((!0) << ((int64_t)numberofbits -1))  | result;
  }else{
    result = (!((!0) << ((int64_t)numberofbits -1))) & result;
  }
  return result;
}



