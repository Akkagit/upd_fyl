
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  This file is part of the G-Nut C++ library.
 
-*/

#include <iostream>

#include "gutils/gcommon.h"
#include "gutils/gtimesync.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {


// sampling synchronization filter for epochs (return true if the epoch fits sampling)
// ----------
bool time_sync(const t_gtime& epo, double smp, double scl, t_glog* glog)
{
  gtrace("time_sync");

  // OLD for >=1Hz only
  // if( int(round(epoch.sod()+epoch.dsec()))%int(sampl) != 0) return false;

  // do filtering with sampling rate
  //    but still consider ROUNDING e.g. possible CLOCK DRIFT cases!
  if( smp > 0 )
  {
    // and||sampl >=1Hz (high-rate) requested, assume clock drift well-below the range of decimals of sampling rate
    if( scl > 0 && smp <= 1 ){
      int smpl = (int)round( scl * smp);
      int iepo = (int)(round(epo.sod() *scl
			                     + epo.dsec()*scl)); // synced to .0 day i.e >=1Hz suggests .0 day epochs at least!

      if( smpl == 0 ) return false; // to be save if mixed high/low-rate sampling occurs
      int resi = iepo%smpl;

#ifdef DEBUG
      if( glog ) glog->comment(1,"epo_sync",epo.str_ymdhms("high-rate data: ")+dbl2str(epo.dsec())
                                            + " i:" + int2str(iepo) + " s:" + int2str(smpl)
			                    + " m:" + int2str(resi) + " d:" + int2str(scl));
#endif
      if( resi != 0 ) return false;
    }

    // inp && sampl < 1Hz, i.e. low-rate data only!
    else  if( int(round(epo.sod()+epo.dsec()))%int(smp) != 0 ){
#ifdef DEBUG   
        if( glog ) glog->comment(1,"epo_sync",epo.str_ymdhms("low-rate data: ")+dbl2str(epo.dsec()));
#endif
      return false;
    }

#ifdef DEBUG
    if( glog ) glog->comment(1,"epo_sync",epo.str_ymdhms("ok, epoch fits sampling: ") + dbl2str(epo.dsec()));
#endif
  }
  return true;
}


} // namespace
