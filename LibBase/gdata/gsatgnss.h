
#ifndef  GSATGNSS_H
#define  GSATGNSS_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  Purpose: statistical function (1D)
  Version: $ Rev: $

  2014-04-18 /PV: created
  2018-09-28 /JD: revised

-*/

#include <string>
#include <vector>

#include "newmat/newmat.h"
#include "gdata/gdata.h"
#include "gdata/gobsgnss.h"
#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gall/gallnav.h"

using namespace std;

namespace gnut {

class t_gsatgnss : public t_gdata {
   
 public:
   t_gsatgnss();
// t_gsatgnss(string site, string sat, const t_gtime& t);  // TOTO ODSTRANIT ?!
   t_gsatgnss(t_gobsgnss* obs);
   virtual ~t_gsatgnss();
   
   void addcrd( const t_gtriple& crd );    // add satellite position
   void addclk( double d );                // add satellite clocks at the transmision time
   void addecl(t_gtriple& sat, t_gtime& epoch); // determine wheather eclipsed or not   
   
   int  addprd(t_gallnav* gnav, bool msk_health = true );  // add satellite pos, clk and ecl ( // true to support INP:chk_health settings only
   
   void addele( double d );                // add satellite elevation                           
   void addazi( double d );                // add satellite azimuth                             
   void addrho( double d );                // add satellite rho-vector
   
   t_gtriple satcrd();                     // get satellite position                       
   double  clk();			   // get satellite clocks at the transmision time 
   double dclk();
   double  ele();			   // get satellite elevation
   double  ele_deg();			   // get satellite elevation [deg]
   double  azi();			   // get satellite azimuth                        
   double  rho();			   // get satellite rho-vector
   bool    ecl();                          // get eclipsing

   void clear();
   bool valid();
   
   t_gobsgnss* gobs;                       // pointer na gobsgnss

 private:
   virtual void _clear();
   virtual bool _valid() const;

   t_gtriple         _satcrd;  // satellite position (X,Y,Z)
   double            _clk;     // satellite clocks (precise, time of transmision) (meters)
   double            _dclk;    // satellite clocks drift(precise, time of transmision) (meters/s)
   double            _ele;     // satellite elevation
   double            _azi;     // satellite azimuth
   double            _rho;     // satellite-station geometrical distance
   bool              _eclipse; // eclipsed satellite
   t_gtime           _lastEcl; // last eclipse time
};

} // namespace

#endif
