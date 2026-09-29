
#ifndef GDATA_H
#define GDATA_H 
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  Purpose: 
  Version: $Rev:$

  2011-03-25 /JD: created

-*/

#include <sstream>
#include <iostream>

#include "gdata/gmonit.h"
#include "gutils/gmutex.h"
#include "gio/glog.h"
#include "gio/gnote.h"
#include "gutils/gcommon.h"

using namespace std;

namespace gnut {

class BASE_LIBRARY_EXPORT t_gdata : public t_gmonit {

 public:
  t_gdata();
  t_gdata( const t_gdata& data );
  virtual ~t_gdata();

  t_gdata& operator=( const t_gdata& data );
     
  enum ID_TYPE {
    NONE,    ///< = 0,  none
    OBJ,     ///< = 1,  object
    TRN,     ///< = 2,  transmitter
    REC,     ///< = 3,  receiver
	REC_LEO, ///		receiver in LEO
    FIL,     ///< = 4,  file
       
    OBS,     ///< = 10, obseravation base
    OBSGNSS, ///< = 11, gnss observations
    SATDATA, ///< = 12, gnss observations + satellite data
	       
    QCDATA,  ///< = xx, data quality control
    QCPROD,  ///< = xx, prod quality control
    
    EPH,     ///< = 20, navigation base
    EPHGPS,  ///< = 21, navigation
    EPHGLO,  ///< = 22, navigation
    EPHGAL,  ///< = 23, navigation
    EPHQZS,  ///< = 24, navigation
    EPHBDS,  ///< = 24, navigation       
    EPHSBS,  ///< = 25, navigation 
    EPHIRN,  ///< = 26, navigation        
    EPHPREC, ///< = 27, sp3/clocks
    EPHRTCM, ///< = 28, navigation + RTCM

    GRID,    ///< = xx, regular data grid
    GEOBASE, ///< = xx,   geoid data grid
    NWMBASE, ///< = xx, surface data grid
    NWMSURF, ///< = xx, surface data grid
    NWMPROF, ///< = xx, profile data grid

	RESOBS,
	RESPAR,

    ALLGIO,      ///< =   , all files
    ALLNAV,      ///< = 28, all navigation all
    ALLPREC,     ///< = 29, all sp3 + rinexc
	ALLORB,
    ALLRTCM,     ///< = 30, all rtcm + nav
    ALLOBS,      ///< = 31, all observations
    ALLOBJ,      ///< = 32, all objects
    ALLPCV,      ///< = 33, all PCV
    ALLOTL,      ///< = 34, all OTL
    ALLSURF,     ///< = xx, all NWM SURF
    ALLPROF,     ///< = xx, all NWM PROF
    ALLPROD,     ///< = xx, all PROD
    ALLBIAS,     ///< = xx, all PROD    
	ALLPANNEL,   ///<
	ALLSOLAR,	 ///<
	ALLPOLEUT1,	 ///<
	AllATTITUDE, ///<
	ALLRECOVER,  ///<
	ALLALBEDO,   ///<
	ALLLRA,
	ALLKBR,
	ALLEOP,
	ALLRB,
	ALLSTA,
	ALLNPT,
	ALLPREOBS,
	ALLPREPROC,

    STRBUFF, ///< = xx, generic product (ASCII string) for strbuff encoder
    POS,     ///< = 35, XYZT position/time
    POST,    ///< = 36, SP3 + CLOCKS (satellite position)    
    MET,     ///< = xx, meteorological parameters
    TRP,     ///< = 37, tropospheric parameters (ztd+gradients)
    TRPSLT,  ///< = 38, tropospheric parameters (slant delays)
    CLK,     ///< = 37, clocks
    ION,     ///< = xx, ionospheric parameters
	IONEX,   ///< = xx, ionospheric delay from tec grid products (GIM)

    PCV,     ///< = 40, PCV model
    OTL,     ///< = 41, ocean loading model
    BIAS,    ///< = 42, code & phase biases
    ERP,     ///< = 43, Earth orientation model

	
    SOL,			///< = 50  solution
	IMUDATA,
	LCI_POS,
	CAMDATA,
	LIDARDATA,
	UPD,            // add for upd
	AMBFLAG,
	AMBFLAG13,
	AMBFLAG14,
	AMBFLAG15,
	AMBUPD,
	//UPD_EPOCH,
	IFCB,
	AMBINP,			// ambinp file for ambfix
	AMBCON,			// amb constraint

	ALLDE,			// all planeteph
	ALLICS,
	ALLICSLEO,
	ALLSUM,
	ALLRESULT,

	AUG,
	
	EGM,			// EGM20008 file
	LEAPSECOND,
	OCEANTIDE,
	DESAISCOPOLECOEF,
	ORB,			// orb file
	SATPARS, 		// sat_parameters_new file    
	NPT,
	SLRF,
	KBROMC,
	LAST
  };

  enum ID_GROUP {
    GRP_NONE,    ///< = 0,  // none
    GRP_OBSERV,  ///< = 10, // observations
    GRP_EPHEM,   ///< = 20, // ephemerides
    GRP_PRODUCT, ///< = 30, // positions/products
    GRP_MODEL,   ///< = 40, // models
    GRP_SOLUT,   ///< = 50, // solutions
    GRP_OBJECT,  ///< = 60, // objects       
    GRP_GRID,    ///< = 70, // grid data
    GRP_GIO,     ///< = 80, // gio
    GRP_QC,      ///< = 90, // quality control
	GRP_IMU,
    GRP_LAST
  };

    virtual string show(int verb);
//  virtual void show(ostringstream& os, int verb);
//  virtual void show( int verb, int repeat );
//  virtual string show( int verb, int repeat, string& buff );
  
  void glog(t_glog* l){ _log = l; }                 // set/get glog pointer
  t_glog* glog()const{ return _log; }

  void gnote(t_gallnote* n){ _gnote = n; }         // set/get gallnote pointer
  t_gallnote* gnote()const{ return _gnote; }

  ID_TYPE  id_type()const{  return _type; }
  ID_GROUP id_group()const{ return _group; }

  string  str_type()const;
  string  str_group()const;

  static string type2str(ID_TYPE type);
      
//  void lock()   const { this->_mutex_data.lock();   }
//  void unlock() const { this->_mutex_data.unlock(); }
//  void lock()   const { this->_mutex.lock();   }  // OK for BOOST !!!
//  void unlock() const { this->_mutex.unlock(); }  // OK for BOOST !!!

  void lock()   const { this->_gmutex.lock();   }
  void unlock() const { this->_gmutex.unlock(); }

 protected:
  int id_type(  ID_TYPE  t );
  int id_group( ID_GROUP g );

//  mutable Mutex  _mutex_data;
//  mutable boost::mutex   _mutex;  // OK for BOOST
  mutable t_gmutex  _gmutex;
  t_glog*           _log;
  t_gallnote*       _gnote;
  ID_TYPE           _type;
  ID_GROUP          _group;
   
 private:

};

} // namespace

#endif
