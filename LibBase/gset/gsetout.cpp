/**
*
* @verbatim
	 (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

	 (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
		Research Institute of Geodesy, Topography and Cartography
		Ondrejov 244, 251 65, Czech Republic
  @endverbatim
*
* @file		gsetout.cpp
* @brief	implements output setting class
* @author   Jan Dousa
* @version	1.0.0
* @date		2012-10-23
*
*/

#include <iomanip>
#include <sstream>
#include <algorithm>

#include "gset/gsetout.h"
#include "gutils/gfileconv.h"

using namespace std;
using namespace pugi;

namespace gnut {


// Convertor for OUT formats
// ----------
OFMT t_gsetout::str2ofmt(const string& s)
{ 
  string tmp = s;
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
  if( tmp == "OUT"     )  return XXX_OUT;
  if( tmp == "TSA"     )  return TSA_OUT;   
  if( tmp == "SUM"     )  return SUM_OUT;
  if( tmp == "LOG"     )  return LOG_OUT;
  if( tmp == "XTR"     )  return XTR_OUT;
  if( tmp == "NAV"     )  return NAV_OUT;
  if( tmp == "XML"     )  return XML_OUT;
  if( tmp == "XQC"     )  return XQC_OUT;
  if( tmp == "PPP"     )  return PPP_OUT;
  if( tmp == "FLT"     )  return FLT_OUT;
  if( tmp == "SMT"     )  return SMT_OUT;
  if( tmp == "LSQ"     )  return LSQ_OUT;
  if( tmp == "RES"     )  return RES_OUT;
  if( tmp == "GRD"     )  return GRD_OUT;
  if( tmp == "SRF"     )  return SRF_OUT;
  if( tmp == "CFG"     )  return CFG_OUT;
  if( tmp == "FIT"     )  return FIT_OUT;
  if( tmp == "RINEXN"  )  return RINEXN_OUT;
  if( tmp == "RINEXN2" )  return RINEXN2_OUT;
  if( tmp == "RINEXO"  )  return RINEXO_OUT;
  if( tmp == "RINEXC"  )  return RINEXC_OUT;
  if( tmp == "SINEX"   )  return SINEX_OUT;
  if( tmp == "TROSINEX")  return TROSINEX_OUT;
  if( tmp == "TROSINEX0") return TROSINEX0_OUT;
  if( tmp == "KML")       return KML_OUT;
  if( tmp == "PRDQC")     return PRDQC_OUT;
  if (tmp == "ENU")		  return ENU_OUT;
  if (tmp == "RESULT")    return RESULT_OUT;
  if (tmp == "PRE")           return PRE_OUT;
  if (tmp == "RB")          return  RB_OUT;
  if (tmp == "STA")           return STA_OUT;
  if (tmp == "KBROMC")           return KBROMC_OUT;   
  if (tmp == "INS")		 return INS_OUT;
  if (tmp == "PREOBS")		 return PREOBS_OUT;
  if (tmp == "INSKF")	 return INSKF_OUT;
  if (tmp == "INSKFPK")	 return INSKFPK_OUT;
  if (tmp == "ORB")      return ORB_OUT;
  if (tmp == "UPD")      return UPD_OUT;
  if (tmp == "SP3")		 return SP3_OUT;
  if (tmp == "ICS")      return ICS_OUT;
  if (tmp == "ICSLEO")   return ICSLEO_OUT;
  if (tmp == "ION")      return ION_OUT;
  if (tmp == "CAR")      return RECCLK_OUT;
  if (tmp == "CAS")	     return SATCLK_OUT;
  if (tmp == "RECCLK")	 return RECCLK_OUT;
  if (tmp == "RECCLK13") return RECCLK13_OUT;// add by xiongyun
  if (tmp == "SATCLK")	 return SATCLK_OUT;
  if (tmp == "SATCLK13") return SATCLK13_OUT;// add by xiongyun
  if (tmp == "AMBUPD")    return AMBUPD_OUT;
  if (tmp == "AMBFLAG_DIR") return AMBFLAG_DIR_OUT; // add by jqwu
  if (tmp == "ORBDIF")    return ORBDIF_OUT;
  if (tmp == "DIF")       return ORBDIF_OUT;
  if (tmp == "ORBFIT")    return ORBFIT_OUT;// added by yqyuan
  if (tmp == "AMBCON_LEO")return AMBCONLEO_OUT;
  if (tmp == "AMBCON")	  return AMBCON_OUT;
  if (tmp == "CON")   	  return AMBCON_OUT;
  if (tmp == "RECOVER")   return RECOVER_OUT;
  if (tmp == "RCV")      return RECOVER_OUT;
  if (tmp == "AUG")      return AUG_OUT;
  if (tmp == "TROP")     return TROP_OUT;
  if (tmp == "IMU")     return IMU_OUT;
  return OFMT(-1);
}

// Convertor for OUT formats
// ----------
string t_gsetout::ofmt2str(const OFMT& f)
{ 
  switch ( f ){        
   case XXX_OUT:     return "OUT";
   case TSA_OUT:     return "TSA";     
   case LOG_OUT:     return "LOG";
   case SUM_OUT:     return "SUM";
   case XTR_OUT:     return "XTR";
   case NAV_OUT:     return "NAV";
   case XML_OUT:     return "XML";
   case XQC_OUT:     return "XQC";
   case PPP_OUT:     return "PPP";
   case LSQ_OUT:     return "LSQ";
   case FLT_OUT:     return "FLT";
   case SMT_OUT:     return "SMT";
   case RES_OUT:     return "RES";
   case GRD_OUT:     return "GRD";
   case SRF_OUT:     return "SRF";
   case PRE_OUT:     return  "PRE";
   case PREOBS_OUT:     return "PREOBS";
   case RB_OUT:      return "RB";
   case KBROMC_OUT:   return "KBROMC";
   case STA_OUT:     return "STA";
   case FIT_OUT:     return "FIT";
   case CFG_OUT:     return "CFG";
   case RINEXN_OUT:  return "RINEXN";
   case RINEXN2_OUT: return "RINEXN2";
   case RINEXO_OUT:  return "RINEXO";
   case RINEXC_OUT:  return "RINEXC";
   case SINEX_OUT:   return "SINEX";
   case TROSINEX_OUT:return "TROSINEX";
   case TROSINEX0_OUT:return "TROSINEX0";
   case KML_OUT:     return "KML";     
   case INS_OUT:	 return "INS";
   case INSKF_OUT:	 return "INSKF";
   case INSKFPK_OUT: return "INSKFPK";
   case UPD_OUT:     return "UPD_OUT";
   case SP3_OUT:	 return "SP3";
   case RECCLK_OUT:	 return "RECCLK";
   case RECCLK13_OUT:return "RECCLK13";// add by xiongyun
   case SATCLK_OUT:	 return "SATCLK";
   case SATCLK13_OUT:return "SATCLK13";// add by xiongyun
   case AMBUPD_OUT:  return "AMBUPD";
   case RESULT_OUT:  return "RESULT";
   case RECOVER_OUT: return "RECOVER";
   case ICS_OUT:     return "ICS";
   case ICSLEO_OUT:  return "ICSLEO";
   case ION_OUT:     return "ION";
   case AUG_OUT:     return "AUG";
   case TROP_OUT:    return "TROP";
   case IMU_OUT:     return "IMU";
   default:          return "UNDEF";
  }
  return "UNDEF";
}


// Constructor
// ----------
t_gsetout::t_gsetout()
  : t_gsetbase(),
     _append(false),
     _verb(0),
     _upd(0),
     _len(0),
     _smp(0)
{
   _set.insert(XMLKEY_OUT);
}


// Destructor
// ----------
t_gsetout::~t_gsetout()
{}


// Get formats output size
// ----------
int t_gsetout::output_size(const string& fmt)
{
  _gmutex.lock();
   
  int tmp = _outputs(fmt).size();

  _gmutex.unlock(); return tmp;
}

/* ======================
 * REPLACED BY _UPD_LOG!!!
 * ======================
 *
// Read glog from settings
int t_gsetout::glog(t_glog* glog)
{
  if( !glog ) return -1;

  _log = glog;
  string logstr( this->outputs("log") );

  if( _log && !logstr.empty() ){
    _log->mask(logstr);
    _log->append( this->append() );
  }

  int verb = this->verb();
  if( verb > _log->verb() ) _log->verb( verb ); // only if higher verbosity

  return 0;
}
*/


// Return value
// ----------
int t_gsetout::verb()
{
  _gmutex.lock();

  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).attribute("verb").as_int();

  _gmutex.unlock(); return tmp;
}


// Retrun value
// -----------
bool t_gsetout::append()
{
  _gmutex.lock();
   
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).attribute("append").as_bool();

  _gmutex.unlock(); return tmp;
}

   
// Get string outputs
// ----------
string t_gsetout::outputs(const string& fmt)
{ 
  _gmutex.lock();
  
  string tmp = _outputs(fmt);

  _gmutex.unlock(); return tmp;
}

int t_gsetout::sp3_obslimit()
{
	_gmutex.lock();
	int obs;
	obs = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).attribute("obslimit").as_int();
	
	_gmutex.unlock(); return obs;
}
   
// get string output version
// ---------
string t_gsetout::version(const string& fmt)
{
  _gmutex.lock();
   
  string ver = DEFAULT_FILE_VER;
  xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).child(fmt.c_str());
  if( ! fmt.empty() &&
      ! node.attribute("ver").empty() )
  {
    ver = node.attribute("ver").as_string();
  }

  _gmutex.unlock(); return ver;
}

   
// get time offset [min] for the output file name
// ---------
int t_gsetout::file_toff(const string& fmt)
{
  _gmutex.lock();
   
  int upd = DEFAULT_FILE_OFF;
  xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).child(fmt.c_str());
  if( ! fmt.empty() &&
      ! node.attribute("toff").empty() )
  {
    upd = node.attribute("toff").as_int();
  }

  _gmutex.unlock(); return upd;
}


// set time system for the output file name
// ---------
t_gtime::t_tsys t_gsetout::file_tsys(const string& fmt)
{
  _gmutex.lock();

  t_gtime::t_tsys upd = t_gtime::str2tsys( DEFAULT_FILE_SYS );
  xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).child(fmt.c_str());
  if( ! fmt.empty() &&
      ! node.attribute("tsys").empty() )
  {
    upd = t_gtime::str2tsys( node.attribute("tsys").as_string() );
  }

  _gmutex.unlock(); return upd;
}


   
// get update [min] for output file update
// ---------
int t_gsetout::out_update(const string& fmt)
{
  _gmutex.lock();
   
  int upd = DEFAULT_FILE_UPD;
  xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).child(fmt.c_str());
  if( ! fmt.empty() &&
      ! node.attribute("upd").empty() )
  {
    upd = node.attribute("upd").as_int();
  }

  _gmutex.unlock(); return upd;
}

   
// get period [min] for output file content
// ----------
int t_gsetout::out_length(const string& fmt)
{
  _gmutex.lock();
   
  int len = DEFAULT_FILE_LEN;
  xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).child(fmt.c_str());

  if( ! fmt.empty() &&
      ! node.attribute("len").empty() )
  {
    len = node.attribute("len").as_int();
  }

  _gmutex.unlock(); return len;
}

   
// get sample [sec] for output file data sampling
// ---------
float t_gsetout::out_sample(const string& fmt)
{
  _gmutex.lock();
   
  float smp = DEFAULT_FILE_SMP;
  xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).child(fmt.c_str());
  if( ! fmt.empty() &&
      ! node.attribute("smp").empty() )
  {
    smp = node.attribute("smp").as_float();
  }
   
  _gmutex.unlock(); return smp;
}

string t_gsetout::output_port(const string & fmt)
{
	string port("");
	xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).child(fmt.c_str());
	if (!fmt.empty() &&
		!node.attribute("port").empty())
	{
		port = node.attribute("port").as_string();
	}
	return port;
}

// Get output formats
// ----------
set<string> t_gsetout::oformats()
{
  return _oformats();
}


// Get output formats
// ----------
set<string> t_gsetout::_oformats()
{
  set<string> tmp;
  for( xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).first_child(); node; node = node.next_sibling() ){
    tmp.insert( node.name() );
  }
  return tmp;
}


// Get formats outputs
// ----------
string t_gsetout::_outputs(const string& fmt)
{
  string str;

  //// jdhuang for debug
  //for (xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).first_child(); node; node = node.next_sibling())
  //{
	 // cout << node.name() << endl;
  //}

  for( xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).first_child(); node; node = node.next_sibling() ){
    if( node.name() == fmt ){
      istringstream is( node.child_value() );
      while( is >> str && !is.fail() ){
        if( str.find("://") == string::npos ) str = GFILE_PREFIX + str;
        return str;
      }       
    }
  }
  return "";
}

// settings check
// ----------
void t_gsetout::check()
{
  this->_upd_glog(); // update log!
   
  _gmutex.lock();
   
  // check existence of nodes/attributes
  xml_node parent = _doc.child(XMLKEY_ROOT);
  xml_node node = _default_node(parent, XMLKEY_OUT);

  // check existence of attributes
  _default_attr(node,"append", _append);

  // check supported input formats (see OFMT enum !)
  set<string> ofmt = _oformats();
  set<string>::const_iterator itFMT = ofmt.begin();
  while( itFMT != ofmt.end() ){
    string fmt = *itFMT;
    OFMT  ofmt = str2ofmt( fmt );
    if( ofmt < 0 ){
      _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).remove_child( node.child(fmt.c_str()) );
      if( _log ) _log->comment(0, "Warning: " + fmt + " out format not implemented [gsetout::check()]!");
      else                cerr << "Warning: " + fmt + " out format not implemented [gsetout::check()]!\n";
      itFMT++;
      continue;
    }
  
    // check application-specific output format
    if( _OFMT_supported.find( ofmt ) == _OFMT_supported.end() ){
      _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).remove_child( node.child(fmt.c_str()) );
      if( _log ) _log->comment(0, "Warning: " + fmt + " out format not supported by this application!");
      else                cerr << "Warning: " + fmt + " out format not supported by this application!\n";       
    }     
    itFMT++;            
  }  
   
  _gmutex.unlock(); return;
   
  this->_upd_glog(); // update log!   
}


// settings help
// ----------
void t_gsetout::help()
{
  _gmutex.lock();
       
  cerr << " <outputs append=\""  << _append << "\" verb=\"" << _verb << "\" >\n"
       << "   <flt> file://dir/name </flt>    \t\t <!-- filter output encoder -->\n"
       << " </outputs>\n";
   
  cerr << "\t<!-- outputs description:\n"
       << "\t <encoder> path </encoder>\n"
       << "\t ... \n"
       << "\t where path contains [file,tcp,ntrip]:// depending on the application\n"
       << "\t -->\n\n";

   _gmutex.unlock(); return;
}


// Read glog from settings
void t_gsetout::_upd_glog()
{
  if( _log ){ // only if exists
    string logstr( this->outputs("log") );

    if( !logstr.empty() ){
      _log->mask(logstr);
      _log->append( this->append() );
      _log->verb( this->verb() );
    }
  //  int verb = this->verb();
  //  if( verb > _log->verb() ) _log->verb( verb ); // only if higher verbosity
  }   
}


} // namespace
