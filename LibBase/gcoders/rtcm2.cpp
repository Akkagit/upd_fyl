/*
 * rtcm2.cpp
 *
 *  Created on: Jan 6, 2015
 *      Author: jakub
 */

#include <iostream>
#include <memory>

#include "gdata/gnav.h"
#include "gdata/geph.h"
#include "gdata/gobsgnss.h"
#include "gcoders/rtcm2.h"
#include "gall/gallobj.h"
#include "gall/gallnav.h"
#include "gall/gallobs.h"
#include "gutils/gtime.h"
using namespace gnut;

t_rtcm2::t_rtcm2(){
  for(int i = 0; i<64; i++){//all non-defined messages set to unsupported
    _messagetemplates[i] = &t_rtcm2::msguspt;
  }
  //...
  _messagetemplates[17] = &t_rtcm2::msg17;//GPS ephemeris
  _messagetemplates[18] = &t_rtcm2::msg18;//CarrierPhases uncorrected
  _messagetemplates[19] = &t_rtcm2::msg19;//PseudoRanges uncorrected
  //...

  clearstate();

}

t_rtcm2::~t_rtcm2() {
}
/*
void t_rtcm2::printbin(uint64_t data){
  std::cout << (((data&(1<<31))!=0)?1:0) << (((data&(1<<30))!=0)?1:0) << " ";
  for(int i=0;i<24;i++){std::cout << (((data&(1<<(29-i)))!=0)?1:0);}
  std::cout << " ";
  for(int i=0;i<6;i++){std::cout << (((data&(1<<(5-i)))!=0)?1:0);}
  //std::cout << " ";
  std::cout << std::endl;
}
*/
void t_rtcm2::clearstate(){
  _synced = false;
  _lastparity = 0;
  for(int i=0;i<33;i++){
    _words[i] = 0;
  }
  _obsepochdata.clear();
  _stationcnt.clear();

  t_rtcm_base::clearstate();
}

//worddata value contains following data MSB to LSB: word data (D1..D24); word parity (D25..D30). Returns -1 if fail
int t_rtcm2::getwordbits(int bitnumber, uint64_t & worddata){
  worddata = 0;
  int resultsz = 0;
  int i = bitnumber%8;
  while ( true ){
    if ( (bitnumber/8)+1 >= *_buffsize ) {std::cout << "word reading failed, not enough data" << std::endl;return -1;}//return if not enough data in the buffer
    while ( i < 6 ){//reading rolled 6bit-byte
      if ( resultsz >= 30 ) return resultsz;
      worddata = worddata << 1;
      worddata += ( ((*_buffer)[bitnumber/8] & BITSELECTOR[7-i]) != 0 ) ? 1 : 0;
      i++;
      resultsz++;
    }
    bitnumber += 8; //move up one byte
    i=0;
  }
  return 0;//unreachable code...
}

//read data field as an unsigned long from a single messages data-field; returns 0 when called with invalid arguments
uint64_t t_rtcm2::getnumber(int wordnumber, int startbitindex, int length){
  wordnumber--;
  if((length + startbitindex) > 24 || length < 0 || startbitindex < 0 || wordnumber > 32 || wordnumber < 0)return 0;
  uint64_t aword = ((_words[wordnumber] & 0x40000000) != 0) ? (_words[wordnumber] ^ 0x3FFFFFC0) : _words[wordnumber];
  uint64_t bitmask = 0;
  for (int i=0; i<length; i++){
    bitmask = bitmask<<1;
    bitmask |= 1;
  }
  bitmask = bitmask << (30-startbitindex-length);
  aword &= bitmask;
  return(aword >> (30-startbitindex-length));
}

//Called on a position in the buffer determines whether this bit is a start of a new rtcm_v2 message
//Returns:
//Not message start                       ...  < 0
//Suspected message start (missing data)  ... == 0
//Confirmed message start                 ...  > 0 (== number of bits in the message); message type is filled properly
int t_rtcm2::messagefound(int bitnumber, int & messagetype){
  //get the potential word and test if the message begins here
  uint64_t theword = 0;
  if ( getwordbits(bitnumber, theword) <= 0 ) {
    std::cout << "Not enough data in the buffer." << std::endl;
    return 0;
  }
  _words[0] = theword;
  if( getnumber(1,0,8) != 0x66 || (bitnumber%8) >= 6 ) {//check bits in preamble position for sequence 0110 0110
    //std::cout << bitnumber <<" Preamble mismatch." << std::endl;
    return -1;//does not match preamble
  } else {

	//possible message start detected...
    //std::cout << "Possible message start detected." << std::endl;

    //first header word
	if(!check_word_parity(theword))return -1;
	//_words[0] = theword;
	//std::cout << "Header B1 parity OK." << std::endl;
    messagetype = (int)getnumber(1,8,6);

    //second header word
    if( getwordbits(bitnumber+40, theword) <= 0 ) {
      //std::cout << "Not enough data in the buffer." << std::endl;
      return 0;
    }
    if(!check_word_parity(theword))return -1; //parity failed
    _words[1] = theword;
    //std::cout << "Header B2 parity OK." << std::endl;

    //print header words:
    /*std::cout << "Header B1: ";
    printbin(_words[0]);
    std::cout << "Header B2: ";
    printbin(_words[1]);
    */

    //data words
    int numberofdatawords = (int)getnumber(2,16,5);
    //std::cout << "Message length reads: "<< numberofdatawords << " + 2 words." << std::endl;

    for(int i = 0; i < numberofdatawords; i++){
      if( getwordbits(bitnumber + (40*(2 + i)), theword) <= 0 ) {
        //std::cout << "Not enough data in the buffer." << std::endl;
        return 0;
      }
      if(!check_word_parity(theword)){
        //std::cout << "Word #" << (i+3) << " parity Failed." << std::endl;
        return -1;
      } //parity failed
      _words[2 + i] = theword;
    }

    //message found, report and return
    //std::cout << "Complete message detected!" << std::endl;
    if ( messagetype == 0 ) { messagetype = 64; } //exception: binary 000000 means message type 64
    //if (messagetype >= 60)return -1;
    return (80+(40*numberofdatawords));
  }
}


//uses internal object data to determine if parity is valid, adjusts bits 31&30(MSb) to valid *D29 *D30 (previous parity)
bool t_rtcm2::check_word_parity (uint64_t & aword) {
  bool result = false;
  if(_synced){
    aword |= ((_lastparity & 0x03) << 30);
    unsigned char parity = count_parity(aword);
    //std::cout << "check: " << (aword & 0x3F) << " ?= " << (parity & 0x3F) << std::endl;
    result = ((parity & 0x3F) == (unsigned char)(aword & 0x3F));
  }else{
    result = resync_word_parity(aword);
  }
  _synced = result;
  _lastparity = aword & 0x03;
  return result;
}

bool t_rtcm2::resync_word_parity (uint64_t & aword) {
  for (uint64_t j = 0; j<4; j++){
	aword &= 0x3FFFFFFF;
    aword |= (j << 30);
    unsigned char parity = count_parity(aword);
    //std::cout << "check: " << (aword & 0x3F) << " ?= " << (parity & 0x3F) << std::endl;
    if((parity & 0x3F) == (unsigned char)(aword & 0x3F)) return true;
  }
  return false;
}

//aword contains following data MSB to LSB: (previous parity: *D29 *D30); word data (D1..D24); word parity (D25..D30).
unsigned char t_rtcm2::count_parity(uint64_t aword){
  aword = ((aword & 0x40000000) != 0) ? (aword ^ 0x3FFFFFC0) : aword;
  unsigned char parity = 0;
  for(int i = 0; i<6; i++){
    parity = parity << 1;
    uint64_t tmpw = aword & parityform[i];
    parity += ptable[(tmpw>>0) & 0xFF]^ptable[(tmpw>>8) & 0xFF]^ptable[(tmpw>>16) & 0xFF]^ptable[(tmpw>>24) & 0xFF];
  }
  return parity;
}

//########################################################################################
//######################  DECODING FUNCTION DECLARATIONS  ################################
//########################################################################################

int t_rtcm2::readmessage(){
  int messagetype = (int)getnumber(1,8,6);
  return (this->*_messagetemplates[messagetype])();
}

int t_rtcm2::msguspt(){
  int messagetype = (int)getnumber(1,8,6);
  std::cout << "Decoding of RTCM2 message #" << ((messagetype == 0)?64:messagetype) << " not supported!" << std::endl;
  return 0;
}

string t_rtcm2::resolvestation(uint64_t stationid){
  if(_stationcnt.find(stationid) != _stationcnt.end()){
    _stationcnt[stationid]++;
  }else{
    if(!_stationcnt.empty()){
      //TODO: generate warning, multiple station IDs in a single stream
    }
    _stationcnt[stationid] = 1;
  }
  char site[10] = "";//reference station cstring
  sprintf( site, "%04lluRTCM2", stationid);//NOTE: STRING DOES NOT MATCH RINEX IDs! RINEX has format of station_id: XXXXMRCCC, RTCM3 has number in range: 0..4095
  return string(site);
}

//NOTE: This message is broadcasted only when IODC does not match IODE in the most recent navigation message.
//      It contains the data from previous one. UNTESTED
int t_rtcm2::msg17(){//GPS Ephemeris
  int messagewords = (int)getnumber(2,16,5)+2;
  if (messagewords != 22){std::cout << "words:  22 != " << messagewords << " Invalid message! (length mismatch)" << std::endl;return 0;}

  t_gnavdata data = { 0 };
  t_gtime epoch;
  string prn;
  char prncstring[10];
  shared_ptr<t_gnav> geph = make_shared<t_gnav>();

  // ## READING INDIVIDUAL DATA FIELDS ##
  data[21] = (double)getnumber( 3, 0,10); // week
  data[19] = (double)decode2int(getnumber( 3,10,14),14) * G_PI; //I DOT
  data[ 3] = (double)getnumber( 4, 0, 8); //IODE
  long toc = bmag(getnumber( 4, 8,16),  4);//t_oc
  data[ 1] = bmag(decode2int(getnumber( 5, 0,16),16),-43);//a_f1; long has 64bits on every system so this should work
  data[ 2] = bmag(decode2int(getnumber( 5,16, 8), 8),-55);//a_f2
  data[ 4] = bmag(decode2int(getnumber( 6, 0,16),16), -5);//C_rs
  data[ 5] = bmag(decode2int((getnumber( 6,16, 8) <<  8) | (getnumber( 7, 0, 8)),16),-43) * G_PI;//delta n
  data[ 7] = bmag(decode2int(getnumber( 7, 8,16),16),-29);//C_uc
  data[ 8] = bmag((getnumber(8,0,24) << 8) | (getnumber(9,0,8)),-33);//eccentricity
  data[ 9] = bmag(decode2int(getnumber( 9, 8,16),16),-29);//C_us
  data[10] = bmag((getnumber(10, 0,24) << 8) | (getnumber(11, 0, 8)),-19);//A^(1/2)
  data[11] = bmag(getnumber(11, 8,16),  4);//t_oe
  data[13] = bmag(decode2int((getnumber(12, 0,24) <<  8) | (getnumber(13, 0, 8)),32),-31) * G_PI;//OMEGA_0
  data[12] = bmag(decode2int(getnumber(13, 8,16),16),-19);//C_ic
  data[15] = bmag(decode2int((getnumber(14, 0,24) <<  8) | (getnumber(15, 0, 8)),32),-31) * G_PI;//i_0
  data[14] = bmag(decode2int(getnumber(15, 8,16),16),-29);//C_is
  data[17] = bmag(decode2int((getnumber(16, 0,24) <<  8) | (getnumber(17, 0, 8)),32),-31) * G_PI;//omega
  data[16] = bmag(decode2int(getnumber(17, 8,16),16), -5);//C_rc
  data[18] = bmag(decode2int(getnumber(18, 0,24),24),-43) * G_PI;//OmegaDot
  data[ 6] = bmag(decode2int((getnumber(19, 0,24) <<  8) | (getnumber(20, 0, 8)),32),-31) * G_PI;//M_0
  data[26] = (double)getnumber(20, 8,10);//IODC
  data[ 0] = bmag(decode2int((getnumber(20,18, 6) << 16) | (getnumber(21, 0,16)),22),-31);//a_f0
  sprintf(prncstring,"%c%02i", 'G', (int)((getnumber(21,16, 5) == 0)? 32 : getnumber(21,16, 5)) );//PRN ID
  // data[FILL] = (double)getnumber(21,21,3);//FILL ...
  data[25] = bmag(decode2int(getnumber(22, 0, 8), 8),-31);//t_GD
  data[20] = (double)getnumber(22, 8, 2);//L2 code
  data[23] = (double)getnumber(22,10, 4);//SV accuracy
  data[24] = (double)getnumber(22,14, 6);//SV HEALTH
  data[22] = (double)getnumber(22,20, 1);//L2 P-code data flag
  // data[FILL] = (double)getnumber(22,21,3);//FILL ...
  //Time Of Transmission is left set to zero (not available from this message)

  //determine time, estimate the epoch
  epoch = estimateepochGPS(toc + (data[21]*604800.0));

  geph = make_shared<t_gnavgps>();
  geph->glog(*_log); epoch.tsys(t_gtime::GPS);
  epoch.from_gws(data[21], toc);
  geph->data2nav( (string)prncstring, epoch, data );

  // skip non-valid messages
  set<string> msg;
	 if( geph->chk(msg) < 0 ){
			if( *_log ) {
				 for(auto it = msg.begin(); it != msg.end(); it++){
						(*_log)->comment(1, "rtcm2", prn + epoch.str_ymdhms(" ") + " " + *it);
				 }
			}
			return 0;
	 }else{  
			if( *_log ) {
				 for(auto it = msg.begin(); it != msg.end(); it++){
						(*_log)->comment(1, "rtcm2", prn + epoch.str_ymdhms(" ") + " " + *it);
				 }
			}						 
	 }
	 
  // fill data
  if( *_log ) (*_log)->comment(2,"rtcm2","add sat nav " + (string)prncstring + epoch.str_ymdhms(" ") );
  map<string,t_gdata*>::iterator it = (*_data).begin();
  while( it != (*_data).end() ){
    if (it->second->id_type()  == t_gdata::ALLNAV ||
        it->second->id_group() == t_gdata::GRP_EPHEM){
      ((t_gallnav*)it->second)->add( geph );
    }
    it++;
  }

  return 1; //(return number of records read, -1 means error has occurred, 0 means message was skipped for some reason)
}

int t_rtcm2::passOBSd(){
  int cnt = 0;
  //NOTE: taken from rinexo2, passing data from single epoch of MSM observations to t_gall
  map<string,t_gdata*>::iterator itDAT = (*_data).begin();
  // adding channel canal to t_gobsgnss from gallobj
  // !!!! must be before adding obs data to t_gallobs
  while( itDAT != (*_data).end() ){
    if( itDAT->second->id_type() != t_gdata::ALLOBJ ){
      itDAT++;
      continue;
    }
    vector<shared_ptr<t_gobsgnss>>::iterator itOBS = _obsepochdata.begin();
    while( itOBS != _obsepochdata.end() ){
      t_gallobj* all_obj = (t_gallobj*)itDAT->second;
      shared_ptr<t_grec> null_obj, one_obj = dynamic_pointer_cast<t_grec>(all_obj->obj( (*itOBS)->sat() ));
      if( one_obj != null_obj && one_obj->id_type() == t_gdata::TRN ){
        int ch = one_obj->channel();
        (*itOBS)->channel(ch);
      }
      itOBS++;
    }
    itDAT++;
  }
  // fill t_gobsgnss
  itDAT = (*_data).begin();
  while( itDAT != (*_data).end() ){

    if( itDAT->second->id_type() != t_gdata::ALLOBS ){
      itDAT++;
      continue;
    }

    vector<shared_ptr<t_gobsgnss>>::const_iterator itOBS = _obsepochdata.begin();
    while( itOBS != _obsepochdata.end() ){

      if( (*_log) && (*_log)->verb() >= 4 )
          (*_log)->comment(1,"rtcm3",(*itOBS)->sat() + " obs filled "
                                  + (*itOBS)->epoch().str_ymdhms()
                            + " " + (*itOBS)->site());

      ((t_gallobs*)itDAT->second)->addobs( *itOBS );
      itOBS++;
    }
    cnt++;
    itDAT++;
  }
  _obsepochdata.clear();
  return cnt;
}

shared_ptr<t_gobsgnss> t_rtcm2::getgobsgnss(t_gtime epoch, uint64_t siteid, string sat){
  //Check if the stream has only a single station ID
  string site = resolvestation(siteid);

  //Search for observation object in the _legepochdata vector
  shared_ptr<t_gobsgnss> nullsptrgobs;
  shared_ptr<t_gobsgnss> resultsptr;
  vector<shared_ptr<t_gobsgnss>>::iterator itOBS = _obsepochdata.begin();
  while(itOBS != _obsepochdata.end()){
    if ((*itOBS)->epoch() == epoch){ //NOTE: this should always resolve as true
      if ((*itOBS)->site().compare(site) == 0){
        if ((*itOBS)->sat().compare(sat) == 0){
          resultsptr = (*itOBS);
        }
      }
    }else{
      //TODO: warning epoch time mismatch within epoch, ignore by default, add as new t_gobsgnss
    }
  }
  //if not found add newly constructed object
  if(resultsptr == nullsptrgobs){
    resultsptr.reset(new t_gobsgnss(site, sat, epoch));
    _obsepochdata.push_back(resultsptr);
  }
  return resultsptr;
}


int t_rtcm2::msg18(){//RTK Uncorrected Carrier Phases
  shared_ptr<t_gobsgnss> gobs;
  int cnt = 0;
  int epochend = 0;

  // ## READING INDIVIDUAL DATA FIELDS ##
  //common data
  int siteid = getnumber( 1, 14,10); //reference station number
  double zcnt = (double)getnumber( 2, 0,13) * 0.6; // Modified Z count [s]
  int messagewords = (int)getnumber(2,16,5)+2;
  int band = (int)getnumber(3,0,2);//Frequency... 0=L1, 1=L2, 2&3...Reserved
  getnumber(3,2,2);//Reserved
  double tus = (double)getnumber(3,4,20);//GNSS Time Of Measurement [us]
  //satellite specific data
  int wordsread = 3;
  while(wordsread < messagewords){
    epochend = getnumber(wordsread+1,0,1);//MULTIPLE MESAGE INDICATOR, if 1: more messages with identical time-tag follow
    int code = getnumber(wordsread+1,1,1);//Code INDICATOR: 0 = C/A Code; 1 = P-Code
    int satsys = getnumber(wordsread+1,2,1);//GNSS system: 0 = GPS; 1 = GLONASS
    int satslot = getnumber(wordsread+1,3,5);//satellite ID ( slot number GLONASS / PRN GPS ); if GPS 0=32
    getnumber(wordsread+1,8,3);//Data quality
    getnumber(wordsread+1,11,5);//Cumulative Loss Of Continuity Indicator
    double phcyls = (double)decode2int((( getnumber(wordsread+1,16,8) <<24 ) + getnumber(wordsread+2,0,24)),32) / 256.0; //CARRIER PHASE [cycles]
    //NOTE: while phase cycles can be easily converted to PhaseRange, the range of this message is smaller than that of the value
    //      the number of roll-overs is, according to the standard, up to a user to determine. +-8,388,608 cycles
    wordsread += 2;

    char prncstring[10];
    t_gtime epoch;
    if (satsys == 0){ //GPS
      sprintf(prncstring,"%c%02i", 'G', satslot);
      epoch = estimateepochGPShour(zcnt + (tus/1000000));
    }else{ //GLONASS
      sprintf(prncstring,"%c%02i", 'R', satslot);
      epoch = estimateepochGLOhour(zcnt + (tus/1000000));
    }

    gobs = getgobsgnss(epoch, siteid, (string)prncstring);
    
    GOBS signaltype;
    if(code == 0) {
      if (band == 0) {signaltype = L1C;}
      else if (band == 1) {signaltype = L1P;}
      else {continue;} //skip it, reserved
    }else{
      if (band == 0) {signaltype = L2C;}
      else if (band == 1) {signaltype = L2P;}
      else {continue;} //skip it, reserved
    }
    
    //compute PhaseRange from phase cycles if L1 Pseudorange is available as reference and add new obs
    double refrange = gobs->getobs(C1C);
    if ((int64_t)(refrange) != 0){
      double cylwl = 1/(gobs->frequency(t_gobs(signaltype).band())); // cycle length [m]
      double phrng = phcyls * cylwl; // phase range modulo overflow limit (has opposite sign to Pseudorange)
      cylwl = 8388608.0 * cylwl * 2; //overflow limit
      int32_t nov = refrange / cylwl;
      nov = (phcyls > 0)? nov + 1 : nov; //number of overflows
      phrng = nov * cylwl - phrng; //total phase range [m]
      gobs->addobs(signaltype, phrng);
    }
  }

  if (epochend == 0) cnt += passOBSd();
  return cnt;//(return number of records read, -1 means error has occurred, 0 means message was skipped for some reason)
}

int t_rtcm2::msg19(){//RTK Uncorrected Pseudoranges
  shared_ptr<t_gobsgnss> gobs;
  int cnt = 0;
  int epochend = 0;

  // ## READING INDIVIDUAL DATA FIELDS ##
  //common data
  int siteid = getnumber( 1, 14,10); //reference station number
  double zcnt = (double)getnumber( 2, 0,13) * 0.6; // Modified Z count [s]
  int messagewords = (int)getnumber(2,16,5)+2;
  int band = (int)getnumber(3,0,2);//Frequency... 0=L1, 1=L2, 2&3...Reserved
  getnumber(3,2,2);//Smoothing Interval
  double tus = getnumber(3,4,20);//GNSS Time Of Measurement [us]
  //satellite specific data
  int wordsread = 3;
  while(wordsread < messagewords){
    epochend = getnumber(wordsread+1,0,1);//MULTIPLE MESAGE INDICATOR, if 1: more messages with identical time-tag follow
    int code = getnumber(wordsread+1,1,1);//Code INDICATOR: 0 = C/A Code; 1 = P-Code
    int satsys = getnumber(wordsread+1,2,1);//GNSS system: 0 = GPS; 1 = GLONASS
    int satslot = getnumber(wordsread+1,3,5);//satellite ID ( slot number GLONASS / PRN GPS ); if GPS 0=32
    getnumber(wordsread+1,8,4);//Data quality
    getnumber(wordsread+1,12,4);//Multipath error
    double psrng = (double)(( getnumber(wordsread+1,16,8) <<24 ) + getnumber(wordsread+2,0,24)) * 0.02; //pseudorange [m]
    wordsread += 2;

    char prncstring[10];
    t_gtime epoch;
    if (satsys == 0){ //GPS
      sprintf(prncstring,"%c%02i", 'G', satslot);
      epoch = estimateepochGPShour(zcnt + (tus/1000000));
    }else{ //GLONASS
      sprintf(prncstring,"%c%02i", 'R', satslot);
      epoch = estimateepochGLOhour(zcnt + (tus/1000000));
    }

    gobs = getgobsgnss(epoch, siteid, (string)prncstring);

    GOBS signaltype;
    if(code == 0) {
      if (band == 0) {signaltype = C1C;}
      else if (band == 1) {signaltype = C1P;}
      else {continue;} //skip it, reserved
    }else{
      if (band == 0) {signaltype = C2C;}
      else if (band == 1) {signaltype = C2P;}
      else {continue;} //skip it, reserved
    }
    
    gobs->addobs(signaltype, psrng);
  }

  if (epochend == 0) cnt += passOBSd();
  return cnt; //(return number of records read, -1 means error has occurred, 0 means message was skipped for some reason)
}

