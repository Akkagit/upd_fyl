/*
 * rtcm3.cpp
 *
 *  Created on: Jan 6, 2015
 *      Author: jakub
 */

#include <iostream>

#include "gdata/geph.h"
#include "gdata/gobsgnss.h"
#include "gcoders/rtcm3.h"
#include "gall/gallobj.h"
#include "gall/gallnav.h"
#include "gall/gallobs.h"
#include "gutils/gtime.h"
using namespace gnut;

t_rtcm3::t_rtcm3() 
{
	for (int i = 0; i < 231; i++) {//all non-defined messages set to unsupported
		_messagetemplates[i] = &t_rtcm3::msguspt;
	}
	_messagetemplates[1] = &t_rtcm3::msg1001;//GPS L1 Observables
	_messagetemplates[2] = &t_rtcm3::msg1002;//GPS L1 Extended Observables
	_messagetemplates[3] = &t_rtcm3::msg1003;//GPS L1 + L2 Observables
	_messagetemplates[4] = &t_rtcm3::msg1004;//GPS L1 + L2 Extended Observables
	//...
	_messagetemplates[9] = &t_rtcm3::msg1009;//GLONASS L1 Observables
	_messagetemplates[10] = &t_rtcm3::msg1010;//GLONASS L1 Extended Observables
	_messagetemplates[11] = &t_rtcm3::msg1011;//GLONASS L1 + L2 Observables
	_messagetemplates[12] = &t_rtcm3::msg1012;//GLONASS L1 + L2 Extended Observables
	//...
	_messagetemplates[19] = &t_rtcm3::msg1019;//GPS ephemeris
	_messagetemplates[20] = &t_rtcm3::msg1020;//GLONASS ephemeris
	_messagetemplates[44] = &t_rtcm3::msg1044;//QZSS ephemeris
	_messagetemplates[45] = &t_rtcm3::msg1045;//Galileo F/NAV ephemeris
	//...
	_messagetemplates[71] = &t_rtcm3::msg1071;// GPS MSM1
	_messagetemplates[72] = &t_rtcm3::msg1072;// GPS MSM2
	_messagetemplates[73] = &t_rtcm3::msg1073;// GPS MSM3
	_messagetemplates[74] = &t_rtcm3::msg1074;// GPS MSM4
	_messagetemplates[75] = &t_rtcm3::msg1075;// GPS MSM5
	_messagetemplates[76] = &t_rtcm3::msg1076;// GPS MSM6
	_messagetemplates[77] = &t_rtcm3::msg1077;// GPS MSM7
	//... MSM Reserved
	_messagetemplates[81] = &t_rtcm3::msg1081;// GLONASS MSM1
	_messagetemplates[82] = &t_rtcm3::msg1082;// GLONASS MSM2
	_messagetemplates[83] = &t_rtcm3::msg1083;// GLONASS MSM3
	_messagetemplates[84] = &t_rtcm3::msg1084;// GLONASS MSM4
	_messagetemplates[85] = &t_rtcm3::msg1085;// GLONASS MSM5
	_messagetemplates[86] = &t_rtcm3::msg1086;// GLONASS MSM6
	_messagetemplates[87] = &t_rtcm3::msg1087;// GLONASS MSM7
	//... MSM Reserved
	_messagetemplates[91] = &t_rtcm3::msg1091;// Galileo MSM1
	_messagetemplates[92] = &t_rtcm3::msg1092;// Galileo MSM2
	_messagetemplates[93] = &t_rtcm3::msg1093;// Galileo MSM3
	_messagetemplates[94] = &t_rtcm3::msg1094;// Galileo MSM4
	_messagetemplates[95] = &t_rtcm3::msg1095;// Galileo MSM5
	_messagetemplates[96] = &t_rtcm3::msg1096;// Galileo MSM6
	_messagetemplates[97] = &t_rtcm3::msg1097;// Galileo MSM7
	//... MSM Reserved
	_messagetemplates[111] = &t_rtcm3::msg1111;// QZSS MSM1
	_messagetemplates[112] = &t_rtcm3::msg1112;// QZSS MSM2
	_messagetemplates[113] = &t_rtcm3::msg1113;// QZSS MSM3
	_messagetemplates[114] = &t_rtcm3::msg1114;// QZSS MSM4
	_messagetemplates[115] = &t_rtcm3::msg1115;// QZSS MSM5
	_messagetemplates[116] = &t_rtcm3::msg1116;// QZSS MSM6
	_messagetemplates[117] = &t_rtcm3::msg1117;// QZSS MSM7
	//... MSM Reserved
	_messagetemplates[121] = &t_rtcm3::msg1121;// BeiDou MSM1
	_messagetemplates[122] = &t_rtcm3::msg1122;// BeiDou MSM2
	_messagetemplates[123] = &t_rtcm3::msg1123;// BeiDou MSM3
	_messagetemplates[124] = &t_rtcm3::msg1124;// BeiDou MSM4
	_messagetemplates[125] = &t_rtcm3::msg1125;// BeiDou MSM5
	_messagetemplates[126] = &t_rtcm3::msg1126;// BeiDou MSM6
	_messagetemplates[127] = &t_rtcm3::msg1127;// BeiDou MSM7
	//... MSM Reserved

	initobjstate();
}

t_rtcm3::~t_rtcm3() {
	delete _legobsepoch;
}

void t_rtcm3::initobjstate() {
	_legobsepoch = new t_gtime();
	clearstate();
}
void t_rtcm3::clearstate() {
	_stmap.clear();//erase all data for LLI detection algorithm
	_msmepochdata.clear();

	_legobsepoch->from_ymdhms(1980, 0, 0, 0, 0, 0);
	_legepochdata.clear();


	_lastephdata.clear();
	_stationcnt.clear();

	// default time estimation tolerance values
	_maxdteph = 20000.0; //NOTE: default value for ephemeris is large to prevent unimplemented GLO time conversion from producing warnings
	_maxdtobs = 1000.0;

	t_rtcm_base::clearstate();
}
//Called on a position in the buffer determines whether this bit is a start of a new rtcm_v3 message
//Returns:
//Not message start                       ...  < 0
//Suspected message start (missing data)  ... == 0
//Confirmed message start                 ...  > 0 (== number of bits in the message); message type is filled properly
int t_rtcm3::messagefound(int bitnumber, int & messagetype) {
	if (bitnumber % 8 != 0) { return -1; } //ensure that location is at a beginning of byte
	if (bitnumber + 8 >= (*_buffsize) * 8) { return 0; } //ensure that buffer has enough data to start checking
	if (getuint(bitnumber, 8) != 0xD3) { return -1; } //does not match preamble
	else {//possible message start detected

	  //std::cout << " possible message start detected..." << std::endl;

	  //message structure: 8b preamble; 6b reserved; 10b size(Bytes); #MESSAGE#DATA#; 24b CRC
		if (bitnumber + 24 >= (*_buffsize) * 8) {

			std::cout << "buffsize: " << *_buffsize << " is too small" << std::endl;

			return 0;
		} //message not complete (impossible to read message length and check CRC)
		else {
			int datalength = (int)getuint(bitnumber + 14, 10);

			//std::cout << "message length reads: 48 + " << (datalength*8) << " bits" << std::endl;

			if (bitnumber + 48 + (datalength * 8) >= (*_buffsize) * 8) {

				std::cout << "buffsize: " << *_buffsize << " is too small" << std::endl;

				return 0;
			} //message not complete (impossible to check CRC)
			else {

				//std::cout << "Initiating CRC check..." << std::endl;
				//std::cout << "CRC from file=  " << getuint(bitnumber+24+(datalength*8), 24) << std::endl;

				unsigned char crcbuffer[1030] = { 0 };
				for (int x = 0; x < (6 + datalength); x++) {
					crcbuffer[x] = getuint(bitnumber + 8 * x, 8);

					//std::cout << "Byte"<< x << " = " << (int)crcbuffer[x] << std::endl;

				}

				long crcchecksum = getCrc24Q(crcbuffer, 3 + datalength); //compute CRC

				//std::cout << "CRC computed =  " << crcchecksum << std::endl;


				if (crcchecksum != (long)getuint(bitnumber + 24 + (datalength * 8), 24)) {//check CRC
				  //std::cout << "CRC mismatch." << std::endl;
					return -1;//failed... false positive start of message detection / data corrupted
				}
				else {
					if (datalength < 2) {

						std::cout << "Message type mismatch." << std::endl;

						return -1; //invalid message format: no message type specified...
					}
					else {

						//std::cout << "Message detected!" << std::endl;

						messagetype = (int)getuint(bitnumber + 24, 12); //read message type
						return 48 + (datalength * 8); //return
					}
				}
			}
		}
	}
}


//CRC24Q check algorithm
const int64_t t_rtcm3::crc24q[256] =
{ 0x00000000, 0x01864CFB, 0x028AD50D, 0x030C99F6, 0x0493E6E1, 0x0515AA1A, 0x061933EC, 0x079F7F17,
0x08A18139, 0x0927CDC2, 0x0A2B5434, 0x0BAD18CF, 0x0C3267D8, 0x0DB42B23, 0x0EB8B2D5, 0x0F3EFE2E,
0x10C54E89, 0x11430272, 0x124F9B84, 0x13C9D77F, 0x1456A868, 0x15D0E493, 0x16DC7D65, 0x175A319E,
0x1864CFB0, 0x19E2834B, 0x1AEE1ABD, 0x1B685646, 0x1CF72951, 0x1D7165AA, 0x1E7DFC5C, 0x1FFBB0A7,
0x200CD1E9, 0x218A9D12, 0x228604E4, 0x2300481F, 0x249F3708, 0x25197BF3, 0x2615E205, 0x2793AEFE,
0x28AD50D0, 0x292B1C2B, 0x2A2785DD, 0x2BA1C926, 0x2C3EB631, 0x2DB8FACA, 0x2EB4633C, 0x2F322FC7,
0x30C99F60, 0x314FD39B, 0x32434A6D, 0x33C50696, 0x345A7981, 0x35DC357A, 0x36D0AC8C, 0x3756E077,
0x38681E59, 0x39EE52A2, 0x3AE2CB54, 0x3B6487AF, 0x3CFBF8B8, 0x3D7DB443, 0x3E712DB5, 0x3FF7614E,
0x4019A3D2, 0x419FEF29, 0x429376DF, 0x43153A24, 0x448A4533, 0x450C09C8, 0x4600903E, 0x4786DCC5,
0x48B822EB, 0x493E6E10, 0x4A32F7E6, 0x4BB4BB1D, 0x4C2BC40A, 0x4DAD88F1, 0x4EA11107, 0x4F275DFC,
0x50DCED5B, 0x515AA1A0, 0x52563856, 0x53D074AD, 0x544F0BBA, 0x55C94741, 0x56C5DEB7, 0x5743924C,
0x587D6C62, 0x59FB2099, 0x5AF7B96F, 0x5B71F594, 0x5CEE8A83, 0x5D68C678, 0x5E645F8E, 0x5FE21375,
0x6015723B, 0x61933EC0, 0x629FA736, 0x6319EBCD, 0x648694DA, 0x6500D821, 0x660C41D7, 0x678A0D2C,
0x68B4F302, 0x6932BFF9, 0x6A3E260F, 0x6BB86AF4, 0x6C2715E3, 0x6DA15918, 0x6EADC0EE, 0x6F2B8C15,
0x70D03CB2, 0x71567049, 0x725AE9BF, 0x73DCA544, 0x7443DA53, 0x75C596A8, 0x76C90F5E, 0x774F43A5,
0x7871BD8B, 0x79F7F170, 0x7AFB6886, 0x7B7D247D, 0x7CE25B6A, 0x7D641791, 0x7E688E67, 0x7FEEC29C,
0x803347A4, 0x81B50B5F, 0x82B992A9, 0x833FDE52, 0x84A0A145, 0x8526EDBE, 0x862A7448, 0x87AC38B3,
0x8892C69D, 0x89148A66, 0x8A181390, 0x8B9E5F6B, 0x8C01207C, 0x8D876C87, 0x8E8BF571, 0x8F0DB98A,
0x90F6092D, 0x917045D6, 0x927CDC20, 0x93FA90DB, 0x9465EFCC, 0x95E3A337, 0x96EF3AC1, 0x9769763A,
0x98578814, 0x99D1C4EF, 0x9ADD5D19, 0x9B5B11E2, 0x9CC46EF5, 0x9D42220E, 0x9E4EBBF8, 0x9FC8F703,
0xA03F964D, 0xA1B9DAB6, 0xA2B54340, 0xA3330FBB, 0xA4AC70AC, 0xA52A3C57, 0xA626A5A1, 0xA7A0E95A,
0xA89E1774, 0xA9185B8F, 0xAA14C279, 0xAB928E82, 0xAC0DF195, 0xAD8BBD6E, 0xAE872498, 0xAF016863,
0xB0FAD8C4, 0xB17C943F, 0xB2700DC9, 0xB3F64132, 0xB4693E25, 0xB5EF72DE, 0xB6E3EB28, 0xB765A7D3,
0xB85B59FD, 0xB9DD1506, 0xBAD18CF0, 0xBB57C00B, 0xBCC8BF1C, 0xBD4EF3E7, 0xBE426A11, 0xBFC426EA,
0xC02AE476, 0xC1ACA88D, 0xC2A0317B, 0xC3267D80, 0xC4B90297, 0xC53F4E6C, 0xC633D79A, 0xC7B59B61,
0xC88B654F, 0xC90D29B4, 0xCA01B042, 0xCB87FCB9, 0xCC1883AE, 0xCD9ECF55, 0xCE9256A3, 0xCF141A58,
0xD0EFAAFF, 0xD169E604, 0xD2657FF2, 0xD3E33309, 0xD47C4C1E, 0xD5FA00E5, 0xD6F69913, 0xD770D5E8,
0xD84E2BC6, 0xD9C8673D, 0xDAC4FECB, 0xDB42B230, 0xDCDDCD27, 0xDD5B81DC, 0xDE57182A, 0xDFD154D1,
0xE026359F, 0xE1A07964, 0xE2ACE092, 0xE32AAC69, 0xE4B5D37E, 0xE5339F85, 0xE63F0673, 0xE7B94A88,
0xE887B4A6, 0xE901F85D, 0xEA0D61AB, 0xEB8B2D50, 0xEC145247, 0xED921EBC, 0xEE9E874A, 0xEF18CBB1,
0xF0E37B16, 0xF16537ED, 0xF269AE1B, 0xF3EFE2E0, 0xF4709DF7, 0xF5F6D10C, 0xF6FA48FA, 0xF77C0401,
0xF842FA2F, 0xF9C4B6D4, 0xFAC82F22, 0xFB4E63D9, 0xFCD11CCE, 0xFD575035, 0xFE5BC9C3, 0xFFDD8538 };

int64_t t_rtcm3::getCrc24Q(const unsigned char *buffer, int length)
{
	long long crc = 0;
	for (int i = 0; i < length; i++)
		crc = (crc << 8) ^ crc24q[(unsigned char)(buffer[i]) ^ ((crc >> 16) & 0xFF)];
	return crc & 0xFFFFFF;
}


//########################################################################################
//######################  DECODING FUNCTION DECLARATIONS  ################################
//########################################################################################
//keep in mind that since RTCM3 starts at the beginning of the byte & every checked byte is consumed, message always starts at position 0 in the buffer
int t_rtcm3::readmessage() {
	int decodeindex = (int)getuint(24, 12);
	decodeindex = ((decodeindex < 1001) | (decodeindex > 1230)) ? 0 : decodeindex - 1000;
	return (this->*_messagetemplates[decodeindex])();
}

string t_rtcm3::resolvestation(uint64_t stationid) {
	if (_stationcnt.find(stationid) != _stationcnt.end()) {
		_stationcnt[stationid]++;
	}
	else {
		if (!_stationcnt.empty()) {
			//TODO: generate warning, multiple station IDs in a single stream
		}
		_stationcnt[stationid] = 1;
	}
	char site[10] = "";//reference station cstring
	sprintf(site, "%04lluRTCM3", stationid);//NOTE: STRING DOES NOT MATCH RINEX IDs! RINEX has format of station_id: XXXXMRCCC, RTCM3 has number in range: 0..4095
	return string(site);
}

int t_rtcm3::msguspt() {
	int messagetype = (int)getuint(24, 12);
	int lastdigit = messagetype % 10;
	int recordsread = 0;
	if (messagetype >= 1070 && messagetype <= 1229 && lastdigit >= 1 && lastdigit <= 7) {//message is reserved for MSM1 - MSM7
		recordsread = passMSMd();
	}
	//TODO: this is info for log
	std::cout << "Decoding of RTCM3 message #" << ((messagetype == 0) ? 64 : messagetype) << " not supported!" << std::endl;
	return recordsread;
}
//########################################################################################
//###########################  DECODING RTK OBS MESSAGES  ################################
int t_rtcm3::gpsRTKh(uint64_t* headdata) { //headdata is expected to be an array with minimum size of 6 elements
	headdata[0] = getuint(36, 12); //reference station ID
	headdata[1] = getuint(48, 30); //TOW
	headdata[2] = getuint(78, 1); //Synchronous GNSS flag*
	headdata[3] = getuint(79, 5); //Number of satellite signals
	headdata[4] = getuint(84, 1); //Divergence free smoothing indicator
	headdata[5] = getuint(85, 3); //Smoothing interval
	//* data for all available satellites should be contained in a single message,
	//  if the SyncGNSSflag is set, data will be processed at first change of epoch, else every message has own epoch
	//  This does not indicate end of an epoch, only whether all systems observations are captured at the same moment or not
	return 88; //starting position to continue reading on the next bit after the header ended
}

int t_rtcm3::gloRTKh(uint64_t* headdata) { //headdata is expected to be an array with minimum size of 6 elements
	headdata[0] = getuint(36, 12); //reference station ID
	headdata[1] = getuint(48, 27); //GLONASS epoch time t_k [milliseconds]
	headdata[2] = getuint(75, 1); //Synchronous GNSS flag*
	headdata[3] = getuint(76, 5); //Number of satellite signals
	headdata[4] = getuint(81, 1); //Divergence free smoothing indicator
	headdata[5] = getuint(82, 3); //Smoothing interval
	//* data for all available satellites should be contained in a single message,
	//  if the SyncGNSSflag is set, data will be processed at first change of epoch, else every message has own epoch
	//  This does not indicate end of an epoch, only whether all systems observations are captured at the same moment or not
	return 85; //starting position to continue reading on the next bit after the header ended
}
shared_ptr<t_gobsgnss> t_rtcm3::leggetgobsgnss(t_gtime epoch, uint64_t siteid, string sat) {
	//Check if the stream has only a single station ID
	string site = resolvestation(siteid);

	//Search for observation object in the _legepochdata vector
	shared_ptr<t_gobsgnss> nullsptrgobs;
	shared_ptr<t_gobsgnss> resultsptr;
	vector<shared_ptr<t_gobsgnss>>::iterator itOBS = _legepochdata.begin();
	while (itOBS != _legepochdata.end()) {
		if ((*itOBS)->epoch() == epoch) { //NOTE: since legacy_obs epoch end is determined by comparison of epoch times this should always resolve as true
			if ((*itOBS)->site().compare(site) == 0) {
				if ((*itOBS)->sat().compare(sat) == 0) {
					resultsptr = (*itOBS);
				}
			}
		}
		itOBS++;
	}
	//if not found add newly constructed object
	if (resultsptr == nullsptrgobs) {
		resultsptr.reset(new t_gobsgnss(site, sat, epoch));
		_legepochdata.push_back(resultsptr);
	}
	return resultsptr;
}
int t_rtcm3::passRTKo() {
	int cnt = 0;

	//NOTE: taken from rinexo2, passing data from single epoch of MSM observations to t_gall
	map<string, t_gdata*>::iterator itDAT = (*_data).begin();
	// adding channel canal to t_gobsgnss from gallobj
	// !!!! must be before adding obs data to t_gallobs
	while (itDAT != (*_data).end()) {
		if (itDAT->second->id_type() != t_gdata::ALLOBJ) {
			itDAT++;
			continue;
		}
		vector<shared_ptr<t_gobsgnss>>::iterator itOBS = _legepochdata.begin();
		while (itOBS != _legepochdata.end()) {
			t_gallobj* all_obj = (t_gallobj*)itDAT->second;
			shared_ptr<t_grec> null_obj, one_obj = dynamic_pointer_cast<t_grec>(all_obj->obj((*itOBS)->sat()));
			if (one_obj != null_obj && one_obj->id_type() == t_gdata::TRN) {
				int ch = one_obj->channel();
				(*itOBS)->channel(ch);
			}
			itOBS++;
		}
		itDAT++;
	}
	// fill t_gobsgnss
	itDAT = (*_data).begin();
	while (itDAT != (*_data).end()) {

		if (itDAT->second->id_type() != t_gdata::ALLOBS) {
			itDAT++;
			continue;
		}

		vector<shared_ptr<t_gobsgnss>>::const_iterator itOBS = _legepochdata.begin();
		while (itOBS != _legepochdata.end()) {

			if ((*_log) && (*_log)->verb() >= 4)
				(*_log)->comment(1, "rtcm3", (*itOBS)->sat() + " obs filled "
					+ (*itOBS)->epoch().str_ymdhms()
					+ " " + (*itOBS)->site());

			((t_gallobs*)itDAT->second)->addobs(*itOBS);
			itOBS++;
		}
		cnt++;
		itDAT++;
	}
	_legepochdata.clear();

	return cnt;
}

int t_rtcm3::msg1001() {//GPS L1 Observables
	int cnt = 0;
	t_gtime epochtime;
	uint64_t headdata[6] = { 0 };
	int position = gpsRTKh(headdata); //read header data
	epochtime = estimateepochGPS(headdata[1] / 1000.0, true);
	//determine if this is start of new epoch, if it is, pass data to t_gall
	if (epochtime != (*_legobsepoch)) {
		cnt += passRTKo();
		(*_legobsepoch) = epochtime;
	}
	//read satellite data
	for (uint64_t i = 0; i < headdata[3]; i++) {
		char prncstring[10];
		sprintf(prncstring, "%c%02i", 'G', (int)getuint(position + 0, 6));//satellite ID

		getuint(position + 6, 1);//GPS L1 Code Indicator
		getuint(position + 7, 24);//GPS L1 Pseudorange
		get2int(position + 31, 20);//GPS L1 PhaseRange - Pseudorange
		getuint(position + 51, 7);//GPS L1 Lock time Indicator

		position += 58; //move position by length of one satellite block
		//NOTE: t_gobsgnss is not created, Pseudorange Modulus Ambiguity not available, impossible to calculate data
	}

	//if station does not take measurement in a synchronous manner, every message has it's own epoch
	if (headdata[2] == 0) cnt += passRTKo(); //NOTE: not sure about epoch times of different legacy messages of the same system though

	return cnt;
}

int t_rtcm3::msg1002() {//GPS L1 Extended Observables
	int cnt = 0;
	t_gtime epochtime;
	uint64_t headdata[6] = { 0 };
	int position = gpsRTKh(headdata); //read header data
	epochtime = estimateepochGPS(headdata[1] / 1000.0, true);
	//determine if this is start of new epoch, if it is, pass data to t_gall
	if (epochtime != (*_legobsepoch)) {
		cnt += passRTKo();
		(*_legobsepoch) = epochtime;
	}
	//read satellite data
	for (uint64_t i = 0; i < headdata[3]; i++) {
		char prncstring[10];
		sprintf(prncstring, "%c%02i", 'G', (int)getuint(position + 0, 6));//satellite ID

		int l1code = getuint(position + 6, 1);//GPS L1 Code Indicator
		double l1ps = getuint(position + 7, 24);//GPS L1 Pseudorange
		double l1ph = get2int(position + 31, 20);//GPS L1 PhaseRange - Pseudorange
		getuint(position + 51, 7);//GPS L1 Lock time Indicator
		double psma = getuint(position + 58, 8);//GPS Integer L1 Pseudorange Modulus Ambiguity
		getuint(position + 66, 8);//GPS L1 CNR

		position += 74; //move position by length of one satellite block

		shared_ptr<t_gobsgnss> gobs = leggetgobsgnss(epochtime, headdata[0], (string)prncstring);
		//TODO: special observation values (out of range/not available)
		l1ps = (0.02*l1ps) + (299792.458*psma);
		l1ph = l1ps + (0.0005*l1ph);

		if (l1code == 0) { //C/A code
			gobs->addobs(C1C, l1ps);
			gobs->addobs(L1C, l1ph);
		}
		else {            //P(Y) code
			gobs->addobs(C1P, l1ps);
			gobs->addobs(L1P, l1ph);
		}

	}

	//if station does not take measurement in a synchronous manner, every message has it's own epoch
	if (headdata[2] == 0) cnt += passRTKo(); //NOTE: not sure about epoch times of different legacy messages of the same system though

	return cnt;
}

int t_rtcm3::msg1003() {//GPS L1 + L2 Observables
	int cnt = 0;
	t_gtime epochtime;
	uint64_t headdata[6] = { 0 };
	int position = gpsRTKh(headdata); //read header data
	epochtime = estimateepochGPS(headdata[1] / 1000.0, true);
	//determine if this is start of new epoch, if it is, pass data to t_gall
	if (epochtime != (*_legobsepoch)) {
		cnt += passRTKo();
		(*_legobsepoch) = epochtime;
	}
	//read satellite data
	for (uint64_t i = 0; i < headdata[3]; i++) {
		char prncstring[10];
		sprintf(prncstring, "%c%02i", 'G', (int)getuint(position + 0, 6));//satellite ID

		getuint(position + 6, 1);//GPS L1 Code Indicator
		getuint(position + 7, 24);//GPS L1 Pseudorange
		get2int(position + 31, 20);//GPS L1 PhaseRange - L1 Pseudorange
		getuint(position + 51, 7);//GPS L1 Lock time Indicator

		getuint(position + 58, 2);//GPS L2 Code Indicator
		get2int(position + 60, 14);//GPS L2-L1 Pseudorange Difference
		get2int(position + 74, 20);//GPS L2 PhaseRange - L1 Pseudorange
		getuint(position + 94, 7);//GPS L2 Lock Time Indicator

		position += 101; //move position by length of one satellite block
		//NOTE: t_gobsgnss is not created, Pseudorange Modulus Ambiguity not available, impossible to calculate data
	}

	//if station does not take measurement in a synchronous manner, every message has it's own epoch
	if (headdata[2] == 0) cnt += passRTKo(); //NOTE: not sure about epoch times of different legacy messages of the same system though

	return cnt;
}

int t_rtcm3::msg1004() {//GPS L1 + L2 Extended Observables
	int cnt = 0;
	t_gtime epochtime;
	uint64_t headdata[6] = { 0 };
	int position = gpsRTKh(headdata); //read header data
	epochtime = estimateepochGPS(headdata[1] / 1000.0, true);
	//determine if this is start of new epoch, if it is, pass data to t_gall
	if (epochtime != (*_legobsepoch)) {
		cnt += passRTKo();
		(*_legobsepoch) = epochtime;
	}
	//read satellite data
	for (uint64_t i = 0; i < headdata[3]; i++) {
		char prncstring[10];
		sprintf(prncstring, "%c%02i", 'G', (int)getuint(position + 0, 6));//satellite ID

		int l1code = getuint(position + 6, 1);//GPS L1 Code Indicator
		double l1ps = getuint(position + 7, 24);//GPS L1 Pseudorange
		double l1ph = get2int(position + 31, 20);//GPS L1 PhaseRange - L1 Pseudorange
		getuint(position + 51, 7);//GPS L1 Lock time Indicator
		double psma = getuint(position + 58, 8);//GPS Integer L1 Pseudorange Modulus Ambiguity
		getuint(position + 66, 8);//GPS L1 CNR

		int l2code = getuint(position + 74, 2);//GPS L2 Code Indicator
		double l2ps = get2int(position + 76, 14);//GPS L2-L1 Pseudorange Difference
		double l2ph = get2int(position + 90, 20);//GPS L2 PhaseRange - L1 Pseudorange
		getuint(position + 110, 7);//GPS L2 Lock Time Indicator
		getuint(position + 117, 8);//GPS L2 CNR

		position += 125; //move position by length of one satellite block

		shared_ptr<t_gobsgnss> gobs = leggetgobsgnss(epochtime, headdata[0], (string)prncstring);
		//TODO: special observation values (out of range/not available)
		l1ps = (0.02*l1ps) + (299792.458*psma);
		l1ph = l1ps + (0.0005*l1ph);
		l2ps = l1ps + (0.02*l2ps);
		l2ph = l1ps + (0.0005*l2ph);

		if (l1code == 0) { //C/A code
			gobs->addobs(C1C, l1ps);
			gobs->addobs(L1C, l1ph);
		}
		else {            //P(Y) code
			gobs->addobs(C1P, l1ps);
			gobs->addobs(L1P, l1ph);
		}
		if (l2code == 0) { //C/A code or L2C code
			gobs->addobs(C2C, l2ps);
			gobs->addobs(L2C, l2ph);
		}
		else if (l2code == 1) { //P(Y) code direct //NOTE: or Y code based on Antispoofing
			gobs->addobs(C2P, l2ps);
			gobs->addobs(L2P, l2ph);
		}
		else if (l2code == 2) { //P(Y) code cross-correlated

		}
		else { // Correlated P/Y

		}

	}

	//if station does not take measurement in a synchronous manner, every message has it's own epoch
	if (headdata[2] == 0) cnt += passRTKo(); //NOTE: not sure about epoch times of different legacy messages of the same system though

	return cnt;
}

int t_rtcm3::msg1009() {//GLONASS L1 Observables
	int cnt = 0;
	t_gtime epochtime;
	uint64_t headdata[6] = { 0 };
	int position = gloRTKh(headdata); //read header data
	epochtime = estimateepochGLOday(headdata[1] / 1000.0);
	//determine if this is start of new epoch, if it is, pass data to t_gall
	if (epochtime != (*_legobsepoch)) {
		cnt += passRTKo();
		(*_legobsepoch) = epochtime;
	}
	//read satellite data
	for (uint64_t i = 0; i < headdata[3]; i++) {
		char prncstring[10];
		sprintf(prncstring, "%c%02i", 'R', (int)getuint(position + 0, 6));//satellite ID
		getuint(position + 7, 5);//GLONASS Satellite Frequency Channel Number

		getuint(position + 6, 1);//GLONASS L1 Code Indicator
		getuint(position + 12, 25);//GLONASS L1 Pseudorange
		get2int(position + 37, 20);//GLONASS L1 PhaseRange - Pseudorange
		getuint(position + 57, 7);//GLONASS L1 Lock time Indicator

		position += 64; //move position by length of one satellite block
		//NOTE: t_gobsgnss is not created, Pseudorange Modulus Ambiguity not available, impossible to calculate data
	}

	//if station does not take measurement in a synchronous manner, every message has it's own epoch
	if (headdata[2] == 0) cnt += passRTKo(); //NOTE: not sure about epoch times of different legacy messages of the same system though

	return cnt;
}

int t_rtcm3::msg1010() {//GLONASS L1 Extended Observables
	int cnt = 0;
	t_gtime epochtime;
	uint64_t headdata[6] = { 0 };
	int position = gloRTKh(headdata); //read header data
	epochtime = estimateepochGLOday(headdata[1] / 1000.0);
	//determine if this is start of new epoch, if it is, pass data to t_gall
	if (epochtime != (*_legobsepoch)) {
		cnt += passRTKo();
		(*_legobsepoch) = epochtime;
	}
	//read satellite data
	for (uint64_t i = 0; i < headdata[3]; i++) {
		char prncstring[10];
		sprintf(prncstring, "%c%02i", 'R', (int)getuint(position + 0, 6));//satellite ID
		int chindicator = getuint(position + 7, 5);//GLONASS Satellite Frequency Channel Number

		int l1code = getuint(position + 6, 1);//GLONASS L1 Code Indicator
		double l1ps = getuint(position + 12, 25);//GLONASS L1 Pseudorange
		double l1ph = get2int(position + 37, 20);//GLONASS L1 PhaseRange - Pseudorange
		getuint(position + 57, 7);//GLONASS L1 Lock time Indicator
		double psma = getuint(position + 64, 7);//GLONASS Integer L1 Pseudorange Modulus Ambiguity
		getuint(position + 71, 8);//GLONASS L1 CNR

		position += 79; //move position by length of one satellite block

		shared_ptr<t_gobsgnss> gobs = leggetgobsgnss(epochtime, headdata[0], (string)prncstring);
		gobs->channel(chindicator - 7);

		//TODO: special observation values (out of range/not available)
		l1ps = (0.02*l1ps) + (599584.916*psma);
		l1ph = l1ps + (0.0005*l1ph);

		if (l1code == 0) { //C/A code
			gobs->addobs(C1C, l1ps);
			gobs->addobs(L1C, l1ph);
		}
		else {            //P(Y) code
			gobs->addobs(C1P, l1ps);
			gobs->addobs(L1P, l1ph);
		}
	}

	//if station does not take measurement in a synchronous manner, every message has it's own epoch
	if (headdata[2] == 0) cnt += passRTKo(); //NOTE: not sure about epoch times of different legacy messages of the same system though

	return cnt;
}

int t_rtcm3::msg1011() {//GLONASS L1 + L2 Observables
	int cnt = 0;
	t_gtime epochtime;
	uint64_t headdata[6] = { 0 };
	int position = gloRTKh(headdata); //read header data
	epochtime = estimateepochGLOday(headdata[1] / 1000.0);
	//determine if this is start of new epoch, if it is, pass data to t_gall
	if (epochtime != (*_legobsepoch)) {
		cnt += passRTKo();
		(*_legobsepoch) = epochtime;
	}
	//read satellite data
	for (uint64_t i = 0; i < headdata[3]; i++) {
		char prncstring[10];
		sprintf(prncstring, "%c%02i", 'R', (int)getuint(position + 0, 6));//satellite ID
		getuint(position + 7, 5);//GLONASS Satellite Frequency Channel Number

		getuint(position + 6, 1);//GLONASS L1 Code Indicator
		getuint(position + 12, 25);//GLONASS L1 Pseudorange
		get2int(position + 37, 20);//GLONASS L1 PhaseRange - Pseudorange
		getuint(position + 57, 7);//GLONASS L1 Lock time Indicator

		getuint(position + 64, 2);//GLONASS L2 Code Indicator
		get2int(position + 66, 14);//GLONASS L2-L1 Pseudorange Difference
		get2int(position + 80, 20);//GLONASS L2 PhaseRange - L1 Pseudorange
		getuint(position + 100, 7);//GLONASS L2 Lock Time Indicator

		position += 107; //move position by length of one satellite block
		//NOTE: t_gobsgnss is not created, Pseudorange Modulus Ambiguity not available, impossible to calculate data
	}

	//if station does not take measurement in a synchronous manner, every message has it's own epoch
	if (headdata[2] == 0) cnt += passRTKo(); //NOTE: not sure about epoch times of different legacy messages of the same system though

	return cnt;
}

int t_rtcm3::msg1012() {//GLONASS L1 + L2 Extended Observables
	int cnt = 0;
	t_gtime epochtime;
	uint64_t headdata[6] = { 0 };
	int position = gloRTKh(headdata); //read header data
	epochtime = estimateepochGLOday(headdata[1] / 1000.0);
	//determine if this is start of new epoch, if it is, pass data to t_gall
	if (epochtime != (*_legobsepoch)) {
		cnt += passRTKo();
		(*_legobsepoch) = epochtime;
	}
	//read satellite data
	for (uint64_t i = 0; i < headdata[3]; i++) {
		char prncstring[10];
		sprintf(prncstring, "%c%02i", 'R', (int)getuint(position + 0, 6));//satellite ID
		int chindicator = getuint(position + 7, 5);//GLONASS Satellite Frequency Channel Number

		int l1code = getuint(position + 6, 1);//GLONASS L1 Code Indicator
		double l1ps = getuint(position + 12, 25);//GLONASS L1 Pseudorange
		double l1ph = get2int(position + 37, 20);//GLONASS L1 PhaseRange - Pseudorange
		getuint(position + 57, 7);//GLONASS L1 Lock time Indicator
		double psma = getuint(position + 64, 7);//GLONASS Integer L1 Pseudorange Modulus Ambiguity
		getuint(position + 71, 8);//GLONASS L1 CNR

		int l2code = getuint(position + 79, 2);//GLONASS L2 Code Indicator
		double l2ps = get2int(position + 81, 14);//GLONASS L2-L1 Pseudorange Difference
		double l2ph = get2int(position + 95, 20);//GLONASS L2 PhaseRange - L1 Pseudorange
		getuint(position + 115, 7);//GLONASS L2 Lock Time Indicator
		getuint(position + 122, 8);//GLONASS L2 CNR

		position += 130; //move position by length of one satellite block

		shared_ptr<t_gobsgnss> gobs = leggetgobsgnss(epochtime, headdata[0], (string)prncstring);
		gobs->channel(chindicator - 7);

		//TODO: special observation values (out of range/not available)
		l1ps = (0.02*l1ps) + (599584.916*psma);
		l1ph = l1ps + (0.0005*l1ph);
		l2ps = l1ps + (0.02*l2ps);
		l2ph = l1ps + (0.0005*l2ph);

		if (l1code == 0) { //C/A code
			gobs->addobs(C1C, l1ps);
			gobs->addobs(L1C, l1ph);
		}
		else {            //P(Y) code
			gobs->addobs(C1P, l1ps);
			gobs->addobs(L1P, l1ph);
		}
		if (l2code == 0) { //C/A code
			gobs->addobs(C2C, l2ps);
			gobs->addobs(L2C, l2ph);
		}
		else if (l2code == 1) { //P code
			gobs->addobs(C2P, l2ps);
			gobs->addobs(L2P, l2ph);
		}
		else { //Reserved
		}
	}

	//if station does not take measurement in a synchronous manner, every message has it's own epoch
	if (headdata[2] == 0) cnt += passRTKo(); //NOTE: not sure about epoch times of different legacy messages of the same system though

	return cnt;
}

//########################################################################################
//#############################  DECODING NAV MESSAGES  ##################################

int t_rtcm3::passEPHd(char gnsssys, string prnstr, t_gnav* newephdata) { // function for handling ephemeris data hand-off to t_gall structure
  // wrap pointer a smart pointer
	shared_ptr<t_gnav> geph;
	geph.reset(newephdata);

	//determine whether this message is a repetition of previous one; if not add it
	std::map<char, shared_ptr<t_gnav>>::iterator itlasteph;
	int recordspassed = 0;
	itlasteph = _lastephdata.find(gnsssys);
	if (itlasteph == _lastephdata.end() || (itlasteph->second->epoch() != newephdata->epoch())) {

		// ADDING EPHEMERIS DATA
		recordspassed = 1;

		// skip non-valid messages
		set<string> msg;

		if (geph->chk(msg) < 0) {
			if (*_log) (*_log)->comment(1, "rtcm3", "skipped ephemeris message due to failed check. ");
			recordspassed = 0;
		}
		else {

			// add to _data
			if (*_log) (*_log)->comment(2, "rtcm3", "add sat nav " + prnstr + " " + newephdata->epoch().str_ymdhms(" "));
			map<string, t_gdata*>::iterator it = (*_data).begin();
			while (it != (*_data).end()) {
				if (it->second->id_type() == t_gdata::ALLNAV || it->second->id_group() == t_gdata::GRP_EPHEM) {
					((t_gallnav*)it->second)->add(geph);

					// if object container and GLONASS ephemeris, add channel
				}
				else if (it->second->id_type() == t_gdata::ALLOBJ) {
					if (geph->id_type() == t_gdata::EPHGLO) {
						t_gallobj* all_obj = (t_gallobj*)it->second;
						shared_ptr<t_gobj>  one_obj = all_obj->obj(geph->sat());
						shared_ptr<t_gobj>  objnull;
						int ch = dynamic_pointer_cast<t_gnavglo>(geph)->channel();
						if (one_obj != objnull) {
							one_obj->channel(ch);
						}
					}
				}
				it++;
			}
		}
		recordspassed = 1;
	}
	else { //skip adding ephemeris, one with identical epoch was added with this stream before
		if (*_log) (*_log)->comment(2, "rtcm3", "skipped sat nav " + prnstr + " " + newephdata->epoch().str_ymdhms(" ") + " due to repeating message.");
	}
	_lastephdata[gnsssys] = geph;
	return recordspassed;
}

int t_rtcm3::msg1019() {//GPS Ephemeris
	t_gnavdata data = { 0 };
	t_gtime epoch;
	string prn;
	char prncstring[10];

	sprintf(prncstring, "%c%02i", 'G', (int)getuint(36, 6));//satellite ID
	data[21] = (double)getuint(42, 10);//week number
	data[23] = (double)getuint(52, 4);//SV accuracy
	data[20] = (double)getuint(56, 2);//code on L2
	data[19] = bmag(get2int(58, 14), -43) * G_PI;//IDOT
	data[3] = (double)getuint(72, 8);//IODE
	long toc = bmag(getuint(80, 16), 4);//t_oc
	data[2] = bmag(get2int(96, 8), -55);//a_f2
	data[1] = bmag(get2int(104, 16), -43);//a_f1
	data[0] = bmag(get2int(120, 22), -31);//a_f0
	data[26] = (double)getuint(142, 10);//IODC
	data[4] = bmag(get2int(152, 16), -5);//C_rs
	data[5] = bmag(get2int(168, 16), -43) * G_PI;//DELTA n
	data[6] = bmag(get2int(184, 32), -31) * G_PI;//M_0
	data[7] = bmag(get2int(216, 16), -29);//C_uc
	data[8] = bmag(getuint(232, 32), -33);//eccentricity
	data[9] = bmag(get2int(264, 16), -29);//C_us
	data[10] = bmag(getuint(280, 32), -19);//A^1/2
	data[11] = bmag(getuint(312, 16), 4);//t_oe
	data[12] = bmag(get2int(328, 16), -29);//C_ic
	data[13] = bmag(get2int(344, 32), -31) * G_PI;//OMEGA_0
	data[14] = bmag(get2int(376, 16), -29);//C_is
	data[15] = bmag(get2int(392, 32), -31) * G_PI;//i_0
	data[16] = bmag(get2int(424, 16), -5);//C_rc
	data[17] = bmag(get2int(440, 32), -31) * G_PI;//omega (argument of perigee)
	data[18] = bmag(get2int(472, 24), -43) * G_PI;//OMEGADOT
	data[25] = bmag(get2int(496, 8), -31);//t_GD
	data[24] = (double)getuint(504, 6);//SV health
	data[22] = (double)getuint(510, 1);//L2P data flag
	data[28] = (double)getuint(511, 1);//fit interval

	//Time estimation
	double wnumber = data[21];
	epoch = estimateepochGPS(toc + (wnumber*604800.0));

	/*
	  std::cout << "PRN: " << prncstring << std::endl;
	  std::cout << "GPSTIME: " << epoch.str_ymdhms() << std::endl;
	  std::cout << "Toc = " << toc << " IODC(7-0) = " << (((uint64_t)data[26])&0xFF) << " IODE = " << data[3] << " Toe = " << data[11] << "\n" << std::endl;
	*/
	//compile into data structure
	t_gnav* geph = new t_gnavgps();
	geph->glog(*_log);
	geph->data2nav((string)prncstring, epoch, data);

	//pass data to t_gall
	return passEPHd('G', string(prncstring), geph); //(return number of records read, -1 means error has occurred, 0 means message was not passed for some reason)
}

int t_rtcm3::msg1020() {//GLONASS Ephemeris
	t_gnavdata data = { 0 };
	t_gtime epoch;
	string prn;
	char prncstring[10];

	sprintf(prncstring, "%c%02i", 'R', (int)getuint(36, 6));//GLONASS satellite ID (slot number)
	data[10] = (double)getuint(42, 5) - 7;//GLONASS Frequency Channel Number
	//getuint(47,1);//almanac health C_n word
	//getuint(48,1);//almanac health availability indicator
	//int p1 = getuint(49,2);//GLONASS P1 // NOTE: there is a discrepancy in the standards involving GLONASS time calculation from t_b based on this parameter
	data[2] = (double)getuint(51, 12);//t_k (bit field)
	data[6] = (double)getuint(63, 1);//MSB B_n word (satellite health status)
	//getuint(64,1);//GLONASS P2
	long t_b = getuint(65, 7);//t_b
	data[4] = bmag(getsint(72, 24), -20);//x_d
	data[3] = bmag(getsint(96, 27), -11);//x
	data[5] = bmag(getsint(123, 5), -30);//x_dd
	data[8] = bmag(getsint(128, 24), -20);//y_d
	data[7] = bmag(getsint(152, 27), -11);//y
	data[9] = bmag(getsint(179, 5), -30);//y_dd
	data[12] = bmag(getsint(184, 24), -20);//z_d
	data[11] = bmag(getsint(208, 27), -11);//z
	data[13] = bmag(getsint(235, 5), -30);//z_dd
	//getuint(240,1);//GLONASS P3
	data[1] = bmag(getsint(241, 11), -30);//gama_n
	//getuint(252,2);//GLONASS-M P
	//getuint(254,1);//GLONASS-M l_n
	data[0] = -bmag(getsint(255, 22), -30);//GLONASS tau_n
	//getsint(277,5);//GLONASS-M DELTA(tau_n)
	data[14] = (double)getuint(282, 5);//E_n (age of operation)
	//getuint(287,1);//GLONASS-M P4
	//getuint(288,4);//GLONASS-M F_t?
	int dayof4yM = getuint(292, 11);//GLONASS-M N_T
	int ismm = getuint(303, 2);//GLONASS-M M, GLONASS M data availability
	//getuint(305,1);//availability of additional data
	int dayof4y = getuint(306, 11);//GLONASS N^A
	//getsint(317,32);//tau_c
	//int epoch4y = getuint(349,5);//GLONASS-M N_4
	//get2int(354,22);//GLONASS-M tau_GPS
	//getuint(376,1);//GLONASS-M l_n (fifth string)
	//getuint(377,7);//reserved

	if (ismm != 0)dayof4y = dayof4yM;

	//Time estimation
	double dnumber = dayof4y;
	epoch = estimateepochGLO((dnumber*86400.0) + ((double)t_b*900.0));

	/*
	  std::cout << "SVID:" << prncstring << std::endl;
	  std::cout << "UTCTIME: " << epoch.str_ymdhms() << std::endl;
	  std::cout << "T_b = " << t_b*0.25 << " Time (t_k)= "<< data[2] << " = " << ( ((((uint64_t)data[2])>>7)&0x1F)-3 )<< "h" << ((((uint64_t)data[2])>>1)&0x3F) << "m" << (((uint64_t)data[2])%2)*30 << " Y4 = " << epoch4y<< " D4 = " << dayof4y << "\n" << std::endl;
	*/
	//compile into data structure
	t_gnav* geph = new t_gnavglo();
	geph->glog(*_log);
	geph->data2nav((string)prncstring, epoch, data);

	//pass data to t_gall
	return passEPHd('R', string(prncstring), geph); //(return number of records read, -1 means error has occurred, 0 means message was not passed for some reason)
}

int t_rtcm3::msg1044() {//QZSS Ephemeris
	t_gnavdata data = { 0 };
	t_gtime epoch;
	string prn;
	char prncstring[10];

	sprintf(prncstring, "%c%02i", 'J', (int)(getuint(36, 4) + 192));//satellite ID
	long toc = bmag(getuint(40, 16), 4);//toc
	data[2] = bmag(get2int(56, 8), -55);//af2
	data[1] = bmag(get2int(64, 16), -43);//af1
	data[0] = bmag(get2int(80, 22), -31);//af0
	data[3] = (double)getuint(102, 8);//IODE
	data[4] = bmag(get2int(110, 16), -5);//crs
	data[5] = bmag(get2int(126, 16), -43) * G_PI;//DELTA n0
	data[6] = bmag(get2int(142, 32), -31) * G_PI;//M0
	data[7] = bmag(get2int(174, 16), -29);//cuc
	data[8] = bmag(getuint(190, 32), -33);//e
	data[9] = bmag(get2int(222, 16), -29);//cus
	data[10] = bmag(getuint(238, 32), -19);//A^1/2
	data[11] = (double)getuint(270, 16) * 16;//toe
	data[12] = bmag(get2int(286, 16), -29);//cic
	data[13] = bmag(get2int(302, 32), -31) * G_PI;//OMEGA0
	data[14] = bmag(get2int(334, 16), -29);//cis
	data[15] = bmag(get2int(350, 32), -31) * G_PI;//i0
	data[16] = bmag(get2int(382, 16), -5);//crc
	data[17] = bmag(get2int(398, 32), -31) * G_PI;//omega n
	data[18] = bmag(get2int(430, 24), -43) * G_PI;//OMEGA DOT
	data[19] = bmag(get2int(454, 14), -43) * G_PI;//i0 DOT
	data[20] = (double)getuint(468, 2);//L2 codes
	data[21] = (double)getuint(470, 10);//week number 0..1023
	data[23] = (double)getuint(480, 4) * 5.73;//URA ~ ACC
	data[24] = (double)getuint(484, 6);//SV health
	data[25] = bmag(get2int(490, 8), -31);//TGD (L1-C/A group delay)
	data[26] = (double)getuint(498, 10);//IODC
	data[28] = (double)getuint(508, 1);//fit interval

	//NOTE: RTCM3 standard states resolution of 'toc' and 'toe' data fields to be 24, while referencing the QZSS ICD, where this value is 16
	//      based on experiment, the value from the ICD has been found correct and is being used

	//Time estimation
	double wnumber = data[21];
	epoch = estimateepochQZS(toc + (wnumber*604800.0));

	/*
	  std::cout << "PRN: " << prncstring << std::endl;
	  std::cout << "GPSTIME: " << epoch.str_ymdhms() << std::endl;
	  std::cout << "Toc = " << toc << " IODC(7-0) = " << (((uint64_t)data[26])&0xFF) << " IODE = " << data[3] << " Toe = " << data[11] << "\n" << std::endl;
	*/
	//compile into data structure
	t_gnav* geph = new t_gnavqzs();
	geph->glog(*_log);
	geph->data2nav((string)prncstring, epoch, data);

	//pass data to t_gall
	return passEPHd('J', string(prncstring), geph); //(return number of records read, -1 means error has occurred, 0 means message was not passed for some reason)
}

int t_rtcm3::msg1045() {//Galileo Ephemeris
	t_gnavdata data = { 0 };
	t_gtime epoch;
	string prn;
	char prncstring[10];

	sprintf(prncstring, "%c%02i", 'E', (int)getuint(36, 6));//satellite ID
	data[21] = (double)getuint(42, 12);//Galileo week number
	data[3] = (double)getuint(54, 10);//Galileo IODnav ?= iode
	double sisa = (double)getuint(64, 8);//Galileo SV SISA
	data[23] = (sisa < 50) ? (0.01*sisa) : ((sisa < 75) ? (0.5 + 0.02*(sisa - 50)) : ((sisa < 100) ? (1 + 0.04*(sisa - 75)) : ((sisa < 126) ? (2 + 0.16*(sisa - 100)) : 0)));// SISA -> (acc)[m]
	data[19] = bmag(get2int(72, 14), -43) * G_PI;// IDOT
	long toc = (long)getuint(86, 14) * 60;//toc
	data[2] = bmag(get2int(100, 6), -59);//af2
	data[1] = bmag(get2int(106, 21), -46);//af1
	data[0] = bmag(get2int(127, 31), -34);//af0
	data[4] = bmag(get2int(158, 16), -5);//crs
	data[5] = bmag(get2int(174, 16), -43) * G_PI;//DELTA n
	data[6] = bmag(get2int(190, 32), -31) * G_PI;//M0
	data[7] = bmag(get2int(222, 16), -29);//cuc
	data[8] = bmag(getuint(238, 32), -33);//e
	data[9] = bmag(get2int(270, 16), -29);//cus
	data[10] = bmag(getuint(286, 32), -19);//a^1/2
	data[11] = (double)getuint(318, 14) * 60;//toe
	data[12] = bmag(get2int(332, 16), -29);//cic
	data[13] = bmag(get2int(348, 32), -31) * G_PI;//OMEGA0
	data[14] = bmag(get2int(380, 16), -29);//cis
	data[15] = bmag(get2int(396, 32), -31) * G_PI;//i0
	data[16] = bmag(get2int(428, 16), -5);//crc
	data[17] = bmag(get2int(444, 32), -31) * G_PI;//omega
	data[18] = bmag(get2int(476, 24), -43) * G_PI;//OMEGADOT
	data[25] = bmag(get2int(500, 10), -32);//BGD_E5a (group delay)
	//getuint(510, 2);//Galileo NAV SIGNAL HEALTH STATUS (OSHS)
	//getuint(512,1);//Galileo NAV DATA VALIDITY STATUS (OSDVS)
	//getuint(513,7);//reserved

	//Time estimation
	double wnumber = data[21];
	epoch = estimateepochGAL(toc + (wnumber*604800.0));

	/*
	  std::cout << "PRN: " << prncstring << std::endl;
	  std::cout << "GALTIME: " << epoch.str_ymdhms() << std::endl;
	  std::cout << "Toc = " << toc << " IODE = " << data[3] << " Toe = " << data[11] << "\n" << std::endl;
	*/
	//compile into data structure
	t_gnav* geph = new t_gnavgal();
	geph->glog(*_log);
	geph->data2nav((string)prncstring, epoch, data);

	//pass data to t_gall
	return passEPHd('E', string(prncstring), geph); //(return number of records read, -1 means error has occurred, 0 means message was not passed for some reason)
}


//########################################################################################
//#############################  DECODING MSM MESSAGES  ##################################

int t_rtcm3::passMSMd() {
	if (getuint(78, 1) == 1) { //Multiple Message Bit (55th bit of message +24bit header) @ RTCM3.2 section-3.5.15.3.8
		return 0; //epoch continues
	}
	int cnt = 0; //counter of records passed to G-NUT

	//NOTE: taken from rinexo2, passing data from single epoch of MSM observations to t_gall
	map<string, t_gdata*>::iterator itDAT = (*_data).begin();
	// adding channel canal to t_gobsgnss from gallobj
	// !!!! must be before adding obs data to t_gallobs
	while (itDAT != (*_data).end()) {
		if (itDAT->second->id_type() != t_gdata::ALLOBJ) {
			itDAT++;
			continue;
		}
		vector<shared_ptr<t_gobsgnss>>::iterator itOBS = _msmepochdata.begin();
		while (itOBS != _msmepochdata.end()) {
			t_gallobj* all_obj = (t_gallobj*)itDAT->second;
			shared_ptr<t_grec> null_obj, one_obj = dynamic_pointer_cast<t_grec>(all_obj->obj((*itOBS)->sat()));
			if (one_obj != null_obj && one_obj->id_type() == t_gdata::TRN) {
				int ch = one_obj->channel();
				(*itOBS)->channel(ch);
			}
			itOBS++;
		}
		itDAT++;
	}
	// fill t_gobsgnss
	itDAT = (*_data).begin();
	while (itDAT != (*_data).end()) {

		if (itDAT->second->id_type() != t_gdata::ALLOBS) {
			itDAT++;
			continue;
		}

		vector<shared_ptr<t_gobsgnss>>::const_iterator itOBS = _msmepochdata.begin();
		while (itOBS != _msmepochdata.end()) {

			if ((*_log) && (*_log)->verb() >= 4)
				(*_log)->comment(1, "rtcm3", (*itOBS)->sat() + " obs filled "
					+ (*itOBS)->epoch().str_ymdhms()
					+ " " + (*itOBS)->site());

			((t_gallobs*)itDAT->second)->addobs(*itOBS);
			itOBS++;
		}
		cnt++;
		itDAT++;
	}

	_msmepochdata.clear();
	return cnt;
}

//MSM head
int t_rtcm3::msgMSMh(uint64_t* datafield, int& nsat, int& nsig) {
	datafield[0] = getuint(36, 12); //reference station ID
	datafield[1] = getuint(48, 30); //GNSS epoch time
	datafield[2] = getuint(78, 1); //Multiple Message Bit
	datafield[3] = getuint(79, 3); //IODS
	//datafield[ 4] = getuint( 82, 7); //reserved
	datafield[5] = getuint(89, 2); //clock steering indicator
	datafield[6] = getuint(91, 2); //external clock indicator
	datafield[7] = getuint(93, 1); //divergence-free smoothing indicator
	datafield[8] = getuint(94, 3); //smoothing interval
	datafield[9] = getuint(97, 64); //GNSS Satellite Mask
	datafield[10] = getuint(161, 32); //GNSS Signal Mask
	//determine the size of the 'GNSS Cell Mask' field
	nsat = 0;
	nsig = 0;
	for (int i = 0; i < 64; i++) {
		nsat = (((datafield[9] >> i) & (uint64_t)1) != 0) ? nsat + 1 : nsat;
		nsig = (((datafield[10] >> i) & (uint64_t)1) != 0) ? nsig + 1 : nsig;
	}
	if (nsat*nsig > 64)return -1;
	datafield[11] = getuint(193, nsat*nsig);//GNSS Cell Mask
	return (193 + nsat*nsig);//starting position for reading the next field
}

//satellite id decoding
bool t_rtcm3::msmgetsatid(char* satidcstr, int msmsatmaskid, char gnsssys) {
	bool result = false;
	msmsatmaskid++;//converts index (0..63) into satellite slot/PRN numbering (1..64)
	if (gnsssys == 'G') {        //GPS
		sprintf(satidcstr, "%c%02i", gnsssys, msmsatmaskid);
		result = ((msmsatmaskid < 64) && (msmsatmaskid > 0)) ? true : false;
	}
	else if (gnsssys == 'R') {  //GLO
		sprintf(satidcstr, "%c%02i", gnsssys, msmsatmaskid);
		result = ((msmsatmaskid < 25) && (msmsatmaskid > 0)) ? true : false;
	}
	else if (gnsssys == 'E') {  //GAL
		sprintf(satidcstr, "%c%02i", gnsssys, msmsatmaskid);
		result = ((msmsatmaskid < 51) && (msmsatmaskid > 0)) ? true : false; //NOTE: number 51 & 52 are GIOVE-A/B, since they are not described in RINEX they are ignored as non-valid...
	}
	else if (gnsssys == 'J') {  //QZS
		sprintf(satidcstr, "%c%02i", gnsssys, /*192+*/msmsatmaskid);
		result = ((msmsatmaskid < 11) && (msmsatmaskid > 0)) ? true : false;
	}
	else if (gnsssys == 'C') {  //BDS
		sprintf(satidcstr, "%c%02i", gnsssys, msmsatmaskid);
		result = ((msmsatmaskid < 38) && (msmsatmaskid > 0)) ? true : false;
		/*}else if(gnsssys == 'S'){  //SBAS
			sprintf( satidcstr, "%c%02i", gnsssys, msmsatmaskid );
			result = ( (msmsatmaskid < 64) && (msmsatmaskid > 0) )? true : false; */
	}
	else {
		result = false; //unrecognized satellite system
	}
	return result;
}

//observation results' object setup/lookup for MSM reading functions; resolves epoch, station, satellite, signal code lookup tables
shared_ptr<t_gobsgnss> t_rtcm3::msmgetgobsgnss(uint64_t* msmheaddata, int satmaskindex, char gnsssys, const GOBS* sigmapC, const GOBS* sigmapL, const GOBS* sigmapD/*, const GOBS* sigmapS*/) {
	//Check if the stream has only a single station ID
	string stid = resolvestation(msmheaddata[0]); //reference STation ID

	//resolve satellite and station IDs (convert to RINEX style)
	char svid[10] = "";//Satellite Vehicle ID
	t_gtime obstime;
	shared_ptr<t_gobsgnss> nullsptrgobs;
	if (!msmgetsatid(svid, satmaskindex, gnsssys)) {
		return nullsptrgobs; //TODO: produce a warning here? unrecognized satellite
	}

	// use gnsssys to resolve epoch and assign proper signal code tables
	double tow; //time of week, used as temporary, scale differs based on system
	switch (gnsssys) {
	case 'G':
		tow = msmheaddata[1];
		obstime = estimateepochGPS(tow / 1000.0, true);
		sigmapC = gpsCtable;
		sigmapL = gpsLtable;
		sigmapD = gpsDtable;
		//sigmapS=gpsStable;
		break;
	case 'R':break;
		tow = (((msmheaddata[1] >> 27) & 7) * 86400.0) + ((double)(msmheaddata[1] & 0x7FFFFFF) / 1000.0); // (seconds in weeks previous days) + (seconds of day)
		obstime = estimateepochGLO(tow, true);
		sigmapC = gloCtable;
		sigmapL = gloLtable;
		sigmapD = gloDtable;
		//sigmapS=gloStable;
		break;
	case 'E':break;
		tow = msmheaddata[1];
		obstime = estimateepochGAL(tow / 1000.0, true);
		sigmapC = galCtable;
		sigmapL = galLtable;
		sigmapD = galDtable;
		//sigmapS=galStable;
		break;
	case 'J':break;
		tow = msmheaddata[1];
		obstime = estimateepochQZS(tow / 1000.0, true);
		sigmapC = qzsCtable;
		sigmapL = qzsLtable;
		sigmapD = qzsDtable;
		//sigmapS=qzsStable;
		break;
	case 'C':break;
		tow = msmheaddata[1];
		obstime = estimateepochBDS(tow / 1000.0, true);
		sigmapC = bdsCtable;
		sigmapL = bdsLtable;
		sigmapD = bdsDtable;
		//sigmapS=bdsStable;
		break;
		//default:
	}

	//look for existing t_gobsgnss object
	shared_ptr<t_gobsgnss> targetgobs;
	for (uint32_t i = 0; i < _msmepochdata.size(); i++) {
		if (_msmepochdata[i]->epoch() == obstime) { //epoch-check, times should match
			if (_msmepochdata[i]->site().compare(stid) == 0) {
				if (_msmepochdata[i]->sat().compare((string)svid) == 0) {
					targetgobs = _msmepochdata[i];
				}
			}
		}
		else { // epoch-check failed
		 //TODO: warning, time difference within epoch
			return nullsptrgobs;//TODO: what if the time does not match?(discard:merge:duplicate)? discard by default
		}
	}

	//generate new t_gobsgnss object if data is not discarded and object does not yet exist
	if (targetgobs == nullsptrgobs) {
		targetgobs.reset(new t_gobsgnss(stid, svid, obstime));
		_msmepochdata.push_back(targetgobs);
	}
	return targetgobs;
}

//    > DATA RECONSTRUCTION <
//PSeudo-range reconstruction
double t_rtcm3::msmgetPS(int64_t n_ms, int32_t roughdata, int32_t finedata, bool highprecision) {
	int finemag = (highprecision) ? -29 : -24;
	double resultdata = CLIGHT / 1000 * (n_ms + bmag(roughdata, -10)*bmag(finedata, finemag));//pseudo-range [m]
	return resultdata;
}

//PHase-cycle / PHase-range reconstruction; n_ms should must be specified (> 0) in order to get meaningful result
double t_rtcm3::msmgetPH(int64_t n_ms, int32_t roughdata, int32_t finedata, double frequency, bool highprecision) {
	if (n_ms < 0) n_ms = 0;
	int finemag = (highprecision) ? -31 : -29;
	double resultdata = CLIGHT / 1000 * (n_ms + bmag(roughdata, -10)*bmag(finedata, finemag));//phase-range [m]
	resultdata = ((resultdata*frequency / CLIGHT));//phase-cycle [full cycles]
	return resultdata;
}

//Doppler-Shift / phase-range-rate reconstruction
double t_rtcm3::msmgetDS(int32_t roughdata, int32_t finedata, double frequency, bool highprecision) {
	//reconstruction of high precision and standard precision PhaseRangeRates is identical for RTCM3 MSMs
	double resultdata = (double)roughdata + (double)finedata / 10000;//Phase-Range-Rate [m/s]
	resultdata = frequency*(CLIGHT / (CLIGHT - resultdata) - 1);//RINEX-like DOPPLER [Hz]... = f_dshifted - f_base
	return resultdata;
}

//Signal-Strength reconstruction
//double t_rtcm3::msmgetSS(int64_t n_ms, int32_t roughdata, int32_t finedata, bool highprecision=false){}


t_rtcm3::obslockdata t_rtcm3::msmgetMLT(uint32_t locktimeindicator, bool highprecision) {//return minimum lock time [ms]
	obslockdata lockdata;
	if (highprecision) {
		if (locktimeindicator < 64) { lockdata = { 0,(locktimeindicator),       1 }; }
		else if (locktimeindicator < 96) { lockdata = { 0,(2 * locktimeindicator - 64),       2 }; }
		else if (locktimeindicator < 128) { lockdata = { 0,(4 * locktimeindicator - 256),       4 }; }
		else if (locktimeindicator < 160) { lockdata = { 0,(8 * locktimeindicator - 768),       8 }; }
		else if (locktimeindicator < 192) { lockdata = { 0,(16 * locktimeindicator - 2048),      16 }; }
		else if (locktimeindicator < 224) { lockdata = { 0,(32 * locktimeindicator - 5120),      32 }; }
		else if (locktimeindicator < 256) { lockdata = { 0,(64 * locktimeindicator - 12288),      64 }; }
		else if (locktimeindicator < 288) { lockdata = { 0,(128 * locktimeindicator - 28672),     128 }; }
		else if (locktimeindicator < 320) { lockdata = { 0,(256 * locktimeindicator - 65536),     256 }; }
		else if (locktimeindicator < 352) { lockdata = { 0,(512 * locktimeindicator - 147456),     512 }; }
		else if (locktimeindicator < 384) { lockdata = { 0,(1024 * locktimeindicator - 327680),    1024 }; }
		else if (locktimeindicator < 416) { lockdata = { 0,(2048 * locktimeindicator - 720896),    2048 }; }
		else if (locktimeindicator < 448) { lockdata = { 0,(4096 * locktimeindicator - 1572864),    4096 }; }
		else if (locktimeindicator < 480) { lockdata = { 0,(8192 * locktimeindicator - 3407872),    8192 }; }
		else if (locktimeindicator < 512) { lockdata = { 0,(16384 * locktimeindicator - 7340032),   16384 }; }
		else if (locktimeindicator < 544) { lockdata = { 0,(32768 * locktimeindicator - 15728640),   32768 }; }
		else if (locktimeindicator < 576) { lockdata = { 0,(65536 * locktimeindicator - 33554432),   65536 }; }
		else if (locktimeindicator < 608) { lockdata = { 0,(131072 * locktimeindicator - 71303168),  131072 }; }
		else if (locktimeindicator < 640) { lockdata = { 0,(262144 * locktimeindicator - 150994944),  262144 }; }
		else if (locktimeindicator < 672) { lockdata = { 0,(524288 * locktimeindicator - 318767104),  524288 }; }
		else if (locktimeindicator < 704) { lockdata = { 0,(1048576 * locktimeindicator - 671088640), 1048576 }; }
		else if (locktimeindicator == 704) { lockdata = { 0,(2097152 * locktimeindicator - 1409286144), 2097152 }; }
		else { lockdata = { 0,                         0              ,       0 }; }//reserved or invalid
	}
	else {
		lockdata = { 0, 0, 0 }; // suppcoeff == 0 => low precision
		lockdata.locktime = (((locktimeindicator > 16) || (locktimeindicator <= 0)) ? 0 : (1 << (4 + locktimeindicator)));
	}
	return lockdata;
}

bool t_rtcm3::msmlligetprevtimes(obslockdata& oldlockdata, obslockdata newlockdata, uint32_t stationID, string satelliteID, GOBS signaltype) {
	bool found = _stmap[stationID][satelliteID].find(signaltype) != _stmap[stationID][satelliteID].end();
	if (found) {
		//extract data
		oldlockdata.epochtime = _stmap[stationID][satelliteID][signaltype].epochtime;
		oldlockdata.locktime = _stmap[stationID][satelliteID][signaltype].locktime;
		oldlockdata.suppcoeff = _stmap[stationID][satelliteID][signaltype].suppcoeff;
		_stmap[stationID][satelliteID].erase(signaltype);//delete read data from the map
	}
	_stmap[stationID][satelliteID][signaltype] = { newlockdata.epochtime, newlockdata.locktime, newlockdata.suppcoeff };//store new data
	return found;
}

//Loss-of-Lock-Indicator computation, using RTCMv3.2 standard-suggested algorithm (as epoch-time pass number of milliseconds in current week)
int t_rtcm3::msmgetLLI(uint32_t stationID, string satelliteID, GOBS signaltype, uint32_t epochtime, obslockdata lockdata) {
	lockdata.epochtime = epochtime;
	obslockdata oldlockdata = { 0, 0, 0 };
	if (!msmlligetprevtimes(oldlockdata, lockdata, stationID, (string)satelliteID, signaltype))return 1;//no past data => unknown => cycle slip expected by default

	// conversion of times to signed longs
	int64_t setime = epochtime;
	int64_t soetime = oldlockdata.epochtime;
	//variable naming @ RTCM 3.2 Amends-1-2 standard, section-3.1.2.1.1
	int64_t p = oldlockdata.locktime;
	int64_t n = lockdata.locktime;
	int64_t dt = ((setime - soetime) < 0) ? (setime - soetime) + 604800000 : (setime - soetime);//[ms]
	int64_t a = oldlockdata.suppcoeff;
	int64_t b = lockdata.suppcoeff;

	int lli = 1;//Loss of Lock presumed by default
	if (b > 0 && a > 0) {
		//high precision
		if (((p == n) && (dt < a)) ||
			((p < n) && (b > p) && (dt <= n)) || //this & next condition <=> ((p<n)&&(dt<=n)), left like this to improve transparency of standard application
			((p < n) && (b <= p) && (dt <= n))) lli = 0;//tracking continuous
	}
	else {
		//low precision
		if (((p == n) && (dt < p)) ||
			((p < n) && (dt <= n))) lli = 0;//tracking continuous
	}
	return lli;//LOC
}

// get satellite section size (datafields) for MSM type
int t_rtcm3::msmsatsz(int msmtype) { //msmtype is in range (1 .. 7)
	int sz = 1;
	switch (msmtype) {
	case 1:
		sz = 1; break;
	case 2:
		sz = 1; break;
	case 3:
		sz = 1; break;
	case 4:
		sz = 2; break;
	case 5:
		sz = 4; break;
	case 6:
		sz = 2; break;
	case 7:
		sz = 4; break;
	}
	return sz;
}
// get signal section size (datafields) for MSM type
int t_rtcm3::msmsigsz(int msmtype) { //msmtype is in range (1 .. 7)
	int sz = 1;
	switch (msmtype) {
	case 1:
		sz = 1; break;
	case 2:
		sz = 3; break;
	case 3:
		sz = 4; break;
	case 4:
		sz = 5; break;
	case 5:
		sz = 6; break;
	case 6:
		sz = 5; break;
	case 7:
		sz = 6; break;
	}
	return sz;
}

//field-reading function for InternalLooping
int t_rtcm3::msmgetiloopfield(int32_t* resultdatafield, int& srcstartbit, int srcfieldnumber, int loopsize, int srcfieldsize, bool issigned) {
	for (int i = 0; i < loopsize; i++) {
		resultdatafield[i + srcfieldnumber*loopsize] = (issigned) ? getuint(srcstartbit, srcfieldsize) : getuint(srcstartbit, srcfieldsize);
		srcstartbit += srcfieldsize;
	}
	return srcstartbit;
}

//MSM satellite data
int t_rtcm3::msgMSMsat123(int32_t* datafield, int& nsat, int startbit) {
	msmgetiloopfield(datafield, startbit, 0, nsat, 10);//GNSS Satellite rough ranges modulo 1ms
	return startbit;
}

int t_rtcm3::msgMSMsat46(int32_t* datafield, int& nsat, int startbit) {
	msmgetiloopfield(datafield, startbit, 0, nsat, 8);//number of integer milliseconds in GNSS satellite rough ranges
	msmgetiloopfield(datafield, startbit, 1, nsat, 10);//GNSS Satellite rough ranges modulo 1ms
	return startbit;
}

int t_rtcm3::msgMSMsat57(int32_t* datafield, int& nsat, int startbit) {
	msmgetiloopfield(datafield, startbit, 0, nsat, 8);//number of integer milliseconds in GNSS satellite rough ranges
	msmgetiloopfield(datafield, startbit, 1, nsat, 4);//Extended Satellite Information
	msmgetiloopfield(datafield, startbit, 2, nsat, 10);//GNSS Satellite rough ranges modulo 1ms
	msmgetiloopfield(datafield, startbit, 3, nsat, 14, true);//GNSS Satellite rough PhaseRangeRates
	return startbit;
}

//MSM signal data
int t_rtcm3::msgMSMsig1(int32_t* datafield, int& ncell, int startbit) {
	msmgetiloopfield(datafield, startbit, 0, ncell, 15, true);//GNSS Satellite fine pseudorange
	return startbit;
}

int t_rtcm3::msgMSMsig2(int32_t* datafield, int& ncell, int startbit) {
	msmgetiloopfield(datafield, startbit, 0, ncell, 22, true);//GNSS fine phase range data
	msmgetiloopfield(datafield, startbit, 1, ncell, 4);//PhaseRange Lock Time Indicator
	msmgetiloopfield(datafield, startbit, 2, ncell, 1);//Half-cycle ambiguity indicator
	return startbit;
}

int t_rtcm3::msgMSMsig3(int32_t* datafield, int& ncell, int startbit) {
	msmgetiloopfield(datafield, startbit, 0, ncell, 15, true);//GNSS Satellite fine pseudorange
	msmgetiloopfield(datafield, startbit, 1, ncell, 22, true);//GNSS fine phase range data
	msmgetiloopfield(datafield, startbit, 2, ncell, 4);//PhaseRange Lock Time Indicator
	msmgetiloopfield(datafield, startbit, 3, ncell, 1);//Half-cycle ambiguity indicator
	return startbit;
}

int t_rtcm3::msgMSMsig4(int32_t* datafield, int& ncell, int startbit) {
	msmgetiloopfield(datafield, startbit, 0, ncell, 15, true);//GNSS Satellite fine pseudorange
	msmgetiloopfield(datafield, startbit, 1, ncell, 22, true);//GNSS fine phase range data
	msmgetiloopfield(datafield, startbit, 2, ncell, 4);//PhaseRange Lock Time Indicator
	msmgetiloopfield(datafield, startbit, 3, ncell, 1);//Half-cycle ambiguity indicator
	msmgetiloopfield(datafield, startbit, 4, ncell, 6);//GNSS signal CNRs
	return startbit;
}

int t_rtcm3::msgMSMsig5(int32_t* datafield, int& ncell, int startbit) {
	msmgetiloopfield(datafield, startbit, 0, ncell, 15, true);//GNSS Satellite fine pseudorange
	msmgetiloopfield(datafield, startbit, 1, ncell, 22, true);//GNSS fine phase range data
	msmgetiloopfield(datafield, startbit, 2, ncell, 4);//PhaseRange Lock Time Indicator
	msmgetiloopfield(datafield, startbit, 3, ncell, 1);//Half-cycle ambiguity indicator
	msmgetiloopfield(datafield, startbit, 4, ncell, 6);//GNSS signal CNRs
	msmgetiloopfield(datafield, startbit, 5, ncell, 15, true);//GNSS signal fine PhaseRangeRates
	return startbit;
}

int t_rtcm3::msgMSMsig6(int32_t* datafield, int& ncell, int startbit) {
	msmgetiloopfield(datafield, startbit, 0, ncell, 20, true);//GNSS extended resolution pseudorange
	msmgetiloopfield(datafield, startbit, 1, ncell, 24, true);//GNSS extended resolution PhaseRange data
	msmgetiloopfield(datafield, startbit, 2, ncell, 10);//PhaseRange Lock Time Indicator with extended range and resolution
	msmgetiloopfield(datafield, startbit, 3, ncell, 1);//Half-cycle ambiguity indicator
	msmgetiloopfield(datafield, startbit, 4, ncell, 10);//GNSS signal CNRs with extended resolution
	return startbit;
}

int t_rtcm3::msgMSMsig7(int32_t* datafield, int& ncell, int startbit) {
	msmgetiloopfield(datafield, startbit, 0, ncell, 20, true);//GNSS extended resolution pseudorange
	msmgetiloopfield(datafield, startbit, 1, ncell, 24, true);//GNSS extended resolution PhaseRange data
	msmgetiloopfield(datafield, startbit, 2, ncell, 10);//PhaseRange Lock Time Indicator with extended range and resolution
	msmgetiloopfield(datafield, startbit, 3, ncell, 1);//Half-cycle ambiguity indicator
	msmgetiloopfield(datafield, startbit, 4, ncell, 10);//GNSS signal CNRs with extended resolution
	msmgetiloopfield(datafield, startbit, 5, ncell, 15, true);//GNSS signal fine PhaseRangeRates
	return startbit;
}

//MSM type decoding
int t_rtcm3::msgMSM1(uint64_t* headdata, int& nsat, char gnsssys) {
	int msmtype = 1;
	//read head, and calculate amount of data
	int nsig = 0;//number of signals for each satellite
	int position = msgMSMh(headdata, nsat, nsig);
	if (position < 0) { return 0; }//invalid cell-mask size
	int ncell = 0;//number of cells
	for (int i = 0; i < 64; i++) {
		ncell = (((headdata[11] >> i)&(uint64_t)1) != 0) ? ncell + 1 : ncell;
	}

	//read satellite data
	int32_t* satdata = new int32_t[nsat*msmsatsz(msmtype)];
	position = msgMSMsat123(satdata, nsat, position);

	//read signal data
	int32_t* sigdata = new int32_t[ncell*msmsigsz(msmtype)];
	position = msgMSMsig1(sigdata, ncell, position);

	//_______________________________________________________
	//interpret all available data, 1 satellite  => 1 t_gobsgnss record

	int satmaskslot = 0;//indicates position of tested bit within satellite-mask, determines satellite ID
	int sigmaskslot = 0;//indicates position of tested bit within signal-mask, determines signal type
	int cellmaskslot = 0;//indicates position of tested bit within cell-mask

	int obsid = 0;//indicates position of satellite information within satdata field
	int cellid = 0;//indicates position of signal information within sigdata field

	//create an observation data structure for every satellite
	for (satmaskslot = 0; satmaskslot < 64; satmaskslot++) {//satmaskslot is an index of satellite in RTCM3 MSM satellite mask
		if ((headdata[9] & ((uint64_t)1 << (63 - satmaskslot))) != 0) {//if n-th position of satellite-mask is set (MSB is first..LSB last)...

		  //create/find t_gobsgnss and setup signal-type tables for appropriate satellite system
			const GOBS* signalmapC = NULL;
			const GOBS* signalmapL = NULL;
			const GOBS* signalmapD = NULL;
			//const GOBS* signalmapS = NULL;
			shared_ptr<t_gobsgnss> nullsptrgobs;
			shared_ptr<t_gobsgnss> gobs = msmgetgobsgnss(headdata, satmaskslot, gnsssys, signalmapC, signalmapL, signalmapD/*, signalmapS*/);
			if (gobs == nullsptrgobs) {//unrecognized satellite
				obsid++;// this is important when counting data fields properly later
				continue;//NOTE: unknown satellites are skipped
			}

			//add all related signal data
			for (sigmaskslot = 0; sigmaskslot < 32; sigmaskslot++) {//sigmaskslot is an index of this signal from RTCM3 signal mask
				if ((headdata[10] & ((uint32_t)1 << (31 - sigmaskslot))) != 0) {//if signal is available in this message
					if ((headdata[11] & ((uint64_t)1 << (63 - cellmaskslot))) != 0) {//if signal is available for this satellite
					  /*
					  //process signals' data and add results into t_gobsgnss structure
					  GOBS signalcode;
					  double obsvalue;

					  //__ADD__ALL__OBSERVABLES__AVAILABLE__FOR__THE__SIGNAL_____________________________________________________

					  //add an observable
					  signalcode = signalmapC[sigmaskslot];
					  int roughrangefieldnumber = obsid*msmsatsz(msmtype)+0;//obsID*numberOfFieldsPerSatellite + numberOfRequiredFieldForSatellite
					  int finerangefieldnumber = cellid*msmsigsz(msmtype)+0;//cellID*numberOfFieldsPerSignal + numberOfRequiredFieldForSignal
					  obsvalue = msmgetPS(0,satdata[roughrangefieldnumber],sigdata[finerangefieldnumber]);//NOTE: the zero is a number of milliseconds since the beginning of the epoch
					  gobs->addobs(signalcode, obsvalue);
					  */ //NOTE: MSM1, MSM2 & MSM3 messages do contain only modulo 1ms values for observation times
					  //       "=> absolute value cannot be calculated and reading message is dropped
				   //__ALL__OBSERVABLES__FOR__SIGNAL__ADDED___________________________________________________________________

						cellid++;//move on to the next field for reading a signal
					}
					cellmaskslot++;//move on to resolving the next signal
				}
			}
			obsid++;//move to record of the next satellite
		}
	}

	//cleanup
	delete[] satdata;
	delete[] sigdata;

	return 1;
}

int t_rtcm3::msgMSM2(uint64_t* headdata, int& nsat, char gnsssys) {
	int msmtype = 2;
	//read head, and calculate amount of data
	int nsig = 0;
	int position = msgMSMh(headdata, nsat, nsig);
	if (position < 0) { return 0; }//invalid cell-mask size
	int ncell = 0;
	for (int i = 0; i < 64; i++) {
		ncell = (((headdata[11] >> i)&(uint64_t)1) != 0) ? ncell + 1 : ncell;
	}

	//read satellite data
	int32_t* satdata = new int32_t[nsat*msmsatsz(msmtype)];
	position = msgMSMsat123(satdata, nsat, position);

	//read signal data
	int32_t* sigdata = new int32_t[ncell*msmsigsz(msmtype)];
	position = msgMSMsig2(sigdata, ncell, position);

	//_______________________________________________________
	//interpret all available data, 1 satellite  => 1 t_gobsgnss record

	int satmaskslot = 0;//indicates position of tested bit within satellite-mask, determines satellite ID
	int sigmaskslot = 0;//indicates position of tested bit within signal-mask, determines signal type
	int cellmaskslot = 0;//indicates position of tested bit within cell-mask

	int obsid = 0;//indicates position of satellite information within satdata field
	int cellid = 0;//indicates position of signal information within sigdata field

	//create an observation data structure for every satellite
	for (satmaskslot = 0; satmaskslot < 64; satmaskslot++) {//sat_i is an index of satellite in RTCM3 MSM satellite mask
		if ((headdata[9] & ((uint64_t)1 << (63 - satmaskslot))) != 0) {//if n-th position of satellite-mask is set (MSB is first..LSB last)...

		  //create t_gobsgnss and setup signal-type tables for appropriate satellite system
			const GOBS* signalmapC = NULL;
			const GOBS* signalmapL = NULL;
			const GOBS* signalmapD = NULL;
			//const GOBS* signalmapS;
			shared_ptr<t_gobsgnss> nullsptrgobs;
			shared_ptr<t_gobsgnss> gobs = msmgetgobsgnss(headdata, satmaskslot, gnsssys, signalmapC, signalmapL, signalmapD/*, signalmapS*/);
			if (gobs == nullsptrgobs) {//unrecognized satellite
				obsid++;// this is important when counting data fields properly later
				continue;//NOTE: unknown satellites are skipped
			}

			//add all related signal data
			for (sigmaskslot = 0; sigmaskslot < 32; sigmaskslot++) {//sig_i is an index of this signal from RTCM3 signal mask
				if ((headdata[10] & ((uint32_t)1 << (31 - sigmaskslot))) != 0) {//if signal is available in this message
					if ((headdata[11] & ((uint64_t)1 << (63 - cellmaskslot))) != 0) {//if signal is available for this satellite
					  /*
					  //process signals' data and add results into t_gobsgnss structure
					  GOBS signalcode;
					  double obsvalue;

					  //__ADD__ALL__OBSERVABLES__AVAILABLE__FOR__THE__SIGNAL_____________________________________________________

					  //add an observable
					  signalcode = signalmapL[sigmaskslot];
					  obsvalue = msmgetPH(0,satdata[obsid*msmsatsz(msmtype)+0],sigdata[cellid*msmsigsz(msmtype)+0],gobs->frequency(t_gobs(signalcode).band()));
					  gobs->addobs(signalcode, obsvalue);

					  //add LLI
					  //compute number of milliseconds elapsed in this epoch from the head data to use for LLI detection algorithm
					  uint32_t epochtime_ms;
					  if(gnsssys != 'R'){ //other systems (30bit epoch time in header represents number of milliseconds from it's beginning)
						epochtime_ms = (double)headdata[1];
					  }else{ //GLONASS, (MSB> 3bit day of week; 27bit milliseconds of the day <LSB)
						epochtime_ms = ((double)((headdata[1] >> 27) & 7) * 86400000.0) + (double)(headdata[1] & 0x7FFFFFF);
					  }
					  obsvalue = msmgetLLI(headdata[0],gobs->sat(),signalcode,epochtime_ms,msmgetMLT(sigdata[cellid*msmsigsz(msmtype)+1]));
					  obsvalue = ((int)sigdata[cellid*msmsigsz(msmtype)+2] == 0)? obsvalue : 1;//NOTE: consider half-cycle ambiguity flag; if this is set consider it LossOfLock
					  gobs->addlli(signalcode,(int)obsvalue);
					  */ //NOTE: MSM1, MSM2 & MSM3 messages do contain only modulo 1ms values for observation times
					  //__ALL__OBSERVABLES__FOR__SIGNAL__ADDED___________________________________________________________________

						cellid++;//move on to the next field for reading a signal
					}
					cellmaskslot++;//move on to resolving the next signal
				}
			}
			obsid++;//move to record of the next satellite
		}
	}

	//cleanup
	delete[] satdata;
	delete[] sigdata;

	return 1;
}

int t_rtcm3::msgMSM3(uint64_t* headdata, int& nsat, char gnsssys) {
	int msmtype = 3;
	//read head, and calculate amount of data
	int nsig = 0;
	int position = msgMSMh(headdata, nsat, nsig);
	if (position < 0) { return 0; }//invalid cell-mask size
	int ncell = 0;
	for (int i = 0; i < 64; i++) {
		ncell = (((headdata[11] >> i)&(uint64_t)1) != 0) ? ncell + 1 : ncell;
	}

	//read satellite data
	int32_t* satdata = new int32_t[nsat*msmsatsz(msmtype)];
	position = msgMSMsat123(satdata, nsat, position);

	//read signal data
	int32_t* sigdata = new int32_t[ncell*msmsigsz(msmtype)];
	position = msgMSMsig3(sigdata, ncell, position);

	//_______________________________________________________
	//interpret all available data, 1 satellite  => 1 t_gobsgnss record

	int satmaskslot = 0;//indicates position of tested bit within satellite-mask, determines satellite ID
	int sigmaskslot = 0;//indicates position of tested bit within signal-mask, determines signal type
	int cellmaskslot = 0;//indicates position of tested bit within cell-mask

	int obsid = 0;//indicates position of satellite information within satdata field, also determines position of result data within the datafield array
	int cellid = 0;//indicates position of signal information within sigdata field

	//create an observation data structure for every satellite
	for (satmaskslot = 0; satmaskslot < 64; satmaskslot++) {//sat_i is an index of satellite in RTCM3 MSM satellite mask
		if ((headdata[9] & ((uint64_t)1 << (63 - satmaskslot))) != 0) {//if n-th position of satellite-mask is set (MSB is first..LSB last)...

		  //create t_gobsgnss and setup signal-type tables for appropriate satellite system
			const GOBS* signalmapC = NULL;
			const GOBS* signalmapL = NULL;
			const GOBS* signalmapD = NULL;
			//const GOBS* signalmapS;
			shared_ptr<t_gobsgnss> nullsptrgobs;
			shared_ptr<t_gobsgnss> gobs = msmgetgobsgnss(headdata, satmaskslot, gnsssys, signalmapC, signalmapL, signalmapD/*, signalmapS*/);
			if (gobs == nullsptrgobs) {//unrecognized satellite
				obsid++;// this is important when counting data fields properly later
				continue;//NOTE: unknown satellites are skipped
			}

			//add all related signal data
			for (sigmaskslot = 0; sigmaskslot < 32; sigmaskslot++) {//sig_i is an index of this signal from RTCM3 signal mask
				if ((headdata[10] & ((uint32_t)1 << (31 - sigmaskslot))) != 0) {//if signal is available in this message
					if ((headdata[11] & ((uint64_t)1 << (63 - cellmaskslot))) != 0) {//if signal is available for this satellite
					  /*
					  //process signals' data and add results into t_gobsgnss structure
					  GOBS signalcode;
					  double obsvalue;

					  //__ADD__ALL__OBSERVABLES__AVAILABLE__FOR__THE__SIGNAL_____________________________________________________

					  //add an observable (Pseudorange)
					  signalcode = signalmapC[sigmaskslot];
					  int roughrangefieldnumber = obsid*msmsatsz(msmtype)+0;//obsID*numberOfFieldsPerSatellite + numberOfRequiredFieldForSatellite
					  int finerangefieldnumber = cellid*msmsigsz(msmtype)+0;//cellID*numberOfFieldsPerSignal + numberOfRequiredFieldForSignal
					  obsvalue = msmgetPS(0,satdata[roughrangefieldnumber],sigdata[finerangefieldnumber]);
					  gobs->addobs(signalcode, obsvalue);

					  //add an observable (Phaserange)
					  signalcode = signalmapL[sigmaskslot];
					  obsvalue = msmgetPH(0,satdata[obsid*msmsatsz(msmtype)+0],sigdata[cellid*msmsigsz(msmtype)+1],gobs->frequency(t_gobs(signalcode).band()));
					  gobs->addobs(signalcode, obsvalue);

					  //add LLI
					  //compute number of milliseconds elapsed in this epoch from the head data to use for LLI detection algorithm
					  uint32_t epochtime_ms;
					  if(gnsssys != 'R'){ //other systems (30bit epoch time in header represents number of milliseconds from it's beginning)
						epochtime_ms = (double)headdata[1];
					  }else{ //GLONASS, (MSB> 3bit day of week; 27bit milliseconds of the day <LSB)
						epochtime_ms = ((double)((headdata[1] >> 27) & 7) * 86400000.0) + (double)(headdata[1] & 0x7FFFFFF);
					  }
					  obsvalue = msmgetLLI(headdata[0],gobs->sat(),signalcode,epochtime_ms,msmgetMLT(sigdata[cellid*msmsigsz(msmtype)+2]));
					  obsvalue = ((int)sigdata[cellid*msmsigsz(msmtype)+3] == 0)? obsvalue : 1;//NOTE: consider half-cycle ambiguity flag; if this is set consider it LossOfLock
					  gobs->addlli(signalcode,(int)obsvalue);
					  */ //NOTE: MSM1, MSM2 & MSM3 messages do contain only modulo 1ms values for observation times
					  //__ALL__OBSERVABLES__FOR__SIGNAL__ADDED___________________________________________________________________

						cellid++;//move on to the next field for reading a signal
					}
					cellmaskslot++;//move on to resolving the next signal
				}
			}
			obsid++;//move to record of the next satellite
		}
	}

	//cleanup
	delete[] satdata;
	delete[] sigdata;

	return 1;
}

int t_rtcm3::msgMSM4(uint64_t* headdata, int& nsat, char gnsssys) {
	int msmtype = 4;
	//read head, and calculate amount of data
	int nsig = 0;//number of signals for each satellite
	int position = msgMSMh(headdata, nsat, nsig);
	if (position < 0) { return 0; }//invalid cell-mask size
	int ncell = 0;//number of cells
	for (int i = 0; i < 64; i++) {
		ncell = (((headdata[11] >> i)&(uint64_t)1) != 0) ? ncell + 1 : ncell;
	}

	//read satellite data
	int32_t* satdata = new int32_t[nsat*msmsatsz(msmtype)];
	position = msgMSMsat46(satdata, nsat, position);

	//read signal data
	int32_t* sigdata = new int32_t[ncell*msmsigsz(msmtype)];
	position = msgMSMsig4(sigdata, ncell, position);

	//_______________________________________________________
	//interpret all available data, 1 satellite  => 1 t_gobsgnss record

	int satmaskslot = 0;//indicates position of tested bit within satellite-mask, determines satellite ID
	int sigmaskslot = 0;//indicates position of tested bit within signal-mask, determines signal type
	int cellmaskslot = 0;//indicates position of tested bit within cell-mask

	int obsid = 0;//indicates position of satellite information within satdata field
	int cellid = 0;//indicates position of signal information within sigdata field

	//create an observation data structure for every satellite
	for (satmaskslot = 0; satmaskslot < 64; satmaskslot++) {//satmaskslot is an index of satellite in RTCM3 MSM satellite mask
		if ((headdata[9] & ((uint64_t)1 << (63 - satmaskslot))) != 0) {//if n-th position of satellite-mask is set (MSB is first..LSB last)...

		  //create/find t_gobsgnss and setup signal-type tables for appropriate satellite system
			const GOBS* signalmapC = NULL;
			const GOBS* signalmapL = NULL;
			const GOBS* signalmapD = NULL;
			//const GOBS* signalmapS = NULL;
			shared_ptr<t_gobsgnss> nullsptrgobs;
			shared_ptr<t_gobsgnss> gobs = msmgetgobsgnss(headdata, satmaskslot, gnsssys, signalmapC, signalmapL, signalmapD/*, signalmapS*/);
			if (gobs == nullsptrgobs) {//unrecognized satellite
				obsid++;// this is important when counting data fields properly later
				continue;//NOTE: unknown satellites are skipped
			}

			//add all related signal data
			for (sigmaskslot = 0; sigmaskslot < 32; sigmaskslot++) {//sigmaskslot is an index of this signal from RTCM3 signal mask
				if ((headdata[10] & ((uint32_t)1 << (31 - sigmaskslot))) != 0) {//if signal is available in this message
					if ((headdata[11] & ((uint64_t)1 << (63 - cellmaskslot))) != 0) {//if signal is available for this satellite
					  //process signals' data and add results into t_gobsgnss structure
						GOBS signalcode;
						double obsvalue;

						//__ADD__ALL__OBSERVABLES__AVAILABLE__FOR__THE__SIGNAL_____________________________________________________

						//add an observable (Pseudorange)
						signalcode = signalmapC[sigmaskslot];
						int obsmsfieldnumber = obsid*msmsatsz(msmtype) + 0;     //integer number of milliseconds in observation
						int roughrangefieldnumber = obsid*msmsatsz(msmtype) + 1;//obsID*numberOfFieldsPerSatellite + numberOfRequiredFieldForSatellite
						int finerangefieldnumber = cellid*msmsigsz(msmtype) + 0;//cellID*numberOfFieldsPerSignal + numberOfRequiredFieldForSignal
						obsvalue = msmgetPS(satdata[obsmsfieldnumber], satdata[roughrangefieldnumber], sigdata[finerangefieldnumber]);
						if ((int64_t)(gobs->getobs(signalcode)) == 0)gobs->addobs(signalcode, obsvalue);// add this observation if it's (int64_t)gobs_value == 0 (unknown)

						//add an observable (Phaserange)
						signalcode = signalmapL[sigmaskslot];
						int finephasefieldnumber = cellid*msmsigsz(msmtype) + 1;//cellID*numberOfFieldsPerSignal + numberOfRequiredFieldForSignal
						obsvalue = msmgetPH(satdata[obsmsfieldnumber], satdata[roughrangefieldnumber], sigdata[finephasefieldnumber], gobs->frequency(t_gobs(signalcode).band()));
						if ((int64_t)(gobs->getobs(signalcode)) == 0)gobs->addobs(signalcode, obsvalue);// add this observation if it's (int64_t)gobs_value == 0 (unknown)

						//add LLI
						//compute number of milliseconds elapsed in this epoch from the head data to use for LLI detection algorithm
						uint32_t epochtime_ms;
						if (gnsssys != 'R') { //other systems (30bit epoch time in header represents number of milliseconds from it's beginning)
							epochtime_ms = (double)headdata[1];
						}
						else { //GLONASS, (MSB> 3bit day of week; 27bit milliseconds of the day <LSB)
							epochtime_ms = ((double)((headdata[1] >> 27) & 7) * 86400000.0) + (double)(headdata[1] & 0x7FFFFFF);
						}
						obsvalue = msmgetLLI(headdata[0], gobs->sat(), signalcode, epochtime_ms, msmgetMLT(sigdata[cellid*msmsigsz(msmtype) + 2]));
						obsvalue = ((int)sigdata[cellid*msmsigsz(msmtype) + 3] == 0) ? obsvalue : 1;//NOTE: consider half-cycle ambiguity flag; if this is set consider it LossOfLock
						gobs->addlli(signalcode, (int)obsvalue);// add this observation, NOTE: no decent way to distinguish between HIprecision and normal result

						//add CNRs
						//NOTE: data decoded, signal strength hand-over to t_gobsgnss to be implemented here

						//__ALL__OBSERVABLES__FOR__SIGNAL__ADDED___________________________________________________________________

						cellid++;//move on to the next field for reading a signal
					}
					cellmaskslot++;//move on to resolving the next signal
				}
			}
			obsid++;//move to record of the next satellite
		}
	}

	//cleanup
	delete[] satdata;
	delete[] sigdata;

	return 1;
}

int t_rtcm3::msgMSM5(uint64_t* headdata, int& nsat, char gnsssys) {
	int msmtype = 5;
	//read head, and calculate amount of data
	int nsig = 0;//number of signals for each satellite
	int position = msgMSMh(headdata, nsat, nsig);
	if (position < 0) { return 0; }//invalid cell-mask size
	int ncell = 0;//number of cells
	for (int i = 0; i < 64; i++) {
		ncell = (((headdata[11] >> i)&(uint64_t)1) != 0) ? ncell + 1 : ncell;
	}

	//read satellite data
	int32_t* satdata = new int32_t[nsat*msmsatsz(msmtype)];
	position = msgMSMsat57(satdata, nsat, position);

	//read signal data
	int32_t* sigdata = new int32_t[ncell*msmsigsz(msmtype)];
	position = msgMSMsig5(sigdata, ncell, position);

	//_______________________________________________________
	//interpret all available data, 1 satellite  => 1 t_gobsgnss record

	int satmaskslot = 0;//indicates position of tested bit within satellite-mask, determines satellite ID
	int sigmaskslot = 0;//indicates position of tested bit within signal-mask, determines signal type
	int cellmaskslot = 0;//indicates position of tested bit within cell-mask

	int obsid = 0;//indicates position of satellite information within satdata field
	int cellid = 0;//indicates position of signal information within sigdata field

	//create an observation data structure for every satellite
	for (satmaskslot = 0; satmaskslot < 64; satmaskslot++) {//satmaskslot is an index of satellite in RTCM3 MSM satellite mask
		if ((headdata[9] & ((uint64_t)1 << (63 - satmaskslot))) != 0) {//if n-th position of satellite-mask is set (MSB is first..LSB last)...

		  //create/find t_gobsgnss and setup signal-type tables for appropriate satellite system
			const GOBS* signalmapC = NULL;
			const GOBS* signalmapL = NULL;
			const GOBS* signalmapD = NULL;
			//const GOBS* signalmapS = NULL;
			shared_ptr<t_gobsgnss> nullsptrgobs;
			shared_ptr<t_gobsgnss> gobs = msmgetgobsgnss(headdata, satmaskslot, gnsssys, signalmapC, signalmapL, signalmapD/*, signalmapS*/);
			if (gobs == nullsptrgobs) {//unrecognized satellite
				obsid++;// this is important when counting data fields properly later
				continue;//NOTE: unknown satellites are skipped
			}

			// add channel for GLONASS sat
			if (gnsssys == 'R') {
				uint32_t chindicator = satdata[obsid*msmsatsz(msmtype) + 1];
				if (chindicator < 14) {
					gobs->channel(chindicator - 7);//available only for MSM5 & MSM7
				}//NOTE: 14 is reserved, 15 is unknown or invalid value
			}

			//add all related signal data
			for (sigmaskslot = 0; sigmaskslot < 32; sigmaskslot++) {//sigmaskslot is an index of this signal from RTCM3 signal mask
				if ((headdata[10] & ((uint32_t)1 << (31 - sigmaskslot))) != 0) {//if signal is available in this message
					if ((headdata[11] & ((uint64_t)1 << (63 - cellmaskslot))) != 0) {//if signal is available for this satellite
					  //process signals' data and add results into t_gobsgnss structure
						GOBS signalcode;
						double obsvalue;

						//__ADD__ALL__OBSERVABLES__AVAILABLE__FOR__THE__SIGNAL_____________________________________________________

						//add an observable (Pseudorange)
						signalcode = signalmapC[sigmaskslot];
						int obsmsfieldnumber = obsid*msmsatsz(msmtype) + 0;     //integer number of milliseconds in observation
						int roughrangefieldnumber = obsid*msmsatsz(msmtype) + 2;//obsID*numberOfFieldsPerSatellite + numberOfRequiredFieldForSatellite
						int finerangefieldnumber = cellid*msmsigsz(msmtype) + 0;//cellID*numberOfFieldsPerSignal + numberOfRequiredFieldForSignal
						obsvalue = msmgetPS(satdata[obsmsfieldnumber], satdata[roughrangefieldnumber], sigdata[finerangefieldnumber]);
						if ((int64_t)(gobs->getobs(signalcode)) == 0)gobs->addobs(signalcode, obsvalue);// add this observation if it's (int64_t)gobs_value == 0 (unknown)

						//add an observable (Phaserange)
						signalcode = signalmapL[sigmaskslot];
						int finephasefieldnumber = cellid*msmsigsz(msmtype) + 1;//cellID*numberOfFieldsPerSignal + numberOfRequiredFieldForSignal
						obsvalue = msmgetPH(satdata[obsmsfieldnumber], satdata[roughrangefieldnumber], sigdata[finephasefieldnumber], gobs->frequency(t_gobs(signalcode).band()));
						if ((int64_t)(gobs->getobs(signalcode)) == 0)gobs->addobs(signalcode, obsvalue);// add this observation if it's (int64_t)gobs_value == 0 (unknown)

						//add LLI
						//compute number of milliseconds elapsed in this epoch from the head data to use for LLI detection algorithm
						uint32_t epochtime_ms;
						if (gnsssys != 'R') { //other systems (30bit epoch time in header represents number of milliseconds from it's beginning)
							epochtime_ms = (double)headdata[1];
						}
						else { //GLONASS, (MSB> 3bit day of week; 27bit milliseconds of the day <LSB)
							epochtime_ms = ((double)((headdata[1] >> 27) & 7) * 86400000.0) + (double)(headdata[1] & 0x7FFFFFF);
						}
						obsvalue = msmgetLLI(headdata[0], gobs->sat(), signalcode, epochtime_ms, msmgetMLT(sigdata[cellid*msmsigsz(msmtype) + 2]));
						obsvalue = ((int)sigdata[cellid*msmsigsz(msmtype) + 3] == 0) ? obsvalue : 1;//NOTE: consider half-cycle ambiguity flag; if this is set consider it LossOfLock
						gobs->addlli(signalcode, (int)obsvalue);// add this observation, NOTE: no decent way to distinguish between HIprecision and normal result

						//add CNRs
						//NOTE: data decoded, signal strength hand-over to t_gobsgnss to be implemented here

						//add an observable (PhaseRangeRate/DopplerShift)
						signalcode = signalmapD[sigmaskslot];
						int roughdopplerfieldnumber = obsid*msmsatsz(msmtype) + 3;
						int finedopplerfieldnumber = cellid*msmsigsz(msmtype) + 5;
						obsvalue = msmgetDS(satdata[roughdopplerfieldnumber], sigdata[finedopplerfieldnumber], gobs->frequency(t_gobs(signalcode).band()));
						if ((int64_t)(gobs->getobs(signalcode)) == 0)gobs->addobs(signalcode, obsvalue);// add this observation if it's (int64_t)gobs_value == 0 (unknown)

						//__ALL__OBSERVABLES__FOR__SIGNAL__ADDED___________________________________________________________________

						cellid++;//move on to the next field for reading a signal
					}
					cellmaskslot++;//move on to resolving the next signal
				}
			}
			obsid++;//move to record of the next satellite
		}
	}

	//cleanup
	delete[] satdata;
	delete[] sigdata;

	return 1;
}

int t_rtcm3::msgMSM6(uint64_t* headdata, int& nsat, char gnsssys) {
	int msmtype = 6;
	//read head, and calculate amount of data
	int nsig = 0;//number of signals for each satellite
	int position = msgMSMh(headdata, nsat, nsig);
	if (position < 0) { return 0; }//invalid cell-mask size
	int ncell = 0;//number of cells
	for (int i = 0; i < 64; i++) {
		ncell = (((headdata[11] >> i)&(uint64_t)1) != 0) ? ncell + 1 : ncell;
	}

	//read satellite data
	int32_t* satdata = new int32_t[nsat*msmsatsz(msmtype)];
	position = msgMSMsat46(satdata, nsat, position);

	//read signal data
	int32_t* sigdata = new int32_t[ncell*msmsigsz(msmtype)];
	position = msgMSMsig6(sigdata, ncell, position);

	//_______________________________________________________
	//interpret all available data, 1 satellite  => 1 t_gobsgnss record

	int satmaskslot = 0;//indicates position of tested bit within satellite-mask, determines satellite ID
	int sigmaskslot = 0;//indicates position of tested bit within signal-mask, determines signal type
	int cellmaskslot = 0;//indicates position of tested bit within cell-mask

	int obsid = 0;//indicates position of satellite information within satdata field
	int cellid = 0;//indicates position of signal information within sigdata field

	//create an observation data structure for every satellite
	for (satmaskslot = 0; satmaskslot < 64; satmaskslot++) {//satmaskslot is an index of satellite in RTCM3 MSM satellite mask
		if ((headdata[9] & ((uint64_t)1 << (63 - satmaskslot))) != 0) {//if n-th position of satellite-mask is set (MSB is first..LSB last)...

		  //create/find t_gobsgnss and setup signal-type tables for appropriate satellite system
			const GOBS* signalmapC = NULL;
			const GOBS* signalmapL = NULL;
			const GOBS* signalmapD = NULL;
			//const GOBS* signalmapS = NULL;
			shared_ptr<t_gobsgnss> nullsptrgobs;
			shared_ptr<t_gobsgnss> gobs = msmgetgobsgnss(headdata, satmaskslot, gnsssys, signalmapC, signalmapL, signalmapD/*, signalmapS*/);
			if (gobs == nullsptrgobs) {//unrecognized satellite
				obsid++;// this is important when counting data fields properly later
				continue;//NOTE: unknown satellites are skipped
			}

			//add all related signal data
			for (sigmaskslot = 0; sigmaskslot < 32; sigmaskslot++) {//sigmaskslot is an index of this signal from RTCM3 signal mask
				if ((headdata[10] & ((uint32_t)1 << (31 - sigmaskslot))) != 0) {//if signal is available in this message
					if ((headdata[11] & ((uint64_t)1 << (63 - cellmaskslot))) != 0) {//if signal is available for this satellite
					  //process signals' data and add results into t_gobsgnss structure
						GOBS signalcode;
						double obsvalue;

						//__ADD__ALL__OBSERVABLES__AVAILABLE__FOR__THE__SIGNAL_____________________________________________________

						//add an observable (Pseudorange)
						signalcode = signalmapC[sigmaskslot];
						int obsmsfieldnumber = obsid*msmsatsz(msmtype) + 0;     //integer number of milliseconds in observation
						int roughrangefieldnumber = obsid*msmsatsz(msmtype) + 1;//obsID*numberOfFieldsPerSatellite + numberOfRequiredFieldForSatellite
						int finerangefieldnumber = cellid*msmsigsz(msmtype) + 0;//cellID*numberOfFieldsPerSignal + numberOfRequiredFieldForSignal
						obsvalue = msmgetPS(satdata[obsmsfieldnumber], satdata[roughrangefieldnumber], sigdata[finerangefieldnumber], true);
						gobs->addobs(signalcode, obsvalue);// add this observation, HiPrecision

						//add an observable (Phaserange)
						signalcode = signalmapL[sigmaskslot];
						int finephasefieldnumber = cellid*msmsigsz(msmtype) + 1;//cellID*numberOfFieldsPerSignal + numberOfRequiredFieldForSignal
						obsvalue = msmgetPH(satdata[obsmsfieldnumber], satdata[roughrangefieldnumber], sigdata[finephasefieldnumber], gobs->frequency(t_gobs(signalcode).band()), true);
						gobs->addobs(signalcode, obsvalue);// add this observation, HiPrecision

						//add LLI
						//compute number of milliseconds elapsed in this epoch from the head data to use for LLI detection algorithm
						uint32_t epochtime_ms;
						if (gnsssys != 'R') { //other systems (30bit epoch time in header represents number of milliseconds from it's beginning)
							epochtime_ms = (double)headdata[1];
						}
						else { //GLONASS, (MSB> 3bit day of week; 27bit milliseconds of the day <LSB)
							epochtime_ms = ((double)((headdata[1] >> 27) & 7) * 86400000.0) + (double)(headdata[1] & 0x7FFFFFF);
						}
						obsvalue = msmgetLLI(headdata[0], gobs->sat(), signalcode, epochtime_ms, msmgetMLT(sigdata[cellid*msmsigsz(msmtype) + 2], true));
						obsvalue = ((int)sigdata[cellid*msmsigsz(msmtype) + 3] == 0) ? obsvalue : 1;//NOTE: consider half-cycle ambiguity flag; if this is set consider it LossOfLock
						gobs->addlli(signalcode, (int)obsvalue);// add this observation, HiPrecision

						//add CNRs
						//NOTE: data decoded, signal strength hand-over to t_gobsgnss to be implemented here

						//__ALL__OBSERVABLES__FOR__SIGNAL__ADDED___________________________________________________________________

						cellid++;//move on to the next field for reading a signal
					}
					cellmaskslot++;//move on to resolving the next signal
				}
			}
			obsid++;//move to record of the next satellite
		}
	}

	//cleanup
	delete[] satdata;
	delete[] sigdata;

	return 1;
}

int t_rtcm3::msgMSM7(uint64_t* headdata, int& nsat, char gnsssys) {
	int msmtype = 7;
	//read head, and calculate amount of data
	int nsig = 0;//number of signals for each satellite
	int position = msgMSMh(headdata, nsat, nsig);
	if (position < 0) { return 0; }//invalid cell-mask size
	int ncell = 0;//number of cells
	for (int i = 0; i < 64; i++) {
		ncell = (((headdata[11] >> i)&(uint64_t)1) != 0) ? ncell + 1 : ncell;
	}

	//read satellite data
	int32_t* satdata = new int32_t[nsat*msmsatsz(msmtype)];
	position = msgMSMsat57(satdata, nsat, position);

	//read signal data
	int32_t* sigdata = new int32_t[ncell*msmsigsz(msmtype)];
	position = msgMSMsig7(sigdata, ncell, position);

	//_______________________________________________________
	//interpret all available data, 1 satellite  => 1 t_gobsgnss record

	int satmaskslot = 0;//indicates position of tested bit within satellite-mask, determines satellite ID
	int sigmaskslot = 0;//indicates position of tested bit within signal-mask, determines signal type
	int cellmaskslot = 0;//indicates position of tested bit within cell-mask

	int obsid = 0;//indicates position of satellite information within satdata field
	int cellid = 0;//indicates position of signal information within sigdata field

	//create an observation data structure for every satellite
	for (satmaskslot = 0; satmaskslot < 64; satmaskslot++) {//satmaskslot is an index of satellite in RTCM3 MSM satellite mask
		if ((headdata[9] & ((uint64_t)1 << (63 - satmaskslot))) != 0) {//if n-th position of satellite-mask is set (MSB is first..LSB last)...

		  //create/find t_gobsgnss and setup signal-type tables for appropriate satellite system
			const GOBS* signalmapC = NULL;
			const GOBS* signalmapL = NULL;
			const GOBS* signalmapD = NULL;
			//const GOBS* signalmapS = NULL;
			shared_ptr<t_gobsgnss> nullsptrgobs;
			shared_ptr<t_gobsgnss> gobs = msmgetgobsgnss(headdata, satmaskslot, gnsssys, signalmapC, signalmapL, signalmapD/*, signalmapS*/);
			if (gobs == nullsptrgobs) {//unrecognized satellite
				obsid++;// this is important when counting data fields properly later
				continue;//NOTE: unknown satellites are skipped
			}

			// add channel for GLONASS sat
			if (gnsssys == 'R') {
				uint32_t chindicator = satdata[obsid*msmsatsz(msmtype) + 1];
				if (chindicator < 14) {
					gobs->channel(chindicator - 7);//available only for MSM5 & MSM7
				}//NOTE: 14 is reserved, 15 is unknown or invalid value
			}

			//add all related signal data
			for (sigmaskslot = 0; sigmaskslot < 32; sigmaskslot++) {//sigmaskslot is an index of this signal from RTCM3 signal mask
				if ((headdata[10] & ((uint32_t)1 << (31 - sigmaskslot))) != 0) {//if signal is available in this message
					if ((headdata[11] & ((uint64_t)1 << (63 - cellmaskslot))) != 0) {//if signal is available for this satellite
					  //process signals' data and add results into t_gobsgnss structure
						GOBS signalcode;
						double obsvalue;

						//__ADD__ALL__OBSERVABLES__AVAILABLE__FOR__THE__SIGNAL_____________________________________________________

						//add an observable (Pseudorange)
						signalcode = signalmapC[sigmaskslot];
						int obsmsfieldnumber = obsid*msmsatsz(msmtype) + 0;     //integer number of milliseconds in observation
						int roughrangefieldnumber = obsid*msmsatsz(msmtype) + 2;//obsID*numberOfFieldsPerSatellite + numberOfRequiredFieldForSatellite
						int finerangefieldnumber = cellid*msmsigsz(msmtype) + 0;//cellID*numberOfFieldsPerSignal + numberOfRequiredFieldForSignal
						obsvalue = msmgetPS(satdata[obsmsfieldnumber], satdata[roughrangefieldnumber], sigdata[finerangefieldnumber], true);
						gobs->addobs(signalcode, obsvalue);// add this observation, HiPrecision

						//add an observable (Phaserange)
						signalcode = signalmapL[sigmaskslot];
						int finephasefieldnumber = cellid*msmsigsz(msmtype) + 1;//cellID*numberOfFieldsPerSignal + numberOfRequiredFieldForSignal
						obsvalue = msmgetPH(satdata[obsmsfieldnumber], satdata[roughrangefieldnumber], sigdata[finephasefieldnumber], gobs->frequency(t_gobs(signalcode).band()), true);
						gobs->addobs(signalcode, obsvalue);// add this observation, HiPrecision

						//add LLI
						//compute number of milliseconds elapsed in this epoch from the head data to use for LLI detection algorithm
						uint32_t epochtime_ms;
						if (gnsssys != 'R') { //other systems (30bit epoch time in header represents number of milliseconds from it's beginning)
							epochtime_ms = (double)headdata[1];
						}
						else { //GLONASS, (MSB> 3bit day of week; 27bit milliseconds of the day <LSB)
							epochtime_ms = ((double)((headdata[1] >> 27) & 7) * 86400000.0) + (double)(headdata[1] & 0x7FFFFFF);
						}
						obsvalue = msmgetLLI(headdata[0], gobs->sat(), signalcode, epochtime_ms, msmgetMLT(sigdata[cellid*msmsigsz(msmtype) + 2], true));
						obsvalue = ((int)sigdata[cellid*msmsigsz(msmtype) + 3] == 0) ? obsvalue : 1;//NOTE: consider half-cycle ambiguity flag; if this is set consider it LossOfLock
						gobs->addlli(signalcode, (int)obsvalue);// add this observation, HiPrecision

						//add CNRs
						//NOTE: data decoded, signal strength hand-over to t_gobsgnss to be implemented here

						//add an observable (PhaseRangeRate/DopplerShift)
						signalcode = signalmapD[sigmaskslot];
						int roughdopplerfieldnumber = obsid*msmsatsz(msmtype) + 3;
						int finedopplerfieldnumber = cellid*msmsigsz(msmtype) + 5;
						obsvalue = msmgetDS(satdata[roughdopplerfieldnumber], sigdata[finedopplerfieldnumber], gobs->frequency(t_gobs(signalcode).band()), true);
						gobs->addobs(signalcode, obsvalue);// add this observation, HiPrecision

						//__ALL__OBSERVABLES__FOR__SIGNAL__ADDED___________________________________________________________________

						cellid++;//move on to the next field for reading a signal
					}
					cellmaskslot++;//move on to resolving the next signal
				}
			}
			obsid++;//move to record of the next satellite
		}
	}

	//cleanup
	delete[] satdata;
	delete[] sigdata;

	return 1;
}


// GPS MSM
int t_rtcm3::msg1071() {//GPS MSM1 - Compact_ pseudoranges
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM1((uint64_t*)headdata, numberofsatellites, 'G');
	return passMSMd();
}
int t_rtcm3::msg1072() {//GPS MSM2 - Compact_ phase-ranges
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM2((uint64_t*)headdata, numberofsatellites, 'G');
	return passMSMd();
}
int t_rtcm3::msg1073() {//GPS MSM3 - Compact_ pseudoranges & phase-ranges
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM3((uint64_t*)headdata, numberofsatellites, 'G');
	return passMSMd();
}
int t_rtcm3::msg1074() {//GPS MSM4 - Full_ pseudoranges & phase-ranges + CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM4((uint64_t*)headdata, numberofsatellites, 'G');
	return passMSMd();
}
int t_rtcm3::msg1075() {//GPS MSM5 - Full_ pseudoranges, phase-ranges, phase-range rate & CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM5((uint64_t*)headdata, numberofsatellites, 'G');
	return passMSMd();
}
int t_rtcm3::msg1076() {//GPS MSM6 - High Resolution_ pseudoranges & phase-ranges + CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM6((uint64_t*)headdata, numberofsatellites, 'G');
	return passMSMd();
}
int t_rtcm3::msg1077() {//GPS MSM7 - High Resolution_ pseudoranges, phase-ranges, phase-range rate & CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM7((uint64_t*)headdata, numberofsatellites, 'G');
	return passMSMd();
}

// GLONASS MSM
int t_rtcm3::msg1081() {//GLONASS MSM1 - Compact_ pseudoranges
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM1((uint64_t*)headdata, numberofsatellites, 'R');
	return passMSMd();
}
int t_rtcm3::msg1082() {//GLONASS MSM2 - Compact_ phase-ranges
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM2((uint64_t*)headdata, numberofsatellites, 'R');
	return passMSMd();
}
int t_rtcm3::msg1083() {//GLONASS MSM3 - Compact_ pseudoranges & phase-ranges
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM3((uint64_t*)headdata, numberofsatellites, 'R');
	return passMSMd();
}
int t_rtcm3::msg1084() {//GLONASS MSM4 - Full_ pseudoranges & phase-ranges + CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM4((uint64_t*)headdata, numberofsatellites, 'R');
	return passMSMd();
}
int t_rtcm3::msg1085() {//GLONASS MSM5 - Full_ pseudoranges, phase-ranges, phase-range rate & CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM5((uint64_t*)headdata, numberofsatellites, 'R');
	return passMSMd();
}
int t_rtcm3::msg1086() {//GLONASS MSM6 - High Resolution_ pseudoranges & phase-ranges + CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM6((uint64_t*)headdata, numberofsatellites, 'R');
	return passMSMd();
}
int t_rtcm3::msg1087() {//GLONASS MSM7 - High Resolution_ pseudoranges, phase-ranges, phase-range rate & CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM7((uint64_t*)headdata, numberofsatellites, 'R');
	return passMSMd();
}

// Galileo MSM
int t_rtcm3::msg1091() {//Galileo MSM1 - Compact_ pseudoranges
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM1((uint64_t*)headdata, numberofsatellites, 'E');
	return passMSMd();
}
int t_rtcm3::msg1092() {//Galileo MSM2 - Compact_ phase-ranges
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM2((uint64_t*)headdata, numberofsatellites, 'E');
	return passMSMd();
}
int t_rtcm3::msg1093() {//Galileo MSM3 - Compact_ pseudoranges & phase-ranges
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM3((uint64_t*)headdata, numberofsatellites, 'E');
	return passMSMd();
}
int t_rtcm3::msg1094() {//Galileo MSM4 - Full_ pseudoranges & phase-ranges + CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM4((uint64_t*)headdata, numberofsatellites, 'E');
	return passMSMd();
}
int t_rtcm3::msg1095() {//Galileo MSM5 - Full_ pseudoranges, phase-ranges, phase-range rate & CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM5((uint64_t*)headdata, numberofsatellites, 'E');
	return passMSMd();
}
int t_rtcm3::msg1096() {//Galileo MSM6 - High Resolution_ pseudoranges & phase-ranges + CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM6((uint64_t*)headdata, numberofsatellites, 'E');
	return passMSMd();
}
int t_rtcm3::msg1097() {//Galileo MSM7 - High Resolution_ pseudoranges, phase-ranges, phase-range rate & CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM7((uint64_t*)headdata, numberofsatellites, 'E');
	return passMSMd();
}

// QZSS MSM
int t_rtcm3::msg1111() {//QZSS MSM1 - Compact_ pseudoranges
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM1((uint64_t*)headdata, numberofsatellites, 'J');
	return passMSMd();
}
int t_rtcm3::msg1112() {//QZSS MSM2 - Compact_ phase-ranges
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM2((uint64_t*)headdata, numberofsatellites, 'J');
	return passMSMd();
}
int t_rtcm3::msg1113() {//QZSS MSM3 - Compact_ pseudoranges & phase-ranges
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM3((uint64_t*)headdata, numberofsatellites, 'J');
	return passMSMd();
}
int t_rtcm3::msg1114() {//QZSS MSM4 - Full_ pseudoranges & phase-ranges + CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM4((uint64_t*)headdata, numberofsatellites, 'J');
	return passMSMd();
}
int t_rtcm3::msg1115() {//QZSS MSM5 - Full_ pseudoranges, phase-ranges, phase-range rate & CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM5((uint64_t*)headdata, numberofsatellites, 'J');
	return passMSMd();
}
int t_rtcm3::msg1116() {//QZSS MSM6 - High Resolution_ pseudoranges & phase-ranges + CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM6((uint64_t*)headdata, numberofsatellites, 'J');
	return passMSMd();
}
int t_rtcm3::msg1117() {//QZSS MSM7 - High Resolution_ pseudoranges, phase-ranges, phase-range rate & CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM7((uint64_t*)headdata, numberofsatellites, 'J');
	return passMSMd();
}

// BeiDou MSM
int t_rtcm3::msg1121() {//BeiDou MSM1 - Compact_ pseudoranges
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM1((uint64_t*)headdata, numberofsatellites, 'C');
	return passMSMd();
}
int t_rtcm3::msg1122() {//BeiDou MSM2 - Compact_ phase-ranges
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM2((uint64_t*)headdata, numberofsatellites, 'C');
	return passMSMd();
}
int t_rtcm3::msg1123() {//BeiDou MSM3 - Compact_ pseudoranges & phase-ranges
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM3((uint64_t*)headdata, numberofsatellites, 'C');
	return passMSMd();
}
int t_rtcm3::msg1124() {//BeiDou MSM4 - Full_ pseudoranges & phase-ranges + CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM4((uint64_t*)headdata, numberofsatellites, 'C');
	return passMSMd();
}
int t_rtcm3::msg1125() {//BeiDou MSM5 - Full_ pseudoranges, phase-ranges, phase-range rate & CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM5((uint64_t*)headdata, numberofsatellites, 'C');
	return passMSMd();
}
int t_rtcm3::msg1126() {//BeiDou MSM6 - High Resolution_ pseudoranges & phase-ranges + CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM6((uint64_t*)headdata, numberofsatellites, 'C');
	return passMSMd();
}
int t_rtcm3::msg1127() {//BeiDou MSM7 - High Resolution_ pseudoranges, phase-ranges, phase-range rate & CNR
	uint64_t headdata[15] = { 0 };
	int numberofsatellites = 0;
	msgMSM7((uint64_t*)headdata, numberofsatellites, 'C');
	return passMSMd();
}


