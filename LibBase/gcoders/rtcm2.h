/*
 * rtcm2.h
 *
 *  Created on: Jan 6, 2015
 *      Author: jakub
 */

#ifndef RTCM2_H_
#define RTCM2_H_

#include <vector>
#include <memory>

#include "gcoders/rtcm_base.h"
#include "gdata/gobsgnss.h"
#include "gdata/gnav.h"

namespace gnut {

	static const unsigned char ptable[] = {
	  0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,
	  1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,
	  1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,
	  0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,
	  1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,
	  0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,
	  0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,
	  1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0
	};

	static const unsigned long parityform[6] = {
			0xBB1F3480l,
			0x5D8F9A40l,
			0xAEC7CD00l,
			0x5763E680l,
			0x6BB1F340l,
			0x8B7A89C0l };

	class t_rtcm2 : public t_rtcm_base {
	public:
		t_rtcm2();
		virtual ~t_rtcm2();

		virtual int messagefound(int bitnumber, int & messagetype);//returns: message-over-size..{-2}; N..{-1}; Need-more-Data..{0}; Y..{message-length [bits]}
		virtual void clearstate();
		virtual int readmessage();//RTCM3 starts with the byte (no need for bitnumber), RTCM2 keeps words in memory (no need for bitnumber)

	private:

		std::vector<shared_ptr<t_gobsgnss>> _obsepochdata; //Observations data vector for a single epoch (sent to t_gall when Multiple Message Indicator is not set)
		std::map<int32_t, int32_t> _stationcnt; //station messages counter, <rtcm-stationID, number of messages decoded> (raise warning when a station is being added and map is not empty)

		uint64_t getnumber(int messageindex, int startbitnumber, int length);
		bool check_word_parity(uint64_t & aword);
		bool resync_word_parity(uint64_t & aword);
		int trymessagefound(int bitnumber, int & messagetype);
		int getwordbits(int bitnumber, uint64_t & worddata);
		void printbin(uint64_t data);

		unsigned char count_parity(uint64_t aword);

		unsigned char _lastparity;
		bool _synced;
		uint64_t _words[33];

		string resolvestation(uint64_t stationid); //convert station ID into RINEX-like string, generate warning if new station (except the first one) starts appearing in the stream

		shared_ptr<t_gobsgnss> getgobsgnss(t_gtime epoch, uint64_t siteid, string sat);//get pointer to proper t_gobsgnss object in a _obsepochdata (create it if necessary)
		int passOBSd();

		//specific message decoding
		int(t_rtcm2::*_messagetemplates[64])(); //message decoding functions field, RTCM2 ~ msg types 0..63 (type 0 => type 64)
		//specific message decoding functions
		int msguspt();//unsupported message type, (default)
		//...
		int msg17();//GPS ephemeris
		int msg18();//RTK Uncorrected Carrier Phase Measurements
		int msg19();//RTK Uncorrected Pseudorange Measurements
		//...



	};
}
#endif /* RTCM2_H_ */
