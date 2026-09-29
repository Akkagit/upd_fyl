/*
 * rtcm_base.h
 *
 *  Created on: Jan 5, 2015
 *      Author: Jakub
 */

#ifndef RTCM_BASE_H_
#define RTCM_BASE_H_

#include <map>
#include <vector>
#include <string>
#include <cstring>
#include "gdata/gdata.h"
#include "gutils/gmutex.h"
#include "gutils/gtime.h"
#include "gset/gsetbase.h"  // TEMPORARY ?

static const char BITSELECTOR[8] = { -0x80,0x40,0x20,0x10,0x08,0x04,0x02,0x01 };//MSB ~ #0, LSB ~ #7
namespace gnut {
	class t_rtcm_base {
	public:
		t_rtcm_base();
		virtual ~t_rtcm_base();

		void rtcminit(char** targetbuffer, int* buffsize, map<string, t_gdata*>* datapointer, t_glog** logpointer);//fill pointers to buffer and other variables appropriately

		//message detection
		virtual void clearstate() {}//clear internal state of this utility class ( should be done before every call of messagefound()... )
		virtual int messagefound(int bitnumber, int& messagetype) { return -1; }//returns: message-over-size..{-2}; N..{-1}; Need-more-Data..{0}; Y..{message-length}
		virtual int readmessage() { return 0; }//RTCM3 starts with the byte (no need for bitnumber), RTCM2 keeps words in memory (no need for bitnumber)
		//read message is to be called right after successful message detection (before consuming or otherwise changing buffer contents or internal rtcm_object state)

		// time estimation arguments
		void setmaxdteph(double maxdteph);
		double getmaxdteph();
		void setmaxdtobs(double maxdtobs);
		double getmaxdtobs();

		void settref(t_gtime tref);
		void settrefrealtime();
		const t_gtime* gettref() const;

	protected:

		char** _buffer;								//pointer to buffer pointer/array
		int* _buffsize;								// pointer to buffer size
		map<string, t_gdata*>*   _data;                // pointer to data pointer
		t_glog**                _log;                 // pointer to log pointer

		//double bmag(double number, int exp);//binary magnitude bmag(number, exp) = number * 2^exp

		//binary buffer to C++style type conversion functions
		virtual bool getonebit(int bitnumber);//return pointer to array containing number of valid bits from specified place in thr buffer
		virtual int64_t get2int(int startbit, int numberofbits);//return value of 2's complement integer of specific length, read from buffer at specified location
		virtual uint64_t getuint(int startbit, int numberofbits);//return value of unsigned integer of specific length, read from buffer at specified location
		virtual int64_t getsint(int startbit, int numberofbits);//return value of sign-magnitude integer of specific length, read from buffer at specified location
		virtual int64_t decode2int(uint64_t unumber, int numberofbits);//take a given number of LSB bits of a unsigned long and return number they represent if perceived as 2's complement

		//utility
		double bmag(double number, int exp);//binary magnitude bmag(number, exp) = number * 2^exp

		//GNSS system time estimation
		// >> general <<
		// convert relative time within an epoch (tepoch) to an absolute time record (all doubles are in seconds)
		// maxdelta is the maximum allowed deviation of the estimated value from the reference time, produces a warning if exceeded, ! set to 0 if >= (epochsz/2)
		t_gtime estimateabstime(t_gtime tref, double tepoch, double maxdelta, double epochsz, t_gtime tstart, t_gtime::t_tsys ts);

		//Following functions estimate absolute GNSSs system times
		// tepoch   : time of epoch [sec] (to be converted to absolute value)
		// tref     : reference time, current system time is used if not explicitly specified
		// maxdelta : maximum deviation of estimated value from reference value, if exceeded a warning is generated in log
		//            maxdelta should not exceed length of half of the estimated epoch; if it does, it is modified to match it
		// weekepoch: a switch between false->the long epoch(GNSS system specific) and true->the week epoch (for observations)
		// >> GPS <<
		t_gtime estimateepochGPS(double tepoch, t_gtime tref, double maxdelta, bool weekepoch = false);
		t_gtime estimateepochGPS(double tepoch, t_gtime tref, bool weekepoch = false);
		t_gtime estimateepochGPS(double tepoch, double maxdelta, bool weekepoch = false);
		t_gtime estimateepochGPS(double tepoch, bool weekepoch = false);
		t_gtime estimateepochGPShour(double tepoch);
		// >> GLO <<
		t_gtime estimateepochGLO(double tepoch, t_gtime tref, double maxdelta, bool weekepoch = false);
		t_gtime estimateepochGLO(double tepoch, t_gtime tref, bool weekepoch = false);
		t_gtime estimateepochGLO(double tepoch, double maxdelta, bool weekepoch = false);
		t_gtime estimateepochGLO(double tepoch, bool weekepoch = false);
		t_gtime estimateepochGLOday(double tepoch);
		t_gtime estimateepochGLOhour(double tepoch);
		// >> GAL <<
		t_gtime estimateepochGAL(double tepoch, t_gtime tref, double maxdelta, bool weekepoch = false);
		t_gtime estimateepochGAL(double tepoch, t_gtime tref, bool weekepoch = false);
		t_gtime estimateepochGAL(double tepoch, double maxdelta, bool weekepoch = false);
		t_gtime estimateepochGAL(double tepoch, bool weekepoch = false);
		// >> QZS <<
		t_gtime estimateepochQZS(double tepoch, t_gtime tref, double maxdelta, bool weekepoch = false);
		t_gtime estimateepochQZS(double tepoch, t_gtime tref, bool weekepoch = false);
		t_gtime estimateepochQZS(double tepoch, double maxdelta, bool weekepoch = false);
		t_gtime estimateepochQZS(double tepoch, bool weekepoch = false);
		// >> BDS <<
		t_gtime estimateepochBDS(double tepoch, t_gtime tref, double maxdelta, bool weekepoch = false);
		t_gtime estimateepochBDS(double tepoch, t_gtime tref, bool weekepoch = false);
		t_gtime estimateepochBDS(double tepoch, double maxdelta, bool weekepoch = false);
		t_gtime estimateepochBDS(double tepoch, bool weekepoch = false);

		t_gtime _tref;
		t_gtime* _trefptr;
		double _maxdteph; //holds maximum deviation of estimated times of ephemeris messages, from given relative times
		double _maxdtobs;
	};
}
#endif /* RTCM_BASE_H_ */
