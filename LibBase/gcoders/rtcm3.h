/*
 * rtcm3.h
 *
 *  Created on: Jan 6, 2015
 *      Author: jakub
 */

#ifndef RTCM3_H_
#define RTCM3_H_

#include <vector>
#include <memory>
#include "gcoders/rtcm_base.h"
#include "gdata/gobsgnss.h"
#include "gdata/gnav.h"

namespace gnut {

	class t_rtcm3 : public t_rtcm_base {
	public:
		t_rtcm3();
		virtual ~t_rtcm3();

		virtual int messagefound(int bitnumber, int & messagetype);//returns: message-over-size..{-2}; N..{-1}; Need-more-Data..{0}; Y..{message-length [bits]}
		virtual int readmessage();//RTCM3 starts with the byte (no need for bitnumber), RTCM2 keeps words in memory (no need for bitnumber)
		virtual void clearstate();

	private:
		void initobjstate();

		static const int64_t crc24q[256];
		static int64_t getCrc24Q(const unsigned char *buffer, int length);

		//message decoding variables
		int(t_rtcm3::*_messagetemplates[230 + 1])(); //message decoding functions field, RTCM3 ~ msg types 1001..1230 unsupported messages routed to function[0]; filled in constructor

		// Loss of Lock Indicator computation variables for MSM messages
		struct obslockdata {
			uint32_t epochtime; // epoch time [ms]
			uint32_t locktime;  // minimum lock time [ms]
			uint32_t suppcoeff; // supplementary coefficient k (RTCM 3.2 Amends-1-2 : Table 3.5-75); 0 if not high precision
		};
		typedef std::map< GOBS, obslockdata > t_mapllisig; // LLI detection algorithm data for signals of satellite;                  signal-type : {epoch-milliseconds, minimum-lock-time-milliseconds}
		typedef std::map< std::string, t_mapllisig >  t_mapllisv; // LLI detection algorithm data for satellites by RINEX-like ID string;   satellite-ID : map-of-signals'-data
		typedef std::map< uint32_t, t_mapllisv >  t_mapllist; // LLI detection algorithm data (map of stations);                          station-ID : map-of-satellites'-data
		t_mapllist   _stmap;   //==> t_mapllisv _svmap; ==> t_mapllisig _sigmap;


		// MSM SIGNAL MASKS:     |-01-|-02-|-03-|-04-|-05-|-06-|-07-|-08-|-09-|-10-|-11-|-12-|-13-|-14-|-15-|-16-|-17-|-18-|-19-|-20-|-21-|-22-|-23-|-24-|-25-|-26-|-27-|-28-|-29-|-30-|-31-|-32-|
		//GPS MSM observation codes
		const GOBS gpsCtable[32]{ X, C1C, C1P, C1W,   X,   X,   X, C2C, C2P, C2W,   X,   X,   X,   X, C2S, C2L, C2X,   X,   X,   X,   X, C5I, C5Q, C5X,   X,   X,   X,   X,   X,   X,   X,   X };
		const GOBS gpsLtable[32]{ X, L1C, L1P, L1W,   X,   X,   X, L2C, L2P, L2W,   X,   X,   X,   X, L2S, L2L, L2C,   X,   X,   X,   X, L5I, L5Q, L5X,   X,   X,   X,   X,   X,   X,   X,   X };
		const GOBS gpsDtable[32]{ X, D1C, D1P, D1W,   X,   X,   X, D2C, D2P, D2W,   X,   X,   X,   X, D2S, D2L, D2X,   X,   X,   X,   X, D5I, D5Q, D5X,   X,   X,   X,   X,   X,   X,   X,   X };
		//const GOBS gpsStable[32] {   X, S1C, S1P, S1W,   X,   X,   X, S2C, S2P, S2W,   X,   X,   X,   X, S2S, S2L, S2X,   X,   X,   X,   X, S5I, S5Q, S5X,   X,   X,   X,   X,   X,   X,   X,   X};
		  //GLO MSM observation codes
		const GOBS gloCtable[32]{ X, C1C, C1P,   X,   X,   X,   X, C2C, C2P,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X };
		const GOBS gloLtable[32]{ X, L1C, L1P,   X,   X,   X,   X, L2C, L2P,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X };
		const GOBS gloDtable[32]{ X, D1C, D1P,   X,   X,   X,   X, D2C, D2P,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X };
		//const GOBS gloStable[32] {   X, S1C, S1P,   X,   X,   X,   X, S2C, S2P,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X};
		  //GAL MSM observation codes
		const GOBS galCtable[32]{ X, C1C, C1A, C1B, C1X, C1Z,   X, C6C, C6A, C6B, C6X, C6Z,   X, C7I, C7Q, C7X,   X, C8I, C8Q, C8X,   X, C5I, C5Q, C5X,   X,   X,   X,   X,   X,   X,   X,   X };
		const GOBS galLtable[32]{ X, L1C, L1A, L1B, L1X, L1Z,   X, L6C, L6A, L6B, L6X, L6Z,   X, L7I, L7Q, L7X,   X, L8I, L8Q, L8X,   X, L5I, L5Q, L5X,   X,   X,   X,   X,   X,   X,   X,   X };
		const GOBS galDtable[32]{ X, D1C, D1A, D1B, D1X, D1Z,   X, D6C, D6A, D6B, D6X, D6Z,   X, D7I, D7Q, D7X,   X, D8I, D8Q, D8X,   X, D5I, D5Q, D5X,   X,   X,   X,   X,   X,   X,   X,   X };
		//const GOBS galStable[32] {   X, S1C, S1A, S1B, S1X, S1Z,   X, S6C, S6A, S6B, S6X, S6Z,   X, S7I, S7Q, S7X,   X, S8I, S8Q, S8X,   X, S5I, S5Q, S5X,   X,   X,   X,   X,   X,   X,   X,   X};
		  //QZS MSM observation codes
		const GOBS qzsCtable[32]{ X, C1C,   X,   X,   X,   X,   X,   X, C6S, C6L, C6X,   X,   X,   X, C2S, C2L, C2X,   X,   X,   X,   X, C5I, C5Q, C5X,   X,   X,   X,   X,   X, C1S, C1L, C1X };
		const GOBS qzsLtable[32]{ X, L1C,   X,   X,   X,   X,   X,   X, L6S, L6L, L6X,   X,   X,   X, L2S, L2L, L2X,   X,   X,   X,   X, L5I, L5Q, L5X,   X,   X,   X,   X,   X, L1S, L1L, L1X };
		const GOBS qzsDtable[32]{ X, D1C,   X,   X,   X,   X,   X,   X, D6S, D6L, D6X,   X,   X,   X, D2S, D2L, D2X,   X,   X,   X,   X, D5I, D5Q, D5X,   X,   X,   X,   X,   X, D1S, D1L, D1X };
		//const GOBS qzsStable[32] {   X, S1C,   X,   X,   X,   X,   X,   X, S6S, S6L, S6X,   X,   X,   X, S2S, S2L, S2X,   X,   X,   X,   X, S5I, S5Q, S5X,   X,   X,   X,   X,   X, S1S, S1L, S1X};
		  //BDS MSM observation codes
		const GOBS bdsCtable[32]{ X, C1I, C1Q, C1X,   X,   X,   X, C6I, C6Q, C6X,   X,   X,   X, C7I, C7Q, C7X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X };
		const GOBS bdsLtable[32]{ X, L1I, L1Q, L1X,   X,   X,   X, L6I, L6Q, L6X,   X,   X,   X, L7I, L7Q, L7X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X };
		const GOBS bdsDtable[32]{ X, D1I, D1Q, D1X,   X,   X,   X, D6I, D6Q, D6X,   X,   X,   X, D7I, D7Q, D7X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X };
		//const GOBS bdsStable[32] {   X, S1I, S1Q, S1X,   X,   X,   X, S6I, S6Q, S6X,   X,   X,   X, S7I, S7Q, S7X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X};

		  //Epoch grouping variables
			// "=>> Observation
		  //TODO: function to produce warnings if time inconsistencies are found within observation vector??
		  //double _msmmaxdt... //maximum allowed deviation caused by type conversions etc. in a single epoch of MSM observations [s]
		  //double _legmaxdt... //maximum allowed deviation caused by type conversions etc. in a single epoch of MSM observations [s]
		std::vector<shared_ptr<t_gobsgnss>> _msmepochdata; //MSM observations data vector for a single epoch (sent to t_gall when end of epoch flag is caught)

		t_gtime* _legobsepoch;
		std::vector<shared_ptr<t_gobsgnss>> _legepochdata; //Legacy observations data vector for a single epoch (sent to t_gall when Multiple Message Indicator is not set)
		//NOTE: the values of observation epoch data should not be combined, when being added to gall
			// "=>> Ephemeris
		std::map<char, shared_ptr<t_gnav>>  _lastephdata; //<GNSS system identifier character, last ephemeris data for system>; do not send data to t_gall if time matches
		//NOTE: there is no way to handle different sets of ephemeris for Galileo in this implementation

		std::map<int32_t, int32_t> _stationcnt; //station messages counter, <rtcm-stationID, number of messages decoded> (raise warning when a station is being added and map is not empty)

		string resolvestation(uint64_t stationid); //convert station ID into RINEX-like string, generate warning if new station (except the first one) starts appearing in the stream

		//message decoding functions
		int msguspt();//unsupported message type, (default)

		// ################ RTK Observation Data #################
		shared_ptr<t_gobsgnss> leggetgobsgnss(t_gtime epoch, uint64_t siteid, string sat);//get pointer to proper t_gobsgnss object in a _legepochdata (create it if necessary)
		int passRTKo(); // pass legacy observation data to t_gall if epoch end

		// ############ GPS RTK Observation Data #################
		int gpsRTKh(uint64_t* headdata);//GPS RTK Observations header
		int msg1001();//GPS L1 Observables
		int msg1002();//GPS L1 Extended Observables
		int msg1003();//GPS L1 + L2 Observables
		int msg1004();//GPS L1 + L2 Extended Observables

		//...

		// ########## GLONASS RTK Observation Data ###############
		int gloRTKh(uint64_t* headdata);//GLONASS RTK Observations header
		int msg1009();//GLONASS L1 Observables
		int msg1010();//GLONASS L1 Extended Observables
		int msg1011();//GLONASS L1 + L2 Observables
		int msg1012();//GLONASS L1 + L2 Extended Observables

		// ...
		// ################ Ephemeris Data #######################
		int passEPHd(char gnsssys, string prnstr, t_gnav* newephdata); // function for releasing eph_data to t_gall structure, if newephdata has identical t_oe/t_oc as oldeph, data is dropped
		// expects a pointer to object constructed using 'new' as newephdata

		int msg1019();//GPS ephemeris
		int msg1020();//GLONASS ephemeris
		int msg1044();//QZSS ephemeris
		int msg1045();//Galileo F/NAV ephemeris

		// ################ MSM DATA DECODING ####################`
		  //general MSM decoding functions

		  //field-reading function for InternalLooping
		int msmgetiloopfield(int32_t* resultdatafield, int& srcstartbit, int srcfieldnumber, int loopsize, int srcfieldsize, bool issigned = false);
		//G-nut structure integration functions
		bool msmgetsatid(char* satidcstr, int msmsatmaskid, char gnsssys);
		shared_ptr<t_gobsgnss> msmgetgobsgnss(uint64_t* msmheaddata, int satmaskindex, char gnsssys, const GOBS* sigmapC, const GOBS* sigmapL, const GOBS* sigmapD/*, const GOBS* sigmapS*/);


		//MSM data field processing into observables
		double msmgetPS(int64_t n_ms, int32_t roughdata, int32_t finedata, bool highprecision = false);//PSeudo-range reconstruction
		double msmgetPH(int64_t n_ms, int32_t roughdata, int32_t finedata, double frequency, bool highprecision = false);//PHase-cycle / PHase-range reconstruction
		double msmgetDS(int32_t roughdata, int32_t finedata, double frequency, bool highprecision = false);//Doppler-Shift / phase-range-rate reconstruction
		//double msmgetSS();//Signal-Strength reconstruction

			//Loss of Lock Indicator computation functions
		obslockdata msmgetMLT(uint32_t locktimeindicator, bool highprecision = false);//Minimum-Lock-Time reconstruction
		bool msmlligetprevtimes(obslockdata& oldlockdata, obslockdata newlockdata, uint32_t stationID, string satelliteID, GOBS signaltype);//search module's internal data for record of last phase observation, returns false if references weren't filled (data not found)
		int msmgetLLI(uint32_t stationID, string satelliteID, GOBS signaltype, uint32_t epochtime, obslockdata lockdata);//Loss-of-Lock-Indicator computation

		//default MSM message handling function (NOTE: catches the end of epoch when it ends with unimplemented MSM message)
		int passMSMd();  //returns number of records passed to G-NUT, unsupported function handler calls this for all reserved MSM1-7s
						//NOTE: all supported MSMs MUST call this function when returning

			// MSM sizes of satellite and signal data sections
		int msmsatsz(int msmtype); //get number of fields describing a satellite for this type of MSM
		int msmsigsz(int msmtype); //get number of fields describing a signal for this type of MSM

			// MSM section reading functions
		int msgMSMh(uint64_t* datafield, int& nsat, int& nsig);//general MSM header*
		int msgMSMsat123(int32_t* datafield, int& nsat, int startbit);//satellite MSM data for MSM1, MSM2 & MSM3 *
		int msgMSMsat46(int32_t* datafield, int& nsat, int startbit);//satellite MSM data for MSM4 & MSM6 *
		int msgMSMsat57(int32_t* datafield, int& nsat, int startbit);//satellite MSM data for MSM5 & MSM7 *
		int msgMSMsig1(int32_t* datafield, int& ncell, int startbit);//signal MSM data for MSM1 *
		int msgMSMsig2(int32_t* datafield, int& ncell, int startbit);//signal MSM data for MSM2 *
		int msgMSMsig3(int32_t* datafield, int& ncell, int startbit);//signal MSM data for MSM3 *
		int msgMSMsig4(int32_t* datafield, int& ncell, int startbit);//signal MSM data for MSM4 *
		int msgMSMsig5(int32_t* datafield, int& ncell, int startbit);//signal MSM data for MSM5 *
		int msgMSMsig6(int32_t* datafield, int& ncell, int startbit);//signal MSM data for MSM6 *
		int msgMSMsig7(int32_t* datafield, int& ncell, int startbit);//signal MSM data for MSM7 *
			// * - the datafield should be filled with zeroes and be sufficiently large for maximum possible data size (13 for header, ??? for sat, ??? for sig)
			//???- the size of datafields required for sig and sat data extraction varies and therefore will be dynamically allocated within parsers

		int msgMSM1(uint64_t* headdata, int& nsat, char gnsssys);//MSM1 field parser
		int msgMSM2(uint64_t* headdata, int& nsat, char gnsssys);//MSM2 field parser
		int msgMSM3(uint64_t* headdata, int& nsat, char gnsssys);//MSM3 field parser
		int msgMSM4(uint64_t* headdata, int& nsat, char gnsssys);//MSM4 field parser
		int msgMSM5(uint64_t* headdata, int& nsat, char gnsssys);//MSM5 field parser
		int msgMSM6(uint64_t* headdata, int& nsat, char gnsssys);//MSM6 field parser
		int msgMSM7(uint64_t* headdata, int& nsat, char gnsssys);//MSM7 field parser

		  //specific MSM decoding functions
		int msg1071();//GPS MSM1 - Compact_ pseudoranges
		int msg1072();//GPS MSM2 - Compact_ phase-ranges
		int msg1073();//GPS MSM3 - Compact_ pseudoranges & phase-ranges
		int msg1074();//GPS MSM4 - Full_ pseudoranges & phase-ranges + CNR
		int msg1075();//GPS MSM5 - Full_ pseudoranges, phase-ranges, phase-range rate & CNR
		int msg1076();//GPS MSM6 - High Resolution_ pseudoranges & phase-ranges + CNR
		int msg1077();//GPS MSM7 - High Resolution_ pseudoranges, phase-ranges, phase-range rate & CNR

		int msg1081();//GLONASS MSM1 - Compact_ pseudoranges
		int msg1082();//GLONASS MSM2 - Compact_ phase-ranges
		int msg1083();//GLONASS MSM3 - Compact_ pseudoranges & phase-ranges
		int msg1084();//GLONASS MSM4 - Full_ pseudoranges & phase-ranges + CNR
		int msg1085();//GLONASS MSM5 - Full_ pseudoranges, phase-ranges, phase-range rate & CNR
		int msg1086();//GLONASS MSM6 - High Resolution_ pseudoranges & phase-ranges + CNR
		int msg1087();//GLONASS MSM7 - High Resolution_ pseudoranges, phase-ranges, phase-range rate & CNR

		int msg1091();//Galileo MSM1 - Compact_ pseudoranges
		int msg1092();//Galileo MSM2 - Compact_ phase-ranges
		int msg1093();//Galileo MSM3 - Compact_ pseudoranges & phase-ranges
		int msg1094();//Galileo MSM4 - Full_ pseudoranges & phase-ranges + CNR
		int msg1095();//Galileo MSM5 - Full_ pseudoranges, phase-ranges, phase-range rate & CNR
		int msg1096();//Galileo MSM6 - High Resolution_ pseudoranges & phase-ranges + CNR
		int msg1097();//Galileo MSM7 - High Resolution_ pseudoranges, phase-ranges, phase-range rate & CNR

		int msg1111();//QZSS MSM1 - Compact_ pseudoranges
		int msg1112();//QZSS MSM2 - Compact_ phase-ranges
		int msg1113();//QZSS MSM3 - Compact_ pseudoranges & phase-ranges
		int msg1114();//QZSS MSM4 - Full_ pseudoranges & phase-ranges + CNR
		int msg1115();//QZSS MSM5 - Full_ pseudoranges, phase-ranges, phase-range rate & CNR
		int msg1116();//QZSS MSM6 - High Resolution_ pseudoranges & phase-ranges + CNR
		int msg1117();//QZSS MSM7 - High Resolution_ pseudoranges, phase-ranges, phase-range rate & CNR

		int msg1121();//BeiDou MSM1 - Compact_ pseudoranges
		int msg1122();//BeiDou MSM2 - Compact_ phase-ranges
		int msg1123();//BeiDou MSM3 - Compact_ pseudoranges & phase-ranges
		int msg1124();//BeiDou MSM4 - Full_ pseudoranges & phase-ranges + CNR
		int msg1125();//BeiDou MSM5 - Full_ pseudoranges, phase-ranges, phase-range rate & CNR
		int msg1126();//BeiDou MSM6 - High Resolution_ pseudoranges & phase-ranges + CNR
		int msg1127();//BeiDou MSM7 - High Resolution_ pseudoranges, phase-ranges, phase-range rate & CNR
		// ...

	};

}
#endif /* RTCM3_H_ */
