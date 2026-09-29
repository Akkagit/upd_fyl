
#ifndef GNSS_H
#define GNSS_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: definition of GNSS data
  Version: $ Rev: $

  2012-09-26 /JD: created
  2019-11-19 /glfeng: add BDS-3 frequency information
		     // BeiDou C02(IQX) C07(7-IQX) C06(IQX) C05(5-DPX)  C09(7-DPZ)  C08(8-DPX)  C01(1-DPX)
		     //        B1I      B2I(BDS-2) B3I      B2a(BDS-3)  B2b(BDS-3)  B2(BDS-3)   B1C(BDS-3)
             

-*/

#include <map>
#include <set>
#include <vector>
#include <string>

#include "gutils/gtriple.h"
#include "gutils/gcommon.h"   // pragma
#include "gutils/ExportBase.h"
using namespace std;

// ------------------------------------------------------------------------------------------------------
// ENUMS
// ------------------------------------------------------------------------------------------------------

namespace gnut
{

	// GNSS systems and augmentations
	// ------------------------------
	enum GSYS 
	{ // GXX = -1,
		GPS, GAL, GLO, BDS, QZS, SBS, IRN, LEO, GNS
	};

	// GNSS freq Sequence ID
	// ---------------------
	enum FREQ_SEQ 
	{
		FREQ_1 = 1, FREQ_2 = 2, FREQ_3 = 3, FREQ_4 = 4, FREQ_5 = 5, FREQ_6 = 6, FREQ_7 = 7,
		FREQ_X = 999
	};

	// GNSS frequencies
	// ----------------
	enum GFRQ 
	{ // FXX = -1,
		G01 = 10, G02 = 11, G05 = 12,     // GPS
		R01 = 20, R02 = 21,               // GLONASS FDMA
		R01_CDMA = 30, R02_CDMA = 31,
		R03_CDMA = 32, R05_CDMA = 33,     // GLONASS CDMA
		E01 = 50, E05 = 51, E07 = 52,
		E08 = 53, E06 = 54,     // Galileo
		// BeiDou C02(IQX) C07(7-IQX) C06(IQX) C05(5-DPX)  C09(7-DPZ)  C08(8-DPX)  C01(1-DPX)
		//        B1I      B2I(BDS-2) B3I      B2a(BDS-3)  B2b(BDS-3)  B2(BDS-3)   B1C(BDS-3)
		C02 = 60, C07 = 61, C06 = 62, C05 = 63, C09 = 64, C08 = 65, C01 = 66,
		J01 = 70, J02 = 71, J05 = 72,
		J06 = 73,     // QZSS
		S01 = 80, S05 = 81,               // SBAS
		I05 = 90,  // I09,                // IRNSS
		LAST_GFRQ = 999
	};

	// GNSS receiver types
	// -------------------
	enum RECTYPE {
		P1P2,    // receiver providing C1, P1, P2
		C1X2,    // cross-correlation
		C1P2     // modern receivers providing C1, P2
	};

	// Broadcast messages types
	// ------------------------
	// Broadcast messages types
	enum GNAVTYPE { FNAV, INAV, INAV_E01, INAV_E07, CNAV, NAV };

	// enum GNAVTYPE {  NAV_G01, CNAV_G02, CNAV_G05,
	//                  NAV_R01,
	//                  NAV_C02,
	//                  FNAV_E05, INAV_E01, INAV_E07, CNAV_E06, GNAV_E01, GNAV_E06,
	//                  NAV_DEF
	//              };

	 // GNSS type/band/attr definitions
	 // -------------------------------
	enum GOBSTYPE {
		TYPE_C = 1, TYPE_L = 2, TYPE_D = 3, TYPE_S = 4,
		TYPE_P = 101, // only for P-code!
		TYPE_SLR,
		TYPE = 999  // ""  UNKNOWN
	};
	enum GOBSBAND {
		BAND_1 = 1, BAND_2 = 2, BAND_3 = 3, BAND_5 = 5,
		BAND_6 = 6, BAND_7 = 7, BAND_8 = 8, BAND_9 = 9, // Band_9 is set for BDS-3 B2b
		BAND_A = 101, BAND_B = 102, BAND_C = 103, BAND_D = 104,
		BAND_SLR,
		BAND = 999  // ""  UNKNOWN
	};
	enum GOBSATTR {
		ATTR_A, ATTR_B, ATTR_C, ATTR_D, ATTR_I, ATTR_L, ATTR_M, ATTR_N,
		ATTR_P, ATTR_Q, ATTR_S, ATTR_W, ATTR_X, ATTR_Y, ATTR_Z,
		ATTR_NULL,    // " " 2CHAR code
		ATTR = 999    // ""  UNKNOWN
	};

	// GNSS observations
	// -----------------
	enum GOBS {

		// psedorange [in meters] (RINEX 3.x)
		C1A = 0, C1B, C1C, C1D, C1I, C1L, C1M, C1P, C1S, C1Q, C1W, C1X, C1Y, C1Z,
		C2C, C2D, C2I, C2L, C2M, C2P, C2S, C2Q, C2W, C2X, C2Y,
		C3I, C3Q, C3X,
		C5A, C5B, C5C, C5D, C5I, C5P, C5Q, C5X,
		C6A, C6B, C6C, C6I, C6L, C6S, C6Q, C6X, C6Z,
		C7I, C7Q, C7X,
		C8D, C8I, C8P, C8Q, C8X,
		C9D, C9P, C9Z,  // BDS-3 B2b

		// carrier phase [in whole cycles] (RINEX 3.x)
		L1A = 100, L1B, L1C, L1D, L1I, L1L, L1M, L1N, L1P, L1S, L1Q, L1W, L1X, L1Y, L1Z,
		L2C, L2D, L2I, L2L, L2M, L2N, L2P, L2S, L2Q, L2W, L2X, L2Y,
		L3I, L3Q, L3X,
		L5A, L5B, L5C, L5D, L5I, L5P, L5Q, L5X,
		L6A, L6B, L6C, L6I, L6L, L6S, L6Q, L6X, L6Z,
		L7I, L7Q, L7X,
		L8D, L8I, L8P, L8Q, L8X,
		L9D, L9P, L9Z,  // BDS-3 B2b

		// doppler [cycles/sec] (RINEX 3.x)
		D1A = 200, D1B, D1C, D1D, D1I, D1L, D1M, D1N, D1P, D1S, D1Q, D1W, D1X, D1Y, D1Z,
		D2C, D2D, D2I, D2L, D2M, D2N, D2P, D2S, D2Q, D2W, D2X, D2Y,
		D3I, D3Q, D3X,
		D5A, D5B, D5C, D5D, D5I, D5P, D5Q, D5X,
		D6A, D6B, D6C, D6I, D6L, D6S, D6Q, D6X, D6Z,
		D7I, D7Q, D7X,
		D8D, D8I, D8P, D8Q, D8X,
		D9D, D9P, D9Z,  // BDS-3 B2b


		// signal strength [DBHZ] (RINEX 3.x)
		S1A = 300, S1B, S1C, S1D, S1I, S1L, S1M, S1N, S1P, S1S, S1Q, S1W, S1X, S1Y, S1Z,
		S2C, S2D, S2I, S2L, S2M, S2N, S2P, S2S, S2Q, S2W, S2X, S2Y,
		S3I, S3Q, S3X,
		S5A, S5B, S5C, S5D, S5I, S5P, S5Q, S5X,
		S6A, S6B, S6C, S6I, S6L, S6S, S6Q, S6X, S6Z,
		S7I, S7Q, S7X,
		S8D, S8I, S8P, S8Q, S8X,
		S9D, S9P, S9Z,  // BDS-3 B2b

		// special cases: v2.x or unknown tracking modes
		P1 = 1000, P2, P5, C1, C2, C5, C6, C7, C8, CA, CB, CC, CD,
		L1 = 1100, L2, L5, L6, L7, L8, LA, LB, LC, LD,
		D1 = 1200, D2, D5, D6, D7, D8, DA, DB, DC, DD,
		S1 = 1300, S2, S5, S6, S7, S8, SA, SB, SC, SD,

		X // LAST_GOBS
	};

	enum GOBS_LC { LC_UNDEF = 0, LC_L1 = 1, LC_L2 = 2, LC_L3 = 3, LC_L4 = 4, LC_L5 = 5, LC_IF, LC_MW, LC_NL, LC_WL, LC_GF };


	// ------------------------------------------------------------------------------------------------------
	// TYPEDEF
	// ------------------------------------------------------------------------------------------------------
	typedef vector< GOBSATTR >          t_vec_attr;
	typedef vector< GOBSBAND >          t_vec_band;
	typedef vector< GFRQ     >          t_vec_freq;
	typedef map< GOBSTYPE, t_vec_attr > t_map_attr;
	typedef map< GOBSBAND, t_map_attr > t_map_type;
	typedef map< GSYS, set<string> >    t_map_sats;
	typedef map< GSYS, set<string> >    t_map_gnav;
	typedef map< GSYS, t_map_type  >    t_map_gnss;
	typedef map< GSYS, t_vec_band  >    t_map_band;
	typedef map< GSYS, t_vec_freq  >    t_map_freq;

	typedef map< GOBSBAND, t_gtriple >  t_map_pcos;  // triple: ATX  NORTH / EAST / UP
	typedef map< GSYS, t_map_pcos >  t_map_offs;

	// ------------------------------------------------------------------------------------------------------
	// GLOBAL FUNCTIONS
	// ------------------------------------------------------------------------------------------------------
	BASE_LIBRARY_EXPORT GOBSATTR str2gobsattr(string s);           // get GOBSATTR enum from gobs string
	BASE_LIBRARY_EXPORT GOBSBAND str2gobsband(string s);           // get GOBSBAND enum from gobs string
	BASE_LIBRARY_EXPORT GOBSTYPE str2gobstype(string s);           // get GOBSTYPE enum from gobs string
	BASE_LIBRARY_EXPORT GNAVTYPE str2gnavtype(string s);           // get GNAVTYPE enum from gobs string
	BASE_LIBRARY_EXPORT FREQ_SEQ str2sysfreq(string s);           // get FREQ_SEQ enum from string (add by xiongyun)
	BASE_LIBRARY_EXPORT FREQ_SEQ str2gnssfreq(string s);           // get FREQ_SEQ enum from gobs string

   // GOBSATTR gobs2gobsattr( GOBS o );            // get GOBSATTR enum from gobs string
   // GOBSBAND gobs2gobsband( GOBS o );            // get GOBSBAND enum from gobs string
   // GOBSTYPE gobs2gobstype( GOBS o );            // get GOBSTYPE enum from gobs string

	BASE_LIBRARY_EXPORT GOBSATTR char2gobsattr(char c);            // get GOBSATTR enum from char
	BASE_LIBRARY_EXPORT GOBSBAND char2gobsband(char c);            // get GOBSBAND enum from char
	BASE_LIBRARY_EXPORT GOBSBAND  int2gobsband(int  c);            // get GOBSBAND enum from char   
	BASE_LIBRARY_EXPORT GOBSTYPE char2gobstype(char c);            // get GOBSTYPE enum from char
	BASE_LIBRARY_EXPORT FREQ_SEQ char2gnssfreq(char c);
	BASE_LIBRARY_EXPORT GOBS_LC   int2gobsfreq(int  c);

	BASE_LIBRARY_EXPORT string gobsattr2str(GOBSATTR e);           // get string enum from GOBSATTR 
	BASE_LIBRARY_EXPORT string gobsband2str(GOBSBAND e);           // get string enum from GOBSBAND
	BASE_LIBRARY_EXPORT string gobstype2str(GOBSTYPE e);           // get string enum from GOBSTYPE 
	BASE_LIBRARY_EXPORT string gnavtype2str(GNAVTYPE e);           // get string enum from GNAVTYPE

	BASE_LIBRARY_EXPORT string gobs2str(GOBS);                     // get string from GOBS enum
	BASE_LIBRARY_EXPORT GOBS str2gobs(string s);                   // get GOBS enum from string
	BASE_LIBRARY_EXPORT GOBS tba2gobs(GOBSTYPE t, GOBSBAND b, GOBSATTR a); // get GOBS from type, band, and attribute
	BASE_LIBRARY_EXPORT string gfreqseq2str(FREQ_SEQ f);           // convert FREQ_SEQ to string

	BASE_LIBRARY_EXPORT int gobs2band(GOBS o);                     // get band from GOBS enum

	BASE_LIBRARY_EXPORT GOBS pha2snr(GOBS o);                      // get GOBS enum (pha->snr)
	BASE_LIBRARY_EXPORT GOBS pl2snr(GOBS o);                       // get GOBS enum (pha or code->snr) add wh
	BASE_LIBRARY_EXPORT bool gobs_code(GOBS o);                    // get true for code obs
	BASE_LIBRARY_EXPORT bool gobs_phase(GOBS o);                   // get true for phase obs
	BASE_LIBRARY_EXPORT bool gobs_doppler(GOBS o);                   // get true for doppler obs
	BASE_LIBRARY_EXPORT bool gobs_snr(GOBS o);                    // get true for snr obs 

	BASE_LIBRARY_EXPORT t_map_sats GNSS_SATS();                      // static map of default GNSS satellites                       
	BASE_LIBRARY_EXPORT t_map_gnav GNSS_GNAV();                      // static map of default GNSS navigation types
	//BASE_LIBRARY_EXPORT t_map_freq GNSS_FREQ_PRIORITY();	            // static map of default GNSS freq priorities
	//BASE_LIBRARY_EXPORT t_map_band GNSS_BAND_PRIORITY(); 	          // static map of default GNSS band priorities
	BASE_LIBRARY_EXPORT t_map_gnss GNSS_DATA_PRIORITY(); 	          // static map of default GNSS data types/bands/attrs priorities

	t_map_band GNSS_BAND_SORTED(); 	            // static map of sorted GNSS band w.r.t. wavelength
	vector<GOBSBAND> sort_band(GSYS gs, set<GOBSBAND>& bands);        // sort set of bands w.r.t. wavelength

	t_map_offs GNSS_PCO_OFFSETS(); 	      // static map of default GNSS PCO offsets

	set<GSYS> GNSS_SUPPORTED();                  // supported GNSS


	// const t_map_sats GNSS_SATS = {
	//	{ GPS, { (string)"G01",(string)"G02",(string)"G03",(string)"G04",(string)"G05",(string)"G06",(string)"G07",(string)"G08",(string)"G09",(string)"G10",
	//		 (string)"G11",(string)"G12",(string)"G13",(string)"G14",(string)"G15",(string)"G16",(string)"G17",(string)"G18",(string)"G19",(string)"G20",
	//		 (string)"G21",(string)"G22",(string)"G23",(string)"G24",(string)"G25",(string)"G26",(string)"G27",(string)"G28",(string)"G29",(string)"G30",
	//		 (string)"G31",(string)"G32" } },

	//	{GLO, { (string)"R01",(string)"R02",(string)"R03",(string)"R04",(string)"R05",(string)"R06",(string)"R07",(string)"R08",(string)"R09",(string)"R10",
	//		 (string)"R11",(string)"R12",(string)"R13",(string)"R14",(string)"R15",(string)"R16",(string)"R17",(string)"R18",(string)"R19",(string)"R20",
	//		 (string)"R21",(string)"R22",(string)"R23",(string)"R24" }},

	//	{GAL, { (string)"E01",(string)"E02",(string)"E03",(string)"E04",(string)"E05",(string)"E06",(string)"E07",(string)"E08",(string)"E09",(string)"E10",
	//		 (string)"E11",(string)"E12",(string)"E13",(string)"E14",(string)"E15",(string)"E16",(string)"E17",(string)"E18",(string)"E19",(string)"E20",
	//		 (string)"E21",(string)"E22",(string)"E23",(string)"E24",(string)"E25",(string)"E26",(string)"E27",(string)"E28",(string)"E29",(string)"E30",
	//		 (string)"E31",(string)"E32",(string)"E33",(string)"E34",(string)"E35",(string)"E36" }},

	//	{BDS, { (string)"C01",(string)"C02",(string)"C03",(string)"C04",(string)"C05",(string)"C06",(string)"C07",(string)"C08",(string)"C09",(string)"C10",
	//		 (string)"C11",(string)"C12",(string)"C13",(string)"C14",(string)"C15",(string)"C16",(string)"C17",(string)"C18",(string)"C19",(string)"C20",
	//		 (string)"C21",(string)"C22",(string)"C23",(string)"C24",(string)"C25",(string)"C26",(string)"C27",(string)"C28",(string)"C29",(string)"C30",
	//		 (string)"C31",(string)"C32",(string)"C33",(string)"C34",(string)"C35",(string)"C36",(string)"C37",(string)"C38",(string)"C39" } },

	//	{SBS, { (string)"S20",(string)"S24",(string)"S25",(string)"S26",(string)"S27",(string)"S28",(string)"S29",
	//		 (string)"S33",(string)"S35",(string)"S36",(string)"S37",(string)"S38",
	//		 (string)"S40",
	//		 (string)"S83" }},

	//	{QZS, { (string)"J01",(string)"J02",(string)"J03",(string)"J07" }},

	//	{LEO, { (string)"201",(string)"202",(string)"203",(string)"204",(string)"205",(string)"206",(string)"207",(string)"208",
	//		 (string)"209",(string)"210",(string)"211",(string)"212",(string)"213",(string)"214",(string)"215",(string)"216",(string)"217",(string)"218",
	//		 (string)"219",(string)"220",(string)"221",(string)"222",(string)"223",(string)"224",(string)"225",(string)"226",(string)"227",(string)"228",
	//		 (string)"229",(string)"230",(string)"231",(string)"232",(string)"233",(string)"234",(string)"235",(string)"236",(string)"237",(string)"238",
	//		 (string)"239",(string)"240",(string)"241",(string)"242",(string)"243",(string)"244",(string)"245",(string)"246",(string)"247",(string)"248",
	//		 (string)"249",(string)"250",(string)"251",(string)"252",(string)"253",(string)"254",(string)"255",(string)"256",(string)"257",(string)"258",
	//		 (string)"259",(string)"260" }}, //HONG YAN LEO SATS PRN: 201 - 260 add leo}
	//	};                     // static map of default GNSS satellites              


	// const t_map_gnav GNSS_GNAV = {
	//		{GPS,{ (string)"NAV",(string)"CNAV"}},
	//		{GLO,{ (string)"NAV" }},
	//		{GAL,{ (string)"NAV",(string)"INAV_E01",(string)"INAV_E01",(string)"INAV",(string)"FNAV" }},
	//		{BDS,{ (string)"NAV" }},
	//		{SBS,{ (string)"NAV" }},
	//		{QZS,{ (string)"NAV" }},
	//};                      // static map of default GNSS navigation types

	 const t_map_freq GNSS_FREQ_PRIORITY = {
		{GPS , { LAST_GFRQ, G01,      G02,      G05}},
		{GLO , { LAST_GFRQ, R01,      R02,      R03_CDMA, R05_CDMA}},
		{GAL , { LAST_GFRQ, E01,      E05,      E07,      E08,      E06     }},
		{BDS , { LAST_GFRQ, C02,      C07,      C06,      C05,      C09,       C08,       C01     }},
		{QZS , { LAST_GFRQ, J01,      J02,      J05,      J06               }},
		{SBS , { LAST_GFRQ, S01,      S05                                   }},
		{GNS , {                                                            }},
	};	            // static map of default GNSS freq priorities
	
	 const t_map_band GNSS_BAND_PRIORITY= {
		{GPS,{ BAND, BAND_1, BAND_2, BAND_5                 }},
		{GLO,{ BAND, BAND_1, BAND_2, BAND_3, BAND_5			}},
		{GAL,{ BAND, BAND_1, BAND_5, BAND_7, BAND_8, BAND_6 }},
		{BDS,{ BAND, BAND_2, BAND_7, BAND_6, BAND_5, BAND_9, BAND_8, BAND_1}},
		{QZS,{ BAND, BAND_1, BAND_2, BAND_5, BAND_6         }},
		{SBS,{ BAND, BAND_1, BAND_5							}},
		{GNS,{												}},
	};	          // static map of default GNSS band priorities
	//BASE_LIBRARY_EXPORT const t_map_gnss GNSS_DATA_PRIORITY; 	          // static map of default GNSS data types/bands/attrs priorities



   // ------------------------------------------------------------------------------------------------------
   // STATIC MAPS
	// ------------------------------------------------------------------------------------------------------

   /*
	// c++11 (problem with Debian 5.0)
	// -------------------------------
	static t_map_gnss GNSS_DATA_PRIORITY = {
	 {GPS, { {BAND_1,{ {TYPE_C, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M,         ATTR_NULL }},
				   {TYPE_L, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_N, ATTR_NULL }},
				   {TYPE_D, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_N, ATTR_NULL }},
				   {TYPE_S, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_N, ATTR_NULL }},
					   {TYPE_P, { ATTR_NULL }}
			 }},
			 {BAND_2,{ {TYPE_C, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_D,         ATTR_NULL }},
				   {TYPE_L, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_D, ATTR_N, ATTR_NULL }},
				   {TYPE_D, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_D, ATTR_N, ATTR_NULL }},
				   {TYPE_S, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_D, ATTR_N, ATTR_NULL }},
					   {TYPE_P, { ATTR_NULL }}
			 }},
			 {BAND_5,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}
		 }},
			 {BAND_A,{ {TYPE_C, { ATTR_NULL }},
				   {TYPE_L, { ATTR_NULL }},
				   {TYPE_D, { ATTR_NULL }},
				   {TYPE_S, { ATTR_NULL }}
		 }},
			 {BAND_B,{ {TYPE_C, { ATTR_NULL }},
				   {TYPE_L, { ATTR_NULL }},
				   {TYPE_D, { ATTR_NULL }},
				   {TYPE_S, { ATTR_NULL }}
		 }}
	 }},

	 {GLO, { {BAND_1,{ {TYPE_C, { ATTR_C, ATTR_P, ATTR_NULL }},
				   {TYPE_L, { ATTR_C, ATTR_P, ATTR_NULL }},
				   {TYPE_D, { ATTR_C, ATTR_P, ATTR_NULL }},
				   {TYPE_S, { ATTR_C, ATTR_P, ATTR_NULL }}
					   {TYPE_P, { ATTR_NULL }}
			 }},
			 {BAND_2,{ {TYPE_C, { ATTR_C, ATTR_P, ATTR_NULL }},
				   {TYPE_L, { ATTR_C, ATTR_P, ATTR_NULL }},
				   {TYPE_D, { ATTR_C, ATTR_P, ATTR_NULL }},
				   {TYPE_S, { ATTR_C, ATTR_P, ATTR_NULL }}
					   {TYPE_P, { ATTR_NULL }}
			 }},
			 {BAND_3,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}
		 }},
			 {BAND_A,{ {TYPE_C, { ATTR_NULL }},
				   {TYPE_L, { ATTR_NULL }},
				   {TYPE_D, { ATTR_NULL }},
				   {TYPE_S, { ATTR_NULL }}
		 }},
			 {BAND_B,{ {TYPE_C, { ATTR_NULL }},
				   {TYPE_L, { ATTR_NULL }},
				   {TYPE_D, { ATTR_NULL }},
				   {TYPE_S, { ATTR_NULL }}
		 }}
	 }},

	 {GAL, { {BAND_1,{ {TYPE_C, { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z, ATTR_NULL }},
				   {TYPE_L, { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z, ATTR_NULL }},
				   {TYPE_D, { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z, ATTR_NULL }},
				   {TYPE_S, { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z, ATTR_NULL }}
			 }},
			 {BAND_5,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}
		 }},
			 {BAND_6,{ {TYPE_C, { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z, ATTR_NULL }},
				   {TYPE_L, { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z, ATTR_NULL }},
				   {TYPE_D, { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z, ATTR_NULL }},
				   {TYPE_S, { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z, ATTR_NULL }}
		 }},
			 {BAND_7,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}
		 }},
			 {BAND_8,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}
		 }},
			 {BAND_A,{ {TYPE_C, { ATTR_NULL }},
				   {TYPE_L, { ATTR_NULL }},
				   {TYPE_D, { ATTR_NULL }},
				   {TYPE_S, { ATTR_NULL }}
		 }},
			 {BAND_B,{ {TYPE_C, { ATTR_NULL }},
				   {TYPE_L, { ATTR_NULL }},
				   {TYPE_D, { ATTR_NULL }},
				   {TYPE_S, { ATTR_NULL }}
		 }}
	 }},

	 {BDS, { {BAND_1,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}, // WILL BE DELETED AFTER CONVERSION IN GCODERS
				   {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}, // WILL BE DELETED AFTER CONVERSION IN GCODERS
				   {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}, // WILL BE DELETED AFTER CONVERSION IN GCODERS
				   {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}  // WILL BE DELETED AFTER CONVERSION IN GCODERS
			 }},
			 {BAND_2,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}
			 }},
			 {BAND_6,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}
		 }},
			 {BAND_7,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
				   {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}
		 }}
	 }},

	 {SBS, { {BAND_1,{ {TYPE_C, { ATTR_C }},
				   {TYPE_L, { ATTR_C }},
				   {TYPE_D, { ATTR_C }},
				   {TYPE_S, { ATTR_C }},
					   {TYPE_P, { ATTR_NULL }}
			 }},
			 {BAND_5,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X }},
				   {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X }},
				   {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X }},
				   {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X }}
		 }}
	 }},

	 {QZS, { {BAND_1,{ {TYPE_C, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_Z }},
				   {TYPE_L, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_Z }},
				   {TYPE_D, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_Z }},
				   {TYPE_S, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_Z }}
			 }},
			 {BAND_2,{ {TYPE_C, { ATTR_S, ATTR_L, ATTR_X }},
				   {TYPE_L, { ATTR_S, ATTR_L, ATTR_X }},
				   {TYPE_D, { ATTR_S, ATTR_L, ATTR_X }},
				   {TYPE_S, { ATTR_S, ATTR_L, ATTR_X }}
		 }},
			 {BAND_5,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X }},
				   {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X }},
				   {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X }},
				   {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X }}
		 }},
			 {BAND_6,{ {TYPE_C, { ATTR_S, ATTR_L, ATTR_X }},
				   {TYPE_L, { ATTR_S, ATTR_L, ATTR_X }},
				   {TYPE_D, { ATTR_S, ATTR_L, ATTR_X }},
				   {TYPE_S, { ATTR_S, ATTR_L, ATTR_X }}
		 }}
	 }}
	};

   // c++11 (problem with Debian 5.0)
   // -------------------------------
   static map<GSYS, set<string> > GNSS_SATS =
   {
	 { GPS, { "G01","G02","G03","G04","G05","G06","G07","G08","G09","G10",
			  "G11","G12","G13","G14","G15","G16","G17","G18","G19","G20",
			  "G21","G22","G23","G24","G25","G26","G27","G28","G29","G30",
			  "G31","G32" } },

	 { GLO, { "R01","R02","R03","R04","R05","R06","R07","R08","R09","R10",
			  "R11","R12","R13","R14","R15","R16","R17","R18","R19","R20",
			  "R21","R22","R23","R24" } },

	 { GAL, { "E01","E02","E03","E04","E05","E06","E07","E08","E09","E10",
			  "E11","E12","E13","E14","E15","E16","E17","E18","E19","E20",
			  "E21","E22","E23","E24","E25","E26","E27","E28","E29","E30" } },

	 { BDS, { "C01","C02","C03","C04","C05","C06","C07","C08","C09","C10",
			  "C11","C12","C13","C14","C15","C16","C17","C18","C19","C20",
			  "C21","C22","C23","C24","C25","C26","C27","C28","C29","C30" } },

	 { SBS, { "S20","S24","S25","S26","S27","S28","S29",
			  "S33","S35","S36","S37","S38",
			  "S40",
			  "S83" } },

	 { QZS, { "J01" } }
   };
   */

} // namespace

#endif // GOBS_H

