/**
*
* @verbatim
	 (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  @endverbatim
*
* @file		gsetinp.h
* @brief	implements input setting class
* @author   Jan Dousa
* @version	1.0.0
* @date		2012-10-23
*
*/


#ifndef GSETINP_H
#define GSETINP_H

#define XMLKEY_INP "inputs"  ///< The defination of input node in XML

#include <map>
#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gutils/gtypeconv.h"
#include "gset/gsetbase.h"

using namespace std;
using namespace pugi;

namespace gnut
{

	/// input file order is important here!
	enum IFMT
	{
		FMT_INP,      
		IONEX_INP,
		SP3_INP,      
		SP3_REF,
		SP3_PRD,      
		ATX_INP,      ///< atx file, need to be before all OBSERVATIONS to create objects
		RINEXC_INP,   
		RINEXC_REF,
		RINEXC_PRD,
		RINEXM_INP,
		RINEXN_INP,
		RINEXN_REF,
		RINEXN_PRD,
		RINEXO_INP,     ///< rinexo file, need to be after RINEXN (for GLONASS chanels)!
		BNCOBS_INP,     ///< need to be after RINEXN (for GLONASS chanels)!
		BNCCORR_INP,	///< added by zhShen on 20191025
		BNCBRDC_INP,	///< added by zhShen on 20191025
		BNCSEND_INP,    ///< added by xiongyun
		TROSINEX_INP,
		TROSINEX0_INP,
		SINEX_INP,
		RTCM_INP,
		TSMET_INP,
		BLQ_INP,
		BNCRTCM_INP,    // need to be after RINEXN to fill the corrections !
		BNCRTCM_REF,
		BNCRTCM_PRD,
		TROGRID_INP,
		FCB_INP,

		GRIB_CHMI_INP,
		GRIB_HARM_INP,
		GRIB_ERA_INP,
		NCD3_ICS_INP,

		BIABERN_INP,
		BIASINEX_INP,

		RAO_TEXT_INP,
		RAO_IGRA_INP,
		RAO_BADC_INP,

		IMU_INP,			///< imu file for inertial navigation system (by zhshen)
		IMU_POS_INP,		///< GNSS positon and velocity file for loose-coupled integration (by zhshen)
		IMAGE_INP,
		//UPD_EPOCH_INP,
		LIDAR_INP,
		UPD_EPOCH_INP,
		ION_INP,            ///< ion file for gns
		ICS_INP,			///< ics file for gns
		ICSLEO_INP,         ///< ics file for leo
		DE_INP,             ///< DE405 file
//      OCEANLOAD_INP,	    ///< oceanload file (jdhuang,blq == oceanload)
		OCEANTIDE_INP,      ///< oceantide file
		POLEUT1_INP,        ///< poleut1 file
		LEAPSECOND_INP,     ///< leapsecond file
		EGM_INP,            ///< EGM file
		PANNEL_INP,         ///< pannel file for leo
		ORB_INP,            ///< orb file for gns
		ORBLEO_INP,         ///< orb file for leo 
		LEO_PANNEL_INP,     ///< pannel file for leo
		SOLAR_INP,          ///< solar parameters file
		ATTITUDE_INP,       ///< attitude file
		RB_INP,
		STA_INP,
		SATPARS_INP,        ///< sat parameters file 
		LRA_INP,
		NPT_INP,
		KBR_INP,
		SLRF_INP,
		ALBEDO_INP,
		AMBCON_INP,			/// ambcon file
		AMBCONLEO_INP,			/// ambcon file
		AMBINP_INP,			/// ambinp file
		AMBUPD_INP,          ///< ambupd file
		AMBFLAG12_INP,         ///< ambflag file
		AMBFLAG13_INP,       ///< ambflag13 file
		AMBFLAG14_INP,       ///< ambflag14 file
		AMBFLAG15_INP,       ///< ambflag15 file
		UPD_INP ,             ///< upd file
		DESAISCOPOLECOEF_INP,
		RECOVER_INP,       /// RESFILE
		PREOBS_INP,
		IFCB_INP,
		PANDA_RECOVER_INP,
		AUG_INP,
		UNDEF


	};

	/// The class for input settings
	class BASE_LIBRARY_EXPORT t_gsetinp : public virtual t_gsetbase
	{
	public:
		/// constructor
		t_gsetinp();
		/// destructor
		~t_gsetinp();

		/// settings check
		void check();      
		/// settings help
		void help();                                  

		/**
		 * @brief change from string to IFMT
		 * @param[in] s file format
		 * @return IFMT : file format
		 */
		static IFMT   str2ifmt(const string& s);

		/**
		 * @brief change from IFMT to string
		 * @param[in] f file format
		 * @return string : file format
		 */
		static string ifmt2str(const IFMT&   f);

		/**
		 * @brief get format input size
		 * @param[in] fmt file format
		 * @return int : format input size
		 */
		int input_size(const string& fmt);

		// jdhuang : add for check input
		bool check_input(const string& fmt);

		/**
		 * @brief get format inputs (all in multimap)
		 * @return multimap<IFMT, string> : format inputs (all in multimap)
		 */
		multimap<IFMT, string> inputs_all();  

		// jdhuang
		multimap<IFMT, string> inputs_all(t_gtime ref);        
		multimap<IFMT, string> replace_all(t_gtime ref, const set<string>& sats, const set<string>& recs, const multimap<IFMT, string>& inputs);


		/**
		 * @brief get format inputs (ordered)
		 * @param[in] fmt file format
		 * @return vector<string> : format inputs (ordered)
		 */
		vector<string> inputs(const string& fmt); 
		vector<string> inputs(const IFMT&   ifmt); 

		/**
		 * @brief get formats       (ordered)
		 * @return set<string> :    formats (ordered)
		 */
		set<string> iformats();   

		/**
		 * @brief get rtcm correction stream
		 * @return string : rtcm correction stream
		 */
		string corrStream(); 

		string rinexo();

		string bncobs();

		string bncbrdc();

		string bnccorr();

		/// checking navigation Rinex
		bool chkNavig();            
		/// checking satellite healthy status
		bool chkHealth();                             

	protected:

		/**
		 * @brief get all the file name of ftm
		 * @param[in] ftm file format
		 * @return vector<string> : all the file name of ftm
		 */
		vector<string> _inputs(const string& fmt);

		// jdhuang :
		// add for get fmt file
		vector<string> _inputs(const   IFMT& fmt);

		/**
		 * @brief get all the ftm in input node
		 * @return set<string> : all the ftm in input node
		 */
		set<string> _iformats();

	protected:
		set<IFMT> _IFMT_supported;                    ///< vector of supported IFMTs (app-specific)

		bool      _chkNavig;                          ///< check navigation data or not
		bool      _chkHealth;                         ///< check data health or not

		string    _corrStream;                        ///< error message string

	private:

	};

} // namespace

#endif
