
/**
*
* @verbatim
	 (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  @endverbatim
*
* @file		gsetout.h
* @brief	implements output setting class
* @author   Jan Dousa
* @version	1.0.0
* @date		2012-10-23
*
*/


#ifndef GSETOUT_H
#define GSETOUT_H

#define XMLKEY_OUT "outputs"        ///< The defination of outputs module in xml file

#define DEFAULT_FILE_VER ""         ///< default version for file format
#define DEFAULT_FILE_UPD     0      ///< default auto update [min] for file saving
#define DEFAULT_FILE_LEN     0      ///< default auto length [min] for file saving
#define DEFAULT_FILE_SMP     0      ///< default auto sample [sec] for file saving
#define DEFAULT_FILE_OFF     0      ///< default file offset [min] for file saving
#define DEFAULT_FILE_SYS  "UTC"     ///< default file system       for file saving

#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"
#include "gsetbase.h"

using namespace std;
using namespace pugi;

namespace gnut
{

	///< the order is important here!
	enum OFMT
	{
		XXX_OUT,
		TSA_OUT,
		MPT_OUT,
		LOG_OUT,
		SUM_OUT,
		XTR_OUT,
		NAV_OUT,
		XML_OUT,
		XQC_OUT,
		PPP_OUT,
		FLT_OUT,
		SMT_OUT,
		LSQ_OUT,
		RES_OUT,
		PREOBS_OUT,
		GRD_OUT,
		SRF_OUT,
		PRE_OUT,
		RB_OUT,
		STA_OUT,
		KBROMC_OUT,
		FIT_OUT,
		CFG_OUT,
		RINEXN_OUT,
		RINEXN2_OUT,
		RINEXO_OUT,
		RINEXC_OUT,
		RTCM_OUT,
		BNCOBS_OUT,
		BNCRTCM_OUT,
		SINEX_OUT,
		TROSINEX_OUT,
		TROSINEX0_OUT,
		KML_OUT,
		INS_OUT,
		INSKF_OUT,
		ENU_OUT,
		INSKFPK_OUT,
		PRDQC_OUT,
		ORB_OUT,         ///< orb out
		SP3_OUT,         ///< sp3 out
		RECCLK_OUT,      ///< rec clk out
		RECCLK13_OUT,    ///< rec clk out for freq 13(add by xiongyun)
		SATCLK_OUT,      ///< sat clk out
		SATCLK13_OUT,    ///< sat clk out for freq 13(add by xiongyun)
		UPD_OUT,         ///< upd out
		AMBUPD_OUT,
		AMBCONLEO_OUT,  ///< ambupd out
		RESULT_OUT,
		ICS_OUT,         ///< ics out
		ICSLEO_OUT,      ///< icsleo out
		ION_OUT,
		AMBCON_OUT,		 ///< ambcon out
		ORBDIF_OUT,		 ///< orbit differ out
		// Added by yqyuan
		ORBFIT_OUT,		 ///< orbit fit out						  
		RECOVER_OUT,
		AUG_OUT,
		TROP_OUT,
		AMBFLAG_DIR_OUT,
		IMU_OUT
	};

	/// The class for settings of output
	class BASE_LIBRARY_EXPORT t_gsetout : public virtual t_gsetbase
	{
	public:
		/// constructor
		t_gsetout();
		/// destructor
		~t_gsetout();

		/**
		 * @brief change from string to OFMT
		 * @param[in] s file format
		 * @return OFMT : file format
		 */
		static OFMT   str2ofmt(const string& s);

		/**
		 * @brief change from OFMT to string
		 * @param[in] f file format
		 * @return string : file format
		 */
		static string ofmt2str(const OFMT&   f);

		/// settings check
		void check();         
		/// settings help
		void help();                               

		// attributes
		/**
		 * @brief  get verbosity attribute
		 * @return int : verbosity attribute
		 */
		int verb();   

		/**
		 * @brief  get append request
		 * @return bool : append request
		 */
		bool append();       

		int sp3_obslimit();

		/// update glog (mask,verb)
		virtual void _upd_glog();                  

		// elements
		/**
		 * @brief  get format output size
		 * @param[in] fmt file format
		 * @return int : format output size
		 */
		int output_size(const string& fmt);    

		/**
		 * @brief  get string outputs
		 * @param[in] fmt file format
		 * @return string : string outputs
		 */
		string outputs(const string& fmt);


		/**
		 * output info via tcp port (added by zhshen)
		 * return the tcp port
		 */
		string output_port(const string& fmt);

		/**
		 * @brief  get formats
		 * @return set<string> : all the outputs
		 */
		set<string> oformats();                     

		/**
		 * @brief  get string output version
		 * @param[in] fmt file format
		 * @return string : string output version
		 */
		string version(const string& fmt);  

		/**
		 * @brief  get time offset [min] for output filename
		 * @param[in] fmt file format
		 * @return int : time offset [min] for output filename
		 */
		int             file_toff(const string& fmt);

		/**
		 * @brief  get time system for output filename
		 * @param[in] fmt file format
		 * @return t_tsys : time system  for output filename
		 */
		t_gtime::t_tsys file_tsys(const string& fmt);

		/**
		 * @brief get update [min] for output file update
		 * @param[in] fmt file format
		 * @return int : update [min] for output file update
		 */
		int   out_update(const string& fmt);       
		
		/**
		 * @brief  get length [min] for output file content
		 * @param[in] fmt file format
		 * @return int : length [min] for output file content
		 */
		int   out_length(const string& fmt);     
		
		/**
		 * @brief  get sample [sec] for output file data
		 * @param[in] fmt file format
		 * @return float : sample [sec] for output file data
		 */
		float out_sample(const string& fmt);   

	protected:

		/**
		 * @brief  get string output file 
		 * @param[in] fmt file format
		 * @return string : string output file 
		 */
		string _outputs(const string& fmt);

		/**
		 * @brief  get all string output file
		 * @return set<string> : all string output file
		 */
		set<string> _oformats();

		set<OFMT> _OFMT_supported;                  ///< vector of supported OFMTs (app-specific)
		bool  _append;                              ///< append mode
		int   _verb;                                ///< output verbosity   
		int   _upd;                                 ///< update [min] for output file update
		int   _len;                                 ///< length [min] for output file content
		float _smp;                                 ///< sample [sec] for output file data

	private:

	};

} // namespace


#endif
