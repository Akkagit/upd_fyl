/**
*
* @verbatim
History
-1.0 BoWong  2019-02-22  creat the file.
-1.1 glfeng  2019-03-02  Adjusted some functions and variables.
-1.2 BoWong  2019-04-09  Adding Doxygen Style Code Remarks.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gupdlsq.h
* @brief	compute Wide-Lane, Ex-Wide-Lane and Narrow-Lane upd.
*
* @author   BoWong, Wuhan University
* @version	1.0.0
* @date		2019-04-09
*
*/
#ifndef GUPDLSQ_H
#define GUPDLSQ_H

#include <algorithm>
#include "gutils/ExportLSQ.h"
#include "gutils/gmutex.h"
#include "gset/gsetbase.h"
#include "gall/gallobs.h"
#include "gall/gallambupd.h"
#include "gall/gallambflag.h"
#include "gdata/gupd.h"
#include "gall/gallbias.h"
#include "gmodels/gprecisemodel.h"
#include "gdata/gobsgnss.h"
#include "gambfix/gambcommon.h"
#include "gset/gcfg_upd.h"
#include "gcoders/upd.h"
#include "gcoders/ambflag.h"
#include "gio/gfile.h"
#include "gutils/ginfolog.h"
#include "gproc/gqualitycontrol.h"

namespace great
{
	/**
	* @brief class for upd calculate main function.
	*/
	class LSQ_LIBRARY_EXPORT t_gupdlsq 
	{
	public:

		/**
		* @brief default constructor, init some internal variables.
		* @param[in]  mark	upd mode.
		* @param[in]  gset	upd setting of xml.
		*/
		t_gupdlsq (UPDTYPE updmode, bool carrier_range, t_gsetbase* set, t_gallproc* data, t_glog* log);

		/** @brief default destructor. */
		virtual ~t_gupdlsq();

		/** @brief output UPD/IFCB File. */
		bool GenerateProduct();

		/**
		* @brief batch processing loop over epoch.
		* @param[in] beg   begin epoch.
		* @param[in] end   end epoch.
		* @param[in] intv  sampling interval.
		* @return	 true - successfully; false - failure
		*/
		bool processBatch(const t_gtime& beg, const t_gtime& end,const double &intv);

		///**
		//* @brief UPD mode to string .
		//* @param[in] type   upd mode(UPDTYPE).
		//* @return	 WL , EWL , NL change to string successfully
		//*            UNDEF - failure
		//*/
		//string UPDType2String(UPDTYPE type);

		void set_sites(set<string> sites) { _sites = sites; };
	protected:

#ifdef BMUTEX
		boost::mutex  _mutex;
#endif
		t_gmutex            _gmutex;

		/** @brief UPD mode:EWL/WL/NL. */
		UPDTYPE             _updmode;

		/** @brief whether write carrier-range. */
		bool                _carrier_range;

		/** @brief Genereal log output. */
		t_glog*				_glog;         

		 /** @brief Base setting. */
		t_gsetbase*			_set; 

		/** @brief Data used to calulate. */
		t_gallproc*         _data = nullptr;           ///< all data
		t_gallobs*          _allobs = nullptr;         ///< obs data
		t_gallambflag*      _allambflag = nullptr;     ///< ambflag data
		t_gallambflag*      _allambflag13 = nullptr;   ///< ambflag13 data for EWL
		t_gallambflag*      _allambflag14 = nullptr;   ///< ambflag14 data for EWL24
		t_gallambflag*      _allambflag15 = nullptr;   ///< ambflag15 data for EWL25
		t_gallambupd*       _allambupd = nullptr;      ///< ambupd data
		t_gallambupd*       _allfixedamb = nullptr;     ///< all fixed ambiguity
		t_gallbias*			_allbias = nullptr;        ///< dcb data
		t_gupd*             _allWLUPD = nullptr;       ///< WL upd, used for nl estimate
		t_gallnav*          _allnav = nullptr;        ///< used in glonass process
		/** @brief Result. */
		t_gupd*             _cal_UPD = nullptr;           ///< upd computed after updlsq
		t_gifcb*            _gifcb = nullptr;

		/** @brief Output for amb res files. */
		t_epoch_amb_res _AMB_RES;          ///< ambiguity residual data
		t_giof*	        _fambres;          ///< output file of ambres

		OBSCOMBIN _obs_type;               ///< IF/UDUC

		/** @brief Combined wavelength. */	
		map<string, int> _glo_freq_num;
		map<string, map<string, double>> _sys_wavelen;  

		/** @brief Sites & Sats & ParList */
		set<string>     _sites;           ///< sites list
		set<string>     _sats;            ///< satellites list
		string          _sys;             ///< system list(string)
		vector<string>  _parlist;         ///< params list

		/** @brief reference site . */
		string _iref_site;     ///< site reference
		string _jref_site;     ///< site reference
		
		/** @brief ambiguity of satellite-station pair == xx */
		map<string, map<string, double>>  _sat_site_amb;

		/** @brief use for WL/EWL estimation: record last sat of OW per site. */
		string _lastSat;  

        /** @brief use for judge whether ambupd file is end. */        
		map<string, t_gtime> _site_ambupd_end;

		//map<GSYS, vector<GOBSBAND>>   _SYS_bands;
		map< GSYS, map<FREQ_SEQ, GOBSBAND> > _band_index;

		// Quality_control
		t_gbds_codebias_cor   _bds_codebias_cor;
		t_goutliers_process   _outliers_proc;

	protected:

		bool _processWL(const t_gtime& beg, const t_gtime& end, const double &intv);
		bool _processNL(const t_gtime& beg, const t_gtime& end, const double &intv);
		bool _processIFCB(const t_gtime& beg, const t_gtime& end, const double &intv);

		/**
		* @brief combine wavelengths of WL/EWL/NL .
		* @param[in] 
		* @return    true combine wavelength successfully, false failure
		*/
		bool _setWavelength();

		/**
		* @brief Init data.
		* @param[in]
		* @return    return false if no data
		*/
		bool _initUPDProcData(t_gallproc * data);

		/**
		* @brief set sites and satellites and init params according upd mode.
		* @return	 true - successful, false failure
		*/
		bool _initSatSite();

		/**
		* @brief set output file
		* @return	 true - successful, false failure
		*/
		void _setOut();

		/**
		* @brief combine observations of the EWL/WL.
		* @param[in] site          site name.
		* @param[in] T_beg         start epoch.
		* @param[in] T_end         end epoch.
		* @param[in] intv          sampling rate in seconds.
		* @param[in] one_site_amb  widelane combination observations.
		* @param[in] OW            oneway ambiguities.
		* @param[out] one_site_amb
		* @param[out] OW
		* @return    true form (ex-)widelane observations successfully, false failure
		*/
		bool _combineObsIF(string site, const t_gtime& T_beg, const t_gtime& T_end, double intv,
			t_epo_sat_amb& one_site_amb, t_all_sats_OW& OW);
		bool _combineObsALL(string site, const t_gtime& T_beg, const t_gtime& T_end, double intv,
			t_epo_sat_amb& one_site_amb, t_all_sats_OW& OW);
		bool _combineObsGFIF(string site, const t_gtime& T_beg, const t_gtime& T_end, double intv,
			t_epo_sat_amb& one_site_GFIF, t_all_sats_OW& OW);

		/**
		* @brief calculate wide-lane/extra-wide-lane ambiguities for UPD estimation.
		* @param[in] intv          sampling rate in seconds.
		* @param[in] one_site_amb  widelane combination observations.
		* @param[in] OW            oneway ambiguities.
		* @param[out] OW
		* @return    true calculate successfully, false failure
		*/
		bool _calculateAmbWL(t_epo_sat_amb& one_site_amb, t_all_sats_OW& OW);

		/**
		* @brief apply satellite wide-lane UPDs to the computed wide-lane integer ambiguities,
		*	    then convert LC ambiguities to narrow-lane ambiguties with the integer wide-lane ambiguities.
		* @param[in] site   site name.
		* @param[in] t      current epoch.
		* @param[in] AMB    ambupd value from gambupd data.
		* @return    true calculate successfully, false failure
		*/
		bool _calculateAmbNL(string site, const t_gtime& t, t_epoch_amb& AMB);

		/**
		* @brief connect WL/EWL ambiguities, make each satellite-station pair have only one ambiguity.
		* @param[in] site   site name.
		* @param[in] OW     oneway ambiguities..
		* @param[in] AMB    ambupd value from gambupd data.
		* @return    true calculate successfully, false failure
		*/
		bool _connectAmb(string site, t_all_sats_OW& OW);

		/**
		* @brief Select the station with the most observation satellites as the reference station.
		* @param[in] t   current epoch.
		* @return    true calculate successfully, false failure
		*/
		bool _setRefStation(const t_gtime& t);   

		/**
		* @brief Get initial values of UPDs from UD-ambiguities based on a given reference station.
		* @param[in] epoch   current epoch.
		* @return    true calculate successfully, false failure
		*/
		bool _initUPD(const t_gtime& epoch);   

		/**
		* @brief estimate upd of each satellite and station.
		* @param[in] t   current epoch.
		* @return    true calculate successfully, false failure
		*/
		bool _updEstimate(const t_gtime& t);

		/**
		* @brief Output AMB_RES files in NL UPD Estimation.
		* @param[in] t   current epoch.
		* @return    true write successfully, false failure
		*/
		bool _writeAmbResidual(const t_gtime& t);
		void _readOfile(set<string>::iterator it_site, t_gallobs& obs);

		/**
		* @brief Write carrier-range in NL UPD Estimation.
		* @param[in] t   current epoch.
		* @return    true write successfully, false failure
		*/
		bool _writeCarrierRange(const t_gtime &t_beg, const double &intv);
	private:

	};


}
#endif //
