/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-10-24	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		npt.h
 * @brief		The base class used to decode npt file information.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-10-24
 */
#ifndef NPT_H
#define NPT_H

#include "gcoders/gcoder.h"
#include "gutils/ExportUtl.h"
#include "gdata/gnpt.h"

using namespace std;
using namespace gnut;
namespace great
{

	class UTL_LIBRARY_EXPORT t_npt : public t_gcoder
	{
	public:
		/** @brief default constructor. */
		t_npt();

		/**
		 * @brief default constructor.
		 *
		 * @param[in]  s        setbase control
		 * @param[in]  version  version of the gcoder
		 * @param[in]  sz       size of the buffer
		 */
		t_npt(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_npt();

		/**
		 * @brief decode the data body of the npt data file.
		 *
		 * decode data body of npt file, all the data read will store in the gnpt.h
		 *
		 * @param[in]  buff        buffer of the data
		 * @param[in]  sz          buffer size of the data
		 * @param[in]  errmsg      error message of the data decoding
		 * @return
			@retval >=0 consume size of body decoding
			@retval <0  finish reading
		 */
		virtual int decode_head(char* buff, int sz, vector<string>& errmsg);

		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

		void setBeg(t_gtime& beg);

		void setEnd(t_gtime& end);

		void setType(string ifLEO);

		void setSat(string sat);

		string NptName2LEO(string name);
		// added by yqyuan
		string NptName2Prn(string name);

	protected:
		string _crd_literal;     ///< Consolidated Ranging Data format
		int _format_version;     ///< Format Version
		int _prod_year;          ///< Year of file production
		int _prod_mon ;           ///< Month of file production
		int _prod_day;           ///< Day of file production
		int _prod_hour;         ///< Hour of file production (UTC)

		string _stn_name;        ///< Station name from official list
		int	_cdp_pad_id;         ///< Crustal Dynamics Project Pad Identifier
		int _cdp_sys_num;        ///< Crustal Dynamics Project 2-digit system number
		int _cdp_occ_num;        ///< Crustal Dynamics Project 2-digit occupancy sequence number
		int _stn_timescale;       ///< Station Epoch Time Scale - indicates the time scale reference.

		string _target_name;     ///< Target name from official list
		int _ilrs_id;            ///< ILRS Satellite Identifier (Based on COSPAR ID)
		int _sic;                ///< SIC (Satellite Identification Code)
		int _norad;              ///< NORAD ID
		int _SC_timescale;       ///< Spacecraft Epoch Time Scale (transponders only)
		int _target_type;        ///< Target type
		t_gtime _beg;            ///< begin time in xml
		t_gtime _end;            ///< end time in xml
		string  _ifLEO;

		int _data_release;        ///< A flag to indicate the data release:
		int _refraction_app_ind;    ///< Tropospheric refraction correction applied indicator
		int _CofM_app_ind;          ///< Center of mass correction applied indicator
		int _xcv_amp_app_ind;       ///< Receive amplitude correction applied indicator
		int _stn_sysdelay_app_ind;    ///< Station system delay applied indicator
		int _SC_sysdelay_app_ind;      ///< Spacecraft system delay applied (transponders) indicator
		int _range_type_ind;          ///<  Range type indicator
		int _data_qual_alert_ind;         ///< Data quality alert indicator

		double _wavelength;

		int _cal_type_ind;
		int _cal_shift_type_ind;
		int _d40_detector_channel;

		double _d40_sec_of_day;
		int _type_of_data;
		string _temp_sysconfig_id;
		int _num_points_recorded;
		int _num_point_used;
		double _one_way_target_dist;
		double _cal_sys_delay;
		double _cal_delay_shift;
		double _cal_rms;
		double _cal_skew;
		double _cal_kurtosis;
		double _cal_PmM;

		double _d20_sec_of_day;
		double _pressure;
		double _temperature;
		double _humidity;
		double _value_origin;

		vector<double> _t_mete;
		vector<double> _Atm_pressure;
		vector<double> _Dry_temperature;
		vector<double> _Atm_humidity;

		vector<double> _dsec;
		vector<double> _two_way_time;
		vector<double> _Sigma_obs;
		vector<int> _Code_n_point_window;

		double _d11_sec_of_day;                 ///< Seconds of day
		double _d11_time_of_flight;             ///< Time of flight in seconds
		int _d11_epoch_event;                   ///< Epoch Event - indicates the time event reference.
		double _np_window_length;               ///< Normal point window length (seconds)
		int _num_ranges;                        ///< Number of raw ranges (after editing) compressed into the normal point.
		double _bin_rms;                        ///< Bin RMS from the mean of raw accepted time of flight values minus the trend function (ps)
		double _bin_skew;                       ///< Bin skew from the mean of raw accepted time of flight values minus the trend function.
		double _bin_kurtosis;                   ///< Bin kurtosis from the mean of raw accepted time of flight values minus the trend function.
		double _bin_PmM;                        ///< Bin peak ¨C mean value (ps)
		double _return_rate;                    ///< Return rate (%) for SLR or signal to noise ratio for LLR.
		int tot_epoch = 0;
		t_gnpt gnpt;
		int tag = 1;
		int connect = 1;

		int n_mete = 0;
		int n_data = 0;

		int start_year, start_mon, start_day, start_hour, start_min, start_sec;
		int end_year, end_mon, end_day, end_hour, end_min, end_sec;
		
		string _sat;
	};
}

#endif