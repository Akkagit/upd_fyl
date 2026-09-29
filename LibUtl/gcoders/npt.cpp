/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-10-24	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		npt.cpp
 * @brief		The base class used to decode npt file information.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-10-24
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <vector>


#include "gcoders/npt.h"
#include "gdata/gnpt.h"
#include "gall/gallnpt.h"


using namespace std;
namespace great
{
	t_npt::t_npt()
	{
	}

	t_npt::t_npt(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_npt::constructor");

	}

	t_npt::~t_npt()
	{

	}

	int t_npt::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		if (t_gcoder::_add2buffer(buff, sz) == 0) { return 0; };
		return -1;
	}

	int t_npt::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

#ifdef DEBUG
		cout << " BUFFER : \n" << _buffer << "\n size = " << sz << " END OF BUFFER \n\n"; cout.flush();
#endif

		int sitsize = 0;
		int tmpsize = 0;
		string line;

		int data_type;
		
		t_gtime _beg_epo;
		t_gtime _end_epo;

		string tmp;
		int tmpi;

		double temp_dsec = 0.0;

		bool lformer_first;
		bool llatter_first;
		double det_t = 0.0;
		double det_t_former = 0.0;
		double det_t_latter = 0.0;
		int i_former;
		int i_latter;

		_beg_epo.tsys(t_gtime::UTC);
		_end_epo.tsys(t_gtime::UTC);

		try
		{
			while ((tmpsize = t_gcoder::_getline(line, sitsize)) >= 0)
			{
				sitsize += tmpsize;

				if ((line.substr(0, 2)).compare("00") == 0) continue;

				if ((line.substr(0, 2)).compare("H1") == 0 || (line.substr(0, 2)).compare("h1") == 0)
				{
					connect = 0;
					istringstream istr(line);
					istr >> tmp;
					istr >> _crd_literal >> _format_version >> _prod_year >> _prod_mon >> _prod_day >> _prod_hour;
					continue;

				}
				//while ((tmpsize = t_gcoder::_getline(line, sitsize)) >= 0)
				//{
					//sitsize += tmpsize;
				if ((line.substr(0, 2)).compare("H8") == 0 || (line.substr(0, 2)).compare("h8") == 0)
				{
					tag = 0;
					//break;
				}
				else if ((line.substr(0, 2)).compare("00") == 0) continue;
				else if ((line.substr(0, 2)).compare("H2") == 0 || (line.substr(0, 2)).compare("h2") == 0)
				{
					istringstream istr(line);
					istr >> tmp;
					istr >> _stn_name;
					/*if (_stn_name.compare("YARL") == 0)
					{
						int test = 0;
					}*/

					if (tot_epoch == 198)
					{
						int test = 0;
					}
					istr >> _cdp_pad_id >> _cdp_sys_num >> _cdp_occ_num >> _stn_timescale;
					gnpt.setName(_stn_name);
					gnpt.setID(_cdp_pad_id);
					gnpt.setTimeScale(_stn_timescale);
					continue;
				}
				else if ((line.substr(0, 2)).compare("H3") == 0 || (line.substr(0, 2)).compare("h3") == 0)
				{
					istringstream istr(line);
					istr >> tmp;
					istr >> _target_name;
					istr >> _ilrs_id >> _sic >> _norad >> _SC_timescale >> _target_type;
					gnpt.setLEOid(_ilrs_id);
					continue;

				}
				else if ((line.substr(0, 2)).compare("H4") == 0 || (line.substr(0, 2)).compare("h4") == 0)
				{
					istringstream istr(line);
					istr >> tmp;
					istr >> data_type >> start_year >> start_mon >> start_day
						>> start_hour >> start_min >> start_sec
						>> end_year >> end_mon >> end_day >> end_hour >> end_min >> end_sec;

					istr >> _data_release >> _refraction_app_ind >> _CofM_app_ind
						>> _xcv_amp_app_ind >> _stn_sysdelay_app_ind
						>> _SC_sysdelay_app_ind >> _range_type_ind >> _data_qual_alert_ind;

					_beg_epo.from_ymdhms(start_year, start_mon, start_day, start_hour, start_min, start_sec,false);
					_end_epo.from_ymdhms(end_year, end_mon, end_day, end_hour, end_min, end_sec,false);

					if (_beg_epo.sod() == 36051)
					{
						int bcd = 1;
					}
					if (end_year == -1)
					{
						_end_epo.from_dmjd(_beg.dmjd() - 1);
					}
					if ((_beg_epo.dmjd() > _end.dmjd()) || ((_end_epo.dmjd() + 1.0 / 86400.0) < _beg.dmjd()))
					{
						while ((tmpsize = t_gcoder::_getline(line, sitsize)) >= 0)
						{
							sitsize += tmpsize;
							if ((line.substr(0, 2)).compare("H8") == 0 || (line.substr(0, 2)).compare("h8") == 0)
							{
								tag = 1;
								break;
							}

						}
						tag = 1;
						connect = 1;
						continue;
					}
					else
					{
						gnpt.now = _beg_epo;
						gnpt.setTrop(_refraction_app_ind);
						gnpt.setMass(_CofM_app_ind);

						temp_dsec = 0.0;

						//tag = 0;

						n_mete = 0;
						_t_mete.clear();
						_Atm_pressure.clear();
						_Dry_temperature.clear();
						_Atm_humidity.clear();

						n_data = 0;
						_dsec.clear();
						_two_way_time.clear();
						_Sigma_obs.clear();
						_Code_n_point_window.clear();
						continue;
					}

				}
				else if ((line.substr(0, 2)).compare("C0") == 0 || (line.substr(0, 2)).compare("c0") == 0)
				{
					istringstream istr(line);
					istr >> tmp;
					istr >> tmpi;
					istr >> _wavelength;
					gnpt.setWaveLength(_wavelength*0.001);   //nm => micormeter
					continue;
				}
				else if ((line.substr(0, 2)).compare("40") == 0)
				{
					istringstream istr(line);
					istr >> tmpi;
					istr >> _d40_sec_of_day >>
						_type_of_data >> _temp_sysconfig_id >>
						_num_points_recorded >> _num_point_used >>
						_one_way_target_dist >> _cal_sys_delay >> _cal_delay_shift >>
						_cal_rms >> _cal_skew >> _cal_kurtosis >> _cal_PmM;
					istr >> _cal_type_ind >>
						_cal_shift_type_ind >> _d40_detector_channel;

					gnpt.setSystemDelay(_cal_sys_delay);
					gnpt.setSigmaCalibration(_cal_rms);
					continue;
				}
				else if ((line.substr(0, 2)).compare("20") == 0)
				{
					n_mete += 1;
					istringstream istr(line);
					istr >> tmpi;
					istr >> _d20_sec_of_day >> _pressure >>
						_temperature >> _humidity >> _value_origin;
					_t_mete.push_back(_d20_sec_of_day);
					_Atm_pressure.push_back(_pressure);
					_Dry_temperature.push_back(_temperature);
					_Atm_humidity.push_back(_humidity * 0.01);
					continue;
				}
				else if ((line.substr(0, 2)).compare("11") == 0)
				{
					n_data += 1;
					istringstream istr(line);
					istr >> tmpi;
					istr >> _d11_sec_of_day >> _d11_time_of_flight >>
						_temp_sysconfig_id >> _d11_epoch_event >> _np_window_length >>
						_num_ranges >> _bin_rms >> _bin_skew >> _bin_kurtosis >> _bin_PmM >>
						_return_rate;

					_dsec.push_back(_d11_sec_of_day);
					_two_way_time.push_back(_d11_time_of_flight * 1e12);
					_Sigma_obs.push_back(_bin_rms);

					///< for Normal Point widow indicator
					if (abs(_np_window_length - 0.0) < 1e-8)
						_Code_n_point_window.push_back(0);
					else if (abs(_np_window_length - 5.0) < 1e-8)
						_Code_n_point_window.push_back(1);
					else if (abs(_np_window_length - 15.0) < 1e-8)
						_Code_n_point_window.push_back(3);
					else if (abs(_np_window_length - 20.0) < 1e-8)
						_Code_n_point_window.push_back(4);
					else if (abs(_np_window_length - 30.0) < 1e-8)
						_Code_n_point_window.push_back(5);
					else if (abs(_np_window_length - 60.0) < 1e-8)
						_Code_n_point_window.push_back(6);
					else if (abs(_np_window_length - 120.0) < 1e-8)
						_Code_n_point_window.push_back(7);
					else if (abs(_np_window_length - 180.0) < 1e-8)
						_Code_n_point_window.push_back(8);
					else if (abs(_np_window_length - 300.0) < 1e-8)
						_Code_n_point_window.push_back(9);
					else
					{
						_log->comment(0, "Error(get_slr_data_npt):Unknown Code_n_point_window");
					}
					continue;
				}
				else if ((line.substr(0, 2)).compare("10") == 0)
				{
					n_data += 1;

				}
				//}

				/*if (tmpsize == 0)
				{
					connect = 0;
				}*/

				if (tag == 0 && connect == 0)
				{
					for (int i = 0; i <= n_data - 1; i++)
					{
						lformer_first = true;
						llatter_first = true;
						det_t = 0.0;
						det_t_former = 0.0;
						det_t_latter = 0.0;
						i_former = -1;
						i_latter = -1;

						gnpt.setdsec(_dsec[i]);

						if ((_dsec[i] - temp_dsec) < 0)
						{
							//gnpt.now.from_dmjd(gnpt.now.dmjd(false) + 1, false);
							gnpt.now.from_ymd(start_year, start_mon, start_day + 1, 0, 0.0);
						}

						temp_dsec = _dsec[i];

						//gnpt.now.from_ymd(start_year, start_mon, start_day, 0, 0.0);
						//gnpt.now.from_dmjd(gnpt.now.mjd(false) + (gnpt.getdsec() / 86400.0), false);
						// corrected by yqyuan to keep consistent with NPT time system (UTC)
						//gnpt.now.from_mjd(gnpt.now.mjd(), (int)gnpt.getdsec(), gnpt.getdsec() - (int)gnpt.getdsec());
						if (_ifLEO == "NAV")
						{
							gnpt.now.from_mjd(gnpt.now.mjd(), (int)gnpt.getdsec(), gnpt.getdsec() - (int)gnpt.getdsec(), false);
						}
						else
						{
							gnpt.now.from_mjd(gnpt.now.mjd(), (int)gnpt.getdsec(), gnpt.getdsec() - (int)gnpt.getdsec());
						}
						

						double time22 = gnpt.now.dmjd();

						gnpt.setTwoWayTime(_two_way_time[i]);
						gnpt.setSigmaObs(_Sigma_obs[i]);
						gnpt.setSigma(0.02);
						gnpt.setCode(_Code_n_point_window[i]);

						if (n_mete == 1)
						{
							gnpt.setAtmPressure(_Atm_pressure[n_mete - 1]);
							gnpt.setDryTem(_Dry_temperature[n_mete - 1]);
							gnpt.setAtmHumidity(_Atm_humidity[n_mete - 1]);
						}
						else if (n_mete > 1)
						{

							for (int i_m = 0; i_m <= n_mete - 1; i_m++)
							{
								det_t = gnpt.getdsec() - _t_mete[i_m];
								if (det_t > 0.0)
								{
									if (lformer_first)
									{
										lformer_first = false;
										det_t_former = det_t;
										i_former = i_m;
									}
									else if (det_t < det_t_former)
									{
										det_t_former = det_t;
										i_former = i_m;
									}
								}
								else if (det_t < 0.0)
								{
									if (llatter_first)
									{
										llatter_first = false;
										det_t_latter = det_t;
										i_latter = i_m;
									}
									else if (det_t > det_t_latter)
									{
										det_t_latter = det_t;
										i_latter = i_m;
									}
								}
							}

							if (i_former != -1 && i_latter != -1)
							{
								gnpt.setAtmPressure(((_Atm_pressure[i_latter] - _Atm_pressure[i_former]) / (_t_mete[i_latter] - _t_mete[i_former]))*(gnpt.getdsec() - _t_mete[i_former]) + _Atm_pressure[i_former]);

								gnpt.setDryTem(((_Dry_temperature[i_latter] - _Dry_temperature[i_former]) / (_t_mete[i_latter] - _t_mete[i_former]))*(gnpt.getdsec() - _t_mete[i_former]) + _Dry_temperature[i_former]);

								gnpt.setAtmHumidity(((_Atm_humidity[i_latter] - _Atm_humidity[i_former]) / (_t_mete[i_latter] - _t_mete[i_former]))*(gnpt.getdsec() - _t_mete[i_former]) + _Atm_humidity[i_former]);

							}
							else if (i_former == -1 && i_latter != -1)
							{
								gnpt.setAtmPressure(_Atm_pressure[i_latter]);
								gnpt.setDryTem(_Dry_temperature[i_latter]);
								gnpt.setAtmHumidity(_Atm_humidity[i_latter]);

							}
							else if (i_former != -1 && i_latter == -1)
							{
								gnpt.setAtmPressure(_Atm_pressure[i_former]);
								gnpt.setDryTem(_Dry_temperature[i_former]);
								gnpt.setAtmHumidity(_Atm_humidity[i_former]);

							}
							else
							{
								_log->comment(0, "Error(get_slr_data_npt):11 There is no suitable meteorological data in Range Record");
							}
						}
						else
						{
							_log->comment(0, "Error(get_slr_data_npt):11 There is no suitable meteorological data in Range Record");
							continue;
						}

						if ((gnpt.now.dmjd() < _end.dmjd()) && (gnpt.now.dmjd() > _beg.dmjd()))
						{

							map<string, t_gdata*>::iterator it = _data.begin();

							if (tot_epoch == 198)
							{
								int abc = 0;
							}

							while (it != _data.end()) {

								if (it->second->id_type() == t_gdata::ALLNPT) {
									if (_ifLEO == "LEO")
									{
										int sizenpt = ((t_gallnpt*)it->second)->getSize(_sat);
										
										tot_epoch = sizenpt;
										
									}
									
									((t_gallnpt*)it->second)->setAllNpt(_sat,tot_epoch, gnpt);
								}

								it++;
							}

							tot_epoch = tot_epoch + 1;
						}
					}
					tag = 1;
				}


			}

		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_npt::decode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_npt::decode_head throw exception" << endl;
			}
			_mutex.unlock();
			return -1;
		}

		t_gcoder::_consume(sitsize);
		_mutex.unlock();
		return sitsize;
	}

	void t_npt::setBeg(t_gtime& beg)
	{
		_beg = beg;
	}

	void t_npt::setEnd(t_gtime& end)
	{
		_end = end;
	}

	void t_npt::setSat(string sat)
	{
		_sat = sat;
	}

	void t_npt::setType(string ifLEO)
	{
		_ifLEO = ifLEO;
	}
	
	string t_npt::NptName2LEO(string name)
	{
		string leo = "";
		if (name.find("swarma") != string::npos) {
			leo = "SWARM-A";
		}
		else if (name.find("swarmb") != string::npos) {
			leo = "SWARM-B";
		}
		else if (name.find("swarmc") != string::npos) {
			leo = "SWARM-C";
		}
		else if (name.find("gracea") != string::npos) {
			leo = "GRACE-A";
		}
		else if (name.find("graceb") != string::npos) {
			leo = "GRACE-B";
		}
		else if (name.find("gracec") != string::npos) {
			leo = "GRACE-C";
		}
		else if (name.find("graced") != string::npos) {
			leo = "GRACE-D";
		}
		else if (name.find("sentinel3a") != string::npos) {
			leo = "SENTINEL-3A";
		}
		else if (name.find("sentinel3b") != string::npos) {
			leo = "SENTINEL-3B";
		}
		else if (name.find("jason3") != string::npos) {
			leo = "JASON-3";
		}
		else if (name.find("jason2") != string::npos) {
			leo = "JASON-2";
		}
		else if (name.find("haiyang2a") != string::npos) {
			leo = "HAIYANG-2A";
		}
		else if (name.find("tiangong2") != string::npos) {
			leo = "TIANGONG-2";
		}
		return leo;
	}

	// added by yqyuan
	string t_npt::NptName2Prn(string name)
	{
		string prn = "";
		if (name.find("gps") != string::npos) {
			if (name.find("gps35") != string::npos) {
				// to be added
			}
			else if(name.find("gps36") != string::npos) {
				// to be added
			}
		}
		else if (name.find("glonass") != string::npos) {
			// to be added
		}
		else if (name.find("galileo") != string::npos) {
			if (name.find("galileo101") != string::npos || name.find("galileo-101") != string::npos) {
				prn = "E11";
			}
			else if (name.find("galileo102") != string::npos || name.find("galileo-102") != string::npos) {
				prn = "E12";
			}
			else if (name.find("galileo103") != string::npos || name.find("galileo-103") != string::npos) {
				prn = "E19";
			}
			else if (name.find("galileo104") != string::npos || name.find("galileo-104") != string::npos) {
				prn = "E20";
			}
			else if (name.find("galileo201") != string::npos || name.find("galileo-201") != string::npos) {
				prn = "E18";
			}
			else if (name.find("galileo202") != string::npos || name.find("galileo-202") != string::npos) {
				prn = "E14";
			}
			else if (name.find("galileo203") != string::npos || name.find("galileo-203") != string::npos) {
				prn = "E26";
			}
			else if (name.find("galileo204") != string::npos || name.find("galileo-204") != string::npos) {
				prn = "E22";
			}
			else if (name.find("galileo205") != string::npos || name.find("galileo-205") != string::npos) {
				prn = "E24";
			}
			else if (name.find("galileo206") != string::npos || name.find("galileo-206") != string::npos) {
				prn = "E30";
			}
			else if (name.find("galileo207") != string::npos || name.find("galileo-207") != string::npos) {
				prn = "E07";
			}
			else if (name.find("galileo208") != string::npos || name.find("galileo-208") != string::npos) {
				prn = "E08";
			}
			else if (name.find("galileo209") != string::npos || name.find("galileo-209") != string::npos) {
				prn = "E09";
			}
			else if (name.find("galileo210") != string::npos || name.find("galileo-210") != string::npos) {
				prn = "E01";
			}
			else if (name.find("galileo211") != string::npos || name.find("galileo-211") != string::npos) {
				prn = "E02";
			}
			else if (name.find("galileo212") != string::npos || name.find("galileo-212") != string::npos) {
				prn = "E03";
			}
			else if (name.find("galileo213") != string::npos || name.find("galileo-213") != string::npos) {
				prn = "E04";
			}
			else if (name.find("galileo214") != string::npos || name.find("galileo-214") != string::npos) {
				prn = "E05";
			}
			else if (name.find("galileo215") != string::npos || name.find("galileo-215") != string::npos) {
				prn = "E21";
			}
			else if (name.find("galileo216") != string::npos || name.find("galileo-216") != string::npos) {
				prn = "E25";
			}
			else if (name.find("galileo217") != string::npos || name.find("galileo-217") != string::npos) {
				prn = "E27";
			}
			else if (name.find("galileo218") != string::npos || name.find("galileo-218") != string::npos) {
				prn = "E31";
			}
			else if (name.find("galileo219") != string::npos || name.find("galileo-219") != string::npos) {
				prn = "E36";
			}
			else if (name.find("galileo220") != string::npos || name.find("galileo-220") != string::npos) {
				prn = "E13";
			}
			else if (name.find("galileo221") != string::npos || name.find("galileo-221") != string::npos) {
				prn = "E15";
			}
			else if (name.find("galileo222") != string::npos || name.find("galileo-222") != string::npos) {
				prn = "E33";
			}
			else {
				;
			}
		}
		else if (name.find("compass") != string::npos) {
			if (name.find("compassg1") != string::npos) {
				prn = "C01";
			}
			else if (name.find("compassi3") != string::npos) {
				prn = "C08";
			}
			else if (name.find("compassi5") != string::npos) {
				prn = "C10";
			}
			else if (name.find("compassm3") != string::npos) {
				prn = "C11";
			}
			else if (name.find("compassi6b") != string::npos) {
				prn = "C13";
			}
			else {
				;
			}
		}
		else if (name.find("beidou3") != string::npos) {
			if (name.find("beidou3m2") != string::npos) {
				prn = "C20";
			}
			else if (name.find("beidou3m3") != string::npos) {
				prn = "C21";
			}
			else if (name.find("beidou3m9") != string::npos) {
				prn = "C29";
			}
			else if (name.find("beidou3m10") != string::npos) {
				prn = "C30";
			}
			else {
				;
			}
		}
		else if (name.find("qzs") != string::npos) {
			if (name.find("qzs1") != string::npos) {
				prn = "J01";
			}
			else if (name.find("qzs2") != string::npos) {
				prn = "J02";
			}
			else if (name.find("qzs3") != string::npos) {
				prn = "J07";
			}
			else if (name.find("qzs4") != string::npos) {
				prn = "J03";
			}
			else {
				;
			}
		}
		else {
			;
		}

		return prn;
	}
}