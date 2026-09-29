/**
*
* @verbatim
	History
	 -1.0 Jiande Huang  2019-07-18  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gxmlfactory.h
* @brief		main control class for get xml
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-18
*
*/


#ifndef GXMLFACTORY_H
#define GXMLFACTORY_H

#include "gutils/ExportUtl.h"
#include "gset/gsetbase.h"
#include "gset/gsetinp.h"
#include "gio/glog.h"
#include <stdint.h>

using namespace gnut;

namespace great
{

	/**
	* @brief main control class for caculate clk
	*/
	class UTL_LIBRARY_EXPORT t_gxmlfactory
	{
	public:
		
		t_gxmlfactory();
		t_gxmlfactory(t_gsetbase* set, t_glog* log);

		virtual ~t_gxmlfactory();

		bool ProcessBatch();
	protected:
		t_gsetbase* _set;
		t_glog*     _log;
		set<string> _rec_list;
		set<string> _sat_list;
		set<string> _sys_list;
		t_gtime		_beg_time;
		t_gtime		_end_time;

		map<string, string>    _proj_list;
		map<string, string>    _xml_list;
		multimap<IFMT, string> _inp_list;

		string _xml_eg;
		string _app_type;
		string _sp3_type;
		string _dcb_root;
		string _sp3_root;
		string _brd_root;
		string _rinexc_type;
		string _sinex_root;
		string _rinexc_root;
		string _rinexo_root;
		string _sysfile_root;
		string _project_root;
		bool _scan_sysfile = false;
		bool _scan_rinexo = false;
		double _time_interval = 0.0;

		bool _isWIN = false;
		bool _isLINUX = false;

		bool _init_proj_dirs(string root, const t_gtime& beg_time, const t_gtime& end_time);
		bool _init_xml_file(string path, const t_gtime& crt_time);

		bool _init_podIF_xml_file(const t_gtime& crt_time, pugi::xml_node& parent);
		bool _init_pceIF_xml_file(const t_gtime& crt_time, pugi::xml_node& parent);//add by xiongyun
		bool _insert_gen(const t_gtime& crt_time,pugi::xml_node& parent);
		bool _insert_inp_POD_IF(const t_gtime& crt_time,pugi::xml_node& parent);
		bool _insert_inp_PCE_IF(const t_gtime& crt_time, pugi::xml_node& parent);//add by xiongyun

		bool _insert_inp_rnxo(const t_gtime& crt_time, pugi::xml_node& parent);
		bool _insert_inp_rnxc(const t_gtime& crt_time, pugi::xml_node& parent);
		bool _insert_inp_ics(const t_gtime& crt_time, pugi::xml_node& parent);
		bool _insert_inp_orb(const t_gtime& crt_time, pugi::xml_node& parent);
		bool _insert_inp_snx(const t_gtime& crt_time, pugi::xml_node& parent);
		bool _insert_inp_dcb(const t_gtime& crt_time, pugi::xml_node& parent);
		bool _insert_inp_sp3(const t_gtime& crt_time, pugi::xml_node& parent);
		bool _insert_inp_brd(const t_gtime& crt_time, pugi::xml_node& parent);
		bool _insert_inp_sys(pugi::xml_node& parent);

		bool _insert_par(const t_gtime& crt_time, pugi::xml_node& parent);
		bool _insert_par_PCE(const t_gtime& crt_time, pugi::xml_node& parent);//add by xiongyun
		bool _insert_rec(const t_gtime& crt_time, pugi::xml_node& parent);
		bool _insert_gnss(pugi::xml_node& parent);
		bool _insert_process(pugi::xml_node& parent);
		bool _insert_root();

		int32_t createDirectory(const std::string &directoryPath);
		//string  recs2rnxos(const t_gtime& crt_time, const set<string>& recs);
		//string  three_day_rnxcs(const t_gtime& crt_time);
		string  sgg_file(const t_gtime& crt_time,const string & type);

		string _list_atx	= "igs_absolute_14.atx";
		string _list_blq	= "oceanload";
		string _list_de_win = "jpleph_de405_win";
		string _list_de_lux = "jpleph_de405_lux";
		string _list_satpars = "sat_parameters";
		string _list_poleut1 = "poleut1";

		string _replace_YYYY(const t_gtime& crt_time, const string& path);
		string _replace_DOY(const t_gtime& crt_time, const string& path);

	};
	
}

#endif /*GCLKLSQ_H*/