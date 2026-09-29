/**
*
* @verbatim
	History
	 -1.0 GREAT	    2019-01-04 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			  gturboedit.h
* @brief		  CRS coordinate transform to ACR coordinate
*
* @author         GREAT, Wuhan University
* @version		  1.0.0
* @date			  2019-01-04
*
*/

#ifndef GTURBOEDIT_H
#define GTURBOEDIT_H

#include <string>
#include <map>
#include <unordered_map>
#include "gutils/gcycleslip.h"
#include "gdata/gsatdata.h"
#include "gall/gallambflag.h"
#include "gset/gsetturboedit.h"
#include "gutils/ExportBase.h"

namespace great
{
	class BASE_LIBRARY_EXPORT t_gturboedit : public t_gcycleslip
	{
	public:
		t_gturboedit();
		t_gturboedit(t_gsetbase * gset, t_glog * glog, int index = 2);
		virtual ~t_gturboedit();
		
		void set_amb_flag(const string& rec, const string& sat, int flag);
		int  get_amb_flag(const string& rec, const string& sat);
		int  get_active_amb(string site);
		void set_active_amb(string site, int active_num);                         // glfeng
		int  num_of_amb_arc(const string& site, const string& prn, const t_gtime& time);
		bool use_of_obs(const string& site, const string& prn, const t_gtime& time);
		bool is_carrier_range(const string& site, const string& prn, const t_gtime& time);
		bool is_carrier_range(const string& site, const string& prn, const t_gtime& time, double &c1, double &c2);
		bool cycle_slip(t_gsatdata& obsdata, const t_gtime& time);
		bool cycle_slip123(t_gsatdata& obsdata, const t_gtime& time);

		set<string> get_sitelist_of_logfile() const;
		bool new_amb(const string& rec, const string& sat);
		void set_new_amb(const string& rec, const string& sat, bool isNew);
		t_gtime get_crt_amb_end(const string& rec, const string& sat);
		virtual map<string, int  >		get_active_amb() { return _active_amb; }
		void add_ambflag(string site, string sat, string description, t_gtime beg, t_gtime end);  // glfeng

		map<string, bool>& logfile_exist() { return _amb_info_file_exist; }
		void merge_logfile_exist(const map<string, bool>& logfile);

	protected:
		/** @brief read log file and record cycle slip info [panda file format]*/
		void _read_logfie(const set<string>& rec, const t_gtime& epoch, int index);
		void _read_logfile(const set<string>& rec, const t_gtime& epoch, int index);
		map<string, int  > _active_amb;              ///< Max ambc in one epoch  glfeng
		map<string, bool > _amb_info_file_exist;		     ///< recording exitence of panda logfile
		unordered_map <string, unordered_map <string, int> > _amb_flag;	 ///< record the site-sat amb arc
		unordered_map <string, unordered_map< string, vector<pair<t_gtime, t_gtime> > > > _cycle_flag;				///< recording cycle slip info in panda log file
		unordered_map <string, unordered_map< string, vector<pair<t_gtime, t_gtime> > > > _cycle_flag_unused;		///< recording unused cycle slip info in panda log file
		unordered_map <string, unordered_map <string, bool> > _new_amb;

		t_gsetbase* _gset = nullptr;
		t_glog*     _glog = nullptr;
		shared_ptr<t_gallambflag> _gambflag = nullptr;
		bool _apply_carrier_range = false;
	};
}

#endif // !GTURBOEDIT_H
