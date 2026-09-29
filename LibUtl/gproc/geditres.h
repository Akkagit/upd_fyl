/**
*
* @verbatim
	History
	 -1.0 Hongjie Zheng		2019-09-06  creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    geditres.h
* @brief		Edit residuals function
* @author       Hongjie Zheng, Wuhan University
* @version		1.0.0
* @date		    2019-09-06
*
*/
#ifndef GEDITRES_H
#define GEDITRES_H

#include <gutils/gtime.h>
#include <gio/glog.h>
#include <gutils/ginfolog.h>
#include <gall/gallambflag.h>
#include <gall/gallrecover.h>
#include <gset/gsetbase.h>
#include <gutils/ExportUtl.h>
#include <gall/gallsum.h>
#include <gcoders/sum.h>
#include <gcoders/preproc.h>
#include <gdata/gpreprocslr.h>
#include <gall/gallpreproc.h>

#define TEST "##TEST##"
#define TEST_INFO "##TEST_INFO##"


namespace great
{
	class UTL_LIBRARY_EXPORT t_geditres
	{
	public:
		const static int GOOD;
		const static int OLDAMB;
		const static int NEWAMB;
		const static int LOWELE;
		const static int NOTALL;
		const static int DELSHRT;
		const static int DELBAD;
		const static int DELRMS;
		const static int NEWBAD;
		const static int NODATA;

	public:

		t_geditres(t_glog* log=nullptr);
		t_geditres(t_gsetbase* set,t_glog* log=nullptr);
		virtual ~t_geditres();

		void set_recover_data(t_gallrecover* recover_data);
		void set_recover_data_by_great_resfile(string resfile_name);
		void set_recover_data_by_panda_resfile(string resfile_name);
		void set_preOBS(string filename);

		void set_ambflag_data(t_gallambflag* ambflag_data);
		void set_ambflag_data(const vector<string>& logfiles);

		void set_logfiles(const vector<string>& logfiles);

		t_gallambflag get_gambflag_data() const;

		void processBatch();

		void update_logfiles();

	protected:

		void _rms_table();

		void _slr_res_edit();
		void _check_resid();

		int _find_flag(int istart, int istop, const vector<int>& flg, string flag_description, string direction);

		void  _find_flag_upd(int istart, int istop, const vector<int>& flg, string flag_desription, string direction, int& iepo, int i0, int RLscale, int &j1, int& j2);

		bool _check_jump(vector<int>& flag_data, const vector<double>& val_data);

		bool _remove_short(vector<int>& flag_data);

		void _update_ambflagdata();

		void _generate_pre_obs();

	public:

		// setting
		t_gtime beg_time;
		t_gtime end_time;
		int interval;
		string freq;
		string mode; //L12,L13
		int bad_elisp;  //mm
		int jump_elisp; //mm
		int short_elisp; //seconds
		double res_max;
		double rec_num;

	protected:


		// flag for point if deleted true for deleted false fro not deleted
		// TODO:: Change by shared ptr
		bool _amb_ptrflag;
		bool _recover_ptrflag;

		vector<string> _logfiles;

		// data
		t_gallambflag* _all_ambflag = nullptr;
		t_gallrecover* _all_recover = nullptr;
		t_gallpreproc* _all_preproc = nullptr;
		t_gallsum*     _all_sum = nullptr;

		//_log;
		t_glog* _log = nullptr;

		t_gsetbase* _gset = nullptr;	///< setting from xml file

		// record the iamb 
		map<string, map <string, vector< int > > > _res_ambflag;

		map<int, map<int, int>> _nepotbl;
		map<int, map<int, int>> _nepotbl_slr;
		map<string, map<int, map<int, double>>>        _rmstable_slr;
		map<t_gobscombtype, map<int, map<int,double>>> _rmstable;


		
	public:

		static bool istrue(int flag_code, string flag_type);
	};
}

#endif // !GEDITRES_H
