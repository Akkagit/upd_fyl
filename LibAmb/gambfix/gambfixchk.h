

#ifndef GAMBCHK_H
#define GAMBCHK_H 

#include "gutils/ExportAMB.h"
#include "gambfix/gambfix.h"
#include "gset/gsetamb.h"
#include "gdata/gambcon.h"

using namespace std;

namespace great
{

	/**
	*@ brief class for LAMADA search finding fixing ambiguity.
	*/
	class AMB_LIBRARY_EXPORT t_gambchk
	{
	public:
		t_gambchk();
		t_gambchk(t_glog* log, t_gsetbase* set);
		
		virtual ~t_gambchk();

		
		// 确定beg到end时间内的单差模糊度
		bool ProcessBatch();
		bool add_data(t_gallrecover* res1, t_gallrecover* res2, t_gambcon* ambcon);


	protected:
		bool get_setting();
		void write_ambcon();
		
	private:
		string         _class_id = "t_gambchk";
		double         _intv;
		t_gtime        _beg;
		t_gtime        _end;
		double         _bias;
		t_glog*        _glog = nullptr;	///< log file for process
		t_gsetbase*    _gset = nullptr;	///< setting from xml file
		t_gallrecover* _grecover1 = nullptr;
		t_gallrecover* _grecover2 = nullptr;
		/*t_map_time_equ& _res1;
		t_map_time_equ& _res2;*/
		map<string, map<string, map<int,int>>>_flag;
		t_gambcon*     _gambcon = nullptr;
		int            _ambnum;
		set<string>    _sites;
		set<string>    _sats;

	};
}

#endif