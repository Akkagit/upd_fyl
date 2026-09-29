/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gallcoder_amb.cpp
* @brief	function for ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#include "gcoders/gallcoder_amb.h"
#include "gcoders/gcoder.h"
#include "gcoders/ambcon.h"
#include "gcoders/ambinp.h"
#include "gcoders/ambupd.h"
#include "gcoders/upd.h"

#include "gio/gio.h"
#include "gio/glog.h"
#include "gio/gfile.h"

#include "gall/gallsatparam.h"
#include "gall/gallprec.h"
#include "gall/gallorb.h"
#include "gdata/gleapsecond.h"

using namespace gnut;
namespace great
{
	t_gallcoder_amb::t_gallcoder_amb()
	{
		_base_support.insert(AMBCON_INP);
		_base_support.insert(AMBUPD_INP);
		_base_support.insert(AMBINP_INP);
		_base_support.insert(UPD_INP);
	}

	t_gallcoder_amb::~t_gallcoder_amb()
	{
	}

	bool t_gallcoder_amb::support(const IFMT& ifmt)
	{
		if (find(_base_support.begin(), _base_support.end(), ifmt) == std::end(_base_support)) return false;
		return true;
	}

	void t_gallcoder_amb::coder(const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj)
	{
		if (!t_gallcoder_amb::support(ifmt)) return;

		t_gcoder* gcoder = nullptr;
		switch (ifmt)
		{
		case AMBCON_INP: { gcoder = new t_ambcon(this->_gset, "", 10240);         break; }
		case AMBINP_INP: { gcoder = new t_ambinp(this->_gset, "", 4096);         break; }
		case AMBUPD_INP: { gcoder = new t_ambupd(this->_gset, "", 4096);		 break; }
		case UPD_INP: { gcoder = new t_upd(this->_gset, "", 4096);      break; }
		default:
		{
			if (_glog)
			{
				_glog->comment(t_glog::LOG_LV::LOG_ERROR, _class_id, "coder", "not support for this format " + t_gsetinp::ifmt2str(ifmt));
			}
			return;
		}
		}

		if (gcoder)
		{
			// prepare gio for the file
			t_gio* gio = new t_gfile;
			gio->glog(this->_glog);
			gio->path(path);

			// Put the file into gcoder
			gcoder->clear();
			gcoder->path(path);
			gcoder->glog(this->_glog);

			// Put the data container into gcoder
			_data_id++;
			gcoder->add_data("ID"+int2str(_data_id), data);
			gcoder->add_data("OBJ", obj);

			// Put the gcoder into the gio
			// Note, gcoder contain the gdata and gio contain the gcoder
			gio->coder(gcoder);

			// Read the data from file here
			t_gtime runepoch = t_gtime::current_time(t_gtime::GPS);
			gio->run_read();
			t_gtime lstepoch = t_gtime::current_time(t_gtime::GPS);
			this->_glog->comment(t_glog::LOG_LV::LOG_INFO, "coder", "main", "READ: " + path + " time: "
				+ dbl2str(lstepoch.diff(runepoch)) + " sec");
			if (gio) { delete gio;       gio = nullptr; }
			if (gcoder) { delete gcoder; gcoder = nullptr; }
		}
	}

	void t_gallcoder_amb::creat(const IFMT& ifmt, map<t_gdata::ID_TYPE, t_gdata*>& map_gdata)
	{
		if (!t_gallcoder_amb::support(ifmt)) return;

		auto type = ifmt2data(ifmt);
		if (map_gdata.count(type) > 0) return;

		switch (ifmt)
		{
		case AMBCON_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gambcon();
			dynamic_cast<t_gambcon*>(map_gdata[type])->glog(_glog);
			break;
		}
		case AMBINP_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gambinp();
			dynamic_cast<t_gambinp*>(map_gdata[type])->glog(_glog);
			break;
		}
		case AMBUPD_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gallambupd();
			dynamic_cast<t_gallambupd*>(map_gdata[type])->glog(_glog);
			break;
		}
		case UPD_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gupd();
			dynamic_cast<t_gupd*>(map_gdata[type])->glog(_glog);
			break;
		}
		default:
		{
			map_gdata[type] = nullptr;
			if (_glog)
			{
				_glog->comment(t_glog::LOG_LV::LOG_ERROR, _class_id, "creat", "not support for this format : " + t_gsetinp::ifmt2str(ifmt));
			}
			return;
		}
		}
	}
}