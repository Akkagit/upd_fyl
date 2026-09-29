#include "gcoders/gallcoder_utl.h"
#include "gcoders/gcoder.h"

#include "gio/gio.h"
#include "gio/glog.h"
#include "gio/gfile.h"

#include "gcoders/lra.h"
#include "gcoders/npt.h"
#include "gcoders/rb.h"
#include "gcoders/slrf.h"
#include "gcoders/preproc.h"

#include "gdata/gslrf.h"
#include "gall/gallnpt.h"
#include "gall/gallrb.h"
#include "gall/galllra.h"
#include "gall/gallpreproc.h"

using namespace great;
namespace great
{
	t_gallcoder_utl::t_gallcoder_utl():
		_class_id("t_gallcoder_utl"),
		_data_id(0)
	{
		_base_support.insert(RB_INP);
		_base_support.insert(PREOBS_INP);
		_base_support.insert(LRA_INP);
		_base_support.insert(SLRF_INP);
	}

	t_gallcoder_utl::~t_gallcoder_utl()
	{
	}

	bool t_gallcoder_utl::support(const IFMT& ifmt)
	{
		if (find(_base_support.begin(), _base_support.end(), ifmt) == std::end(_base_support)) return false;
		return true;
	}

	void t_gallcoder_utl::coder(const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj)
	{
		if (!t_gallcoder_utl::support(ifmt)) return;

		t_gcoder* gcoder = nullptr;
		switch (ifmt)
		{
		case RB_INP: { gcoder = new t_rb(this->_gset, "", 4096);      break; }
		case PREOBS_INP: { gcoder = new t_preproc(this->_gset, "", 4096);      break; }
		case LRA_INP: { gcoder = new t_lra(this->_gset, "", 4096);         break; }
		case SLRF_INP: { gcoder = new t_slrf(this->_gset, "", 4096);         break; }
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

	void t_gallcoder_utl::creat(const IFMT& ifmt, map<t_gdata::ID_TYPE, t_gdata*>& map_gdata)
	{
		if (!t_gallcoder_utl::support(ifmt)) return;

		auto type = ifmt2data(ifmt);
		if (map_gdata.count(type) > 0) return;

		switch (ifmt)
		{
		case RB_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gallrb();
			dynamic_cast<t_gallrb*>(map_gdata[type])->glog(_glog);
			break;
		}
		case PREOBS_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gallpreproc();
			dynamic_cast<t_gallpreproc*>(map_gdata[type])->glog(_glog);
			break;
		}
		case LRA_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_galllra();
			dynamic_cast<t_galllra*>(map_gdata[type])->glog(_glog);
			break;
		}
		case SLRF_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gslrf();
			dynamic_cast<t_gslrf*>(map_gdata[type])->glog(_glog);
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