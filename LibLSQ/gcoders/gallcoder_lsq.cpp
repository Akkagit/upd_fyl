/**
*
* @verbatim
History
	-1.0 Jiande Huang  2019-10-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gallcoder_lsq.cpp
* @brief	main decode function for all files used in podlsq
*
*
* @author   Jiande Huang, Wuhan University
* @version	1.0.0
* @date		2019-10-25
*
*/

#include "gcoders/gallcoder.h"
#include "gcoders/gallcoder_lsq.h"
#include "gio/gio.h"
#include "gio/glog.h"
#include "gio/gfile.h"
namespace great
{

	t_gallcoder_lsq::t_gallcoder_lsq()
	{
	}

	t_gallcoder_lsq::~t_gallcoder_lsq()
	{
	}


	bool t_gallcoder_lsq::support(const IFMT& ifmt)
	{
		if (find(_lsq_support.begin(), _lsq_support.end(), ifmt) == std::end(_lsq_support)) return false;
		return true;
	}

	void t_gallcoder_lsq::coder(const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj)
	{
		if (!t_gallcoder_lsq::support(ifmt)) return;

		t_gcoder* gcoder = nullptr;
		switch (ifmt)
		{
		case IFMT::UNDEF:
		default:
		{
			cout << "not support for this format " << endl;
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
			this->_glog->comment(0, "main", "READ: " + path + " time: "
				+ dbl2str(lstepoch.diff(runepoch)) + " sec");

			if (gio) { delete gio;       gio = nullptr; }
			if (gcoder) { delete gcoder; gcoder = nullptr; }
		}
	}

	void t_gallcoder_lsq::creat(const IFMT& ifmt, map<t_gdata::ID_TYPE, t_gdata*>& map_gdata)
	{
		if (!t_gallcoder_lsq::support(ifmt)) return;

		auto type = ifmt2data(ifmt);
		if (map_gdata.count(type) > 0) return;
		map_gdata[type] = nullptr;
		switch (ifmt)
		{
		case IFMT::UNDEF:
		default:
		{
			map_gdata[type] = nullptr;
			cout << "not support for this format " << endl;
			return;
		}
		}

	}
}
