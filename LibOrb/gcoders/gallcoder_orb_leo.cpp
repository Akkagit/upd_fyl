#include "gcoders/gallcoder.h"
#include "gcoders/gallcoder_orb_leo.h"
#include "gio/gio.h"
#include "gio/glog.h"
#include "gio/gfile.h"
#include "gcoders/icsleo.h"
#include "gcoders/leapsecond.h"
#include "gcoders/satparameters.h"
#include "gcoders/poleut1.h"
#include "gcoders/dvpteph405.h"
#include "gcoders/oceantide.h"
#include "gcoders/egm2008.h"

#include "gcoders/attitude.h"
#include "gcoders/albedo.h"
#include "gcoders/solargeomag.h"
#include "gcoders/pannel.h"

#include "gdata/gsatparam.h"
#include "gdata/gnavde.h"
#include "gall/gallsatparam.h"
#include "gdata/gegmcoef.h"
#include "gdata/gotdata.h"

#include "gall/gallattitude.h"
#include "gall/gallpannel.h"
#include "gdata/gsolargeomag.h"
#include "gdata/galbedo.h"

namespace great
{

	t_gallcoder_orb_leo::t_gallcoder_orb_leo()
	{
		_orb_support.insert(ICSLEO_INP);
		_orb_support.insert(SP3_INP);
		_orb_support.insert(POLEUT1_INP);
		_orb_support.insert(LEAPSECOND_INP);
		_orb_support.insert(EGM_INP);
		_orb_support.insert(OCEANTIDE_INP);
		_orb_support.insert(SATPARS_INP);
		_orb_support.insert(DE_INP);

		_orb_support.insert(ATTITUDE_INP);
		_orb_support.insert(SOLAR_INP);
		_orb_support.insert(PANNEL_INP);
		_orb_support.insert(ALBEDO_INP);

	}

	t_gallcoder_orb_leo::~t_gallcoder_orb_leo()
	{
	}


	bool t_gallcoder_orb_leo::support(const IFMT& ifmt)
	{
		if (find(_orb_support.begin(), _orb_support.end(), ifmt) == std::end(_orb_support)) return false;
		return true;
	}

	void t_gallcoder_orb_leo::coder(const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj)
	{
		if (!t_gallcoder_orb_leo::support(ifmt)) return;

		t_gcoder* gcoder = nullptr;
		switch (ifmt)
		{
		case ICSLEO_INP:        {gcoder = new t_ics_leo(this->_gset, "", 20480);    break; }
		case POLEUT1_INP:    {gcoder = new t_poleut1(this->_gset, "", 4096);    break; }
		case LEAPSECOND_INP: {gcoder = new t_leapsecond(this->_gset, "", 4096); break; }
		case OCEANTIDE_INP:  {gcoder = new t_oceantide(this->_gset, "", 4096);  break; }
		case EGM_INP:        {gcoder = new t_egm2008(this->_gset, "", 4096);    break; }
		case SATPARS_INP:    {gcoder = new t_satparameters(this->_gset, "", 4096); break; }
		case DE_INP: {gcoder = new t_dvpteph405(this->_gset, "", 4096); break; }
		case ATTITUDE_INP: {gcoder = new t_attitude(this->_gset, "", 4096); break; }
		case PANNEL_INP: {gcoder = new t_pannel(this->_gset, "", 4096); break; }
		case SOLAR_INP: {gcoder = new t_solargeomag(this->_gset, "", 4096); break; }
		case ALBEDO_INP: {gcoder = new t_albedo(this->_gset, "", 4096); break; }
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

	void t_gallcoder_orb_leo::creat(const IFMT& ifmt, map<t_gdata::ID_TYPE, t_gdata*>& map_gdata)
	{
		if (!t_gallcoder_orb_leo::support(ifmt)) return;

		auto type = ifmt2data(ifmt);
		if (map_gdata.count(type) > 0) return;
		map_gdata[type] = nullptr;
		switch (ifmt)
		{
		case ICSLEO_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gallics_leo();
			dynamic_cast<t_gallics_leo*>(map_gdata[type])->glog(_glog);
			break;
		}
		case EGM_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gegmcoef();
			dynamic_cast<t_gegmcoef*>(map_gdata[type])->glog(_glog);
			break;
		}
		case OCEANTIDE_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gocean_tide();
			dynamic_cast<t_gocean_tide*>(map_gdata[type])->glog(_glog);
			break;
		}
		case POLEUT1_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gpoleut1();
			dynamic_cast<t_gpoleut1*>(map_gdata[type])->glog(_glog);
			break;
		}
		case LEAPSECOND_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gleapsecond();
			dynamic_cast<t_gleapsecond*>(map_gdata[type])->glog(_glog);
			break;
		}
		case SATPARS_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gallsatinfo();
			dynamic_cast<t_gallsatinfo*>(map_gdata[type])->glog(_glog);
			break;
		}
		case DE_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gnavde();
			dynamic_cast<t_gnavde*>(map_gdata[type])->glog(_glog);
			break;
		}

		case ATTITUDE_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gallatt();
			dynamic_cast<t_gallatt*>(map_gdata[type])->glog(_glog);
			break;
		}
		case PANNEL_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gallpannel();
			dynamic_cast<t_gallpannel*>(map_gdata[type])->glog(_glog);
			break;
		}
		case SOLAR_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gsolargeomag();
			dynamic_cast<t_gsolargeomag*>(map_gdata[type])->glog(_glog);
			break;
		}
		case ALBEDO_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_galbedo();
			dynamic_cast<t_galbedo*>(map_gdata[type])->glog(_glog);
			break;
		}
		
		default:
		{
			map_gdata[type] = nullptr;
			cout << "not support for this format " << endl;
			return;
		}
		}

	}
}
