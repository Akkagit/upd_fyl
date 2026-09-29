#include "gprod/gprodics.h"
namespace great
{
	t_gprodics::t_gprodics(const t_gtime& t, shared_ptr<t_gobj> pt)
		: t_gprod(t,pt)
	{
	}

	t_gprodics::~t_gprodics()
	{
	}

	void t_gprodics::set_sat(const string& sat)
	{
		_sat = sat;
	}

	void t_gprodics::set_sys(const string& sys)
	{
		_sys = sys;
	}

	void t_gprodics::set_record(const string& type, const double& value, const double& rms, const t_gtime& beg, const t_gtime& end)
	{
		sat_record record;
		record.par_type = type;
		record.value = value;
		record.rms = rms;
		record.beg_mjd_day = (int)(beg.dmjd());
		record.beg_mjd_sod = (beg.dmjd() - record.beg_mjd_day) * 86400.0;
		record.end_mjd_day = (int)(end.dmjd());
		record.end_mjd_sod = (end.dmjd() - record.end_mjd_day) * 86400.0;

		_value.insert(pair<ICS_PAR, sat_record>(str2ics(type), record));
	}

	double t_gprodics::get_rms(const string& type, const string& sat, const t_gtime& beg, const t_gtime& end)
	{
		return 0.0;
	}

	double t_gprodics::get_value(const string& type, const string& sat, const t_gtime& beg, const t_gtime& end)
	{
		return 0.0;
	}
}


