
#ifndef GPRODICS_H
#define GPRODICS_H 

#include "gprod/gprod.h"
#include "gutils/gtriple.h"
#include "gforce/gfcommon.h"
#include "gutils/ExportOI.h"
using namespace std;
using namespace gnut;
namespace great
{
	class OI_LIBRARY_EXPORT t_gprodics : public t_gprod
	{

	public:
		struct sat_record
		{
			string par_type;
			double value;
			double rms;
			int    beg_mjd_day;
			double beg_mjd_sod;

			int    end_mjd_day;
			double end_mjd_sod;
		};

		t_gprodics(const t_gtime& t, shared_ptr<t_gobj> pt = nullobj);
		virtual ~t_gprodics();

		void set_sat(const string& sat);
		void set_sys(const string& sys);
		void set_record(const string& type, const double& value, const double& rms, const t_gtime& beg, const t_gtime& end);

		double get_rms(const string& type, const string& sat, const t_gtime& beg, const t_gtime& end);
		double get_value(const string& type, const string& sat, const t_gtime& beg, const t_gtime& end);

	protected:
		string  _sat;
		string  _sys;
		map<ICS_PAR, sat_record> _value;
	private:

	};

} // namespace

#endif
