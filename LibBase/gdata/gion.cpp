#include "gdata/gion.h"
#include "gcoders/ion.h"

namespace great
{
	great::t_gion::t_gion()
	{
		id_type(t_gdata::ION);
	}

	great::t_gion::~t_gion()
	{
	}

	void great::t_gion::set_header(const t_gtime& beg, const t_gtime& end, const double& intv, const int& num_amb, const set<string>& sats, const set<string>& recs)
	{
		this->_beg = beg;
		this->_end = end;
		this->_intv = intv;
		this->_sats = sats;
		this->_recs = recs;
		this->_numb = num_amb;
	}

	void great::t_gion::get_header(t_gtime& beg, t_gtime& end, double& intv, set<string>& sats, set<string>& recs)
	{
		beg = this->_beg;
		end = this->_end;
		intv = this->_intv;
		sats = this->_sats;
		recs = this->_recs;
	}

	double t_gion::get_value(const string& rec, const string& sat, const t_gtime& beg, const t_gtime& end)
	{
		auto key = make_pair(rec, sat);
		if (_ions.count(key) >0)
		{
			auto ions = _ions[key];

			auto iter_find = _ions[key].lower_bound(beg);
			if (iter_find->second.beg <= beg && iter_find->second.end >= end) return iter_find->second.value;
			else return 0.0;
		}
		return 0.0;
	}

	void great::t_gion::add_record(const t_gion_record& record)
	{
		string rec = record.rec;
		string sat = record.sat;
		t_gtime beg = record.beg;
		t_gtime end = record.end;
		auto key = make_pair(rec, sat);
		this->_ions[key][beg] = record;

		if (this->_beg > beg) this->_beg = beg;
		if (this->_end < end) this->_end = end;
		if (this->_recs.count(rec) == 0) this->_recs.insert(rec);
		if (this->_sats.count(sat) == 0) this->_sats.insert(sat);
		this->_numb++;

	}

	map<t_gtime, t_gion_record>& t_gion::get_records(string rec, string sat)
	{
		// TODO: 在此处插入 return 语句
		auto key = make_pair(rec, sat);
		auto idx = _ions.count(key);
		map<t_gtime, t_gion_record> tmp;
		if (idx > 0) return _ions[key];
		else return tmp;
	}

	t_gion_record::t_gion_record()
	{
	}

	great::t_gion_record::t_gion_record(const string& type, const string& rec, const string& sat, const double& value, const double& sigma, const t_gtime& beg, const t_gtime& end)
	{
		this->ion_type = type; // SION or VION
		this->value = value;
		this->sigma = sigma;
		this->sat = sat;
		this->rec = rec;
		this->beg = beg;
		this->end = end;
	}

	great::t_gion_record::~t_gion_record()
	{
	}

}