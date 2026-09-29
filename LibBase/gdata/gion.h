/**
*
* @verbatim
History
	-1.0 Jiande huang  2019-11-14  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gion.h
* @brief	Storage the ion files' data
*
* @author   Jiande huang, Wuhan University
* @version	1.0.0
* @date		2019-11-14
*
*/

#ifndef GION_H
#define GION_H

#include "gutils/ExportBase.h"
#include "gdata/gdata.h"

using namespace gnut;

namespace great
{
	class BASE_LIBRARY_EXPORT t_gion_record
	{
	public:
		t_gion_record();
		t_gion_record(const string& type, const string& rec, const string& sat,
			const double& value, const double& sigma,
			const t_gtime& beg, const t_gtime& end);
		virtual ~t_gion_record();

		string ion_type; // SION or VION
		double value;
		double sigma;
		string sat;
		string rec;
		t_gtime beg;
		t_gtime end;
	};

	/**
	*@brief	 Class for storaging ion file data
	*/
	class BASE_LIBRARY_EXPORT t_gion : public t_gdata
	{
	public:

		/** @brief default constructor. */
		t_gion();

		/** @brief default destructor. */
		virtual ~t_gion();

		void    set_header(const t_gtime& beg, const t_gtime& end, const double& intv, const int& num_amb, const set<string>& sats, const set<string>& recs);
		void    get_header(t_gtime& beg,t_gtime& end, double& intv, set<string>& sats, set<string>& recs);
		double  get_value(const string& rec, const string& sat, const t_gtime& beg, const t_gtime& end);
		void   add_record(const t_gion_record& record);
		map<t_gtime, t_gion_record>& get_records(string rec, string sat);
	protected:
		map< pair<string, string>, map<t_gtime, t_gion_record> > _ions; // rec_sat_beg
	private:
		t_gtime		_beg;
		t_gtime     _end;
		double		_intv = 300.0;
		int			_numb = 0;
		set<string> _sats;
		set<string> _recs;
	};
}

#endif // !GION_H
