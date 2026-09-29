/**
*
* @verbatim
History
-1.0 ZhengHJ  2019-09-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		  grecoverdata.h
* @brief	  Storage all recover observ equation/ parameter information
*
* @author     ZhengHJ, Wuhan University
* @version	  1.0.0
* @date		  2019-09-25
*
*/


#ifndef GRECOVERDATA_H
#define GRECOVERDATA_H

#include <set>
#include "gutils/gtime.h"
#include "gmodels/gpar.h"
#include "gutils/ExportBase.h"

using namespace gnut;

namespace great
{
	class BASE_LIBRARY_EXPORT t_grcover_head
	{
	public:
		t_grcover_head();
		~t_grcover_head();

		t_gtime get_beg_time() const;
		t_gtime get_end_time() const;

	public:

		double interval;
		set<string>  sat_list;
		set<string>  site_list;
		set<t_gtime> time_list;
	};

	class BASE_LIBRARY_EXPORT t_grecover_data : public t_gdata
	{

	public:
		t_grecover_data();
		virtual ~t_grecover_data();

	public:
		virtual t_gtime get_recover_time() const = 0;
		virtual string convert2strline() const = 0;
	};

	class BASE_LIBRARY_EXPORT t_grecover_equation : public t_grecover_data
	{
	public:
		t_grecover_equation(const t_gtime& time, const string site, string sat);
		t_grecover_equation(const t_grecover_equation& other);
		~t_grecover_equation();

		void operator =(const t_grecover_equation& other);

		t_gtime get_recover_time() const override;
		string convert2strline() const override;

		void set_recover_equation(const t_gobscombtype& obstype, const pair<double, double>& resinfo, int is_newamb);

	public:
		t_gtime time;
		string sat_name;
		string site_name;
		t_gobscombtype obstype;
		double weight;
		double resuidal;
		int is_newamb;
	};


	class BASE_LIBRARY_EXPORT t_grecover_par : public t_grecover_data
	{
	public:

		t_grecover_par(const t_grecover_par& ohter);
		t_grecover_par(const t_gpar& par, double correct_value);
		~t_grecover_par();

		t_gtime get_recover_time() const override;
		string convert2strline() const override;

		bool operator<(const t_grecover_par&) const;

	public:
		t_gpar par;
		// Correction
		double correct_value;
	};

}

#endif // !GRECOVERDATA_H