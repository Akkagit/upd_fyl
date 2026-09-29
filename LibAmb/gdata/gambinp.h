/**
*
* @verbatim
History
	-1.0 Jie Li  2019-10-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambinp.h
* @brief	Storage the ambinp files' data include par_amb and neq
*
* @author   Jie Li, Wuhan University
* @version	1.0.0
* @date		2019-10-25
*
*/

#ifndef GAMBINP_H
#define GAMBINP_H

#include "gdata/gdata.h"
#include "gmodels/gpar.h"
#include "gutils/ExportAMB.h"

using namespace gnut;

namespace great
{
	class AMB_LIBRARY_EXPORT t_gambinp_head
	{
	public:
		t_gambinp_head();

		~t_gambinp_head();

	public:
		t_gtime beg;
		t_gtime end;
		int		 num_amb;
		double  interval;
		string fix_mode;
		set<string>	sat;
		set<string> site;

	};

	/**
	*@brief	 Class for storaging ambinp file data
	*/
	class AMB_LIBRARY_EXPORT t_gambinp : public t_gdata
	{
	public:

		/** @brief default constructor. */
		t_gambinp();

		/** @brief default destructor. */
		virtual ~t_gambinp();


		void add_hdinfo(t_gambinp_head* hd);
		void add_ambpar(t_gpar ambpar);
		void add_neq(string parname, vector<double>& q);

		void get_hdinfo(t_gambinp_head* hd);
		void get_allambpar(vector<t_gpar>& ambpar);
		void get_neq(vector<string>& parlist, vector<vector<double> >& neq);
	protected:
		vector<t_gpar>			_amb;
		vector<vector<double>>	_neq;
		vector<string>			_par_list;		//ipt for neq(parname)
		t_gambinp_head*			_ambinp_hd = nullptr;
		
	private:

	};

}
#endif // !GAMBINP_H