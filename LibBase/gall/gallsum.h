/**
*
* @verbatim
	History
	 -1.0 Hongmin Zhang 2020-01-14 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gallsum.h
* @brief		    Storage the sum of residuals.
*
* @author       Hongmin Zhang, Wuhan University
* @version		1.0.0
* @date		    2020-01-14
*
*/

#ifndef GALLSUM_H
#define GALLSUM_H

#include "gdata/gdata.h"
#include "gset/gsetproc.h"

using namespace std;
using namespace gnut;

namespace great
{
	/**
	*@brief	   Class for all LEO ICS data storaging
	*/
	class BASE_LIBRARY_EXPORT t_gallsum : public t_gdata
	{
	public:
		t_gallsum();
		~t_gallsum();

		void setSat(vector<string>& satlist);
		void setSite(vector<string>& sitelist);
		void setRMS(map<t_gobscombtype, map<int, map<int, double>>>& rms);
		void setNEPO(map<int, map<int, int>>& nepo);
		void setNEPO_SLR(map<int, map<int, int>>& nepo);
		void setRMS_SLR(map<string, map<int, map<int, double>>>& rms);

		vector<string>& getSat();
		vector<string>& getSite();
		map<t_gobscombtype, map<int, map<int, double>>>& getRMS();
		map<int, map<int, int>>& getNEPO();
		map<string, map<int, map<int, double>>>& getRMS_SLR();
		map<int, map<int, int>>& getNEPO_SLR();

	protected:
		vector<string> _sat_list;
		vector<string> _site_list;
		map<t_gobscombtype, map<int, map<int, double>>> _rmstable;
		map<string, map<int, map<int, double>>> _rmstable_slr;
		map<int, map<int, int>> _nepotbl;
		map<int, map<int, int>> _nepotbl_slr;
	};
}


#endif