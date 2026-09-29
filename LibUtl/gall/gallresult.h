#ifndef GALLRESULT_H
#define GALLRESULT_H

#include <iostream>
#include <string.h>
#include "gutils/ExportUtl.h"
#include "gutils/gmutex.h"
#include "gdata/gdata.h"
#include "gdata/gresult.h"
#include <vector>

using namespace std;
using namespace gnut;

namespace great
{

	class UTL_LIBRARY_EXPORT t_gresult : public t_gdata
	{
	public:
		t_gresult();
		virtual ~t_gresult();


		void setOneResult(int i, t_gresult_one* oneresult);
		void setOnesat(string sat, map<int, t_gresult_one>& mapresult);
		map<int, t_gresult_one>* getMap();
		void clearmap();
		map<string, map<int, t_gresult_one>>* getMapSat();


	protected:
		map<int, t_gresult_one> _map_result;
		map<string, map<int, t_gresult_one>> _map_result_sat;

	};




}

#endif