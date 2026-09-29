
#ifndef GSLRF_H
#define GSLRF_H

#include <iostream>
#include <string.h>
#include "gdata/gdata.h"
#include <vector>

using namespace std;
using namespace gnut;

namespace great
{
	class t_gslrf_one
	{
	public:
		int _code;
		string _DOmes;
		string _name_gps;
		double _xsta[6];
		double _blh[3];
		t_gtime now;
		double _offset[3];

	};

	class BASE_LIBRARY_EXPORT  t_gslrf : public t_gdata
	{
	public:
		t_gslrf();
		virtual ~t_gslrf();

		t_gslrf_one& getOnegslrf(int i);

		void setCode(int i,int code);

		void setDomes(int i, string domes);

		void setNameGPS(int i, string name);

		void setXsta(int i, int num, double xv);

		int getCode(int i);

		void setGSLRF(t_gslrf& gslrf);

		void setREF(int i,t_gtime now);

		map<int, t_gslrf_one>& getgslrf();

		int getSize();
		
		
	protected:
		map<int, t_gslrf_one> _gslrf;

	};
}

#endif