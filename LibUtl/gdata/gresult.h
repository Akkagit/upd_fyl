/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-11-06	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gresult.h
 * @brief		The base class used to store slr caliboration result.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-11-06
 */

#ifndef GRESULT_H
#define GRESULT_H

#include <iostream>
#include <string.h>
#include "gutils/ExportUtl.h"
#include "gutils/gmutex.h"
#include "gdata/gdata.h"
#include <vector>

using namespace std;
using namespace gnut;

namespace great
{

	class UTL_LIBRARY_EXPORT t_gresult_one
	{
	public:

		t_gresult_one();
		virtual ~t_gresult_one();

		void setTime( t_gtime& time);
		void setCode( int code);
		void setName( string name);
		void setNameL( string name);
		void setOMC( double omc);
		void setElev( double elev);
		void setDLRA( double dlra);

		t_gtime _now;
		int _code;
		string _name;
		string _name_long;
		double _omc;
		double _elev;
		double _dlra;

	protected:
		

	};



	

	
}

#endif