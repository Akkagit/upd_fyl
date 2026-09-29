/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2020-03-06	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gallsta.h
 * @brief		The base class used to store range biase information.
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2020-03-06
 */

#ifndef GALLSTA_H
#define GALLSTA_H

#include <iostream>
#include <string.h>
#include "gdata/gdata.h"

using namespace std;
using namespace gnut;

namespace great
{
	class BASE_LIBRARY_EXPORT sta_cor
	{
	public:
		
		double dx;
		double dy;
		double dz;
	};

	class BASE_LIBRARY_EXPORT t_gallsta : public t_gdata
	{
	public:
		t_gallsta();
		virtual ~t_gallsta();

		void addsta(string name, double dx,double dy,double dz);

		map<string, sta_cor>& getsta();

	private:
		map<string, sta_cor> _Sta_Cor;

	};

}

#endif