/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-10-27	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		galllra.h
 * @brief		The base class used to store lra file information.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-10-27
 */

#ifndef GALLLRA_H
#define GALLLRA_H

#include <iostream>
#include <string.h>
#include "gdata/gdata.h"

using namespace std;
using namespace gnut;

namespace great
{
	class BASE_LIBRARY_EXPORT lra_leo
	{
	public:
		double         _azi1;
		double         _azi2;
		double         _dazi;       // azimuth sampling
		double         _dzen;       // zenith sampling
		double         _zen1;       // zenith1
		double         _zen2;       // zenith2
		double         _wavelength;
		void clearall();
	};

	class BASE_LIBRARY_EXPORT t_galllra : public t_gdata
	{
	public:
		t_galllra();
		virtual ~t_galllra();

		typedef map<int, double>    t_map_Z;      // map of ZEN-dependence
		typedef map<int, t_map_Z>   t_map_A;      // map of AZI-dependence & ZEN-dependence
		typedef map<vector<string>, t_map_A>  t_map_ALL;      
		typedef map<vector<string>, lra_leo>   t_map_head;

		void   dazi(vector<string>& sats, double d);          // set/get azimuth sampling
		double dazi(string sat);
		void   azi1(vector<string>& sats, double d);        // set/get zenith1
		double azi1(string sat);
		void   azi2(vector<string>& sats, double d);          // set/get zenith2
		double azi2(string sat);
		void   dzen(vector<string>& sats, double d);          // set/get zenith sampling
		double dzen(string sat);
		void   zen1(vector<string>& sats, double d);       // set/get zenith1
		double zen1(string sat);
		void   zen2(vector<string>& sats, double d);         // set/get zenith2
		double zen2(string sat);
		void   wave(vector<string>& sats, double d);
		double wave(string sat);

		void addlra(vector<string>& sats,map<int, t_map_Z>& mapA);

		double getlra(string sat, int zen, int azi);
		int lrasize();
		string getleo(string sat);

	private:
		t_map_head     _maplra;
		//t_map_A        _mapA;
		t_map_ALL      _mapALL;
	};
}

#endif