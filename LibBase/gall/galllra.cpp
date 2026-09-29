/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-10-27	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		galllra.cpp
 * @brief		The base class used to store lra file information.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-10-27
 */

#include "gall/galllra.h"

using namespace std;

namespace great
{
	/** default constructor */
	t_galllra::t_galllra()
	{
		_mapALL.clear();
		id_type(t_gdata::ALLLRA);
	}

	/** default destructor */
	t_galllra::~t_galllra()
	{

	}
	void lra_leo::clearall()
	{
		_azi1 = 0.0;
		_azi2 = 0.0;
		_dazi = 0.0;
		_dzen = 0.0;
		_zen1 = 0.0;
		_zen2 = 0.0;
		_wavelength = 0.0;

	}

	void  t_galllra::dazi(vector<string>& sats, double d)
	{
		_maplra[sats]._dazi = d;
	}

	double t_galllra::dazi(string sat)
	{
		double dazi = 0.0;
		string leo = getleo(sat);
		for (auto iter = _maplra.begin(); iter != _maplra.end(); iter++)
		{
			auto tab = find((iter->first).begin(), (iter->first).end(), leo);
			if (tab != (iter->first).end())
			{
				dazi = iter->second._dazi;
			}
		}
		return dazi; 
	}

	void  t_galllra::azi1(vector<string>& sats, double d)
	{ 
		_maplra[sats]._azi1 = d;
	}

	double t_galllra::azi1(string sat)
	{
		double azi1;
		string leo = getleo(sat);
		for (auto iter = _maplra.begin(); iter != _maplra.end(); iter++)
		{
			auto tab = find((iter->first).begin(), (iter->first).end(), leo);
			if (tab != (iter->first).end())
			{
				azi1 = iter->second._azi1;
			}
		}
		return azi1; 
	}

	void  t_galllra::azi2(vector<string>& sats, double d)
	{ 
		_maplra[sats]._azi2 = d;
	}

	double t_galllra::azi2(string sat)
	{
		double azi2;
		string leo = getleo(sat);
		for (auto iter = _maplra.begin(); iter != _maplra.end(); iter++)
		{
			auto tab = find((iter->first).begin(), (iter->first).end(), leo);
			if (tab != (iter->first).end())
			{
				azi2 = iter->second._azi2;
			}
		}
		return azi2;
	}

	void  t_galllra::dzen(vector<string>& sats, double d)
	{ 
		_maplra[sats]._dzen = d;
	}          // set/get zenith sampling

	double t_galllra::dzen(string sat)
	{ 
		double dzen = 0.0;
		string leo = getleo(sat);
		for (auto iter = _maplra.begin(); iter != _maplra.end(); iter++)
		{
			auto tab = find((iter->first).begin(), (iter->first).end(), leo);
			if (tab != (iter->first).end())
			{
				dzen = iter->second._dzen;
			}
		}
		return dzen; 
	}

	void  t_galllra::zen1(vector<string>& sats, double d)
	{ 
		_maplra[sats]._zen1 = d;
	}          // set/get zenith1

	double t_galllra::zen1(string sat)
	{ 
		double zen1;
		string leo = getleo(sat);
		for (auto iter = _maplra.begin(); iter != _maplra.end(); iter++)
		{
			auto tab = find((iter->first).begin(), (iter->first).end(), leo);
			if (tab != (iter->first).end())
			{
				zen1 = iter->second._zen1;
			}
		}
		return zen1; 
	}

	void  t_galllra::zen2(vector<string>& sats, double d)
	{ 
		_maplra[sats]._zen2 = d;
	}          // set/get zenith2

	double t_galllra::zen2(string sat)
	{ 
		double zen2;
		string leo = getleo(sat);
		for (auto iter = _maplra.begin(); iter != _maplra.end(); iter++)
		{
			auto tab = find((iter->first).begin(), (iter->first).end(), leo);
			if (tab != (iter->first).end())
			{
				zen2 = iter->second._zen2;
			}
		}
		return zen2; 
	}

	void t_galllra::wave(vector<string>& sats, double d)
	{
		_maplra[sats]._wavelength = d;
	}

	double t_galllra::wave(string sat)
	{
		double wave;
		string leo = getleo(sat);
		for (auto iter = _maplra.begin(); iter != _maplra.end(); iter++)
		{
			auto tab = find((iter->first).begin(), (iter->first).end(), leo);
			if (tab != (iter->first).end())
			{
				wave = iter->second._wavelength;
			}
		}
		return wave;
	}

	void t_galllra::addlra(vector<string>& sats,map<int, t_map_Z>& mapA)
	{
		_mapALL[sats] = mapA;
	}

	int t_galllra::lrasize()
	{
		return _mapALL.size();
	}

	double t_galllra::getlra(string sat, int zen, int azi)
	{
		t_map_Z mapZ;
		double cor;
		t_map_A _mapA;
		string leo;
		leo = getleo(sat);
		for (auto iter = _mapALL.begin(); iter != _mapALL.end(); iter++)
		{
			auto tab = find((iter->first).begin(), (iter->first).end(), leo);
			if (tab != (iter->first).end())
			{
				_mapA = iter->second;
			}
		}
		mapZ = _mapA[azi];
		cor = mapZ[zen];
		return cor;
	}

	string t_galllra::getleo(string sat)
	{
		string leo;
		if (sat.find("GRACE") != string::npos)
		{
			leo = "GRACE";
		}
		else if (sat.find("SWARM") != string::npos)
		{
			leo = "SWARM";
		}
		else if (sat.find("CHAMP") != string::npos)
		{
			leo = "CHAMP";
		}
		else if (sat.find("GOCE") != string::npos)
		{
			leo = "GOCE";
		}
		else if (sat.find("SENTINEL") != string::npos)
		{
			leo = "SENTINEL";
		}
		else
		{
			leo = sat;
		}
		return leo;
	}

}