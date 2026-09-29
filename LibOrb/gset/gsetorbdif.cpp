/**
* @verbatim
	History
	 -1.0 Jie Li	 2019-09-05 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gsetintgra.cpp
* @brief		orbit differ parameter set for satellite orbit differ from XML
*
* @author       Jie Li, Wuhan University
* @version		1.0.0
* @date			2019-09-05
*
*/

#include "gset/gsetorbdif.h"
#include <algorithm>
#include <sstream>

//#define DEBUG
using namespace std;

namespace great {

	t_gsetorbdif::t_gsetorbdif()
	{

	}

	t_gsetorbdif::~t_gsetorbdif()
	{
	}

	string t_gsetorbdif::trans()
	{
		string trans = _doc.child(XMLKEY_ROOT).child(XMLKEY_ORBDIF).child_value(XMLKEY_ORBDIF_TRANS);
		transform(trans.begin(), trans.end(), trans.begin(), ::toupper);
		return trans;
	}

	// added by yqyuan for orbfit
	string t_gsetorbdif::trans_orbfit()
	{
		string trans = _doc.child(XMLKEY_ROOT).child(XMLKEY_ORBFIT).child_value(XMLKEY_ORBDIF_TRANS);
		transform(trans.begin(), trans.end(), trans.begin(), ::toupper);
		return trans;
	}					   
	set<string> t_gsetorbdif::excsys()
	{
		set<string> excsys;
		string sys;
		istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_ORBDIF).child_value(XMLKEY_ORBDIF_EXCSYS));
		while (is >> sys) {
			transform(sys.begin(), sys.end(), sys.begin(), ::toupper);
			excsys.insert(sys);
		}
		return excsys;
	}
	// added by yqyuan
	set<string> t_gsetorbdif::excsys_orbfit()
	{
		set<string> excsys;
		string sys;
		istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_ORBFIT).child_value(XMLKEY_ORBDIF_EXCSYS));
		while (is >> sys) {
			transform(sys.begin(), sys.end(), sys.begin(), ::toupper);
			excsys.insert(sys);
		}
		return excsys;
	}

	set<string> t_gsetorbdif::excsat()
	{
		set<string> excsat;
		string sat;
		istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_ORBDIF).child_value(XMLKEY_ORBDIF_EXCSAT));
		while (is >> sat) {
			transform(sat.begin(), sat.end(), sat.begin(), ::toupper);
			excsat.insert(sat);
		}
		return excsat;
	}
	// added by yqyuan for orbfit
	set<string> t_gsetorbdif::excsat_orbfit()
	{
		set<string> excsat;
		string sat;
		istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_ORBFIT).child_value(XMLKEY_ORBDIF_EXCSAT));
		while (is >> sat) {
			transform(sat.begin(), sat.end(), sat.begin(), ::toupper);
			excsat.insert(sat);
		}
		return excsat;
	}
			  
	void t_gsetorbdif::check()
	{
	}

	void t_gsetorbdif::help()
	{
		cerr << "<orbdif>" << endl
			<< "  <excsat>  < / excsat> " << endl
			<< "  <excsys>  < / excsys> " << endl
			<< "  <trans> STRD < / trans>" << endl
			<< "< / orbdif>" << endl;
	}
	// added by yqyuan
	void t_gsetorbdif::help_orbfit()
	{
		cerr << "<orbfit>" << endl
			<< "  <excsat>  < / excsat> " << endl
			<< "  <excsys>  < / excsys> " << endl
			<< "  <trans> STRD < / trans>" << endl
			<< "< / orbfit>" << endl;
	}			

}

