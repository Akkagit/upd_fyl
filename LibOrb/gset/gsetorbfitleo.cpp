

#include "gset/gsetorbfitleo.h"
#include "gforce/gfcommon.h"
#include <algorithm>
#include <sstream>

//#define DEBUG
using namespace std;

namespace great {

	t_gsetorbfitleo::t_gsetorbfitleo()
	{

	}

	t_gsetorbfitleo::~t_gsetorbfitleo()
	{
	}

	string t_gsetorbfitleo::trans()
	{
		string trans = _doc.child(XMLKEY_ROOT).child(XMLKEY_ORBDIF).child_value(XMLKEY_ORBDIF_TRANS);
		transform(trans.begin(), trans.end(), trans.begin(), ::toupper);
		if (!trans.empty())
		{
			trans.erase(0, trans.find_first_not_of(" "));
			trans.erase(trans.find_last_not_of(" ") + 1);
		}
		return trans;
	}

	

	set<string> t_gsetorbfitleo::leo()
	{
		set<string> excsat;
		string sat;
		istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_ORBDIF).child_value(XMLKEY_ORBDIF_LEO));
		while (is >> sat) {
			transform(sat.begin(), sat.end(), sat.begin(), ::toupper);
			excsat.insert(sat);
		}
		return excsat;
	}

	double t_gsetorbfitleo::unit()
	{
		istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_ORBDIF).child_value(XMLKEY_ORBDIF_UNIT));
		string unitstr;
		double unit = 1.0;
		while (is >> unitstr) {
			transform(unitstr.begin(), unitstr.end(), unitstr.begin(), ::tolower);
			if (unitstr == "mm" || unitstr == "")
			{
				unit = 1.0;
			}
			else if (unitstr == "cm")
			{
				unit = 1E-1;
			}
			else if (unitstr == "dm")
			{
				unit = 1E-2;
			}
			else if (unitstr == "m")
			{
				unit = 1E-3;
			}
			else
			{
				cout << "orbfitleo:no unit like this: " << unitstr << endl;
			}
		}
		
		
		return unit;
	}
	bool t_gsetorbfitleo::isfit()
	{
		string fit;
		bool isfit;
		istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_ORBDIF).child_value(XMLKEY_ORBDIF_FIT));
		while (is >> fit) {
			isfit = str2bool(fit);
		}
		return isfit;
	}

	bool t_gsetorbfitleo::isinter()
	{
		string fit;
		bool isinter = true;;
		istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_ORBDIF).child_value(XMLKEY_ORBDIF_INTER));
		while (is >> fit) {
			isinter = str2bool(fit);
		}
		return isinter;
	}


	void t_gsetorbfitleo::check()
	{
	}

	void t_gsetorbfitleo::help()
	{
		cerr << "<orbdifleo>" << endl
			<< "  <leo>  < / leo> " << endl
			<< "  <trans> STRD < / trans>" << endl
			<< "< / orbdifleo>" << endl;
	}

}

