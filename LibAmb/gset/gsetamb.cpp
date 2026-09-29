/**
*
* @verbatim
History
-1.0 glFeng   2019-02-27 creat the file.
-1.1 BoWong   2019-04-10 Adding Doxygen Style Code Remarks
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		   gsetamb.cpp
* @brief	   control set from XML
*
* @author      glFeng, Wuhan University
* @version	   1.0.0
* @date		   2019-04-10
*
*/

#include "gset/gsetamb.h"
#include <sstream>

using namespace std;
using namespace pugi;

namespace great
{
	/** @brief default constructor. */
	t_gsetamb::t_gsetamb()
	 : t_gsetbase()
	{
		_gmutex.lock();

		_set.insert(XMLKEY_AMBIGUITY);

		//_upd_type = UPD_TYPE::NONE;
		//_fix_mode = FIX_MODE::NO;
		//_ratio = 3.0;                 
		//_addleo = false;
		//_independ_baseline = false;
		//_min_common_time = 900.0;
		//_wl_interval = 30.0;
		//_max_baseline_length = 5000.0;
		//_map_amb_decision[UPD_TYPE::WL]["maxdev"] = 0.15;
		//_map_amb_decision[UPD_TYPE::WL]["maxsig"] = 0.15;
		//_map_amb_decision[UPD_TYPE::WL]["alpha"] = 0.001;
		//_map_amb_decision[UPD_TYPE::NL]["maxdev"] = 0.15;
		//_map_amb_decision[UPD_TYPE::NL]["maxsig"] = 0.15;
		//_map_amb_decision[UPD_TYPE::NL]["alpha"] = 0.001;

		_gmutex.unlock();
	}

	/** @brief default destructor. */
	 t_gsetamb::~t_gsetamb()
	{

	}

	 /**
	 * @brief settings check.
	 * @return	  void
	 */
	 void t_gsetamb::check()
	 {
		 _gmutex.lock();
		 _gmutex.unlock();
		 return;
	 }

	 /**
	 * @brief settings help.
	 * @return	  void
	 */
	 void t_gsetamb::help()
	 {
		 _gmutex.lock();
		 cerr << "<ambiguity>\n"
			 << "<upd> EWL/WL/NL/NONE </upd>\n"
			 << "<fix_mode> ROUND/SEARCH/NO </fix_mode>\n"
			 << "<ratio> 3.0 </ratio>\n"
			 << "<add_leo> NO </add_leo>\n"
			 << "<all_baselines> NO </all_baselines>\n"
			 << "<min_common_time> 30 </min_common_time>\n"
			 << "<baseline_length_limit> 3500 </baseline_length_limit>\n"
			 << "<widelane_interval> 30 </widelane_interval>\n"
			 << "<widelane_decision     maxdev = \"0.15\" maxsig = \"0.10\" alpha = \"1000\"/>\n"
			 << "<narrowlane_decision   maxdev = \"0.15\" maxsig = \"0.10\" alpha = \"1000\"/>\n";

		 _gmutex.unlock();
		 return;
	 }

	 UPDTYPE t_gsetamb::get_updtype()
	 {
		 _gmutex.lock();
		
		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("upd");

		 _gmutex.unlock(); return str2updmode(trim(tmp));
	 }

	 t_gsetamb::FIX_MODE t_gsetamb::fix_mode()
	 {
		 if (!_doc) str2fixmode(string());
		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("fix_mode");
		 return str2fixmode(trim(tmp));
	 }

	 t_gsetamb::UPD_MODE t_gsetamb::upd_mode() {
		 if (!_doc) str2fixmode(string());
		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("upd_mode");
		 return str2upd_mode(trim(tmp));
	 }

	 DD_MODEL t_gsetamb::dd_mode()
	 {
		 if (!_doc) str2fixmode(string());
		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("dd_mode");

		 DD_MODEL mode;
		      if (tmp.find("IF_CB_WN")  != string::npos) mode = DD_MODEL::IF_CB_WN;
		 else if (tmp.find("RAW_CB_WN") != string::npos) mode = DD_MODEL::RAW_CB_WN;
		 else if (tmp.find("RAW_CB")    != string::npos) mode = DD_MODEL::RAW_CB;
		 else if (tmp.find("NONE" )     != string::npos) mode = DD_MODEL::NONE;
		 else
		 {
			 throw logic_error("wrong type of dd mode in node amb.");
		 }
			 
		 return mode;
	 }

	 double t_gsetamb::lambda_ratio()
	 {
		 _gmutex.lock();

		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("ratio");

		 _gmutex.unlock(); return str2dbl(tmp);
	 }

	 double t_gsetamb::bootstrapping()
	 {
		 _gmutex.lock();

		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("boot");

		 if (tmp.empty())
		 {
			 _gmutex.unlock();
			 return -0.001;
		 }

		 _gmutex.unlock(); return str2dbl(tmp);
	 }

	 bool t_gsetamb::addleo()
	 {
		 _gmutex.lock();

		 istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("add_leo"));
		 string tmp;
		 is >> tmp;
		 bool isaddleo = (tmp == "YES" || tmp == "yes");

		 _gmutex.unlock(); return isaddleo;
	 }

	 bool t_gsetamb::independ_baseline()
	 {
		 _gmutex.lock();

		 istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("all_baselines"));
		 string tmp;
		 is >> tmp;
		 bool isindepend_baseline = (tmp == "NO" || tmp == "no");

		 _gmutex.unlock(); return isindepend_baseline;
	 }

	 double t_gsetamb::min_common_time()
	 {
		 _gmutex.lock();

		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("min_common_time");

		 _gmutex.unlock(); return str2dbl(tmp);
	 }

	 double t_gsetamb::max_baseline_length()
	 {
		 _gmutex.lock();

		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("baseline_length_limit");

		 _gmutex.unlock(); return str2dbl(tmp);
	 }

	 double t_gsetamb::wl_interval() 
	 {
		 _gmutex.lock();

		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("widelane_interval");

		 _gmutex.unlock(); return str2dbl(tmp);
	 }

	 map<string, double> t_gsetamb::get_amb_decision(string type)
	 {
		 _gmutex.lock();

		 xml_node tmp_set;
		 map<string, double> amb_decision;
		 if (type == "EWL") {
			 tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child("extra_widelane_decision");
			 string tmp = tmp_set.attribute("maxdev").value();
			 amb_decision["maxdev"] = tmp.empty() ? 0.07 : str2dbl(tmp);
			 tmp = tmp_set.attribute("maxsig").value();
			 amb_decision["maxsig"] = tmp.empty() ? 0.10 : str2dbl(tmp);
			 tmp = tmp_set.attribute("alpha").value();
			 amb_decision["alpha"] = tmp.empty() ? 1000 : str2dbl(tmp);
		 }
		 else if (type == "WL") {
			 tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child("widelane_decision");
			 string tmp = tmp_set.attribute("maxdev").value();
			 amb_decision["maxdev"] = tmp.empty() ? 0.25 : str2dbl(tmp);
			 tmp = tmp_set.attribute("maxsig").value();
			 amb_decision["maxsig"] = tmp.empty() ? 0.10 : str2dbl(tmp);
			 tmp = tmp_set.attribute("alpha").value();
			 amb_decision["alpha"] = tmp.empty() ? 1000 : str2dbl(tmp);
		 }
		 else if (type == "NL") {
			 tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child("narrowlane_decision");
			 string tmp = tmp_set.attribute("maxdev").value();
			 amb_decision["maxdev"] = tmp.empty() ? 0.25 : str2dbl(tmp);
			 tmp = tmp_set.attribute("maxsig").value();
			 amb_decision["maxsig"] = tmp.empty() ? 0.10 : str2dbl(tmp);
			 tmp = tmp_set.attribute("alpha").value();
			 amb_decision["alpha"] = tmp.empty() ? 1000 : str2dbl(tmp);
		 }

		 _gmutex.unlock(); return amb_decision;
	 }

	 bool t_gsetamb::part_ambfix()
	 {
		 _gmutex.lock();

		 istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("part_fix"));
		 string tmp;
		 is >> tmp;
		 bool is_part_fix = (tmp == "YES" || tmp == "yes");

		 _gmutex.unlock(); return is_part_fix;
	 }

	 bool t_gsetamb::carrier_range_out()
	 {
		 _gmutex.lock();

		 istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("carrier_range_out"));
		 string tmp;
		 is >> tmp;
		 bool is_out = (tmp == "YES" || tmp == "yes");

		 _gmutex.unlock(); return is_out;
	 }

	 bool t_gsetamb::apply_carrier_range()
	 {
		 _gmutex.lock();

		 istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("apply_carrier_range"));
		 string tmp = is.str();
		 bool is_apply = (tmp == "YES" || tmp == "yes");

		 _gmutex.unlock(); return is_apply;
	 }

	 int t_gsetamb::part_ambfix_num()
	 {
		 _gmutex.lock();

		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("part_fix_num");

		 if (tmp.empty())
		 {
			 _gmutex.unlock();
			 return 2;
		 }
		 _gmutex.unlock(); return str2dbl(tmp);
	 }

	 //bool t_gsetamb::apply_irc()
	 //{
		// _gmutex.lock();

		// istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("apply_irc"));
		// string tmp;
		// is >> tmp;
		// bool is_apply_irc = (tmp == "YES" || tmp == "yes");

		// _gmutex.unlock(); return is_apply_irc;
	 //}


	 /**
	  * @brief change string to FIXMODE.
	  * @param[in]  str   fix mode string
	  * @return	  FIXMODE
	  */
	 t_gsetamb::FIX_MODE t_gsetamb::str2fixmode(string str)
	 {
		 if (str == "NO") 
		 {
			 return t_gsetamb::FIX_MODE::NO;
		 }
		 else if (str == "ROUND") 
		 {
			 return t_gsetamb::FIX_MODE::ROUND;
		 }
		 else if (str == "SEARCH") 
		 {
			 return t_gsetamb::FIX_MODE::SEARCH;
		 }
		 else if (str == "HOLD")
		 {
			 return t_gsetamb::FIX_MODE::HOLD;
		 }
		 else
		 {
			 cout << "*** warning: not defined ambiguity fixing mode [" << str << "]\n"; cout.flush();
		 }
		 return t_gsetamb::FIX_MODE::NO;
	 }

	 /**
	 * @brief change string to UPD_MODE.
	 * @param[in]  str   upd mode string
	 * @return	  UPD_MODE
	 */
	 t_gsetamb::UPD_MODE t_gsetamb::str2upd_mode(string str) {
		 if (str == "UPD" || str == "upd") {
			 return t_gsetamb::UPD_MODE::UPD;
		 }
		 else if (str == "IRC" || str == "irc") {
			 return t_gsetamb::UPD_MODE::IRC;
		 }
		 else if (str == "OSB" || str == "osb") {
			 return t_gsetamb::UPD_MODE::OSB;
		 }
		 else {
			 cout << "*** warning: not defined upd mode [" << str << "]\n"; cout.flush();
			 return t_gsetamb::UPD_MODE::UPD;
		 }
	 }

	 UPDTYPE t_gsetamb::str2updmode(string str)
	 {
		 str = trim(str);
		 if (str == "wl" || str == "WL")
		 {
			 return UPDTYPE::WL;
		 }
		 else if (str == "ewl" || str == "EWL")
		 {
			 return UPDTYPE::EWL;
		 }
		 else if (str == "ewl24" || str == "EWL24")
		 {
			 return UPDTYPE::EWL24;
		 }
		 else if (str == "ewl25" || str == "EWL25")
		 {
			 return UPDTYPE::EWL25;
		 }
		 else if (str == "ewl_epoch" || str == "EWL_epoch")
		 {
			 return UPDTYPE::EWL_epoch;
		 }
		 else if (str == "nl" || str == "NL")
		 {
			 return UPDTYPE::NL;
		 }
		 else if (str == "ifcb" || str == "IFCB")
		 {
			 return UPDTYPE::IFCB;
		 }
		 else
		 {
			 cout << "*** warning: not defined upd mode [" << str << "]\n"; cout.flush();
		 }
		 return UPDTYPE::NONE;
	 }

	 bool t_gsetamb::isSetRefSat()
	 {
		 _gmutex.lock();

		 istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("set_refsat"));
		 string tmp;
		 is >> tmp;
		 bool isSetRefsat = (tmp == "YES" || tmp == "yes");

		 _gmutex.unlock(); return isSetRefsat;
	 }


}
