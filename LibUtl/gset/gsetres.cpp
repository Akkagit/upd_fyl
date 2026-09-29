/*
* @verbatim
	================================================================<!-->
	<xml>
		 <editres>
			<freq>        LC  </freq>		LC/L1/L2/L3/L4/L5
			<short_elisp> 600 </short_elisp>
			<jump_elisp>  40  </jump_elisp> 80/40
			<bad_elisp>   40  </bad_elisp>
		 </editres>

	</xml>
 @endverbatim

* @verbatim
	History
	 -1.0 Hongjie Zheng		2019-09-06  creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gsetres.cpp
* @brief		the xml set for edit residual
*
* @author       Hongjie Zheng, Wuhan University
* @version		1.0.0
* @date			2019-09-06
*
*/

#include <sstream>
#include "gset/gsetres.h"

using namespace std;
using namespace pugi;

namespace great
{

	 t_gsetres::t_gsetres():
		 _default_freq("LC"),
		 _default_mode("L12"),
		 _default_jump_elisp(80),
		_default_bad_elisp(80),
		_default_short_elisp(600),
		 _default_res_max(0.2),
		 _default_sta_num(10)
	 {
		 _set.insert(XMLKEY_RES);
	 }
	 t_gsetres::~t_gsetres()
	 {
	 }
	 void t_gsetres::check()
	 {
		 _gmutex.lock();
		 // check existence of nodes/attributes
		 xml_node parent = _doc.child(XMLKEY_ROOT);
		 xml_node node = _default_node(parent, XMLKEY_RES);

		 _gmutex.unlock();
		 return;
	 }
	 void t_gsetres::help()
	 {
		 _gmutex.lock();
		 cerr << "<editres>\n"
			 << "<short_elisp> 600 </short_elisp> \n"
			 << "<jump_elisp>  50  </jump_elisp> \n"
			 << "<bad_elisp>   50 </bad_elisp> \n";

		 _gmutex.unlock();
		 return;
	 }
	 double t_gsetres::get_jump_elisp()
	 {
		 _gmutex.lock();
		 //xml_node parent = _doc.child(XMLKEY_ROOT).child(XMLKEY_RES);
		 //xml_node node = _default_node(parent, "jump_elisp");
		 //string tmp = node.value();
		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_RES).child_value("jump_elisp");
		 _gmutex.unlock(); 

		 if (tmp == "") {
			 return _default_jump_elisp;
		 }
		 else {
			 return str2dbl(tmp);
		 }
	 }
	 string t_gsetres::get_mode()
	 {
		 _gmutex.lock();
		 //xml_node parent = _doc.child(XMLKEY_ROOT).child(XMLKEY_RES);
		 //xml_node node = _default_node(parent, "jump_elisp");
		 //string tmp = node.value();
		 string tmp = _doc.child(XMLKEY_ROOT).child_value("mode");
		 _gmutex.unlock();

		 if (tmp == "") {
			 return "UNDEF";
		 }
		 else {
			 return tmp;
		 }
	 }
	 double t_gsetres::get_sta_num()
	 {
		 _gmutex.lock();
		 //xml_node parent = _doc.child(XMLKEY_ROOT).child(XMLKEY_RES);
		 //xml_node node = _default_node(parent, "bad_elisp");
		 //string tmp = node.value();
		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_RES).child_value("sta_num");
		 _gmutex.unlock();

		 if (tmp == "") {
			 return _default_sta_num;
		 }
		 else {
			 return str2dbl(tmp);
		 }
	 }
	 double t_gsetres::get_res_max()
	 {
		 _gmutex.lock();
		 //xml_node parent = _doc.child(XMLKEY_ROOT).child(XMLKEY_RES);
		 //xml_node node = _default_node(parent, "bad_elisp");
		 //string tmp = node.value();
		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_RES).child_value("res_max");
		 _gmutex.unlock();

		 if (tmp == "") {
			 return _default_res_max;
		 }
		 else {
			 return str2dbl(tmp);
		 }
	 }
	 double t_gsetres::get_bad_elisp()
	 {
		 _gmutex.lock();
		 //xml_node parent = _doc.child(XMLKEY_ROOT).child(XMLKEY_RES);
		 //xml_node node = _default_node(parent, "bad_elisp");
		 //string tmp = node.value();
		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_RES).child_value("bad_elisp");
		 _gmutex.unlock();

		 if (tmp == "") {
			 return _default_bad_elisp;
		 }
		 else {
			 return str2dbl(tmp);
		 }
	 }
	 double t_gsetres::get_short_elisp()
	 {
		 _gmutex.lock();
		 //xml_node parent = _doc.child(XMLKEY_ROOT).child(XMLKEY_RES);
		 //xml_node node = _default_node(parent, "short_elisp");
		 //string tmp = node.value();
		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_RES).child_value("short_elisp");
		 _gmutex.unlock();

		 if (tmp == "") {
			 return _default_short_elisp;
		 }
		 else {
			 return str2dbl(tmp);
		 }
	 }

	 string t_gsetres::get_freq()
	 {
		 _gmutex.lock();
		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_RES).child_value("freq");
		 _gmutex.unlock();

		 if (tmp.empty()) {
			 return _default_freq;
		 }
		 else {
			 return trim(tmp);
		 }
	 }

	 string t_gsetres::get_freq_mode()
	 {
		 _gmutex.lock();
		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_RES).child_value("mode");
		 _gmutex.unlock();

		 if (tmp.empty())
		 {
			 return _default_mode;
		 }
		 else 
		 {
			 return trim(tmp);
		 }
	 }
}
