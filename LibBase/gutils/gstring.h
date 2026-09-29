/**
* Some tools for string
* @verbatim
	History
	 -1.0 jdhuang	 2019-04-07 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gstring.h
* @brief		for string
*
* @author       jdhuang, Wuhan University
* @version		1.0.0
* @date			2019-04-07
*
*/

#include <string>
#include <vector>
#include <set>
#include "gutils/ExportBase.h"
using namespace std;

namespace great
{
	BASE_LIBRARY_EXPORT std::string trim(std::string s);
	BASE_LIBRARY_EXPORT std::string clean(std::string s);
	BASE_LIBRARY_EXPORT void split(const std::string& s, std::string delim, std::vector< std::string >& ret);
	
	BASE_LIBRARY_EXPORT std::string set2str(std::set<std::string> &s);
	BASE_LIBRARY_EXPORT std::string str2upper(const std::string& str);
	BASE_LIBRARY_EXPORT std::string str2lower(const std::string& str);
	BASE_LIBRARY_EXPORT std::string do_fraction(long double value, int digits);
	BASE_LIBRARY_EXPORT std::string format(const char* fmt, ...);

	// 去除多余空格
	BASE_LIBRARY_EXPORT std::string remove_surplus_spaces(const std::string& src);

}
