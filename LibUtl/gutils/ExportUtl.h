/**
*
* @verbatim
	History
	 -1.0 Huang Jiande	 2019-11-06 creat the file.
  @endverbatim
*
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			ExportUtl.h
* @brief		for dll on windows
*
* @author		Huang Jiande
* @version		1.0.0
* @date			2019-11-06
*
*/

#ifndef UTL_EXPORT_H
#define UTL_EXPORT_H

// ½ûÓÃVS¾¯¸æ
#if defined(_MSC_VER)
#pragma warning( disable : 4244 )
#pragma warning( disable : 4251 )
#pragma warning( disable : 4275 )
#pragma warning( disable : 4512 )
#pragma warning( disable : 4267 )
#pragma warning( disable : 4702 )
#pragma warning( disable : 4511 )
#pragma warning( disable : 4996 )
#endif

#if defined(_MSC_VER) || defined(__CYGWIN__) || defined(__MINGW32__) || defined( __BCPLUSPLUS__)  || defined( __MWERKS__)
#  if defined(GREAT_UTL_LIBRARY )
#    define UTL_LIBRARY_EXPORT  __declspec(dllexport)
#  else
#    define UTL_LIBRARY_EXPORT  __declspec(dllimport)
#  endif
#else
#  define UTL_LIBRARY_EXPORT
#endif


#ifdef _MSC_VER
# if (_MSC_VER >= 1300)
#  define __STL_MEMBER_TEMPLATES
# endif
#endif

// for check whether output the log file or not 
//#define DEBUG_LOG
//#define DEBUG_B_LOG
//#define DEBUG_P_LOG
//#define DEBUG_L_LOG
//#define DEBUG_RES_LOG
//#define DEBUG_LOG
//#define DEBUG_LOG

#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>
#include <exception>
#include <iostream>

// for debug
// #include <vld.h>
using namespace std;

#endif 

