/**
*
* @verbatim
	History
	 -1.0 zhanglong	 2019-01-06 creat the file.
  @endverbatim
*
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			ExportOI.h
* @brief		for dll on windows
*
* @author		zhanglong
* @version		1.0.0
* @date			2019-01-06
*
*/

#ifndef OI_EXPORT_H
#define OI_EXPORT_H

// for debug
//#include <vld.h>

// 禁用VS警告
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
#  if defined(GREAT_OI_LIBRARY )
#    define OI_LIBRARY_EXPORT  __declspec(dllexport)
#  else
#    define OI_LIBRARY_EXPORT  __declspec(dllimport)
#  endif
#else
#  define OI_LIBRARY_EXPORT
#endif


#ifdef _MSC_VER
# if (_MSC_VER >= 1300)
#  define __STL_MEMBER_TEMPLATES
# endif
#endif


#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
using namespace std;

//#include <vld.h>

#endif 

