/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-01-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		ExportBase.h
* @brief	for the use of Lib Base.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*
*/

#ifndef BASE_EXPORT_H
#define BASE_EXPORT_H


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
#  if defined(GREAT_BASE_LIBRARY )
#    define BASE_LIBRARY_EXPORT  __declspec(dllexport)
#  else
#    define BASE_LIBRARY_EXPORT  __declspec(dllimport)
#  endif
#else
#  define BASE_LIBRARY_EXPORT
#endif


#ifdef _MSC_VER
# if (_MSC_VER >= 1300)
#  define __STL_MEMBER_TEMPLATES
# endif
#endif



#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <string.h>
#include <stdio.h>
#include <iostream>

using namespace std;


#ifndef SIZE_INT
#define SIZE_INT sizeof(int)
#endif

#ifndef SIZE_DBL
#define SIZE_DBL sizeof(double)
#endif // !SIZE_DBL


//#include "vld.h"

#endif //BASE_EXPORT_H

