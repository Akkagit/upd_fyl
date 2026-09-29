/**
*
* @verbatim
	History
	  -1.0 jdhuang 2019-05-17  creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gsetlsq.h
* @brief	The setting of lsq estimator
*
* @author   glfeng BoWong xjhan, Wuhan University
* @version	1.0.0
* @date		2019-04-08
*
*/

#ifndef GSETFLT_H
#define GSETFLT_H

#define XMLKEY_FLT "filter"

#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gutils/gtypeconv.h"
#include "gset/gsetbase.h"
#include "gutils/ExportBase.h"
using namespace std;

namespace gnut
{
	/// The class of settings for lsq estimator
	class BASE_LIBRARY_EXPORT t_gsetlsq : public virtual t_gsetbase
	{
	public:
		t_gsetlsq();
		virtual ~t_gsetlsq();

		void check();                                  // settings check
		void help();                                   // settings help

	};

} // namespace

#endif
