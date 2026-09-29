/**
*
* @verbatim
	History
	 -1.0 GREAT	    2019-01-04 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			  gcycleslip.h
* @brief		  CRS coordinate transform to ACR coordinate
*
* @author         GREAT, Wuhan University
* @version		  1.0.0
* @date			  2019-01-04
*
*/

#ifndef GINFOLOG_H
#define GINFOLOG_H

#include "gio/glog.h"
#include "gutils/ExportBase.h"
/*
	ERROR 0 : error information
	NOTE  0 : important note information
	WARN  0 : important warnning information
	WARN  1 : warnning information that not important
	NOTE  2 :
*/

#define LV1 "ERROR  "
#define LV2 "WANNING"
#define LV3 "NOTE   "
#define LV4 "OUTPUT "

namespace gnut
{
	enum class LOG_LV { LOG_ERROR, LOG_WARN, LOG_NOTE, LOG_OUTPUT, LOG_DEF};

	void   BASE_LIBRARY_EXPORT write_log(t_glog* log, const LOG_LV& LV, const string& class_id, const string& function_id, const string& description);

	void   BASE_LIBRARY_EXPORT write_log(t_glog* log, int l, string description);
	void   BASE_LIBRARY_EXPORT write_log(t_glog* log, int l, string ID, string description);

	void   BASE_LIBRARY_EXPORT writeLogInfo(t_glog* log, int l, string description);
	void   BASE_LIBRARY_EXPORT writeLogInfo(t_glog* log, int l, string ID, string description);

	void   BASE_LIBRARY_EXPORT write_log_info(t_glog* log, int l, string description);
	void   BASE_LIBRARY_EXPORT write_log_info(t_glog* log, int l, string ID, string description);



}
#endif // !GINFOLOG_H
