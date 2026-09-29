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

#include "ginfolog.h"
#include <iostream>
#include <iomanip>

void BASE_LIBRARY_EXPORT gnut::write_log(t_glog* log, const LOG_LV& LV, const string& class_id, const string& function_id, const string& description)
{

}

void gnut::write_log(t_glog * log, int l, string description)
{
	write_log_info(log, l, description);
}

void gnut::write_log(t_glog * log, int l, string ID, string description)
{
	if (log)
	{
		log->comment(l, ID, description);
	}
	else
	{
		std::cout << std::setw(8) << ID + " : " << std::setw(88) << description << endl;
	}
}

void gnut::writeLogInfo(t_glog * log, int l, string description)
{
	if (log)
	{
		log->comment(l, description);
	}
	else
	{
		std::cout << setw(88) << description << endl;
	}
}

void gnut::writeLogInfo(t_glog * log, int l, string ID, string description)
{
	if (log)
	{
		log->comment(l, ID, description);
	}
	else
	{
		std::cout << std::setw(8) << ID + " : " << std::setw(88) << description << endl;
	}
}

void gnut::write_log_info(t_glog * log, int l, string description)
{
	if (log)
	{
		log->comment(l, description);
	}
	else
	{
		std::cout << setw(88) << description << endl;
	}
}

void gnut::write_log_info(t_glog * log, int l, string ID, string description)
{
	if (log)
	{
		log->comment(l, ID, description);
	}
	else
	{
		std::cout << std::setw(8) << ID + " : " << std::setw(88) << description << endl;
	}
}
