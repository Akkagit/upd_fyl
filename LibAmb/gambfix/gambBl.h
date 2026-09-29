/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambBl.h
* @brief	baseline of SD and DD ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#ifndef GAMBBL_H
#define GAMBBL_H 

#include "gutils/ExportAMB.h"
#include "gio/glog.h"
#include "gutils/gnss.h"

using namespace std;
using namespace gnut;

namespace great
{
	// 双差模糊度
	class AMB_LIBRARY_EXPORT t_gamb_baseline
	{
	public:
		// 默认构造函数
		t_gamb_baseline();
		t_gamb_baseline(string obj1, string obj2, double distance, int nsys);
		virtual ~t_gamb_baseline();

		pair<string, string> bl_pair() const;


		string obj1; 
		string obj2; 
		int    nsys; 
		double distance; 

	public:
		bool isSame(const string& obj_1st, const string& obj_2nd);
	};
}

#endif