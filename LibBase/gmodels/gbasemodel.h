/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2020-05-15 created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    t_gabasemode.h
* @brief	    base interface model  
* @mainpage     C++ version 
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-05-18
*
*/
#ifndef GBASEMODEL_H
#define GBASEMODEL_H

#include <vector>
#include "gutils/gtime.h"
#include "gall/gallpar.h"
#include "gdata/gsatdata.h"

using namespace std;
using namespace gnut;

namespace great
{
    class BASE_LIBRARY_EXPORT t_gbaseEquation
    {
    public:
        t_gbaseEquation();
        virtual ~t_gbaseEquation();

        t_gbaseEquation &operator+(t_gbaseEquation &Other);

        vector<vector<pair<int, double>>> B;
        vector<double> P;
        vector<double> l;
    };

    // interface class for all mode 
    class BASE_LIBRARY_EXPORT t_gbasemodel
    {
    public:
        t_gbasemodel();
        ~t_gbasemodel();

        virtual bool cmb_equ(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,t_gbaseEquation& result)=0;
    };
}

#endif