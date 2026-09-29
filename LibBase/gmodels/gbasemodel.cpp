/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2020-05-15 created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    t_gabasemode.cpp
* @brief	    base interface model
* @mainpage     C++ version
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-05-18
*
*/

#include "gmodels/gbasemodel.h"


namespace great
{
    t_gbaseEquation::t_gbaseEquation()
    {
    }

    t_gbaseEquation::~t_gbaseEquation()
    {

    }

    t_gbaseEquation& t_gbaseEquation::operator+(t_gbaseEquation &Other)
    {
        this->B.insert(this->B.begin(), Other.B.begin(), Other.B.end());
        this->P.insert(this->P.begin(), Other.P.begin(), Other.P.end());
        this->l.insert(this->l.begin(), Other.l.begin(), Other.l.end());
        return *this;
    }

    t_gbasemodel::t_gbasemodel()
    {
    }

    t_gbasemodel::~t_gbasemodel()
    {
    }
}