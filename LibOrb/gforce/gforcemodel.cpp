/**
 *
 * @verbatim
	 History
	  -1.0 Jiande Huang 2019-01-14 creat the file.
	  -1.1 Jie Li  2019-04-01 add code in doxygen style.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gforcemodel.cpp
 * @brief		The container of all the force models used in space craft orbit intergration.
 *
 * @author      Jiande Huang, Wuhan University
 * @version		1.0.0
 * @date		2019-01-14
 */
#include "gforcemodel.h"


 /** @brief clear the force model and delete the point.
 * @param[in]  ptr	 force pointer need to be delete
 */
void great::t_gforce_model::destoryForceModel(t_gforce * ptr)
{
	delete ptr;
	ptr = NULL;
}
