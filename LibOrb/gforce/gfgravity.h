/**
 *
 * @verbatim
	 History
	  -1.0 Jiande Huang 2019-01-14 creat the file.
      -1.1 Jiande Huang 2019-03-28 Add the construtor and destructor.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gfgravity.h
 * @brief		The base class for gravity model, following could be EGM08_12 or EGM08_150 and so on.
 *
 * @author      Jiande Huang, Wuhan University
 * @version		1.0.0
 * @date		2019-01-14
 */

#ifndef GFGRAVITY_H
#define GFGRAVITY_H

#include "gutils/ExportOI.h"
#include "gio/glog.h"
#include "gforce/gforce.h"

namespace great
{
	/** 
	 * @brief This is the base class for gravity model, 
		 other gravity field models will inherit this class.
	 */	
	class OI_LIBRARY_EXPORT t_gfgravity : public t_gforce
	{
	public:
		/** 
		* @brief constructor. 
		* @param [in]  log		log file
		*/
		t_gfgravity(t_glog* log);

		/** @brief default destructor. */
		virtual ~t_gfgravity() {};

		/** 
		 * @brief The function will show you the information of this gravity model, 
		 *	such as its order and model type. 
		 */
		virtual void show() = 0;
	};

}
#endif // !GFGRAVITY