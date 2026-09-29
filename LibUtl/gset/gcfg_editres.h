/**
*
* @verbatim
	History
	 -1.0 Hongjie Zheng		2019-09-06  creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gcfg_editres.h
* @brief		edit residual control set from XML
*
* @author       Hongjie Zheng, Wuhan University
* @version		1.0.0
* @date		    2019-09-06
*
*/

#ifndef GCFG_EDITRES_H
#define GCFG_EDITRES_H

#include "gutils/ExportUtl.h"
#include "gset/gsetinp.h"
#include "gset/gsetout.h"
#include "gset/gsetres.h"
#include "gset/gsetgen.h"
#include "gset/gsetgnss.h"

namespace great
{
	/**
	*@brief	   Class for UPD contrl setting
	*/
	class UTL_LIBRARY_EXPORT t_gcfg_editres :
		public t_gsetinp,
		public t_gsetout,
		public t_gsetgen,
		public t_gsetgnss,
		public t_gsetres
	{
	public:

		/** @brief default constructor. */
		t_gcfg_editres();

		/** @brief default destructor. */
		virtual ~t_gcfg_editres();

		/** 
		* @brief settings check. 
		* @return void
		*/
		void check();                                 // settings check

		/** 
		* @brief settings help. 
		* @return void
		*/
		void help();                                  // settings help

	protected:

	private:
	};

}

#endif 