/**
*
* @verbatim
	History
	 -1.0 liyuhao	 2020-07-28 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gcfg_clkdif.h
* @brief		clock differ control set from XML
*
* @author		liyuhao, Wuhan University
* @version		1.0.0
* @date			2020-07-28
*
*/

#ifndef GCFG_CLKDIF_H
#define GCFG_CLKDIF_H

#include "gutils/ExportUtl.h"
#include "gset/gsetinp.h"
#include "gset/gsetout.h"
#include "gset/gsetres.h"
#include "gset/gsetgen.h"
#include "gset/gsetgnss.h"

namespace great
{
	class UTL_LIBRARY_EXPORT t_gcfg_clkdif :
		public t_gsetinp,
		public t_gsetout,
		public t_gsetgen,
		public t_gsetgnss
	{
	public:

		/** @brief default constructor. */
		t_gcfg_clkdif();

		/** @brief default destructor. */
		virtual ~t_gcfg_clkdif();

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


