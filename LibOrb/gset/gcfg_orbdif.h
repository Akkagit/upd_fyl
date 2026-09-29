/**
*
* @verbatim
	History
	 -1.0 Jie Li	 2019-09-05 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gcfg_orbdif.h
* @brief		orbit differ control set from XML
*
* @author		Jie Li, Wuhan University
* @version		1.0.0
* @date			2019-09-05
*
*/


#ifndef GCFG_ORBDIF
#define GCFG_ORBDIF

#include "gset/gsetgen.h"
#include "gset/gsetgnss.h"
#include "gset/gsetinp.h"
#include "gset/gsetout.h"
#include "gset/gsetrec.h"
#include "gset/gsetorbdif.h"


namespace great
{
	/**
	*@brief	   Class for GNS orbit differ contrl setting
	*/
	class OI_LIBRARY_EXPORT t_gcfg_orbdif : public t_gsetgen,
		public t_gsetinp,
		public t_gsetout,
		public t_gsetgnss,
		public t_gsetrec,
		public t_gsetorbdif
	{
	public:
		/** @brief  constructor. */
		t_gcfg_orbdif();

		/** @brief  destructor. */
		virtual ~t_gcfg_orbdif();

		/** @brief settings check. */
		void check();

		/** @brief settings help. */
		void help();

	protected:

	private:
	};

}
#endif // !GCFG_OI_H
