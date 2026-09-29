/**
*
* @verbatim
	History
	 -1.0 Hongmin Zhang     2019-11-20 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			    gsetslr.h
* @brief		    part of SLR settings from XML file
*
* @author		    Hongmin Zhang, Wuhan University
* @version		    1.0.0
* @date			    2019-11-20
*
*/


#ifndef GSETSLR_H
#define GSETSLR_H

#define XMLKEY_SLR "slr"

#include <vector>
#include <map>
#include <string>
#include <set>
#include "gutils/ExportUtl.h"
#include "gset/gsetbase.h"

using namespace std;
using namespace pugi;
using namespace gnut;
namespace great
{
	/**
	*@brief	   Class for LEO settings
	*/
	class UTL_LIBRARY_EXPORT t_gsetslr : public virtual t_gsetbase
	{
	public:
		/** @brief constructor. */
		t_gsetslr();

		/** @brief destructor. */
		virtual ~t_gsetslr();

		/** @brief settings check. */
		void check();

		/** @brief settings help. */
		void help();

		
		std::string getSLRXmlSet();

	protected:

	};
} // namespace

#endif