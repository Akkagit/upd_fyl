/**
*
* @verbatim
    History
     -1.0 Yujie Qin     2019-01-06 creat the file.
     -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
	 -1.2 Yujie Qin     2019-07-03 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			    gsetleo.h
* @brief		    part of LEO settings from XML file
*
* @author		    Yujie Qin, Wuhan University
* @version		    1.0.0
* @date			    2019-03-27
*
*/


#ifndef GSETLEO_H
#define GSETLEO_H

#define XMLKEY_LEO "leo"

#include <vector>
#include <map>
#include <string>
#include <set>
#include "gutils/ExportOI.h"
#include "gset/gsetbase.h"
#include "gutils/gleo.h"

using namespace std;
using namespace pugi;
using namespace gnut;
namespace great
{
	/**
	*@brief	   Class for LEO settings
	*/
	class OI_LIBRARY_EXPORT t_gsetleo : public virtual t_gsetbase
	{
	public:
		/** @brief constructor. */
		t_gsetleo();

		/** @brief destructor. */
		virtual ~t_gsetleo();

		/** @brief settings check. */
		void check(); 

		/** @brief settings help. */
		void help();

		/**
	    * @brief return all LEO satellites in config file(string)
	    * @return  all LEO satellites in config file(string)
	    */
		std::set<std::string> sat();          

		/**
	    * @brief return all LEO satellites in config file(GLEOLONG enum)
	    * @return  all LEO satellites in config file(GLEOLONG enum)
	    */
		std::vector<GLEOLONG> esat();    

	protected:

	};
} // namespace

#endif