/**
*
* The format :
* @verbatim
		<orbdif>
			<excsat> exclude_sat </excsat>
			<excsys> exclude_sat </excsys>
			<trans> STRD </trans>
		</orbdif>
  @endverbatim

* @verbatim
	History
	 -1.0 Jie Li	 2019-09-05 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gsetintgra.h
* @brief		orbit differ parameter  set for satellite orbit differ from XML
*
* @author       Jie Li, Wuhan University
* @version		1.0.0
* @date			2019-09-05
*
*/

#ifndef GSETORBDIF_H
#define GSETORBDIF_H

#define XMLKEY_ORBDIF		  "orbdif"
#define XMLKEY_ORBFIT		  "orbfit"								 
#define XMLKEY_ORBDIF_EXCSAT  "excsat"
#define XMLKEY_ORBDIF_EXCSYS  "excsys"
#define XMLKEY_ORBDIF_TRANS   "trans"


#include "gset/gsetbase.h"
#include "gutils/ExportOI.h"
using namespace gnut;

namespace great
{
	class OI_LIBRARY_EXPORT t_gsetorbdif : public virtual t_gsetbase
	{
	public:
		/** @brief constructor. */
		t_gsetorbdif();

		/** @brief destructor. */
		virtual ~t_gsetorbdif();

		/**
		* @brief get the trans parameteres of orbit differ.
		* @return the trans parameteres(STRD)
		*/
		string trans();
		// added by yqyuan for orbfit
		string trans_orbfit();				   

		/**
		* @brief get the satellite systems that needn't differ.
		* @return the exclude system
		*/
		set<string> excsys();
		// added by yqyuan for orbfit
		set<string> excsys_orbfit();				   


		/**
		* @brief get the satellites that needn't differ.
		* @return the exclude sat
		*/
		set<string> excsat();
		// added by yqyuan for orbfit
		set<string> excsat_orbfit();				   


		/** @brief settings check. */
		void check();

		/** @brief settings help. */
		void help();
		// added by yqyuan for orbfit
		void help_orbfit();				   


	private:
		
	};
}
#endif
