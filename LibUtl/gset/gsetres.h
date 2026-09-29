/*

*/
/**
*

* @verbatim
	================================================================<!-->
	<xml>
		 <editres>
			<freq>        LC  </freq>		LC/L1/L2/L3/L4/L5
			<short_elisp> 600 </short_elisp>
			<jump_elisp>  40  </jump_elisp> 80/40
			<bad_elisp>   40  </bad_elisp>
         </editres>	

	</xml>
 @endverbatim

* @verbatim
	History
	 -1.0 Hongjie Zheng		2019-09-06  creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gsetres.h
* @brief		the xml set for edit residual
*
* @author       Hongjie Zheng, Wuhan University
* @version		1.0.0
* @date			2019-09-06
*
*/
#ifndef GSETRES_H
#define GSETRES_H

#define XMLKEY_RES "editres"

#include "gset/gsetbase.h"
#include "gutils/ExportUtl.h"

using namespace std;
using namespace gnut;
namespace great
{

	class UTL_LIBRARY_EXPORT t_gsetres : public virtual t_gsetbase
	{
	public:

		/** @brief default constructor. */
		t_gsetres();

		/** @brief default destructor. */
		virtual ~t_gsetres();

		/**
		* @brief settings check.
		* @return void
		*/
		void check();

		/**
		* @brief settings help.
		* @return void
		*/
		void help();

		double get_jump_elisp();
		double get_bad_elisp();
		double get_short_elisp();
		double get_sta_num();
		double get_res_max();
		string get_mode();
		string get_freq();
		string get_freq_mode();
	protected:
		const string _default_freq;
		const string _default_mode;
		const double _default_jump_elisp;
		const double _default_bad_elisp;
		const double _default_short_elisp;
		const double _default_sta_num;
		const double _default_res_max;


	private:

	};



}
#endif
