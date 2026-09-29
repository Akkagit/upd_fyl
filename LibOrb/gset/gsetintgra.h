/**
*
* The format : 
* @verbatim
        <integrator>
        	<rkf    order = "7"  rkf_num  = "30"                 />
        	<adams  order = "11" inp_step = "60" out_step = "300"/>
        </integrator>
  @endverbatim

* @verbatim
	History
	 -1.0 Jiande Huang	 2019-04-07 creat the file.
	 -1.1 Yujie Qin      2019-07-04 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gsetintgra.h
* @brief		number integrator  set for satellite orbit integral from XML
*
* @author       Jiande Huang, Wuhan University
* @version		1.0.0
* @date			2019-04-07
*
*/

#ifndef GSETINTGR_H
#define GSETINTGR_H

#define XMLKEY_INTEGRATOR       "integrators"
#define XMLKEY_INTEGRATOR_RKF   "rkf"
#define XMLKEY_INTEGRATOR_ADAMS "adams"

#include "gset/gsetbase.h"
#include "gutils/ExportOI.h"

using namespace gnut;

namespace great
{
	class OI_LIBRARY_EXPORT t_gsetintgra : public virtual t_gsetbase
	{
	public:
		/** @brief constructor. */
		t_gsetintgra();

		/** @brief destructor. */
		virtual ~t_gsetintgra();

		/** 
		* @brief get the name of rkf model(including the order of rkf). 
		* @return the name of rkf
		*/
		string rkf_model();

		/** 
		* @brief get the name of adams model(including the order of adams). 
		* @return the name of adams
		*/
		string adams_model();

		/** 
		* @brief get the number of rkf points. 
		* @return the number of rkf points
		*/
		int rkf_num();

		/** 
		* @brief get the order of rkf. 
		* @return the order of rkf
		*/
		int rkf_order();

		/** 
		* @brief get the order of adams. 
		* @return the order of adams
		*/
		int adams_order();

		/** 
		* @brief get the input step of adams. 
		* @return the input step of adams
		*/
		double adams_inp_step();

		/** 
		* @brief get the output step of adams. 
		* @return the output step of adams
		*/
		double adams_out_step();

		/** @brief settings check. */
		void check();

		/** @brief settings help. */
		void help();

	//protected:
	//	int _rkf_num;
	//	int _rkf_order;
	//	int _adams_order;
	//	int _adams_inp_step;
	//	int _adams_out_step;

	private:
		// some defalut num, list here for the change convenience
		const int _default_rkf_num        = 30;       ///< the default number of rkf points
		const int _default_rkf_order      = 7;        ///< the default order of rkf
		const int _default_adams_order    = 11;       ///< the default order of adams
		const double _default_adams_inp_step = 60;    ///< the default input step of adams(unit:seconds)
		const double _default_adams_out_step = 300;   ///< the default output step of adams(unit:seconds)
	};
}
#endif
