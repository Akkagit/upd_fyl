/**
*
* @verbatim
	History
	 -1.0 zhangwei  2019-03-18  creat
	 -1.1 Qinyujie  2019-03-20  optim
	 -1.2 Qinyujie  2019-04-11  add Doxygen style code
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gorblsq_leo.h
* @brief		lsq processing for leo
*
* @author       zhangwei, Wuhan University
* @version		1.0.0
* @date		    2019-03-18
*
*/
#ifndef _GORBLSQLEO_H
#define _GROBLSQLEO_H

#include "gutils/ExportLSQ.h"

#include <algorithm>

#include "gproc/glsq.h"
#include "gmodels/glsqprocIF.h"
#include "gmodels/gprecisemodel.h"
#include "gall/gallics_leo.h"
#include "gdata/gorb.h"
#include "gdata/gpoleut1.h"
#include "gdata/gleapsecond.h"
#include "gall/gallattitude.h"
using namespace std;

namespace great
{
	/**
	*@brief	   Class for LEO LSQ processing.
	*/
	class LSQ_LIBRARY_EXPORT t_gorblsq_leo : public t_glsqprocIF
	{
	public:
		/** @brief default constructor. */
		t_gorblsq_leo(t_gsetbase* base, t_gallproc* data, t_glog* log);
		
		/** @brief default destructor. */
		virtual ~t_gorblsq_leo();

	public:

		virtual bool GenerateProduct() override;

	protected:

		t_gallics_leo* _icsleo;        ///< ics data of the leo
		set<string> _satleo;          ///< long name of leo list to be processed
		t_giof*	_result;              /// output file of result

		/** @brief init the lsq parameter. */
		virtual bool _InitLsq(t_gallproc* data) ;

		/** @brief init the output file. */
		void _setOut();

		/** @brief init LEO orbit parameters. */
		void _init_ORB_parameter();

		/** @brief init other parameters except orbit parms. */
		void _init_parameter();

		/** @brief upsate LEO ics information. */
		void _update_icsleo();

	};
}

#endif /* _GORBLSQ_H */