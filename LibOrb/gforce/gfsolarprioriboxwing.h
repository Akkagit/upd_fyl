/**
 *
 * @verbatim
     History
      -1.0 Yongqiang Yuan     2019-11-02 Creat the file.
   @endverbatim
 * Copyright (c) 2019, Wuhan University.All rights reserved.
 *
 * @file		      gfsolarprioriboxwing.h
 * @brief	          general solar radiation model, calculate as sum of accelerations all pannels.
 *
 * @author            Yongqiang Yuan, Wuhan University
 * @version	          1.0.0
 * @date		      2019-11-02
 *
 */
#ifndef GFSOLARPRIORIBOXWING_H
#define GFSOLARPRIORIBOXWING_H

#include <string>
#include "gforce/gfsolar.h"
#include "gforce/gfarg_init.h"
#include "gio/glog.h"
#include "gutils/gconst_oi.h"
#include "gutils/ExportOI.h"

using namespace std;

namespace great
{
	/** @brief class for for BOXWING force model. */
	class OI_LIBRARY_EXPORT t_gfsolar_pbw : public t_gfsolar
	{
	public:

		/** @brief constructor. */
		t_gfsolar_pbw(t_init_arg* arg, t_glog* log);

		/** @brief  destructor. */
		virtual ~t_gfsolar_pbw();

		/** @brief  show the type of model. */
		void show() override;

		/** @brief  check some important pars in the model. */
		bool check() override;

		/** 
		* @brief calculate the acceleration and the corresponding partial matrix.
	    *@param[in] argument Parameters to be used in calculation.
	    */
		bool calculate(t_gforce_calarg* argument) override;
		
		
	protected:
		std::map<ICS_PAR, double>	_pbw_pars;	       ///< 9 solar pars of ECOM
		t_gpannel*            _pannel = nullptr;                 ///< satellite pannel information.
		const double          _ps        = 4.5604e-6;
		const double          _aunit     = 1.4959787061e8;
		
	};
} // namespace

#endif