/**
 *
 * @verbatim
     History
      -1.0 Wei Zhang     2019-01-20 Creat the file.
      -1.1 Yujie Qin     2019-03-04 edit for final version.
      -1.2 Yujie Qin     2019-03-15 Code comments are added in Doxygen format.
	  -1.3 Jiande Huang   2019-03-28 fix some proble with linux make up, such as change #ifndef GFSOLARBOXW_H to GFSOLARBOXWING_H and so on.
	  -1.4 Yujie Qin     2019-07-03 Adding Doxygen Style Code Remarks
   @endverbatim
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		      gfsolarboxwing.h
 * @brief	          general solar radiation model, calculate as sum of accelerations all pannels.
 *
 * @author            Wei Zhang, Wuhan University
 * @version	          1.0.0
 * @date		      2019-01-20
 *
 */
#ifndef GFSOLARBOXWING_H
#define GFSOLARBOXWING_H

#include "gutils/ExportOI.h"

#include <string>
#include "gforce/gfsolar.h"
#include "gforce/gfarg_init.h"
#include "gio/glog.h"
#include "gutils/gconst_oi.h"

using namespace std;

namespace great
{
	/** @brief class for for BOXWING force model. */
	class OI_LIBRARY_EXPORT t_gfsolarbow : public t_gfsolar
	{
	public:

		/** @brief constructor. */
		t_gfsolarbow(t_init_arg* arg, t_glog* log);

		/** @brief  destructor. */
		virtual ~t_gfsolarbow();

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
		map<ICS_PAR, double>  _boxw_par;               ///< may be have a solar par of BOXWING.
		t_gpannel*            _pannel = nullptr;       ///< satellite pannel information.
		bool                  _ispartial = false;      ///< calculate the partial or not.
		double                _mass      = 0.0;        ///< the mass of sat.
		const double          _ps        = 4.5604e-6;
		const double          _aunit     = 1.4959787061e8;
		const ICS_PAR         _parname   = ICS_PAR::SR_scale;   ///< the name of solar par.
		
	};
} // namespace

#endif