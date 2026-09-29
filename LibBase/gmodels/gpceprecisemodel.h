/**
*
* @verbatim
	History
	 -1.0 Xiongyun  2019-12-10  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gpceprecisemodel.h
* @brief		precise model for computing  correction
*
* @author       Xiongyun, Wuhan University
* @version		1.0.0
* @date		    2019-12-10
*
*/


#ifndef GPCEPRECISEMODEL_H
#define	GPCEPRECISEMODEL_H

#include "gmodels/gprecisemodel.h"

//#define DEBUG_LOG

namespace great
{
	class BASE_LIBRARY_EXPORT t_gpceprecisemodel : public t_gprecisemodel
	{
	public:

		t_gpceprecisemodel(t_gallproc* data, t_glog* log, t_gsetbase* setting);

		~t_gpceprecisemodel();

		bool _cmb_equ(const t_gtime& epoch, t_gsatdata& crt_obs, t_gallpar& pars, t_glsqEquationMatrix& equ);
		bool _cmb_equ_IF(const t_gtime& epoch, t_gallpar& par, t_glsqEquationMatrix& equ);
		bool _cmb_equ_ALL(const t_gtime& epoch, t_gallpar& par, t_glsqEquationMatrix& equ);
		bool _omc_obs_IF(const t_gtime& crt_epo, t_gobs& gobs, t_gallpar& pars, double& omc);
		bool _omc_obs_IF13(const t_gtime& crt_epo, t_gobs& gobs, t_gallpar& pars, double& omc);
		bool _omc_obs_ALL(const t_gtime& crt_epo, t_gobs& gobs, t_gallpar& pars, double& omc);
		bool _wgt_obs(t_gobs& gobs, double factorP, double& wgt, bool isTriple = false);
		bool _update_obs_info(const t_gtime& epoch, t_gsatdata& obsdata, t_gallpar& pars);
		bool _update_obj_clk(const string& obj, const t_gtime& epo, t_gallpar& par,
			double& clk, map<string, pair<t_gtime, double> >& obj_clk, GSYS& ctr_sys);

		double cmpObs(t_gtime& epoch, t_gallpar& param, t_gsatdata& obsdata, t_gobs& gobs, bool com);

	protected:
		bool _carrier_range; ///< use carrier range data [true=1, false=0]
		bool   _IRC;         ///< estimate integer recovery clock [true=1, false=0]

	};
}

#endif
