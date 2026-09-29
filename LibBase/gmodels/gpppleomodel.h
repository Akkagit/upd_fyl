/**
*
* @verbatim
History
-1.0 xjhan  2019-11-22  creat the file.

@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gpppleomodel.h
* @brief	precise model for ppp with leo
*
* @author   xjhan , Wuhan University
* @version	1.0.0
* @date		2019-11-22
*
*/

#ifndef GPPPLEOMODEL_H
#define	GPPPLEOMODEL_H

#include "gmodels/gprecisemodel.h"

// add leo
//#define DEBUG

namespace gnut {

	class BASE_LIBRARY_EXPORT t_gpppleomodel :public t_gprecisemodel
	{
	public:
		t_gpppleomodel(t_gallproc* data, t_glog* log, t_gsetbase* setting);
		virtual ~t_gpppleomodel();

		// override for wind up correction
		double windUp(t_gsatdata& satdata, const ColumnVector& rRec) override;

		//double get_raw_dcb(string sat, GOBS obs) { return _simu_raw_dcb[sat][obs]; }
		void apply_GaussNoiseP(string sat, double& obs);
		void apply_GaussNoiseL(string sat, double& obs);

	protected:
		double _noiseP_next;
		double _noiseL_next;
		map<string, bool> _sat_noise_generateP;
		map<string, bool> _sat_noise_generateL;
		//map<string, map<GOBS, double> >   _simu_raw_dcb;
		//map<string, int>                  _glofrq_num;
	};

	bool BASE_LIBRARY_EXPORT raw_CODE_bias(const t_gtime epo, t_gallbias* allbias, map<string, int>& R_freq_ID, map<string, map<GOBS, double>>& raw_bia);
}


#endif /* GPPPLEOMODEL_H */
