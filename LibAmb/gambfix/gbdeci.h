/**
*
* @verbatim
History
-1.0 glfeng  2019-03-05  creat the file.
-1.2 BoWong  2019-04-09  Adding Doxygen Style Code Remarks.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gbdeci.h
* @brief	 Judge whether ambiguities meet the requirement of fixing. 
			 That is judge whether the variance of Amb
			 which have remove upd is less than the threshold.
*
* @author   glfeng BoWong , Wuhan University
* @version	1.0.0
* @date		2019-04-09
*
*/


#ifndef GBDECI_H
#define GBDECI_H 

#include "gutils/ExportAMB.h"
using namespace std;

namespace great
{
	class AMB_LIBRARY_EXPORT t_gbdeci
	{
	public:
		t_gbdeci();
		t_gbdeci(const double& cutdev, const double& cutsig,
			const double& cutprob, const double& cutdeci);
		virtual ~t_gbdeci() {};

		bool bdeci(double dEst, double dSigma, int iIh);
		void bdeci(double dEst, double dSigma, int iIh, double dCutdev, double dCutsig, double& dProb, double& dDeci);

	private:
		double _erf(double dX);
		double _erfc(double dX);
		double _gammp(double dA, double dX);
		double _gammq(double dA, double dX);
		double _gammln(double dXx);
		void _gcf(double& dGammcf, double dA, double dX, double& dGln);
		void _gser(double& dGamser, double dA, double dX, double& dGln);

	private:
		double _cutdev;
		double _cutsig;
		double _cutprob;
		double _cutdeci;
	};
}

#endif
