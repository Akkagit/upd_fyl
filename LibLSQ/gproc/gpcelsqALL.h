/**
*
* @verbatim
	History
	 -1.0 Jiande Huang  2019-07-08  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gpcelsqIF.h
* @brief		main control class for ppp estimate
*
* @author       Jiande Huang, Wuhan University
* @version		1.0.0
* @date		    2019-04-18
*
*/


#ifndef GPCELSQALL_H
#define GPCELSQALL_H

#include "gmodels/glsqprocALL.h"
#include "gutils/ExportLSQ.h"

namespace great
{

	/**
	* @brief main control class for PPP
	*/
	class LSQ_LIBRARY_EXPORT t_gpcelsqALL :public t_glsqproc
	{
	public:

		t_gpcelsqALL();
		t_gpcelsqALL(t_gsetbase* set, t_gallproc* data, t_glog* log);
		virtual ~t_gpcelsqALL();

		bool ProcessBatch(t_gallproc* data, const t_gtime&beg, const t_gtime& end) override;

	protected:

		t_gpceprecisemodel* _bias_pcemodel = nullptr;

		bool _initLsqProdData(t_gallprod* data) override;
		bool _initLsqProcPars(t_glsq* lsq)      override;
		bool _initLsqProcData(t_gallproc* data) override;

		//bool _processOneRecOneSat(const t_gtime& crt_epoch,const std::string& rec,const std::string& sat,
		//	t_gsatdata& crd_obs,int& obs_num) override;
		bool _processOneEpoch(const t_gtime& crt_epoch, std::vector<t_gsatdata>& crt_obs) override;
		bool _processOneRec(const t_gtime& crt_epoch, const std::string& crt_rec, std::vector<t_gsatdata>& crt_obs) override;
		bool _processOneRecOneSat(const t_gtime& crt_epoch, const std::string& rec,
			const std::string& sat,
			t_gsatdata& crd_obs) override;

		bool GenerateProduct() override;

		void _print_clk_solved(t_glsq* lsq, t_gtime& crt_epoch, bool isTriple = false);
	};




}

#endif /*  GPPPLSQ_H */