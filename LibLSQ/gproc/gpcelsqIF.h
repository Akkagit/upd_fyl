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


#ifndef GPCELSQIF_H
#define GPCELSQIF_H

#include "gmodels/glsqprocIF.h"
#include "gutils/ExportLSQ.h"
#include "gset/gsetturboedit.h"
#include "gproc/gexeturboedit.h"
#include "gio/gtcp.h"
#include "gcoders/sp3.h"

namespace great
{

	/**
	* @brief main control class for PPP
	*/
	class LSQ_LIBRARY_EXPORT t_gpcelsqIF :public t_glsqproc
	{
	public:

		t_gpcelsqIF(t_gsetbase* set, t_gallproc* data, t_glog* log);
		virtual ~t_gpcelsqIF();

		bool ProcessBatch(t_gallproc* data, const t_gtime&beg, const t_gtime& end) override;

	protected:

		t_gpceprecisemodel* _bias_pcemodel =nullptr;
		map<string,shared_ptr<t_gexeturboedit> > _gturbo;
		t_gio* _gioout = nullptr;

		bool   _check_out = false;
		bool   _upload = false;     //added by xiongyun
		string _sp3dir;             //added by xiongyun
		string _sp3name;            //added by xiongyun
		map<string, t_glsqEquationMatrix> _map_all_equ;
		vector<string> _invalid_site;
		vector<pair<string, string> > _outlier_data;

		bool _initLsqProdData(t_gallprod* data) override;
		bool _initLsqProcPars(t_glsq* lsq)      override;
		bool _initLsqProcData(t_gallproc* data) override;
		bool _update_ultrasp3(int time_gap);

		bool _processOneRecOneSat(const t_gtime& crt_epoch,const std::string& rec,
			const std::string& sat,
			t_gsatdata& crd_obs) override;

		void _print_clk_solved(t_glsq* lsq, t_gtime& crt_epoch, bool isTriple = false);
		void _print_rtnet_clk(t_glsq* lsq, t_gtime& crt_epoch);
		void set_recover_data_by_resfile(string resfile_name);
		bool outlierDetect(t_glsq& lsq, map<string, t_glsqEquationMatrix>& all_equ);
		void outilerDelete();

		bool _processOneEpoch(const t_gtime& crt_epoch,
			std::vector<t_gsatdata>& crt_obs) override;

		bool _processOneRec(const t_gtime& crt_epoch,
			const std::string& crt_rec,
			std::vector<t_gsatdata>& crt_obs) override;

		bool _processOneRec_thread_safe(const t_gtime& crt_epoch,
			const std::string& crt_rec,
			const std::vector<t_gsatdata>& crt_obs,t_glsqEquationMatrix& equ_result);

		bool GenerateProduct() override;
	};




}

#endif /*  GPPPLSQ_H */