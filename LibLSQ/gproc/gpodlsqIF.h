/**
*
* @verbatim
	History
	 -1.0 zhangwei     2019-03-18  creat
	 -1.1 Qinyujie     2019-03-20  optim
	 -1.2 Qinyujie     2019-04-11  add Doxygen style code
	 -1.3 Huang Jiande 2019-07-12  add UC model
	 -1.4 Huang Jiande 2020-03-14  fix UC model bug
	 -1.5 Huang Jiande 2020-09-11  format the function
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
#ifndef GPODLSQIF_H
#define GPODLSQIF_H

#include "gutils/ExportLSQ.h"
#include "gambfix/gallamb.h"
#include "gdata/gambcon.h"

#include "gmodels/glsqproc.h"

using namespace std;

namespace great
{
	/**
	*@brief	   Class for LEO LSQ processing.
	*/
	class LSQ_LIBRARY_EXPORT t_gpodlsqIF : public t_glsqproc
	{
	public:

		t_gpodlsqIF();

		/** @brief default constructor. */
		t_gpodlsqIF(t_gsetbase* base, t_gallproc* data, t_glog* log);
		
		/** @brief default destructor. */
		virtual ~t_gpodlsqIF();

		bool ProcessBatch(t_gallproc* data, const t_gtime&beg, const t_gtime& end) override;
		bool GenerateProduct() override;

		bool initLsqProcData(t_gallproc* data) { return _initLsqProcData(data); };
		bool initLsqProcPars(t_glsq* lsq) { return _initLsqProcPars(lsq); };
		void get_epoch_solution() { _get_epoch_solution(); };
		void get_correlation_matrix() { _get_correlation_matrix(); };
		void change_update_mode(t_gallproc* eop_data, t_gallpar& pars){ _change_update_mode(eop_data, pars); };

	protected:

		bool _initLsqProdData(t_gallprod* data) override;
		bool _initLsqProcPars(t_glsq* lsq)      override;
		bool _initLsqProcData(t_gallproc* data) override;

		void _update_nav_ics();
		void _update_leo_ics();

		void _get_epoch_solution();
		void _get_correlation_matrix();

		void _print_out_header();
		void _print_out_content(const t_glsq& lsqestimator);
		void _readOfiles(int step,t_gallobs*& obs);

		void _set_correlation_out();
		void _change_update_mode(t_gallproc* eop_data, t_gallpar& pars);
		void _extract_corr_map(SymmetricMatrix& Qnew, vector<string>& name);

		t_giof*	           _result = nullptr;		///< outputs every epoch result file
		t_giof*	           _corr_map = nullptr;

	};
}

#endif /* _GPODLSQIF_LEO_H */