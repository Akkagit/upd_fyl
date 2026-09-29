/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2019-10-21  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gpodlsqIF_leokin.h
* @brief		main control class for geometric leo 
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-18
*
*/


#ifndef GPODLSQIF_LEOKIN_H
#define GPODLSQIF_LEOKIN_H

#include "gproc/gppplsq.h"
#include "gutils/ExportLSQ.h"
#include "gcoders/sp3.h"

namespace great
{

	/**
	* @brief main control class for PPP
	*/
	class LSQ_LIBRARY_EXPORT t_gpodlsqIF_leokin :public t_gppplsq
	{
	public:

		t_gpodlsqIF_leokin(string site,t_gsetbase* set, t_gallproc* data, t_glog* log);
		virtual ~t_gpodlsqIF_leokin();

		/// for gengerator kin product file
		bool GenerateProduct() override;

	protected:

		virtual void _setOut() override;
		bool _read_kin();
		virtual bool _initLsqProcData(t_gallproc* data) override;
		virtual int  ProcessOneEpoch(const t_gtime& now, vector<t_gsatdata>* obs = NULL) override;
		bool _init_rec_crd_pars(t_glsq* lsq, bool isPPP) override;

		virtual bool _processOneRec(const t_gtime& crt_epoch,

			const std::string& crt_rec,

			std::vector<t_gsatdata>& crt_obs) override;

		void _extract_kin();
		void _write_kin_header();
		void _write_kin_data(t_gtime epo, map<string, double> kin_map);
		void _write_kin_end();
		bool _init_sp3(t_gallprec* sp3);
		void _add_spp_sp3(t_gallprec* sp3);
		bool _add_ppp_sp3(t_gallprec* sp3);
		bool _generate_spp_kin();
		

	protected:

		double _sampling;
		int    _obs_limit;
		map<t_gtime, int> _obs_num_spp;
		t_gallprec* _newsp3 = nullptr;
		t_gallprec* _sppsp3 = nullptr;
		t_giof*	_result_kin;
		vector<string> _prn3;
		//t_gallkin* _allkin;   /// storage for the kin data

	};




}

#endif 
/*  GGEOMETRIC_LEO_H */