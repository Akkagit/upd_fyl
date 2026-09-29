/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2020-05-15 created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    t_gcombmodel.h
* @brief	    base combine biase model
* @mainpage     C++ version 
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-05-18
*
*/

#ifndef GCOMBMODEL_H
#define GCOMBMODEL_H

#include "gmodels/gbasemodel.h"
#include "gmodels/gbiasmodel.h"
#include "gmodels/gprecisemodel.h"
#include "gmodels/giono.h"
#include "gall/gallbias.h"
#include "gall/gallproc.h"
#include "gutils/gcommon.h"
#include <map>
#include <memory>

using namespace gnut;
using namespace std;

namespace great
{
    class BASE_LIBRARY_EXPORT t_gcombmodel:virtual public t_gbasemodel
    {
    public:
        t_gcombmodel(t_gsetbase* setting,t_glog* log,shared_ptr<t_gbiasmodel> bias_model,t_gallproc* data);
        virtual ~t_gcombmodel();

		map< GSYS, map<FREQ_SEQ, GOBSBAND> > get_band_index();
		map< GSYS, map<GOBSBAND, FREQ_SEQ> > get_freq_index();
    protected:

        bool _apply_dcb_in_obsdata(t_gsatdata& obsdata);

        map< GSYS, map<FREQ_SEQ, GOBSBAND> > _band_index;
		map< GSYS, map<GOBSBAND, FREQ_SEQ> > _freq_index;


        // baise model
        shared_ptr<t_gbiasmodel> _bias_model;

        t_gprecisemodel* _gprecisemodel = nullptr;
        t_gallbias*      _gallbias      = nullptr;
        t_glog*          _log           = nullptr;

        int             _frequency;
        IFB_MODEL       _ifb_model;
		IFCB_MODEL      _ifcb_model;
        IONMODEL        _ion_model;
		bool            _update_amb_lite;
		
    };

    class BASE_LIBRARY_EXPORT t_gcombIF:virtual public t_gcombmodel
    {
    public:
        t_gcombIF(t_gsetbase *setting,t_glog* log,shared_ptr<t_gbiasmodel> bias_model,t_gallproc* data);
        ~t_gcombIF();

        bool cmb_equ(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,t_gbaseEquation& result) override;
        bool cmb_equ_IF(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,GOBSBAND b1,GOBSBAND b2, t_gbaseEquation& result) ;

    };

    class BASE_LIBRARY_EXPORT t_gcombIF13:virtual public t_gcombmodel
    {
    public:
        t_gcombIF13(t_gsetbase *setting,t_glog* log,shared_ptr<t_gbiasmodel> bias_model,t_gallproc* data);
        ~t_gcombIF13();

        bool cmb_equ(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,t_gbaseEquation& result) override;

    protected:
        t_gcombIF _cmb_IF;
    };

    class BASE_LIBRARY_EXPORT t_gcombALL:virtual public t_gcombmodel
    {
    public:
        t_gcombALL(t_gsetbase *setting,t_glog* log,shared_ptr<t_gbiasmodel> bias_model,t_gallproc* data);
        ~t_gcombALL();

        bool cmb_equ(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, t_gbaseEquation& result) override;

    private:
        map<FREQ_SEQ, t_gpar::t_type> ambtype_list;
		
	};

	class BASE_LIBRARY_EXPORT t_gcombWL :public t_gcombmodel
	{
	public:
		t_gcombWL(t_gsetbase* setting, t_glog* log, shared_ptr<t_gbiasmodel> bias_model, t_gallproc* data);
		~t_gcombWL();

		bool cmb_equ(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, t_gbaseEquation& result) override;
		bool cmb_equ_WL(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, GOBSBAND b1, GOBSBAND b2, t_gbaseEquation& result);
	};

	class BASE_LIBRARY_EXPORT t_gcombDD :virtual public t_gcombALL, virtual public t_gcombIF
	{
	public:
		t_gcombDD(t_gsetbase *setting, t_glog* log, shared_ptr<t_gbiasmodel> bias_model, t_gallproc* data);
		~t_gcombDD();

		bool cmb_equ(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, t_gbaseEquation& result) override;

		void set_observ(OBSCOMBIN observ);
		void set_base_data(vector<t_gsatdata>* data_base);
		void set_site(string site, string site_base);
		void set_rec_info(t_gtriple xyz_base, double clk_rover, double clk_base);
		bool cmb_equ_test(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, t_gbaseEquation& result);


	protected:
		bool _temp_params(t_gallpar& params, t_gallpar& params_temp);
		OBSCOMBIN _observ;
		vector<t_gsatdata> *_data_base;
		string _site, _site_base;
		double _clk_rover, _clk_base;
        t_gtriple _crd_base;
		
	};

	class BASE_LIBRARY_EXPORT t_gcombMIX : virtual public t_gcombIF
	{
	public:
		t_gcombMIX(t_gsetbase *setting, t_glog* log, shared_ptr<t_gbiasmodel> bias_model, t_gallproc* data);
		~t_gcombMIX();
		bool cmb_equ(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, t_gbaseEquation& result) override;

	protected:
		OBSCOMBIN  _observ;
		t_gionex*  _gdata_ionex = nullptr;
		bool _cmb_equ_IF_P1(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, t_gbaseEquation& result);
	};

	
}
#endif