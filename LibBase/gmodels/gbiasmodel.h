/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2020-05-15 created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    t_gbiasemodel.h
* @brief	    mainly about how to cacultae B P l single
* @mainpage     C++ version 
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-05-18
*
*/

#ifndef GBIASEMODEL_H
#define GBIASEMODEL_H

#include "gutils/gtime.h"
#include "gall/gallpar.h"
#include "gdata/gsatdata.h"
#include "gdata/gobsgnss.h"
#include "gmodels/gprecisemodel.h"
#include "gmodels/gpceprecisemodel.h"
#include "gmodels/gpppleomodel.h"
#include "gmodels/gbasemodel.h"

using namespace std;
using namespace gnut;

namespace great
{
    class BASE_LIBRARY_EXPORT t_gbiasmodel
    {
    public:
        t_gbiasmodel();
        virtual ~t_gbiasmodel();

        virtual bool cmb_equ(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,t_gobs& gobs,t_gbaseEquation& result)=0;
        virtual void update_obj_clk(const string& obj, const t_gtime& epo, double clk) = 0;
        virtual double get_rec_clk(const string& obj) = 0; //add by xiongyun

        virtual t_gmodel* precisemodel() { return nullptr; }
    };

    class BASE_LIBRARY_EXPORT t_gprecisebias:public t_gbiasmodel
    {
    public:
        t_gprecisebias(t_gallproc *data, t_glog *log, t_gsetbase *setting);
        ~t_gprecisebias();

		void set_multi_thread(const set<string>& sites);

        bool cmb_equ(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,t_gobs& gobs,t_gbaseEquation& result) override;
        void update_obj_clk(const string& obj, const t_gtime& epo, double clk) override;
        double get_rec_clk(const string& obj) override;

        t_gmodel* precisemodel() { return &gmodel; };
    protected:
        t_gprecisemodel gmodel;
		tuple<string, string, t_gtime> _rec_sat_before;

		bool _multi_thread_flag;
		map<string, shared_ptr<t_gprecisemodel> > _map_site_model;
		map < string, tuple<string, string, t_gtime> > _map_flag;
		t_gmutex _map_mtx;

    private:
        t_glog *_log;
        t_gsetbase* _setting;//test by lvhb
    };

    class BASE_LIBRARY_EXPORT t_gsppbias:public t_gbiasmodel
    {

    };

    class BASE_LIBRARY_EXPORT t_gIRCbias:public t_gbiasmodel
    {
    public:
        t_gIRCbias(t_gallproc* data, t_glog* log, t_gsetbase* setting);
        ~t_gIRCbias();

        bool cmb_equ(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, t_gobs& gobs, t_gbaseEquation& result) override;
        void update_obj_clk(const string& obj, const t_gtime& epo, double clk) override;
        double get_rec_clk(const string& obj) override;

    protected:

        t_gprecisemodel gmodel;

    private:
        t_glog* _log;
    };


	class BASE_LIBRARY_EXPORT t_simuleobias :public t_gbiasmodel
	{
	public:
		t_simuleobias(t_gallproc* data, t_glog* log, t_gsetbase* setting);
		~t_simuleobias();

		bool cmb_equ(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, t_gobs& gobs, t_gbaseEquation& result) override;
		void update_obj_clk(const string& obj, const t_gtime& epo, double clk) override;
		double get_rec_clk(const string& obj) override;

		t_gmodel* precisemodel() { return &gmodel; };

	protected:

		t_gpppleomodel gmodel;
		pair<string, string> _rec_sat_before;
		double _noiseP_crt;
		double _noiseL_crt;

	private:
		t_glog* _log;
		t_gallbias* _allbias;
	};

} // namespace great

#endif