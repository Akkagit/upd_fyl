/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2019-12-23  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gupdatepar.cpp
* @brief		define update parameter base model
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-12-23
*
*/


#include "gupdatepar.h"

namespace great
{

	t_gupdateparinfo::t_gupdateparinfo()
	{
	}
	t_gupdateparinfo::~t_gupdateparinfo()
	{
	}

	bool t_gupdateparinfo::exist(int id)
	{
		return _remove_id.count(id)!=0 ? true : false;
	}

	void t_gupdateparinfo::add(int id)
	{
		if (_remove_id.count(id) != 0) 
		{
			throw std::logic_error("Remove id have repeated!");
		}

		_remove_id.insert(id);

	}

	void t_gupdateparinfo::add(t_gpar newpar)
	{
		_new_parlist.push_back(newpar);
	}

	void t_gupdateparinfo::add(int id, t_gpar newpar)
	{
		this->add(id);
		_new_parlist.push_back(newpar);
	}

	void t_gupdateparinfo::add(vector<pair<int, t_gpar> > update_par, t_glsqEquationMatrix state_equ)
	{
		for (auto iter : update_par) {
			//this->add(iter.first, iter.second);
			this->add(iter.first);
			_equ_parlist.push_back(iter.second);
		}
		_state_equ.add_equ(state_equ);
	}

	void t_gupdateparinfo::get(vector<int>& remove_id)
	{
		remove_id = vector<int>(_remove_id.begin(),_remove_id.end());
	}

	void t_gupdateparinfo::get(vector<t_gpar>& newparlist)
	{
		newparlist = _new_parlist;
	}

	void t_gupdateparinfo::get(vector<t_gpar>& update_newpar, t_glsqEquationMatrix& equ)
	{
		update_newpar = _equ_parlist;
		equ = _state_equ;
	}

	void t_gupdateparinfo::get(vector<int>& remove_id, vector<t_gpar>& newparlist, vector<t_gpar>& equ_parlist, t_glsqEquationMatrix & equ)
	{
		get(remove_id);
		get(newparlist);
		get(equ_parlist, equ);
	}



	t_gupdatepar::t_gupdatepar()
	{
	}

	t_gupdatepar::~t_gupdatepar()
	{
	}

	t_gupdatepar::t_gupdatepar(const t_gupdatepar & Other):
		_state_mode(Other._state_mode),
		_cycleslip(Other._cycleslip)
	{
		
	}

	void t_gupdatepar::set_interval(double interval)
	{
		_intv = interval;
	}

	void t_gupdatepar::set_sig_ion(double sigion)
	{
		_sig_ion = sigion;
	}

	void t_gupdatepar::set_cycleslip(shared_ptr<t_gcycleslip> cycleslip)
	{
		_cycleslip = cycleslip;
	}
	
	void t_gupdatepar::set_amb_update_way(bool way)
	{
	    _lite_update_amb = way;
	}

	void t_gupdatepar::set_par_state_mode(t_gpar::t_type type, int order, double dt, double noise)
	{
		_state_mode[type] = t_statemode(order, dt, noise);
	}

	void t_gupdatepar::update_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info)
	{
		this->_update_amb_pars(epoch, allpars, obsdata, update_info);
	}

	t_gupdateparinfo t_gupdatepar::get_all_update_parameters(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata)
	{
		t_gupdateparinfo update_info;
		// update AMB
		_update_amb_pars(epoch,allpars,obsdata, update_info);

		// update Other pars
		_update_process_pars(epoch, allpars, update_info);

		return update_info;
	}

	void t_gupdatepar::_update_process_pars(const t_gtime& epoch, t_gallpar& allpars,t_gupdateparinfo & update_info)
	{

		// update process pars
		vector<t_gpar> parnew;
		for (unsigned int ipar = 0; ipar < allpars.parNumber(); ipar++)
		{
			// according to time,lremove,update_info
			if (allpars[ipar].end < epoch && allpars[ipar].lremove && !update_info.exist(ipar+1))
			{
				
				// judge state mode
				if (_state_mode.find(allpars[ipar].parType) == _state_mode.end())
				{
					
					_update_process_par(epoch, ipar, allpars, update_info);
				}
				else {
					_update_state_par(epoch, ipar, allpars, update_info);
				}
			}
		}
	}

	void t_gupdatepar::_update_state_par(const t_gtime & epoch, int update_id, t_gallpar & allpars, t_gupdateparinfo& update_info)
	{
		// record the new add par location
		vector<pair<int,t_gpar> > parnewlist;
		t_glsqEquationMatrix state_equ;
		t_glsqEquationMatrix state_equ_tmp;
		vector<t_gpar> equ_parlist;

		update_info.get(equ_parlist, state_equ_tmp);
		int newpar_location = allpars.parNumber() + 1 + equ_parlist.size();
		t_gpar::t_type par_type = allpars[update_id].parType;

		for (int i = 0; i < _state_mode[par_type].order; i++) {
			t_gpar par_temp = allpars[update_id + i];
			double intv = par_temp.end - par_temp.beg;
			par_temp.setTime(epoch, epoch + intv);
			parnewlist.push_back(make_pair(update_id+i+1,par_temp));
		}

		// create the state_equation
		for (int Row = 1; Row <= _state_mode[par_type].order; Row++) {
			vector<pair<int, double> > B;
			double P;
			// for B
			//B.push_back(make_pair(update_id + Row - 1, 1.0));
			B.push_back(make_pair(update_id + Row, 1.0));
			for (int Col = Row; Col <= _state_mode[par_type].order; Col++) {
				B.push_back(make_pair(newpar_location + Col - 1, -1.0));
			}

			// for P
			P = _state_mode[par_type].P(Row, Row);

			t_gobscombtype type;
			state_equ.add_equ(B, P, 0.0, "", "", type, false);
		}

		update_info.add(parnewlist, state_equ);

	}

	void t_gupdatepar::_update_process_par(const t_gtime & epoch, int update_id, t_gallpar & allpars, t_gupdateparinfo& update_info)
	{
		t_gpar par_tmp = allpars[update_id];
		double intv = par_tmp.end - par_tmp.beg;
		par_tmp.setTime(epoch, epoch + intv);
		update_info.add(update_id + 1, par_tmp);
	}

	void t_gupdatepar::_udpate_gps_rec_ifb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo & update_info)
	{
		// for debug
		double my_mjd = epoch.mjd();
		double my_sod = epoch.sod();


		std::set<string> mapRec;
		for (auto& data : obsdata)
		{
			auto  rec = data.site();
			if (mapRec.count(rec) > 0) continue;
			mapRec.insert(rec);
			// add new ifb par
			int idx_gps_ifb = allpars.getParam(rec, t_gpar::GPS_REC_IFB_C3, "");
			if (idx_gps_ifb < 0)
			{
				int idx_max = allpars[allpars.parNumber() - 1].index;
				t_gpar par_gps_rec_ifb = t_gpar(rec, t_gpar::GPS_REC_IFB_C3, idx_max + 1, "");
				par_gps_rec_ifb.value(0.0);
				par_gps_rec_ifb.setTime(epoch, epoch);
				par_gps_rec_ifb.apriori(9000);
				update_info.add(par_gps_rec_ifb);
			}
		}

		// remove old ifb par
		// jdhuang : change mapPRN to mapRec
		for (unsigned int i = 0; i < allpars.parNumber(); i++)
		{
			if (allpars[i].parType == t_gpar::GPS_REC_IFB_C3 && mapRec.count(allpars[i].site) == 0)
			{
				update_info.add(i + 1);
			}
		}
	}



}
