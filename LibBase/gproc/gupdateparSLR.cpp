#include "gupdateparSLR.h"

namespace great
{

	t_gupdateparSLR::t_gupdateparSLR(t_gallproc* data)
	{
		_alldata = data;
	}

	t_gupdateparSLR::~t_gupdateparSLR()
	{
	}

	bool t_gupdateparSLR::_update_erp_par(const t_gtime & epoch, int update_id, t_gallpar & allpars, t_gupdateparinfo& update_info)
	{
		t_gpar par_tmp = allpars[update_id];
		double intv = par_tmp.end - par_tmp.beg;
		t_gtime epoch_beg;
		epoch_beg.from_mjd(epoch.mjd(), 0, 0.0);
		par_tmp.setTime(epoch_beg, epoch_beg + intv);

		if (!_allerp) _allerp = dynamic_cast<t_gpoleut1*>((*_alldata)[t_gdata::ALLPOLEUT1]);
		if (!_allerp) return false;
		_allerp->getEopData(epoch_beg, data1);
		_allerp->getEopData(epoch_beg+intv, data2);
		
		switch (par_tmp.parType)
		{
		case t_gpar::XPOLE:
		{
			par_tmp.value(data1["XPOLE"]);
			break;
		}
		case t_gpar::YPOLE:
		{
			par_tmp.value(data1["YPOLE"]);
			break;
		}
		case t_gpar::DXPOLE:
		{
			par_tmp.value(data2["XPOLE"] - data1["XPOLE"]);
			break;
		}
		case t_gpar::DYPOLE:
		{
			par_tmp.value(data2["YPOLE"] - data1["YPOLE"]);
			break;
		}
		case t_gpar::UT1:
		{
			par_tmp.value(data1["UT1-TAI"]);
			break;
		}
		case t_gpar::DUT1:
		{
			par_tmp.value(data2["UT1-TAI"] - data1["UT1-TAI"]);
			//par_tmp.value(data1["LOD"]);
                        break;
		}

		}

		update_info.add(update_id + 1, par_tmp);

	}

	t_gupdateparinfo t_gupdateparSLR::get_all_update_parameters(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata)
	{
		t_gupdateparinfo update_info;
		
		// update Other pars
		_update_process_pars(epoch, allpars, update_info);

		return update_info;
	}

	void t_gupdateparSLR::_update_process_pars(const t_gtime& epoch, t_gallpar& allpars, t_gupdateparinfo & update_info)
	{

		// update process pars
		vector<t_gpar> parnew;
		for (unsigned int ipar = 0; ipar < allpars.parNumber(); ipar++)
		{
			// according to time,lremove,update_info
			if (allpars[ipar].end < epoch && allpars[ipar].lremove && !update_info.exist(ipar + 1))
			{

				// judge state mode
				if (_state_mode.find(allpars[ipar].parType) == _state_mode.end())
				{
					_update_erp_par(epoch, ipar, allpars, update_info);
				}
				else {
					_update_state_par(epoch, ipar, allpars, update_info);
				}
			}
		}
	}

	void t_gupdateparSLR::_update_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo & update_info)
	{

	}

}
