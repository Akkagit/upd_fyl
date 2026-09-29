// history:
// jdhuang 2019-05-13 add the ud uc case;

#include "gmodels/glsqprocALL.h"
#include "gutils/ginfolog.h"
#include "gdata/gobsgnss.h"
#include "gset/gsetpar.h"
#include "gdata/gion.h"
#include "gcoders/ion.h"
#include "gio/gfile.h"
#include "gproc/gupdateparALL.h"
using namespace std;

namespace great
{
	t_glsqprocALL::t_glsqprocALL()
	{

	}

	t_glsqprocALL::t_glsqprocALL(t_gsetbase* set, t_gallproc* data, t_glog* log) :
		t_glsqproc(set, data, log)
	{
		// Nothing but cout a string here
		gtrace("t_glsqprocALL::t_glsqprocALL");

		// get ion model
		_ion_model  = dynamic_cast<t_gsetproc*>(set)->ion_model();

		// jdhuang: 3 freq
		if (_frequency == 3)
		{
			auto logfile_12 = dynamic_cast<t_gturboedit*>(_slip12.get())->logfile_exist();
			auto logfile_13 = dynamic_cast<t_gturboedit*>(_slip13.get())->logfile_exist();
			dynamic_cast<t_gturboedit*>(_slip12.get()  )->merge_logfile_exist(logfile_12);
			dynamic_cast<t_gturboedit*>(_slip13.get())->merge_logfile_exist(logfile_13);

			shared_ptr<t_gupdatepar> update_all(new t_gupdateparALL(_slip12, _slip13, _band_index, _iono_order));
			update_all->set_amb_update_way(/*_realtime*/_lite_turboedit);
			update_all->set_interval(_obs_intv);
			_lsq->set_update_par(update_all);
		}
		// 2 freq
		else
		{
			shared_ptr<t_gupdatepar> update_all(new t_gupdateparALL(_slip12, _band_index, _iono_order));
			update_all->set_amb_update_way(/*_realtime*/_lite_turboedit);
			update_all->set_interval(_obs_intv);
			_lsq->set_update_par(update_all);
		}

		shared_ptr<t_gbiasmodel> precise_bias(new t_gprecisebias(data,log,set));
		_base_model = shared_ptr<t_gbasemodel>(new t_gcombALL(set, log, precise_bias,data));

	
	}

	t_glsqprocALL::~t_glsqprocALL()
	{

	}


	//bool t_glsqprocALL::_add_sion_pars_value(const t_gtime& beg, const t_gtime& end, t_glsq* lsq, t_gion* gion)
	//{

	//	for (const auto& it_rec : _rec_list)
	//	{
	//		for (const auto& it_sat : _sat_list)
	//		{
	//			int idx = _lsq->_x_solve.getParam(it_rec, t_gpar::SION, it_sat);
	//			if (idx < 0) continue;
	//			else
	//			{
	//				double value = gion ? gion->get_value(it_rec, it_sat, beg, end) : 0.0;
	//				_lsq->_x_solve[idx].value(value);
	//			}
	//		}
	//	}
	//	return false;
	//}

	//void t_glsqprocALL::_extract_ion_file(t_gallrecover& recover_data)
	//{
	//	t_gion ion_data;
	//	vector<t_tuple_ion> ions;
	//	recover_data.get_iondata(ions);

	//	for (const auto& iter : ions)
	//	{
	//		t_gion_record tmp(std::get<0>(iter),
	//			std::get<1>(iter), std::get<2>(iter),
	//			std::get<3>(iter), std::get<4>(iter),
	//			std::get<5>(iter), std::get<6>(iter));

	//		ion_data.add_record(tmp);
	//	}

	//	// encoder data
	//	t_gfile gout;
	//	string ion_path = dynamic_cast<t_gsetout*>(_gset)->outputs("ion");

	//	if (ion_path.empty())
	//	{
	//		ion_path = "grt$(date).ion";
	//	}
	//	substitute(ion_path, "$(date)", _beg_time.str_yyyydoy(), false);

	//	gout.path(ion_path);
	//	gout.glog(_glog);

	//	t_ion coder(_gset);
	//	coder.add_data("ID1", &ion_data);
	//	gout.coder(&coder);

	//	gout.run_write();
	//	return;
	//}

}