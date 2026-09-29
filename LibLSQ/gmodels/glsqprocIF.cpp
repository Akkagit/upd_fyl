// history:
// jdhuang 2019-05-13 add the ud uc case;

#include "gmodels/glsqprocIF.h"
#include "gproc/gupdateparIF.h"
#include "gutils/ginfolog.h"

using namespace std;

namespace great
{
	t_glsqprocIF::t_glsqprocIF()
	{

	}

	t_glsqprocIF::t_glsqprocIF(t_gsetbase* set, t_gallproc* data, t_glog* log) :
		t_glsqproc(set, data, log)
	{
		// Nothing but cout a string here
		if (_frequency == 2)
		{
			shared_ptr<t_gupdatepar> updateIF(new t_gupdateparIF(_slip12, _band_index));
			updateIF->set_amb_update_way(/*_realtime*/_lite_turboedit);
			updateIF->set_interval(_obs_intv);
			_lsq->set_update_par(updateIF);

			shared_ptr<t_gbiasmodel> precise_bias(new t_gprecisebias(data, log, set));
			_base_model = shared_ptr<t_gbasemodel>(new t_gcombIF(set, log, precise_bias,data));
		}
		else if (_frequency == 3)
		{
			shared_ptr<t_gupdatepar> update_if13(new t_gupdateparIF13(_slip12, _slip13, _band_index));
			update_if13->set_amb_update_way(/*_realtime*/_lite_turboedit);
			update_if13->set_interval(_obs_intv);
			_lsq->set_update_par(update_if13);

			shared_ptr<t_gbiasmodel> precise_bias(new t_gprecisebias(data, log, set));
			_base_model = shared_ptr<t_gbasemodel>(new t_gcombIF13(set, log, precise_bias,data));
		}
		else
		{
			throw logic_error("you should use DOUBLE/TRIPLE in xml");
		}


		gtrace("t_glsqprocIF::t_glsqprocIF");
	}

	t_glsqprocIF::~t_glsqprocIF()
	{
		
	}


}
