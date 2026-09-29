
#ifndef _GORBLSQ_H
#define _GROBLSQ_H

#include "gutils/ExportLSQ.h"

#include "gproc/gclklsq.h"
#include "gdata/gorb.h"
#include "gdata/gpoleut1.h"
#include "gall/gallics_nav.h"
#include "gall/gallics_leo.h"
#include "gdata/gleapsecond.h"
#include "gmodels/glsqprocIF.h"


using namespace std;

namespace great
{
	class LSQ_LIBRARY_EXPORT t_gorblsq_com : public t_glsqprocIF
	{
	public:
		// construct
		t_gorblsq_com(t_gsetbase* base, t_glog* obj);
		virtual ~t_gorblsq_com();

	public:

		virtual bool GenerateProduct() override;

	protected:

		virtual bool _InitLsq(t_gallproc* data) ;

		// ref clk
		string          _site_ref;
		string          _sat_ref;
		double			_sig_refclk;
		CONSTRPAR		_crd_est;

		// recclk file (out)
		//t_giof*	_result;              /// output file of result
		string  _recclk_file;
		string  _leoclk_file;
		string  _satclk_file;

		//// needed data
		t_gpoleut1*    _poledata;
		t_gallics_nav* _icsdata;
		t_gallics_leo* _icsleo;        ///< ics data of the leo

		set<string> _satleo;          ///< long name of leo list to be processed

		// set recclk file
		void _setOut() ;

		// init orb parameter
		bool _init_EOP_parameter();
		bool _init_SAT_parameter();
		/** @brief init LEO orbit parameters. */
		void _init_ORB_parameter();
		// write ics file
		void _update_icsfile();
		/** @brief upsate LEO ics information. */
		void _update_icsleo();
	};
}

#endif /* _GORBLSQ_H */