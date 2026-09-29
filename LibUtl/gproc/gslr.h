/**
*
* @verbatim
	History
	 -1.0 Hongmin Zhang		2019-10-25  creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gslr.h
* @brief		Satellite laser ranging
*
* @author       Hongmin Zhang, Wuhan University
* @version		1.0.0
* @date		    2019-10-25
*
*/

#ifndef G_SLR_H
#define G_SLR_H

#include "gio/glog.h"
#include "gio/giof.h"

#include "gall/gallics_nav.h"
#include "gall/gallprec.h"
#include "gall/gallorb.h"
#include "gall/gallobj.h"
#include "gall/gallpreproc.h"
#include "gdata/gattitude.h"
#include "gdata/gresult.h"
#include "gall/gallsta.h"

#include "gdata/gics_nav.h"
#include "gdata/gorb.h"
#include "gdata/gpoleut1.h"
#include "gdata/gslrf.h"
#include "gdata/gnavde.h"
#include "gall/galllra.h"
#include "gall/gallnpt.h"
#include "gall/gallotl.h"
#include "gall/gallresult.h"

#include "gproc/glsq.h"

#include "newmat/newmat.h"
#include "gutils/ExportUtl.h"
#include "gutils/gtrs2crs.h"
#include "gutils/gscf2crsleo.h"
#include "gutils/gcrs2acr.h"
#include "gutils/gscf2crs.h"

#include "gdata/gpreprocslr.h"

#include "gmodels/gtideIERS.h"

using namespace std;

namespace great
{
	/**
	*@brief	   Class for slr validation
	*/
	class UTL_LIBRARY_EXPORT t_gslr
	{
	public:
		/** @brief constructor. */
		t_gslr();

		/**
		* @brief  constructor.
		* @param[in]	set    setting information
		* @param[in]	log    log file
		*/
		t_gslr(t_gsetbase* set, t_glog* log);

		void setDataUsed(string sat);

		/** @brief destructor. */
		virtual ~t_gslr();

		bool processBatch();

		void slrmod(double e, double n, double u);

		bool _init_data();

		void setNPT(t_gallnpt* gnpt);

		void setLRA(t_galllra* glra);

		void setTime(t_gtime now);

		void setOBJ(t_gallobj* gobj);

		void setOTL(t_gallotl* otl);

		void setSLRF(t_gslrf* gslrf);

		void setDE(t_gnavde* gde);

		void setORB(t_gallorb* gorb);

		void setATT(t_gattitude*    attitude);

		void setPoleut(t_gpoleut1* poleut1);

		void enu2trs(double lat, double lon, Matrix& rotmat);

		void get_lracor(double zeni0, double azim0, double* dlra);

		void troposphere_delay_slr(double elev, double* delay, double* partial,int i,int ista);

		void trop_slr_fculzd_hpa(int i, int ista, double elev, double* FCUL_ZTD, double* FCUL_ZHD, double* FCUL_ZWD,double WVP);
		//void _update_rot_matrix(const t_gtime& epoch);
		double trop_slr_fcul_a(double elev, int i, int ista);

		bool _apply_rec_tides(const t_gtime& epoch, t_gtriple& rec);

		bool _get_setting();

		void _set_out();

		void setifLEO(string str);

		void result_reserve(t_gresult& result);

		void setResult(t_gresult& result);

		void setPreObs(t_gallpreproc& preobs);

		void setRB(t_gallrb* rb);

		void setSTA(t_gallsta* sta);

		void preresult_reserve(t_gallpreproc& preresult);

		void _update_rot_matrix(const t_gtime& epoch);

		void _rotate_posvel(ColumnVector& x_f, ColumnVector& v_f, ColumnVector& x_j, ColumnVector* v_j, Matrix& rotmat, Matrix& drotmat);

		void _calculate_sta_pos(t_gtime time, ColumnVector* x_f, ColumnVector* v_f, ColumnVector* x_j, ColumnVector* v_j);

	protected:
		// For log
		t_glog*				_glog;             ///< Genereal   log output. 

		t_gtide*            _tide = nullptr;

		t_gpreprocslr *     _gpreproc;

		t_gallpreproc*      _gpreOBS;

		t_gallrb*         _gallrb;

		t_gallsta*       _gallsta;

		t_gtime            _time;
		// For set
		t_gsetbase*			_set;              ///< Base setting. 

		string              _sat;             ///< LEO

		double              _slre;            ///< the term S : the vector from the center of mass to the reference point

		double              _slrn;            ///< slre and slrn are always negative

		double              _slru;

		t_gallnpt*          _gnpt;

		t_galllra*          _glra;

		string              _ifLEO;

		t_gslrf*            _gslrf;

		t_gnavde*           _gde;

		t_gallorb*          _gorb;

		t_gtrs2crs*  _trs2crs_cord;              ///< class of TRS to CRS

		t_gscf2crsleo*	 _scf2crsleo_cord;       ///< class of SCF to CRS 

		t_gscf2crs*      _scf2crs_cord;

		t_gcrs2acr*     _crs2acr;

		t_gtrs2crs*  _trs2crs_time;              ///< time of TRS to CRS

		t_gtime             _beg;

		t_gtime             _end;

		int                 _sampling;

		t_giof*				_slr;

		t_gtrs2crs*    _trs2crs_2000 = nullptr;

		t_gpoleut1*     _allerp = nullptr;

		shared_ptr<t_gobj>	_crt_obj;

		t_gallobj*          _gallobj;

		t_gpoleut1*         _gpoleut1;         ///< poleut1 data

		t_gresult_one*         _result;

		t_gresult*            _gresult;

		t_gattitude*        _gattitude;        ///< LEO attitude data

	};
}

#endif