/**
*
* @verbatim
	History
	 -1.0 Hongmin Zhang  2019-12-31  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gslrmodel.h
* @brief	    model of constructing B,P,L for slr observations
* @author       Hongmin Zhang, Wuhan University
* @version		1.0.0
* @date		    2019-12-31
*
*/

#ifndef GSLRMODEL_H
#define	GSLRMODEL_H

#include "gmodels/gpppmodel.h"
#include "gmodels/gtideIERS.h"
#include "gproc/glsqmatrix.h"
#include "gproc/glsq.h"
#include "gall/galllra.h"
#include "gall/gallproc.h"
#include "gall/gallorb.h"
#include "gall/gallattitude.h"

#include "gall/gallpreproc.h"
#include "gdata/gslrf.h"
#include "gdata/gorb.h"
#include "gdata/gpoleut1.h"
#include "gutils/gtrs2crs.h"
#include "gutils/gscf2crsleo.h"
#include "gdata/gnavde.h"
#include "gdata/gifcb.h"
#include "gdata/gion.h"

//#define DEBUG_LOG

namespace great
{
	class BASE_LIBRARY_EXPORT t_gslrmodel : public t_gpppmodel
	{
	public:

		t_gslrmodel(t_gallproc* data, t_glog* log, t_gsetbase* setting);

		virtual ~t_gslrmodel();

		void setOTL(t_gallotl* otl);

		bool apply_erp(const t_gtime & epoch,t_gallpar& pars, string site);
		bool _cmb_equ(const t_gtime & epoch, t_gallpar& pars, string sat_leo, t_glsqEquationMatrix& equ, string slr_site);
		//void setOBJ(t_gallobj* gobj);
		bool _prt_obs(const t_gtime & epoch, t_gallpar& pars, vector<pair<int, double>>& coeff, string slr_site,string sat);

		bool _omc_obs(const t_gtime & epoch,  double& omc);

		void _calculate_sta_pos(t_gtime time, ColumnVector* x_f, ColumnVector* v_f, ColumnVector* x_j, ColumnVector* v_j);

		void _rotate_posvel(ColumnVector& x_f, ColumnVector& v_f, ColumnVector& x_j, ColumnVector* v_j, Matrix& rotmat, Matrix& drotmat);

		void _machine_offset(ColumnVector& xyz);

		bool _apply_rec_tides(const t_gtime& epoch, t_gtriple& rec);

		void get_lracor(string sat,double zeni0, double azim0, double* dlra);

		void _update_rot_matrix(const t_gtime& epoch);

		void troposphere_delay_slr(double elev, double* delay, double* partial, t_gtime now,int ista);

		void trop_slr_fculzd_hpa(t_gtime now, int ista, double elev, double* FCUL_ZTD, double* FCUL_ZHD, double* FCUL_ZWD, double WVP);

		double trop_slr_fcul_a(double elev, t_gtime now, int ista);

		double _Partial_Orb_Leo(t_gpar& par, t_gtime epoch);

		double _Partial_Erp(t_gpar& par);

		double _Partial(const t_gtime& epoch, t_gpar& par, string slr_site, string slr_sat);

		void enu2trs(double lat, double lon, Matrix& rotmat);

		t_gslrf*        _allslrf = nullptr;

		t_gallpreproc*      _allpreproc = nullptr;

		t_gallrb*      _allrb = nullptr;


	protected:
		t_gtideIERS* _tide = nullptr;			///< tide correction model

	    
		t_galllra*      _alllra = nullptr;
		
		t_gallorb*      _allorb = nullptr;		  ///< all orb data
		t_gpoleut1*     _allerp = nullptr;  	  ///< all poleut1 data
		t_gnavde *      _allde = nullptr;	 	  ///< all panetnav info
		t_gallatt*      _allatt = nullptr;	      ///< all attitude data for leo sat

		t_gtrs2crs*    _trs2crs_2000 = nullptr;	  ///< trs2crs matrix
		t_gscf2crsleo* _scf2crs_leo = nullptr;	  ///< ant2crs matrix for leo sat

		string				_crt_rec;
		string				_crt_sat;
		
		string             _sat_leo;

		double             _two_way;

		int _siteIndex;

		Matrix _partialSite;

		Matrix _siteRotmat, _drdxpole, _drdypole, _drdut1, _drderp;
		t_gtriple _crs_rec_pco;
		t_gtriple _crs_sta_crd;
		t_gtriple _trs_rec_crd;
		t_gtriple _crs_rec_crd;
		double    _slr_rho;
		
	};
}

#endif