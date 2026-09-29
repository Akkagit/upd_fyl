/**
*
* @verbatim
	History
	 -1.0 Wei Zhang     2019-01-01 creat the file.
	 -1.1 Jie Li		2019-03-31 add doxygen style code
	 -1.2 Jiande Huang  2020-10-10 Optimize code efficiency and structure
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			    gtrs2crs.h
* @brief		    calculate the rotation matrix from TRS to CRS and the corresponding partial derivation matrix
*
* @author		    Wei Zhang, Wuhan University
* @version		    1.0.0
* @date			    2019-01-01
*
*/

#ifndef GTRS2CRS_H
#define GTRS2CRS_H

#include <vector>

#include "gutils/gtime.h"
#include "gdata/gpoleut1.h"
#include "gutils/gconst.h"
#include "gutils/gtriple.h"
#include "gmodels/ginterp.h"
#include "gdata/gleapsecond.h"
#include "gutils/gmatrixconv.h"

#include "gutils/ExportBase.h"
#include "newmat/newmat.h"
#include "newmat/newmatio.h"

using namespace std;

namespace great
{
	/**
	*@brief	   Class for storaging pole and ut1 data for one day
	*/
	class t_pudaily
	{
	public:
		/** @brief default constructor. */
		t_pudaily();

		/** @brief default destructor. */
		virtual ~t_pudaily() {};

		t_gtime  time;    ///< epoch time
		double   xpole;   ///< x pole
		double   ypole;   ///< y pole
		double   ut1;     ///< dut1

	};

	//todo 命名
	/**
	 *@brief	   Class for nutation for one epoch
	 */
	class t_EpochNutation
	{
	public:
		/** @brief default constructor. */
		t_EpochNutation();

		/** @brief default destructor. */
		virtual ~t_EpochNutation() {};

		double T;		///< epoch
		double Pi;		///< nutation angle of the longitude(radian)
		double Ep;		///< nutation angle of the obliquity(radian)
	};
	/**
	 *@brief	   Class for zonaltide data
	 */
	class t_zonaltide
	{
	public:
		/** @brief default constructor. */
		t_zonaltide();

		/** @brief default destructor. */
		virtual ~t_zonaltide() {};

		t_gtime  time;    ///< epoch time
		double   ut1;     ///< dut1
		double   lod;     ///< dlod 
		double   omega;
	};

	/**
	*@brief	   Class for calculating the rotation matrix from TRS to CRS and the corresponding partial derivation matrix
	*/
	class BASE_LIBRARY_EXPORT t_gtrs2crs
	{
	public:
		/** @brief default constructor. */
		t_gtrs2crs();

		/** @brief constructor. */
		t_gtrs2crs(bool erptab, t_gpoleut1* poleut1);

		t_gtrs2crs& operator=(const t_gtrs2crs& Other);

		/** @brief default destructor. */
		virtual ~t_gtrs2crs() {};

		/** @brief calculte process matrix. */
		void calcProcMat(const bool& partial, const int& axis,const double& angle, vector<Matrix>& rot);
		/** @brief calculate the rotation matrix. */
		void calcRotMat(const t_gtime& epoch, const bool& ldxdpole, const bool& ldydpole, const bool& ldudpole);

		//todo const??
		/** @brief return rotation matrix. */
		Matrix& getRotMat();

		/** @brief return du matrix. */
		Matrix& getMatDu();

		/** @brief return dx matrix. */
		Matrix& getMatDx();

		/** @brief return dy matrix. */
		Matrix& getMatDy();

		/** @brief return xpole. */
		double getXpole();
		/** @brief return ypole. */
		double getYpole();
		/** @brief return gmst. */
		double getGmst();
		/** @brief return now epoch. */
		t_gtime getCurtEpoch();

	protected:
		/** @brief interpolate short term tidal corrections(fortran: polut1_ray_interpolation). */
		void	 _tide_corrections(t_gtime& t, t_gtriple& xyu);

		/** @brief calculate short term tidal corrections(fortran: polut1_ray_calc). */
		t_pudaily& _tideCor1Cal(t_pudaily& b);

		/** @brief calculate the diurnal and semi-diurnal variations in Earth Orientation Parameters (x,y, UT1) from ocean tides(fortran: ORTHO_EOP). */
		void   _ORTHO_EOP(t_gtime& t, t_gtriple& eop);

		/** @brief calculate the partials of the tidal variations to the orthoweights(fortran: CNMTX). */
		void   _CNMTX(t_gtime& t, double * h);

		/** @brief evaluate the model of polar motion for a nonrigid Earth due to tidal gravitation. */
		void   _PMSDNUT2(t_gtime& t, double* pm);

		/** @brief evaluate the model of subdiurnal libration in the axial component of rotation, expressed by UT1 and LOD.. */
		void   _UTLIBR(t_gtime& t, double* temp);

		/** @brief calculate zonaltide correction on ut1(fortran: ut1_zonaltide_interpolation). */
		double     _tideCor2(const double& dRmjd);

		 /** @brief evaluate the effects of zonal Earth tides on the rotation of the Earth. */
		void      _RG_ZONT2(const double& dT, double *DUT, double *DLOD, double *DOMEGA);

		/** @brief calculate fund arg. */
		void      _FUNDARG(const double& T, double *L, double *LP, double *F, double *D, double *OM);

		/** @brief calculate _pudata(fortran: polut1_interpolation). */
		void      _calPoleut1(t_gtime& t, double *x, t_gpoleut1* poleut1, string type = "");

		/** @brief interpolate the nutation(fortran: nutation_interpolation). */
		void   _nutInt(const double& dRmjd, double* dpsi, double* deps, const double& step);

	   /** @brief calculate nutation according to IAU 2000A model(fortran: nu2000a). */
		void   _nutCal(const double& dATE1, const double& dATE2, double *dPSI, double *dEPS);

		/** @brief interpolation. */
		double    _interpolation(const int& iOrder, const int& iPoint, double *pdX, double *pdY, const double& dXin);

		/** @brief calculate rotation matrix of Precession and nutation. */
		void      _process2000(const double& dRmjd, const double& dpsi, const double& deps);

	   /** @brief calculate the angle (unit: radian) according to the current epoch and J2000(fortran: sp2000). */
		double    _sp2000(const double& dDATE1, const double& dDATE2);

		/** @brief calculate the earth rotation angle (unit: radian)(fortran: era2000). */
		double    _era2000(const double& dJ1, const double& dJ2);

	   /** @brief calculate Greenwich sidereal time, IAU 2000. */
		double    _gst2000(const double& dRmjd, const double& era, const double& dpsi);

	   /** @brief fmod to (0,2pi). */
		double    _iau_anp(const double& dA);

		/** @brief Equation of the equinoxes complementary terms,
		*	consistent with IAU 2000 resolutions(return complementary terms (radians)). */
		double    _eect2000(const double& dRmjd);

	   /** @brief fmod. */
		double	 _iau_anpm(const double& dA);

		//void rotation_matrix(bool ldot, int iaxis, double angle, double *protmat, double *pdrotmat);  //generate rotation matrix and its derivation given the rotating angle and axis(bool ldot : true - with derivation false - without derivation)

		// inner variable
		bool      _erptab;    ///< if ERP table is used or not
		t_gtime   _epo;       ///< cureent epoch
		t_gtriple _pudata;    ///< pole and ut1 data after interpolation(xpole, ypole and ut1-tai)
		t_gtime   _tdt;		  ///< dynamic time
		double    _arg[5];    ///< fundamental arguements of nutation in radius
		Matrix    _qmat;      //
		double    _epsa;      //

		// out variable
		double    _xpole;     ///< x pole of pole shift
		double    _ypole;     ///< y pole of pole shift
		double    _gmst;      ///< Greenwich Mean Sidereal Time
		Matrix    _rotmat;    ///< rotation matrix from TRS to CRS
		Matrix    _rotdu;     ///< partial of rotmat wrt to ut1.
		Matrix    _rotdx;     ///< partial of rotmat wrt to xpole
		Matrix    _rotdy;     ///< partial of rotmat wrt to ypole

		t_gpoleut1*		  _poleut1;			///< poleut1 data

		static t_pudaily   _tbp[2];  //
		static t_zonaltide _tbz[3];  //

		// constant values
		const double _rr = 1.00273781191135448e0;   ///< r = gmst/ut1 = 1.00273781191135448D0 (ERA2000.f90, IERS20xx)
		const double _djc = 36525.0;				///< days per julian century


		// 在类函数中使用的各种常数变量

		// === void t_gtrs2crs::_process2000(const double& dRmjd, const double& dpsi, const double& deps) ===				
		const double dEps0 = 84381.448 / (RAD2SEC);    //! J2000 obliquity (Lieske et al. 1977)
		const double dRa0 = -0.0146 / (RAD2SEC);       //! The ICRS RA of the J2000 equinox (Chapront et al., 2002)
		const double dPrecor = -0.29965 / (RAD2SEC);   //! The precession and obliquity corrections (radians per century) (page 43, equ(27))
		const double dOblcor = -0.02524 / (RAD2SEC);
		const double dPsibi = -0.041775 / (RAD2SEC);   //! The frame bias corrections in longitude and obliquity (page 43, equ(28) ?)
		const double dEpsbi = -0.0068192 / (RAD2SEC);
		
		// === double t_gtrs2crs::_era2000(const double& dJ1, const double& dJ2) ===
		const double dPI = 3.141592653589793238462643;
		const double d2PI = 6.283185307179586476925287;  //DOUBLE PRECISION DJ0
		const double dJ0 = 51544.5;                      //Reference epoch (J2000), JD

		// === double t_gtrs2crs::_eect2000(const double& dRmjd) ===		
		const double dAS2R = 4.848136811095359935899141 * 0.000001;  //    	 *  Arcseconds to radians
		const double dJC = 36525;                                    //    	 *  Days per Julian century

		// == void  t_gtrs2crs::_nutCal(const double& dATE1, const double& dATE2, double *dPSI, double *dEPS) ===
		const double dDAS2R = 4.848136811095359935899141e-6;
		const double dDMAS2R = dDAS2R / 1e3;
		const double dTURNAS = 1296000e0;
		const double dU2R = dDAS2R / 1e7;
		const double dDJ0 = 51544.5e0;
		const double dDJC = 36525e0;

		// === void t_gtrs2crs::_PMSDNUT2(t_gtime& t,double* pm) === 
		const double _sec2rad = 4.848136811095359935899141e-6;
		const double _turnas = 1296000.0;
		const double rmjd0 = 51544.5;
		const double rad2sec = 86400.0 / d2PI;

	};
} // namespace

#endif