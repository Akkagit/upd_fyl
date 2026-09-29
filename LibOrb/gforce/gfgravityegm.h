/**
 *
 * @verbatim
	 History
	  -1.0 Jie Li & Yujie Qin   2019-02-21 Creat the file.
	  -1.1 Jiande Huang         2019-03-14 Code comments are added in Doxygen format.
   @endverbatim
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		gfgravityegm.h
 * @brief		The class for EGM08 gravity model, including 12th-order EGM08, 120th-order EGM08 and so on.
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-02-21
 *
 */

#ifndef GFGRAVITYEGM_H
#define GFGRAVITYEGM_H

#include "gforce/gfgravity.h"
#include "gdata/gegmcoef.h"
#include "gutils/gtriple.h"
#include "gforce/gfarg_init.h"
#include "gall/gallotl.h"
#include "gdata/gotdata.h"
#include "gdata/gdeadata.h"
using namespace std;

namespace great
{
	/** 
	* @brief EGM 2008 gravity model, 
	*	the solid tide and ocean tide were considered in the Cnm and Snm coefficient.
	* 
	* Different orders of gravity field will be distinguished by different subclasses.
	* The ocean and solid tide are considered here.
	*/
	class t_gfgravity_egm08 : public t_gforce
	{
	public:
		/**
		 * @brief constructor.
		 * @remarks The max ocean degree is 35.
		 * @param [in]  arg		initial argument
		 * @param [in]  log		log file
		 */
		t_gfgravity_egm08(t_init_arg* arg, t_glog* log);
		
		/** @brief defult destructor. */
		virtual ~t_gfgravity_egm08() {};

		/** 
		 * @brief show the information of model.
		 */
		void show();

		/** 
		 * @brief check some important pars in the model.
		 * @return
			@retval true	data has no problem 
			@retval false   some wrong with data
		 */
		bool check();

		/** 
		* @brief calculate the acceleration and partial matrix e.
		* @param [in] argument		Parameters to be used in calculation.
		* @return
			@retval true	calulate sucessfully
			@retval false   calulate unsucessfully
		*/
		bool calculate(t_gforce_calarg* argument);


	protected:
		/**
		* @brief calculate the Cnm coefficents by using the EGM08 coefficents read from file
		* @remarks Specific formulas can be referred to Geodetic books.
		* @return
			@retval true	get Cnm sucessfully
			@retval false   get Cnm unsucessfully
		*/
		bool _isGetCnm();		

		/** 
		* @brief calculate the Snm coefficents by using the EGM08 coefficents read from file
		* @remarks Specific formulas can be referred to Geodetic books.
		* @return
			@retval true	get Snm sucessfully
			@retval false   get Snm unsucessfully
		*/
		bool _isGetSnm();		

		/** 
		* @brief calculate the Anm coefficents by using the EGM08 coefficents read from file
		* @remarks Specific formulas can be referred to Geodetic books.
		* @return
			@retval true	calulate Anm sucessfully
			@retval false   calulate Anm unsucessfully
		*/
		bool _cal_Anm_coeff();	

		/** 
		* @brief calculate gravity acceleration in J2000 by using the coefficients.
		* @remarks Specific formulas can be referred to Geodetic books.And the solid tide and
		*		ocean tide were considered.
		* @param [in] ej2		The rotation matrix from TRS to J2000.
		* @return
			@retval true	calulate acceleration sucessfully
			@retval false   calulate acceleration unsucessfully
		*/
		bool _cal_acc_mat(const Matrix& e2j);		

		/** 
		* @brief calculate gravity acceleration partial matrix in J2000 by using the coefficients.
		* @remarks Specific formulas can be referred to Geodetic books.And the solid tide and
		*		ocean tide were considered.
		* @param [in] ej2		The rotation matrix from TRS to J2000.
		* @return 
			@retval true	calculate amat sucessfully
			@retval false   calculate amat unsucessfully
		*/
		bool _cal_Amat_coeff(const Matrix& e2j);	

		/// get mean pole according to mjd (IERS2010).
		/**
		* @brief calculating and evaluating pole shift point location.Note that the default coordinate system is TRS
		* @param [in] mjd		Dynamic time.
		* @param [out] xpm		The x pos of pole point.
		* @param [out] ypm		The y pos of pole point.
		*/
		void _get_mean_pole(const double& mjd, double& xmp, double& ymp);		

		/// calculate dC and dS for solid tide. 
		/** 
		* @brief calculate solid tide coefficients.
		* @remarks Specific formulas can be referred to Geodetic books.
		* @param [in] tdt Dynamic time.
		* @param [in] gst Greenwich sidereal time.
		* @return
			@retval true	calculate solid tide sucessfully
			@retval false   calculate solid tide unsucessfully
		*/
		bool _calculate_solid_tide(const double& tdt, const double& gst);

		/// calculate dC and dS for ocean tide. 
		/** 
		* @brief calculate ocean tide coefficients.
		* @remarks Specific formulas can be referred to Geodetic books.
		* @param [in] tdt Dynamic time.
		* @param [in] gst Greenwich sidereal time.
		* @return
			@retval true	calculate ocean tide sucessfully
			@retval false   calculate ocean tide unsucessfully
		*/
		bool _calculate_ocean_tide(const double& tdt, const double& gst); 

		bool _calculate_ocean_pole_tide(const double& m1, const double& m2, const int& ndegree);
		/// calculate normalization coefficients for Lengendre functions. 
		/** 
		* @brief calculate normalization coefficient.
		* @param [in] mode Type of planetary ephemeris.
		* @param [in] ndegree The ocean degree n that we used.
		* @param [in] morder  The ocean degree m taht we used.
		* @param [in] ndim    Max ocean degree n.
		* @param [in] ndim    Max ocean degree m.
		* @param [out] norm   Normalize coefficients.
		* @return
			@retval true	calculate normalization coefficient sucessfully
			@retval false   calculate normalization coefficient unsucessfully
		*/
		bool _normalize_coeff(const string& mode, const int& ndegree, const int& morder, const int& ndim, const int& mdim, Matrix& norm);

		/// calculate the doodson angel 
		/** caculate beta angle.
		* @param [in] tdt Dynamic time.
		* @param [in] gst Greenwich sidereal time.
		* @param [out] beta angle.
		*/
		void _cal_beta_arg(const double& tdt, const double& gst, ColumnVector& beta);

		/// calculate coefficients for Lengendre functions and get _lon, _rho and _factor. 
		/** caculate Pnm coefficients.
		* @remarks Specific formulas can be referred to Geodetic books.
		* @param [in] tdt Dynamic time.
		* @param [in] name Planet name.
		* @param [out] Pnm coefficients.
		* @return
			@retval true	calculate Pnm coefficient sucessfully
			@retval false   calculate Pnm coefficient unsucessfully
		*/
		bool _cal_pnm(const string& name, const double& tdt, Matrix& pnm);

		t_gegmcoef*    _gracoef  = nullptr;	///< the coefficients of EGM model.
		t_gallotl*     _gotl     = nullptr;		///< the ocean_load data .
		t_gocean_tide* _gotdata  = nullptr;	///< the ocean_tide data.
		t_gdeadata*    _gdeadata = nullptr;
		t_gnavde*      _gplanets = nullptr;

		int			   _max_constits = 2000;   ///< Max constits from oceantide data.
		int            _ocean_degree = 30;	   ///< degree of ocean tide
		int            _gravi_degree = 12;     ///< degree of gravity model
		int            _pole_degree = 30;
		bool		   _isfirst = true;	       ///< Mark for wherther to caculate the initial coefficients.
		double		   _distance   = 0.0;      ///< Distance from satellite to geocentric.
		double         _lon        = 0.0;	   ///< The longitude of the planet.
		double         _rho        = 0.0;      ///< The distance between the satellite and the planet. 
		double         _factor     = 0.0;      ///< The factor for dC and dS calculation.
		const int	   _max_ocean_degree = 35; ///< Max ocean degree.
		const int      _max_gravi_degree = 180;///< Max gravity degree.

		bool    _poletide = false;
		Matrix  _Anm;             ///< the matrix A.
		Matrix	_Cnm;		      ///< the matrix Cnm.
		Matrix	_Snm;		      ///< the matrix Snm.
		Matrix  _SUM_F;		      ///< for calculat amat.
		Matrix  _Cnm_save;        ///< legendre coeeficients for EGM.
		Matrix  _Snm_save;        ///< legendre coeeficients for EGM.
		Matrix  _solid_tide_Cnm;  ///< solid tide C correction(spherical harmonic coefficients).
		Matrix  _solid_tide_Snm;  ///< solid tide S correction(spherical harmonic coefficients).
		Matrix  _ocean_tide_Cnm;  ///< ocean tide C correction(spherical harmonic coefficients).
		Matrix  _ocean_tide_Snm;  ///< ocean tide S correction(spherical harmonic coefficients).
		Matrix  _ocean_pole_tide_Cnm;
		Matrix  _ocean_pole_tide_Snm;
		ColumnVector  _sat_pos_norm; ///< Unit Vector of Solar Position.
		ColumnVector  _acc_colum;	 ///<
		
		set<string>                   _planet_list;  ///< force model for solide tide calculation(list: sun moom...).
		vector<t_octrec>              _ocean_tide;   ///< ocean tide data read from file.
		vector<double>                _ocean_load;   ///< ocean load data read from file.
		map<string, ColumnVector> _planet_crd_crs;   ///< The planet pos in J2000.
		vector<double>       _load;
		map <int, map< int, double >> _AR;
		map <int, map< int, double >> _BR;
		map <int, map< int, double >> _AI;
		map <int, map< int, double >> _BI;
	
	};
} // namespace
#endif // !GFGRAVITYEGM_H