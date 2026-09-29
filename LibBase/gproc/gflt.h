
#ifndef FLT_H
#define FLT_H 

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: implements Kalman filter class
  Version: $ Rev: $

  2011-01-10 /PV: created

-*/

#include "newmat/newmat.h"
#include "newmat/newmatap.h"

#include "gall/gallpar.h"
#include "gutils/ExportBase.h"

using namespace std;

namespace gnut
{

	// Base class for all filtering technique
	// -------------------------------------------------------

	class BASE_LIBRARY_EXPORT t_gflt
	{
	public:
		t_gflt();
		virtual ~t_gflt();

		virtual void update() {}
		virtual void update(const Matrix& A, const DiagonalMatrix& P, const ColumnVector& l, ColumnVector& dx, SymmetricMatrix& Q) {};
		virtual void update(const Matrix& A, const SymmetricMatrix& P, const ColumnVector& l, ColumnVector& dx, SymmetricMatrix& Q) {};
		virtual void update_AmbEl(const Matrix& A, const DiagonalMatrix& P, const ColumnVector& l, ColumnVector& dx, const SymmetricMatrix& Q, SymmetricMatrix& Q_el) {};
		virtual void smooth(t_gallpar& Xp, t_gallpar& Xu, SymmetricMatrix& Qp, SymmetricMatrix& Qu, t_gallpar& Xsm, SymmetricMatrix& Qsm, DiagonalMatrix& Noise) {};
		//added by zhshen
		virtual void add_data(const t_gallpar& param, const ColumnVector& dx, const SymmetricMatrix& Qx, const double& sigma0, const SymmetricMatrix& Qx0);
		virtual void add_data(const Matrix& A, const SymmetricMatrix& P, const ColumnVector& l);
		virtual void add_data(const double& vtpv, const int& nobs_total, const int& npar_number);
		virtual void add_virtual_obs(const Matrix& A, const SymmetricMatrix& P, const ColumnVector& l);
		virtual void get_virtual_obs(Matrix& A, SymmetricMatrix& P, ColumnVector& l);
		virtual bool outlierDetect(const Matrix& A, const SymmetricMatrix& Q, const ColumnVector& dx, const ColumnVector& l, SymmetricMatrix& P);
		virtual void resetQ();


		/** @brief set Qx Matrix specified location */
		void change_Qx(int row, int col, double xx);
		/** @brief set dx  specified location */
		void change_dx(int n, double xx);
		/** @brief set stdx specified location */
		void change_stdx(int n, double xx);
		/** @brief set sigma0 by  specified sigma */
		void change_sigma0(double sigma);
		/** @brief set vtpv by  specified value */
		void change_vtpv(double xx);


		ColumnVector dx() { return _dx; }
		SymmetricMatrix Qx() { return _Qx; }
		t_gallpar param() { return _param; }
		double sigma0() { return _sigma0; }
		double vtpv() { return _vtpv; }
		int nobs_total() { return _nobs_total; }
		int npar_number() { return _npar_number; }
		void amb(bool state) { _amb = state; }
		bool amb() { return _amb; }
	protected:

		ColumnVector _dx;
		SymmetricMatrix _Qx;
		//LX added for WL constraint
		SymmetricMatrix _Qx_tmp;
		SymmetricMatrix _Qx0;
		Matrix _A, _A_virtual;
		SymmetricMatrix _P, _P_virtual;
		ColumnVector _l, _l_virtual;
		t_gallpar _param;
		double _sigma0, _vtpv;
		int _nobs_total, _npar_number;
		bool _amb = false;
	};

	// --------------------------------------------------------


	// Classical formule for Kalman filter
	// -------------------------------------------------------
	class BASE_LIBRARY_EXPORT t_kalman : public t_gflt
	{

	public:

		t_kalman();
		~t_kalman();

		virtual void update();
		void update(const Matrix& A, const DiagonalMatrix& P, const ColumnVector& l, ColumnVector& dx, SymmetricMatrix& Q);
		void update(const Matrix& A, const SymmetricMatrix& P, const ColumnVector& l, ColumnVector& dx, SymmetricMatrix& Q);
	};
	//-------------------------------------------------------

	// Square root covariance filter
	// ------------------------------------------------------
	class BASE_LIBRARY_EXPORT t_SRF : public t_gflt
	{
	public:
		t_SRF();
		~t_SRF();

		virtual void update();
		void update(const Matrix& A, const DiagonalMatrix& P, const ColumnVector& l, ColumnVector& dx, SymmetricMatrix& Q);
		void update(const Matrix& A, const SymmetricMatrix& P, const ColumnVector& l, ColumnVector& dx, SymmetricMatrix& Q);

	};
	// ------------------------------------------------------

	// Square root information filter
	// ------------------------------------------------------
	class BASE_LIBRARY_EXPORT t_SRIF : public t_gflt
	{
	public:
		t_SRIF();
		~t_SRIF();

		void update(const Matrix& A, const DiagonalMatrix& P, const ColumnVector& l, ColumnVector& dx, SymmetricMatrix& Q);

	};

} // namespace

#endif

