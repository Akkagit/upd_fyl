
#ifndef GSYSCONV_H
#define GSYSCONV_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: reference system conversion utilities
  Version: $ Rev: $

  2011-04-26 /PV: created
  2012-04-06 /JD: extracted coordinate system conversion utilities only

-*/

#include <string>

#include "gutils/gtriple.h"
#include "gdata/gsatdata.h"

#include "newmat/newmat.h"
#include "newmat/newmatio.h"
#include "gutils/ExportBase.h"


using namespace std;

namespace gnut {

	BASE_LIBRARY_EXPORT int ell2xyz(const double* Ell, double* XYZ, bool degrees);           // True/False:  LatLon in Degrees/Radians!
	BASE_LIBRARY_EXPORT int xyz2ell(const double* XYZ, double* Ell, bool degrees);           // True/False:  LatLon in Degrees/Radians!
	BASE_LIBRARY_EXPORT int xyz2neu(const double* XYZ, const double* XYZ_Ref, double* neu);
	BASE_LIBRARY_EXPORT int neu2xyz(const double* Ell, const double* neu, double* xyz);  // ! Radians only: Ell[0], Ell[1]

	BASE_LIBRARY_EXPORT int rao2xyz_rot(const ColumnVector& pos, const ColumnVector& vel, Matrix& R);

	BASE_LIBRARY_EXPORT int rao2xyz(const ColumnVector& pos, const ColumnVector& vel,
		const ColumnVector& rao, ColumnVector& xyz);
	BASE_LIBRARY_EXPORT int xyz2rao(const ColumnVector& pos, const ColumnVector& vel,
		ColumnVector& xyz, ColumnVector& rao);

	BASE_LIBRARY_EXPORT int  ell2xyz(const t_gtriple& ell, t_gtriple& xyz, bool degrees);
	BASE_LIBRARY_EXPORT int  xyz2ell(const t_gtriple& crd, t_gtriple& ell, bool degrees);    // True/False:  LatLon in Degrees/Radians!
	BASE_LIBRARY_EXPORT int  xyz2ell2(t_gtriple& crd, t_gtriple& ell, bool degrees);         // True/False:  LatLon in Degrees/Radians!
	BASE_LIBRARY_EXPORT void xyz2neu(t_gtriple& ell, t_gtriple& xyz, t_gtriple& neu);        // Radians only: Ell[0], Ell[1]
	BASE_LIBRARY_EXPORT void neu2xyz(t_gtriple& ell, t_gtriple& neu, t_gtriple& xyz);        // Radians only: Ell[0], Ell[1]

	BASE_LIBRARY_EXPORT int xyz2neu(t_gtriple& xyz, SymmetricMatrix& Q_xyz, SymmetricMatrix& Q_neu);

	BASE_LIBRARY_EXPORT int ell2ipp(t_gsatdata& satdata, t_gtriple& ell_site, t_gtriple& ell_ipp);
	BASE_LIBRARY_EXPORT int ell2ipp(t_gsatdata& satdata, t_gtriple& ell_site, double radius, double ion_hgt, t_gtriple& ell_ipp);

} // namespace

#endif // GSYSCONV_H
