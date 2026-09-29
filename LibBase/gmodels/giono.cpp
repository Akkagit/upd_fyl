
#include <cmath>
#include <iomanip>

#include "gmodels/giono.h"
#include "gutils/gconst.h"

using namespace std;

namespace great {   

t_giono::t_giono()
{}


t_giono::~t_giono()
{}


// get ZHD
// ----------
double t_giono::ionMapFunc(IONMPFUNC mapfunc, double ele) // ell v RADIANECH !! TREBA SJEDNOTIT
{   
	if (mapfunc == ICOSZ)
	{
		double mf = 1 / sqrt(1.0 - pow(R_SPHERE / (R_SPHERE + 450000.0) * sin(G_PI / 2.0 - ele), 2));
		return mf;
	}
	else if (mapfunc == QFAC)
	{
		return 0.0;
	}
	else
	{
		return 0.0;
	}
}

double t_giono_brdc::getIonoDelay(t_gallnav* nav, t_gsatdata& satdata, const t_gtime& epo, const t_gtriple& site_ell)
{
	// pre-Judge
	if (site_ell[2] < -1000 || satdata.ele() <= 0) return 0.0;

	t_iono_corr gps_alpha = nav->get_iono_corr(IONO_CORR::IO_GPSA);
	t_iono_corr gps_beta = nav->get_iono_corr(IONO_CORR::IO_GPSB);

	/* earth centered angle (semi-circle) */
	double psi = 0.0137 / (satdata.ele() / G_PI + 0.11) - 0.022;

	/* subionospheric latitude/longitude (semi-circle) */
	double phi = site_ell[0] / G_PI + psi * cos(satdata.azi());
	if (phi > 0.416) phi = 0.416;
	else if (phi < -0.416) phi = -0.416;
	double lam = site_ell[1] / G_PI + psi * sin(satdata.azi()) / cos(phi*G_PI);

	/* geomagnetic latitude (semi-circle) */
	phi += 0.064*cos((lam - 1.617)*G_PI);

	/* local time (s) */
	double tt = 43200.0*lam + epo.sow() + epo.dsec();
	tt -= floor(tt / 86400.0)*86400.0; /* 0<=tt<86400 */

	/* slant factor */
	double f = 1.0 + 16.0*pow(0.53 - satdata.ele() / G_PI, 3.0);

	/* ionospheric delay */
	double amp = gps_alpha.x0 + phi * (gps_alpha.x1 + phi * (gps_alpha.x2 + phi * gps_alpha.x3));
	double per = gps_beta.x0 + phi * (gps_beta.x1 + phi * (gps_beta.x2 + phi * gps_beta.x3));
	amp = amp < 0.0 ? 0.0 : amp;
	per = per < 72000.0 ? 72000.0 : per;
	double x = 2.0*G_PI*(tt - 50400.0) / per;

	return CLIGHT * f*(fabs(x) < 1.57 ? 5e-9 + amp*(1.0 + x * x*(-0.5 + x * x / 24.0)) : 5e-9);
}

bool t_giono_tecgrid::getIonoDelay(t_gionex* ionexdata, t_gsatdata& satdata, const t_gtime& epo, t_gtriple& site_pos, double& value, double& rms)
{
	value = 0.0;
	rms = 0.0;

	t_gtriple site_ell;
	xyz2ell(site_pos, site_ell, false);

	/*elevation and site height check*/
	if (site_ell[2] < -1000 || satdata.ele() <= 0) return false;

	/*Slant total electron content by tec grid data*/
	if (!ionexdata->getSTEC(epo, satdata, site_ell, value, rms)) return false;

	/*Slant ionospheric delay(L1) (m)*/	
	double f1 = satdata.frequency(t_gsys::band_priority(satdata.gsys(), FREQ_1));
	double fact = 40.30 * 1e16 / f1 / f1;
	value *= fact;
	rms *= SQR(fact);

	return true;
}

} // namespace
