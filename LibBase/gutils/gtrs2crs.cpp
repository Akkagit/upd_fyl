
/**
*
* @verbatim
	History
	 -1.0 Qin Yujie         2019-01-06 creat the file.
	 -1.1 Zhang Hongmin		2019-01-07 change sth.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			    gtrs2crs.cpp
* @brief		    calculate the rotation matrix from TRS to CRS and the corresponding partial derivation matrix
*
* @author		    Qin Yujie, Wuhan University
* @version		    1.0.0
* @date			    2019-01-06
*
*/
#include "gutils/gtrs2crs.h"



using namespace std;
namespace great
{

	t_pudaily::t_pudaily()
	{
		xpole = 0.0;
		ypole = 0.0;
		ut1 = 0.0;
	}

	t_EpochNutation::t_EpochNutation()
	{
		T = 0;
		Pi = 0;
		Ep = 0;
	}

	t_zonaltide::t_zonaltide()
	{
		ut1 = 0.0;
		lod = 0.0;
		omega = 0.0;
	}
	t_gtrs2crs::t_gtrs2crs()
	{

	}

	t_gtrs2crs::t_gtrs2crs(bool erptab,t_gpoleut1 * poleut1)
	{
		_erptab = erptab;
		_poleut1 = poleut1;
	}

	t_gtrs2crs& t_gtrs2crs::operator=(const t_gtrs2crs & Other)
	{
		_erptab = Other._erptab;
		_poleut1 = Other._poleut1;
		_erptab = Other._erptab;
		_epo = Other._epo;
		_pudata = Other._pudata;
		_tdt = Other._tdt;
		_qmat = Other._qmat;
		_epsa = Other._epsa;
		_xpole = Other._xpole;
		_ypole = Other._ypole;
		_gmst = Other._gmst;
		_rotmat = Other._rotmat;
		_rotdu = Other._rotdu;
		_rotdx = Other._rotdx;
		_rotdy = Other._rotdy;

		for (int i = 0; i < 5; i++) {
			_arg[i] = Other._arg[i];
		}
		return *this;
	}


	//return the values needed by other functions
	Matrix& t_gtrs2crs::getRotMat() { return _rotmat; };
	Matrix& t_gtrs2crs::getMatDu()  {return _rotdu;};
	Matrix& t_gtrs2crs::getMatDx() { return _rotdx; };  
	Matrix& t_gtrs2crs::getMatDy() { return _rotdy; };  
	double t_gtrs2crs::getXpole() { return _xpole; };  
	double t_gtrs2crs::getYpole() { return _ypole; }; 
	double t_gtrs2crs::getGmst() { return _gmst; };
	t_gtime t_gtrs2crs::getCurtEpoch() { return _tdt; };

	//main
	void t_gtrs2crs::calcRotMat(const t_gtime& epoch, const bool& ldxdpole, const bool& ldydpole, const bool& ldudpole)
	{
		const double dRr = 1.00273781191135448;
		double era;           //earth rotation angle in radian
		double sp;
		double psi, eps, gast;

		//double dRmjd;
		double dUt1_tai;
		string strTemp1;
		string strTemp2;
		double dXhelp[6] = { 0 };
		double dPsi_cor;
		double dEps_cor;

		_tdt = epoch;
		_xpole = 0.0;
		_ypole = 0.0;
		dUt1_tai = 0.0;

		if (!_erptab)
		{
			strTemp1 = "    ";
			strTemp2 = " ";
			_calPoleut1(_tdt, dXhelp, _poleut1, strTemp2);

			_xpole = dXhelp[0] / RAD2SEC;
			_ypole = dXhelp[1] / RAD2SEC;
			dUt1_tai = dXhelp[2];
			dPsi_cor = dXhelp[3] / RAD2SEC;
			dEps_cor = dXhelp[4] / RAD2SEC;
		}

		t_gtime sTUT1;
		sTUT1.from_mjd(_tdt.mjd(), int(_tdt.sod() + _tdt.dsec() + (dUt1_tai - 32.184)), (_tdt.sod() + _tdt.dsec() + (dUt1_tai - 32.184)) - int(_tdt.sod() + _tdt.dsec() + (dUt1_tai - 32.184)));

		sp = _sp2000(_tdt.mjd(), (_tdt.sod()+ _tdt.dsec()) / 86400.0);

		vector<Matrix> roty;
		calcProcMat(ldydpole, 1, _ypole,roty);
		vector<Matrix> rotx;
		calcProcMat(ldxdpole, 2, _xpole,rotx);
		vector<Matrix> rotsp;
		calcProcMat(false, 3, -sp, rotsp);

		era = _era2000(sTUT1.mjd()*1.0, (sTUT1.sod()+sTUT1.dsec()) / 86400.0);

		_nutInt(_tdt.dmjd(), &psi, &eps, 0.0625);

		psi = psi + dPsi_cor;
		eps = eps + dEps_cor;

		_process2000(_tdt.dmjd(), psi, eps);
		gast = _gst2000(_tdt.dmjd(), era, psi);

		vector<Matrix> rotu;
		calcProcMat(ldudpole, 3, -gast,rotu);
		Matrix mathlp;
		mathlp = _qmat * rotu[0];
		if (ldudpole)
		{
			rotu[1] = _qmat * rotu[1];
		}
		mathlp = mathlp * rotsp[0];

		//need test
		if (ldudpole)
		{
			rotu[1] = rotu[1] * rotsp[0];
		}
		if (ldxdpole)
		{
			rotx[1] = mathlp * rotx[1];
		}
		if (ldudpole)
		{
			rotu[1] = rotu[1] * rotx[0];
		}

		mathlp = mathlp * rotx[0];
		_rotmat = mathlp * roty[0];

		if (ldudpole)
		{
			rotu[1] = rotu[1] * roty[0];
			rotu[1] = rotu[1] * dRr;
			
		}
		if (ldxdpole)
		{
			rotx[1] = rotx[1] * roty[0];
		}
		if (ldydpole)
		{
			roty[1] = mathlp * roty[1];
		}
		_xpole = _xpole * RAD2SEC;
		_ypole = _ypole *  RAD2SEC;
		_rotdx = rotx[1];
		_rotdy = roty[1];
		_rotdu = rotu[1];
	}

	void t_gtrs2crs::calcProcMat(const bool& partial, const int& axis, const double& angle, vector<Matrix>& rot)
	{
	    rot.clear();
		Matrix rotmat(0.0, 3, 3);
		Matrix drotmat(0.0,3, 3);

		int i, j;
		double sina, cosa;
		cosa = cos(angle);
		sina = sin(angle);
		rotmat(axis, axis) = 1.0;
		i = axis + 1;
		j = axis + 2;
		if (i > 3)  i = i - 3;
		if (j > 3)  j = j - 3;

		rotmat(i, i) = cosa;
		rotmat(j, j) = cosa;
		rotmat(i, j) = sina;
		rotmat(j, i) = -sina;
		if (partial)
		{
			drotmat(i, i) = -sina;
			drotmat(j, j) = -sina;
			drotmat(i, j) = cosa;
			drotmat(j, i) = -cosa;
		}
		rot.push_back(rotmat);
		rot.push_back(drotmat);
		return;
	}

	void  t_gtrs2crs::_calPoleut1(t_gtime& t, double *x, t_gpoleut1* poleut1, string type)
	{
		int i;
		double dx[6] = { 0.0 }, dutlzonaltide, alpha;
		t_gtriple xyu;
		t_gtime t0;
		double dut1_zonaltide = 0.0;
		t0.from_mjd(_tdt.mjd(), 0, 0.0);
		t_gtime t1;
		t1.from_mjd(_tdt.mjd() + 1, 0, 0.0);
		vector<t_pudaily> rec;
		t_pudaily rec0, rec1;
		rec0.time = t0;
		rec1.time = t1;
		map<string, double> data0, data1;
		data0 = (*(poleut1->getPoleUt1DataMap()))[rec0.time];
		rec0.xpole = data0["XPOLE"];
		rec0.ypole = data0["YPOLE"];
		rec0.ut1 = data0["UT1-TAI"];
		rec.push_back(rec0);
		data1 = (*(poleut1->getPoleUt1DataMap()))[rec1.time];
		rec1.xpole = data1["XPOLE"];
		rec1.ypole = data1["YPOLE"];
		rec1.ut1 = data1["UT1-TAI"];
		rec.push_back(rec1);

		for (i = 0; i < 2; i++)
		{
			if ((*poleut1).getUt1Mode() != "UT1R")
			{
				_tide_corrections(rec[i].time, xyu);
				dutlzonaltide = _tideCor2(rec[i].time.mjd());
			}
			rec[i].xpole = rec[i].xpole - xyu[0];
			rec[i].ypole = rec[i].ypole - xyu[1];
			rec[i].ut1 = rec[i].ut1 - xyu[2] - dut1_zonaltide;
		}
		x[0] = 0;
		x[1] = 0;
		x[2] = 0;
		alpha = (t.dmjd() - rec[0].time.mjd()) / (*poleut1).getIntv();
		x[0] = rec[0].xpole + alpha * (rec[1].xpole - rec[0].xpole);
		x[1] = rec[0].ypole + alpha * (rec[1].ypole - rec[0].ypole);
		x[2] = rec[0].ut1 + alpha * (rec[1].ut1 - rec[0].ut1);
		x[3] = data0["DPSI"] + alpha * (data1["DPSI"] - data0["DPSI"]);
		x[4] = data0["DEPSI"] + alpha * (data1["DEPSI"] - data0["DEPSI"]);

		if (type != "UT1R")
		{
			_tide_corrections(t, xyu);
			dutlzonaltide = _tideCor2(t.dmjd());
			x[0] = x[0] + xyu[0];
			x[1] = x[1] + xyu[1];
			x[2] = x[2] + xyu[2] + dutlzonaltide;
		}


	}

	void t_gtrs2crs::_tide_corrections(t_gtime& t, t_gtriple& xyu)
	{
		double stepsize = 0.015 * 86400.0;//unit sec
		double rmjd = t.dmjd();
		static bool isfirst = true;
		static t_pudaily tb0, tb1;
		if (isfirst)
		{
			tb0.time = LAST_TIME;
			tb1.time = FIRST_TIME;
			isfirst = false;
		}
	
		while (t > tb1.time || t < tb0.time)
		{
			if (t > (tb1.time + stepsize) || t < (tb0.time - stepsize))
			{

				tb0.time = t - stepsize;
				tb1.time = t + stepsize;
				
				tb0 = _tideCor1Cal(tb0);
				tb1 = _tideCor1Cal(tb1);
			}

			else if (t < tb0.time)
			{
				tb1 = tb0;
				tb0.time = tb0.time - stepsize;
				tb0 = _tideCor1Cal(tb0);
			}
			else
			{
				tb0 = tb1;
				tb1.time = tb1.time + stepsize;
				tb1 = _tideCor1Cal(tb1);
			}

		}


		double mjds[2]   = { tb0.time.dmjd(), tb1.time.dmjd() };
		double xpoles[2] = { tb0.xpole,	tb1.xpole };
		double ypoles[2] = { tb0.ypole, tb1.ypole };
		double ut1s[2]   = { tb0.ut1, tb1.ut1 };
		
		xyu[0] = _interpolation(1, 2, mjds, xpoles, rmjd);
		xyu[1] = _interpolation(1, 2, mjds, ypoles, rmjd);
		xyu[2] = _interpolation(1, 2, mjds, ut1s,   rmjd);

		return;
	}


	t_pudaily& t_gtrs2crs::_tideCor1Cal(t_pudaily& b)
	{
		t_gtriple eop;
		//t_gtime time1;
		//time1=b.time
		_ORTHO_EOP(b.time, eop);
		b.xpole = eop[0];
		b.ypole = eop[1];
		b.ut1 = eop[2];

		double xtemp[2];
		_PMSDNUT2(b.time,xtemp);
		b.xpole = b.xpole + xtemp[0];
		b.ypole = b.ypole + xtemp[1];

		//double temp[2];
	   _UTLIBR(b.time,xtemp);
		b.ut1 = b.ut1 + xtemp[0];

		b.xpole = b.xpole*1e-6;
		b.ypole = b.ypole*1e-6;
		b.ut1 =	  b.ut1*1e-6;

		return b;

	}
	void t_gtrs2crs::_ORTHO_EOP(t_gtime& t, t_gtriple& eop)
	{
		
		int k, j;
		double orthow[12][3] = {
			-6.77832, 14.86283, -1.76335,
			-14.86323, -6.77846, 1.03364,
			0.47884, 1.45234, -0.27553,
			-1.45303, 0.47888, 0.34569,
			0.16406, -0.42056, -0.12343,
			0.42030, 0.16469, -0.10146,
			0.09398, 15.30276, -0.47119,
			25.73054, -4.30615, 1.28997,
			-4.77974, 0.07564, -0.19336,
			0.28080, 2.28321, 0.02724,
			1.94539, -0.45717, 0.08955,
			-0.73089, -1.62010, 0.04726 };
		double h[12];
	    _CNMTX(t,h);
		for (k = 0; k < 3; k++)
		{
			eop[k] = 0.0;
			for (j = 0; j < 12; j++)
			{
				eop[k] = eop[k] + h[j] * orthow[j][k];
			}

		}
		return;

	}
	void t_gtrs2crs::_CNMTX(t_gtime& t,double *h)
	{
		int i, j, k, m, n, nmax;
		const int nlines = 71;
		double dt60, d1960, dt;
		double nj[nlines], mj[nlines];
		double anm[2][4][3], bnm[2][4][3];
		double ap, am, bp, bm, pinm, alpha;
		double p[3][2], q[3][3];
		//double h[12];
		double sp[6][2] = { { 0.0298, 0.0200 }, { 0.1408, 0.0905 }, { 0.0805, 0.0638 },
		{ 0.6002, 0.3476 }, { 0.3025, 0.1645 }, { 0.1517, 0.0923 } };
		dt = 2.0;
		nmax = 2;
		d1960 = 37076.5;
		for (i = 0; i < 41; i++)
		{
			nj[i] = 2.0;
			mj[i] = 1.0;
		}
		for (i = 41; i < 71; i++)
		{
			nj[i] = 2.0;
			mj[i] = 2.0;
		}
		double hs[71] =
		{
			-1.94, -1.25, -6.64, -1.51, -8.02, -9.47, -50.20, -1.80, -9.54, 1.52, -49.45, -262.21, 1.70, 3.43, 1.94,
			1.37, 7.41, 20.62, 4.14, 3.94, -7.14, 1.37, -122.03, 1.02, 2.89, -7.3, 368.78, 50.01, -1.08, 2.93,
			5.25, 3.95, 20.62, 4.09, 3.42, 1.69, 11.29, 7.23, 1.51, 2.16, 1.38, 1.80, 4.67, 16.01, 19.32,
			1.30, -1.02, -4.51, 120.99, 1.13, 22.98, 1.06, -1.9, -2.18, -23.58, 631.92, 1.92, -4.66, -17.86, 4.47,
			1.97, 17.2, 294.00, -2.46, -1.02, 79.96, 23.83, 2.59, 4.47, 1.95, 1.17
		};
		double phase[71] =
		{
			9.0899831, 8.8234208, 12.1189598, 1.4425700, 4.7381090, 4.4715466, 7.7670857, -2.9093042, 0.3862349, -3.1758666, 0.1196725, 3.4152116, 12.8946194, 5.5137686, 6.4441883,
			-4.2322016, -0.9366625, 8.5427453, 11.8382843, 1.1618945, 5.9693878, -1.2032249, 2.0923141, -1.7847596, 8.0679449, 0.8953321, 4.1908712, 7.4864102, 10.7819493, 0.3137975,
			6.2894282, 7.2198478, -0.161003, 3.1345361, 2.8679737, -4.5128771, 4.9665307, 8.2620698, 11.5576089, 0.6146566, 3.9101957, 20.6617051, 13.2808543, 16.309831, 8.9289802,
			5.0519065, 15.8350306, 8.6624178, 11.9579569, 8.0808832, 4.5771061, 0.7000324, 14.9869335, 11.4831564, 4.3105437, 7.6060827, 3.729009, 10.6350594, 3.2542086, 12.7336164,
			16.0291555, 10.160259, 6.2831853, 2.4061116, 5.0862033, 8.3817423, 11.6772814, 14.9728205, 4.0298682, 7.3254073, 9.1574019
		};
		double freq[71] =
		{
			5.18688050, 5.38346657, 5.38439079, 5.41398343, 5.41490765, 5.61149372, 5.61241794, 5.64201057, 5.64293479, 5.83859664, 5.83952086, 5.84044508, 5.84433381, 5.87485066, 6.03795537,
			6.06754801, 6.06847223, 6.07236095, 6.07328517, 6.10287781, 6.24878055, 6.26505830, 6.26598252, 6.28318449, 6.28318613, 6.29946388, 6.3003881, 6.30131232, 6.30223654, 6.31759007,
			6.33479368, 6.49789839, 6.52841524, 6.52933946, 6.72592553, 6.75644239, 6.76033111, 6.76125533, 6.76217955, 6.98835826, 6.98928248, 11.45675174, 11.4872686, 11.68477889, 11.71529575,
			11.73249771, 11.89560406, 11.91188181, 11.91280603, 11.93000800, 11.94332289, 11.96052486, 12.11031632, 12.12363121, 12.13990896, 12.14083318, 12.15803515, 12.33834347, 12.36886033, 12.37274905,
			12.37367327, 12.54916865, 12.56637061, 12.58357258, 12.59985198, 12.6007762, 12.60170041, 12.60262463, 12.82880334, 12.82972756, 13.06071921
		};
		string numarg[71] =
		{
			"117.655", "125.745", "125.755", "127.545", "127.555", "135.645", "135.655", "137.445", "137.455", "145.535", "145.545", "145.555", "145.755", "147.555", "153.655",
			"155.445", "155.455", "155.655", "155.665", "157.455", "162.556", "163.545", "163.555", "164.554", "164.556", "165.545", "165.555", "165.565", "165.575", "166.554",
			"167.555", "173.655", "175.455", "175.465", "183.555", "185.355", "185.555", "185.565", "185.575", "195.455", "195.465", "225.855", "227.655", "235.755", "237.555",
			"238.554", "244.656", "245.645", "245.655", "246.654", "247.455", "248.454", "253.755", "254.556", "255.545", "255.555", "256.554", "263.655", "265.455", "265.655",
			"265.665", "272.556", "273.555", "274.554", "275.545", "275.555", "275.565", "275.575", "285.455", "285.465", "295.555"
		};


		for (k = -1; k <= 1; k++)
		{
			dt60 = (t.dmjd() - k * dt) - d1960;
			anm[0][1][k + 1] = 0;
			anm[0][2][k + 1] = 0;
			bnm[0][1][k + 1] = 0;
			bnm[0][2][k + 1] = 0;
			for (j = 0; j < nlines; j++)
			{
				n = nj[j];
				m = mj[j];
				pinm = ((n + m) % 2)*d2PI / 4;
				alpha = fmod(phase[j] - pinm, d2PI) + fmod(freq[j] * dt60, d2PI);
				anm[n - 2][m][k + 1] = anm[n - 2][m][k + 1] + hs[j] * cos(alpha);
				bnm[n - 2][m][k + 1] = bnm[n - 2][m][k + 1] - hs[j] * sin(alpha);
			}
		}
		for (m = 1; m <= 2; m++)
		{
			ap = anm[0][m][2] + anm[0][m][0];
			am = anm[0][m][2] - anm[0][m][0];
			bp = bnm[0][m][2] + bnm[0][m][0];
			bm = bnm[0][m][2] - bnm[0][m][0];
			p[0][m - 1] = sp[0][m - 1] * anm[0][m][1];
			p[1][m - 1] = sp[1][m - 1] * anm[0][m][1] - sp[2][m - 1] * ap;
			p[2][m - 1] = sp[3][m - 1] * anm[0][m][1] - sp[4][m - 1] * ap + sp[5][m - 1] * bm;
			q[0][m - 1] = sp[0][m - 1] * bnm[0][m][1];
			q[1][m - 1] = sp[1][m - 1] * bnm[0][m][1] - sp[2][m - 1] * bp;
			q[2][m - 1] = sp[3][m - 1] * bnm[0][m][1] - sp[4][m - 1] * bp - sp[5][m - 1] * am;
			anm[0][m][0] = p[0][m - 1];
			anm[0][m][1] = p[1][m - 1];
			anm[0][m][2] = p[2][m - 1];
			bnm[0][m][0] = q[0][m - 1];
			bnm[0][m][1] = q[1][m - 1];
			bnm[0][m][2] = q[2][m - 1];
		}

		j = 0;
		for (n = 2; n <= nmax; n++)
		{
			for (m = 1; m <= n; m++)
			{
				for (k = -1; k <= 1; k++)
				{
					h[j] = anm[n - 2][m][k + 1];
					h[j + 1] = bnm[n - 2][m][k + 1];
					j = j + 2;
				}
			}
		}
		return ;
	}
	void t_gtrs2crs::_PMSDNUT2(t_gtime& t,double* pm)
	{
		int i, j, jstart;
		//double pm[2];

		//todo 一些常量可以使用外部？？
		double tt, gmst, l, lp, f, d, om;
		double angle, xrate, yrate;
		const int iband = 1;

		double iarg[6][25] =
		{
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			0, -1, -1, -1, 1, 1, 0, 1, 0, 0, 0, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, 0, 0, 0, 1,
			0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 1, 1, 1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 3, 3, -2, -2, -2, -2, -2, 0, -2, 0, 0, 0,
			0, 0, 0, 0, 0, 0, -1, -2, 0, 0, 0, 2, 0, 0, 0, 0, 0, -2, 0, 0, 0, 2, 0, 0, 0,
			-1, 2, 1, 0, 0, -1, 1, 1, 2, 1, 0, 1, 1, 3, 2, -1, -2, -2, -1, -2, 0, -2, 0, -1, 0
		};
		double per[25] =
		{
			6798.3837, 6159.1355, 3231.4956, 2190.3501, 438.3599, 411.80661, 365.24219, 193.55971, 27.431826, 27.321582, 27.212221, 14.698136, 13.718786,
			9.1071941, 9.0950103, 1.1196992, 1.1195149, 1.1134606, 1.0759762, 1.0758059, 1.0347187, 1.0027454, 0.9972696, 0.9971233, 0.9624365
		};
		double xs[25] =
		{
			0, 1.5, -28.5, -4.7, -0.7, 1, 1.2, 1.3, -0.1, 0.9, 0.1, 0, -0.1, -0.1, -0.1, -0.4, -2.3, -0.4, -2.1, -11.4, 0.8, -4.8, 14.3, 1.9, 0.8
		};
		double xc[25] =
		{
			0.6, 0.0, -0.2, -0.1, 0.2, 0.3, 0.2, 0.4, -0.2, 4.0, 0.6, 0.1, 0.3, 0.1, 0.1, 0.3, 1.3, 0.3, 1.2, 6.5, -0.5, 2.7, -8.2, -1.1, -0.4
		};
		double ys[25] =
		{
			-0.1, -0.2, 3.4, 0.6, -0.2, -0.3, -0.2, -0.2, 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, -0.3, -1.3, -0.3, -1.2, -6.5, 0.5, -2.7, 8.2, 1.1, 0.4
		};
		double yc[25] =
		{
			-0.1, 0.1, -3.9, -0.9, -0.7, 1.0, 1.4, 2.9, -1.7, 32.4, 5.1, 0.6, 2.7, 0.9, 0.6, -0.4, -2.3, -0.4, -2.1, -11.4, 0.8, -4.8, 14.3, 1.9, 0.8
		};
		xrate = -3.8;
		yrate = -4.3;

		pm[0] = 0;
		pm[1] = 0;
		tt = (t.dmjd() - rmjd0) / 36525.0;
		gmst = fmod(67310.54841 +
			tt * ((8640184.812866 + 3155760000.0) +
				tt * (0.093104 +
					tt * (-0.0000062))), 86400.0);

		_FUNDARG(tt, &l, &lp, &f, &d, &om);
		double arg[6];
		arg[0] = gmst / rad2sec + dPI;
		arg[0] = fmod(arg[0], d2PI);
		arg[1] = l;
		arg[2] = lp;
		arg[3] = f;
		arg[4] = d;
		arg[5] = om;

		if (iband == 1)
		{
			jstart = 16;
		}
		else
		{
			jstart = 1;
		}
		for (j = jstart - 1; j < 25; j++)
		{
			angle = 0.0;
			for (i = 0; i < 6; i++)
			{
				angle = angle + iarg[i][j] * arg[i];
			}
			angle = fmod(angle, d2PI);

			pm[0] = pm[0] + xs[j] * sin(angle) + xc[j] * cos(angle);
			pm[1] = pm[1] + ys[j] * sin(angle) + yc[j] * cos(angle);
		
		}


		if (iband == 1)
		{
			return;
		}
		pm[0] = pm[0] + xrate * (t.dmjd() - rmjd0) / 365.25;
		pm[1] = pm[1] + yrate * (t.dmjd() - rmjd0) / 365.25;
		return;
	}

	void   t_gtrs2crs::_UTLIBR(t_gtime& t,double* temp)
	{
		int i, j, iarg[6][11];
		//double temp[2];
		double tt;
		double gmst, l, lp, f, d, om, arg[6], per[11], dut1s[11], dut1c[11], dlods[11], dlodc[11], angle;
		int iall[11][6] =
		{
			2, -2, 0, -2, 0, -2,
			2, 0, 0, -2, -2, -2,
			2, -1, 0, -2, 0, -2,
			2, 1, 0, -2, -2, -2,
			2, 0, 0, -2, 0, -1,
			2, 0, 0, -2, 0, -2,
			2, 1, 0, -2, 0, -2,
			2, 0, -1, -2, 2, -2,
			2, 0, 0, -2, 2, -2,
			2, 0, 0, 0, 0, 0,
			2, 0, 0, 0, 0, -1
		};
		for (i = 0; i <= 10; i++)
		{
			for (j = 0; j <= 5; j++)
			{
				iarg[j][i] = iall[i][j];
			}
		}
		double dall[11][5] =
		{
			0.5377239, 0.05, -0.03, -0.3, -0.6,
			0.5363232, 0.06, -0.03, -0.4, -0.7,
			0.5274312, 0.35, -0.20, -2.4, -4.1,
			0.5260835, 0.07, -0.04, -0.5, -0.8,
			0.5175645, -0.07, 0.04, 0.5, 0.8,
			0.5175251, 1.75, -1.01, -12.2, -21.3,
			0.5079842, -0.05, 0.03, 0.3, 0.6,
			0.5006854, 0.04, -0.03, -0.3, -0.6,
			0.5000000, 0.76, -0.44, -5.5, -9.6,
			0.4986348, 0.21, -0.12, -1.5, -2.6,
			0.4985982, 0.06, -0.04, -0.4, -0.8
		};
		int ii = 0;
		for (ii = 0; ii <= 10; ii++)
		{
			per[ii] = dall[ii][0];
			dut1s[ii] = dall[ii][1];
			dut1c[ii] = dall[ii][2];
			dlods[ii] = dall[ii][3];
			dlodc[ii] = dall[ii][4];
		}
		tt = (t.dmjd() - rmjd0) / 36525;
		gmst = fmod((67310.54841 + tt * ((8640184.812866 + 3155760000) + tt * (0.093104 + tt * (-0.0000062)))), 86400.0);
		_FUNDARG(tt, &l, &lp, &f, &d, &om);
		arg[0] = gmst / rad2sec + dPI;
		arg[0] = fmod(arg[0], d2PI);
		arg[1] = l;
		arg[2] = lp;
		arg[3] = f;
		arg[4] = d;
		arg[5] = om;
		temp[0] = 0;
		temp[1] = 0;
		for (j = 0; j <= 10; j++)
		{
			angle = 0;
			for (i = 0; i <= 5; i++)
			{
				angle = angle + iarg[i][j] * arg[i];
			}

			angle = fmod(angle, d2PI);

			temp[0] = temp[0] + dut1s[j] * sin(angle) + dut1c[j] * cos(angle);
			temp[1] = temp[1] + dlods[j] * sin(angle) + dlodc[j] * cos(angle);

		}
		return;
	}

	

	double t_gtrs2crs::_tideCor2(const double& dRmjd)
	{
		static t_zonaltide sTZB[3];
		static double gdStepsize = 0.05; 
		static bool isFirst = true; 
		//int i;
		double dT;
		double pdUt1;


		if (isFirst)
		{
			sTZB[0].time = LAST_TIME;
			sTZB[2].time = FIRST_TIME;
			isFirst = false;
		}

		while ((dRmjd > sTZB[2].time.dmjd()) || (dRmjd < sTZB[0].time.dmjd()))
		{
			if ((dRmjd > sTZB[2].time.dmjd() + gdStepsize) || (dRmjd < sTZB[0].time.dmjd() - gdStepsize))
			{

				sTZB[0].time.from_mjd((int)(dRmjd - gdStepsize), (int)((dRmjd - gdStepsize-(int)(dRmjd - gdStepsize))*86400), ((dRmjd - gdStepsize - (int)(dRmjd - gdStepsize)) * 86400)-(int)(((dRmjd - gdStepsize - (int)(dRmjd - gdStepsize)) * 86400)));
				//sTZB[1].dT = dRmjd;
				sTZB[1].time.from_mjd((int)(dRmjd), (int)((dRmjd-(int)dRmjd)*86400), ((dRmjd - (int)dRmjd) * 86400)-(int)((dRmjd - (int)dRmjd) * 86400));
				//sTZB[2].dT = dRmjd + gdStepsize;
				sTZB[2].time.from_mjd((int)(dRmjd + gdStepsize), (int)((dRmjd + gdStepsize - (int)(dRmjd + gdStepsize)) * 86400), ((dRmjd + gdStepsize - (int)(dRmjd + gdStepsize)) * 86400) - (int)(((dRmjd + gdStepsize - (int)(dRmjd + gdStepsize)) * 86400)));
				for (int i = 0; i < 3; i++)
				{
					dT = (sTZB[i].time.dmjd() - 51544.5) / 36525.0;
					_RG_ZONT2(dT, &sTZB[i].ut1, &sTZB[i].lod, &sTZB[i].omega);
				}

			}

			else if (dRmjd < sTZB[0].time.dmjd())
			{
				sTZB[2] = sTZB[1];
				sTZB[1] = sTZB[0];
				sTZB[0].time.from_mjd((int)(sTZB[0].time.dmjd() - gdStepsize), (int)(((sTZB[0].time.dmjd() - gdStepsize) - (int)(sTZB[0].time.dmjd() - gdStepsize)) * 86400), (((sTZB[0].time.dmjd() - gdStepsize) - (int)(sTZB[0].time.dmjd() - gdStepsize)) * 86400) - (int)(((sTZB[0].time.dmjd() - gdStepsize) - (int)(sTZB[0].time.dmjd() - gdStepsize)) * 86400));
				dT = (sTZB[0].time.dmjd() - 51544.5) / 36525.0;
				_RG_ZONT2(dT, &sTZB[0].ut1, &sTZB[0].lod, &sTZB[0].omega);
			}

			else
			{
				sTZB[0] = sTZB[1];
				sTZB[1] = sTZB[2];
				//sTZB[2].dT = sTZB[2].time.dmjd() + gdStepsize;
				sTZB[2].time.from_mjd((int)(sTZB[2].time.dmjd() + gdStepsize), (int)(((sTZB[2].time.dmjd() + gdStepsize) - (int)(sTZB[2].time.dmjd() + gdStepsize)) * 86400), (((sTZB[2].time.dmjd() + gdStepsize) - (int)(sTZB[2].time.dmjd() + gdStepsize)) * 86400) - (int)(((sTZB[2].time.dmjd() + gdStepsize) - (int)(sTZB[2].time.dmjd() + gdStepsize)) * 86400));
				dT = (sTZB[2].time.dmjd() - 51544.5) / 36525.0;
				_RG_ZONT2(dT, &sTZB[2].ut1, &sTZB[2].lod, &sTZB[2].omega);
			}

		}

		//    	 dut1   = r_interpolation(2,3,TB(:)%t,TB(:)%dut1,rmjd)

		double dTemp1[3];
		double dTemp2[3];
		for (int i = 0; i < 3; i++)
		{
			dTemp1[i] = sTZB[i].time.dmjd();
			dTemp2[i] = sTZB[i].ut1;
		}

		pdUt1 = _interpolation(2, 3, dTemp1, dTemp2, dRmjd);
		return pdUt1;

	}

	void  t_gtrs2crs::_RG_ZONT2(const double& dT, double *DUT, double *DLOD, double *DOMEGA)
	{
		int I;
		//int J;
		double dL;
		double dLP;
		double dF;
		double dD;
		double dOM;
		double dARG;

		/*  ----------------------
		*  Zonal Earth tide model
		*  ----------------------*/
		//*  Number of terms in the zonal Earth tide model
		//  INTEGER NZONT
		const int iNZONT = 62;
		//*  Coefficients for the fundamental arguments
		int iNFUND[5][iNZONT];
		//*  Zonal tide term coefficients
		double dTIDE[6][iNZONT];
		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		*  --------------------------------------------------
		*  Tables of multiples of arguments and coefficients
		*  --------------------------------------------------
		*  Luni-Solar argument multipliers*/
		
		int iNFUNDSUB1[20][5] =
		{
			1, 0, 2, 2, 2,
			2, 0, 2, 0, 1,
			2, 0, 2, 0, 2,
			0, 0, 2, 2, 1,
			0, 0, 2, 2, 2,
			1, 0, 2, 0, 0,
			1, 0, 2, 0, 1,
			1, 0, 2, 0, 2,
			3, 0, 0, 0, 0,
			-1, 0, 2, 2, 1,
			-1, 0, 2, 2, 2,
			1, 0, 0, 2, 0,
			2, 0, 2, -2, 2,
			0, 1, 2, 0, 2,
			0, 0, 2, 0, 0,
			0, 0, 2, 0, 1,
			0, 0, 2, 0, 2,
			2, 0, 0, 0, -1,
			2, 0, 0, 0, 0,
			2, 0, 0, 0, 1
		};
		int i = 0;
		int j = 0;
		for (i = 0; i <= 19; i++)
		{
			for (j = 0; j <= 4; j++)
			{
				iNFUND[j][i] = iNFUNDSUB1[i][j];
			}
		}
		
		int iNFUNDSUB2[20][5] =
		{
			0, -1, 2, 0, 2,
			0, 0, 0, 2, -1,
			0, 0, 0, 2, 0,
			0, 0, 0, 2, 1,
			0, -1, 0, 2, 0,
			1, 0, 2, -2, 1,
			1, 0, 2, -2, 2,
			1, 1, 0, 0, 0,
			-1, 0, 2, 0, 0,
			-1, 0, 2, 0, 1,
			-1, 0, 2, 0, 2,
			1, 0, 0, 0, -1,
			1, 0, 0, 0, 0,
			1, 0, 0, 0, 1,
			0, 0, 0, 1, 0,
			1, -1, 0, 0, 0,
			-1, 0, 0, 2, -1,
			-1, 0, 0, 2, 0,
			-1, 0, 0, 2, 1,
			1, 0, -2, 2, -1
		};
		for (i = 0; i <= 19; i++)
		{
			for (j = 0; j <= 4; j++)
			{
				iNFUND[j][i + 20] = iNFUNDSUB2[i][j];
			}
		}
		
		int iNFUNDSUB3[20][5] =
		{
			-1, -1, 0, 2, 0,
			0, 2, 2, -2, 2,
			0, 1, 2, -2, 1,
			0, 1, 2, -2, 2,
			0, 0, 2, -2, 0,
			0, 0, 2, -2, 1,
			0, 0, 2, -2, 2,
			0, 2, 0, 0, 0,
			2, 0, 0, -2, -1,
			2, 0, 0, -2, 0,
			2, 0, 0, -2, 1,
			0, -1, 2, -2, 1,
			0, 1, 0, 0, -1,
			0, -1, 2, -2, 2,
			0, 1, 0, 0, 0,
			0, 1, 0, 0, 1,
			1, 0, 0, -1, 0,
			2, 0, -2, 0, 0,
			-2, 0, 2, 0, 1,
			-1, 1, 0, 1, 0
		};
		for (i = 0; i <= 19; i++)
		{
			for (j = 0; j <= 4; j++)
			{
				iNFUND[j][i + 40] = iNFUNDSUB3[i][j];
			}
		}
	
		int iNFUNDSUB4[2][5] =
		{
			0, 0, 0, 0, 2,
			0, 0, 0, 0, 1
		};
		for (i = 0; i <= 1; i++)
		{
			for (j = 0; j <= 4; j++)
			{
				iNFUND[j][i + 60] = iNFUNDSUB4[i][j];
			}
		}
		
		double dTIDESUB1[20][6] =
		{
			-0.0235, 0.0000, 0.2617, 0.0000, -0.2209, 0.0000,
			-0.0404, 0.0000, 0.3706, 0.0000, -0.3128, 0.0000,
			-0.0987, 0.0000, 0.9041, 0.0000, -0.7630, 0.0000,
			-0.0508, 0.0000, 0.4499, 0.0000, -0.3797, 0.0000,
			-0.1231, 0.0000, 1.0904, 0.0000, -0.9203, 0.0000,
			-0.0385, 0.0000, 0.2659, 0.0000, -0.2244, 0.0000,
			-0.4108, 0.0000, 2.8298, 0.0000, -2.3884, 0.0000,
			-0.9926, 0.0000, 6.8291, 0.0000, -5.7637, 0.0000,
			-0.0179, 0.0000, 0.1222, 0.0000, -0.1031, 0.0000,
			-0.0818, 0.0000, 0.5384, 0.0000, -0.4544, 0.0000,
			-0.1974, 0.0000, 1.2978, 0.0000, -1.0953, 0.0000,
			-0.0761, 0.0000, 0.4976, 0.0000, -0.4200, 0.0000,
			0.0216, 0.0000, -0.1060, 0.0000, 0.0895, 0.0000,
			0.0254, 0.0000, -0.1211, 0.0000, 0.1022, 0.0000,
			-0.2989, 0.0000, 1.3804, 0.0000, -1.1650, 0.0000,
			-3.1873, 0.2010, 14.6890, 0.9266, -12.3974, -0.7820,
			-7.8468, 0.5320, 36.0910, 2.4469, -30.4606, -2.0652,
			0.0216, 0.0000, -0.0988, 0.0000, 0.0834, 0.0000,
			-0.3384, 0.0000, 1.5433, 0.0000, -1.3025, 0.0000,
			0.0179, 0.0000, -0.0813, 0.0000, 0.0686, 0.0000
		};
		for (i = 0; i <= 19; i++)
		{
			for (j = 0; j <= 5; j++)
			{
				dTIDE[j][i] = dTIDESUB1[i][j];
			}
		}

		double dTIDESUB2[20][6] =
		{
			-0.0244, 0.0000, 0.1082, 0.0000, -0.0913, 0.0000,
			0.0470, 0.0000, -0.2004, 0.0000, 0.1692, 0.0000,
			-0.7341, 0.0000, 3.1240, 0.0000, -2.6367, 0.0000,
			-0.0526, 0.0000, 0.2235, 0.0000, -0.1886, 0.0000,
			-0.0508, 0.0000, 0.2073, 0.0000, -0.1749, 0.0000,
			0.0498, 0.0000, -0.1312, 0.0000, 0.1107, 0.0000,
			0.1006, 0.0000, -0.2640, 0.0000, 0.2228, 0.0000,
			0.0395, 0.0000, -0.0968, 0.0000, 0.0817, 0.0000,
			0.0470, 0.0000, -0.1099, 0.0000, 0.0927, 0.0000,
			0.1767, 0.0000, -0.4115, 0.0000, 0.3473, 0.0000,
			0.4352, 0.0000, -1.0093, 0.0000, 0.8519, 0.0000,
			0.5339, 0.0000, -1.2224, 0.0000, 1.0317, 0.0000,
			-8.4046, 0.2500, 19.1647, 0.5701, -16.1749, -0.4811,
			0.5443, 0.0000, -1.2360, 0.0000, 1.0432, 0.0000,
			0.0470, 0.0000, -0.1000, 0.0000, 0.0844, 0.0000,
			-0.0555, 0.0000, 0.1169, 0.0000, -0.0987, 0.0000,
			0.1175, 0.0000, -0.2332, 0.0000, 0.1968, 0.0000,
			-1.8236, 0.0000, 3.6018, 0.0000, -3.0399, 0.0000,
			0.1316, 0.0000, -0.2587, 0.0000, 0.2183, 0.0000,
			0.0179, 0.0000, -0.0344, 0.0000, 0.0290, 0.0000

		};
		for (i = 0; i <= 19; i++)
		{
			for (j = 0; j <= 5; j++)
			{
				dTIDE[j][i + 20] = dTIDESUB2[i][j];
			}
		}

		double dTIDESUB3[20][6] =
		{
			-0.0855, 0.0000, 0.1542, 0.0000, -0.1302, 0.0000,
			-0.0573, 0.0000, 0.0395, 0.0000, -0.0333, 0.0000,
			0.0329, 0.0000, -0.0173, 0.0000, 0.0146, 0.0000,
			-1.8847, 0.0000, 0.9726, 0.0000, -0.8209, 0.0000,
			0.2510, 0.0000, -0.0910, 0.0000, 0.0768, 0.0000,
			1.1703, 0.0000, -0.4135, 0.0000, 0.3490, 0.0000,
			-49.7174, 0.4330, 17.1056, 0.1490, -14.4370, -0.1257,
			-0.1936, 0.0000, 0.0666, 0.0000, -0.0562, 0.0000,
			0.0489, 0.0000, -0.0154, 0.0000, 0.0130, 0.0000,
			-0.5471, 0.0000, 0.1670, 0.0000, -0.1409, 0.0000,
			0.0367, 0.0000, -0.0108, 0.0000, 0.0092, 0.0000,
			-0.0451, 0.0000, 0.0082, 0.0000, -0.0069, 0.0000,
			0.0921, 0.0000, -0.0167, 0.0000, 0.0141, 0.0000,
			0.8281, 0.0000, -0.1425, 0.0000, 0.1202, 0.0000,
			-15.8887, 0.1530, 2.7332, 0.0267, -2.3068, -0.0222,
			-0.1382, 0.0000, 0.0225, 0.0000, -0.0190, 0.0000,
			0.0348, 0.0000, -0.0053, 0.0000, 0.0045, 0.0000,
			-0.1372, 0.0000, -0.0079, 0.0000, 0.0066, 0.0000,
			0.4211, 0.0000, -0.0203, 0.0000, 0.0171, 0.0000,
			-0.0404, 0.0000, 0.0008, 0.0000, -0.0007, 0.0000
		};
		for (i = 0; i <= 19; i++)
		{
			for (j = 0; j <= 5; j++)
			{
				dTIDE[j][i + 40] = dTIDESUB3[i][j];
			}
		}

		double dTIDESUB4[2][6] =
		{
			7.8998, 0.0000, 0.1460, 0.0000, -0.1232, 0.0000,
			-1617.2681, 0.0000, -14.9471, 0.0000, 12.6153, 0.0000
		};
		for (i = 0; i <= 1; i++)
		{
			for (j = 0; j <= 5; j++)
			{
				dTIDE[j][i + 60] = dTIDESUB4[i][j];
			}
		}
		/*  -------------------------------------
		 *   Computation of fundamental arguments
		 *  -------------------------------------*/

		_FUNDARG(dT, &dL, &dLP, &dF, &dD, &dOM);

		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		// Set initial values to zero.
		*DUT = 0;
		*DLOD = 0;
		*DOMEGA = 0;

		//*  Sum zonal tide terms.
		for (I = 0; I <= 40; I++)
		{
			//*     Formation of multiples of arguments.
			dARG = fmod((iNFUND[0][I] * dL + iNFUND[1][I] * dLP + iNFUND[2][I] * dF + iNFUND[3][I] * dD + iNFUND[4][I] * dOM), d2PI);

			if (dARG < 0)
			{
				dARG = dARG + d2PI;
			}
			//*     Evaluate zonal tidal terms.
			*DUT = *DUT + dTIDE[0][I] * sin(dARG) + dTIDE[1][I] * cos(dARG);

			*DLOD = *DLOD + dTIDE[2][I] * cos(dARG) + dTIDE[3][I] * sin(dARG);

			*DOMEGA = *DOMEGA + dTIDE[4][I] * cos(dARG) + dTIDE[5][I] * sin(dARG);
		}
		//*  Rescale corrections so that they are in units involving seconds.

		*DUT = (*DUT) * 1.0E-4;
		*DLOD = (*DLOD) * 1.0E-5;
		*DOMEGA = (*DOMEGA) * 1.0E-14;
		return;
	}

	void t_gtrs2crs::_FUNDARG(const double& T, double *L, double *LP, double *F, double *D, double *OM)
	{
		if (L == NULL || LP == NULL || F == NULL || D == NULL || OM == NULL)
		{
			cout<< "参数不正确";
		}
		
		const double DAS2R = 4.848136811095359935899141e-6;
		const double TURNAS = 1296000;
		const double D2PI = 6.283185307179586476925287;
	
		*L = fmod(485868.249036 +
			T * (1717915923.2178 +
				T * (31.8792 +
					T * (0.051635 +
						T * (-0.00024470)))), TURNAS) * DAS2R;
		*LP = fmod(1287104.79305 +
			T * (129596581.0481 +
				T * (-0.5532 +
					T * (0.000136 +
						T * (-0.00001149)))), TURNAS)*DAS2R;

		*F = fmod(335779.526232 +
			T * (1739527262.8478 +
				T * (-12.7512 +
					T * (-0.001037 +
						T * (0.00000417)))), TURNAS)*DAS2R;
		*D = fmod(1072260.70369 +
			T * (1602961601.2090 +
				T * (-6.3706 +
					T * (0.006593 +
						T * (-0.00003169)))), TURNAS)*DAS2R;
		*OM = fmod(450160.398036 +
			T * (-6962890.5431 +
				T * (7.4722 +
					T * (0.007702 +
						T * (-0.00005939)))), TURNAS)*DAS2R;
		return;
	}


	/***********************
	*FunctionName：nutation_interpolation
		*Function：
		*InPut：  	dRmjd modified JD
		dStepsize stepsize in days, 0.d0 use the default 0.125
		*OutPut：dpi nutation angle of the longitude(radian)
		dep nutation angle of the obliquity(radian)
		*Return：
		*Other：
	***********************/
	void  t_gtrs2crs::_nutInt(const double& dRmjd, double* dpsi, double* deps, const double& step)
	{
		static bool isFirst = true;  
		static double gdStepsize_used = 0.125;  
		static t_EpochNutation sTB1[3];  
		
		if ((isFirst) || fabs(gdStepsize_used - step) > pow(10.0, -10))
		{
			isFirst = false;
			sTB1[0].T = 1.0e20;
			sTB1[2].T = -1.0e-20;
			if (step != 0.0)
			{
				gdStepsize_used = step;
			}
		}
		
		if ((dRmjd > sTB1[2].T) || (dRmjd < sTB1[0].T))
		{
			if ((dRmjd > (sTB1[2].T + gdStepsize_used)) || (dRmjd < sTB1[0].T - gdStepsize_used))
			{
				sTB1[0].T = dRmjd - gdStepsize_used;
				sTB1[1].T = dRmjd;
				sTB1[2].T = dRmjd + gdStepsize_used;
				_nutCal(sTB1[0].T, 0.0, &sTB1[0].Pi, &sTB1[0].Ep);
				_nutCal(sTB1[1].T, 0.0, &sTB1[1].Pi, &sTB1[1].Ep);
				_nutCal(sTB1[2].T, 0.0, &sTB1[2].Pi, &sTB1[2].Ep);
			}
			else if (dRmjd < sTB1[0].T)
			{
				sTB1[2] = sTB1[1];
				sTB1[1] = sTB1[0];
				sTB1[0].T = sTB1[0].T - gdStepsize_used;
				_nutCal(sTB1[0].T, 0.0, &sTB1[0].Pi, &sTB1[0].Ep);
			}
			else
			{
				sTB1[0] = sTB1[1];
				sTB1[1] = sTB1[2];
				sTB1[2].T = sTB1[2].T + gdStepsize_used;
				_nutCal(sTB1[2].T, 0, &sTB1[2].Pi, &sTB1[2].Ep);
			}

		}
		double dTemp1[3];
		for (int i = 0; i < 3; i++)
		{
			dTemp1[i] = sTB1[i].T;
		}
		double dTemp2[3];
		for (int i = 0; i < 3; i++)
		{
			dTemp2[i] = sTB1[i].Pi;
		}
		double dTemp3[3];
		for (int i = 0; i < 3; i++)
		{
			dTemp3[i] = sTB1[i].Ep;
		}

		*dpsi = _interpolation(2, 3, dTemp1, dTemp2, dRmjd);
		*deps = _interpolation(2, 3, dTemp1, dTemp3, dRmjd);

	}

	/***********************
	*FunctionName：r_interpolation
	*Function：interpolation of linear of 2-order polynomials
	*InPut： iOrder 1 or 2 for linear or quadratic
	iPoint number of points in the input array (x,y) in case of long table (npoint > norder+1)：
	*OutPut：pdX
	pdX
	dXin
	*Return：
	*Other：
	***********************/
	double t_gtrs2crs::_interpolation(const int& iOrder, const int& iPoint, double *pdX, double *pdY, const double& dXin)
	{
		const string strProgname = "_interpolation";

		int i1;
		int i2;
		int i3;
		double dX;
		double dIntv;
		string strTemp1;
		string strTemp2;

		dIntv = (*(pdX + 1)) - (*(pdX));
		//qyj
		if ((dXin < (*(pdX)-0.1*dIntv)) || (dXin > (*(pdX + iPoint - 1) + 0.1*dIntv)))
		{
			cout << " ERROR: input variable out of the table (tbeg,tend,tinput): "
				<< setw(10) << (*pdX) << *(pdX + iPoint - 1) << dXin << endl;
			
		}
		
		if (iPoint == iOrder + 1)
		{
			i1 = 1;
		}
		else if (iPoint > iOrder + 1)
		{
			i1 = (int)((dXin - *pdX) / dIntv) + 1;
			if (i1 < 1)
			{
				i1 = 1;
			}
			if ((iPoint - iOrder) < i1)
			{
				i1 = iPoint - iOrder;
			}

		}
		else
		{
			cout<< " ERROR: table is enough long for interpolation: "
				<< " npoint = " << iPoint
				<< " norder = " << iOrder
				<< endl;
			throw(" ERROR: table is enough long for interpolation");
		}

		
		double dResult;
		if (iOrder == 1)
		{
			i2 = i1 + 1;
			dX = (dXin - *(pdX + i1 - 1)) / (*(pdX + i2 - 1) - *(pdX + i1 - 1));
			dResult = *(pdY + i1 - 1) + dX * (*(pdY + i2 - 1) - *(pdY + i1 - 1));
			return dResult;
		}
		else if (iOrder == 2)
		{
			i2 = i1 + 1;
			i3 = i2 + 1;
			dX = (dXin - *(pdX + i2 - 1)) / (*(pdX + i2 - 1) - *(pdX + i1 - 1));
			dResult = *(pdY + i2 - 1) + dX * (*(pdY + i3 - 1) - *(pdY + i1 - 1)) * 0.5
				+ dX * dX * ((*(pdY + i3 - 1) + *(pdY + i1 - 1)) * 0.5 - *(pdY + i2 - 1));
			return dResult;
		}
		else
		{
			cout << "***ERROR: norder > 2 is not supported ";
			throw("***ERROR: norder > 2 is not supported ");
		}
	}

	// calculate rotation matrix of Precession and nutation 
	void t_gtrs2crs::_process2000(const double& dRmjd, const double& dpsi, const double& deps)
	{

		Matrix dMathlp(0.0, 3, 3);
		Matrix dNmat(0.0, 3, 3);
		Matrix dPmat(0.0, 3, 3);

		vector<Matrix> rot; rot.reserve(5); // add by zhenghj

		//! Interval between fundamental epoch J2000.0 and given date in centuries
		double dT = (dRmjd - 51544.5) / 36525.e0;

		//! IAU 1980 mean obliquity of date and  Precession angles (Lieske et al. 1977) ( equ(30))
		double dEpsa80 = dEps0 + (-46.8150 + (-0.00059 + (0.001813) * dT) * dT) * dT / RAD2SEC;
		double dPsia77 = (5038.7784 + (-1.07259 + (-0.001147) * dT) * dT) * dT / RAD2SEC;
		double dOma77 = dEps0 + ((0.05127 + (-0.007726) * dT) * dT) * dT / RAD2SEC;
		double dChia = (10.5526 + (-2.38064 + (-0.001125) * dT) * dT) * dT / RAD2SEC;

		//    	 ! Apply IAU 2000A precession corrections.
		double dPsia = dPsia77 + dPrecor * dT;
		double dOma = dOma77 + dOblcor * dT;
		_epsa = dEpsa80 + dOblcor * dT;

		calcProcMat(false, 1, _epsa + deps, rot);
		dNmat << rot[0];

		calcProcMat(false, 3, dpsi, rot);
		dMathlp << rot[0];
		dNmat = dMathlp * dNmat;

		calcProcMat(false, 1, -_epsa, rot);
		dMathlp << rot[0];	
		dNmat = dMathlp * dNmat;

		calcProcMat(false, 3, -dChia, rot);
		dPmat << rot[0];

		calcProcMat(false, 1, dOma, rot);
		dMathlp << rot[0];		
		dPmat = dMathlp * dPmat;

		calcProcMat(false, 3, dPsia, rot);
		dMathlp << rot[0];	
		dPmat = dMathlp * dPmat;

		calcProcMat(false, 1, -dEps0, rot);
		dMathlp << rot[0];
		dPmat = dMathlp * dPmat;

		_qmat = dPmat * dNmat;
		calcProcMat(false, 1, dEpsbi, rot);
		dMathlp << rot[0];
		_qmat = dMathlp * _qmat;

		calcProcMat(false, 2, -dPsibi * sin(dEps0), rot);
		dMathlp << rot[0];

		_qmat = dMathlp * _qmat;

		calcProcMat(false, 3, -dRa0, rot);
		dMathlp << rot[0];
		_qmat = dMathlp * _qmat;
	}

	double t_gtrs2crs::_sp2000(const double& dDATE1, const double& dDATE2)
	{

		double dT = ((dDATE1 - dJ0) + dDATE2) / dJC;
		double dTemp = (-47) * 1.0e-6 * dT * dAS2R;
		return dTemp;

	}

	//The result is the Earth Rotation Angle (radians), in the range 0 to 2pi.
	double t_gtrs2crs::_era2000(const double& dJ1, const double& dJ2)
	{	
		double d1 = (dJ1 < dJ2) ? dJ1 : dJ2;
		double d2 = (dJ1 < dJ2) ? dJ2 : dJ1;
		double dT = d1 + (d2 - dJ0);
		double dF = fmod((d1 + 0.5), 1.0) + fmod(d2, 1.0);
		
		return _iau_anp(d2PI * (dF + 0.7790572732640 + 0.00273781191135448 * dT));
	}

	double t_gtrs2crs::_iau_anp(const double& dA)
	{
		double dW = fmod(dA, d2PI);
		if (dW < 0)
		{
			dW = dW + d2PI;
		}
		return dW;
	}

	//calculate Greenwich sidereal time, IAU 2000
	double t_gtrs2crs::_gst2000(const double& dRmjd, const double& era, const double& dpsi)
	{
		double dT = (dRmjd - 51544.5) / 36525.0;
		_gmst = era + (0.014506 + (4612.15739966 + (1.39667721 +
			(-0.00009344 + 0.00001882 * dT) * dT) * dT) * dT) / RAD2SEC;
		double dTemp = _gmst + dpsi * cos(_epsa) + _eect2000(dRmjd);
		return dTemp;
	}

	double t_gtrs2crs::_iau_anpm(const double& dA)
	{
		double dW = fmod(dA, d2PI);
		if (fabs(dW) >= dPI)
		{
			if (dA > 0)
			{
				dW = dW - d2PI;
			}
			else
			{
				dW = dW + d2PI;
			}
		}

		return dW;
	}

	//Equation of the equinoxes complementary terms, consistent with IAU 2000 resolutions(return complementary terms (radians))
	double t_gtrs2crs::_eect2000(const double& dRmjd)
	{
		//    	 *  Time since J2000, in Julian centuries
		double dT;
		//    	 *  Miscellaneous
		int I;
		int J;
		double dA;
		double dS0;
		double dS1;
		//    	 *  Fundamental arguments
		double dFA[14];
		//    	 *  -----------------------------------------
		//    	 *  The series for the EE complementary terms
		//    	 *  -----------------------------------------
		//
		//    	 *  Number of terms in the series
		//    	       INTEGER NE0, NE1
		//    	       PARAMETER ( NE0=  33, NE1=  1 )
		int const iNE0 = 33;
		int const iNE1 = 1;

		//    	 *  Coefficients of l,l',F,D,Om,LMe,LVe,LE,LMa,LJu,LSa,LU,LN,pA
		//    	       INTEGER KE0 ( 14, NE0 ),
		//    	      :        KE1 ( 14, NE1 )
		int iKE0[14][iNE0];

		//
		//    	 *  Sine and cosine coefficients
		//    	       DOUBLE PRECISION SE0 ( 2, NE0 ),
		//    	      :                 SE1 ( 2, NE1 )
		double dSE0[2][iNE0];

		//    	 *  Argument coefficients for t^0
		int iTemp1[10][14] =
		{
			0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 2, -2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 2, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 2, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0
		};
		for (int i = 0; i < 14; i++)
		{
			for (int j = 0; j < 10; j++)
			{
				iKE0[i][j] = iTemp1[j][i];
			}
		}

		int iTemp2[10][14] =
		{
			1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 1, 2, -2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 1, 2, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 4, -4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 1, -1, 1, 0, -8, 12, 0, 0, 0, 0, 0, 0,
			0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			1, 0, 2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			1, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0
		};
		for (int i = 0; i < 14; i++)
		{
			for (int j = 0; j < 10; j++)
			{
				iKE0[i][j + 10] = iTemp2[j][i];
			}
		}

		int iTemp3[10][14] =
		{
			0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 1, -2, 2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 1, -2, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 8, -13, 0, 0, 0, 0, 0, -1,
			0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			2, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			1, 0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 1, 2, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			1, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 4, -2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0
		};
		for (int i = 0; i < 14; i++)
		{
			for (int j = 0; j < 10; j++)
			{
				iKE0[i][j + 20] = iTemp3[j][i];
			}
		}


		int iTemp4[3][14] =
		{
			0, 0, 2, -2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			1, 0, -2, 0, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			1, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0
		};
		for (int i = 0; i < 14; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				iKE0[i][j + 30] = iTemp4[j][i];
			}
		}

		//    	 *  Argument coefficients for t^1
		int iKE1[14][iNE1] = { 0, 0, 0, 0, 1, 0, 0,
			0, 0, 0, 0, 0, 0, 0 };
		//
		//    	 *  Sine and cosine coefficients for t^0
		double dTemp1[10][2] =
		{

			+2640.96e-6, -0.39e-6,
			+63.52e-6, -0.02e-6,
			+11.75e-6, +0.01e-6,
			+11.21e-6, +0.01e-6,
			-4.55e-6, +0.00e-6,
			+2.02e-6, +0.00e-6,
			+1.98e-6, +0.00e-6,
			-1.72e-6, +0.00e-6,
			-1.41e-6, -0.01e-6,
			-1.26e-6, -0.01e-6
		};

		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 10; j++)
			{
				dSE0[i][j] = dTemp1[j][i];
			}
		}

		double dTemp2[10][2] =
		{

			-0.63e-6, +0.00e-6,
			-0.63e-6, +0.00e-6,
			+0.46e-6, +0.00e-6,
			+0.45e-6, +0.00e-6,
			+0.36e-6, +0.00e-6,
			-0.24e-6, -0.12e-6,
			+0.32e-6, +0.00e-6,
			+0.28e-6, +0.00e-6,
			+0.27e-6, +0.00e-6,
			+0.26e-6, +0.00e-6
		};

		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 10; j++)
			{
				dSE0[i][j + 10] = dTemp2[j][i];
			}
		}

		double dTemp3[10][2] =
		{

			-0.21e-6, +0.00e-6,
			+0.19e-6, +0.00e-6,
			+0.18e-6, +0.00e-6,
			-0.10e-6, +0.05e-6,
			+0.15e-6, +0.00e-6,
			-0.14e-6, +0.00e-6,
			+0.14e-6, +0.00e-6,
			-0.14e-6, +0.00e-6,
			+0.14e-6, +0.00e-6,
			+0.13e-6, +0.00e-6
		};

		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 10; j++)
			{
				dSE0[i][j + 20] = dTemp3[j][i];
			}
		}

		double dTemp4[3][2] =
		{

			-0.11e-6, +0.00e-6,
			+0.11e-6, +0.00e-6,
			+0.11e-6, +0.00e-6
		};
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				dSE0[i][j + 30] = dTemp4[j][i];
			}
		}
		//    	 *  Sine and cosine coefficients for t^1
		double dSE1[2][iNE1] = { -0.87*pow(10, -6), +0.00*pow(10, -6) };

		//    	 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		//
		//    	 *  Interval between fundamental epoch J2000.0 and current date (JC).
		dT = (dRmjd - 51544.5) / dJC;
		//    	 *  Fundamental Arguments (from IERS Conventions 2000)
		//
		//    	 *  Mean Anomaly of the Moon.

		dFA[0] = _iau_anpm((485868.249036 +
			(715923.2178 +
			(31.8792 +
				(0.051635 +
				(-0.00024470)
					* dT) * dT) * dT) * dT) * dAS2R
			+ fmod(1325 * dT, 1.0) * d2PI);
		//    	 *  Mean Anomaly of the Sun.
		dFA[1] = _iau_anpm((1287104.793048 +
			(1292581.0481 +
			(-0.5532 +
				(+0.000136 +
				(-0.00001149)
					* dT) * dT) * dT) * dT) * dAS2R
			+ fmod(99 * dT, 1.0) * d2PI);
		//    	 *  Mean Longitude of the Moon minus Mean Longitude of the Ascending
		//    	 *  Node of the Moon.

		dFA[2] = _iau_anpm((335779.526232 +
			(295262.8478 +
			(-12.7512 +
				(-0.001037 +
				(0.00000417)
					* dT) * dT) * dT) * dT) * dAS2R
			+ fmod(1342 * dT, 1.0) * d2PI);


		dFA[3] = _iau_anpm((1072260.703692 +
			(1105601.2090 +
			(-6.3706 +
				(0.006593 +
				(-0.00003169)
					* dT) * dT) * dT) * dT) * dAS2R
			+ fmod(1236 * dT, 1.0) * d2PI);

		//    	 *  Mean Longitude of the Ascending Node of the Moon.

		dFA[4] = _iau_anpm((450160.398036 +
			(-482890.5431 +
			(7.4722 +
				(0.007702 +
				(-0.00005939)
					* dT) * dT) * dT) * dT) * dAS2R
			+ fmod(-5 * dT, 1.0) * d2PI);


		dFA[5] = _iau_anpm(4.402608842 + 2608.7903141574 * dT);
		dFA[6] = _iau_anpm(3.176146697 + 1021.3285546211 * dT);
		dFA[7] = _iau_anpm(1.753470314 + 628.3075849991 * dT);
		dFA[8] = _iau_anpm(6.203480913 + 334.0612426700 * dT);
		dFA[9] = _iau_anpm(0.599546497 + 52.9690962641 * dT);
		dFA[10] = _iau_anpm(0.874016757 + 21.3299104960 * dT);
		dFA[11] = _iau_anpm(5.481293872 + 7.4781598567 * dT);
		dFA[12] = _iau_anpm(5.311886287 + 3.8133035638 * dT);
		dFA[13] = (0.024381750 + 0.00000538691 * dT) * dT;

		//
		//    	 *  Evaluate the EE complementary terms.

		dS0 = 0.0;
		dS1 = 0.0;


		for (I = iNE0 - 1; I >= 0; I--)
		{
			dA = 0.0;
			for (J = 0; J < 14; J++)
			{
				dA = dA + iKE0[J][I] * dFA[J];
			}
			dS0 = dS0 + (dSE0[0][I] * sin(dA) + dSE0[1][I] * cos(dA));
		}

		for (I = iNE1 - 1; I >= 0; I--)
		{
			dA = 0.0;
			for (J = 0; J < 14; J++)
			{
				dA = dA + iKE1[J][I] * dFA[J];
			}
			dS1 = dS1 + (dSE1[0][I] * sin(dA) + dSE1[1][I] * cos(dA));
		}


		return (dS0 + dS1 * dT) * dAS2R;
	}

	void  t_gtrs2crs::_nutCal(const double& dATE1, const double& dATE2, double *dPSI, double *dEPS)
	{
		//todo 变量是否考虑维护成本 ？？
		double dT, dEL, dELP, dF, dD, dOM, dARG, dDP, dDE, dSARG, dCARG,
			dDPSILS, dDEPSLS,
			dAL, dALSU, dAF, dAD, dAOM, dALME, dALVE, dALEA, dALMA,
			dALJU, dALSA, dALUR, dALNE, dAPA, dDPSIPL, dDEPSPL;
		int I, J;

		/*** Luni-Solar nutation model***/

		const int iNLS = 678;
		//int i;
		int(*iNALS)[iNLS] = new int[5][iNLS];
		double(*dCLS)[iNLS] = new double[6][iNLS];//modified by zhuyt
		//int iNALS[5][iNLS];
		//double dCLS[6][iNLS];


		/*** Planetary terms***/

		const int iNPL = 687;
		int(*iNAPL)[iNPL] = new int[14][iNPL];
		int(*iICPL)[iNPL] = new int[4][iNPL];//modified by zhuyt
		//int iNAPL[14][iNPL];
		//int iICPL[4][iNPL];

		/***Tables of argument and term coefficients***/
		/*  Luni-Solar argument multipliers:
			L     L'    F     D     Om  */
			//NALST is the transposition of NALS
		int iNALST[iNLS][5] = {
				{ 0, 0, 0, 0, 1 },
				{ 0, 0, 2, -2, 2 },
				{ 0, 0, 2, 0, 2 },
				{ 0, 0, 0, 0, 2 },
				{ 0, 1, 0, 0, 0 },
				{ 0, 1, 2, -2, 2 },
				{ 1, 0, 0, 0, 0 },
				{ 0, 0, 2, 0, 1 },
				{ 1, 0, 2, 0, 2 },
				{ 0, -1, 2, -2, 2 },
				{ 0, 0, 2, -2, 1 },
				{ -1, 0, 2, 0, 2 },
				{ -1, 0, 0, 2, 0 },
				{ 1, 0, 0, 0, 1 },
				{ -1, 0, 0, 0, 1 },
				{ -1, 0, 2, 2, 2 },
				{ 1, 0, 2, 0, 1 },
				{ -2, 0, 2, 0, 1 },
				{ 0, 0, 0, 2, 0 },
				{ 0, 0, 2, 2, 2 },
				{ 0, -2, 2, -2, 2 },
				{ -2, 0, 0, 2, 0 },
				{ 2, 0, 2, 0, 2 },
				{ 1, 0, 2, -2, 2 },
				{ -1, 0, 2, 0, 1 },
				{ 2, 0, 0, 0, 0 },
				{ 0, 0, 2, 0, 0 },
				{ 0, 1, 0, 0, 1 },
				{ -1, 0, 0, 2, 1 },
				{ 0, 2, 2, -2, 2 },
				{ 0, 0, -2, 2, 0 },
				{ 1, 0, 0, -2, 1 },
				{ 0, -1, 0, 0, 1 },
				{ -1, 0, 2, 2, 1 },
				{ 0, 2, 0, 0, 0 },
				{ 1, 0, 2, 2, 2 },
				{ -2, 0, 2, 0, 0 },
				{ 0, 1, 2, 0, 2 },
				{ 0, 0, 2, 2, 1 },
				{ 0, -1, 2, 0, 2 },
				{ 0, 0, 0, 2, 1 },
				{ 1, 0, 2, -2, 1 },
				{ 2, 0, 2, -2, 2 },
				{ -2, 0, 0, 2, 1 },
				{ 2, 0, 2, 0, 1 },
				{ 0, -1, 2, -2, 1 },
				{ 0, 0, 0, -2, 1 },
				{ -1, -1, 0, 2, 0 },
				{ 2, 0, 0, -2, 1 },
				{ 1, 0, 0, 2, 0 },
				{ 0, 1, 2, -2, 1 },
				{ 1, -1, 0, 0, 0 },
				{ -2, 0, 2, 0, 2 },
				{ 3, 0, 2, 0, 2 },
				{ 0, -1, 0, 2, 0 },
				{ 1, -1, 2, 0, 2 },
				{ 0, 0, 0, 1, 0 },
				{ -1, -1, 2, 2, 2 },
				{ -1, 0, 2, 0, 0 },
				{ 0, -1, 2, 2, 2 },
				{ -2, 0, 0, 0, 1 },
				{ 1, 1, 2, 0, 2 },
				{ 2, 0, 0, 0, 1 },
				{ -1, 1, 0, 1, 0 },
				{ 1, 1, 0, 0, 0 },
				{ 1, 0, 2, 0, 0 },
				{ -1, 0, 2, -2, 1 },
				{ 1, 0, 0, 0, 2 },
				{ -1, 0, 0, 1, 0 },
				{ 0, 0, 2, 1, 2 },
				{ -1, 0, 2, 4, 2 },
				{ -1, 1, 0, 1, 1 },
				{ 0, -2, 2, -2, 1 },
				{ 1, 0, 2, 2, 1 },
				{ -2, 0, 2, 2, 2 },
				{ -1, 0, 0, 0, 2 },
				{ 1, 1, 2, -2, 2 },
				{ -2, 0, 2, 4, 2 },
				{ -1, 0, 4, 0, 2 },
				{ 2, 0, 2, -2, 1 },
				{ 2, 0, 2, 2, 2 },
				{ 1, 0, 0, 2, 1 },
				{ 3, 0, 0, 0, 0 },
				{ 3, 0, 2, -2, 2 },
				{ 0, 0, 4, -2, 2 },
				{ 0, 1, 2, 0, 1 },
				{ 0, 0, -2, 2, 1 },
				{ 0, 0, 2, -2, 3 },
				{ -1, 0, 0, 4, 0 },
				{ 2, 0, -2, 0, 1 },
				{ -2, 0, 0, 4, 0 },
				{ -1, -1, 0, 2, 1 },
				{ -1, 0, 0, 1, 1 },
				{ 0, 1, 0, 0, 2 },
				{ 0, 0, -2, 0, 1 },
				{ 0, -1, 2, 0, 1 },
				{ 0, 0, 2, -1, 2 },
				{ 0, 0, 2, 4, 2 },
				{ -2, -1, 0, 2, 0 },
				{ 1, 1, 0, -2, 1 },
				{ -1, 1, 0, 2, 0 },
				{ -1, 1, 0, 1, 2 },
				{ 1, -1, 0, 0, 1 },
				{ 1, -1, 2, 2, 2 },
				{ -1, 1, 2, 2, 2 },
				{ 3, 0, 2, 0, 1 },
				{ 0, 1, -2, 2, 0 },
				{ -1, 0, 0, -2, 1 },
				{ 0, 1, 2, 2, 2 },
				{ -1, -1, 2, 2, 1 },
				{ 0, -1, 0, 0, 2 },
				{ 1, 0, 2, -4, 1 },
				{ -1, 0, -2, 2, 0 },
				{ 0, -1, 2, 2, 1 },
				{ 2, -1, 2, 0, 2 },
				{ 0, 0, 0, 2, 2 },
				{ 1, -1, 2, 0, 1 },
				{ -1, 1, 2, 0, 2 },
				{ 0, 1, 0, 2, 0 },
				{ 0, -1, -2, 2, 0 },
				{ 0, 3, 2, -2, 2 },
				{ 0, 0, 0, 1, 1 },
				{ -1, 0, 2, 2, 0 },
				{ 2, 1, 2, 0, 2 },
				{ 1, 1, 0, 0, 1 },
				{ 1, 1, 2, 0, 1 },
				{ 2, 0, 0, 2, 0 },
				{ 1, 0, -2, 2, 0 },
				{ -1, 0, 0, 2, 2 },
				{ 0, 1, 0, 1, 0 },
				{ 0, 1, 0, -2, 1 },
				{ -1, 0, 2, -2, 2 },
				{ 0, 0, 0, -1, 1 },
				{ -1, 1, 0, 0, 1 },
				{ 1, 0, 2, -1, 2 },
				{ 1, -1, 0, 2, 0 },
				{ 0, 0, 0, 4, 0 },
				{ 1, 0, 2, 1, 2 },
				{ 0, 0, 2, 1, 1 },
				{ 1, 0, 0, -2, 2 },
				{ -1, 0, 2, 4, 1 },
				{ 1, 0, -2, 0, 1 },
				{ 1, 1, 2, -2, 1 },
				{ 0, 0, 2, 2, 0 },
				{ -1, 0, 2, -1, 1 },
				{ -2, 0, 2, 2, 1 },
				{ 4, 0, 2, 0, 2 },
				{ 2, -1, 0, 0, 0 },
				{ 2, 1, 2, -2, 2 },
				{ 0, 1, 2, 1, 2 },
				{ 1, 0, 4, -2, 2 },
				{ -1, -1, 0, 0, 1 },
				{ 0, 1, 0, 2, 1 },
				{ -2, 0, 2, 4, 1 },
				{ 2, 0, 2, 0, 0 },
				{ 1, 0, 0, 1, 0 },
				{ -1, 0, 0, 4, 1 },
				{ -1, 0, 4, 0, 1 },
				{ 2, 0, 2, 2, 1 },
				{ 0, 0, 2, -3, 2 },
				{ -1, -2, 0, 2, 0 },
				{ 2, 1, 0, 0, 0 },
				{ 0, 0, 4, 0, 2 },
				{ 0, 0, 0, 0, 3 },
				{ 0, 3, 0, 0, 0 },
				{ 0, 0, 2, -4, 1 },
				{ 0, -1, 0, 2, 1 },
				{ 0, 0, 0, 4, 1 },
				{ -1, -1, 2, 4, 2 },
				{ 1, 0, 2, 4, 2 },
				{ -2, 2, 0, 2, 0 },
				{ -2, -1, 2, 0, 1 },
				{ -2, 0, 0, 2, 2 },
				{ -1, -1, 2, 0, 2 },
				{ 0, 0, 4, -2, 1 },
				{ 3, 0, 2, -2, 1 },
				{ -2, -1, 0, 2, 1 },
				{ 1, 0, 0, -1, 1 },
				{ 0, -2, 0, 2, 0 },
				{ -2, 0, 0, 4, 1 },
				{ -3, 0, 0, 0, 1 },
				{ 1, 1, 2, 2, 2 },
				{ 0, 0, 2, 4, 1 },
				{ 3, 0, 2, 2, 2 },
				{ -1, 1, 2, -2, 1 },
				{ 2, 0, 0, -4, 1 },
				{ 0, 0, 0, -2, 2 },
				{ 2, 0, 2, -4, 1 },
				{ -1, 1, 0, 2, 1 },
				{ 0, 0, 2, -1, 1 },
				{ 0, -2, 2, 2, 2 },
				{ 2, 0, 0, 2, 1 },
				{ 4, 0, 2, -2, 2 },
				{ 2, 0, 0, -2, 2 },
				{ 0, 2, 0, 0, 1 },
				{ 1, 0, 0, -4, 1 },
				{ 0, 2, 2, -2, 1 },
				{ -3, 0, 0, 4, 0 },
				{ -1, 1, 2, 0, 1 },
				{ -1, -1, 0, 4, 0 },
				{ -1, -2, 2, 2, 2 },
				{ -2, -1, 2, 4, 2 },
				{ 1, -1, 2, 2, 1 },
				{ -2, 1, 0, 2, 0 },
				{ -2, 1, 2, 0, 1 },
				{ 2, 1, 0, -2, 1 },
				{ -3, 0, 2, 0, 1 },
				{ -2, 0, 2, -2, 1 },
				{ -1, 1, 0, 2, 2 },
				{ 0, -1, 2, -1, 2 },
				{ -1, 0, 4, -2, 2 },
				{ 0, -2, 2, 0, 2 },
				{ -1, 0, 2, 1, 2 },
				{ 2, 0, 0, 0, 2 },
				{ 0, 0, 2, 0, 3 },
				{ -2, 0, 4, 0, 2 },
				{ -1, 0, -2, 0, 1 },
				{ -1, 1, 2, 2, 1 },
				{ 3, 0, 0, 0, 1 },
				{ -1, 0, 2, 3, 2 },
				{ 2, -1, 2, 0, 1 },
				{ 0, 1, 2, 2, 1 },
				{ 0, -1, 2, 4, 2 },
				{ 2, -1, 2, 2, 2 },
				{ 0, 2, -2, 2, 0 },
				{ -1, -1, 2, -1, 1 },
				{ 0, -2, 0, 0, 1 },
				{ 1, 0, 2, -4, 2 },
				{ 1, -1, 0, -2, 1 },
				{ -1, -1, 2, 0, 1 },
				{ 1, -1, 2, -2, 2 },
				{ -2, -1, 0, 4, 0 },
				{ -1, 0, 0, 3, 0 },
				{ -2, -1, 2, 2, 2 },
				{ 0, 2, 2, 0, 2 },
				{ 1, 1, 0, 2, 0 },
				{ 2, 0, 2, -1, 2 },
				{ 1, 0, 2, 1, 1 },
				{ 4, 0, 0, 0, 0 },
				{ 2, 1, 2, 0, 1 },
				{ 3, -1, 2, 0, 2 },
				{ -2, 2, 0, 2, 1 },
				{ 1, 0, 2, -3, 1 },
				{ 1, 1, 2, -4, 1 },
				{ -1, -1, 2, -2, 1 },
				{ 0, -1, 0, -1, 1 },
				{ 0, -1, 0, -2, 1 },
				{ -2, 0, 0, 0, 2 },
				{ -2, 0, -2, 2, 0 },
				{ -1, 0, -2, 4, 0 },
				{ 1, -2, 0, 0, 0 },
				{ 0, 1, 0, 1, 1 },
				{ -1, 2, 0, 2, 0 },
				{ 1, -1, 2, -2, 1 },
				{ 1, 2, 2, -2, 2 },
				{ 2, -1, 2, -2, 2 },
				{ 1, 0, 2, -1, 1 },
				{ 2, 1, 2, -2, 1 },
				{ -2, 0, 0, -2, 1 },
				{ 1, -2, 2, 0, 2 },
				{ 0, 1, 2, 1, 1 },
				{ 1, 0, 4, -2, 1 },
				{ -2, 0, 4, 2, 2 },
				{ 1, 1, 2, 1, 2 },
				{ 1, 0, 0, 4, 0 },
				{ 1, 0, 2, 2, 0 },
				{ 2, 0, 2, 1, 2 },
				{ 3, 1, 2, 0, 2 },
				{ 4, 0, 2, 0, 1 },
				{ -2, -1, 2, 0, 0 },
				{ 0, 1, -2, 2, 1 },
				{ 1, 0, -2, 1, 0 },
				{ 0, -1, -2, 2, 1 },
				{ 2, -1, 0, -2, 1 },
				{ -1, 0, 2, -1, 2 },
				{ 1, 0, 2, -3, 2 },
				{ 0, 1, 2, -2, 3 },
				{ 0, 0, 2, -3, 1 },
				{ -1, 0, -2, 2, 1 },
				{ 0, 0, 2, -4, 2 },
				{ -2, 1, 0, 0, 1 },
				{ -1, 0, 0, -1, 1 },
				{ 2, 0, 2, -4, 2 },
				{ 0, 0, 4, -4, 4 },
				{ 0, 0, 4, -4, 2 },
				{ -1, -2, 0, 2, 1 },
				{ -2, 0, 0, 3, 0 },
				{ 1, 0, -2, 2, 1 },
				{ -3, 0, 2, 2, 2 },
				{ -3, 0, 2, 2, 1 },
				{ -2, 0, 2, 2, 0 },
				{ 2, -1, 0, 0, 1 },
				{ -2, 1, 2, 2, 2 },
				{ 1, 1, 0, 1, 0 },
				{ 0, 1, 4, -2, 2 },
				{ -1, 1, 0, -2, 1 },
				{ 0, 0, 0, -4, 1 },
				{ 1, -1, 0, 2, 1 },
				{ 1, 1, 0, 2, 1 },
				{ -1, 2, 2, 2, 2 },
				{ 3, 1, 2, -2, 2 },
				{ 0, -1, 0, 4, 0 },
				{ 2, -1, 0, 2, 0 },
				{ 0, 0, 4, 0, 1 },
				{ 2, 0, 4, -2, 2 },
				{ -1, -1, 2, 4, 1 },
				{ 1, 0, 0, 4, 1 },
				{ 1, -2, 2, 2, 2 },
				{ 0, 0, 2, 3, 2 },
				{ -1, 1, 2, 4, 2 },
				{ 3, 0, 0, 2, 0 },
				{ -1, 0, 4, 2, 2 },
				{ 1, 1, 2, 2, 1 },
				{ -2, 0, 2, 6, 2 },
				{ 2, 1, 2, 2, 2 },
				{ -1, 0, 2, 6, 2 },
				{ 1, 0, 2, 4, 1 },
				{ 2, 0, 2, 4, 2 },
				{ 1, 1, -2, 1, 0 },
				{ -3, 1, 2, 1, 2 },
				{ 2, 0, -2, 0, 2 },
				{ -1, 0, 0, 1, 2 },
				{ -4, 0, 2, 2, 1 },
				{ -1, -1, 0, 1, 0 },
				{ 0, 0, -2, 2, 2 },
				{ 1, 0, 0, -1, 2 },
				{ 0, -1, 2, -2, 3 },
				{ -2, 1, 2, 0, 0 },
				{ 0, 0, 2, -2, 4 },
				{ -2, -2, 0, 2, 0 },
				{ -2, 0, -2, 4, 0 },
				{ 0, -2, -2, 2, 0 },
				{ 1, 2, 0, -2, 1 },
				{ 3, 0, 0, -4, 1 },
				{ -1, 1, 2, -2, 2 },
				{ 1, -1, 2, -4, 1 },
				{ 1, 1, 0, -2, 2 },
				{ -3, 0, 2, 0, 0 },
				{ -3, 0, 2, 0, 2 },
				{ -2, 0, 0, 1, 0 },
				{ 0, 0, -2, 1, 0 },
				{ -3, 0, 0, 2, 1 },
				{ -1, -1, -2, 2, 0 },
				{ 0, 1, 2, -4, 1 },
				{ 2, 1, 0, -4, 1 },
				{ 0, 2, 0, -2, 1 },
				{ 1, 0, 0, -3, 1 },
				{ -2, 0, 2, -2, 2 },
				{ -2, -1, 0, 0, 1 },
				{ -4, 0, 0, 2, 0 },
				{ 1, 1, 0, -4, 1 },
				{ -1, 0, 2, -4, 1 },
				{ 0, 0, 4, -4, 1 },
				{ 0, 3, 2, -2, 2 },
				{ -3, -1, 0, 4, 0 },
				{ -3, 0, 0, 4, 1 },
				{ 1, -1, -2, 2, 0 },
				{ -1, -1, 0, 2, 2 },
				{ 1, -2, 0, 0, 1 },
				{ 1, -1, 0, 0, 2 },
				{ 0, 0, 0, 1, 2 },
				{ -1, -1, 2, 0, 0 },
				{ 1, -2, 2, -2, 2 },
				{ 0, -1, 2, -1, 1 },
				{ -1, 0, 2, 0, 3 },
				{ 1, 1, 0, 0, 2 },
				{ -1, 1, 2, 0, 0 },
				{ 1, 2, 0, 0, 0 },
				{ -1, 2, 2, 0, 2 },
				{ -1, 0, 4, -2, 1 },
				{ 3, 0, 2, -4, 2 },
				{ 1, 2, 2, -2, 1 },
				{ 1, 0, 4, -4, 2 },
				{ -2, -1, 0, 4, 1 },
				{ 0, -1, 0, 2, 2 },
				{ -2, 1, 0, 4, 0 },
				{ -2, -1, 2, 2, 1 },
				{ 2, 0, -2, 2, 0 },
				{ 1, 0, 0, 1, 1 },
				{ 0, 1, 0, 2, 2 },
				{ 1, -1, 2, -1, 2 },
				{ -2, 0, 4, 0, 1 },
				{ 2, 1, 0, 0, 1 },
				{ 0, 1, 2, 0, 0 },
				{ 0, -1, 4, -2, 2 },
				{ 0, 0, 4, -2, 4 },
				{ 0, 2, 2, 0, 1 },
				{ -3, 0, 0, 6, 0 },
				{ -1, -1, 0, 4, 1 },
				{ 1, -2, 0, 2, 0 },
				{ -1, 0, 0, 4, 2 },
				{ -1, -2, 2, 2, 1 },
				{ -1, 0, 0, -2, 2 },
				{ 1, 0, -2, -2, 1 },
				{ 0, 0, -2, -2, 1 },
				{ -2, 0, -2, 0, 1 },
				{ 0, 0, 0, 3, 1 },
				{ 0, 0, 0, 3, 0 },
				{ -1, 1, 0, 4, 0 },
				{ -1, -1, 2, 2, 0 },
				{ -2, 0, 2, 3, 2 },
				{ 1, 0, 0, 2, 2 },
				{ 0, -1, 2, 1, 2 },
				{ 3, -1, 0, 0, 0 },
				{ 2, 0, 0, 1, 0 },
				{ 1, -1, 2, 0, 0 },
				{ 0, 0, 2, 1, 0 },
				{ 1, 0, 2, 0, 3 },
				{ 3, 1, 0, 0, 0 },
				{ 3, -1, 2, -2, 2 },
				{ 2, 0, 2, -1, 1 },
				{ 1, 1, 2, 0, 0 },
				{ 0, 0, 4, -1, 2 },
				{ 1, 2, 2, 0, 2 },
				{ -2, 0, 0, 6, 0 },
				{ 0, -1, 0, 4, 1 },
				{ -2, -1, 2, 4, 1 },
				{ 0, -2, 2, 2, 1 },
				{ 0, -1, 2, 2, 0 },
				{ -1, 0, 2, 3, 1 },
				{ -2, 1, 2, 4, 2 },
				{ 2, 0, 0, 2, 2 },
				{ 2, -2, 2, 0, 2 },
				{ -1, 1, 2, 3, 2 },
				{ 3, 0, 2, -1, 2 },
				{ 4, 0, 2, -2, 1 },
				{ -1, 0, 0, 6, 0 },
				{ -1, -2, 2, 4, 2 },
				{ -3, 0, 2, 6, 2 },
				{ -1, 0, 2, 4, 0 },
				{ 3, 0, 0, 2, 1 },
				{ 3, -1, 2, 0, 1 },
				{ 3, 0, 2, 0, 0 },
				{ 1, 0, 4, 0, 2 },
				{ 5, 0, 2, -2, 2 },
				{ 0, -1, 2, 4, 1 },
				{ 2, -1, 2, 2, 1 },
				{ 0, 1, 2, 4, 2 },
				{ 1, -1, 2, 4, 2 },
				{ 3, -1, 2, 2, 2 },
				{ 3, 0, 2, 2, 1 },
				{ 5, 0, 2, 0, 2 },
				{ 0, 0, 2, 6, 2 },
				{ 4, 0, 2, 2, 2 },
				{ 0, -1, 1, -1, 1 },
				{ -1, 0, 1, 0, 3 },
				{ 0, -2, 2, -2, 3 },
				{ 1, 0, -1, 0, 1 },
				{ 2, -2, 0, -2, 1 },
				{ -1, 0, 1, 0, 2 },
				{ -1, 0, 1, 0, 1 },
				{ -1, -1, 2, -1, 2 },
				{ -2, 2, 0, 2, 2 },
				{ -1, 0, 1, 0, 0 },
				{ -4, 1, 2, 2, 2 },
				{ -3, 0, 2, 1, 1 },
				{ -2, -1, 2, 0, 2 },
				{ 1, 0, -2, 1, 1 },
				{ 2, -1, -2, 0, 1 },
				{ -4, 0, 2, 2, 0 },
				{ -3, 1, 0, 3, 0 },
				{ -1, 0, -1, 2, 0 },
				{ 0, -2, 0, 0, 2 },
				{ 0, -2, 0, 0, 2 },
				{ -3, 0, 0, 3, 0 },
				{ -2, -1, 0, 2, 2 },
				{ -1, 0, -2, 3, 0 },
				{ -4, 0, 0, 4, 0 },
				{ 2, 1, -2, 0, 1 },
				{ 2, -1, 0, -2, 2 },
				{ 0, 0, 1, -1, 0 },
				{ -1, 2, 0, 1, 0 },
				{ -2, 1, 2, 0, 2 },
				{ 1, 1, 0, -1, 1 },
				{ 1, 0, 1, -2, 1 },
				{ 0, 2, 0, 0, 2 },
				{ 1, -1, 2, -3, 1 },
				{ -1, 1, 2, -1, 1 },
				{ -2, 0, 4, -2, 2 },
				{ -2, 0, 4, -2, 1 },
				{ -2, -2, 0, 2, 1 },
				{ -2, 0, -2, 4, 0 },
				{ 1, 2, 2, -4, 1 },
				{ 1, 1, 2, -4, 2 },
				{ -1, 2, 2, -2, 1 },
				{ 2, 0, 0, -3, 1 },
				{ -1, 2, 0, 0, 1 },
				{ 0, 0, 0, -2, 0 },
				{ -1, -1, 2, -2, 2 },
				{ -1, 1, 0, 0, 2 },
				{ 0, 0, 0, -1, 2 },
				{ -2, 1, 0, 1, 0 },
				{ 1, -2, 0, -2, 1 },
				{ 1, 0, -2, 0, 2 },
				{ -3, 1, 0, 2, 0 },
				{ -1, 1, -2, 2, 0 },
				{ -1, -1, 0, 0, 2 },
				{ -3, 0, 0, 2, 0 },
				{ -3, -1, 0, 2, 0 },
				{ 2, 0, 2, -6, 1 },
				{ 0, 1, 2, -4, 2 },
				{ 2, 0, 0, -4, 2 },
				{ -2, 1, 2, -2, 1 },
				{ 0, -1, 2, -4, 1 },
				{ 0, 1, 0, -2, 2 },
				{ -1, 0, 0, -2, 0 },
				{ 2, 0, -2, -2, 1 },
				{ -4, 0, 2, 0, 1 },
				{ -1, -1, 0, -1, 1 },
				{ 0, 0, -2, 0, 2 },
				{ -3, 0, 0, 1, 0 },
				{ -1, 0, -2, 1, 0 },
				{ -2, 0, -2, 2, 1 },
				{ 0, 0, -4, 2, 0 },
				{ -2, -1, -2, 2, 0 },
				{ 1, 0, 2, -6, 1 },
				{ -1, 0, 2, -4, 2 },
				{ 1, 0, 0, -4, 2 },
				{ 2, 1, 2, -4, 2 },
				{ 2, 1, 2, -4, 1 },
				{ 0, 1, 4, -4, 4 },
				{ 0, 1, 4, -4, 2 },
				{ -1, -1, -2, 4, 0 },
				{ -1, -3, 0, 2, 0 },
				{ -1, 0, -2, 4, 1 },
				{ -2, -1, 0, 3, 0 },
				{ 0, 0, -2, 3, 0 },
				{ -2, 0, 0, 3, 1 },
				{ 0, -1, 0, 1, 0 },
				{ -3, 0, 2, 2, 0 },
				{ 1, 1, -2, 2, 0 },
				{ -1, 1, 0, 2, 2 },
				{ 1, -2, 2, -2, 1 },
				{ 0, 0, 1, 0, 2 },
				{ 0, 0, 1, 0, 1 },
				{ 0, 0, 1, 0, 0 },
				{ -1, 2, 0, 2, 1 },
				{ 0, 0, 2, 0, 2 },
				{ -2, 0, 2, 0, 2 },
				{ 2, 0, 0, -1, 1 },
				{ 3, 0, 0, -2, 1 },
				{ 1, 0, 2, -2, 3 },
				{ 1, 2, 0, 0, 1 },
				{ 2, 0, 2, -3, 2 },
				{ -1, 1, 4, -2, 2 },
				{ -2, -2, 0, 4, 0 },
				{ 0, -3, 0, 2, 0 },
				{ 0, 0, -2, 4, 0 },
				{ -1, -1, 0, 3, 0 },
				{ -2, 0, 0, 4, 2 },
				{ -1, 0, 0, 3, 1 },
				{ 2, -2, 0, 0, 0 },
				{ 1, -1, 0, 1, 0 },
				{ -1, 0, 0, 2, 0 },
				{ 0, -2, 2, 0, 1 },
				{ -1, 0, 1, 2, 1 },
				{ -1, 1, 0, 3, 0 },
				{ -1, -1, 2, 1, 2 },
				{ 0, -1, 2, 0, 0 },
				{ -2, 1, 2, 2, 1 },
				{ 2, -2, 2, -2, 2 },
				{ 1, 1, 0, 1, 1 },
				{ 1, 0, 1, 0, 1 },
				{ 1, 0, 1, 0, 0 },
				{ 0, 2, 0, 2, 0 },
				{ 2, -1, 2, -2, 1 },
				{ 0, -1, 4, -2, 1 },
				{ 0, 0, 4, -2, 3 },
				{ 0, 1, 4, -2, 1 },
				{ 4, 0, 2, -4, 2 },
				{ 2, 2, 2, -2, 2 },
				{ 2, 0, 4, -4, 2 },
				{ -1, -2, 0, 4, 0 },
				{ -1, -3, 2, 2, 2 },
				{ -3, 0, 2, 4, 2 },
				{ -3, 0, 2, -2, 1 },
				{ -1, -1, 0, -2, 1 },
				{ -3, 0, 0, 0, 2 },
				{ -3, 0, -2, 2, 0 },
				{ 0, 1, 0, -4, 1 },
				{ -2, 1, 0, -2, 1 },
				{ -4, 0, 0, 0, 1 },
				{ -1, 0, 0, -4, 1 },
				{ -3, 0, 0, -2, 1 },
				{ 0, 0, 0, 3, 2 },
				{ -1, 1, 0, 4, 1 },
				{ 1, -2, 2, 0, 1 },
				{ 0, 1, 0, 3, 0 },
				{ -1, 0, 2, 2, 3 },
				{ 0, 0, 2, 2, 2 },
				{ -2, 0, 2, 2, 2 },
				{ -1, 1, 2, 2, 0 },
				{ 3, 0, 0, 0, 2 },
				{ 2, 1, 0, 1, 0 },
				{ 2, -1, 2, -1, 2 },
				{ 0, 0, 2, 0, 1 },
				{ 0, 0, 3, 0, 3 },
				{ 0, 0, 3, 0, 2 },
				{ -1, 2, 2, 2, 1 },
				{ -1, 0, 4, 0, 0 },
				{ 1, 2, 2, 0, 1 },
				{ 3, 1, 2, -2, 1 },
				{ 1, 1, 4, -2, 2 },
				{ -2, -1, 0, 6, 0 },
				{ 0, -2, 0, 4, 0 },
				{ -2, 0, 0, 6, 1 },
				{ -2, -2, 2, 4, 2 },
				{ 0, -3, 2, 2, 2 },
				{ 0, 0, 0, 4, 2 },
				{ -1, -1, 2, 3, 2 },
				{ -2, 0, 2, 4, 0 },
				{ 2, -1, 0, 2, 1 },
				{ 1, 0, 0, 3, 0 },
				{ 0, 1, 0, 4, 1 },
				{ 0, 1, 0, 4, 0 },
				{ 1, -1, 2, 1, 2 },
				{ 0, 0, 2, 2, 3 },
				{ 1, 0, 2, 2, 2 },
				{ -1, 0, 2, 2, 2 },
				{ -2, 0, 4, 2, 1 },
				{ 2, 1, 0, 2, 1 },
				{ 2, 1, 0, 2, 0 },
				{ 2, -1, 2, 0, 0 },
				{ 1, 0, 2, 1, 0 },
				{ 0, 1, 2, 2, 0 },
				{ 2, 0, 2, 0, 3 },
				{ 3, 0, 2, 0, 2 },
				{ 1, 0, 2, 0, 2 },
				{ 1, 0, 3, 0, 3 },
				{ 1, 1, 2, 1, 1 },
				{ 0, 2, 2, 2, 2 },
				{ 2, 1, 2, 0, 0 },
				{ 2, 0, 4, -2, 1 },
				{ 4, 1, 2, -2, 2 },
				{ -1, -1, 0, 6, 0 },
				{ -3, -1, 2, 6, 2 },
				{ -1, 0, 0, 6, 1 },
				{ -3, 0, 2, 6, 1 },
				{ 1, -1, 0, 4, 1 },
				{ 1, -1, 0, 4, 0 },
				{ -2, 0, 2, 5, 2 },
				{ 1, -2, 2, 2, 1 },
				{ 3, -1, 0, 2, 0 },
				{ 1, -1, 2, 2, 0 },
				{ 0, 0, 2, 3, 1 },
				{ -1, 1, 2, 4, 1 },
				{ 0, 1, 2, 3, 2 },
				{ -1, 0, 4, 2, 1 },
				{ 2, 0, 2, 1, 1 },
				{ 5, 0, 0, 0, 0 },
				{ 2, 1, 2, 1, 2 },
				{ 1, 0, 4, 0, 1 },
				{ 3, 1, 2, 0, 1 },
				{ 3, 0, 4, -2, 2 },
				{ -2, -1, 2, 6, 2 },
				{ 0, 0, 0, 6, 0 },
				{ 0, -2, 2, 4, 2 },
				{ -2, 0, 2, 6, 1 },
				{ 2, 0, 0, 4, 1 },
				{ 2, 0, 0, 4, 0 },
				{ 2, -2, 2, 2, 2 },
				{ 0, 0, 2, 4, 0 },
				{ 1, 0, 2, 3, 2 },
				{ 4, 0, 0, 2, 0 },
				{ 2, 0, 2, 2, 0 },
				{ 0, 0, 4, 2, 2 },
				{ 4, -1, 2, 0, 2 },
				{ 3, 0, 2, 1, 2 },
				{ 2, 1, 2, 2, 1 },
				{ 4, 1, 2, 0, 2 },
				{ -1, -1, 2, 6, 2 },
				{ -1, 0, 2, 6, 1 },
				{ 1, -1, 2, 4, 1 },
				{ 1, 1, 2, 4, 2 },
				{ 3, 1, 2, 2, 2 },
				{ 5, 0, 2, 0, 1 },
				{ 2, -1, 2, 4, 2 },
				{ 2, 0, 2, 4, 1 }
		};
		for (I = 0; I < iNLS; I++)
			for (J = 0; J < 5; J++)
			{
				iNALS[J][I] = iNALST[I][J];
			}
		//CLST is the transposition of CLS
		double dCLST[iNLS][6] = {
				{ -172064161e0, -174666e0, 33386e0, 92052331e0, 9086e0, 15377e0 },
				{ -13170906e0, -1675e0, -13696e0, 5730336e0, -3015e0, -4587e0 },
				{ -2276413e0, -234e0, 2796e0, 978459e0, -485e0, 1374e0 },
				{ 2074554e0, 207e0, -698e0, -897492e0, 470e0, -291e0 },
				{ 1475877e0, -3633e0, 11817e0, 73871e0, -184e0, -1924e0 },
				{ -516821e0, 1226e0, -524e0, 224386e0, -677e0, -174e0 },
				{ 711159e0, 73e0, -872e0, -6750e0, 0e0, 358e0 },
				{ -387298e0, -367e0, 380e0, 200728e0, 18e0, 318e0 },
				{ -301461e0, -36e0, 816e0, 129025e0, -63e0, 367e0 },
				{ 215829e0, -494e0, 111e0, -95929e0, 299e0, 132e0 },
				{ 128227e0, 137e0, 181e0, -68982e0, -9e0, 39e0 },
				{ 123457e0, 11e0, 19e0, -53311e0, 32e0, -4e0 },
				{ 156994e0, 10e0, -168e0, -1235e0, 0e0, 82e0 },
				{ 63110e0, 63e0, 27e0, -33228e0, 0e0, -9e0 },
				{ -57976e0, -63e0, -189e0, 31429e0, 0e0, -75e0 },
				{ -59641e0, -11e0, 149e0, 25543e0, -11e0, 66e0 },
				{ -51613e0, -42e0, 129e0, 26366e0, 0e0, 78e0 },
				{ 45893e0, 50e0, 31e0, -24236e0, -10e0, 20e0 },
				{ 63384e0, 11e0, -150e0, -1220e0, 0e0, 29e0 },
				{ -38571e0, -1e0, 158e0, 16452e0, -11e0, 68e0 },
				{ 32481e0, 0e0, 0e0, -13870e0, 0e0, 0e0 },
				{ -47722e0, 0e0, -18e0, 477e0, 0e0, -25e0 },
				{ -31046e0, -1e0, 131e0, 13238e0, -11e0, 59e0 },
				{ 28593e0, 0e0, -1e0, -12338e0, 10e0, -3e0 },
				{ 20441e0, 21e0, 10e0, -10758e0, 0e0, -3e0 },
				{ 29243e0, 0e0, -74e0, -609e0, 0e0, 13e0 },
				{ 25887e0, 0e0, -66e0, -550e0, 0e0, 11e0 },
				{ -14053e0, -25e0, 79e0, 8551e0, -2e0, -45e0 },
				{ 15164e0, 10e0, 11e0, -8001e0, 0e0, -1e0 },
				{ -15794e0, 72e0, -16e0, 6850e0, -42e0, -5e0 },
				{ 21783e0, 0e0, 13e0, -167e0, 0e0, 13e0 },
				{ -12873e0, -10e0, -37e0, 6953e0, 0e0, -14e0 },
				{ -12654e0, 11e0, 63e0, 6415e0, 0e0, 26e0 },
				{ -10204e0, 0e0, 25e0, 5222e0, 0e0, 15e0 },
				{ 16707e0, -85e0, -10e0, 168e0, -1e0, 10e0 },
				{ -7691e0, 0e0, 44e0, 3268e0, 0e0, 19e0 },
				{ -11024e0, 0e0, -14e0, 104e0, 0e0, 2e0 },
				{ 7566e0, -21e0, -11e0, -3250e0, 0e0, -5e0 },
				{ -6637e0, -11e0, 25e0, 3353e0, 0e0, 14e0 },
				{ -7141e0, 21e0, 8e0, 3070e0, 0e0, 4e0 },
				{ -6302e0, -11e0, 2e0, 3272e0, 0e0, 4e0 },
				{ 5800e0, 10e0, 2e0, -3045e0, 0e0, -1e0 },
				{ 6443e0, 0e0, -7e0, -2768e0, 0e0, -4e0 },
				{ -5774e0, -11e0, -15e0, 3041e0, 0e0, -5e0 },
				{ -5350e0, 0e0, 21e0, 2695e0, 0e0, 12e0 },
				{ -4752e0, -11e0, -3e0, 2719e0, 0e0, -3e0 },
				{ -4940e0, -11e0, -21e0, 2720e0, 0e0, -9e0 },
				{ 7350e0, 0e0, -8e0, -51e0, 0e0, 4e0 },
				{ 4065e0, 0e0, 6e0, -2206e0, 0e0, 1e0 },
				{ 6579e0, 0e0, -24e0, -199e0, 0e0, 2e0 },
				{ 3579e0, 0e0, 5e0, -1900e0, 0e0, 1e0 },
				{ 4725e0, 0e0, -6e0, -41e0, 0e0, 3e0 },
				{ -3075e0, 0e0, -2e0, 1313e0, 0e0, -1e0 },
				{ -2904e0, 0e0, 15e0, 1233e0, 0e0, 7e0 },
				{ 4348e0, 0e0, -10e0, -81e0, 0e0, 2e0 },
				{ -2878e0, 0e0, 8e0, 1232e0, 0e0, 4e0 },
				{ -4230e0, 0e0, 5e0, -20e0, 0e0, -2e0 },
				{ -2819e0, 0e0, 7e0, 1207e0, 0e0, 3e0 },
				{ -4056e0, 0e0, 5e0, 40e0, 0e0, -2e0 },
				{ -2647e0, 0e0, 11e0, 1129e0, 0e0, 5e0 },
				{ -2294e0, 0e0, -10e0, 1266e0, 0e0, -4e0 },
				{ 2481e0, 0e0, -7e0, -1062e0, 0e0, -3e0 },
				{ 2179e0, 0e0, -2e0, -1129e0, 0e0, -2e0 },
				{ 3276e0, 0e0, 1e0, -9e0, 0e0, 0e0 },
				{ -3389e0, 0e0, 5e0, 35e0, 0e0, -2e0 },
				{ 3339e0, 0e0, -13e0, -107e0, 0e0, 1e0 },
				{ -1987e0, 0e0, -6e0, 1073e0, 0e0, -2e0 },
				{ -1981e0, 0e0, 0e0, 854e0, 0e0, 0e0 },
				{ 4026e0, 0e0, -353e0, -553e0, 0e0, -139e0 },
				{ 1660e0, 0e0, -5e0, -710e0, 0e0, -2e0 },
				{ -1521e0, 0e0, 9e0, 647e0, 0e0, 4e0 },
				{ 1314e0, 0e0, 0e0, -700e0, 0e0, 0e0 },
				{ -1283e0, 0e0, 0e0, 672e0, 0e0, 0e0 },
				{ -1331e0, 0e0, 8e0, 663e0, 0e0, 4e0 },
				{ 1383e0, 0e0, -2e0, -594e0, 0e0, -2e0 },
				{ 1405e0, 0e0, 4e0, -610e0, 0e0, 2e0 },
				{ 1290e0, 0e0, 0e0, -556e0, 0e0, 0e0 },
				{ -1214e0, 0e0, 5e0, 518e0, 0e0, 2e0 },
				{ 1146e0, 0e0, -3e0, -490e0, 0e0, -1e0 },
				{ 1019e0, 0e0, -1e0, -527e0, 0e0, -1e0 },
				{ -1100e0, 0e0, 9e0, 465e0, 0e0, 4e0 },
				{ -970e0, 0e0, 2e0, 496e0, 0e0, 1e0 },
				{ 1575e0, 0e0, -6e0, -50e0, 0e0, 0e0 },
				{ 934e0, 0e0, -3e0, -399e0, 0e0, -1e0 },
				{ 922e0, 0e0, -1e0, -395e0, 0e0, -1e0 },
				{ 815e0, 0e0, -1e0, -422e0, 0e0, -1e0 },
				{ 834e0, 0e0, 2e0, -440e0, 0e0, 1e0 },
				{ 1248e0, 0e0, 0e0, -170e0, 0e0, 1e0 },
				{ 1338e0, 0e0, -5e0, -39e0, 0e0, 0e0 },
				{ 716e0, 0e0, -2e0, -389e0, 0e0, -1e0 },
				{ 1282e0, 0e0, -3e0, -23e0, 0e0, 1e0 },
				{ 742e0, 0e0, 1e0, -391e0, 0e0, 0e0 },
				{ 1020e0, 0e0, -25e0, -495e0, 0e0, -10e0 },
				{ 715e0, 0e0, -4e0, -326e0, 0e0, 2e0 },
				{ -666e0, 0e0, -3e0, 369e0, 0e0, -1e0 },
				{ -667e0, 0e0, 1e0, 346e0, 0e0, 1e0 },
				{ -704e0, 0e0, 0e0, 304e0, 0e0, 0e0 },
				{ -694e0, 0e0, 5e0, 294e0, 0e0, 2e0 },
				{ -1014e0, 0e0, -1e0, 4e0, 0e0, -1e0 },
				{ -585e0, 0e0, -2e0, 316e0, 0e0, -1e0 },
				{ -949e0, 0e0, 1e0, 8e0, 0e0, -1e0 },
				{ -595e0, 0e0, 0e0, 258e0, 0e0, 0e0 },
				{ 528e0, 0e0, 0e0, -279e0, 0e0, 0e0 },
				{ -590e0, 0e0, 4e0, 252e0, 0e0, 2e0 },
				{ 570e0, 0e0, -2e0, -244e0, 0e0, -1e0 },
				{ -502e0, 0e0, 3e0, 250e0, 0e0, 2e0 },
				{ -875e0, 0e0, 1e0, 29e0, 0e0, 0e0 },
				{ -492e0, 0e0, -3e0, 275e0, 0e0, -1e0 },
				{ 535e0, 0e0, -2e0, -228e0, 0e0, -1e0 },
				{ -467e0, 0e0, 1e0, 240e0, 0e0, 1e0 },
				{ 591e0, 0e0, 0e0, -253e0, 0e0, 0e0 },
				{ -453e0, 0e0, -1e0, 244e0, 0e0, -1e0 },
				{ 766e0, 0e0, 1e0, 9e0, 0e0, 0e0 },
				{ -446e0, 0e0, 2e0, 225e0, 0e0, 1e0 },
				{ -488e0, 0e0, 2e0, 207e0, 0e0, 1e0 },
				{ -468e0, 0e0, 0e0, 201e0, 0e0, 0e0 },
				{ -421e0, 0e0, 1e0, 216e0, 0e0, 1e0 },
				{ 463e0, 0e0, 0e0, -200e0, 0e0, 0e0 },
				{ -673e0, 0e0, 2e0, 14e0, 0e0, 0e0 },
				{ 658e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -438e0, 0e0, 0e0, 188e0, 0e0, 0e0 },
				{ -390e0, 0e0, 0e0, 205e0, 0e0, 0e0 },
				{ 639e0, -11e0, -2e0, -19e0, 0e0, 0e0 },
				{ 412e0, 0e0, -2e0, -176e0, 0e0, -1e0 },
				{ -361e0, 0e0, 0e0, 189e0, 0e0, 0e0 },
				{ 360e0, 0e0, -1e0, -185e0, 0e0, -1e0 },
				{ 588e0, 0e0, -3e0, -24e0, 0e0, 0e0 },
				{ -578e0, 0e0, 1e0, 5e0, 0e0, 0e0 },
				{ -396e0, 0e0, 0e0, 171e0, 0e0, 0e0 },
				{ 565e0, 0e0, -1e0, -6e0, 0e0, 0e0 },
				{ -335e0, 0e0, -1e0, 184e0, 0e0, -1e0 },
				{ 357e0, 0e0, 1e0, -154e0, 0e0, 0e0 },
				{ 321e0, 0e0, 1e0, -174e0, 0e0, 0e0 },
				{ -301e0, 0e0, -1e0, 162e0, 0e0, 0e0 },
				{ -334e0, 0e0, 0e0, 144e0, 0e0, 0e0 },
				{ 493e0, 0e0, -2e0, -15e0, 0e0, 0e0 },
				{ 494e0, 0e0, -2e0, -19e0, 0e0, 0e0 },
				{ 337e0, 0e0, -1e0, -143e0, 0e0, -1e0 },
				{ 280e0, 0e0, -1e0, -144e0, 0e0, 0e0 },
				{ 309e0, 0e0, 1e0, -134e0, 0e0, 0e0 },
				{ -263e0, 0e0, 2e0, 131e0, 0e0, 1e0 },
				{ 253e0, 0e0, 1e0, -138e0, 0e0, 0e0 },
				{ 245e0, 0e0, 0e0, -128e0, 0e0, 0e0 },
				{ 416e0, 0e0, -2e0, -17e0, 0e0, 0e0 },
				{ -229e0, 0e0, 0e0, 128e0, 0e0, 0e0 },
				{ 231e0, 0e0, 0e0, -120e0, 0e0, 0e0 },
				{ -259e0, 0e0, 2e0, 109e0, 0e0, 1e0 },
				{ 375e0, 0e0, -1e0, -8e0, 0e0, 0e0 },
				{ 252e0, 0e0, 0e0, -108e0, 0e0, 0e0 },
				{ -245e0, 0e0, 1e0, 104e0, 0e0, 0e0 },
				{ 243e0, 0e0, -1e0, -104e0, 0e0, 0e0 },
				{ 208e0, 0e0, 1e0, -112e0, 0e0, 0e0 },
				{ 199e0, 0e0, 0e0, -102e0, 0e0, 0e0 },
				{ -208e0, 0e0, 1e0, 105e0, 0e0, 0e0 },
				{ 335e0, 0e0, -2e0, -14e0, 0e0, 0e0 },
				{ -325e0, 0e0, 1e0, 7e0, 0e0, 0e0 },
				{ -187e0, 0e0, 0e0, 96e0, 0e0, 0e0 },
				{ 197e0, 0e0, -1e0, -100e0, 0e0, 0e0 },
				{ -192e0, 0e0, 2e0, 94e0, 0e0, 1e0 },
				{ -188e0, 0e0, 0e0, 83e0, 0e0, 0e0 },
				{ 276e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -286e0, 0e0, 1e0, 6e0, 0e0, 0e0 },
				{ 186e0, 0e0, -1e0, -79e0, 0e0, 0e0 },
				{ -219e0, 0e0, 0e0, 43e0, 0e0, 0e0 },
				{ 276e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -153e0, 0e0, -1e0, 84e0, 0e0, 0e0 },
				{ -156e0, 0e0, 0e0, 81e0, 0e0, 0e0 },
				{ -154e0, 0e0, 1e0, 78e0, 0e0, 0e0 },
				{ -174e0, 0e0, 1e0, 75e0, 0e0, 0e0 },
				{ -163e0, 0e0, 2e0, 69e0, 0e0, 1e0 },
				{ -228e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ 91e0, 0e0, -4e0, -54e0, 0e0, -2e0 },
				{ 175e0, 0e0, 0e0, -75e0, 0e0, 0e0 },
				{ -159e0, 0e0, 0e0, 69e0, 0e0, 0e0 },
				{ 141e0, 0e0, 0e0, -72e0, 0e0, 0e0 },
				{ 147e0, 0e0, 0e0, -75e0, 0e0, 0e0 },
				{ -132e0, 0e0, 0e0, 69e0, 0e0, 0e0 },
				{ 159e0, 0e0, -28e0, -54e0, 0e0, 11e0 },
				{ 213e0, 0e0, 0e0, -4e0, 0e0, 0e0 },
				{ 123e0, 0e0, 0e0, -64e0, 0e0, 0e0 },
				{ -118e0, 0e0, -1e0, 66e0, 0e0, 0e0 },
				{ 144e0, 0e0, -1e0, -61e0, 0e0, 0e0 },
				{ -121e0, 0e0, 1e0, 60e0, 0e0, 0e0 },
				{ -134e0, 0e0, 1e0, 56e0, 0e0, 1e0 },
				{ -105e0, 0e0, 0e0, 57e0, 0e0, 0e0 },
				{ -102e0, 0e0, 0e0, 56e0, 0e0, 0e0 },
				{ 120e0, 0e0, 0e0, -52e0, 0e0, 0e0 },
				{ 101e0, 0e0, 0e0, -54e0, 0e0, 0e0 },
				{ -113e0, 0e0, 0e0, 59e0, 0e0, 0e0 },
				{ -106e0, 0e0, 0e0, 61e0, 0e0, 0e0 },
				{ -129e0, 0e0, 1e0, 55e0, 0e0, 0e0 },
				{ -114e0, 0e0, 0e0, 57e0, 0e0, 0e0 },
				{ 113e0, 0e0, -1e0, -49e0, 0e0, 0e0 },
				{ -102e0, 0e0, 0e0, 44e0, 0e0, 0e0 },
				{ -94e0, 0e0, 0e0, 51e0, 0e0, 0e0 },
				{ -100e0, 0e0, -1e0, 56e0, 0e0, 0e0 },
				{ 87e0, 0e0, 0e0, -47e0, 0e0, 0e0 },
				{ 161e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ 96e0, 0e0, 0e0, -50e0, 0e0, 0e0 },
				{ 151e0, 0e0, -1e0, -5e0, 0e0, 0e0 },
				{ -104e0, 0e0, 0e0, 44e0, 0e0, 0e0 },
				{ -110e0, 0e0, 0e0, 48e0, 0e0, 0e0 },
				{ -100e0, 0e0, 1e0, 50e0, 0e0, 0e0 },
				{ 92e0, 0e0, -5e0, 12e0, 0e0, -2e0 },
				{ 82e0, 0e0, 0e0, -45e0, 0e0, 0e0 },
				{ 82e0, 0e0, 0e0, -45e0, 0e0, 0e0 },
				{ -78e0, 0e0, 0e0, 41e0, 0e0, 0e0 },
				{ -77e0, 0e0, 0e0, 43e0, 0e0, 0e0 },
				{ 2e0, 0e0, 0e0, 54e0, 0e0, 0e0 },
				{ 94e0, 0e0, 0e0, -40e0, 0e0, 0e0 },
				{ -93e0, 0e0, 0e0, 40e0, 0e0, 0e0 },
				{ -83e0, 0e0, 10e0, 40e0, 0e0, -2e0 },
				{ 83e0, 0e0, 0e0, -36e0, 0e0, 0e0 },
				{ -91e0, 0e0, 0e0, 39e0, 0e0, 0e0 },
				{ 128e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ -79e0, 0e0, 0e0, 34e0, 0e0, 0e0 },
				{ -83e0, 0e0, 0e0, 47e0, 0e0, 0e0 },
				{ 84e0, 0e0, 0e0, -44e0, 0e0, 0e0 },
				{ 83e0, 0e0, 0e0, -43e0, 0e0, 0e0 },
				{ 91e0, 0e0, 0e0, -39e0, 0e0, 0e0 },
				{ -77e0, 0e0, 0e0, 39e0, 0e0, 0e0 },
				{ 84e0, 0e0, 0e0, -43e0, 0e0, 0e0 },
				{ -92e0, 0e0, 1e0, 39e0, 0e0, 0e0 },
				{ -92e0, 0e0, 1e0, 39e0, 0e0, 0e0 },
				{ -94e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 68e0, 0e0, 0e0, -36e0, 0e0, 0e0 },
				{ -61e0, 0e0, 0e0, 32e0, 0e0, 0e0 },
				{ 71e0, 0e0, 0e0, -31e0, 0e0, 0e0 },
				{ 62e0, 0e0, 0e0, -34e0, 0e0, 0e0 },
				{ -63e0, 0e0, 0e0, 33e0, 0e0, 0e0 },
				{ -73e0, 0e0, 0e0, 32e0, 0e0, 0e0 },
				{ 115e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -103e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 63e0, 0e0, 0e0, -28e0, 0e0, 0e0 },
				{ 74e0, 0e0, 0e0, -32e0, 0e0, 0e0 },
				{ -103e0, 0e0, -3e0, 3e0, 0e0, -1e0 },
				{ -69e0, 0e0, 0e0, 30e0, 0e0, 0e0 },
				{ 57e0, 0e0, 0e0, -29e0, 0e0, 0e0 },
				{ 94e0, 0e0, 0e0, -4e0, 0e0, 0e0 },
				{ 64e0, 0e0, 0e0, -33e0, 0e0, 0e0 },
				{ -63e0, 0e0, 0e0, 26e0, 0e0, 0e0 },
				{ -38e0, 0e0, 0e0, 20e0, 0e0, 0e0 },
				{ -43e0, 0e0, 0e0, 24e0, 0e0, 0e0 },
				{ -45e0, 0e0, 0e0, 23e0, 0e0, 0e0 },
				{ 47e0, 0e0, 0e0, -24e0, 0e0, 0e0 },
				{ -48e0, 0e0, 0e0, 25e0, 0e0, 0e0 },
				{ 45e0, 0e0, 0e0, -26e0, 0e0, 0e0 },
				{ 56e0, 0e0, 0e0, -25e0, 0e0, 0e0 },
				{ 88e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -75e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 85e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 49e0, 0e0, 0e0, -26e0, 0e0, 0e0 },
				{ -74e0, 0e0, -3e0, -1e0, 0e0, -1e0 },
				{ -39e0, 0e0, 0e0, 21e0, 0e0, 0e0 },
				{ 45e0, 0e0, 0e0, -20e0, 0e0, 0e0 },
				{ 51e0, 0e0, 0e0, -22e0, 0e0, 0e0 },
				{ -40e0, 0e0, 0e0, 21e0, 0e0, 0e0 },
				{ 41e0, 0e0, 0e0, -21e0, 0e0, 0e0 },
				{ -42e0, 0e0, 0e0, 24e0, 0e0, 0e0 },
				{ -51e0, 0e0, 0e0, 22e0, 0e0, 0e0 },
				{ -42e0, 0e0, 0e0, 22e0, 0e0, 0e0 },
				{ 39e0, 0e0, 0e0, -21e0, 0e0, 0e0 },
				{ 46e0, 0e0, 0e0, -18e0, 0e0, 0e0 },
				{ -53e0, 0e0, 0e0, 22e0, 0e0, 0e0 },
				{ 82e0, 0e0, 0e0, -4e0, 0e0, 0e0 },
				{ 81e0, 0e0, -1e0, -4e0, 0e0, 0e0 },
				{ 47e0, 0e0, 0e0, -19e0, 0e0, 0e0 },
				{ 53e0, 0e0, 0e0, -23e0, 0e0, 0e0 },
				{ -45e0, 0e0, 0e0, 22e0, 0e0, 0e0 },
				{ -44e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -33e0, 0e0, 0e0, 16e0, 0e0, 0e0 },
				{ -61e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ 28e0, 0e0, 0e0, -15e0, 0e0, 0e0 },
				{ -38e0, 0e0, 0e0, 19e0, 0e0, 0e0 },
				{ -33e0, 0e0, 0e0, 21e0, 0e0, 0e0 },
				{ -60e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 48e0, 0e0, 0e0, -10e0, 0e0, 0e0 },
				{ 27e0, 0e0, 0e0, -14e0, 0e0, 0e0 },
				{ 38e0, 0e0, 0e0, -20e0, 0e0, 0e0 },
				{ 31e0, 0e0, 0e0, -13e0, 0e0, 0e0 },
				{ -29e0, 0e0, 0e0, 15e0, 0e0, 0e0 },
				{ 28e0, 0e0, 0e0, -15e0, 0e0, 0e0 },
				{ -32e0, 0e0, 0e0, 15e0, 0e0, 0e0 },
				{ 45e0, 0e0, 0e0, -8e0, 0e0, 0e0 },
				{ -44e0, 0e0, 0e0, 19e0, 0e0, 0e0 },
				{ 28e0, 0e0, 0e0, -15e0, 0e0, 0e0 },
				{ -51e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -36e0, 0e0, 0e0, 20e0, 0e0, 0e0 },
				{ 44e0, 0e0, 0e0, -19e0, 0e0, 0e0 },
				{ 26e0, 0e0, 0e0, -14e0, 0e0, 0e0 },
				{ -60e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 35e0, 0e0, 0e0, -18e0, 0e0, 0e0 },
				{ -27e0, 0e0, 0e0, 11e0, 0e0, 0e0 },
				{ 47e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ 36e0, 0e0, 0e0, -15e0, 0e0, 0e0 },
				{ -36e0, 0e0, 0e0, 20e0, 0e0, 0e0 },
				{ -35e0, 0e0, 0e0, 19e0, 0e0, 0e0 },
				{ -37e0, 0e0, 0e0, 19e0, 0e0, 0e0 },
				{ 32e0, 0e0, 0e0, -16e0, 0e0, 0e0 },
				{ 35e0, 0e0, 0e0, -14e0, 0e0, 0e0 },
				{ 32e0, 0e0, 0e0, -13e0, 0e0, 0e0 },
				{ 65e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 47e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ 32e0, 0e0, 0e0, -16e0, 0e0, 0e0 },
				{ 37e0, 0e0, 0e0, -16e0, 0e0, 0e0 },
				{ -30e0, 0e0, 0e0, 15e0, 0e0, 0e0 },
				{ -32e0, 0e0, 0e0, 16e0, 0e0, 0e0 },
				{ -31e0, 0e0, 0e0, 13e0, 0e0, 0e0 },
				{ 37e0, 0e0, 0e0, -16e0, 0e0, 0e0 },
				{ 31e0, 0e0, 0e0, -13e0, 0e0, 0e0 },
				{ 49e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 32e0, 0e0, 0e0, -13e0, 0e0, 0e0 },
				{ 23e0, 0e0, 0e0, -12e0, 0e0, 0e0 },
				{ -43e0, 0e0, 0e0, 18e0, 0e0, 0e0 },
				{ 26e0, 0e0, 0e0, -11e0, 0e0, 0e0 },
				{ -32e0, 0e0, 0e0, 14e0, 0e0, 0e0 },
				{ -29e0, 0e0, 0e0, 14e0, 0e0, 0e0 },
				{ -27e0, 0e0, 0e0, 12e0, 0e0, 0e0 },
				{ 30e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -11e0, 0e0, 0e0, 5e0, 0e0, 0e0 },
				{ -21e0, 0e0, 0e0, 10e0, 0e0, 0e0 },
				{ -34e0, 0e0, 0e0, 15e0, 0e0, 0e0 },
				{ -10e0, 0e0, 0e0, 6e0, 0e0, 0e0 },
				{ -36e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -9e0, 0e0, 0e0, 4e0, 0e0, 0e0 },
				{ -12e0, 0e0, 0e0, 5e0, 0e0, 0e0 },
				{ -21e0, 0e0, 0e0, 5e0, 0e0, 0e0 },
				{ -29e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ -15e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ -20e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 28e0, 0e0, 0e0, 0e0, 0e0, -2e0 },
				{ 17e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -22e0, 0e0, 0e0, 12e0, 0e0, 0e0 },
				{ -14e0, 0e0, 0e0, 7e0, 0e0, 0e0 },
				{ 24e0, 0e0, 0e0, -11e0, 0e0, 0e0 },
				{ 11e0, 0e0, 0e0, -6e0, 0e0, 0e0 },
				{ 14e0, 0e0, 0e0, -6e0, 0e0, 0e0 },
				{ 24e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 18e0, 0e0, 0e0, -8e0, 0e0, 0e0 },
				{ -38e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -31e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -16e0, 0e0, 0e0, 8e0, 0e0, 0e0 },
				{ 29e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -18e0, 0e0, 0e0, 10e0, 0e0, 0e0 },
				{ -10e0, 0e0, 0e0, 5e0, 0e0, 0e0 },
				{ -17e0, 0e0, 0e0, 10e0, 0e0, 0e0 },
				{ 9e0, 0e0, 0e0, -4e0, 0e0, 0e0 },
				{ 16e0, 0e0, 0e0, -6e0, 0e0, 0e0 },
				{ 22e0, 0e0, 0e0, -12e0, 0e0, 0e0 },
				{ 20e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -13e0, 0e0, 0e0, 6e0, 0e0, 0e0 },
				{ -17e0, 0e0, 0e0, 9e0, 0e0, 0e0 },
				{ -14e0, 0e0, 0e0, 8e0, 0e0, 0e0 },
				{ 0e0, 0e0, 0e0, -7e0, 0e0, 0e0 },
				{ 14e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 19e0, 0e0, 0e0, -10e0, 0e0, 0e0 },
				{ -34e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -20e0, 0e0, 0e0, 8e0, 0e0, 0e0 },
				{ 9e0, 0e0, 0e0, -5e0, 0e0, 0e0 },
				{ -18e0, 0e0, 0e0, 7e0, 0e0, 0e0 },
				{ 13e0, 0e0, 0e0, -6e0, 0e0, 0e0 },
				{ 17e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -12e0, 0e0, 0e0, 5e0, 0e0, 0e0 },
				{ 15e0, 0e0, 0e0, -8e0, 0e0, 0e0 },
				{ -11e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ 13e0, 0e0, 0e0, -5e0, 0e0, 0e0 },
				{ -18e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -35e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 9e0, 0e0, 0e0, -4e0, 0e0, 0e0 },
				{ -19e0, 0e0, 0e0, 10e0, 0e0, 0e0 },
				{ -26e0, 0e0, 0e0, 11e0, 0e0, 0e0 },
				{ 8e0, 0e0, 0e0, -4e0, 0e0, 0e0 },
				{ -10e0, 0e0, 0e0, 4e0, 0e0, 0e0 },
				{ 10e0, 0e0, 0e0, -6e0, 0e0, 0e0 },
				{ -21e0, 0e0, 0e0, 9e0, 0e0, 0e0 },
				{ -15e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 9e0, 0e0, 0e0, -5e0, 0e0, 0e0 },
				{ -29e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -19e0, 0e0, 0e0, 10e0, 0e0, 0e0 },
				{ 12e0, 0e0, 0e0, -5e0, 0e0, 0e0 },
				{ 22e0, 0e0, 0e0, -9e0, 0e0, 0e0 },
				{ -10e0, 0e0, 0e0, 5e0, 0e0, 0e0 },
				{ -20e0, 0e0, 0e0, 11e0, 0e0, 0e0 },
				{ -20e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -17e0, 0e0, 0e0, 7e0, 0e0, 0e0 },
				{ 15e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ 8e0, 0e0, 0e0, -4e0, 0e0, 0e0 },
				{ 14e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -12e0, 0e0, 0e0, 6e0, 0e0, 0e0 },
				{ 25e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -13e0, 0e0, 0e0, 6e0, 0e0, 0e0 },
				{ -14e0, 0e0, 0e0, 8e0, 0e0, 0e0 },
				{ 13e0, 0e0, 0e0, -5e0, 0e0, 0e0 },
				{ -17e0, 0e0, 0e0, 9e0, 0e0, 0e0 },
				{ -12e0, 0e0, 0e0, 6e0, 0e0, 0e0 },
				{ -10e0, 0e0, 0e0, 5e0, 0e0, 0e0 },
				{ 10e0, 0e0, 0e0, -6e0, 0e0, 0e0 },
				{ -15e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -22e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 28e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ 15e0, 0e0, 0e0, -7e0, 0e0, 0e0 },
				{ 23e0, 0e0, 0e0, -10e0, 0e0, 0e0 },
				{ 12e0, 0e0, 0e0, -5e0, 0e0, 0e0 },
				{ 29e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ -25e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ 22e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -18e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 15e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ -23e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 12e0, 0e0, 0e0, -5e0, 0e0, 0e0 },
				{ -8e0, 0e0, 0e0, 4e0, 0e0, 0e0 },
				{ -19e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -10e0, 0e0, 0e0, 4e0, 0e0, 0e0 },
				{ 21e0, 0e0, 0e0, -9e0, 0e0, 0e0 },
				{ 23e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ -16e0, 0e0, 0e0, 8e0, 0e0, 0e0 },
				{ -19e0, 0e0, 0e0, 9e0, 0e0, 0e0 },
				{ -22e0, 0e0, 0e0, 10e0, 0e0, 0e0 },
				{ 27e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ 16e0, 0e0, 0e0, -8e0, 0e0, 0e0 },
				{ 19e0, 0e0, 0e0, -8e0, 0e0, 0e0 },
				{ 9e0, 0e0, 0e0, -4e0, 0e0, 0e0 },
				{ -9e0, 0e0, 0e0, 4e0, 0e0, 0e0 },
				{ -9e0, 0e0, 0e0, 4e0, 0e0, 0e0 },
				{ -8e0, 0e0, 0e0, 4e0, 0e0, 0e0 },
				{ 18e0, 0e0, 0e0, -9e0, 0e0, 0e0 },
				{ 16e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ -10e0, 0e0, 0e0, 4e0, 0e0, 0e0 },
				{ -23e0, 0e0, 0e0, 9e0, 0e0, 0e0 },
				{ 16e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ -12e0, 0e0, 0e0, 6e0, 0e0, 0e0 },
				{ -8e0, 0e0, 0e0, 4e0, 0e0, 0e0 },
				{ 30e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 24e0, 0e0, 0e0, -10e0, 0e0, 0e0 },
				{ 10e0, 0e0, 0e0, -4e0, 0e0, 0e0 },
				{ -16e0, 0e0, 0e0, 7e0, 0e0, 0e0 },
				{ -16e0, 0e0, 0e0, 7e0, 0e0, 0e0 },
				{ 17e0, 0e0, 0e0, -7e0, 0e0, 0e0 },
				{ -24e0, 0e0, 0e0, 10e0, 0e0, 0e0 },
				{ -12e0, 0e0, 0e0, 5e0, 0e0, 0e0 },
				{ -24e0, 0e0, 0e0, 11e0, 0e0, 0e0 },
				{ -23e0, 0e0, 0e0, 9e0, 0e0, 0e0 },
				{ -13e0, 0e0, 0e0, 5e0, 0e0, 0e0 },
				{ -15e0, 0e0, 0e0, 7e0, 0e0, 0e0 },
				{ 0e0, 0e0, -1988e0, 0e0, 0e0, -1679e0 },
				{ 0e0, 0e0, -63e0, 0e0, 0e0, -27e0 },
				{ -4e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 0e0, 0e0, 5e0, 0e0, 0e0, 4e0 },
				{ 5e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ 0e0, 0e0, 364e0, 0e0, 0e0, 176e0 },
				{ 0e0, 0e0, -1044e0, 0e0, 0e0, -891e0 },
				{ -3e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 0e0, 0e0, 330e0, 0e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ -5e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 0e0, 0e0, 5e0, 0e0, 0e0, 0e0 },
				{ 0e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 6e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -7e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -12e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ -5e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -7e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ 7e0, 0e0, 0e0, -4e0, 0e0, 0e0 },
				{ 0e0, 0e0, -12e0, 0e0, 0e0, -10e0 },
				{ 4e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -7e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ 0e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ 7e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -5e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -5e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -8e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ 9e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 6e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ -5e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -7e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ -5e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 9e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 9e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 8e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ 6e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ -7e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 9e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -5e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ -13e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -7e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 10e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ 10e0, 0e0, 13e0, 6e0, 0e0, -5e0 },
				{ 0e0, 0e0, 30e0, 0e0, 0e0, 14e0 },
				{ 0e0, 0e0, -162e0, 0e0, 0e0, -138e0 },
				{ 0e0, 0e0, 75e0, 0e0, 0e0, 0e0 },
				{ -7e0, 0e0, 0e0, 4e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -5e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 6e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 9e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -7e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 7e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -6e0, 0e0, -3e0, 3e0, 0e0, 1e0 },
				{ 0e0, 0e0, -3e0, 0e0, 0e0, -2e0 },
				{ 11e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ 11e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -1e0, 0e0, 3e0, 3e0, 0e0, -1e0 },
				{ 4e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 0e0, 0e0, -13e0, 0e0, 0e0, -11e0 },
				{ 3e0, 0e0, 6e0, 0e0, 0e0, 0e0 },
				{ -7e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ -7e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ 8e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 11e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 8e0, 0e0, 0e0, -4e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ 11e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -6e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -8e0, 0e0, 0e0, 4e0, 0e0, 0e0 },
				{ -7e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ 6e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ -6e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ 6e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 6e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -5e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 6e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 0e0, 0e0, -26e0, 0e0, 0e0, -11e0 },
				{ 0e0, 0e0, -10e0, 0e0, 0e0, -5e0 },
				{ 5e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ -13e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 7e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -6e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -5e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -7e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 13e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -11e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 6e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -12e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ 0e0, 0e0, -5e0, 0e0, 0e0, -2e0 },
				{ -7e0, 0e0, 0e0, 4e0, 0e0, 0e0 },
				{ 6e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ -5e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 12e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 6e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -6e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 6e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ 6e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -6e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 7e0, 0e0, 0e0, -4e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -5e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -6e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ -6e0, 0e0, 0e0, 3e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 10e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 7e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 7e0, 0e0, 0e0, -3e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 11e0, 0e0, 0e0, 0e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -6e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 5e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -4e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 2e0, 0e0, 0e0 },
				{ 4e0, 0e0, 0e0, -2e0, 0e0, 0e0 },
				{ 3e0, 0e0, 0e0, -1e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 1e0, 0e0, 0e0 },
				{ -3e0, 0e0, 0e0, 2e0, 0e0, 0e0 }
		};
		for (I = 0; I < iNLS; I++)
			for (J = 0; J < 6; J++)
			{
				dCLS[J][I] = dCLST[I][J];
			}
		//NAPLT is the transposition of NAPL
		int iNAPLT[iNPL][14] = {
				{ 0, 0, 0, 0, 0, 0, 0, 8, -16, 4, 5, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -8, 16, -4, -5, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 8, -16, 4, 5, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 2, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -4, 8, -1, -5, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, -8, 3, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, 0, 3, -8, 3, 0, 0, 0, 0 },
				{ -1, 0, 0, 0, 0, 0, 10, -3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 6, -3, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, -8, 3, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 0, -5, 8, -3, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -4, 8, -3, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, -8, 1, 5, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -5, 6, 4, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -5, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -5, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, 2, -5, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -5, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, -2, 5, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 5, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 5, 0, 0, 2 },
				{ 2, 0, -1, -1, 0, 0, 0, 3, -7, 0, 0, 0, 0, 0 },
				{ 1, 0, 0, -2, 0, 0, 19, -21, 3, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 2, -4, 0, -3, 0, 0, 0, 0 },
				{ 1, 0, 0, -1, 1, 0, 0, -1, 0, 2, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, -4, 10, 0, 0, 0 },
				{ -2, 0, 0, 2, 1, 0, 0, 2, 0, 0, -5, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 3, -7, 4, 0, 0, 0, 0, 0 },
				{ 0, 0, -1, 1, 0, 0, 0, 1, 0, 1, -1, 0, 0, 0 },
				{ -2, 0, 0, 2, 1, 0, 0, 2, 0, -2, 0, 0, 0, 0 },
				{ -1, 0, 0, 0, 0, 0, 18, -16, 0, 0, 0, 0, 0, 0 },
				{ -2, 0, 1, 1, 2, 0, 0, 1, 0, -2, 0, 0, 0, 0 },
				{ -1, 0, 1, -1, 1, 0, 18, -17, 0, 0, 0, 0, 0, 0 },
				{ -1, 0, 0, 1, 1, 0, 0, 2, -2, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -8, 13, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 2, -2, 2, 0, -8, 11, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -8, 13, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, -8, 12, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 8, -13, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 8, -14, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 8, -13, 0, 0, 0, 0, 0, 1 },
				{ -2, 0, 0, 2, 1, 0, 0, 2, 0, -4, 5, 0, 0, 0 },
				{ -2, 0, 0, 2, 2, 0, 3, -3, 0, 0, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 0, 0, 0, 2, 0, -3, 1, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 3, -5, 0, 2, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 0, 0, 0, 2, 0, -4, 3, 0, 0, 0 },
				{ 0, 0, -1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, -1, 2, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 2, 0, 0, -2, 2, 0, 0, 0, 0, 0 },
				{ -1, 0, 1, 0, 1, 0, 3, -5, 0, 0, 0, 0, 0, 0 },
				{ -1, 0, 0, 1, 0, 0, 3, -4, 0, 0, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 0, 0, 0, 2, 0, -2, -2, 0, 0, 0 },
				{ -2, 0, 2, 0, 2, 0, 0, -5, 9, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, 0, 0, -1, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, 0, 0, 0, 2, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2 },
				{ -1, 0, 0, 1, 0, 0, 0, 3, -4, 0, 0, 0, 0, 0 },
				{ 0, 0, -1, 1, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0 },
				{ 0, 0, 1, -1, 2, 0, 0, -1, 0, 0, 2, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, -9, 17, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 2, 0, -3, 5, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, -1, 2, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0 },
				{ 1, 0, 0, -2, 0, 0, 17, -16, 0, -2, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, 1, -3, 0, 0, 0 },
				{ -2, 0, 0, 2, 1, 0, 0, 5, -6, 0, 0, 0, 0, 0 },
				{ 0, 0, -2, 2, 0, 0, 0, 9, -13, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 2, 0, 0, -1, 0, 0, 1, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0 },
				{ 0, 0, -1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 },
				{ 0, 0, -2, 2, 0, 0, 5, -6, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, -1, 1, 1, 0, 5, -7, 0, 0, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 0, 0, 6, -8, 0, 0, 0, 0, 0, 0 },
				{ 2, 0, 1, -3, 1, 0, -6, 7, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0 },
				{ 0, 0, -1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, 0, 0, 2, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -8, 15, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -8, 15, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, 0, -9, 15, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 8, -15, 0, 0, 0, 0, 0 },
				{ 1, 0, -1, -1, 0, 0, 0, 8, -15, 0, 0, 0, 0, 0 },
				{ 2, 0, 0, -2, 0, 0, 2, -5, 0, 0, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 0, 0, 0, 2, 0, -5, 5, 0, 0, 0 },
				{ 2, 0, 0, -2, 1, 0, 0, -6, 8, 0, 0, 0, 0, 0 },
				{ 2, 0, 0, -2, 1, 0, 0, -2, 0, 3, 0, 0, 0, 0 },
				{ -2, 0, 1, 1, 0, 0, 0, 1, 0, -3, 0, 0, 0, 0 },
				{ -2, 0, 1, 1, 1, 0, 0, 1, 0, -3, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 0, 0, 0, 2, 0, -3, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 0, 0, 0, 6, -8, 0, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 0, 0, 0, 2, 0, -1, -5, 0, 0, 0 },
				{ -1, 0, 0, 1, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0 },
				{ -1, 0, 1, 1, 1, 0, -20, 20, 0, 0, 0, 0, 0, 0 },
				{ 1, 0, 0, -2, 0, 0, 20, -21, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, 8, -15, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, 0, -10, 15, 0, 0, 0, 0, 0 },
				{ 0, 0, -1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 2, 0, 0, -1, 0, 1, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, -2, 4, 0, 0, 0 },
				{ 2, 0, 0, -2, 1, 0, -6, 8, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, -2, 2, 1, 0, 5, -6, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, 0, -1, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, 0, 1, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2 },
				{ 0, 0, 2, -2, 1, 0, 0, -9, 13, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, 7, -13, 0, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 0, 0, 0, 5, -6, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 9, -17, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -9, 17, 0, 0, 0, 0, 2 },
				{ 1, 0, 0, -1, 1, 0, 0, -3, 4, 0, 0, 0, 0, 0 },
				{ 1, 0, 0, -1, 1, 0, -3, 4, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 2, 0, 0, -1, 2, 0, 0, 0, 0, 0 },
				{ 0, 0, -1, 1, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0 },
				{ 0, 0, -2, 2, 0, 1, 0, -2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 3, -5, 0, 2, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 1, 0, 0, 2, 0, -3, 1, 0, 0, 0 },
				{ -2, 0, 0, 2, 1, 0, 3, -3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 8, -13, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, -1, 1, 0, 0, 8, -12, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, -8, 11, 0, 0, 0, 0, 0, 0 },
				{ -1, 0, 0, 1, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0 },
				{ -1, 0, 0, 0, 1, 0, 18, -16, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, -1, 1, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 3, -7, 4, 0, 0, 0, 0, 0 },
				{ -2, 0, 1, 1, 1, 0, 0, -3, 7, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 2, 0, 0, -1, 0, -2, 5, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, 0, 0, -2, 5, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, -4, 8, -3, 0, 0, 0, 0 },
				{ 1, 0, 0, 0, 1, 0, -10, 3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, 0, -2, 0, 0, 0, 0, 0, 0 },
				{ -1, 0, 0, 0, 1, 0, 10, -3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, 4, -8, 3, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, -5, 0, 0, 0 },
				{ 0, 0, -1, 1, 0, 0, 0, 1, 0, 2, -5, 0, 0, 0 },
				{ 2, 0, -1, -1, 1, 0, 0, 3, -7, 0, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 0, 0, 0, 2, 0, 0, -5, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, -3, 7, -4, 0, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0 },
				{ 1, 0, 0, 0, 1, 0, -18, 16, 0, 0, 0, 0, 0, 0 },
				{ -2, 0, 1, 1, 1, 0, 0, 1, 0, -2, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 2, 0, -8, 12, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, -8, 13, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, 0, 0, -2, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 0, -2, 2, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -1, 2, 0, 0, 0, 0, 1 },
				{ -1, 0, 0, 1, 1, 0, 3, -4, 0, 0, 0, 0, 0, 0 },
				{ -1, 0, 0, 1, 1, 0, 0, 3, -4, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, 0, -2, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, 0, 2, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2 },
				{ 0, 0, 1, -1, 0, 0, 3, -6, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, -3, 5, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 2, 0, -3, 4, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, -2, 4, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, -5, 6, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, -1, 1, 0, 0, 5, -7, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 5, -8, 0, 0, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 1, 0, 6, -8, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, -8, 15, 0, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 1, 0, 0, 2, 0, -3, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 1, 0, 0, 6, -8, 0, 0, 0, 0, 0 },
				{ 1, 0, 0, -1, 1, 0, 0, -1, 0, 1, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, -5, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, -1, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, 1, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2 },
				{ 0, 0, 1, -1, 2, 0, 0, -1, 0, 0, -1, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0 },
				{ 0, 0, -1, 1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -7, 13, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 7, -13, 0, 0, 0, 0, 0 },
				{ 2, 0, 0, -2, 1, 0, 0, -5, 6, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, 0, -8, 11, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, -1, 0, 2, 0, 0, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 0, 0, 0, 4, -4, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, 0, 3, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 2 },
				{ -2, 0, 0, 2, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 2, 0, 0, -4, 8, -3, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 2, 0, 0, 4, -8, 3, 0, 0, 0, 0 },
				{ 2, 0, 0, -2, 1, 0, 0, -2, 0, 2, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 2, 0, 0, -1, 0, 2, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 2, 0, 0, 0, -2, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, 1, -2, 0, 0, 0, 0, 0 },
				{ 0, 0, -1, 1, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0 },
				{ 0, 0, -1, 1, 0, 0, 0, 1, 0, 0, -2, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, 0, -2, 0, 0, 2, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, 3, -6, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 3, -5, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 3, -5, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 1, 0, -3, 4, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -3, 5, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, -3, 5, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 2, -2, 2, 0, -3, 3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -3, 5, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, 0, 1, -4, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -2, 4, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, 0, -3, 4, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -2, 4, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, -2, 4, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -5, 8, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 2, -2, 2, 0, -5, 6, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -5, 8, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -5, 8, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, -5, 7, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -5, 8, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 5, -8, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 2, 0, 0, -1, 0, -1, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0 },
				{ 0, 0, -1, 1, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, 0, -2, 0, 1, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -6, 11, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 6, -11, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, -1, 0, 4, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 1, 0, -4, 0, 0, 0, 0, 0, 0 },
				{ 2, 0, 0, -2, 1, 0, -3, 3, 0, 0, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, 0, -7, 9, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -5, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, 2, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2 },
				{ 0, 0, 2, -2, 2, 0, 0, -2, 0, 2, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 2 },
				{ 0, 0, 0, 0, 1, 0, 3, -5, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, -1, 1, 0, 0, 3, -4, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, -3, 3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, 2, -4, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, 0, -4, 4, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 2, 0, -5, 7, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, -6, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -3, 6, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, 0, -4, 6, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -3, 6, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, -3, 6, 0, 0, 0, 0, 2 },
				{ 0, 0, -1, 1, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 2, -3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -5, 9, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -5, 9, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 5, -9, 0, 0, 0, 0, 0 },
				{ 0, 0, -1, 1, 0, 0, 0, 1, 0, -2, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, 0, -2, 0, 2, 0, 0, 0, 0 },
				{ -2, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, -2, 2, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -6, 10, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, -6, 10, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -2, 3, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -2, 3, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, -2, 2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 2, -3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 2, -3, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, 0, -1, 0, 3, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, -8, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -4, 8, 0, 0, 0, 0, 2 },
				{ 0, 0, -2, 2, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -4, 7, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -4, 7, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, -7, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, -2, 3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, 0, -2, 0, 3, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -5, 10, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 1, 0, -1, 2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -3, 5, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -3, 5, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, -5, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, 1, -3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -1, 2, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, -1, 2, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -7, 11, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -7, 11, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, -2, 2, 0, 0, 4, -4, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, -3, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, -4, 4, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, -1, 1, 0, 0, 4, -5, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -4, 7, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, -4, 6, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -4, 7, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -4, 6, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -4, 6, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, -4, 5, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -4, 6, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 4, -6, 0, 0, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 },
				{ 0, 0, -1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 1, -1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -1, 0, 5, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, -3, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -1, 3, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -7, 12, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, 1, -2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -2, 5, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -1, 0, 4, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, -4, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, -1, 1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -6, 10, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -6, 10, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, 0, -3, 0, 3, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -3, 7, 0, 0, 0, 0, 2 },
				{ -2, 0, 0, 2, 0, 0, 4, -4, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -5, 8, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 5, -8, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -1, 0, 3, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -1, 0, 3, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, -3, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -2, 4, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, -2, 3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -2, 4, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -6, 9, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -6, 9, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 6, -9, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, 1, 0, -2, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, -2, 2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -4, 6, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, -6, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 3, -4, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -1, 0, 2, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, -2, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, 1, 0, -1, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -5, 9, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, -4, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -3, 4, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -3, 4, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 3, -4, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 3, -4, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 1, 0, 0, 2, -2, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, -1, 0, 2, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -3, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, -5, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, -3, 5, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, -3, 4, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -2, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, -1, 0, 1, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, -2, 2, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -8, 14, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, -5, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 5, -8, 3, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 5, -8, 3, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, -8, 3, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -3, 8, -3, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, -2, 5, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -8, 12, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -8, 12, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, -2, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 2 },
				{ 0, 0, 2, -2, 1, 0, -5, 5, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 3, -6, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -3, 6, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, -3, 6, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -1, 4, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -5, 7, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -5, 7, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, -5, 6, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 5, -7, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, 0, -1, 0, 1, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, -1, 0, 3, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -2, 6, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 1, 0, 2, -2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -6, 9, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 6, -9, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -2, 2, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, -2, 1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -5, 7, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 5, -7, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, -2, 2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, -5, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 1, -3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -1, 3, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, -1, 2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -1, 3, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -7, 10, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -7, 10, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -4, 8, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -4, 5, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -4, 5, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 4, -5, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -2, 0, 5, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -9, 13, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -1, 5, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -2, 0, 4, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, -4, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -2, 7, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, -3, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -2, 5, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, -2, 5, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -6, 8, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -6, 8, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 6, -8, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, 2, 0, -2, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -3, 9, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 5, -6, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 5, -6, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -5, 10, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, -3, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, -5, 13, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, -1, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, -2, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, -2, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, -1, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -6, 15, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -8, 15, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -3, 9, -4, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, -5, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -2, 8, -1, -5, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 6, -8, 3, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 1, -1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -6, 16, -4, -5, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -2, 8, -3, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -2, 8, -3, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 6, -8, 1, 5, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 5, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 3, -5, 4, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -8, 11, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -8, 11, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, -8, 11, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 11, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 3, -3, 0, 2, 0, 0, 0, 2 },
				{ 0, 0, 2, -2, 1, 0, 0, 4, -8, 3, 0, 0, 0, 0 },
				{ 0, 0, 1, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, 0, -4, 8, -3, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -3, 7, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -5, 6, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -5, 6, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 5, -6, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 5, -6, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -1, 6, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 7, -9, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 6, -7, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 5, -5, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -1, 4, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, -1, 4, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -7, 9, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -7, 9, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, -3, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, -1, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -4, 4, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -3, 0, 5, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -9, 12, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, 0, -4, 0, 0, 0, 0 },
				{ 0, 0, 2, -2, 1, 0, 1, -1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 7, -8, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, 0, -3, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, 0, -3, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -2, 6, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -6, 7, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 6, -7, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, 0, -2, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, 0, -2, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 5, -4, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 3, -2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 3, -2, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, 0, -1, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, 0, -1, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, -2, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, -2, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, -1, 0, 0, 2 },
				{ 0, 0, 2, -2, 1, 0, 0, 1, 0, -1, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, -8, 16, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, 0, 2, -5, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 7, -8, 3, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -5, 16, -4, -5, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, -1, 8, -3, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -8, 10, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -8, 10, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, -8, 10, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 3, 0, 1, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -3, 8, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -5, 5, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 5, -5, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 5, -5, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 5, -5, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 7, -7, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 7, -7, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 6, -5, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 7, -8, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 5, -3, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 4, -3, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -9, 11, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -9, 11, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, 0, -4, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, 0, -3, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -6, 6, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, 0, -2, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 6, -4, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 3, -1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 3, -1, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 3, -1, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, 0, -1, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, -2, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 5, -2, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 8, -9, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 5, -4, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, -7, 7, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 7, -7, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 4, -2, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 4, -2, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 4, -2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 4, -2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 5, 0, -4, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 5, 0, -3, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 5, 0, -2, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -8, 8, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 8, -8, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 5, -3, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 5, -3, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, -9, 9, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, -9, 9, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, -9, 9, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 9, -9, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 6, -4, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 2 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2 },
				{ 1, 0, 0, -2, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0 },
				{ 1, 0, 0, -2, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0 },
				{ 1, 0, 0, -2, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0 },
				{ 1, 0, 0, -2, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0 },
				{ -1, 0, 0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0 },
				{ -1, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0 },
				{ -1, 0, 0, 2, 0, 0, 0, 4, -8, 3, 0, 0, 0, 0 },
				{ 1, 0, 0, -2, 0, 0, 0, 4, -8, 3, 0, 0, 0, 0 },
				{ -2, 0, 0, 2, 0, 0, 0, 4, -8, 3, 0, 0, 0, 0 },
				{ -1, 0, 0, 0, 0, 0, 0, 2, 0, -3, 0, 0, 0, 0 },
				{ -1, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0 },
				{ -1, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0 },
				{ -1, 0, 0, 2, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0 },
				{ 1, 0, -1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
				{ -1, 0, 0, 2, 0, 0, 0, 2, 0, -3, 0, 0, 0, 0 },
				{ -2, 0, 0, 0, 0, 0, 0, 2, 0, -3, 0, 0, 0, 0 },
				{ 1, 0, 0, 0, 0, 0, 0, 4, -8, 3, 0, 0, 0, 0 },
				{ -1, 0, 1, -1, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0 },
				{ 1, 0, 1, -1, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0 },
				{ -1, 0, 0, 0, 0, 0, 0, 4, -8, 3, 0, 0, 0, 0 },
				{ -1, 0, 0, 2, 1, 0, 0, 2, 0, -2, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0 },
				{ -1, 0, 0, 2, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0 },
				{ -1, 0, 0, 2, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0 },
				{ 1, 0, 0, -2, 1, 0, 0, -2, 0, 2, 0, 0, 0, 0 },
				{ 1, 0, 2, -2, 2, 0, -3, 3, 0, 0, 0, 0, 0, 0 },
				{ 1, 0, 2, -2, 2, 0, 0, -2, 0, 2, 0, 0, 0, 0 },
				{ 1, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0 },
				{ 1, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0 },
				{ 0, 0, 0, -2, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, -2, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0 },
				{ 0, 0, 2, 0, 2, 0, -2, 2, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, 0, 2, 0, 0, -1, 0, 1, 0, 0, 0, 0 },
				{ 0, 0, 2, 0, 2, 0, -1, 1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, 0, 2, 0, -2, 3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 2, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0 },
				{ 0, 0, 1, 1, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
				{ 1, 0, 2, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
				{ -1, 0, 2, 0, 2, 0, 10, -3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
				{ 1, 0, 2, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, 0, 2, 0, 0, 4, -8, 3, 0, 0, 0, 0 },
				{ 0, 0, 2, 0, 2, 0, 0, -4, 8, -3, 0, 0, 0, 0 },
				{ -1, 0, 2, 0, 2, 0, 0, -4, 8, -3, 0, 0, 0, 0 },
				{ 2, 0, 2, -2, 2, 0, 0, -2, 0, 3, 0, 0, 0, 0 },
				{ 1, 0, 2, 0, 1, 0, 0, -2, 0, 3, 0, 0, 0, 0 },
				{ 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
				{ -1, 0, 2, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
				{ -2, 0, 2, 2, 2, 0, 0, 2, 0, -2, 0, 0, 0, 0 },
				{ 0, 0, 2, 0, 2, 0, 2, -3, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, 0, 2, 0, 1, -1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, 0, 2, 0, 0, 1, 0, -1, 0, 0, 0, 0 },
				{ 0, 0, 2, 0, 2, 0, 2, -2, 0, 0, 0, 0, 0, 0 },
				{ -1, 0, 2, 2, 2, 0, 0, -1, 0, 1, 0, 0, 0, 0 },
				{ 1, 0, 2, 0, 2, 0, -1, 1, 0, 0, 0, 0, 0, 0 },
				{ -1, 0, 2, 2, 2, 0, 0, 2, 0, -3, 0, 0, 0, 0 },
				{ 2, 0, 2, 0, 2, 0, 0, 2, 0, -3, 0, 0, 0, 0 },
				{ 1, 0, 2, 0, 2, 0, 0, -4, 8, -3, 0, 0, 0, 0 },
				{ 1, 0, 2, 0, 2, 0, 0, 4, -8, 3, 0, 0, 0, 0 },
				{ 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
				{ 2, 0, 2, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
				{ -1, 0, 2, 2, 2, 0, 0, 2, 0, -2, 0, 0, 0, 0 },
				{ -1, 0, 2, 2, 2, 0, 3, -3, 0, 0, 0, 0, 0, 0 },
				{ 1, 0, 2, 0, 2, 0, 1, -1, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 2, 2, 2, 0, 0, 2, 0, -2, 0, 0, 0, 0 }
		};
		for (I = 0; I < iNPL; I++)
			for (J = 0; J < 14; J++)
			{
				iNAPL[J][I] = iNAPLT[I][J];
			}

		//ICPLT is the transposition of ICPL
		int iICPLT[iNPL][4] = {
				{ 1440, 0, 0, 0 },
				{ 56, -117, -42, -40 },
				{ 125, -43, 0, -54 },
				{ 0, 5, 0, 0 },
				{ 3, -7, -3, 0 },
				{ 3, 0, 0, -2 },
				{ -114, 0, 0, 61 },
				{ -219, 89, 0, 0 },
				{ -3, 0, 0, 0 },
				{ -462, 1604, 0, 0 },
				{ 99, 0, 0, -53 },
				{ -3, 0, 0, 2 },
				{ 0, 6, 2, 0 },
				{ 3, 0, 0, 0 },
				{ -12, 0, 0, 0 },
				{ 14, -218, 117, 8 },
				{ 31, -481, -257, -17 },
				{ -491, 128, 0, 0 },
				{ -3084, 5123, 2735, 1647 },
				{ -1444, 2409, -1286, -771 },
				{ 11, -24, -11, -9 },
				{ 26, -9, 0, 0 },
				{ 103, -60, 0, 0 },
				{ 0, -13, -7, 0 },
				{ -26, -29, -16, 14 },
				{ 9, -27, -14, -5 },
				{ 12, 0, 0, -6 },
				{ -7, 0, 0, 0 },
				{ 0, 24, 0, 0 },
				{ 284, 0, 0, -151 },
				{ 226, 101, 0, 0 },
				{ 0, -8, -2, 0 },
				{ 0, -6, -3, 0 },
				{ 5, 0, 0, -3 },
				{ -41, 175, 76, 17 },
				{ 0, 15, 6, 0 },
				{ 425, 212, -133, 269 },
				{ 1200, 598, 319, -641 },
				{ 235, 334, 0, 0 },
				{ 11, -12, -7, -6 },
				{ 5, -6, 3, 3 },
				{ -5, 0, 0, 3 },
				{ 6, 0, 0, -3 },
				{ 15, 0, 0, 0 },
				{ 13, 0, 0, -7 },
				{ -6, -9, 0, 0 },
				{ 266, -78, 0, 0 },
				{ -460, -435, -232, 246 },
				{ 0, 15, 7, 0 },
				{ -3, 0, 0, 2 },
				{ 0, 131, 0, 0 },
				{ 4, 0, 0, 0 },
				{ 0, 3, 0, 0 },
				{ 0, 4, 2, 0 },
				{ 0, 3, 0, 0 },
				{ -17, -19, -10, 9 },
				{ -9, -11, 6, -5 },
				{ -6, 0, 0, 3 },
				{ -16, 8, 0, 0 },
				{ 0, 3, 0, 0 },
				{ 11, 24, 11, -5 },
				{ -3, -4, -2, 1 },
				{ 3, 0, 0, -1 },
				{ 0, -8, -4, 0 },
				{ 0, 3, 0, 0 },
				{ 0, 5, 0, 0 },
				{ 0, 3, 2, 0 },
				{ -6, 4, 2, 3 },
				{ -3, -5, 0, 0 },
				{ -5, 0, 0, 2 },
				{ 4, 24, 13, -2 },
				{ -42, 20, 0, 0 },
				{ -10, 233, 0, 0 },
				{ -3, 0, 0, 1 },
				{ 78, -18, 0, 0 },
				{ 0, 3, 1, 0 },
				{ 0, -3, -1, 0 },
				{ 0, -4, -2, 1 },
				{ 0, -8, -4, -1 },
				{ 0, -5, 3, 0 },
				{ -7, 0, 0, 3 },
				{ -14, 8, 3, 6 },
				{ 0, 8, -4, 0 },
				{ 0, 19, 10, 0 },
				{ 45, -22, 0, 0 },
				{ -3, 0, 0, 0 },
				{ 0, -3, 0, 0 },
				{ 0, 3, 0, 0 },
				{ 3, 5, 3, -2 },
				{ 89, -16, -9, -48 },
				{ 0, 3, 0, 0 },
				{ -3, 7, 4, 2 },
				{ -349, -62, 0, 0 },
				{ -15, 22, 0, 0 },
				{ -3, 0, 0, 0 },
				{ -53, 0, 0, 0 },
				{ 5, 0, 0, -3 },
				{ 0, -8, 0, 0 },
				{ 15, -7, -4, -8 },
				{ -3, 0, 0, 1 },
				{ -21, -78, 0, 0 },
				{ 20, -70, -37, -11 },
				{ 0, 6, 3, 0 },
				{ 5, 3, 2, -2 },
				{ -17, -4, -2, 9 },
				{ 0, 6, 3, 0 },
				{ 32, 15, -8, 17 },
				{ 174, 84, 45, -93 },
				{ 11, 56, 0, 0 },
				{ -66, -12, -6, 35 },
				{ 47, 8, 4, -25 },
				{ 0, 8, 4, 0 },
				{ 10, -22, -12, -5 },
				{ -3, 0, 0, 2 },
				{ -24, 12, 0, 0 },
				{ 5, -6, 0, 0 },
				{ 3, 0, 0, -2 },
				{ 4, 3, 1, -2 },
				{ 0, 29, 15, 0 },
				{ -5, -4, -2, 2 },
				{ 8, -3, -1, -5 },
				{ 0, -3, 0, 0 },
				{ 10, 0, 0, 0 },
				{ 3, 0, 0, -2 },
				{ -5, 0, 0, 3 },
				{ 46, 66, 35, -25 },
				{ -14, 7, 0, 0 },
				{ 0, 3, 2, 0 },
				{ -5, 0, 0, 0 },
				{ -68, -34, -18, 36 },
				{ 0, 14, 7, 0 },
				{ 10, -6, -3, -5 },
				{ -5, -4, -2, 3 },
				{ -3, 5, 2, 1 },
				{ 76, 17, 9, -41 },
				{ 84, 298, 159, -45 },
				{ 3, 0, 0, -1 },
				{ -3, 0, 0, 2 },
				{ -3, 0, 0, 1 },
				{ -82, 292, 156, 44 },
				{ -73, 17, 9, 39 },
				{ -9, -16, 0, 0 },
				{ 3, 0, -1, -2 },
				{ -3, 0, 0, 0 },
				{ -9, -5, -3, 5 },
				{ -439, 0, 0, 0 },
				{ 57, -28, -15, -30 },
				{ 0, -6, -3, 0 },
				{ -4, 0, 0, 2 },
				{ -40, 57, 30, 21 },
				{ 23, 7, 3, -13 },
				{ 273, 80, 43, -146 },
				{ -449, 430, 0, 0 },
				{ -8, -47, -25, 4 },
				{ 6, 47, 25, -3 },
				{ 0, 23, 13, 0 },
				{ -3, 0, 0, 2 },
				{ 3, -4, -2, -2 },
				{ -48, -110, -59, 26 },
				{ 51, 114, 61, -27 },
				{ -133, 0, 0, 57 },
				{ 0, 4, 0, 0 },
				{ -21, -6, -3, 11 },
				{ 0, -3, -1, 0 },
				{ -11, -21, -11, 6 },
				{ -18, -436, -233, 9 },
				{ 35, -7, 0, 0 },
				{ 0, 5, 3, 0 },
				{ 11, -3, -1, -6 },
				{ -5, -3, -1, 3 },
				{ -53, -9, -5, 28 },
				{ 0, 3, 2, 1 },
				{ 4, 0, 0, -2 },
				{ 0, -4, 0, 0 },
				{ -50, 194, 103, 27 },
				{ -13, 52, 28, 7 },
				{ -91, 248, 0, 0 },
				{ 6, 49, 26, -3 },
				{ -6, -47, -25, 3 },
				{ 0, 5, 3, 0 },
				{ 52, 23, 10, -23 },
				{ -3, 0, 0, 1 },
				{ 0, 5, 3, 0 },
				{ -4, 0, 0, 0 },
				{ -4, 8, 3, 2 },
				{ 10, 0, 0, 0 },
				{ 3, 0, 0, -2 },
				{ 0, 8, 4, 0 },
				{ 0, 8, 4, 1 },
				{ -4, 0, 0, 0 },
				{ -4, 0, 0, 0 },
				{ -8, 4, 2, 4 },
				{ 8, -4, -2, -4 },
				{ 0, 15, 7, 0 },
				{ -138, 0, 0, 0 },
				{ 0, -7, -3, 0 },
				{ 0, -7, -3, 0 },
				{ 54, 0, 0, -29 },
				{ 0, 10, 4, 0 },
				{ -7, 0, 0, 3 },
				{ -37, 35, 19, 20 },
				{ 0, 4, 0, 0 },
				{ -4, 9, 0, 0 },
				{ 8, 0, 0, -4 },
				{ -9, -14, -8, 5 },
				{ -3, -9, -5, 3 },
				{ -145, 47, 0, 0 },
				{ -10, 40, 21, 5 },
				{ 11, -49, -26, -7 },
				{ -2150, 0, 0, 932 },
				{ -12, 0, 0, 5 },
				{ 85, 0, 0, -37 },
				{ 4, 0, 0, -2 },
				{ 3, 0, 0, -2 },
				{ -86, 153, 0, 0 },
				{ -6, 9, 5, 3 },
				{ 9, -13, -7, -5 },
				{ -8, 12, 6, 4 },
				{ -51, 0, 0, 22 },
				{ -11, -268, -116, 5 },
				{ 0, 12, 5, 0 },
				{ 0, 7, 3, 0 },
				{ 31, 6, 3, -17 },
				{ 140, 27, 14, -75 },
				{ 57, 11, 6, -30 },
				{ -14, -39, 0, 0 },
				{ 0, -6, -2, 0 },
				{ 4, 15, 8, -2 },
				{ 0, 4, 0, 0 },
				{ -3, 0, 0, 1 },
				{ 0, 11, 5, 0 },
				{ 9, 6, 0, 0 },
				{ -4, 10, 4, 2 },
				{ 5, 3, 0, 0 },
				{ 16, 0, 0, -9 },
				{ -3, 0, 0, 0 },
				{ 0, 3, 2, -1 },
				{ 7, 0, 0, -3 },
				{ -25, 22, 0, 0 },
				{ 42, 223, 119, -22 },
				{ -27, -143, -77, 14 },
				{ 9, 49, 26, -5 },
				{ -1166, 0, 0, 505 },
				{ -5, 0, 0, 2 },
				{ -6, 0, 0, 3 },
				{ -8, 0, 1, 4 },
				{ 0, -4, 0, 0 },
				{ 117, 0, 0, -63 },
				{ -4, 8, 4, 2 },
				{ 3, 0, 0, -2 },
				{ -5, 0, 0, 2 },
				{ 0, 31, 0, 0 },
				{ -5, 0, 1, 3 },
				{ 4, 0, 0, -2 },
				{ -4, 0, 0, 2 },
				{ -24, -13, -6, 10 },
				{ 3, 0, 0, 0 },
				{ 0, -32, -17, 0 },
				{ 8, 12, 5, -3 },
				{ 3, 0, 0, -1 },
				{ 7, 13, 0, 0 },
				{ -3, 16, 0, 0 },
				{ 50, 0, 0, -27 },
				{ 0, -5, -3, 0 },
				{ 13, 0, 0, 0 },
				{ 0, 5, 3, 1 },
				{ 24, 5, 2, -11 },
				{ 5, -11, -5, -2 },
				{ 30, -3, -2, -16 },
				{ 18, 0, 0, -9 },
				{ 8, 614, 0, 0 },
				{ 3, -3, -1, -2 },
				{ 6, 17, 9, -3 },
				{ -3, -9, -5, 2 },
				{ 0, 6, 3, -1 },
				{ -127, 21, 9, 55 },
				{ 3, 5, 0, 0 },
				{ -6, -10, -4, 3 },
				{ 5, 0, 0, 0 },
				{ 16, 9, 4, -7 },
				{ 3, 0, 0, -2 },
				{ 0, 22, 0, 0 },
				{ 0, 19, 10, 0 },
				{ 7, 0, 0, -4 },
				{ 0, -5, -2, 0 },
				{ 0, 3, 1, 0 },
				{ -9, 3, 1, 4 },
				{ 17, 0, 0, -7 },
				{ 0, -3, -2, -1 },
				{ -20, 34, 0, 0 },
				{ -10, 0, 1, 5 },
				{ -4, 0, 0, 2 },
				{ 22, -87, 0, 0 },
				{ -4, 0, 0, 2 },
				{ -3, -6, -2, 1 },
				{ -16, -3, -1, 7 },
				{ 0, -3, -2, 0 },
				{ 4, 0, 0, 0 },
				{ -68, 39, 0, 0 },
				{ 27, 0, 0, -14 },
				{ 0, -4, 0, 0 },
				{ -25, 0, 0, 0 },
				{ -12, -3, -2, 6 },
				{ 3, 0, 0, -1 },
				{ 3, 66, 29, -1 },
				{ 490, 0, 0, -213 },
				{ -22, 93, 49, 12 },
				{ -7, 28, 15, 4 },
				{ -3, 13, 7, 2 },
				{ -46, 14, 0, 0 },
				{ -5, 0, 0, 0 },
				{ 2, 1, 0, 0 },
				{ 0, -3, 0, 0 },
				{ -28, 0, 0, 15 },
				{ 5, 0, 0, -2 },
				{ 0, 3, 0, 0 },
				{ -11, 0, 0, 5 },
				{ 0, 3, 1, 0 },
				{ -3, 0, 0, 1 },
				{ 25, 106, 57, -13 },
				{ 5, 21, 11, -3 },
				{ 1485, 0, 0, 0 },
				{ -7, -32, -17, 4 },
				{ 0, 5, 3, 0 },
				{ -6, -3, -2, 3 },
				{ 30, -6, -2, -13 },
				{ -4, 4, 0, 0 },
				{ -19, 0, 0, 10 },
				{ 0, 4, 2, -1 },
				{ 0, 3, 0, 0 },
				{ 4, 0, 0, -2 },
				{ 0, -3, -1, 0 },
				{ -3, 0, 0, 0 },
				{ 5, 3, 1, -2 },
				{ 0, 11, 0, 0 },
				{ 118, 0, 0, -52 },
				{ 0, -5, -3, 0 },
				{ -28, 36, 0, 0 },
				{ 5, -5, 0, 0 },
				{ 14, -59, -31, -8 },
				{ 0, 9, 5, 1 },
				{ -458, 0, 0, 198 },
				{ 0, -45, -20, 0 },
				{ 9, 0, 0, -5 },
				{ 0, -3, 0, 0 },
				{ 0, -4, -2, -1 },
				{ 11, 0, 0, -6 },
				{ 6, 0, 0, -2 },
				{ -16, 23, 0, 0 },
				{ 0, -4, -2, 0 },
				{ -5, 0, 0, 2 },
				{ -166, 269, 0, 0 },
				{ 15, 0, 0, -8 },
				{ 10, 0, 0, -4 },
				{ -78, 45, 0, 0 },
				{ 0, -5, -2, 0 },
				{ 7, 0, 0, -4 },
				{ -5, 328, 0, 0 },
				{ 3, 0, 0, -2 },
				{ 5, 0, 0, -2 },
				{ 0, 3, 1, 0 },
				{ -3, 0, 0, 0 },
				{ -3, 0, 0, 0 },
				{ 0, -4, -2, 0 },
				{ -1223, -26, 0, 0 },
				{ 0, 7, 3, 0 },
				{ 3, 0, 0, 0 },
				{ 0, 3, 2, 0 },
				{ -6, 20, 0, 0 },
				{ -368, 0, 0, 0 },
				{ -75, 0, 0, 0 },
				{ 11, 0, 0, -6 },
				{ 3, 0, 0, -2 },
				{ -3, 0, 0, 1 },
				{ -13, -30, 0, 0 },
				{ 21, 3, 0, 0 },
				{ -3, 0, 0, 1 },
				{ -4, 0, 0, 2 },
				{ 8, -27, 0, 0 },
				{ -19, -11, 0, 0 },
				{ -4, 0, 0, 2 },
				{ 0, 5, 2, 0 },
				{ -6, 0, 0, 2 },
				{ -8, 0, 0, 0 },
				{ -1, 0, 0, 0 },
				{ -14, 0, 0, 6 },
				{ 6, 0, 0, 0 },
				{ -74, 0, 0, 32 },
				{ 0, -3, -1, 0 },
				{ 4, 0, 0, -2 },
				{ 8, 11, 0, 0 },
				{ 0, 3, 2, 0 },
				{ -262, 0, 0, 114 },
				{ 0, -4, 0, 0 },
				{ -7, 0, 0, 4 },
				{ 0, -27, -12, 0 },
				{ -19, -8, -4, 8 },
				{ 202, 0, 0, -87 },
				{ -8, 35, 19, 5 },
				{ 0, 4, 2, 0 },
				{ 16, -5, 0, 0 },
				{ 5, 0, 0, -3 },
				{ 0, -3, 0, 0 },
				{ 1, 0, 0, 0 },
				{ -35, -48, -21, 15 },
				{ -3, -5, -2, 1 },
				{ 6, 0, 0, -3 },
				{ 3, 0, 0, -1 },
				{ 0, -5, 0, 0 },
				{ 12, 55, 29, -6 },
				{ 0, 5, 3, 0 },
				{ -598, 0, 0, 0 },
				{ -3, -13, -7, 1 },
				{ -5, -7, -3, 2 },
				{ 3, 0, 0, -1 },
				{ 5, -7, 0, 0 },
				{ 4, 0, 0, -2 },
				{ 16, -6, 0, 0 },
				{ 8, -3, 0, 0 },
				{ 8, -31, -16, -4 },
				{ 0, 3, 1, 0 },
				{ 113, 0, 0, -49 },
				{ 0, -24, -10, 0 },
				{ 4, 0, 0, -2 },
				{ 27, 0, 0, 0 },
				{ -3, 0, 0, 1 },
				{ 0, -4, -2, 0 },
				{ 5, 0, 0, -2 },
				{ 0, -3, 0, 0 },
				{ -13, 0, 0, 6 },
				{ 5, 0, 0, -2 },
				{ -18, -10, -4, 8 },
				{ -4, -28, 0, 0 },
				{ -5, 6, 3, 2 },
				{ -3, 0, 0, 1 },
				{ -5, -9, -4, 2 },
				{ 17, 0, 0, -7 },
				{ 11, 4, 0, 0 },
				{ 0, -6, -2, 0 },
				{ 83, 15, 0, 0 },
				{ -4, 0, 0, 2 },
				{ 0, -114, -49, 0 },
				{ 117, 0, 0, -51 },
				{ -5, 19, 10, 2 },
				{ -3, 0, 0, 0 },
				{ -3, 0, 0, 2 },
				{ 0, -3, -1, 0 },
				{ 3, 0, 0, 0 },
				{ 0, -6, -2, 0 },
				{ 393, 3, 0, 0 },
				{ -4, 21, 11, 2 },
				{ -6, 0, -1, 3 },
				{ -3, 8, 4, 1 },
				{ 8, 0, 0, 0 },
				{ 18, -29, -13, -8 },
				{ 8, 34, 18, -4 },
				{ 89, 0, 0, 0 },
				{ 3, 12, 6, -1 },
				{ 54, -15, -7, -24 },
				{ 0, 3, 0, 0 },
				{ 3, 0, 0, -1 },
				{ 0, 35, 0, 0 },
				{ -154, -30, -13, 67 },
				{ 15, 0, 0, 0 },
				{ 0, 4, 2, 0 },
				{ 0, 9, 0, 0 },
				{ 80, -71, -31, -35 },
				{ 0, -20, -9, 0 },
				{ 11, 5, 2, -5 },
				{ 61, -96, -42, -27 },
				{ 14, 9, 4, -6 },
				{ -11, -6, -3, 5 },
				{ 0, -3, -1, 0 },
				{ 123, -415, -180, -53 },
				{ 0, 0, 0, -35 },
				{ -5, 0, 0, 0 },
				{ 7, -32, -17, -4 },
				{ 0, -9, -5, 0 },
				{ 0, -4, 2, 0 },
				{ -89, 0, 0, 38 },
				{ 0, -86, -19, -6 },
				{ 0, 0, -19, 6 },
				{ -123, -416, -180, 53 },
				{ 0, -3, -1, 0 },
				{ 12, -6, -3, -5 },
				{ -13, 9, 4, 6 },
				{ 0, -15, -7, 0 },
				{ 3, 0, 0, -1 },
				{ -62, -97, -42, 27 },
				{ -11, 5, 2, 5 },
				{ 0, -19, -8, 0 },
				{ -3, 0, 0, 1 },
				{ 0, 4, 2, 0 },
				{ 0, 3, 0, 0 },
				{ 0, 4, 2, 0 },
				{ -85, -70, -31, 37 },
				{ 163, -12, -5, -72 },
				{ -63, -16, -7, 28 },
				{ -21, -32, -14, 9 },
				{ 0, -3, -1, 0 },
				{ 3, 0, 0, -2 },
				{ 0, 8, 0, 0 },
				{ 3, 10, 4, -1 },
				{ 3, 0, 0, -1 },
				{ 0, -7, -3, 0 },
				{ 0, -4, -2, 0 },
				{ 6, 19, 0, 0 },
				{ 5, -173, -75, -2 },
				{ 0, -7, -3, 0 },
				{ 7, -12, -5, -3 },
				{ -3, 0, 0, 2 },
				{ 3, -4, -2, -1 },
				{ 74, 0, 0, -32 },
				{ -3, 12, 6, 2 },
				{ 26, -14, -6, -11 },
				{ 19, 0, 0, -8 },
				{ 6, 24, 13, -3 },
				{ 83, 0, 0, 0 },
				{ 0, -10, -5, 0 },
				{ 11, -3, -1, -5 },
				{ 3, 0, 1, -1 },
				{ 3, 0, 0, -1 },
				{ -4, 0, 0, 0 },
				{ 5, -23, -12, -3 },
				{ -339, 0, 0, 147 },
				{ 0, -10, -5, 0 },
				{ 5, 0, 0, 0 },
				{ 3, 0, 0, -1 },
				{ 0, -4, -2, 0 },
				{ 18, -3, 0, 0 },
				{ 9, -11, -5, -4 },
				{ -8, 0, 0, 4 },
				{ 3, 0, 0, -1 },
				{ 0, 9, 0, 0 },
				{ 6, -9, -4, -2 },
				{ -4, -12, 0, 0 },
				{ 67, -91, -39, -29 },
				{ 30, -18, -8, -13 },
				{ 0, 0, 0, 0 },
				{ 0, -114, -50, 0 },
				{ 0, 0, 0, 23 },
				{ 517, 16, 7, -224 },
				{ 0, -7, -3, 0 },
				{ 143, -3, -1, -62 },
				{ 29, 0, 0, -13 },
				{ -4, 0, 0, 2 },
				{ -6, 0, 0, 3 },
				{ 5, 12, 5, -2 },
				{ -25, 0, 0, 11 },
				{ -3, 0, 0, 1 },
				{ 0, 4, 2, 0 },
				{ -22, 12, 5, 10 },
				{ 50, 0, 0, -22 },
				{ 0, 7, 4, 0 },
				{ 0, 3, 1, 0 },
				{ -4, 4, 2, 2 },
				{ -5, -11, -5, 2 },
				{ 0, 4, 2, 0 },
				{ 4, 17, 9, -2 },
				{ 59, 0, 0, 0 },
				{ 0, -4, -2, 0 },
				{ -8, 0, 0, 4 },
				{ -3, 0, 0, 0 },
				{ 4, -15, -8, -2 },
				{ 370, -8, 0, -160 },
				{ 0, 0, -3, 0 },
				{ 0, 3, 1, 0 },
				{ -6, 3, 1, 3 },
				{ 0, 6, 0, 0 },
				{ -10, 0, 0, 4 },
				{ 0, 9, 4, 0 },
				{ 4, 17, 7, -2 },
				{ 34, 0, 0, -15 },
				{ 0, 5, 3, 0 },
				{ -5, 0, 0, 2 },
				{ -37, -7, -3, 16 },
				{ 3, 13, 7, -2 },
				{ 40, 0, 0, 0 },
				{ 0, -3, -2, 0 },
				{ -184, -3, -1, 80 },
				{ -3, 0, 0, 1 },
				{ -3, 0, 0, 0 },
				{ 0, -10, -6, -1 },
				{ 31, -6, 0, -13 },
				{ -3, -32, -14, 1 },
				{ -7, 0, 0, 3 },
				{ 0, -8, -4, 0 },
				{ 3, -4, 0, 0 },
				{ 0, 4, 0, 0 },
				{ 0, 3, 1, 0 },
				{ 19, -23, -10, 2 },
				{ 0, 0, 0, -10 },
				{ 0, 3, 2, 0 },
				{ 0, 9, 5, -1 },
				{ 28, 0, 0, 0 },
				{ 0, -7, -4, 0 },
				{ 8, -4, 0, -4 },
				{ 0, 0, -2, 0 },
				{ 0, 3, 0, 0 },
				{ -3, 0, 0, 1 },
				{ -9, 0, 1, 4 },
				{ 3, 12, 5, -1 },
				{ 17, -3, -1, 0 },
				{ 0, 7, 4, 0 },
				{ 19, 0, 0, 0 },
				{ 0, -5, -3, 0 },
				{ 14, -3, 0, -1 },
				{ 0, 0, -1, 0 },
				{ 0, 0, 0, -5 },
				{ 0, 5, 3, 0 },
				{ 13, 0, 0, 0 },
				{ 0, -3, -2, 0 },
				{ 2, 9, 4, 3 },
				{ 0, 0, 0, -4 },
				{ 8, 0, 0, 0 },
				{ 0, 4, 2, 0 },
				{ 6, 0, 0, -3 },
				{ 6, 0, 0, 0 },
				{ 0, 3, 1, 0 },
				{ 5, 0, 0, -2 },
				{ 3, 0, 0, -1 },
				{ -3, 0, 0, 0 },
				{ 6, 0, 0, 0 },
				{ 7, 0, 0, 0 },
				{ -4, 0, 0, 0 },
				{ 4, 0, 0, 0 },
				{ 6, 0, 0, 0 },
				{ 0, -4, 0, 0 },
				{ 0, -4, 0, 0 },
				{ 5, 0, 0, 0 },
				{ -3, 0, 0, 0 },
				{ 4, 0, 0, 0 },
				{ -5, 0, 0, 0 },
				{ 4, 0, 0, 0 },
				{ 0, 3, 0, 0 },
				{ 13, 0, 0, 0 },
				{ 21, 11, 0, 0 },
				{ 0, -5, 0, 0 },
				{ 0, -5, -2, 0 },
				{ 0, 5, 3, 0 },
				{ 0, -5, 0, 0 },
				{ -3, 0, 0, 2 },
				{ 20, 10, 0, 0 },
				{ -34, 0, 0, 0 },
				{ -19, 0, 0, 0 },
				{ 3, 0, 0, -2 },
				{ -3, 0, 0, 1 },
				{ -6, 0, 0, 3 },
				{ -4, 0, 0, 0 },
				{ 3, 0, 0, 0 },
				{ 3, 0, 0, 0 },
				{ 4, 0, 0, 0 },
				{ 3, 0, 0, -1 },
				{ 6, 0, 0, -3 },
				{ -8, 0, 0, 3 },
				{ 0, 3, 1, 0 },
				{ -3, 0, 0, 0 },
				{ 0, -3, -2, 0 },
				{ 126, -63, -27, -55 },
				{ -5, 0, 1, 2 },
				{ -3, 28, 15, 2 },
				{ 5, 0, 1, -2 },
				{ 0, 9, 4, 1 },
				{ 0, 9, 4, -1 },
				{ -126, -63, -27, 55 },
				{ 3, 0, 0, -1 },
				{ 21, -11, -6, -11 },
				{ 0, -4, 0, 0 },
				{ -21, -11, -6, 11 },
				{ -3, 0, 0, 1 },
				{ 0, 3, 1, 0 },
				{ 8, 0, 0, -4 },
				{ -6, 0, 0, 3 },
				{ -3, 0, 0, 1 },
				{ 3, 0, 0, -1 },
				{ -3, 0, 0, 1 },
				{ -5, 0, 0, 2 },
				{ 24, -12, -5, -11 },
				{ 0, 3, 1, 0 },
				{ 0, 3, 1, 0 },
				{ 0, 3, 2, 0 },
				{ -24, -12, -5, 10 },
				{ 4, 0, -1, -2 },
				{ 13, 0, 0, -6 },
				{ 7, 0, 0, -3 },
				{ 3, 0, 0, -1 },
				{ 3, 0, 0, -1 }
		};
		for (I = 0; I < iNPL; I++)
			for (J = 0; J < 4; J++)
			{
				iICPL[J][I] = iICPLT[I][J];
			}
		dT = ((dATE1 - dDJ0) + dATE2) / dDJC;

		/***LUNI-SOLAR NUTATION***/

		dEL = fmod(485868.249036e0 +
			dT * (1717915923.2178e0 +
				dT * (31.8792e0 +
					dT * (0.051635e0 +
						dT * (-0.00024470e0)))), dTURNAS) * dDAS2R;
		dELP = fmod(1287104.79305e0 +
			dT * (129596581.0481e0 +
				dT * (-0.5532e0 +
					dT * (0.000136e0 +
						dT * (-0.00001149e0)))), dTURNAS) * dDAS2R;
		dF = fmod(335779.526232e0 +
			dT * (1739527262.8478e0 +
				dT * (-12.7512e0 +
					dT * (-0.001037e0 +
						dT * (0.00000417e0)))), dTURNAS) * dDAS2R;
		dD = fmod(1072260.70369e0 +
			dT * (1602961601.2090e0 +
				dT * (-6.3706e0 +
					dT * (0.006593e0 +
						dT * (-0.00003169e0)))), dTURNAS) * dDAS2R;
		dOM = fmod(450160.398036e0 +
			dT * (-6962890.5431e0 +
				dT * (7.4722e0 +
					dT * (0.007702e0 +
						dT * (-0.00005939e0)))), dTURNAS) * dDAS2R;
		dDP = 0e0;
		dDE = 0e0;
		for (I = iNLS - 1; I >= 0; I--)
		{
			dARG = fmod(iNALS[0][I] * dEL +
				iNALS[1][I] * dELP +
				iNALS[2][I] * dF +
				iNALS[3][I] * dD +
				iNALS[4][I] * dOM, d2PI);
			dSARG = sin(dARG);
			dCARG = cos(dARG);

			dDP = dDP + (dCLS[0][I] + dCLS[1][I] * dT) * dSARG + dCLS[2][I] * dCARG;
			dDE = dDE + (dCLS[3][I] + dCLS[4][I] * dT) * dCARG + dCLS[5][I] * dSARG;
		}
		dDPSILS = dDP * dU2R;
		dDEPSLS = dDE * dU2R;
		dAL = fmod(2.35555598e0 + 8328.6914269554e0 * dT, d2PI);
		dALSU = fmod(6.24006013e0 + 628.301955e0 * dT, d2PI);
		dAF = fmod(1.627905234e0 + 8433.466158131e0 * dT, d2PI);
		dAD = fmod(5.198466741e0 + 7771.3771468121e0 * dT, d2PI);
		dAOM = fmod(2.18243920e0 - 33.757045e0 * dT, d2PI);
		dAPA = (0.02438175e0 + 0.00000538691e0 * dT) * dT;
		dALME = fmod(4.402608842e0 + 2608.7903141574e0 * dT, d2PI);
		dALVE = fmod(3.176146697e0 + 1021.3285546211e0 * dT, d2PI);
		dALEA = fmod(1.753470314e0 + 628.3075849991e0 * dT, d2PI);
		dALMA = fmod(6.203480913e0 + 334.0612426700e0 * dT, d2PI);
		dALJU = fmod(0.599546497e0 + 52.9690962641e0 * dT, d2PI);
		dALSA = fmod(0.874016757e0 + 21.3299104960e0 * dT, d2PI);
		dALUR = fmod(5.481293871e0 + 7.4781598567e0 * dT, d2PI);
		dALNE = fmod(5.321159000e0 + 3.8127774000e0 * dT, d2PI);

		dDP = 0e0;
		dDE = 0e0;

		for (I = iNPL - 1; I >= 0; I--)
		{
			dARG = fmod(iNAPL[0][I] * dAL +
				iNAPL[1][I] * dALSU +
				iNAPL[2][I] * dAF +
				iNAPL[3][I] * dAD +
				iNAPL[4][I] * dAOM +
				iNAPL[5][I] * dALME +
				iNAPL[6][I] * dALVE +
				iNAPL[7][I] * dALEA +
				iNAPL[8][I] * dALMA +
				iNAPL[9][I] * dALJU +
				iNAPL[10][I] * dALSA +
				iNAPL[11][I] * dALUR +
				iNAPL[12][I] * dALNE +
				iNAPL[13][I] * dAPA, d2PI);
			dSARG = sin(dARG);
			dCARG = cos(dARG);

			dDP = dDP + iICPL[0][I] * dSARG + iICPL[1][I] * dCARG;
			dDE = dDE + iICPL[2][I] * dSARG + iICPL[3][I] * dCARG;
		}
		dDPSIPL = dDP * dU2R;
		dDEPSPL = dDE * dU2R;

		double temp1 = dDPSIPL + dDPSILS;
		*dPSI = temp1;
		double temp2 = dDEPSPL + dDEPSLS;
		*dEPS = temp2;
		delete[] iNALS;
		delete[] dCLS;
		delete[] iNAPL;
		delete[] iICPL;//modified by zhuyt
	}

	
}