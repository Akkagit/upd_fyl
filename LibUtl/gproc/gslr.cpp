/**
*
* @verbatim
	History
	 -1.0 Hongmin Zhang		2019-10-25  creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gslr.cpp
* @brief		Satellite laser ranging
*
* @author       Hongmin Zhang, Wuhan University
* @version		1.0.0
* @date		    2019-10-25
*
*/


#include <algorithm>
#include <set>
#include "gslr.h"

#include "gset/gsetout.h" 
#include "gset/gsetorbdif.h"
#include "gutils/gtrs2crs.h"
#include "gutils/gcrs2acr.h"
#include "gutils/gleo.h"

#include "gutils/gtrs2crs.h"
#include "gdata/gnavde.h"
#include "gdata/gresult.h"
#include "gall/gallsta.h"

using namespace std;

namespace great
{
	t_gslr::t_gslr()
	{
	}

	t_gslr::t_gslr(t_gsetbase* set, t_glog* log) 
	{
		_set = set;
		_glog = log;
		_tide = new t_gtideIERS(log);
		//setOBJ(dynamic_cast<t_gallobj*>((*data)[t_gdata::ALLOBJ]));
		_result = new t_gresult_one();
		_gresult = new t_gresult();
		_gpreproc = new t_gpreprocslr();
		//_result = new t_gresult();
	}

	t_gslr::~t_gslr()
	{

	}

	void t_gslr::setDataUsed(string sat)
	{
		//string temp;
		//temp = strl2strs(sat);
		_sat = sat;

		
	}

	void t_gslr::setTime(t_gtime now)
	{
		_time = now;
	}

	void t_gslr::setOBJ(t_gallobj* gobj)
	{
		_gallobj = gobj;

		string _crt_rec = strs2strl(_sat);
		//_crt_sys = obsdata.gsys();
		_crt_obj = _gallobj->obj(_crt_rec);

	}

	void t_gslr::setifLEO(string str)
	{
		_ifLEO = str;
	}

	void t_gslr::setATT(t_gattitude*    attitude)
	{
		_gattitude = attitude;
	}

	void t_gslr::setPoleut(t_gpoleut1* poleut1)
	{
		_gpoleut1 = poleut1;
		_trs2crs_2000 = new t_gtrs2crs(false, _gpoleut1);
	}

	void t_gslr::setNPT(t_gallnpt* gnpt)
	{
		_gnpt = gnpt;
	}

	void t_gslr::setLRA(t_galllra* glra)
	{
		if (glra->lrasize() == 0)
		{
			_glra = NULL;
			return;
		}
		_glra = glra;
	}

	void t_gslr::setSLRF(t_gslrf* gslrf)
	{
		_gslrf = gslrf;
	}

	void t_gslr::setDE(t_gnavde* gde)
	{
		_gde = gde;
	}

	void t_gslr::setORB(t_gallorb* gorb)
	{
		_gorb = gorb;
	}

	void t_gslr::setOTL(t_gallotl* otl)
	{
		if (_tide)
		{
			_tide->setOTL(otl);
		}
	}


	bool t_gslr::processBatch()
	{

		string function_name = "t_gslr::processBatch";
		int m;

		if (!_init_data())
		{
			_glog->comment(0, function_name, "ERROR : init data fail, check your xml file.");
			return false;
		}

		int islr = 0;

		if (_ifLEO == "NAV")
		{
			Matrix G_svn_ilrs;
			Matrix R_svn_ilrs;
			Matrix E_svn_ilrs;
			Matrix C_svn_ilrs;
			Matrix J_svn_ilrs;

			G_svn_ilrs.ReSize(100, 2);
			R_svn_ilrs.ReSize(100, 2);
			E_svn_ilrs.ReSize(100, 2);
			C_svn_ilrs.ReSize(100, 2);
			J_svn_ilrs.ReSize(100, 2);

			G_svn_ilrs = 0;
			R_svn_ilrs = 0;
			E_svn_ilrs = 0;
			C_svn_ilrs = 0;
			J_svn_ilrs = 0;

			G_svn_ilrs(1, 1) = 35;
			G_svn_ilrs(2, 1) = 36;
			G_svn_ilrs(1, 2) = 35;
			G_svn_ilrs(2, 2) = 36;

			for (int i = 1; i <= 34; i++)
			{
				R_svn_ilrs(i, 1) = 714 + i - 1;
				R_svn_ilrs(i, 2) = 100 + i - 1;
			}

			R_svn_ilrs(14, 1) = 728;
			R_svn_ilrs(18, 1) = 733;
			R_svn_ilrs(19, 1) = 734;
			R_svn_ilrs(20, 1) = 731;
			R_svn_ilrs(21, 1) = 732;
			R_svn_ilrs(22, 1) = 735;
			R_svn_ilrs(23, 1) = 736;
			R_svn_ilrs(24, 1) = 737;
			R_svn_ilrs(25, 1) = 738;
			R_svn_ilrs(26, 1) = 701;
			R_svn_ilrs(27, 1) = 742;
			R_svn_ilrs(28, 1) = 743;
			R_svn_ilrs(29, 1) = 744;
			R_svn_ilrs(30, 1) = 745;
			R_svn_ilrs(31, 1) = 746;
			R_svn_ilrs(32, 1) = 747;
			R_svn_ilrs(33, 1) = 754;
			R_svn_ilrs(34, 1) = 755;

			if (_sat[0] == 'C') 
			{
				if (_sat.compare("C01") == 0) // C003
				{
					islr = 1;
					_slre = -0.6241;
					_slrn = -0.5705;
					_slru = 1.1000;
				}
				else if (_sat.compare("C08") == 0) //C008
				{
					islr = 2;
					_slre = -0.4203;
					_slrn = -0.5735;
					_slru =  1.1000;
				}
				else if (_sat.compare("C10") == 0) //C010
				{
					islr = 3;
					_slre = -0.4267;
					_slrn = -0.5732;
					_slru =  1.0995;
				}
				else if (_sat.compare("C11") == 0) //C012
				{
					islr = 4;
					_slre = -0.4264;
					_slrn = -0.5379;
					_slru =  1.1000;					     
				}
				else if (_sat.compare("C13") == 0) //C017
				{
					islr = 4;
					_slre = -0.4215;
					_slrn = -0.5726;
					_slru =  1.1000;
				}
				else if (_sat.compare("C20") == 0) //C202
				{
					islr = 1;
					_slre = 0.5947;
					_slrn =-0.0846;
					_slru = 1.2644;
				}
				else if (_sat.compare("C21") == 0) //C206
				{
					islr = 1;
					_slre = 0.5986;
					_slrn =-0.0866;
					_slru = 1.2650;
				}
				else if (_sat.compare("C29") == 0) //C207
				{
					islr = 1;
					_slre = 0.6095;
					_slrn = 0.4260;
					_slru = 0.6142;
				}
				else if (_sat.compare("C30") == 0) //C208
				{
					islr = 1;
					_slre = 0.6097;
					_slrn = 0.4273;
					_slru = 0.6153;
				}
				else {
					_glog->comment(0, "BDS SVN to ILRS Number:not found" + _sat);
					islr = 0;
				}
			}			
			else if (_sat[0] == 'G')
			{
				islr = 5;
				_slre = 0.86260;
				_slrn = -0.52450;
				_slru = 0.66950;

				shared_ptr<t_gobj>  sat_obj = _gallobj->obj(_sat);
				shared_ptr<t_gpcv>  sat_pcv;
				if (sat_obj != 0)  sat_pcv = sat_obj->pcv(_time);
				int nsvn = str2int(sat_pcv->svcod());
				for (m = 0; m < 100; m++)
				{
					if (nsvn == G_svn_ilrs(m, 1))
					{
						_glog->comment(0, "GPS SVN to ILRS Number:"+to_string(nsvn)+"  "+_sat);
						break;
					}
					
				}

				if (m > 0)
				{
					_glog->comment(0, "GPS SVN to ILRS Number:not found" + to_string(nsvn) + "  " + _sat);
					islr = 0;
				}

			}
			else if (_sat[0] == 'R')
			{
				islr = 6;
				_slre = 0.13700;
				_slrn = 0.00300;
				_slru = 1.87370;  

				shared_ptr<t_gobj>  sat_obj = _gallobj->obj(_sat);
				shared_ptr<t_gpcv>  sat_pcv;
				if (sat_obj != 0)  sat_pcv = sat_obj->pcv(_time);
				int nsvn = str2int(sat_pcv->svcod());
				for (m = 0; m < 100; m++)
				{
					if (nsvn == G_svn_ilrs(m, 1))
					{
						_glog->comment(0, "GLONASS SVN to ILRS Number:" + to_string(nsvn) + "  " + _sat);
						break;
					}

				}

				if (m > 0)
				{
					_glog->comment(0, "GLONASS SVN to ILRS Number:not found" + to_string(nsvn) + "  " + _sat);
					islr = 0;
				}
			}
			else if (_sat[0] == 'E')
			{
				islr = 7;
				_slre = 0.00;  //Ref.Montenbruck et al. 2015 ASR, IOV and FOC are different
				_slrn = -0.00;
				_slru = 0.00;
				if (_sat.compare("E11") == 0)
				{
					_slre = -1.0922;
					_slrn = 0.0340;
					_slru = 0.6206;
				}
				else if (_sat.compare("E12") == 0)
				{
					_slre = -1.0927;
					_slrn = 0.0338;
					_slru = 0.6226;
				}
				else if (_sat.compare("E19") == 0)
				{
					_slre = -1.0927;
					_slrn = 0.0346;
					_slru = 0.6212;
				}
				else if (_sat.compare("E20") == 0)
				{
					_slre = -1.0927;
					_slrn = 0.0340;
					_slru = 0.6225;
				}
				else if (_sat.compare("E18") == 0)
				{
					_slre = 1.0199;
					_slrn = 0.0140;
					_slru = 0.5586;
				}
				else if (_sat.compare("E14") == 0)
				{
					_slre = 1.0146;
					_slrn = 0.0149;
					_slru = 0.5582;
				}
				else if (_sat.compare("E26") == 0)
				{
					_slre = 0.9625;
					_slrn = 0.0183;
					_slru = 0.5593;
				}
				else if (_sat.compare("E22") == 0)
				{
					_slre = 0.9641;
					_slrn = 0.0183;
					_slru = 0.5593;
				}
				else if (_sat.compare("E24") == 0)
				{
					_slre = 0.9616;
					_slrn = 0.0176;
					_slru = 0.5550;
				}
				else if (_sat.compare("E30") == 0)
				{
					_slre = 0.9622;
					_slrn = 0.0180;
					_slru = 0.5552;
				}
				else if (_sat.compare("E07") == 0)
				{
					_slre = 0.9642;
					_slrn = 0.0180;
					_slru = 0.5552;
				}
				else if (_sat.compare("E08") == 0)
				{
					_slre = 0.9641;
					_slrn = 0.0171;
					_slru = 0.5552;
				}
				else if (_sat.compare("E09") == 0)
				{

					_slre = 0.9645;
					_slrn = 0.0179;
					_slru = 0.5555;
				}
				else if (_sat.compare("E01") == 0)
				{
					_slre = 0.9652;
					_slrn = 0.0179;
					_slru = 0.5552;
				}
				else if (_sat.compare("E02") == 0)
				{
					_slre = 0.9661;
					_slrn = 0.0173;
					_slru = 0.5552;
				}
				else if (_sat.compare("E03") == 0)
				{
					_slre = 0.9602;
					_slrn = 0.0176;
					_slru = 0.5552;
				}
				else if (_sat.compare("E04") == 0)
				{
					_slre = 0.9657;
					_slrn = 0.0172;
					_slru = 0.5552;
				}
				else if (_sat.compare("E05") == 0)
				{
					_slre = 0.9646;
					_slrn = 0.0181;
					_slru = 0.5552;
				}
				else if (_sat.compare("E21") == 0)
				{
					_slre = 0.9698;
					_slrn = 0.0162;
					_slru = 0.5551;
				}
				else if (_sat.compare("E25") == 0)
				{
					_slre = 0.9632;
					_slrn = 0.0173;
					_slru = 0.5552;
				}
				else if (_sat.compare("E27") == 0)
				{
					_slre = 0.9653;
					_slrn = 0.0165;
					_slru = 0.5552;
				}
				else if (_sat.compare("E31") == 0)
				{
					_slre = 0.9642;
					_slrn = 0.0168;
					_slru = 0.5552;
				}
				else if (_sat.compare("E36") == 0)
				{
					_slre = 0.9666;
					_slrn = 0.0176;
					_slru = 0.5552; 
				}
				else if (_sat.compare("E13") == 0)
				{
					_slre = 0.9656;
					_slrn = 0.0177;
					_slru = 0.5552;
				}
				else if (_sat.compare("E15") == 0)
				{
					_slre = 0.9647;
					_slrn = 0.0182;
					_slru = 0.5552;
				}
				else if (_sat.compare("E33") == 0)
				{
					_slre = 0.9648;
					_slrn = 0.0176;
					_slru = 0.5552;
				}
			}
			else if (_sat[0] == 'J')
			{
				islr = 8;
				if (_sat.compare("J01") == 0) {
					_slre = 1.1511;
					_slrn = 0.5484;
					_slru = 2.6793;
				}
				else if (_sat.compare("J02") == 0) {
					_slre = 0.9852;
					_slrn = 0.8627;
					_slru = 2.5993;
				}
				else if (_sat.compare("J03") == 0) {
					_slre = 0.9849;
					_slrn = 0.8622;
					_slru = 2.6033;
				}
				else if (_sat.compare("J07") == 0) {
					_slre =-1.0820;
					_slrn = 0.4598;
					_slru = 2.5893;
				}
			}				 				
			if (islr != 0)
			{
				slrmod(_slre, _slrn, _slru);
			}
										
		}
		else
		{

			if (_sat.compare("SWARM-A") == 0)
			{
				islr = 1;
				_slre = -0.4633;
				_slrn = 0.51979;
				_slru = 0.30281;
			}
			else if (_sat.compare("SWARM-B") == 0)
			{
				islr = 2;
				_slre = -0.46356;
				_slrn = 0.52097;
				_slru = 0.30239;
			}
			else if (_sat.compare("SWARM-C") == 0)
			{
				islr = 3;
				_slre = -0.46592;
				_slrn = 0.52009;
				_slru = 0.30289;
			}
			else if (_sat.compare("GRACE-A") == 0)
			{
				islr = 4;
				_slre = -0.6000;
				_slrn = -0.3275;
				_slru = 0.2178;
			}
			else if (_sat.compare("GRACE-B") == 0)
			{
				islr = 5;
				_slre = -0.6000;
				_slrn = -0.3275;
				_slru = 0.2178;
			}
			else if (_sat.compare("GRACE-C") == 0)
			{
				islr = 6;
				_slre = -0.6000;
				_slrn = -0.3275;
				_slru = 0.2208;
			}
			else if (_sat.compare("GRACE-D") == 0)
			{
				islr = 7;
				_slre = -0.6000;
				_slrn = -0.3275;
				_slru = 0.2208;
			}
			else if (_sat.compare("SENTINEL-3A") == 0)
			{
				islr = 8;
				_slre = -0.3630;
				_slrn = 0.4279;
				_slru = 0.7922;
			}
			else if (_sat.compare("SENTINEL-3B") == 0)
			{
				islr = 8;
				_slre = -0.3530;
				_slrn = 0.4359;
				_slru = 0.7892;
			}
			else if (_sat.compare("HAIYANG-2A") == 0)
			{
				islr = 9;
				_slre = -1.01840;
				_slrn = -0.21930;
				_slru = 0.98730;
			}
			else if (_sat.compare("TIANGONG-2") == 0)
			{
				islr = 10;
				_slre = 0.0195;
				_slrn = -0.0044;
				_slru = 1.6112;

			}
			else if (_sat.compare("JASON-2") == 0)
			{
				islr = 11;
				_slre = 0.2172;
				_slrn = 0.5981;
				_slru = 0.6829;

			}
			else if (_sat.compare("JASON-3") == 0)
			{
				islr = 12;
				_slre = 0.2172;
				_slrn = 0.5981;
				_slru = 0.6829;

			}


			slrmod(_slre, _slrn, _slru);

		}

		return true;
	}

	void t_gslr::result_reserve(t_gresult& result)
	{
		result = *_gresult;

	}
	void t_gslr::preresult_reserve(t_gallpreproc& preresult)
	{
		preresult = *_gpreOBS;

	}

	void t_gslr::setResult(t_gresult& result)
	{
		
		_gresult = &result;
	}

	void t_gslr::setPreObs(t_gallpreproc& preobs)
	{
		_gpreOBS = &preobs;
	}

	void t_gslr::setRB(t_gallrb* rb)
	{
		_gallrb = rb;
	}

	void t_gslr::setSTA(t_gallsta* sta)
	{
		_gallsta = sta;
	}

	void t_gslr::slrmod(double e,double n,double u)
	{
		int jd;
		double sod, sodr, sodb;
		double fjd;
		t_gtime now, rtime, btime, time;
		int ista;
		Matrix rot_l2f;
		Matrix rotmat;
		ColumnVector xyz, xtmp;
		t_gtriple xsun, vsun, xlun, vlun;
		t_gtriple xsat, vsat;
		ColumnVector x_f(0.0,3), v_f(0.0,3), x_j(0.0,3), v_j(0.0,3), r1(0.0, 3), dloudx(0.0, 3), u_r1(0.0, 3);
		double _xsat[3], _vsat[3];
		ColumnVector cxsat(0.0, 3), cvsat(0.0, 3);
		ColumnVector cxsat_km(0.0, 3), cvsat_km(0.0, 3);
		ColumnVector cxsun(0.0, 3), cxsun_km(0.0, 3);
		double dlra_f = 0.0, dlra = 0.0;
		double r1length = 0.0;
		double azim_leo, elev_leo;
		double azim, elev;
		double trpdel, ztdpart;
		double temp;
		double press_wvp;
		double fcul_ztd, fcul_zhd, fcul_zwd;
		double fak_map;
		double std;
		double ztd_dif;
		double reldel0, reldel1;
		double sitrad, satrad;
		Matrix xtmp_leo1;
		ColumnVector xtmp_leo(0.0, 3), xtmp_lrr(0.0, 3);

		xtmp_leo1.ReSize(3, 1);

		//t_gprecisemodel gmodel=new()

		xyz.ReSize(3, 1);
		xtmp.ReSize(3, 1);
		
		rotmat.ReSize(3, 3);

		xyz(1) = e;
		xyz(2) = n;
		xyz(3) = u;

		double r2d = 180.0 / 3.141592653589;

		int size = _gnpt->getSize(_sat);

		int size_sitslr = _gslrf->getSize();

		Matrix rot_scf2crsleo(0.0, 3, 3); //rotation matrix of SCF to CRS, just for leo!
		Matrix rot_scf2crs(0.0, 3, 3);
		Matrix rot_trs2crs(0.0, 3, 3);   // rotation matrix of TRS to CRS
		Matrix rot_du(0.0, 3, 3);
		// prepare the trs2crs and crs2scf and so on;
		_trs2crs_cord = new t_gtrs2crs(false, _gpoleut1);
		_scf2crsleo_cord = new t_gscf2crsleo(_sat);
		_scf2crs_cord = new t_gscf2crs();
		
		_gresult->clearmap();

		for (int j = 0; j < size_sitslr; j++)
		{
			for (int jj = 0; jj < 3; jj++)
			{
				_gslrf->getOnegslrf(j)._xsta[jj] = _gslrf->getOnegslrf(j)._xsta[jj] + _gslrf->getOnegslrf(j)._xsta[jj + 3] * ((_beg.dmjd(false) + _end.dmjd(false)) / 2.0 - _gslrf->getOnegslrf(j).now.dmjd(false)) / 365.25;
			}
			xyz2ell(_gslrf->getOnegslrf(j)._xsta, _gslrf->getOnegslrf(j)._blh, false);
			_gslrf->getOnegslrf(j)._blh[2] = _gslrf->getOnegslrf(j)._blh[2] / 1000;
			for (int jj = 0; jj < 3; jj++)
			{
				_gslrf->getOnegslrf(j)._xsta[jj] = _gslrf->getOnegslrf(j)._xsta[jj] / 1000;    //m->km
			}
			
		}

		map<string,map<string,double>> RB = _gallrb->getrb();
		//map<string correction>* STA_COR = _gallsta->getsta();

		double range_bias = 0.0;

		for (int i = 0; i < size; i++)
		{
			dlra = 0.0;
			dlra_f = 0.0;
			_gnpt->getAllNpt(_sat,i).setDelay((_gnpt->getAllNpt(_sat,i)).getTwoWayTime()*1e-12);

			jd = _gnpt->getAllNpt(_sat, i).now.mjd();
			
			now.tsys(t_gtime::UTC);
			now=_gnpt->getAllNpt(_sat, i).now;
			now.tsys(t_gtime::GPS);

			// corrected by yqyuan from UTC to GPST
			// sod = now.sod()+now.dsec();
			if (_ifLEO == "NAV")
			{
				sod = now.sod(false) + now.dsec(false) + now.leapsec();
			}
			else
			{
				sod = now.sod()+now.dsec();
			}

			sodr = _gnpt->getAllNpt(_sat, i).getdelay()/2.0 + sod;

			sodb = _gnpt->getAllNpt(_sat, i).getdelay() + sod;

			for ( ista = 0; ista < size_sitslr; ista++)
			{
				if (_gnpt->getAllNpt(_sat, i).getCode() == _gslrf->getCode(ista))
				{
					break;
				}
			}
			if (ista >= size_sitslr)
			{
				continue;
			}


			if (RB.empty())
			{
				range_bias = 0.0;
			}
			else
			{
				for (auto iterRB = RB.begin(); iterRB != RB.end(); iterRB++)
				{
					for (auto itersta = iterRB->second.begin(); itersta != iterRB->second.end(); itersta++)
					{
						if (itersta->first == _gnpt->getAllNpt(_sat, i).getName() && iterRB->first == _sat)
						{
							range_bias = RB[_sat][_gnpt->getAllNpt(_sat, i).getName()];
							break;
						}
					}
					
				}
				
			}

			int code = _gnpt->getAllNpt(_sat, i).getCode();


			rot_l2f.ReSize(3, 3);

			enu2trs(_gslrf->getOnegslrf(ista)._blh[0], _gslrf->getOnegslrf(ista)._blh[1], rot_l2f);

			fjd = jd + (sod + 32.1840 + 19.00) / 86400.0;

			_gde->get_pv(fjd, "SUN", xsun, vsun);
			_gde->get_pv(fjd, "MOON", xlun, vlun);

			for (int ij = 0; ij < 3; ij++)
			{
				cxsun(ij + 1) = xsun.crd(ij);
				
			}

			rtime.from_mjd(jd, (int)sodr, sodr-(int)sodr);

			//m
			_gorb->get_crs_sat_crd(rtime, _sat, _xsat[0], _xsat[1], _xsat[2]);
			xsat.set(_xsat);

			//m/s
			_gorb->get_crs_sat_vel(rtime, _sat, _vsat[0], _vsat[1], _vsat[2]);
			vsat.set(_vsat);

			if (_ifLEO == "NAV")
			{
				for (int ij = 0; ij < 3; ij++)
				{
					cxsat(ij + 1) = _xsat[ij];
					cvsat(ij + 1) = _vsat[ij];
				}

				t_gtime tdt_grtime = rtime;
				tdt_grtime.tsys(t_gtime::TT);
				//ChangeTimeSys(GPS2TDT, rtime, tdt_grtime);
				double tdt = tdt_grtime.dmjd();

				_trs2crs_cord->calcRotMat(tdt_grtime, true, true, true);
				rot_trs2crs = _trs2crs_cord->getRotMat();

				for (int mn = 0; mn < 3; mn++)
				{
					cxsat_km(mn + 1) = cxsat(mn + 1)*1e-3;
					cvsat_km(mn + 1) = cvsat(mn + 1)*1e-3;
					cxsun_km(mn + 1) = cxsun(mn + 1)*1e-3;
				}
				//scf2crs
				_scf2crs_cord->calcRotMat(cxsat_km,cvsat_km,cxsun_km,rot_scf2crs);

				for (int k = 0; k < 3; k++)
				{
					cxsat(k + 1) = cxsat(k + 1) + (xyz(1)*rot_scf2crs(k + 1, 1) + xyz(2)*rot_scf2crs(k + 1, 2) + xyz(3)*rot_scf2crs(k + 1, 3));
				}

			}
			else
			{
				for (int ij = 0; ij < 3; ij++)
				{
					cxsat(ij + 1) = _xsat[ij];
					cvsat(ij + 1) = _vsat[ij];
				}

				t_gtime tdt_grtime = rtime;
				tdt_grtime.tsys(t_gtime::TT);
				//ChangeTimeSys(GPS2TDT, rtime, tdt_grtime);
				double tdt = tdt_grtime.dmjd();

				_trs2crs_cord->calcRotMat(tdt_grtime, true, true, true);
				rot_trs2crs = _trs2crs_cord->getRotMat();

				//scf2crs
				_scf2crsleo_cord->calcRotMat(rtime, _gattitude, cxsat, cvsat, rot_trs2crs, rot_scf2crsleo);

				xtmp = rot_scf2crsleo * xyz;

				for (int m = 0; m < 3; m++)
				{
					cxsat(m + 1) = xsat.crd(m) + xtmp(m + 1);
				}
			}

			for (int m = 0; m < 3; m++)
			{
				double ttmp= _gslrf->getOnegslrf(ista)._xsta[m] + _gslrf->getOnegslrf(ista)._offset[m] * 0.001;
				x_f(m + 1) = _gslrf->getOnegslrf(ista)._xsta[m] + _gslrf->getOnegslrf(ista)._offset[m] * 0.001;
				x_f(m + 1) = x_f(m + 1)*1000.0;
				v_f(m + 1) = 0.0;
			}

			//add station corrections
			for (auto iterSTA = _gallsta->getsta().begin(); iterSTA != _gallsta->getsta().end(); iterSTA++)
			{
				if (iterSTA->first == _gnpt->getAllNpt(_sat, i).getName())
				{
					x_f(1) += iterSTA->second.dx;
					x_f(2) += iterSTA->second.dy;
					x_f(3) += iterSTA->second.dz;
					break;
				}
			}

			time.from_mjd(jd, (int)sod, (sod - (int)sod));

			_calculate_sta_pos(time, &x_f, &v_f, &x_j, &v_j);
			//// TRS2CRS 
			//t_gtime tdt_gtime;
			//ChangeTimeSys(GPS2TDT, time, tdt_gtime);
			////double tdt = tdt_gtime.dmjd(false);

			//_trs2crs_cord->calcRotMat(tdt_gtime, true, true, true);
			//rot_trs2crs = _trs2crs_cord->getRotMat();
			//rot_du = _trs2crs_cord->getMatDu();
			//rot_du = rot_du * OMGE_DOT_GAL;

			//_rotate_posvel(x_f, v_f, x_j, &v_j, rot_trs2crs, rot_du);

			//t_gtriple tx_f(x_f);

			//bool tide_valid = _apply_rec_tides(time, tx_f);
			//if (!tide_valid)
			//{
			//	_glog->comment(1, "apply tide failed ");
			//	
			//}

			////tx_f = tx_f + _crt_obj->eccxyz(time);
			//for (int m = 0; m < 3; m++)
			//{
			//	x_f(m + 1) = tx_f.crd(m);
			//}
			//
			//x_j = rot_trs2crs * x_f;
			//v_j = rot_trs2crs * v_f;

			/*for (int m = 0; m < 3; m++)
			{
				_gslrf->getOnegslrf(ista)._xsta[m] = x_f(m + 1)*1e-3;
			}
			xyz2ell(_gslrf->getOnegslrf(ista)._xsta, _gslrf->getOnegslrf(ista)._blh, false);*/

			for (int m = 0; m < 3; m++)
			{
				r1(m + 1) = cxsat(m + 1) - x_j(m + 1);
			}

			r1length = sqrt(r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3));

			/*for (int m = 0; m < 3; m++)
			{
				u_r1(m + 1) = r1(m + 1) / r1length;
			}*/
			/*Matrix r1t;
			r1t.ReSize(1, 3);
			r1t = r1.t();*/
			xtmp_leo1 = r1.t() * rot_scf2crsleo;

			for (int m = 0; m < 3; m++)
			{
				xtmp_leo(m + 1) = xtmp_leo1(1, m + 1);
			}

			for (int m = 0; m < 3; m++)
			{
				xtmp_leo(m+1) = -xtmp_leo(m+1);
			}

			if (_ifLEO == "LEO")
			{
				if (_sat.find("SENTINEL") != string::npos)
				{
					xtmp_leo(1) -= 1.13403;
					xtmp_leo(2) -= 0.637905;
					xtmp_leo(3) -= 0.80118;
					xtmp_lrr(1) = xtmp_leo(2);
					xtmp_lrr(2) = -xtmp_leo(1);
					xtmp_lrr(3) = xtmp_leo(3);
					azim_leo = atan2(xtmp_lrr(2), xtmp_lrr(1));
					elev_leo = atan2(xtmp_lrr(3),sqrt(xtmp_lrr(1)*xtmp_lrr(1) + xtmp_lrr(2)*xtmp_lrr(2)));
				}
				else
				{
					azim_leo = atan2(xtmp_leo(2), xtmp_leo(1));
					elev_leo = atan2(sqrt(xtmp_leo(1)*xtmp_leo(1) + xtmp_leo(2)*xtmp_leo(2)), xtmp_leo(3));
				}
				get_lracor(elev_leo, azim_leo, &dlra_f);
			}
			

			//_crs2acr->calcRotMat(cxsat, cvsat, rotmat);

			// TRS2CRS 
			/*_trs2crs_cord->calcRotMat(time, true, true, true);
			rot_trs2crs = _trs2crs_cord->getRotMat();*/

			//x_j = rot_trs2crs * x_f;

			for (int m = 0; m < 3; m++)
			{
				double ttmp = _gslrf->getOnegslrf(ista)._xsta[m] + _gslrf->getOnegslrf(ista)._offset[m] * 0.001;
				x_f(m + 1) = _gslrf->getOnegslrf(ista)._xsta[m] + _gslrf->getOnegslrf(ista)._offset[m] * 0.001;
				x_f(m + 1) = x_f(m + 1)*1000.0;
				v_f(m + 1) = 0.0;
			}

			for (auto iterSTA = _gallsta->getsta().begin(); iterSTA != _gallsta->getsta().end(); iterSTA++)
			{
				if (iterSTA->first == _gnpt->getAllNpt(_sat, i).getName())
				{
					x_f(1) += iterSTA->second.dx;
					x_f(2) += iterSTA->second.dy;
					x_f(3) += iterSTA->second.dz;
					break;
				}
			}

			btime.from_mjd(jd, (int)sodb, sodb - (int)sodb);

			_calculate_sta_pos(btime, &x_f, &v_f, &x_j, &v_j);
			// TRS2CRS 
		 //   tdt_gtime;
			//ChangeTimeSys(GPS2TDT, btime, tdt_gtime);
			////double tdt = tdt_gtime.dmjd(false);

			//_trs2crs_cord->calcRotMat(tdt_gtime, true, true, true);
			//rot_trs2crs = _trs2crs_cord->getRotMat();
			//rot_du = _trs2crs_cord->getMatDu();
			//rot_du = rot_du * OMGE_DOT_GAL;

			//_rotate_posvel(x_f, v_f, x_j, &v_j, rot_trs2crs, rot_du);

			//t_gtriple tx_f1(x_f);

			//tide_valid = _apply_rec_tides(btime, tx_f1);
			//if (!tide_valid)
			//{
			//	_glog->comment(1, "apply tide failed ");

			//}

			////tx_f = tx_f + _crt_obj->eccxyz(time);
			//for (int m = 0; m < 3; m++)
			//{
			//	x_f(m + 1) = tx_f1.crd(m);
			//}

			//x_j = rot_trs2crs * x_f;
			//v_j = rot_trs2crs * v_f;

			/*for (int m = 0; m < 3; m++)
			{
				_gslrf->getOnegslrf(ista)._xsta[m] = x_f(m + 1)*1e-3;
			}
			xyz2ell(_gslrf->getOnegslrf(ista)._xsta, _gslrf->getOnegslrf(ista)._blh, false);*/

			for (int m = 0; m < 3; m++)
			{
				r1(m + 1) = cxsat(m + 1) - x_j(m + 1);
			}

			for (int m = 0; m < 3; m++)
			{
				dloudx(m + 1) = r1(m + 1) / sqrt(r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3));
			}

			r1length = r1length+sqrt(r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3));

			xtmp_leo1 = r1.t() * rot_scf2crsleo;

			for (int m = 0; m < 3; m++)
			{
				xtmp_leo(m + 1) = xtmp_leo1(1, m + 1);
			}
			for (int m = 0; m < 3; m++)
			{
				xtmp_leo(m + 1) = -xtmp_leo(m + 1);
			}
			
			if (_ifLEO == "LEO")
			{
				if (_sat.find("SENTINEL") != string::npos)
				{
					xtmp_leo(1) -= 1.13403;
					xtmp_leo(2) -= 0.637905;
					xtmp_leo(3) -= 0.80118;
					xtmp_lrr(1) = xtmp_leo(2);
					xtmp_lrr(2) = -xtmp_leo(1);
					xtmp_lrr(3) = xtmp_leo(3);
					azim_leo = atan2(xtmp_lrr(2), xtmp_lrr(1));
					elev_leo = atan2(xtmp_lrr(3),sqrt(xtmp_lrr(1)*xtmp_lrr(1) + xtmp_lrr(2)*xtmp_lrr(2)));
				}
				else
				{
					azim_leo = atan2(xtmp_leo(2), xtmp_leo(1));
					elev_leo = atan2(sqrt(xtmp_leo(1)*xtmp_leo(1) + xtmp_leo(2)*xtmp_leo(2)), xtmp_leo(3));
				}
				get_lracor(elev_leo, azim_leo, &dlra);
			}
			

			Matrix xtmp1;
			xtmp1.ReSize(1, 3);

			xtmp1 = r1.t()*rot_trs2crs;

			for (int m = 0; m < 3; m++)
			{
				xtmp(m + 1) = xtmp1(1, m + 1);
			}

			xtmp1 = xtmp.t()*rot_l2f;

			for (int m = 0; m < 3; m++)
			{
				xtmp(m + 1) = xtmp1(1, m + 1);
			}

			azim = atan2(xtmp(2), xtmp(1));
			elev = atan2(xtmp(3), sqrt(xtmp(1)*xtmp(1) + xtmp(2)*xtmp(2)));

			_gnpt->getAllNpt(_sat, i).setElev(elev*R2D);

			troposphere_delay_slr(elev, &trpdel, &ztdpart, i, ista);

			temp = _gnpt->getAllNpt(_sat, i).getDryTem() - 273.150;
			press_wvp = _gnpt->getAllNpt(_sat, i).getAtmHumidity()*exp(2.30260*(0.78580 + 7.50*temp / (237.30 + temp)));

			trop_slr_fculzd_hpa(i, ista, elev, &fcul_ztd, &fcul_zhd, &fcul_zwd, press_wvp);

			fak_map = trop_slr_fcul_a(elev, i, ista);

			std = fak_map * fcul_ztd;
			ztd_dif = (trpdel - std * 1e-3)*1e3;
			trpdel = std * 1.e-3;

			reldel0 = -2.0*(x_j(1)*v_j(1) + x_j(2)*v_j(2) + x_j(3)*v_j(3)) / (CLIGHT*1e-3);
			reldel0 = reldel0 + 2.0*(cxsat(1)*1e-3*cvsat(1)*1e-3 + cxsat(2)*1e-3*cvsat(2)*1e-3 + cxsat(3)*1e-3*cvsat(3)*1e-3) / (CLIGHT*1e-3);

			sitrad = sqrt(x_j(1)*x_j(1) + x_j(2)*x_j(2) + x_j(3)*x_j(3));
			satrad = sqrt(cxsat(1)*1e-3*cxsat(1)*1e-3 + cxsat(2)*1e-3*cxsat(2)*1e-3 + cxsat(3)*1e-3*cxsat(3)*1e-3);

			reldel1 = 2.0*EARTH_GM / ((CLIGHT*1e-3)*(CLIGHT*1e-3));
			reldel1 = log((sitrad + satrad + (r1length * 1e-3) / 2.0) / (sitrad + satrad - (r1length * 1e-3) / 2.0));
			reldel1 = 2.0*EARTH_GM/((CLIGHT*1e-3)*(CLIGHT*1e-3)) * log((sitrad + satrad + (r1length * 1e-3) / 2.0) / (sitrad + satrad - (r1length * 1e-3) / 2.0));

			double omc = (_gnpt->getAllNpt(_sat, i).getdelay()*(CLIGHT*1e-3) - (r1length*1e-3 + 2.0*trpdel + 2.0*reldel1))*1e3 / 2.0;
			omc = ((dlra / 1e3 + dlra_f / 1e3) / 2.0);
			_gnpt->getAllNpt(_sat, i).setOMC((_gnpt->getAllNpt(_sat, i).getdelay()*(CLIGHT*1e-3) - (r1length*1e-3 + 2.0*trpdel + 2.0*reldel1))*1e3 / 2.0 -range_bias+ ((dlra / 1e3 + dlra_f / 1e3) / 2.0));

			_result->setTime(_gnpt->getAllNpt(_sat, i).now);
			_result->setCode(code);
			_result->setName(_gnpt->getAllNpt(_sat, i).getName());
			
			_result->setOMC(_gnpt->getAllNpt(_sat, i).getOMC());
			_result->setElev( _gnpt->getAllNpt(_sat, i).getElev());
			

			if (_ifLEO == "LEO")
			{
				_result->setDLRA(dlra / 1000.0);
				_result->setNameL(_sat);
			}
			else
			{
				_result->setNameL(_sat);
			}

			_gresult->setOneResult(i, _result);
			_gresult->setOnesat(_sat, *(_gresult->getMap()));

			if (abs(_gnpt->getAllNpt(_sat, i).getOMC()) < 0.2 && _gnpt->getAllNpt(_sat, i).getName()!="AREL")
			{
				_gpreproc->setAtmHumidity(_gnpt->getAllNpt(_sat, i).getAtmHumidity());
				_gpreproc->setAtmPressure(_gnpt->getAllNpt(_sat, i).getAtmPressure());
				_gpreproc->setCode(_gnpt->getAllNpt(_sat, i).getCode());
				_gpreproc->setDryTem(_gnpt->getAllNpt(_sat, i).getDryTem());
				_gpreproc->setTwowayTime(_gnpt->getAllNpt(_sat, i).getdelay());
				_gpreproc->setWave(_gnpt->getAllNpt(_sat, i).getWave());
				_gpreproc->setSta(_gnpt->getAllNpt(_sat, i).getName());
				_gpreOBS->setSLROBS(_sat,now, *_gpreproc);
				//_gpreOBS->setSLROBS(now, _gnpt->getAllNpt(_sat, i).getdelay(), _gnpt->getAllNpt(_sat, i).getCode());//GPS    TWO-WAY TIME
			}
			else
			{
				int defg = 1;
			}
		}


	}

	void t_gslr::_update_rot_matrix(const t_gtime& epoch)
	{
		if (!_trs2crs_2000 && !_allerp)
		{
			if (_glog)
				_glog->comment(1, "t_gorbmodel::_update_rotmatrix ", "have no poleut1 data,cant' compute trs2crs!!");
			throw exception();
		}

		if (!_trs2crs_2000)
		{
			_trs2crs_2000 = new t_gtrs2crs(false, _allerp);
		}

		t_gtime tdt = epoch;
		tdt.tsys(t_gtime::TT);
		if (tdt != _trs2crs_2000->getCurtEpoch())
		{
			_trs2crs_2000->calcRotMat(tdt, true, true, true);
		}

#ifdef DEBUG_TRS
		//jdhuang : check for the trs2crs_2000
		cout << "Mark = trs2crs  : "
			<< setw(25) << tdt.mjd()
			<< setw(25) << tdt.sod()
			<< setw(25) << tdt.dsec() << endl
			<< setw(25) << fixed << setprecision(15)
			<< trs2crs_2000->rotmat() << endl
			<< trs2crs_2000->rotdx() << endl
			<< trs2crs_2000->rotdy() << endl
			<< trs2crs_2000->rotdu() << endl;
#endif // DEBUG
		return;

	}

	void t_gslr::_calculate_sta_pos(t_gtime time, ColumnVector* x_f, ColumnVector* v_f, ColumnVector* x_j, ColumnVector* v_j)
	{
		Matrix rot_scf2crsleo(0.0, 3, 3); //rotation matrix of SCF to CRS, just for leo!
		Matrix rot_trs2crs(0.0, 3, 3);   // rotation matrix of TRS to CRS
		Matrix rot_du(0.0, 3, 3);
		// TRS2CRS 
		t_gtime tdt_gtime = time;
		tdt_gtime.tsys(t_gtime::TT);
		//ChangeTimeSys(GPS2TDT, time, tdt_gtime);
		//double tdt = tdt_gtime.dmjd();

		_trs2crs_cord->calcRotMat(tdt_gtime, true, true, true);
		rot_trs2crs = _trs2crs_cord->getRotMat();
		rot_du = _trs2crs_cord->getMatDu();
		rot_du = rot_du * OMGE_DOT_GAL;

		_rotate_posvel(*x_f, *v_f, *x_j, v_j, rot_trs2crs, rot_du);

		t_gtriple tx_f(*x_f);

		bool tide_valid = _apply_rec_tides(time, tx_f);
		if (!tide_valid)
		{
			_glog->comment(1, "apply tide failed ");

		}

		//tx_f = tx_f + _crt_obj->eccxyz(time);
		for (int m = 0; m < 3; m++)
		{
			(*x_f)(m + 1) = tx_f.crd(m);
		}

		(*x_j) = rot_trs2crs * (*x_f);
	}

	void t_gslr::_rotate_posvel(ColumnVector& x_f, ColumnVector& v_f, ColumnVector& x_j, ColumnVector* v_j, Matrix& rotmat, Matrix& drotmat)
	{
		ColumnVector xhelp(0.0,3);
		
		x_j = rotmat * x_f;
		*v_j = rotmat * v_f;
		xhelp = drotmat * x_f;
		*v_j = *v_j + xhelp;
	}

	bool t_gslr::_apply_rec_tides(const t_gtime& epoch, t_gtriple& rec)
	{
		if (!_tide)
		{
			_glog->comment(0, "ERROR", "The tide ptr is nullptr.");
			return false;
		}

		_update_rot_matrix(epoch);
		double xpole = _trs2crs_2000->getXpole();
		double ypole = _trs2crs_2000->getYpole();
		double gast = _trs2crs_2000->getGmst();
		Matrix rot_trs2crs = _trs2crs_2000->getRotMat();
		string _site = "";

		t_gtriple tide(0.0, 0.0, 0.0);
		try
		{
			// solid tide
			t_gtideIERS* solid_ptr = dynamic_cast<t_gtideIERS*>(_tide);
			if (!solid_ptr)
			{
				_glog->comment(0, "ERROR", "can not get the solid tide ptr.");
				return false;
			}
			t_gtriple solid_earth = solid_ptr->tide_solid(epoch, rec, rot_trs2crs, _gde);

			// ocean load
			t_gtriple load_ocean = _tide->load_ocean(epoch, _site, rec);

			//pole tide
			t_gtideIERS* pole_ptr = dynamic_cast<t_gtideIERS*>(_tide);
			if (!pole_ptr)
			{
				_glog->comment(0, "ERROR", "can not get the pole tide ptr.");
				return false;
			}
			t_gtriple tide_pole = pole_ptr->tide_pole_pod(epoch, xpole, ypole, rec);

			//atmosph
			t_gtriple load_atmosph = _tide->load_atmosph();

			// freq tide
			t_gtideIERS* freq_ptr = dynamic_cast<t_gtideIERS*>(_tide);
			if (!freq_ptr)
			{
				_glog->comment(0, "ERROR", "can not get the freq tide ptr.");
				return false;
			}
			t_gtriple tide_freq = freq_ptr->tide_freq(_site, rec, gast);

			// jdhuang  for test
			tide = solid_earth +
				load_ocean +
				tide_pole +
				load_atmosph +
				tide_freq;

#ifdef DEBUG_PRECISEMODEL
			cout << "Tide" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd()
				<< " " << epoch.sod()
				<< " " << tide[0]
				<< " " << tide[1]
				<< " " << tide[2] << endl;

			cout << "solid_earth" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd()
				<< " " << epoch.sod()
				<< " " << solid_earth[0]
				<< " " << solid_earth[1]
				<< " " << solid_earth[2] << endl;

			cout << "load_ocean" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd()
				<< " " << epoch.sod()
				<< " " << load_ocean[0]
				<< " " << load_ocean[1]
				<< " " << load_ocean[2] << endl;

			cout << "tide_pole" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd()
				<< " " << epoch.sod()
				<< " " << tide_pole[0]
				<< " " << tide_pole[1]
				<< " " << tide_pole[2] << endl;

			cout << "load_atmosph" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd()
				<< " " << epoch.sod()
				<< " " << load_atmosph[0]
				<< " " << load_atmosph[1]
				<< " " << load_atmosph[2] << endl;

			cout << "tide_freq" + _site
				<< fixed << setprecision(15) << setw(20)
				<< " " << epoch.mjd()
				<< " " << epoch.sod()
				<< " " << tide_freq[0]
				<< " " << tide_freq[1]
				<< " " << tide_freq[2] << endl;
#endif // DEBUG

		}
		catch (...)
		{
			_glog->comment(0, "ERROR", "can not get the tide ptr.");
			return false;
		}

		//unit to m
		rec = rec + tide * 1.e3;
		return true;

	}


	double t_gslr::trop_slr_fcul_a(double elev,int i,int ista)
	{
		double EPSILON, SINE, T_C, COSPHI;
		double A10, A11, A12, A13;
		double A20, A21, A22, A23;
		double A30, A31, A32, A33;
		double A1, A2, A3;
		double MAP_ZEN;
		double trop_slr_FCUL_A;

		EPSILON = elev;
		SINE = sin(EPSILON);
		T_C = _gnpt->getAllNpt(_sat, i).getDryTem() - 273.1500;

		COSPHI = cos(_gslrf->getOnegslrf(ista)._blh[0]);

		A10 = 0.121008e-02;
		A11 = 0.17295e-05;
		A12 = 0.3191e-04;
		A13 = -0.18478e-07;
		
		A20 = 0.304965e-02;
		A21 = 0.2346e-05;
		A22 = -0.1035e-03;
		A23 = -0.1856e-07;

		A30 = 0.68777e-01;
		A31 = 0.1972e-04;
		A32 = -0.3458e-02;
		A33 = 0.1060e-06;

		A1 = A10 + A11 * T_C + A12 * COSPHI + A13 * _gslrf->getOnegslrf(ista)._blh[2] * 1000.0;
		A2 = A20 + A21 * T_C + A22 * COSPHI + A23 * _gslrf->getOnegslrf(ista)._blh[0] * 1000.0;
		A3 = A30 + A31 * T_C + A32 * COSPHI + A33 * _gslrf->getOnegslrf(ista)._blh[0] * 1000.0;

		MAP_ZEN = (1.000 + A1 / (1.000 + A2 / (1.000 + A3)));

		trop_slr_FCUL_A = MAP_ZEN / (SINE + A1 / (SINE + A2 / (SINE + A3)));

		return trop_slr_FCUL_A;

	}

	void t_gslr::trop_slr_fculzd_hpa(int i, int ista, double elev,double* FCUL_ZTD,double* FCUL_ZHD,double* FCUL_ZWD,double WVP)
	{
		double k0, k1, k2, k3, XC, CORR;
		double SIGMA, W0, W1, W2, W3, FH, FNH;
		double F;

		XC = 375.00;
		k0 = 238.01850;
		k1 = 19990.9750;
		k2 = 57.3620;
		k3 = 579.551740;

		W0 = 295.2350;
		W1 = 2.64220;
		W2 = -0.0323800;
		W3 = 0.0040280;

		SIGMA = 1000.0 / (_gnpt->getAllNpt(_sat, i).getWave()*1000.0);

		F = 1 - 0.00266*cos(2.0 * _gslrf->getOnegslrf(ista)._blh[0]) - 0.00028e-3 * _gslrf->getOnegslrf(ista)._blh[2] * 1000.0;

		CORR = 1.00 + 0.534e-6 * (XC - 450);

		FH = 0.010*CORR*((k1*(k0 + SIGMA * SIGMA)) / ((k0 - SIGMA * SIGMA)*(k0 - SIGMA * SIGMA)) + k3 * (k2 + SIGMA * SIGMA) / ((k2 - SIGMA * SIGMA)*(k2 - SIGMA * SIGMA)));

		*FCUL_ZHD = 2.416579e-3 * FH*_gnpt->getAllNpt(_sat, i).getAtmPressure() / F;

		FNH = 0.0031010*(W0 + 3.0*W1*SIGMA*SIGMA + 5.0*W2*SIGMA*SIGMA*SIGMA*SIGMA + 7.0*W3*SIGMA*SIGMA*SIGMA*SIGMA*SIGMA*SIGMA);

		*FCUL_ZWD = 1e-4 * (5.316e0*FNH - 3.759*FH)*WVP / F;

		*FCUL_ZTD = *FCUL_ZHD + *FCUL_ZWD;
	}

	void t_gslr::troposphere_delay_slr(double elev, double* delay, double* partial, int i,int ista)
	{
		double T0, fw, es, fpos, fwvl, e0, a, k, b;

		T0 = _gnpt->getAllNpt(_sat, i).getDryTem();
		
		fw = 1.000620 + 3.14e-6 * _gnpt->getAllNpt(_sat, i).getAtmPressure() + 5.6e-7 * (T0 - 273.15)*(T0 - 273.15);
		es = 0.01*exp(1.2378847e-5 * T0*T0 - 1.9121316e-2 * T0 + 33.937110470 - 6.3431645e3 / T0);
		fpos = 1.0 - 0.0026*cos(2.0*_gslrf->getOnegslrf(ista)._blh[0]) - 0.00031*_gslrf->getOnegslrf(ista)._blh[2];
		fwvl = 0.9650 + 0.0164 / _gnpt->getAllNpt(_sat, i).getWave() / _gnpt->getAllNpt(_sat, i).getWave() + 0.000228*pow((_gnpt->getAllNpt(_sat, i).getWave()), -4);
		e0 = _gnpt->getAllNpt(_sat, i).getAtmHumidity() * es*fw * 1e-2;
		a = 0.002357* _gnpt->getAllNpt(_sat, i).getAtmPressure() + 0.000141*e0;
		k = 1.163 - 0.00968*cos(2 * _gslrf->getOnegslrf(ista)._blh[0]) - 0.00104*T0 + 0.00001435*_gnpt->getAllNpt(_sat, i).getAtmPressure();
			
		b = 1.084e-8 * _gnpt->getAllNpt(_sat, i).getAtmPressure()*T0*k + 4.734e-8 * _gnpt->getAllNpt(_sat, i).getAtmPressure()*_gnpt->getAllNpt(_sat, i).getAtmPressure() / T0 * 2.0 / (3.0 - 1.0 / k);
		*partial = 1.0 / (sin(elev) + b / (a + b) / (sin(elev) + 0.010));
		*delay = fwvl / fpos * (a + b)*(*partial) * 1e-3;

	}



	void t_gslr::get_lracor(double zeni0,double azim0,double* dlra)
	{
		if (!_glra)
		{
			*dlra = 0.0;
			return;
		}
		double zeni;
		double azim;
		int ntmp;
		int iazi, izen;
		double x1, x2;
		double alpha;

		zeni = zeni0 * R2D;

		if (zeni > _glra->zen2(_sat))
			zeni = _glra->zen2(_sat) - 1e-5;
		if (zeni < _glra->zen1(_sat))
			zeni = _glra->zen1(_sat);

		if (_glra->dzen(_sat) == 0.0)
		{
			_glog->comment(0, "slr", "ERROR : slr ==> dzen in lra is zero");
			return;
		}

		azim = azim0 * R2D;
		if (azim < 0.0)
			azim = azim + 360.0;
		if (_sat.find("SENTINEL") != string::npos)
		{
			if (azim > 60.0)
			{
				ntmp = int(azim / 60.0);
				azim = azim - 60.0*ntmp;
			}
		}
		else
		{
			if (azim > 90.0)
			{
				ntmp = int(azim / 90.0);
				azim = azim - 90.0*ntmp;
			}
			if (azim > 45.0)
				azim = 90.0 - azim;
		}

		iazi = 0;

		if (_glra->dazi(_sat) != 0)
			iazi = int(azim / _glra->dazi(_sat) - 1e-5);

		if (_glra->dzen(_sat) != 0)
			izen = int((zeni - _glra->zen1(_sat)) / _glra->dzen(_sat));

		x1 = _glra->getlra(_sat,izen, iazi);
		x2 = _glra->getlra(_sat,izen + 1, iazi);

		if (iazi != 0)
		{
			alpha = (azim - _glra->azi1(_sat)) / _glra->dazi(_sat) - iazi;
			x1 = x1 + (_glra->getlra(_sat, izen, iazi + 1) - _glra->getlra(_sat, izen, iazi))*alpha;
			x2 = x2 + (_glra->getlra(_sat, izen + 1, iazi + 1) - _glra->getlra(_sat, izen + 1, iazi))*alpha;
		}

		alpha = (zeni - _glra->zen1(_sat)) / _glra->dzen(_sat) - izen ;
		*dlra = *dlra + x1 + (x2 - x1)*alpha;

		if (_sat.find("SENTINEL") != string::npos)
		{
			*dlra = *dlra*(-1);
		}
	}


	void t_gslr::enu2trs(double lat, double lon, Matrix& rotmat)
	{
		double coslat, sinlat, coslon, sinlon;

		coslat = cos(lat - G_PI / 2.0);
		sinlat = sin(lat - G_PI / 2.0);

		coslon = cos(-G_PI / 2.0 - lon);
		sinlon = sin(-G_PI / 2.0 - lon);

		rotmat(1, 1) = coslon;
		rotmat(1, 2) = sinlon * coslat;
		rotmat(1, 3) = sinlon * sinlat;
		rotmat(2, 1) = -sinlon;
		rotmat(2, 2) = coslon * coslat;
		rotmat(2, 3) = coslon * sinlat;
		rotmat(3, 1) = 0.0;
		rotmat(3, 2) = -sinlat;
		rotmat(3, 3) = coslat;

	}

	bool t_gslr::_init_data()
	{
		if ( !_gnpt || !_gslrf || !_gde || !_gpoleut1)
		{
			_glog->comment(0, "slr", "ERROR : slr ==> t_gslr::no comeplete file !");
			return false;
		}
		// get begin/end time,sampling, sat list and trans mode from xml
		if (!_get_setting())
		{
			_glog->comment(0, "slr", "ERROR : slr ==> t_gslr::there is something wrong witn xml !");
			return false;
		}

		// init lsq
		/*if (_lsqEstimator) {
			delete _lsqEstimator;
		}
		_lsqEstimator = new t_glsq(NULL);*/

		//_init_par();

		_set_out();

		return true;
	}

	bool t_gslr::_get_setting()
	{
		_beg = dynamic_cast<t_gsetgen*>(_set)->beg();
		_end = dynamic_cast<t_gsetgen*>(_set)->end();
		_sampling = dynamic_cast<t_gsetgen*>(_set)->sampling();

		//_end = _beg + ((int)((_end - _beg) / _sampling))*_sampling;

		// get satlist from xml
		//_sat = dynamic_cast<t_gsetorbdifleo*>(_set)->leo();

		//if (_sat.empty())
		//{
		//	_glog->comment(0, "orbdifleo", "ERROR : orbdif ==> No leos, check your xml file !");
		//	return false;
		//}

		//// get trans from xml
		//_trans = dynamic_cast<t_gsetorbdifleo*>(_set)->trans();

		return true;
	}

	void t_gslr::_set_out()
	{
		string tmp = dynamic_cast<t_gsetout*>(_set)->outputs("slr");
		if (tmp.empty())
		{
			tmp = "slr" + int2str(_beg.year()) + int2str(_beg.doy(), 3);
		}
		else
		{
			substitute(tmp, "$(date)", int2str(_beg.year()) + int2str(_beg.doy(), 3), false);
		}

		_slr = new t_giof;
		_slr->mask(tmp);
		_slr->append(false);

	}

}