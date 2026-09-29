/**
*
* @verbatim
History
-1.0 xjhan  2019-11-22  creat the file.

@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gpppleomodel.cpp
* @brief	precise model for ppp with leo
*
* @author   xjhan , Wuhan University
* @version	1.0.0
* @date		2019-11-22
*
*/

#include "cstdlib"
#include "gmodels/gpppleomodel.h"
#include "gall/gallobs.h"
#include "gutils/ginfolog.h"

//#define DEBUG_GLO

namespace gnut {

	t_gpppleomodel::t_gpppleomodel(t_gallproc* data, t_glog* log, t_gsetbase* setting)
		:t_gprecisemodel(data, log, setting)
	{
		
		//if (dynamic_cast<t_gallobs*>((*data)[t_gdata::ALLOBS]))
		//{
		//	_glofrq_num = dynamic_cast<t_gallobs*>((*data)[t_gdata::ALLOBS])->glo_freq_num();

		//	t_gtime beg_time = dynamic_cast<t_gsetgen*>(setting)->beg();

		//	raw_CODE_bias(beg_time, _gallbias, _glofrq_num, _simu_raw_dcb);
		//}
	}

	t_gpppleomodel::~t_gpppleomodel()
	{

	}

	double t_gpppleomodel::windUp(t_gsatdata& satdata, const ColumnVector& rRec)
	{
		gtrace("t_gpppleomodel::windUp");

		return 0.0;
	}

	void t_gpppleomodel::apply_GaussNoiseP(string sat, double& obs)
	{
		double mean = 0.0;
		double sigma = 1;//1m
		double epsilon = 2.2250738585072014e-308;
		double random1 = 0, random2 = 0;

		if (_sat_noise_generateP.find(sat) == _sat_noise_generateP.end())
		{
			_sat_noise_generateP[sat] = false;
		}

		_sat_noise_generateP[sat] = !(_sat_noise_generateP[sat]);

		if (!_sat_noise_generateP[sat])
		{
			obs = obs + _noiseP_next * sigma + mean;
			return;
		}

		while (random1 <= epsilon)
		{
			random1 = rand() % 100 / (float)(100);  // rand() % 100 --> 产生0~99这100个整数中的一个随机整数
			random2 = rand() % 100 / (float)(100);  // rand() % 100 / (float)(100) --> 产生0.00~1.00 中的随机小数（2位）
		}

		double noise = sqrt(-2.0 * log10(random1)) * cos(2.0 * G_PI * random2);
		_noiseP_next = sqrt(-2.0 * log10(random1)) * sin(2.0 * G_PI * random2);
		obs = obs + noise * sigma + mean;
	}

	void t_gpppleomodel::apply_GaussNoiseL(string sat, double& obs)
	{
		double mean = 0.0;
		double sigma = 0.01;//3mm
		double epsilon = 2.2250738585072014e-308;
		double random1 = 0, random2 = 0;

		if (_sat_noise_generateL.find(sat) == _sat_noise_generateL.end())
		{
			_sat_noise_generateL[sat] = false;
		}

		_sat_noise_generateL[sat] = !(_sat_noise_generateL[sat]);

		if (!_sat_noise_generateL[sat])
		{
			obs = obs + _noiseL_next * sigma + mean;
			return;
		}

		while (random1 <= epsilon)
		{
			random1 = rand() % 100 / (float)(100);
			random2 = rand() % 100 / (float)(100);
		}

		double noise = sqrt(-2.0 * log10(random1)) * cos(2.0 * G_PI * random2);
		_noiseL_next = sqrt(-2.0 * log10(random1)) * sin(2.0 * G_PI * random2);
		obs = obs + noise * sigma + mean;
	}

	bool raw_CODE_bias(const t_gtime epo, t_gallbias* allbias, map<string, int>& R_freq_ID, map<string, map<GOBS, double>>& raw_bia)
	{
		double l1, l2;
		double cor[5];
		set<GSYS> system = { GPS, GLO, LEO };
		t_map_sats all_sats = GNSS_SATS();

		for (auto sys : system)
		{
			set<string> sats = all_sats[sys];
			for (auto sat : sats)
			{
				if (sys == GPS || sys == LEO)
				{
					l1 = CLIGHT / double(G01_F);
					l2 = CLIGHT / double(G02_F);
				}
				else
				{
					l1 = CLIGHT / double(R01_F(R_freq_ID[sat]));
					l2 = CLIGHT / double(R02_F(R_freq_ID[sat]));
				}

				double c2 = l1 / l2;
				double c1 = 1.0 / (1.0 - c2 * c2);
				c2 = c2 * c2 * c1;

				double p1p2 = allbias->get(epo, sat, P1, P2, "COD_R");  //单天解没问题，多天会有问题
				double p1c1 = allbias->get(epo, sat, P1, C1, "COD_R");
				double p2c2 = allbias->get(epo, sat, P2, C2, "COD_R");

				cor[0] = c2 * p1p2;
				cor[1] = c1 * p1p2;
				cor[2] = cor[0] + p1c1;
				cor[3] = cor[1] + p2c2;
				cor[4] = cor[1] + p1c1;

				raw_bia[sat][C1C] = cor[2];
				raw_bia[sat][C1P] = cor[0];
				raw_bia[sat][C1Y] = cor[0];
				raw_bia[sat][C1W] = cor[0];
				raw_bia[sat][C2C] = cor[3];
				raw_bia[sat][C2P] = cor[1];
				raw_bia[sat][C2Y] = cor[1];
				raw_bia[sat][C2W] = cor[1];
				raw_bia[sat][C2X] = cor[4];
			}
		}

		return true;
	}

}