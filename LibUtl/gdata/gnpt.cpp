/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-10-27	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gnpt.cpp
 * @brief		The base class used to store npt file information.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-10-27
 */

#include "gnpt.h"

using namespace std;

namespace great
{
	/** default constructor */
	t_gnpt::t_gnpt()
	{
		now.tsys(t_gtime::UTC);
		id_type(t_gdata::NPT);
	}

	/** default destructor */
	t_gnpt::~t_gnpt()
	{

	}

	void t_gnpt::setName(string name)
	{
		_name = name;
	}

	void t_gnpt::setID(int code)
	{
		_code = code;
	}

	void t_gnpt::setTimeScale(int scale)
	{
		_flag_time_scale = scale;
	}

	void t_gnpt::setLEOid(int leo)
	{
		_ileo = leo;
	}

	void t_gnpt::setTrop(int trop)
	{
		_Flag_trop_cor = trop;
	}

	void t_gnpt :: setMass(int mass)
	{
		_Flag_mass_cor = mass;
	}

	void t_gnpt::setWaveLength(double wave)
	{
		_wave_length = wave;
	}

	void t_gnpt::setSystemDelay(double delay)
	{
		_System_delay = delay;
	}

	void t_gnpt::setSigmaCalibration(double sigma)
	{
		_Sigma_Calibration = sigma;
	}

	void t_gnpt::setdsec(double sec)
	{
		_dsec = sec;
	}

	void t_gnpt::setTwoWayTime(double time)
	{
		_two_way_time = time;
	}

	void t_gnpt::setSigmaObs(double obs)
	{
		_Sigma_obs = obs;
	}

	void t_gnpt::setSigma(double Sigma)
	{
		_sigma = Sigma;
	}

	void t_gnpt::setAtmPressure(double pressure)
	{
		_Atm_pressure = pressure;
	}

	void t_gnpt::setDryTem(double degree)
	{
		_Dry_temperature = degree;
	}

	void t_gnpt::setAtmHumidity(double humidity)
	{
		_Atm_humidity = humidity;
	}


	void t_gnpt::setCode(int code)
	{
		_Code_n_point_window = code;
	}

	double t_gnpt::getdsec()
	{
		return _dsec;
	}

	double t_gnpt::getTwoWayTime()
	{
		return _two_way_time;
	}

	void t_gnpt::setOMC(double omc)
	{
		_omc = omc;
	}

	void t_gnpt::setDelay(double delay)
	{
		_delay = delay;
	}

	double t_gnpt::getdelay()
	{
		return _delay;
	}

	int t_gnpt::getCode()
	{
		return _code;
	}

	void t_gnpt::setElev(double elev)
	{
		_elev = elev;
	}

	double t_gnpt::getElev()
	{
		return _elev;
	}

	double t_gnpt::getDryTem()
	{
		return _Dry_temperature;
	}

	double t_gnpt::getAtmPressure()
	{
		return _Atm_pressure;
	}

	double t_gnpt::getWave()
	{
		return _wave_length;
	}

	double t_gnpt::getAtmHumidity()
	{
		return _Atm_humidity;
	}

	double t_gnpt::getOMC()
	{
		return _omc;
	}

	string t_gnpt::getName()
	{
		return _name;
	}
}