
#ifndef GNPT_H
#define GNPT_H

#include <iostream>
#include <string.h>
#include "gutils/ExportUtl.h"
#include "gutils/gmutex.h"
#include "gdata/gdata.h"

using namespace std;
using namespace gnut;

namespace great
{

	class UTL_LIBRARY_EXPORT  t_gnpt : public t_gdata
	{
	public:
		t_gnpt();
		virtual ~t_gnpt();

		void setName(string name);
		void setID(int code);
		void setTimeScale(int scale);
		void setLEOid(int leo);
		void setTrop(int trop);
		void setMass(int mass);
		void setWaveLength(double wave);
		void setSystemDelay(double delay);
		void setSigmaCalibration(double sigma);
		void setdsec(double sec);
		void setTwoWayTime(double time);
		void setSigmaObs(double obs);
		void setSigma(double Sigma);
		void setAtmPressure(double pressure);
		void setDryTem(double degree);
		void setAtmHumidity(double humidity);
		void setCode(int code);
		void setOMC(double omc);
		void setDelay(double delay);
		void setElev(double elev);

		double getdsec();
		double getdelay();
		double getTwoWayTime();
		int    getCode();
		double getElev();
		double getDryTem();
		double getAtmPressure();
		double getWave();
		double getAtmHumidity();
		double getOMC();
		string getName();

		t_gtime now;

	protected:
		string _name;            ///< slr station name
		int _code;               ///< alr station code
		int _flag_time_scale;    ///< Station Epoch Time Scale - indicates the time scale reference.
		int _ileo;            ///< ILRS Satellite Identifier (Based on COSPAR ID)
		int _Flag_trop_cor;      
		int _Flag_mass_cor;
		double _wave_length;
		double _System_delay;
		double _Sigma_Calibration;
		double _dsec;
		double _two_way_time;
		double _Sigma_obs;
		double _sigma;
		int _Code_n_point_window;
		double _Atm_pressure;
		double _Dry_temperature;
		double _Atm_humidity;
		double _omc = 0.0;
		double _delay=0.0;
		double _elev;
		
	};
}

#endif