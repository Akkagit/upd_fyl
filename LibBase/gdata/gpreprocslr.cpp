#include "gpreprocslr.h"
#include "gdata/gdata.h"

using namespace std;

namespace great
{

	/** default constructor */
	t_gpreprocslr::t_gpreprocslr()
	{
		//id_type(t_gdata::GRESULT);
	}

	/** default destructor */
	t_gpreprocslr::~t_gpreprocslr()
	{

	}

	/*void t_gpreprocslr::setTime(t_gtime& time)
	{
		_now = time;
	}
*/
	void t_gpreprocslr::setCode(double code)
	{
		_code = code;
	}

	void t_gpreprocslr::setTwowayTime(double sec)
	{
		_two_way_time = sec;
	}

	void t_gpreprocslr::setAtmPressure(double pressure)
	{
		_atm_pressure = pressure;
	}

	void t_gpreprocslr::setDryTem(double drytem)
	{
		_dry_tem = drytem;
	}

	void t_gpreprocslr::setWave(double wave)
	{
		_wave = wave;
	}

	void t_gpreprocslr::setAtmHumidity(double humidity)
	{
		_atm_humidity = humidity;
	}

	void t_gpreprocslr::setSta(string name)
	{
		_sta_name = name;
	}

	void t_gpreprocslr::setSat(string sat)
	{
		_sat_name = sat;
	}

}