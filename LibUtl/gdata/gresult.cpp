#include "gresult.h"
#include "gdata/gdata.h"

using namespace std;

namespace great
{

	/** default constructor */
	t_gresult_one::t_gresult_one()
	{
		//id_type(t_gdata::GRESULT);
	}

	/** default destructor */
	t_gresult_one::~t_gresult_one()
	{

	}

	void t_gresult_one::setTime(t_gtime& time)
	{
		_now = time;
	}

	void t_gresult_one::setCode(int code)
	{
		_code = code;
	}

	void t_gresult_one::setName(string name)
	{
		_name = name;
	}

	void t_gresult_one::setNameL(string name)
	{
		_name_long = name;
	}

	void t_gresult_one::setOMC( double omc)
	{
		_omc = omc;
	}

	void t_gresult_one::setElev( double elev)
	{
		_elev = elev;
	}

	void t_gresult_one::setDLRA(double dlra)
	{
		_dlra = dlra;
	}

	
}