/**
 *
 * @verbatim
	 History
	  -1.0 Yujie Qin   2019-01-06 Creat the file.
	  -1.1 Yujie Qin   2019-03-05 edit for final version.
	  -1.2 Yujie Qin   2019-03-15 Code comments are added in Doxygen format.
	  -1.3 Jiande Huang 2019-04-02 Adding abnormal throwing mechanism of functions, but need check by someone else focusing on LEO 
   @endverbatim
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		     gfcustomer.cpp
 * @brief		     The class for customer model, get the acceleration caused by customer.
 *
 * @author           Yujie Qin, Wuhan University
 * @version		     1.0.0
 * @date		     2019-01-06
 *
 */

#include "gforce/gfcustomer.h"
#include "gutils/ginfolog.h"
#include "newmat/newmatio.h"

using namespace std;
namespace great
{
	/** init some variables based on the default information.
	*/
	t_gfcustomer::t_gfcustomer()
	{
		for (int i = 0; i < 6; i++)
		{
			_mode[i] = 0.0;
		}
		_period_time = 1;
		_isPartial = false;
		_same_sat = true;
		_cus_par =
		{
			ICS_PAR::EMP_Ac, ICS_PAR::EMP_Bc, ICS_PAR::EMP_Sc1,ICS_PAR::EMP_Cc1, ICS_PAR::EMP_Sc2, ICS_PAR::EMP_Cc2,
			ICS_PAR::EMP_Aa, ICS_PAR::EMP_Ba, ICS_PAR::EMP_Sa1,ICS_PAR::EMP_Ca1, ICS_PAR::EMP_Sa2, ICS_PAR::EMP_Ca2,
			ICS_PAR::EMP_Ar, ICS_PAR::EMP_Br, ICS_PAR::EMP_Sr1,ICS_PAR::EMP_Cr1, ICS_PAR::EMP_Sr2, ICS_PAR::EMP_Cr2
		};
	}

	/** constructor. */
	t_gfcustomer::t_gfcustomer(t_init_arg* arg, t_glog* log) : t_gforce(log)
	{
		string function_name = "t_gfcustomer::t_gfcustomer";

		t_init_customer *cus = NULL;
		cus = dynamic_cast<t_init_customer*>(arg);
		if (!cus)
		{
			write_log(log, 4, function_name, "ERROR : init failed");
			write_log(log, 4, function_name, "ERROR : no init argument");
		}
		else
		{
			for (int i = 0; i < 6; i++)
			{
				_mode[i] = 0.0;
			}

			_period_time = cus->period;
			_isPartial = cus->ispartial;
			_mass = cus->mass;
			_same_sat = true;
			_cus_par =
			{
				ICS_PAR::EMP_Ac,ICS_PAR::EMP_Bc, ICS_PAR::EMP_Sc1, ICS_PAR::EMP_Cc1, ICS_PAR::EMP_Sc2, ICS_PAR::EMP_Cc2,
				ICS_PAR::EMP_Aa,ICS_PAR::EMP_Ba, ICS_PAR::EMP_Sa1, ICS_PAR::EMP_Ca1, ICS_PAR::EMP_Sa2, ICS_PAR::EMP_Ca2,
				ICS_PAR::EMP_Ar,ICS_PAR::EMP_Br, ICS_PAR::EMP_Sr1, ICS_PAR::EMP_Cr1, ICS_PAR::EMP_Sr2, ICS_PAR::EMP_Cr2
			};
		}
	}

	/** destructor. */
	t_gfcustomer::~t_gfcustomer()
	{

	}

	/** calculate the accleration and the corresponding partial matrix.
	*@param [in] argument Parameters to be used in calculation.
	*/
	bool t_gfcustomer::calculate(t_gforce_calarg* argument)
	{
		t_gfcalarg_customer* arg = NULL;
		arg = dynamic_cast<t_gfcalarg_customer*>(argument);
		if (!arg || arg->scf2crs->is_zero() || arg->scf2crs->nrows() != 3 || arg->scf2crs->ncols() != 3)
		{
			cout << "change ERROR" << endl;
			return false;
		}

		try
		{
			int npara = (*(arg->cust_pars)).size();
			_acc = 0.0;
			_partial.ReSize(npara, 3);

			int dinx;
			ICS_PAR parname;
			_customer_pars = *(arg->cust_pars);
			// For one satellite, the index in 18 customer parameters and the corresponding index in the parameter list will not change
			// So it's not necessary to judge the index every time, just for the first time

			//static map<int, int> parindex;
			//static bool isfirst = true;
			map<ICS_PAR, double>::iterator it;
			if (_same_sat)
			{
				_same_sat = false;
				for (int i = 0; i < 18; i++)
				{
					parname = _cus_par[i];
					dinx = 0;
					for (it = _customer_pars.begin(); it != _customer_pars.end(); it++)
					{
						dinx++;
						if (it->first == parname)
						{
							_par_index[i] = dinx;
							break;
						}
					}
				}
			}
			set_mode(*(arg->period), *(arg->time));

			int iacr, iterm;
			map<int, int>::iterator it_ind;

			for (it_ind = _par_index.begin(); it_ind != _par_index.end(); it_ind++)
			{
				parname = _cus_par[it_ind->first];
				iacr = it_ind->first / 6 + 1;
				iterm = it_ind->first - (iacr - 1) * 6;

				_partial.Row(it_ind->second) << (*(arg->scf2crs)).Column(iacr).t() * _mode[iterm];

				if (_customer_pars[parname] == 0.0)
				{
					continue;
				}

				_acc += (*(arg->scf2crs)).Column(iacr) * _mode[iterm] * _customer_pars[parname];

			}

			return true;
		}
		catch (...)
		{
			write_log(_glog, 0, "ERROR : t_gfcustomer::calculate throw exception while calculating Density");
			return false;
		}
	}

	/** set mode for the customer model: constant, linear, one cpr and two cpr.
	*@param [in] period the period of sat.
	*@param [in] dt the second of dynamic time.
	*/
	void t_gfcustomer::set_mode(double period, double dt)
	{
		try
		{
			_mode[0] = 1.0 * 1e-9;
			_mode[1] = dt * 1e-9;
			_mode[2] = 1.0 * cos(period * dt) * 1e-9;
			_mode[3] = 1.0 * sin(period * dt) * 1e-9;
			_mode[4] = 1.0 * cos(period / 2.0 * dt) * 1e-9;
			_mode[5] = 1.0 * sin(period / 2.0 * dt) * 1e-9;
		}
		catch (...)
		{
			write_log(_glog, 0, "ERROR : t_gfcustomer::set_Mode throw exception while calculating Density");
		}
	}
}