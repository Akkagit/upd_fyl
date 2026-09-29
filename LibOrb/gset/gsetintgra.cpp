/**
*
* @verbatim
	History
	 -1.0 Jiande Huang	 2019-04-07 creat the file.
	 -1.1 Yujie Qin      2019-07-04 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gsetintgra.cpp
* @brief		number integrator  set for satellite orbit integral from XML
*
* @author       Jiande Huang, Wuhan University
* @version		1.0.0
* @date			2019-04-07
*
*/

#include "gset/gsetintgra.h"
using namespace gnut;
namespace great
{
	t_gsetintgra::t_gsetintgra() : t_gsetbase()
	{

	}

	t_gsetintgra::~t_gsetintgra()
	{

	}

	string t_gsetintgra::rkf_model()
	{
		//get rkf model node
		string rkf_order = _doc.child(XMLKEY_ROOT).child(XMLKEY_INTEGRATOR).child(XMLKEY_INTEGRATOR_RKF).attribute("order").value();

#ifdef DEBUG
		cout << "The rkf_model is : " << rkf_order + "-th rkf" << endl;
#endif // DEBUG

		if (rkf_order.empty())
		{
			_log->comment(1, "WARNNING : no setting for the rkf order, use the defalut number : " + int2str(_default_rkf_num));
			return "RKF" + int2str(_default_rkf_order - 1) + "(" + int2str(_default_rkf_order) + ")";
		}
		else
		{
			return "RKF" + int2str(str2int(rkf_order) - 1) + "(" + rkf_order + ")";
		}
	}

	string t_gsetintgra::adams_model()
	{
		//get rkf model node
		string adams_order = _doc.child(XMLKEY_ROOT).child(XMLKEY_INTEGRATOR).child(XMLKEY_INTEGRATOR_ADAMS).attribute("order").value();

#ifdef DEBUG
		cout << "The adams model is : " << adams_order + "-th adams" << endl;
#endif // DEBUG

		if (adams_order.empty())
		{
			_log->comment(1, "WARNNING : no setting for the adams order, use the defalut number : " + int2str(_default_rkf_num));
			return "Adams_" + int2str(_default_adams_order);
		}
		else
		{
			return "Adams_" + adams_order;
		}
	}

	int t_gsetintgra::rkf_num()
	{
		//get rkf model node
		string rkf_num = _doc.child(XMLKEY_ROOT).child(XMLKEY_INTEGRATOR).child(XMLKEY_INTEGRATOR_RKF).attribute("rkf_num").value();

#ifdef DEBUG
		cout << "The setting of rkf_num is : " + rkf_num << endl;
#endif // DEBUG

		if (rkf_num.empty())
		{
			_log->comment(1, "WARNNING : no setting for the rkf_num, use the defalut number : " + int2str(_default_rkf_num));
			return _default_rkf_num;
		}
		else
		{
			return str2int(rkf_num);
		}

	}

	int t_gsetintgra::rkf_order()
	{
		//get rkf model node
		string rkf_order = _doc.child(XMLKEY_ROOT).child(XMLKEY_INTEGRATOR).child(XMLKEY_INTEGRATOR_RKF).attribute("order").value();

#ifdef DEBUG
		cout << "The setting of rkf_order is : " + rkf_order << endl;
#endif // DEBUG

		if (rkf_order.empty())
		{
			_log->comment(1, "WARNNING : no setting for the rkf_num, use the defalut number : " + int2str(_default_rkf_order));
			return _default_rkf_order;
		}
		else
		{
			return str2int(rkf_order);
		}
	}

	int t_gsetintgra::adams_order()
	{
		//get rkf model node
		string adams_order = _doc.child(XMLKEY_ROOT).child(XMLKEY_INTEGRATOR).child(XMLKEY_INTEGRATOR_RKF).attribute("order").value();

#ifdef DEBUG
		cout << "The setting of adams_order is : " + adams_order << endl;
#endif // DEBUG

		if (adams_order.empty())
		{
			_log->comment(1, "WARNNING : no setting for the rkf_num, use the defalut number : " + int2str(_default_adams_order));
			return _default_adams_order;
		}
		else
		{
			return str2int(adams_order);
		}
	}

	double t_gsetintgra::adams_inp_step()
	{
		//get rkf model node
		string adams_inp_step = _doc.child(XMLKEY_ROOT).child(XMLKEY_INTEGRATOR).child(XMLKEY_INTEGRATOR_ADAMS).attribute("inp_step").value();

#ifdef DEBUG
		cout << "The setting of adams_inp_step is : " + adams_inp_step << endl;
#endif // DEBUG

		if (adams_inp_step.empty())
		{
			_log->comment(1, "WARNNING : no setting for the rkf_num, use the defalut number : " + dbl2str(_default_adams_inp_step));
			return _default_adams_inp_step;
		}
		else
		{
			return str2dbl(adams_inp_step);
		}
	}

	double t_gsetintgra::adams_out_step()
	{
		//get rkf model node
		string adams_out_step = _doc.child(XMLKEY_ROOT).child(XMLKEY_INTEGRATOR).child(XMLKEY_INTEGRATOR_ADAMS).attribute("out_step").value();

#ifdef DEBUG
		cout << "The setting of adams_out_step is : " + adams_out_step << endl;
#endif // DEBUG

		if (adams_out_step.empty())
		{
			_log->comment(1, "WARNNING : no setting for the rkf_num, use the defalut number : " + dbl2str(_default_adams_out_step));
			return _default_adams_out_step;
		}
		else
		{
			return str2dbl(adams_out_step);
		}
	}

	void t_gsetintgra::check()
	{
	}

	void t_gsetintgra::help()
	{
		_gmutex.lock();

		// jdhuang
		cerr << "<number_integrator>" << endl
			<< "  <rkf    order = \"7\"  rkf_num  = \"30\" / > " << endl
			<< "  <adams  order = \"11\" inp_step = \"60\" out_step = \"300\" / > " << endl
			<< "</number_integrator>" << endl;

		_gmutex.unlock();
		return;
	}
}
