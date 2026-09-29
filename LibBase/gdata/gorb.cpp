/**
 *
 * @verbatim
	 History
	  -1.0	Jie Li		2019-01-06	creat the file.
	  -1.1  Jie Li		2019-04-01	add code in doxygen style.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gorb.cpp
 * @brief		The base class used to storage the orb data.
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 */

#include "gdata/gorb.h"
#include "gmodels/gpoly.h"

using namespace std;

namespace great {

	/** @brief constructor. */
	t_gorb::t_gorb()
	{
		gtrace("t_orb::constructor");
		id_type(t_gdata::ORB);
	}

	/** @brief destructor. */
	t_gorb:: ~t_gorb()
	{

	}

	/**
	* @brief add orb head data.
	* @param[in]   orb_hd   orb head data
	*/
	void t_gorb::add_head(const t_orb_hd& orb_hd)
	{
		gtrace("t_orb::add_head");
		_gmutex.lock();

		_orb_hd = orb_hd;

		_gmutex.unlock();

		return;
	}

	/**
		* @brief add orb sat data(include sat list,sat num).
		* @param[in]   sys   system
		* @param[in]   sats  sat list
		*/
	void t_gorb::add_sat_info(set<string>& sys, set<string>& sats)
	{
		if (sys.size() > 1)
		{
			_orb_hd.sat_group = "GNS";
		}
		else
		{
			_orb_hd.sat_group = *(sys.begin());
		}

		//sat type
		set<string>::iterator itsys = sys.begin();
		for (itsys; itsys != sys.end(); ++itsys)
		{
			char sys_tmp = t_gsys::str2char(*itsys);
			_orb_hd.sat_type.insert(sys_tmp);
		}
		if (IS_SUPPORTED_LEO(*sats.begin()))
		{
			_orb_hd.sat_type.clear();
			_orb_hd.sat_type.insert('L');
			_orb_hd.sat_group = "LEO";
		}

		//sat num
		_orb_hd.sat_num = sats.size();
		//sat list
		_orb_hd.sat_prn = sats;

	}

	/**
	* @brief add orb funct data to map.
	* @param[in]   epo			  epoch(=(epoch-ref)/out_step)
	* @param[in]   sat			  satellite prn
	* @param[in]   one_funct_rec  one record(one line) of funct of a sat in epo
	*/
	void t_gorb::add_funct(int epo, string sat, t_gtriple& one_funct_rec)
	{
		gtrace("t_orb::add_funct");
		_gmutex.lock();

		_funct[sat][epo].push_back(one_funct_rec);

		_gmutex.unlock();

		return;
	}

	/**
	* @brief add orb funct data to map.
	* @param[in]   epo		epoch(=(epoch-ref)/out_step)
	* @param[in]   sat		satellite prn
	* @param[in]   funct	all funct of one sat in epo
	*/
	void t_gorb::add_funct(int epo, string sat, const vector<double>& funct)
	{
		//gtrace("t_orb::add_funct");
		//_gmutex.lock();
		int num = funct.size();
		int npar = num / 3;
		t_sat_funct tmp_sat;
		for (int i = 0; i < npar; i++)
		{
			int index = i * 3;
			t_gtriple tmp(funct[index], funct[index + 1], funct[index + 2]);
			tmp_sat.push_back(tmp);
		}
		_funct[sat][epo] = tmp_sat;
		//_gmutex.unlock();
	}

	/**
	* @brief add orb funct data to map.
	* @param[in]   sat		satellite prn
	* @param[in]   funct	all funct of one sat
	*/
	void t_gorb::add_funct(string sat, map<int, vector< t_gtriple> >& funct)
	{
		_funct[sat] = funct;
	}

	/**
	* @brief add parlist for one sat to map.
	* @param[in]   sat		satellite prn
	* @param[in]   parlist	parlistt of one sat
	*/
	void t_gorb::add_parlist(string sat, vector<string> parlist)
	{
		_par_list[sat] = parlist;
	}

	/** @brief whether the orb data is empty.
	* @return  bool
	*	@retval   true   orb data is empty
	*   @retval   false  orb data is existent
	*/
	bool t_gorb::is_empty()
	{
		if (_funct.size() == 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	/**
	* @brief get orb head data.
	* @param[out]   orb_hd	orb head data
	*/
	void t_gorb::get_head(t_orb_hd& orb_hd)
	{
		// add jdhuang
		// 有的时候，orb文件没有某卫星的数据，但是头文件会显示orb文件存在某颗卫星
		set<string> orb_data_sat_list;
		for (const auto& iter : _funct)
		{
			orb_data_sat_list.insert(iter.first);
		}

		_orb_hd.sat_prn = orb_data_sat_list;
		_orb_hd.sat_num = orb_data_sat_list.size();

		orb_hd = _orb_hd;
	}

	/**
	* @brief get orb funct data.
	* @param[out]   funct	whole orb funct data
	*/
	void t_gorb::get_funct(map<string, map<int, vector< t_gtriple> > >& funct)
	{
		funct = _funct;
	}

	/**
	* @brief get orb funct data of a satellite.
	* @param[in]   sat		satellite prn
	* @param[out]  funct	orb funct data of the sat
	*/
	void t_gorb::get_funct(string sat, map<int, t_sat_funct >& funct)
	{
		funct = _funct[sat];
	}

	/**
	* @brief intergral a satellite's position in specific epoch.
	* @param[in]   epoch	epoch time
	* @param[in]   sat		satellite prn
	* @param[out]  px_sat	x coordinates
	* @param[out]  py_sat	y coordinates
	* @param[out]  pz_sat	z coordinates
	* @return   if intergral successfully
	*/
	bool t_gorb::get_xyz(const t_gtime& epoch, string sat, double& px_sat, double& py_sat, double& pz_sat)
	{

		if (_funct.find(sat) == _funct.end())
		{
			return false;
		}
		// NOTICE::lagrange order = 8
		const int poly_order = 8;
		double epo_now = (epoch - _orb_hd.t_ref) / _orb_hd.admout;
		int epo_left = round(epo_now - (int)(poly_order / 2) + 1.e-5 / _orb_hd.admout);
		int epo_right = epo_left + poly_order;
		

		/*if (_funct[sat].find(epo_left) == _funct[sat].end() || _funct[sat].find(epo_right) == _funct[sat].end())
		{
			return false;
		}*/
		map<int, t_sat_funct> ::iterator itepo = _funct[sat].begin();
		if (epo_left < itepo->first)
		{
			epo_left = itepo->first;
			epo_right = epo_left + 8;
		}
		map<int, t_sat_funct> ::reverse_iterator ritepo = _funct[sat].rbegin();
		if (epo_right > ritepo->first)
		{
			epo_right = ritepo->first;
			epo_left = epo_right - 8;
		}

		// cycle for xyz 
		double	xyz[3];
		try
		{
			for (int crd = 0; crd < 3; crd++) {
				vector<double> Y;
				vector<double> X;
				for (int i = epo_left; i <= epo_right; i++) {
					Y.push_back(_funct[sat][i][0][crd]);
					X.push_back(i);
				}
				xyz[crd] = lagrange_interpolate(X, Y, epo_now);
			}
		}
		catch (const std::exception& e)
		{
			std::cout << e.what() << endl;
			return false;
		}
		// result
		px_sat = xyz[0] * 1000;
		py_sat = xyz[1] * 1000;
		pz_sat = xyz[2] * 1000;

		return true;

	}

	bool t_gorb::get_vel(const t_gtime& epoch, string sat, double& vx_sat, double& vy_sat, double& vz_sat)
	{

		if (_funct.find(sat) == _funct.end())
		{
			return false;
		}
		// NOTICE::lagrange order = 8
		const int poly_order = 8;
		double epo_now = (epoch - _orb_hd.t_ref) / _orb_hd.admout;
		int epo_left = round(epo_now - (int)(poly_order / 2) + 1.e-5 / _orb_hd.admout);
		int epo_right = epo_left + poly_order;

		/*if (_funct[sat].find(epo_left) == _funct[sat].end() || _funct[sat].find(epo_right) == _funct[sat].end())
		{
			return false;
		}*/
		map<int, t_sat_funct> ::iterator itepo = _funct[sat].begin();
		if (epo_left < itepo->first)
		{
			epo_left = itepo->first;
			epo_right = epo_left + 8;
		}
		map<int, t_sat_funct> ::reverse_iterator ritepo = _funct[sat].rbegin();
		if (epo_right > ritepo->first)
		{
			epo_right = ritepo->first;
			epo_left = epo_right - 8;
		}

		// cycle for xyz 
		double	vel[3];
		try
		{
			for (int crd = 0; crd < 3; crd++) {
				vector<double> Y;
				vector<double> X;
				for (int i = epo_left; i <= epo_right; i++) {
					Y.push_back(_funct[sat][i][0][crd]);
					X.push_back(i);
				}
				vel[crd] = lagrange_interpolate(X, Y, epo_now, true);
				vel[crd] /= _orb_hd.admout;
			}
		}
		catch (const std::exception& e)
		{
			std::cout << e.what() << endl;
			return false;
		}
		// result
		vx_sat = vel[0] * 1000;
		vy_sat = vel[1] * 1000;
		vz_sat = vel[2] * 1000;

		return true;

	}
	bool t_gorb::get_partial(const t_gtime& epoch, string sat, Matrix& partial)
	{
		// NOTICE::lagrange order = 8
		const int poly_order = 8;
		double epo_now = (epoch - _orb_hd.t_ref);
		epo_now /= (double)_orb_hd.admout;
		//epo_now = 0.99973255290999 - 1.0;
		int epo_left = floor(epo_now) - (int)(poly_order / 2) + 1;
		int epo_right = floor(epo_now) + (int)(poly_order / 2) + 1;

		if (_funct[sat].find(epo_left) == _funct[sat].end() || _funct[sat].find(epo_right) == _funct[sat].end())
		{
			return false;
		}

		// cycle for funct 
		partial.ReSize(_funct[sat][0].size() - 1, 3);
		try
		{
			for (int crd = 0; crd < 3; crd++) {
				for (int par_num = 1; par_num < _funct[sat][0].size(); par_num++)
				{
					vector<double> Y;
					vector<double> X;

					//cout << "par num" << par_num << endl;
					for (int i = epo_left; i <= epo_right; i++) {
						Y.push_back(_funct[sat][i][par_num][crd]);

						X.push_back(i);
					}
					partial(par_num, crd + 1) = lagrange_interpolate(X, Y, epo_now);
				}
			}
		}
		catch (const std::exception& e)
		{
			std::cout << e.what() << endl;
			return false;
		}
		// result
		return true;
	}

	// for LEO
	bool t_gorb::get_partial_leo(const t_gtime& epoch, string sat, Matrix& partial, vector<t_gpar>& orbpar)
	{
		const int poly_order = 8;
		// current epoch
		int epo_left, epo_right;
		double epo_now;
		if (!get_interp_points(epoch, poly_order, sat, epo_left, epo_right, epo_now))
		{
			return false;
		}
		//
		//vector<t_gpar> orbpar = icspar.getOrbPar();
		int npar = orbpar.size();
		Matrix partialPV_now(npar, 6), partialPV_t1(npar, 6);  // P V partial for current epoch and P V partial for X param begin time
		partialPV_now = { 0.0 };
		partialPV_t1 = { 0.0 };
		// get all partials for current epoch
		if (!get_partial_interp(epo_left, epo_right, epo_now, sat, partialPV_now))
		{
			return false;
		}
		//Debug
		////Debug
		/*cout << "Current epoch partial from interpolation" << endl;
		//cout << "Current epoch partial from interpolation" << endl;
		for (int i = 1; i <= partialPV_now.Nrows(); i++)
		//for (int i = 1; i <= partialPV_now.Nrows(); i++)
		{
			for (int j = 1; j <= partialPV_now.Ncols(); j++)
			{
				cout << setw(25) << scientific << setprecision(16) << partialPV_now(i, j);
			}
			cout << endl;
		}*/
		//{
		//	for (int j = 1; j <= partialPV_now.Ncols(); j++)
		//	{
		//		cout << setw(25) << scientific << setprecision(16) << partialPV_now(i, j);
		//	}
		//	cout << endl;
		//}
		//
		Matrix pvpart_now(6, 6), xpart_now(6, 1);   // [d(P, V)now/d(P, V)0], [d(P, V)now/d(X)i]
		pvpart_now = partialPV_now.submatrix(1, 6, 1, 6).t();
		// get partials before current epoch for all PWC parameters
		//map<int, map<int, Matrix>> pwcpart;  // index(in ICS param list: 13 in all currently)->index(in this kind of PWC param list)->partial[d(P, V)now/d(X)i]
		int indics = -1;  // index of param in all ICS param list(total 13 currently)(start from 0)
		//map<t_gpar::t_type, map<int, int>> activePWC;  // partype->nactive->index(in ICS param list: 13 in all currently)
		for (int ipar = 0; ipar < npar;)
		{
			indics++;
			if (!orbpar[ipar].lPWC)
			{
				for (int i = 1; i <= 3; i++)
				{
					partial(ipar + 1, i) = partialPV_now(indics + 1, i);   // non-PWC params
				}
				ipar++;
				continue;
			}
			xpart_now = partialPV_now.submatrix(indics + 1, indics + 1, 1, 6).t();
			int nactive = 0;
			t_gpar parpwc;
			Matrix xpart_b(6, 1), xpart_f(6, 1)/*, xpart_m(6, 1)*/;   // partial matrix of two adjacent PWC param for one kind of PWC param(back and forward) and xpart_f minus xpart_b
			xpart_b = { 0.0 };
			xpart_f = { 0.0 };

			// search for number of active PWC param for this kind of PWC param in current epoch
			int nPWC = orbpar[ipar].nPWC;
			double duration_beg = epoch - orbpar[ipar].beg;
			double duration_end = epoch - orbpar[ipar + nPWC - 1].end;
			if (duration_beg < 0)
			{
				if ((duration_beg + 31 * 60.0) < 0)
				{
					ipar += nPWC;
					continue;
				}
				else
				{
					nactive = 1;
				}
			}
			else if (duration_end > 0)
			{
				if ((duration_end - 31 * 60.0) >= 0)
				{
					ipar += nPWC;
					continue;
				}
				else
				{
					nactive = nPWC;
				}
			}
			else
			{
				nactive = floor(duration_beg / (orbpar[ipar].end - orbpar[ipar].beg)) + 1;
			}

			// calculate the partial matrix
			for (int iPWC = 0; iPWC < nactive; iPWC++)
			{
				xpart_b = xpart_f;
				parpwc = orbpar[ipar];
				/*double dt1 = epoch - parpwc.beg;
				double dt2 = epoch - parpwc.end;
				if (iPWC == 0)
				{
					dt1 += 31 * 60.0;
				}*/
				if (iPWC == nactive - 1)
				{
					//dt2 -= 31 * 60.0;
					xpart_f = xpart_now;
					/*if (!get_interp_points(parpwc.end, poly_order, sat, epo_left, epo_right, epo_now))
					{
						return false;
					}*/
				}
				else
				{
					if (!get_interp_points(parpwc.end, poly_order, sat, epo_left, epo_right, epo_now))
					{
						return false;
					}

					if (!get_partial_interp(epo_left, epo_right, epo_now, sat, partialPV_t1))
					{
						return false;
					}
					//Debug
					/*cout << "PWC param start epoch partial from interpolation" << endl;
					for (int i = 1; i <= partialPV_t1.Nrows(); i++)
					{
					for (int j = 1; j <= partialPV_t1.Ncols(); j++)
					{
					cout << setw(25) << scientific << setprecision(16) << partialPV_t1(i, j);
					}
					cout << endl;
					}*/
					//
					////Debug
					//cout << "PWC param start epoch partial from interpolation" << endl;
					//for (int i = 1; i <= partialPV_t1.Nrows(); i++)
					//{
					//for (int j = 1; j <= partialPV_t1.Ncols(); j++)
					//{
					//cout << setw(25) << scientific << setprecision(16) << partialPV_t1(i, j);
					//}
					//cout << endl;
					//}
					////
					// i for begin time of current PWC param
					Matrix pvpart(6, 6), xpart(6, 1), xpart0(6, 1);  // [d(P, V)i/d(P, V)0], [d(P, V)i/d(X)i], [d(P, V)0/d(X)j]
					pvpart = partialPV_t1.submatrix(1, 6, 1, 6).t();
					/*cout << "pvpart" << endl;
					for (int i = 1; i <= 6; i++)
					{
						for (int j = 1; j <= 6; j++)
						{
							cout << setw(25) << setprecision(16) << scientific << pvpart(i, j) << endl;
						}
					}*/
					xpart = partialPV_t1.submatrix(indics + 1, indics + 1, 1, 6).t();
					/*cout << "xpart" << endl;
					for (int i = 1; i <= 6; i++)
					{
						cout << setw(25) << setprecision(16) << scientific << xpart(i, 1) << endl;
					}*/
					xpart0 = pvpart.i() * xpart;  // [d(P, V)0/d(X)j] = inv([d(P, V)i/d(P, V)0]) * [d(P, V)i/d(X)i]
					//pwcpart[indics][iPWC] = xpart0;
					/*cout << "xpart0" << endl;
					for (int i = 1; i <= 6; i++)
					{
						cout << setw(25) << setprecision(16) << scientific << xpart0(i, 1) << endl;
					}*/

					xpart_f = pvpart_now * xpart0;  //[d(P, V)now/d(X)i] = [d(P, V)now/d(P, V)0] * [d(P, V)0/d(X)i]
					//xpart_m = xpart_f - xpart_b;
				}
				//Debug
				/*cout << "xpart_f" << endl;
				for (int i = 1; i <= 6; i++)
				{
					cout << setw(25) << setprecision(16) << scientific << xpart_f(i, 1) << endl;
				}*/
				//
				for (int k = 1; k <= 3; k++)
				{
					partial(ipar + 1, k) = xpart_f(k, 1) - xpart_b(k, 1);
				}
				ipar++;
				//pwcpart[indics][iPWC] = xpart1;
				//nactive++;
			}
			ipar += (nPWC - nactive);
			//activePWC[parpwc.parType][nactive] = indics;
		}

		// result
		return true;
	}

	//
	bool t_gorb::get_partial_interp(int epo_left, int epo_right, double epo_now, string sat, Matrix& partial)
	{
		try
		{
			for (int crd = 0; crd < 3; crd++) {
				for (int par_num = 1; par_num < _funct[sat][0].size(); par_num++)
				{
					vector<double> Y;
					vector<double> X;

					//cout << "par num" << par_num << endl;
					for (int i = epo_left; i <= epo_right; i++) {
						Y.push_back(_funct[sat][i][par_num][crd]);
						X.push_back(i);
					}
					partial(par_num, crd + 1) = lagrange_interpolate(X, Y, epo_now);
					partial(par_num, crd + 4) = lagrange_interpolate(X, Y, epo_now, true) / _orb_hd.admout;
				}
			}
		}
		catch (const std::exception& e)
		{
			std::cout << e.what() << endl;
			return false;
		}
		return true;
	}

	//
	bool t_gorb::get_interp_points(const t_gtime& epoch, int norder, string sat, int& epo_left, int& epo_right, double& epo_now)
	{
		epo_now = (epoch - _orb_hd.t_ref) / (double)_orb_hd.admout;
		epo_left = round(epo_now - (int)(norder / 2) + 1.e-5 / _orb_hd.admout);
		epo_right = epo_left + norder;
		if (_funct[sat].find(epo_left) == _funct[sat].end() || _funct[sat].find(epo_right) == _funct[sat].end())
		{
			return false;
		}
		return true;
	}

	void t_gorb::get_parlist(map<string,vector<string> >& parlist)
	{
		parlist = _par_list;
	}
}