/* ----------------------------------------------------------------------
Purpose: save data of orb file
Version: $ Rev: $ 1.0

2019-04-10 /HJ_Zheng: created

-*/
/**
 *
 * @verbatim
	 History
	  -1.0  Hongjie Zheng   2019-04-10 Creat the file.
	  -1.1  Jie Li			2019-07-04 add doxygen style code remarks
   @endverbatim
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		gorb.h
 * @brief		The class for storaging all orb data.
 *
 * @author      Hongjie Zheng, Wuhan University
 * @version		1.0.0
 * @date		2019-04-10
 *
 */

#include "gall/gallorb.h"
#include "gmodels/gpoly.h"
#include "gutils/gleo.h"
#include "gutils/ginfolog.h"

using namespace std;
using namespace gnut;

namespace great
{

	/** @brief default constructor. */
	t_gallorb::t_gallorb() :
		t_gdata()
	{
		gtrace("t_gallorb::constructor");
		id_type(t_gdata::ALLORB);
	}
	/** @brief default destructor. */
	t_gallorb:: ~t_gallorb()
	{

	}
	/**
	* @brief add orb data include header and body of one satellite to map.
	* @param[in]   epo		epoch
	* @param[in]   sat		satellite prn
	* @param[in]   orb_hd   orb head data
	* @param[in]   funct    position and patial matrix
	*/
	void t_gallorb::addData(int epo, string sat, const t_orb_hd& orb_hd, const t_gtriple& funct)
	{
		_orbs[sat].first = orb_hd;
		_orbs[sat].second[epo].push_back(funct);
	}

	/**
	* @brief add orb data include header and body of one satellite to map.
	* @param[in]   epo		epoch
	* @param[in]   sat		satellite prn
	* @param[in]   orb_hd   orb head data
	* @param[in]   funct    position and patial matrix
	*/
	void t_gallorb::addData(int epo, string sat, const t_orb_hd& orb_hd, const vector<t_gtriple>& funct)
	{
		_orbs[sat].first = orb_hd;
		_orbs[sat].second[epo] = funct;
	}

	/**
	* @brief interpolate a satellite's position(unit: m) in specific epoch.
	* @param[in]   epoch	epoch time
	* @param[in]   sat		satellite prn
	* @param[out]  px_sat	x coordinates
	* @param[out]  py_sat	y coordinates
	* @param[out]  pz_sat	z coordinates
	* @return   if interpolate successfully
	*	@retval   true   interpolate successfully
	*   @retval   false  interpolate unsuccessfully
	*/
	bool t_gallorb::get_crs_sat_crd(const t_gtime& epoch, const string& sat, double& px_sat, double& py_sat, double& pz_sat) 
	{

		if (_orbs.find(sat) == _orbs.end())
		{
			return false;
		}
		// NOTICE::lagrange order = 8
		const t_orb_hd& hd_info = _orbs[sat].first;
		map<int, vector<t_gtriple> >& funct = _orbs[sat].second;

		// NOTICE::lagrange order = 8
		const int poly_order = 8;
		int epo_left, epo_right;
		double epo_now;
		if (!getInterpPoints(epoch, poly_order, sat, epo_left, epo_right, epo_now))
		{
				return false;
		}

		map<int, t_sat_funct> ::iterator itepo = funct.begin();
		if (epo_left < itepo->first)
		{
			write_log_info(_log, 0, "ERROR", "Cannot perform the lagrange interpolation: input epoch is beyond the epoch table of orb file");
			return false;
		}
		map<int, t_sat_funct> ::reverse_iterator ritepo = funct.rbegin();
		if (epo_right > ritepo->first)
		{
			write_log_info(_log, 0, "ERROR", "Cannot perform the lagrange interpolation: input epoch is beyond the epoch table of orb file");
			return false;
		}

		// cycle for xyz 
		double	xyz[3];
		try
		{
			vector<double> Y;
			vector<double> X;
			Y.resize(epo_right - epo_left + 1);
			X.resize(epo_right - epo_left + 1);
			for (int crd = 0; crd < 3; crd++) 
            {
				for (int i = epo_left; i <= epo_right; i++) {
					Y[i - epo_left] = funct[i][0][crd];
					X[i - epo_left] = i;
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

	/**
	* @brief interpolate a satellite's velocity(unit: m/s) in specific epoch.
	* @param[in]   epoch	epoch time
	* @param[in]   sat		satellite prn
	* @param[out]  vx_sat	x velocity
	* @param[out]  vy_sat	y velocity
	* @param[out]  vz_sat	z velocity
	* @return   if interpolate successfully
	*	@retval   true   interpolate successfully
	*   @retval   false  interpolate unsuccessfully
	*/
	bool t_gallorb::get_crs_sat_vel(const t_gtime& epoch, string sat, double& vx_sat, double& vy_sat, double& vz_sat) {

		if (_orbs.find(sat) == _orbs.end())
		{
			return false;
		}
		// NOTICE::lagrange order = 8
		t_orb_hd& orb_hd = _orbs[sat].first;
		auto& orb_data = _orbs[sat].second;
		const int poly_order = 8;
		int epo_left, epo_right;
		double epo_now;
		if (!getInterpPoints(epoch, poly_order, sat, epo_left, epo_right, epo_now))
		{
				return false;
		}
		//double epo_now = (epoch - orb_hd.t_ref) / orb_hd.admout;
		//int epo_left = floor(epo_now) - (int)(poly_order / 2) + 1;
		//int epo_right = floor(epo_now) + (int)(poly_order / 2) + 1;

		/*if (_funct[sat].find(epo_left) == _funct[sat].end() || _funct[sat].find(epo_right) == _funct[sat].end())
		{
			return false;
		}*/
		map<int, t_sat_funct> ::iterator itepo = orb_data.begin();
		if (epo_left < itepo->first)
		{
			//epo_left = itepo->first; Commented by zhangwei: lagrange interpolation cannot be performed if the current epoch beyond the epoch range of orb file
			//epo_right = epo_left + 8;
			write_log_info(_log, 0, "ERROR", "Cannot perform the lagrange interpolation: input epoch is beyond the epoch table of orb file");
			return false;
		}
		map<int, t_sat_funct> ::reverse_iterator ritepo = orb_data.rbegin();
		if (epo_right > ritepo->first)
		{
			//epo_right = ritepo->first; Commented by zhangwei: lagrange interpolation cannot be performed if the current epoch beyond the epoch range of orb file
			//epo_left = epo_right - 8;
			write_log_info(_log, 0, "ERROR", "Cannot perform the lagrange interpolation: input epoch is beyond the epoch table of orb file");
			return false;
		}


		// cycle for xyz 
		double	vel[3];
		try
		{
			vector<double> Y;
			vector<double> X;
			Y.resize(epo_right - epo_left + 1);
			X.resize(epo_right - epo_left + 1);
			for (int crd = 0; crd < 3; crd++) {
				//vector<double> Y;
				//vector<double> X;
				for (int i = epo_left; i <= epo_right; i++) {
					//Y.push_back(_allorb[sat].second[i][0][crd]);
					//X.push_back(i);
					Y[i - epo_left] = orb_data[i][0][crd];
					X[i - epo_left] = i;
				}
				vel[crd] = lagrange_interpolate(X, Y, epo_now, true);
				vel[crd] /= orb_hd.admout;
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

	/**
	* @brief interpolate a satellite's partial matrix(only dxdx0..dxdvx0...dxdnp) in specific epoch.
	* @param[in]   epo		epoch time
	* @param[in]   sat		satellite prn
	* @param[out]  partial	partial matrix
	* @return   if interpolate successfully
	*	@retval   true   interpolate successfully
	*   @retval   false  interpolate unsuccessfully
	*/
	bool t_gallorb::getGnsPartial(const t_gtime& epoch, string sat, Matrix& partial) 
	{
		if (_orbs.find(sat) == _orbs.end())
		{
			return false;
		}
		// NOTICE::lagrange order = 8
		const t_orb_hd& hd_info = _orbs[sat].first;
		map<int, vector<t_gtriple> >& funct = _orbs[sat].second;

		const int poly_order = 8;
		int epo_left, epo_right;
		double epo_now;
		if (!getInterpPoints(epoch, poly_order, sat, epo_left, epo_right, epo_now))
		{
		    return false;
		}
		//double epo_now = (epoch - hd_info.t_ref);
		//epo_now /= (double)hd_info.admout;
		//epo_now = 0.99973255290999 - 1.0;
		//int epo_left = floor(epo_now) - (int)(poly_order / 2) + 1;
		//int epo_right = floor(epo_now) + (int)(poly_order / 2) + 1;

		if (funct.find(epo_left) == funct.end() || funct.find(epo_right) == funct.end())
		{
			return false;
		}

		// cycle for funct 
		partial.ReSize(funct[0].size() - 1, 3);
		try
		{
			vector<double> Y;
			vector<double> X;
			Y.resize(epo_right - epo_left + 1);
			X.resize(epo_right - epo_left + 1);
			for (int crd = 0; crd < 3; crd++) 
			{
				for (int par_num = 1; par_num < funct[0].size(); par_num++)
				{
					//vector<double> Y;
					//vector<double> X;

					//cout << "par num" << par_num << endl;
					for (int i = epo_left; i <= epo_right; i++) 
					{
						//Y.push_back(funct[i][par_num][crd]);
						//X.push_back(i);
						Y[i - epo_left] = funct[i][par_num][crd];
						X[i - epo_left] = i;
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

	/**
	* @brief interpolate a LEO satellite's partial matrix(especially for PWC parameters) in specific epoch.
	* @param[in]   epoch	epoch time
	* @param[in]   sat		satellite prn
	* @param[in]   orbpar	par list
	* @param[out]  partial	partial matrix
	* @return   if interpolate successfully
	*	@retval   true   interpolate successfully
	*   @retval   false  interpolate unsuccessfully
	*/
	bool t_gallorb::getLeoPartial(const t_gtime& epoch, string sat, Matrix& partial, const vector<t_gpar>& orbpar)
	{
		const int poly_order = 8;
		// current epoch
		int epo_left = 0, epo_right = 0;
		double epo_now = 0.0;
		static bool first = true;

		static vector<string> leolist;  ///< Modified by zhangwei
		if (leolist.size() == 0)
		{
			leolist.push_back(sat);
		}
		else
		{
			if (std::find(leolist.begin(), leolist.end(), sat) != leolist.end())
			{
				first = false;
			}
			else
			{
				leolist.push_back(sat);
				first = true;
			}
		}

		static map<string, map<t_gpar::t_type, map<t_gtime, Matrix> > > Partial_PV;
		if (first) {
			Matrix Partial_PV_init(orbpar.size(), 6);
			for (auto sat_iter = _orbs.begin(); sat_iter != _orbs.end(); sat_iter++)
			{
				if(sat != sat_iter->first)
				{
                    continue;
				}
				/*if (!IS_SUPPORTED_LEO(sat)) {
					continue;
				}*/

				for (auto par_iter = orbpar.begin(); par_iter != orbpar.end(); par_iter++)
				{
					if (!par_iter->lPWC) {
						continue;
					}

					if (!getInterpPoints(par_iter->end, poly_order, sat, epo_left, epo_right, epo_now))
					{
						continue;
					}

					if (!getPartialInterp(epo_left, epo_right, epo_now, sat, Partial_PV_init))
					{
						continue;
					}

					Partial_PV[sat][par_iter->parType][par_iter->end] = Partial_PV_init;
				}

			}
			//first = false;
		}

		if (!getInterpPoints(epoch, poly_order, sat, epo_left, epo_right, epo_now))
		{
			return false;
		}
		//
		//vector<t_gpar> orbpar = icspar.getOrbPar();
		int npar = orbpar.size();
		Matrix partialPV_now(npar, 6), partialPV_t1(npar, 6);  // P V partial for current epoch and P V partial for X param begin time
		partialPV_now = { 0.0 };
		partialPV_t1 = { 0.0 };
		Matrix pvpart_now(6, 6), xpart_now(6, 1);   // [d(P, V)now/d(P, V)0], [d(P, V)now/d(X)i]
		
		/*if (npar != 0)
		{*/
		// get all partials for current epoch
		if (!getPartialInterp(epo_left, epo_right, epo_now, sat, partialPV_now))
		{
			return false;
		}
		pvpart_now = partialPV_now.submatrix(1, 6, 1, 6).t();
		//}
		
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
		    else if (duration_end >= 0)
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
				nactive = floor(duration_beg / orbpar[ipar].piecewise_time / 60) + 1;
				if (nactive > nPWC)
				{
					nactive = nPWC;
					//todo,yjqin
				}
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
					//if (!get_interp_points(parpwc.end, poly_order, sat, epo_left, epo_right, epo_now))
					//{
					//	return false;
					//}

					//if (!get_partial_interp(epo_left, epo_right, epo_now, sat, partialPV_t1))
					//{
					//	return false;
					//}
					if (Partial_PV[sat][parpwc.parType].find(parpwc.end) == Partial_PV[sat][parpwc.parType].end()) {
						return false;
					}
					partialPV_t1 = Partial_PV[sat][parpwc.parType][parpwc.end];

					// i for begin time of current PWC param
					Matrix pvpart(6, 6), xpart(6, 1), xpart0(6, 1);  // [d(P, V)i/d(P, V)0], [d(P, V)i/d(X)i], [d(P, V)0/d(X)j]
					pvpart = partialPV_t1.submatrix(1, 6, 1, 6).t();

					xpart = partialPV_t1.submatrix(indics + 1, indics + 1, 1, 6).t();

					xpart0 = pvpart.i() * xpart;  // [d(P, V)0/d(X)j] = inv([d(P, V)i/d(P, V)0]) * [d(P, V)i/d(X)i]


					xpart_f = pvpart_now * xpart0;  //[d(P, V)now/d(X)i] = [d(P, V)now/d(P, V)0] * [d(P, V)0/d(X)i]
					//xpart_m = xpart_f - xpart_b;
				}

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

	/**
	* @brief interpolate a satellite's partial matrix(both current pos and vel) in specific epoch.
	* @param[in]   epo_left		first interpolate point
	* @param[in]   epo_right	last interpolate point
	* @param[in]   epo_now		current epoch time
	* @param[in]   sat			satellite prn
	* @param[out]  partial		partial matrix(pos and vel)
	* @return   if interpolate successfully
	*	@retval   true   interpolate successfully
	*   @retval   false  interpolate unsuccessfully
	*/
	bool t_gallorb::getPartialInterp(int epo_left, int epo_right, double epo_now, string sat, Matrix& partial)
	{
		try
		{
			auto& orb_hd = _orbs[sat].first;
			auto& orb_data = _orbs[sat].second;
			int par_size = _orbs[sat].second[0].size();
			vector<double> Y;
			vector<double> X;
			X.resize(epo_right - epo_left + 1);
			Y.resize(epo_right - epo_left + 1);
			vector<t_sat_funct*> funct_data;
			funct_data.resize(epo_right - epo_left + 1);
			for (int i = epo_left; i <= epo_right; i++) {
				funct_data[i - epo_left] = &orb_data[i];
			}

			for (int crd = 0; crd < 3; crd++) {
				/*for (int par_num = 1; par_num < _allorb[sat].second[0].size(); par_num++)*/
				for (int par_num = 1; par_num < par_size; par_num++)
				{
					//vector<double> Y;
					//vector<double> X;

					//cout << "par num" << par_num << endl;
					for (int i = epo_left; i <= epo_right; i++) {
						//Y.push_back(_allorb[sat].second[i][par_num][crd]);
						//X.push_back(i);
						//Y[i - epo_left] = orb_data[i][par_num][crd];
						Y[i - epo_left] = (*funct_data[i - epo_left])[par_num][crd];
						X[i - epo_left] = i;
					}
					partial(par_num, crd + 1) = lagrange_interpolate(X, Y, epo_now);
					//partial(par_num, crd + 4) = lagrange_interpolate(X, Y, epo_now, true) / _allorb[sat].first.admout;
					partial(par_num, crd + 4) = lagrange_interpolate(X, Y, epo_now, true) / orb_hd.admout;
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

	/**
	* @brief get the interpolate point
	* @param[in]   epo			epoch input(int)
	* @param[in]   norder		lagrange interpolate order
	* @param[in]   sat			satellite prn
	* @param[out]  epo_left		first interpolate point
	* @param[out]  epo_right	last interpolate point
	* @param[out]  epo_now		current epoch time(mjd)
	* @return   if get interpolate point successfully
	*	@retval   true   get interpolate poin successfully
	*   @retval   false  get interpolate poin unsuccessfully
	*/
	bool t_gallorb::getInterpPoints(const t_gtime& epoch, const int& norder, const string& sat, int& epo_left, int& epo_right, double& epo_now)
	{
		if (_orbs.find(sat) == _orbs.end())
		{
			return false;
		}
		t_orb_hd& orb_hd = _orbs[sat].first;
		epo_now = (epoch - orb_hd.t_ref) / (double)orb_hd.admout;
		epo_left = round(epo_now - (int)(norder / 2) + 1.e-5 / orb_hd.admout);
		epo_right = epo_left + norder;
		
		const map<int, t_sat_funct >& funct = _orbs[sat].second;
		if (funct.find(epo_left) == funct.end() || funct.find(epo_right) == funct.end())
		{
			return false;
		}
		return true;
	}

	/**
	* @brief get begin, end and reference time of one sat.
	* @param[in]   sat		satellite prn
	* @param[in]   beg		begin time
	* @param[in]   end		end time
	* @param[in]   ref		reference time
	*/
	void t_gallorb::getTime(string sat, t_gtime& beg, t_gtime& end, t_gtime& ref)
	{
		if (sat.empty())
		{
			auto sats = getSatlist();
			if (sats.empty()) return;

			t_gtime max_beg = FIRST_TIME;
			t_gtime min_end = LAST_TIME;
			for (const auto& iter : sats)
			{
				t_gtime tmp_beg = _orbs[iter].first.beg;
				t_gtime tmp_end = _orbs[iter].first.end;

				max_beg = (tmp_beg > max_beg) ? tmp_beg : max_beg;
				min_end = (tmp_end < min_end) ? tmp_end : min_end;
				ref = _orbs[iter].first.t_ref;

			}

			beg = max_beg;
			end = min_end;

		}
		else
		{
			beg = _orbs[sat].first.beg;
			end = _orbs[sat].first.end;
			ref = _orbs[sat].first.t_ref;
		}
	}

	/**
	* @brief get sat list.
	* @param[out]   sats		sat list
	*/
	set<string> t_gallorb::getSatlist()
	{
		set<string> sats;
		for (auto itsat = _orbs.begin(); itsat != _orbs.end(); ++itsat)
		{
			sats.insert((*itsat).first);
		}
		return sats;
	}

}
