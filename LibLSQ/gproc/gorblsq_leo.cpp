/**
*
* @verbatim
History
-1.0 zhangwei  2019-03-18  creat
-1.1 Qinyujie  2019-03-20  optim
-1.2 Qinyujie  2019-04-11  add Doxygen style code
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gorblsq_leo.cpp
* @brief		lsq processing for leo
*
* @author       zhangwei, Wuhan University
* @version		1.0.0
* @date		    2019-03-18
*
*/
#include "gproc/gorblsq_leo.h"
#include "gutils/gleo.h"
#include "gset/gsetout.h"


namespace great
{

	/** @brief default constructor. */
	t_gorblsq_leo::t_gorblsq_leo(t_gsetbase* set, t_gallproc* data, t_glog* log) :
		t_glsqprocIF(set, data, log),
		_result(0)
	{
		try
		{
			//_sites = dynamic_cast<t_gsetgen*>(set)->recskd();
			for (auto iter = _rec_list.begin(); iter != _rec_list.end(); iter++)
			{
				string satleo = strs2strl(*iter);
				_satleo.insert(satleo);
			}
			_setOut();
		}
		catch (...)
		{
			if (_glog)
			{
				_glog->comment(1, "ERROR:  t_gorblsq_leo::t_gorblsq_leo() constructor  unknown mistake");
			}
			else
			{
				cout << "ERROR:  t_gorblsq_leo::t_gorblsq_leo() constructor  unknown mistake" << endl;
			}
		}

	}

	/** @brief default destructor. */
	t_gorblsq_leo::~t_gorblsq_leo()
	{
		// delete _result
		if (_result)
		{
			if (_result->is_open())
			{
				_result->flush();
				_result->close();
			}

			delete _result;
			_result = NULL;
		}
	}


	bool t_gorblsq_leo::GenerateProduct()
	{
		_extract_clkfile(*_gall_rcv, t_gallprec::AR);

		_update_icsleo();

		return true;
	}

	bool t_gorblsq_leo::_InitLsq(t_gallproc* data)
	{
		_icsleo = dynamic_cast<t_gallics_leo*>((*data)[t_gdata::ALLICSLEO]);
		if (!_icsleo)
		{
			if (_glog) {
				_glog->comment(1, "icsleo data not found!!!");
			}
			return false;
		}
		set<string> longsite = _icsleo->get_leo_sat_list();
		set<string> shortsite;
		for (auto site = longsite.begin(); site != longsite.end(); site++)
		{
			shortsite.insert(strl2strs(*site));
		}
		_rec_list = shortsite;
		if (_rec_list.size() == 0) {
			if (_glog) {
				_glog->comment(1, "icsleo file have no leo data!");
			}
			return false;
		}

		_init_parameter();
		_init_ORB_parameter();

		return true;
	}


	/** @brief init LEO orbit parameters. */
	void t_gorblsq_leo::_init_ORB_parameter()
	{
		try
		{
			int ipar_num = _lsq->_x_solve.parNumber();
			// cycle sat list
			for (auto sat_iter = _satleo.begin(); sat_iter != _satleo.end(); sat_iter++)
			{
				auto ics = _icsleo->get_leo_ics(*sat_iter);
				map<ICS_PAR, vector<t_icsleo_record>>* leopar;
				leopar = ics->leoics();
				int npar = 0;
				vector<t_gpar> pars;
				// cycle par list
				for (auto par_iter = leopar->begin(); par_iter != leopar->end(); par_iter++)
				{
					t_gpar::t_type par_type = t_gics_leo::icspar2gpar(par_iter->first);
					vector<t_icsleo_record> onepar = par_iter->second;
					for (int k = 0; k < onepar.size(); k++)
					{
						//int indexPWC = 0;
						int npwc = ceil((_end_time - _beg_time) / (onepar[0].end - onepar[0].beg));
						if (onepar[k].beg >= _beg_time && onepar[k].beg < _end_time || onepar[k].end > _beg_time && onepar[k].end <= _end_time)
						{
							npar++;
							t_gpar par_temp(strl2strs(*sat_iter), par_type, ++ipar_num, "", false);
							par_temp.setTime(onepar[k].beg, onepar[k].end);
							par_temp.value(onepar[k].value);
							par_temp.apriori(onepar[k].sig);
							par_temp.value(onepar[k].value);
							if (double_eq(onepar[k].sig, 0.0))
							{
								switch (par_type)
								{
								case t_gpar::PXSAT:
								case t_gpar::PYSAT:
								case t_gpar::PZSAT:
									par_temp.apriori(100.0); break;
								case t_gpar::VXSAT:
								case t_gpar::VYSAT:
								case t_gpar::VZSAT:
									par_temp.apriori(1.0); break;
								default:
									par_temp.apriori(0.1);
									break;
								}
							}
							else
							{
								par_temp.apriori(onepar[k].sig);
							}
							if (onepar.size() > 1)  // PWC parameters
							{
								par_temp.nPWC = npwc;
								par_temp.lPWC = true;
								//par_temp.indPWC = indexPWC;
								//indexPWC++;
							}
							_lsq->add_parameter(par_temp);
							pars.push_back(par_temp);
						}

					}
					_lsq->_x_solve.setOrbPar(*sat_iter, pars);
				}
			}
		}
		catch (...)
		{
			if (_glog)
			{
				_glog->comment(1, "ERROR:  t_gorblsq_leo::_init_ORB_parameter  unknown mistake");
			}
			else
			{
				cout << "ERROR:  t_gorblsq_leo::_init_ORB_parameter  unknown mistake" << endl;
			}
		}
	}

	/** @brief init other parameters except orbit params. */
	void t_gorblsq_leo::_init_parameter()
	{
		try
		{
			gtrace("t_gorblsq_leo::_init_parameter");

			if (_gall_obj == nullptr) {
				throw exception();
			}

			// count par number
			int ipar_num = _lsq->_x_solve.parNumber();

			// init recclk isb/ifb
			// loop loop over sites from settings
			for (auto it = _rec_list.begin(); it != _rec_list.end(); it++)
			{
				//string site = it->first;
				string site = *it;

				// init the ISB/IFB 
				// TODO  set in estimate equations
				if (_sys_list.size() > 1)
				{
					set<string>::iterator iter_sys = _sys_list.begin();
					// for scene that haven't GPS
					if (_sys_list.find("GPS") == _sys_list.end())
						iter_sys++;
					// set ISB
					for (; iter_sys != _sys_list.end(); iter_sys++)
					{
						t_gpar* par_ISB;
						GSYS sys = t_gsys::str2gsys(*iter_sys);
						switch (sys)
						{
						case gnut::GPS:
							break;
						case gnut::GAL:
							par_ISB = new t_gpar(site, t_gpar::GAL_ISB, ++ipar_num, "");
							par_ISB->value(0.0);
							par_ISB->setTime(_beg_time, _end_time);
							par_ISB->apriori(3000.0);
							_lsq->add_parameter(*par_ISB);
							break;
						case gnut::GLO:
							par_ISB = new t_gpar(site, t_gpar::GLO_ISB, ++ipar_num, "");
							par_ISB->value(0.0);
							par_ISB->setTime(_beg_time, _end_time);
							par_ISB->apriori(3000.0);
							_lsq->add_parameter(*par_ISB);
							break;
						case gnut::BDS:
							par_ISB = new t_gpar(site, t_gpar::BDS_ISB, ++ipar_num, "");
							par_ISB->value(0.0);
							par_ISB->setTime(_beg_time, _end_time);
							par_ISB->apriori(3000.0);
							_lsq->add_parameter(*par_ISB);
							break;
						case gnut::QZS:
							par_ISB = new t_gpar(site, t_gpar::QZS_ISB, ++ipar_num, "");
							par_ISB->value(0.0);
							par_ISB->setTime(_beg_time, _end_time);
							par_ISB->apriori(3000.0);
							_lsq->add_parameter(*par_ISB);
							break;
						default:
							break;
						}
						delete par_ISB;
					}
				}
				// set GLO IFB
				if (_sys_list.find("GLO") != _sys_list.end())
				{
					set<string> glo_sats = GNSS_SATS()[GLO];
					for (auto sat_iter = glo_sats.begin(); sat_iter != glo_sats.end(); sat_iter++)
					{
						t_gpar par_GLOIFB(site, t_gpar::GLO_IFPB, ++ipar_num, *sat_iter);
						par_GLOIFB.value(0.0);
						par_GLOIFB.setTime(_beg_time, _end_time);
						par_GLOIFB.apriori(3000.0);
						_lsq->add_parameter(par_GLOIFB);
					}
				}

				//intit the Clock
				t_gpar par_clk(site, t_gpar::CLK, ++ipar_num, "");
				par_clk.value(0.0);
				par_clk.setTime(_beg_time, _beg_time);
				par_clk.apriori(9000);
				_lsq->add_parameter(par_clk);
			}
		}
		catch (...)
		{
			if (_glog)
			{
				_glog->comment(1, "ERROR:  t_gorblsq_leo::_init_parameter  unknown mistake");
			}
			else
			{
				cout << "ERROR:  t_gorblsq_leo::_init_parameter  unknown mistake" << endl;
			}
		}
	}

	/** @brief init the output file. */
	void t_gorblsq_leo::_setOut()
	{
		try
		{
			// set output file
			string tmp;
			tmp = dynamic_cast<t_gsetout*>(_gset)->outputs("flt");
			if (!tmp.empty()) {
				substitute(tmp, "$(date)", int2str(_beg_time.year()) + int2str(_beg_time.doy()), false);
				_result = new t_giof;
				_result->tsys(t_gtime::GPS);
				_result->mask(tmp);
				_result->append(dynamic_cast<t_gsetout*>(_gset)->append());
			}

		}
		catch (...)
		{
			if (_glog)
			{
				_glog->comment(1, "ERROR:  t_gorblsq_leo::_setOut  unknown mistake");
			}
			else
			{
				cout << "ERROR:  t_gorblsq_leo::_setOut  unknown mistake" << endl;
			}
		}
	}

	/** @brief upsate LEO ics information. */
	void t_gorblsq_leo::_update_icsleo()
	{
		try
		{
			// update
			for (auto sat_iter = _satleo.begin(); sat_iter != _satleo.end(); sat_iter++)
			{

				auto ics = _icsleo->get_leo_ics(*sat_iter);
				map<ICS_PAR, vector<t_icsleo_record>>* par_list = ics->leoics();

				// cycle par list
				for (auto par_iter = par_list->begin(); par_iter != par_list->end(); par_iter++)
				{
					t_gpar::t_type par_type = t_gics_leo::icspar2gpar(par_iter->first);
					vector<t_icsleo_record> &pars = par_iter->second;
					for (int k = 0; k < pars.size(); k++)
					{
						t_gtime begtime = pars[k].beg;
						t_gtime endtime = pars[k].end;
						int par_idx = _lsq->_x_solve.getParam(strl2strs(*sat_iter), par_type, "", begtime, endtime);
						if (par_idx >= 0)
						{
							ColumnVector dx = _lsq->dx();
							ColumnVector stdx = _lsq->stdx();

							// update value && sig
							pars[k].value += dx(par_idx + 1);
							pars[k].sig = stdx(par_idx + 1);

						}

					}

				}

			}
		}
		catch (...)
		{
			if (_glog)
			{
				_glog->comment(1, "ERROR:  t_gorblsq_leo::_update_icsleo  unknown mistake");
			}
			else
			{
				cout << "ERROR:  t_gorblsq_leo::_update_icsleo  unknown mistake" << endl;
			}
		}
	}


}