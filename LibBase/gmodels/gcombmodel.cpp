/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2020-05-15 created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    t_gcombmodel.cpp
* @brief	    base combine biase model
* @mainpage     C++ version
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-05-18
*
*/

#include "gmodels/gcombmodel.h"
#include <math.h>


namespace great
{
	t_gcombmodel::t_gcombmodel(t_gsetbase* setting, t_glog* log, shared_ptr<t_gbiasmodel> bias_model, t_gallproc* data) :
		_bias_model(bias_model),
		_log(log)
	{
		_band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::GPS);
		_band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::GAL);
		_band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::GLO);
		_band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::BDS);
		_band_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::QZS);

		_freq_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::GPS);
		_freq_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::GAL);
		_freq_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::GLO);
		_freq_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::BDS);
		_freq_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::QZS);
		
		

		_frequency = dynamic_cast<t_gsetproc*>(setting)->frequency();
		_ion_model = dynamic_cast<t_gsetproc*>(setting)->ion_model();
		_ifb_model = dynamic_cast<t_gsetproc*>(setting)->ifb_model();
		_ifcb_model = dynamic_cast<t_gsetproc*>(setting)->ifcb_model();
		
		_gallbias = dynamic_cast<t_gallbias*>((*data)[t_gdata::ALLBIAS]);
		_gprecisemodel = dynamic_cast<t_gprecisemodel*>(_bias_model->precisemodel());

		_update_amb_lite = (dynamic_cast<t_gsetturboedit*>(setting)->liteMode() || dynamic_cast<t_gsetproc*>(setting)->realtime());

	}

	t_gcombmodel::~t_gcombmodel()
	{
	}


	map<GSYS, map<FREQ_SEQ, GOBSBAND>> t_gcombmodel::get_band_index()
	{
		return _band_index;
	}

	map<GSYS, map<GOBSBAND, FREQ_SEQ>> t_gcombmodel::get_freq_index()
	{
		return _freq_index;
	}

	t_gcombIF::t_gcombIF(t_gsetbase* setting, t_glog* log, shared_ptr<t_gbiasmodel> bias_model, t_gallproc* data) :
		t_gcombmodel(setting, log, bias_model, data)
	{
	}

	t_gcombIF::~t_gcombIF()
	{
	}

	bool t_gcombIF::cmb_equ(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, t_gbaseEquation& result)
	{
		// check Obs type
		GOBSBAND b1 = _band_index[obsdata.gsys()][FREQ_1];
		GOBSBAND b2 = _band_index[obsdata.gsys()][FREQ_2];


		auto band1 = GNSS_BAND_PRIORITY.at(obsdata.gsys())[1];
		auto band3 = GNSS_BAND_PRIORITY.at(obsdata.gsys())[3];

		//obsdata.apply_dcb(*_gallbias, (b1 == band1 && b2 == band3 && _frequency == 2));
		obsdata.apply_bias(*_gallbias);
		
		return this->cmb_equ_IF(epoch, params, obsdata, b1, b2, result);
		
	}

	bool t_gcombIF::cmb_equ_IF(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, GOBSBAND b1, GOBSBAND b2, t_gbaseEquation& result)
	{

		t_gobs gobsP1 = obsdata.select_range(b1);
		t_gobs gobsP2 = obsdata.select_range(b2);
		t_gobs gobsL1 = obsdata.select_phase(b1);
		t_gobs gobsL2 = obsdata.select_phase(b2);

		vector<pair<t_gobs, t_gobs> > type_list;
		type_list.push_back(make_pair(gobsP1, gobsP2));
		type_list.push_back(make_pair(gobsL1, gobsL2));

		// IF coef
		double coef1, coef2;
		obsdata.coef_ionofree(b1, coef1, b2, coef2);

		t_glsqEquationMatrix equ_IF;

		// cycle [Code and Phase]
		for (auto item : type_list)
		{
			t_gobs gobs1 = item.first;
			t_gobs gobs2 = item.second;

			//combine f1 and f2
			t_gbaseEquation temp_equ;
			if (!_bias_model->cmb_equ(epoch, params, obsdata, gobs1, temp_equ))
			{
				return false;
			}

			if (!_bias_model->cmb_equ(epoch, params, obsdata, gobs2, temp_equ))
			{
				return false;
			}

			vector<pair<int, double> > coef_IF;
			double P_IF = 0.0, l_IF = 0.0;

			if (temp_equ.B[0].size() != temp_equ.B[1].size())
			{
				throw logic_error("coeff size is not equal in f1 and f2");
			}
			for (int i = 0; i < temp_equ.B[0].size(); i++)
			{
				if (temp_equ.B[0][i].first != temp_equ.B[1][i].first) 
				{
					throw logic_error("coeff par is not the same in f1 and f2");
				}
				// combine coeff
				coef_IF.push_back(make_pair(temp_equ.B[0][i].first, coef1 * temp_equ.B[0][i].second + coef2 * temp_equ.B[1][i].second));
			}
			// combine P
			P_IF = 1.0 / (pow(coef1, 2) / temp_equ.P[0] + pow(coef2, 2) / temp_equ.P[1]);
			// combine l
			l_IF = coef1 * temp_equ.l[0] + coef2 * temp_equ.l[1];

			if (_ifcb_model == IFCB_MODEL::COR && _frequency == 3 && gobs1.is_phase() && _freq_index[obsdata.gsys()][b2] == FREQ_3)
			{
				apply_IFCB_GPS(obsdata, _gprecisemodel->_gifcb, l_IF, OBSCOMBIN::IONO_FREE);
			}

			// CORRECT AMB
			if (gobs1.is_phase() && (!obsdata.is_carrier_range(b1) || !obsdata.is_carrier_range(b2)))
			{
				t_gpar::t_type amb_type = t_gpar::AMB_IF;
				if (_freq_index[obsdata.gsys()][b2] == FREQ_3) 
				{
					amb_type = t_gpar::AMB13_IF;
				}

				int idx = params.getParam(obsdata.site(), amb_type, obsdata.sat());

				if (idx < 0)
				{
					return false;
				}

				// if amb is new  then init value with Obs - Modelobs
				//params[idx].value() = 0;
				
				if (double_eq(params[idx].value(), 0.0) || params[idx].beg == epoch)
				{
					double obs_P3 = obsdata.P3(gobsP1, gobsP2);
					double obs_L3 = obsdata.L3(gobsL1, gobsL2);					
				
					params[idx].value(obs_L3 - obs_P3);
				}
				// update B l
				coef_IF.push_back(make_pair(idx + 1, 1.0));
				l_IF -= params[idx].value();
			}

			// CORRECT CLK
			// remove 12 clk and isb
			if (_freq_index[obsdata.gsys()][b2] >= FREQ_3)
			{
				for (auto iter = coef_IF.begin(); iter != coef_IF.end();)
				{
					t_gpar::t_type type = params[iter->first - 1].parType;
					if (type == t_gpar::CLK ||
						type == t_gpar::GAL_ISB ||
						type == t_gpar::BDS_ISB ||
						type == t_gpar::QZS_ISB)
					{
						iter = coef_IF.erase(iter);
						continue;
					}
					iter++;
				}
			}

			// add 13 CLK
			if (_freq_index[obsdata.gsys()][b2] == FREQ_3)
			{
				t_gpar::t_type clk_type;
				if (obsdata.gsys() == GPS)       clk_type = t_gpar::G15CLK;
				else if (obsdata.gsys() == GAL)  clk_type = t_gpar::E15CLK;
				else if (obsdata.gsys() == BDS)  clk_type = t_gpar::C15CLK;
				else if (obsdata.gsys() == QZS)  clk_type = t_gpar::J15CLK;

                int idx = params.getParam(obsdata.site(), clk_type, "");
                if (idx <0)
				{
                    return false;
                }

                //update_value
                int idx_clk12 = params.getParam(obsdata.site(), t_gpar::CLK, "");
				
				params[idx].value( params[idx_clk12].value());

                coef_IF.push_back(make_pair(idx + 1, 1.0 - obsdata.drate()));
            }

			// jdhuang : for debug
			//string tmp,name;
			//string time(epoch.str_ymdhms());
			//char tmp_wgt[18];  sprintf(tmp_wgt, "%16.8f", P_IF);  string wgt(tmp_wgt);
			//char tmp_rsi[18];  sprintf(tmp_rsi, "%16.8f", l_IF);  string rsi(tmp_rsi);
			//if (gobs1.is_code()) name = "OBSP  ";
			//else name = "OBSL  ";
			//tmp = name + time +
			//	"Range rec : " + obsdata.site() +
			//	" satP : " + obsdata.sat() +
			//	" wgtP : " + wgt +
			//	" rsiP : " + rsi; //+ 
#ifdef DEBUG
			cout << tmp << endl;
#endif
			t_gobscombtype type(gobs1, OBSCOMBIN::IONO_FREE);
            equ_IF.add_equ(coef_IF,P_IF,l_IF,obsdata.site(),obsdata.sat(), type, false);
        }

        if (dynamic_cast<t_glsqEquationMatrix*>(&result))
        {
            t_glsqEquationMatrix *lsq_result = dynamic_cast<t_glsqEquationMatrix *>(&result);
            lsq_result->add_equ(equ_IF);
        }
        else
        {
            result = result + equ_IF;
        }
	
        return true;
    }



    t_gcombIF13::t_gcombIF13(t_gsetbase *setting,t_glog* log,shared_ptr<t_gbiasmodel> bias_model,t_gallproc* data):
        t_gcombmodel(setting,log,bias_model,data),
        _cmb_IF(setting,log,bias_model,data)
    {
    }

    t_gcombIF13::~t_gcombIF13()
    {
    }


    bool t_gcombIF13::cmb_equ(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,t_gbaseEquation& result)
    {
		obsdata.apply_bias(*_gallbias);
        // check Obs type
        GOBSBAND b1 = _band_index[obsdata.gsys()][FREQ_1];
        GOBSBAND b2 = _band_index[obsdata.gsys()][FREQ_2];
        GOBSBAND b3 = _band_index[obsdata.gsys()][FREQ_3];

        if (!_cmb_IF.cmb_equ_IF(epoch,params,obsdata,b1,b2,result))
        {
            return false;
        }

        _cmb_IF.cmb_equ_IF(epoch, params, obsdata, b1, b3, result);

        return true;
    }

    t_gcombALL::t_gcombALL(t_gsetbase *setting,t_glog* log,shared_ptr<t_gbiasmodel> bias_model,t_gallproc* data):
        t_gcombmodel(setting,log,bias_model,data)
    {
		ambtype_list[FREQ_1] = t_gpar::AMB_L1;
		ambtype_list[FREQ_2] = t_gpar::AMB_L2;
		ambtype_list[FREQ_3] = t_gpar::AMB_L3;
		ambtype_list[FREQ_4] = t_gpar::AMB_L4;
		ambtype_list[FREQ_5] = t_gpar::AMB_L5;
    }

    t_gcombALL::~t_gcombALL()
    {
    }


    bool t_gcombALL::cmb_equ(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,t_gbaseEquation& result)
    {
		obsdata.apply_bias(*_gallbias);
        // check band_list
        map<FREQ_SEQ, GOBSBAND>& crt_bands=_band_index[obsdata.gsys()];
      

		string grec = obsdata.site();
		string gsat = obsdata.sat();
		GSYS   gsys = obsdata.gsys();

        if (crt_bands.empty())
		{
			_log->logInfo("t_gcombALL", "cmb_equ", "crt bands is empty for ==> " + t_gsys::gsys2str(gsys));
			return false;
        }

        for (const auto &iter : crt_bands)
        {
            t_glsqEquationMatrix equ_ALL;
			vector<pair<int, double> > coefP;
			vector<pair<int, double> > coefL;

			double omcP = 0.0;
			double omcL = 0.0;
			double wgtP = 0.0;
			double wgtL = 0.0;
			GOBSBAND band = iter.second;
			t_gobs  obsP = obsdata.select_range(band, true);
			t_gobs  obsL = obsdata.select_phase(band, true);

			auto freq = iter.first;
			// check Obs Valid 
			if (obsP.type() != TYPE_C && obsP.type() != TYPE_P)
			{
				_log->logInfo("t_gcombALL", "cmb_equ", "check your obs file, we have no range obs for " + gobsband2str(band));
				if (freq == FREQ_1) return false;
				else continue;
				
			}

			if (obsL.type() != TYPE_L)
			{
				_log->logInfo("t_gcombALL", "cmb_equ", "check your obs file, we have no phase obs for " + gobsband2str(band));
				if (freq == FREQ_1) return false;
				else continue;
			}

			// jdhuang add for freq3 check
		
			if (_frequency == 1 && freq >= FREQ_2) continue; // 单频模式，处理 freq1
			if (_frequency == 2 && freq >= FREQ_3) continue; // 双频模式，只处理freq 1 和 freq 2

			// jdhuang
			if (freq == FREQ_1 && !obsdata.tb12() && !obsdata.tb13())
			{
				_log->logInfo("t_gprecisemodel", "_cmb_equ_ALL", "bad arc in log tb for freq 1" + grec + " + " + gsat);
				continue;
			}

			if (!_update_amb_lite && _frequency >= 3 && ((freq == FREQ_2 && !obsdata.tb12())||( freq == FREQ_3 && !obsdata.tb13())||(freq == FREQ_4 && !obsdata.tb14())||(freq == FREQ_5 && !obsdata.tb15())))
			{
				_log->logInfo("t_gprecisemodel", "_cmb_equ_ALL", "bad arc in log tb for freq X" + grec + " + " + gsat);
				continue;
			}

			t_gbaseEquation tempP;
			if (!_bias_model->cmb_equ(epoch, params, obsdata, obsP, tempP))
			{
				// jdhaung
				continue;
			
				//return false;
			}

			t_gbaseEquation tempL;
			if (!_bias_model->cmb_equ(epoch, params, obsdata, obsL, tempL))
			{
				continue;
				
				//return false;
			}

			//cout << "Test Xgg used   " + gsat + "  " <<  freq << "  " + epoch.str_ymdhms() + "  " << gobs2str(obsP.gobs()) << "  " + gobs2str(obsL.gobs()) << endl;

			// add AMB
			int idx = params.getParam(grec, ambtype_list[_freq_index[gsys][band]], gsat);
			if(idx >=0){//lvhb modified for npp testing,20200917
				if (idx < 0)
				{
					_log->logError("t_gcombALL", "cmb_equ", "check your log file, the amb idx < 0 for :" + grec + "_" + gsat);
					return false;
				}
		

				// if amb is new  then init value with Obs - Modelobs
				if ( double_eq(params[idx].value(), 0.0) || params[idx].beg == epoch)
				{
					// jdhuang change for ambfix
					double obs_P = obsdata.obs_C(obsP);
					double obs_L = obsdata.obs_L(obsL);

					//auto cmb_type = OBSCOMBIN::RAW_ALL;
					//double CmbIntraFreq = 0.0;
					//CmbIntraFreq = _gprecisemodel->intraFreqDcb("CAS", _gallbias, obsdata, &obsP);
					//CmbIntraFreq += _gprecisemodel->cmbInterFreqDcb("CAS", cmb_type, _gallbias, obsdata, &obsP, &obsP);
					//obs_P = obs_P - CmbIntraFreq;
		
			
					params[idx].value(obs_L - obs_P);
  
				}
			
				tempL.B[0].push_back(make_pair(idx+1,1.0));
				tempL.l[0] -= params[idx].value();
				//ifcb correction
				if (_ifcb_model == IFCB_MODEL::COR && _frequency >= 3 && freq == FREQ_3) 
				{
					apply_IFCB_GPS(obsdata, _gprecisemodel->_gifcb, tempL.l[0], OBSCOMBIN::RAW_ALL);
				}
			}
			
			// jdhuang
			// ifb corrections
			//	if (_frequency == 3)
			//	{
				//auto cmb_type = OBSCOMBIN::RAW_ALL;
				//double CmbInterFreq = _gprecisemodel->cmbInterFreqDcb("CAS", cmb_type, _gallbias, obsdata, &obsP, &obsP);
				//tempP.l[0] -= CmbInterFreq;
			//	} 
			// jdhuang : for debug
	  //		string tmp;
	  //		string time(epoch.str_ymdhms());
	  //		char tmp_wgt[18];  sprintf(tmp_wgt, "%16.8f", tempP.P[0]);  string wgt(tmp_wgt);
	  //		char tmp_rsi[18];  sprintf(tmp_rsi, "%16.8f", tempP.l[0]);  string rsi(tmp_rsi);
	  //		tmp = time +
	  //			"Range rec : " + obsdata.site() +
	  //			" wgtP : " + wgt +
	  //			" rsiP : " + rsi; //+ 
			//cout << "OBSP  " + tmp << endl;
	//		cout << "obs_P    " << obsdata.sat() <<fixed <<setw(12)<< setprecision(3) << obsdata.obs_C(obsP) << endl;
	//		cout << "obs_L    " << obsdata.sat() << fixed << setw(12) << setprecision(3) << obsdata.obs_L(obsL) << endl;
	//// jdhuang : for debug
	  // 		sprintf(tmp_wgt, "%16.8f", tempL.P[0]);  string wgt_L(tmp_wgt);
	  // 		sprintf(tmp_rsi, "%16.8f", tempL.l[0]);  string rsi_L(tmp_rsi);
	  // 		tmp = time +
	  // 			"Range rec : " + obsdata.site() +
	  // 			" satL : " + obsdata.sat() +
	  // 			" wgtL : " + wgt_L +
	  // 			" rsiL : " + rsi_L;// +
			//cout << "OBSL  " + tmp << endl;

			t_gobscombtype typeP(obsP, OBSCOMBIN::RAW_ALL);
			t_gobscombtype typeL(obsL, OBSCOMBIN::RAW_ALL);
            equ_ALL.add_equ(tempP.B[0], tempP.P[0], tempP.l[0], obsdata.site(), obsdata.sat(), typeP,false);
            equ_ALL.add_equ(tempL.B[0], tempL.P[0], tempL.l[0], obsdata.site(), obsdata.sat(), typeL,false);


            // add result
            if (dynamic_cast<t_glsqEquationMatrix*>(&result))
            {
                t_glsqEquationMatrix *lsq_reuslt = dynamic_cast<t_glsqEquationMatrix *>(&result);
                lsq_reuslt->add_equ(equ_ALL);
            }
            else
            {
                result = result+ equ_ALL;
            }
        }

        return true;
    }
	    t_gcombWL::t_gcombWL(t_gsetbase* setting, t_glog* log, shared_ptr<t_gbiasmodel> bias_model, t_gallproc* data) :
        t_gcombmodel(setting, log, bias_model, data)
    {
    }

	t_gcombWL::~t_gcombWL()
	{
	}

    bool t_gcombWL::cmb_equ(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, t_gbaseEquation& result)
    {
		obsdata.apply_bias(*_gallbias);
        // check Obs type
        GOBSBAND b1 = _band_index[obsdata.gsys()][FREQ_1];
        GOBSBAND b2 = _band_index[obsdata.gsys()][FREQ_2];

        return this->cmb_equ_WL(epoch, params, obsdata, b1, b2, result);
    }


    bool t_gcombWL::cmb_equ_WL(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, GOBSBAND b1, GOBSBAND b2, t_gbaseEquation& result)
    {
        // only for freq1 and freq2
        t_gobs gobsP1 = obsdata.select_range(b1);
        t_gobs gobsP2 = obsdata.select_range(b2);
        t_gobs gobsL1 = obsdata.select_phase(b1);
        t_gobs gobsL2 = obsdata.select_phase(b2);

        vector<pair<t_gobs, t_gobs> > type_list;
        type_list.push_back(make_pair(gobsP1, gobsP2));
        type_list.push_back(make_pair(gobsL1, gobsL2));

        // WL coef
        double coef1, coef2;
        
        t_glsqEquationMatrix equ_WL;

        // cycle [Code and Phase]
        for (auto item : type_list)
        {
            t_gobs gobs1 = item.first;
            t_gobs gobs2 = item.second;

            if (gobs1.is_phase())
            {
                obsdata.coef_widelane(b1, coef1, b2, coef2);
            }
            else
            {
                obsdata.coef_narrlane(b1, coef1, b2, coef2);
            }

            //combine P1 and P2
            t_gbaseEquation temp_equ;
            if (!_bias_model->cmb_equ(epoch, params, obsdata, gobs1, temp_equ))
            {
                return false;
            }
            if (!_bias_model->cmb_equ(epoch, params, obsdata, gobs2, temp_equ))
            {
                return false;
            }

            vector<pair<int, double> > coef_WL;
            double P_WL = 0.0, l_WL = 0.0;
            if (temp_equ.B[0].size() != temp_equ.B[1].size()) {
                throw logic_error("coeff size is not equal in f1 and f2");
            }
            for (int i = 0; i < temp_equ.B[0].size(); i++)
            {
                if (temp_equ.B[0][i].first != temp_equ.B[1][i].first) {
                    throw logic_error("coeff par is not the same in f1 and f2");
                }
                // combine coeff
                coef_WL.push_back(make_pair(temp_equ.B[0][i].first, coef1 * temp_equ.B[0][i].second + coef2 * temp_equ.B[1][i].second));
            }

			// combine P
			P_WL = sqrt(pow(coef1, 2) * pow(temp_equ.P[0], 2) + pow(coef2, 2) * pow(temp_equ.P[1], 2));
			// combine l
			l_WL = coef1 * temp_equ.l[0] + coef2 * temp_equ.l[1];

            // CORRECT AMB
            if (gobs1.is_phase())
            {
                t_gpar::t_type amb_type = t_gpar::AMB_WL;

                int idx = params.getParam(obsdata.site(), amb_type, obsdata.sat());

                if (idx < 0)
                {
                    return false;
                }

                // if amb is new  then init value with Obs - Modelobs
                if (double_eq(params[idx].value(), 0.0) || params[idx].beg == epoch)
                {
                    double MW_value = -obsdata.MW(gobs1, gobs2);
                    //correct_DCB("CAS", _gallbias, obsdata, MW_value, OBSCOMBIN::MW_COMBIN, &gobs1, &gobs2);  // glfeng for Normalize
                    params[idx].value(MW_value);
                }

                // update B l
                coef_WL.push_back(make_pair(idx + 1, 1.0));
                l_WL -= params[idx].value();
            }

            // jdhuang : for debug
            if (!gobs1.is_phase())
            {
                string tmp;
                string time(epoch.str_ymdhms());
                char tmp_wgt[18];  sprintf(tmp_wgt, "%16.8f", P_WL);  string wgt(tmp_wgt);
                char tmp_rsi[18];  sprintf(tmp_rsi, "%16.8f", l_WL);  string rsi(tmp_rsi);
                tmp = time +
                    "Range rec : " +
                    " satP : " + obsdata.sat() +
                    " wgtP : " + wgt +
                    " rsiP : " + rsi; //+ 
                //cout << "OBSP  " + tmp << endl;
            }
            else
            {
                string tmp;
                string time(epoch.str_ymdhms());
                char tmp_wgt[18];  sprintf(tmp_wgt, "%16.8f", P_WL);  string wgt(tmp_wgt);
                char tmp_rsi[18];  sprintf(tmp_rsi, "%16.8f", l_WL);  string rsi(tmp_rsi);
                tmp = time +
                    "Range rec : " +
                    " satL : " + obsdata.sat() +
                    " wgtL : " + wgt +
                    " rsiL : " + rsi; //+ 
               // cout << "OBSL  " + tmp << endl;

            }
			t_gobscombtype type(gobs1, OBSCOMBIN::WL_COMBIN);
            equ_WL.add_equ(coef_WL, P_WL, l_WL, obsdata.site(), obsdata.sat(),type, false);
        }

        if (dynamic_cast<t_glsqEquationMatrix*>(&result))
        {
            t_glsqEquationMatrix* lsq_result = dynamic_cast<t_glsqEquationMatrix*>(&result);
            lsq_result->add_equ(equ_WL);
        }
        else
        {
            result = result + equ_WL;
        }
        return true;
    }

	t_gcombDD::t_gcombDD(t_gsetbase *setting, t_glog* log, shared_ptr<t_gbiasmodel> bias_model, t_gallproc* data) :
		t_gcombmodel(setting, log, bias_model, data),
		t_gcombALL(setting, log, bias_model, data),
		t_gcombIF(setting, log, bias_model, data)
	{
		_data_base = nullptr;
		_site = "";
		_site_base = "";
	}

	t_gcombDD::~t_gcombDD()
	{
	}

	void t_gcombDD::set_observ(OBSCOMBIN observ)
	{
		_observ = observ;
	}

	void t_gcombDD::set_base_data(vector<t_gsatdata>* data_base)
	{
		_data_base = data_base;
	}

	void t_gcombDD::set_site(string site, string site_base)
	{
		_site = site;
		_site_base = site_base;
	}

	void t_gcombDD::set_rec_info(t_gtriple xyz_base, double clk_rover, double clk_base)
	{
		_crd_base = xyz_base;
		_clk_rover = clk_rover;
		_clk_base = clk_base;
	}

	bool t_gcombDD::cmb_equ(t_gtime & epoch, t_gallpar & params, t_gsatdata & obsdata, t_gbaseEquation & result)
	{
		t_gsatdata obsdata_other;
		if (_data_base == nullptr) return false;
		for (int i = 0; i < _data_base->size(); i++)
		{
			if (obsdata.sat() == (*_data_base)[i].sat()) obsdata_other = (*_data_base)[i];
		}
		if (obsdata_other.sat() == "") return false;
		//obsdata.apply_dcb(*_gallbias);
		//obsdata_other.apply_dcb(*_gallbias);
		int npar_orig = params.parNumber();
		t_gallpar params_temp;
		_temp_params(params, params_temp);
		if (_observ == IONO_FREE)
		{
			// check Obs type
			GOBSBAND b1 = _band_index[obsdata.gsys()][FREQ_1];
			GOBSBAND b2 = _band_index[obsdata.gsys()][FREQ_2];

			t_gobs gobsP1 = obsdata.select_range(b1);
			t_gobs gobsP2 = obsdata.select_range(b2);
			t_gobs gobsL1 = obsdata.select_phase(b1);
			t_gobs gobsL2 = obsdata.select_phase(b2);

			t_gobs gobsP1_other = obsdata_other.select_range(b1);
			t_gobs gobsP2_other = obsdata_other.select_range(b2);
			t_gobs gobsL1_other = obsdata_other.select_phase(b1);
			t_gobs gobsL2_other = obsdata_other.select_phase(b2);

			vector<pair<t_gobs, t_gobs> > type_list;
			type_list.push_back(make_pair(gobsP1, gobsP2));
			type_list.push_back(make_pair(gobsL1, gobsL2));

			vector<pair<t_gobs, t_gobs> > type_list_other;
			type_list_other.push_back(make_pair(gobsP1_other, gobsP2_other));
			type_list_other.push_back(make_pair(gobsL1_other, gobsL2_other));

			// IF coef
			double coef1, coef2;
			obsdata.coef_ionofree(b1, coef1, b2, coef2);

			t_glsqEquationMatrix equ_IF;

			// meter [Code and Phase]
			for (int i=0;i<type_list.size();i++)
			{
				t_gobs gobs1 = type_list[i].first;
				t_gobs gobs2 = type_list[i].second;
				t_gobs gobs1_other = type_list_other[i].first;
				t_gobs gobs2_other = type_list_other[i].second;

				//combine P1 and P2
				t_gbaseEquation temp_equ;
				t_gtime crt = obsdata.epoch();
				if (!_bias_model->cmb_equ(crt, params_temp, obsdata, gobs1, temp_equ))
				{
					/*return false;*/
					continue;//Lvhb changed in 202012
				}
				if (!_bias_model->cmb_equ(crt, params_temp, obsdata, gobs2, temp_equ))
				{
					//return false;
					continue;//Lvhb changed in 202012
				}
				t_gbaseEquation temp_equ_other;
				t_gtime crt_other = obsdata_other.epoch();
				if (!_bias_model->cmb_equ(crt_other, params_temp, obsdata_other, gobs1_other, temp_equ_other))
				{
					//return false;
					continue;//Lvhb changed in 202012
				}
				if (!_bias_model->cmb_equ(crt_other, params_temp, obsdata_other, gobs2_other, temp_equ_other))
				{
					//return false;
					continue;//Lvhb changed in 202012
				}
				vector<pair<int, double> > coef_IF;
				double P_IF = 0.0, l_IF = 0.0;

				if (temp_equ.B[0].size() != temp_equ.B[1].size()) {
					throw logic_error("coeff size is not equal in f1 and f2");
				}
				for (int i = 0; i < temp_equ.B[0].size(); i++)
				{
					if (temp_equ.B[0][i].first != temp_equ.B[1][i].first) {
						throw logic_error("coeff par is not the same in f1 and f2");
					}
					// combine coeff
					if (temp_equ.B[0][i].first > npar_orig) continue;
					coef_IF.push_back(make_pair(temp_equ.B[0][i].first, coef1 * temp_equ.B[0][i].second + coef2 * temp_equ.B[1][i].second));
				}

				if (temp_equ_other.B[0].size() != temp_equ_other.B[1].size()) {
					throw logic_error("coeff size is not equal in f1 and f2");
				}
				for (int i = 0; i < temp_equ_other.B[0].size(); i++)
				{
					if (temp_equ_other.B[0][i].first != temp_equ_other.B[1][i].first) {
						throw logic_error("coeff par is not the same in f1 and f2");
					}
					// combine coeff
					if (temp_equ_other.B[0][i].first > npar_orig) continue;
					coef_IF.push_back(make_pair(temp_equ_other.B[0][i].first, -(coef1 * temp_equ_other.B[0][i].second + coef2 * temp_equ_other.B[1][i].second)));
				}

				// combine P
				P_IF = 1.0 / (pow(coef1, 2) / temp_equ.P[0] + pow(coef2, 2) / temp_equ.P[1]);
				// combine l
				l_IF = coef1 * temp_equ.l[0] + coef2 * temp_equ.l[1]
					- (coef1 * temp_equ_other.l[0] + coef2 * temp_equ_other.l[1]);


				// CORRECT AMB
				if (gobs1.is_phase())
				{
					t_gpar::t_type amb_type = t_gpar::AMB_IF;
					if (_freq_index[obsdata.gsys()][b2] == FREQ_3) {
						amb_type == t_gpar::AMB13_IF;
					}

					int idx = params_temp.getParam(obsdata.site(), amb_type, obsdata.sat());

					if (idx < 0)
					{
						return false;
					}

					// if amb is new  then init value with Obs - Modelobs
					if (double_eq(params_temp[idx].value(), 0.0) || params_temp[idx].beg == epoch)
					{
						double obs_P3 = obsdata.P3(gobsP1, gobsP2);
						double obs_L3 = obsdata.L3(gobsL1, gobsL2);
						////////////apply_DCB(_gallbias, obsdata, obs_P3, OBSCOMBIN::IONO_FREE, &gobs1, &gobs2);
						//pars[idx].value(obs - ModelObs);
						params_temp[idx].value(obs_L3 - obs_P3);
					}

					// update B l
					coef_IF.push_back(make_pair(idx + 1, 1.0));
					l_IF -= params_temp[idx].value();
				}

				// CORRECT CLK
				//remove 12 clk and isb
				if (_freq_index[obsdata.gsys()][b2] >= FREQ_3)
				{
					for (auto iter = coef_IF.begin(); iter != coef_IF.end();)
					{
						t_gpar::t_type type = params_temp[iter->first - 1].parType;
						if (type == t_gpar::CLK ||
							type == t_gpar::GAL_ISB ||
							type == t_gpar::BDS_ISB ||
							type == t_gpar::QZS_ISB)
						{
							iter = coef_IF.erase(iter);
							continue;
						}
						iter++;
					}
				}
				// add 13 CLK
				if (_freq_index[obsdata.gsys()][b2] == FREQ_3)
				{
					t_gpar::t_type clk_type;
					if (obsdata.gsys() == GPS)
						clk_type = t_gpar::G15CLK;
					else if (obsdata.gsys() == GAL)
						clk_type = t_gpar::E15CLK;
					else if (obsdata.gsys() == BDS)
						clk_type = t_gpar::C15CLK;
					else if (obsdata.gsys() == QZS)
						clk_type = t_gpar::J15CLK;

					int idx = params_temp.getParam(obsdata.site(), clk_type, "");
					if (idx < 0) {
						return false;
					}
					//update_value
					int idx_clk12 = params_temp.getParam(obsdata.site(), t_gpar::CLK, obsdata.sat());
					params_temp[idx] = params_temp[idx_clk12];

					coef_IF.push_back(make_pair(idx + 1, 1.0 - obsdata.drate()));
				}

				// add 14 15 ... clk
				// ....
				t_gobscombtype type(gobs1, OBSCOMBIN::IONO_FREE);
				equ_IF.add_equ(coef_IF, P_IF, l_IF, obsdata.site(), obsdata.sat(), type, false);
			}

			if (dynamic_cast<t_glsqEquationMatrix*>(&result))
			{
				t_glsqEquationMatrix *lsq_result = dynamic_cast<t_glsqEquationMatrix *>(&result);
				lsq_result->add_equ(equ_IF);
			}
			else
			{
				result = result + equ_IF;
			}
		}
		else if (_observ == RAW_ALL || _observ == RAW_SINGLE || _observ==RAW_MIX)
		{
			// check band_list
			map<FREQ_SEQ, GOBSBAND> crt_bands = _band_index[obsdata.gsys()];
			// added RAW_MIX wh
			if ( ( _observ == RAW_SINGLE || _observ == RAW_MIX) && crt_bands.size() > 1)
			{
				auto first_band = (*crt_bands.begin());
				crt_bands.clear();
				crt_bands.insert(first_band);
			}
			map<FREQ_SEQ, t_gpar::t_type> ambtype_list = {
				{FREQ_1, t_gpar::AMB_L1},
				{FREQ_2, t_gpar::AMB_L2},
				{FREQ_3, t_gpar::AMB_L3},
				{FREQ_4, t_gpar::AMB_L4},
				{FREQ_5, t_gpar::AMB_L5} };

			if (crt_bands.empty()) {
				return false;
			}
			t_glsqEquationMatrix equ_ALL;
			t_gbaseEquation tempP, tempL;
			vector<int> index(10, 0);//added by lvhb in 20201211,10=5freq*2site;

			for (int isite = 0; isite < 2; isite++)
			{
				t_gsatdata *satdata_ptr;
				if (isite == 0) satdata_ptr = &obsdata_other;
				else satdata_ptr = &obsdata;

				for (const auto &iter : crt_bands)
				{
					if (iter.first > _frequency) continue;//modified by lvhb in 20201211

					GOBSBAND band = iter.second;
					t_gobs  obsP = satdata_ptr->select_range(band);
					t_gobs  obsL = satdata_ptr->select_phase(band);

					// check Obs Valid 
					// for P code modified by zhshen 
					if (obsP.type() != TYPE_C && obsP.type() != TYPE_P)
					{
						_log->logInfo("t_gprecisemodel", "_cmb_equ_ALL", "check your obs file, we have no range obs for " + gobsband2str(band));
						continue;//Lvhb changed in 202012
					}

					if (obsL.type() != TYPE_L)
					{
						_log->logInfo("t_gprecisemodel", "_cmb_equ_ALL", "check your obs file, we have no phase obs for " + gobsband2str(band));
						continue;//Lvhb changed in 202012
					}

					// jdhuang add for freq3 check
					auto freq = iter.first;
					t_gtime crt = satdata_ptr->epoch();
					if (!_bias_model->cmb_equ(crt, params_temp, *satdata_ptr, obsP, tempP))
					{
						continue;//Lvhb changed in 202012
					}

					if (!_bias_model->cmb_equ(crt, params_temp, *satdata_ptr, obsL, tempL))
					{
						continue;//Lvhb changed in 202012
					}
					index[isite * crt_bands.size() + freq - 1] = tempL.l.size();
					// add AMB for ROVER
					if (satdata_ptr->site() == _site)
					{
						int idx = params_temp.getParam(satdata_ptr->site(), ambtype_list[_freq_index[satdata_ptr->gsys()][band]], satdata_ptr->sat());

						if (idx < 0)
						{
							return false;
						}
						// if amb is new  then init value with Obs - Modelobs
						if (double_eq(params_temp[idx].value(), 0.0) || params_temp[idx].beg == epoch)
						{
							// jdhuang change for ambfix
							//pars[idx].value(obs - ModelObs);
							double obs_P = satdata_ptr->obs_C(obsP);
							double obs_L = satdata_ptr->obs_L(obsL);
							//apply_DCB(_gallbias, obsdata, obs_P, OBSCOMBIN::RAW_ALL, &obsP);

							// jdhuang : without sion pars, a test
							params_temp[idx].value(obs_L - obs_P);

						}
						tempL.B.back().push_back(make_pair(idx + 1, 1.0));
						tempL.l.back() -= params_temp[idx].value();
					}
				}
			}

			// make differenced equations

			int freq_count = 0;
			for (const auto &iter : crt_bands)
			{
				if(iter.first > _frequency) continue;//modified by lvhb in 20201211

				freq_count++;

				GOBSBAND band = iter.second;
				t_gobs  obsP = obsdata.select_range(band);
				t_gobs  obsL = obsdata.select_phase(band);

				vector<pair<int, double>> B_P, B_L;
				double P_P, P_L, l_P, l_L;
				//int ibase = freq_count - 1;
				//int irover = crt_bands.size()  + freq_count - 1;
				//if(!index[ibase]||!index[irover]) continue;
				int ibase = index[freq_count - 1];
				int irover = index[crt_bands.size() + freq_count - 1];
				if(!ibase || !irover) continue;
				ibase -= 1; irover -= 1;

				for (int i = 0; i < tempP.B[irover].size(); i++)
				{
					if (tempP.B[irover][i].first > npar_orig) continue;
					B_P.push_back(tempP.B[irover][i]);
				}
				for (int i = 0; i < tempP.B[ibase].size(); i++)
				{
					if (tempP.B[ibase][i].first > npar_orig) continue;
					B_P.push_back(make_pair(tempP.B[ibase][i].first, -tempP.B[ibase][i].second));
				}
				for (int i = 0; i < tempL.B[irover].size(); i++)
				{
					if (tempL.B[irover][i].first > npar_orig) continue;
					B_L.push_back(tempL.B[irover][i]);
				}
				for (int i = 0; i < tempL.B[ibase].size(); i++)
				{
					if (tempL.B[ibase][i].first > npar_orig) continue;
					B_L.push_back(make_pair(tempL.B[ibase][i].first, -tempL.B[ibase][i].second));
				}
				
				P_P = 1 / (1 / tempP.P[irover] + 1 / tempP.P[ibase]); l_P = tempP.l[irover] - tempP.l[ibase];
				P_L = 1 / (1 / tempL.P[irover] + 1 / tempL.P[ibase]); l_L = tempL.l[irover] - tempL.l[ibase];

				t_gobscombtype typeP(obsP, OBSCOMBIN::RAW_ALL);
				t_gobscombtype typeL(obsL, OBSCOMBIN::RAW_ALL);
				equ_ALL.add_equ(B_P, P_P, l_P, obsdata.site(), obsdata.sat(), typeP, false);
				equ_ALL.add_equ(B_L, P_L, l_L, obsdata.site(), obsdata.sat(), typeL, false);

			}

			// add result
			if (dynamic_cast<t_glsqEquationMatrix*>(&result))
			{
				t_glsqEquationMatrix *lsq_reuslt = dynamic_cast<t_glsqEquationMatrix *>(&result);
				lsq_reuslt->add_equ(equ_ALL);
			}
			else
			{
				result = result + equ_ALL;
			}
		}
		for (int ipar = 0; ipar < params.parNumber(); ipar++)
			params[ipar] = params_temp[ipar];
		return true;
	}

	bool t_gcombDD::cmb_equ_test(t_gtime & epoch, t_gallpar & params, t_gsatdata & obsdata, t_gbaseEquation & result)
	{
        obsdata.apply_bias(*_gallbias);
		t_gsatdata obsdata_other;
		if (_data_base == nullptr) return false;
		for (int i = 0; i < _data_base->size(); i++)
		{
			if (obsdata.sat() == (*_data_base)[i].sat()) obsdata_other = (*_data_base)[i];
		}
		if (obsdata_other.sat() == "") return false;

		t_glsqEquationMatrix res_base, res_rove;
		// base
		int npar_orig = params.parNumber();
		t_gallpar params_temp;
		if (!_temp_params(params, params_temp)) return false;

		if (_observ == IONO_FREE)
		{
			if (!t_gcombIF::cmb_equ(epoch, params_temp, obsdata, res_rove)
				|| !t_gcombIF::cmb_equ(epoch, params_temp, obsdata_other, res_base))
				return false;
		}
		else if (_observ = RAW_ALL)
		{
			if (!t_gcombALL::cmb_equ(epoch, params_temp, obsdata, res_rove)
				|| !t_gcombALL::cmb_equ(epoch, params_temp, obsdata_other, res_base))
				return false;
		}
		else return false;

		int base_obs_num = res_base.B.size();
		int base_par_num = res_base.B[0].size();
		int rove_obs_num = res_rove.B.size();
		int rove_par_num = res_rove.B[0].size();
		if (base_obs_num != rove_obs_num) {
			_log->logInfo("t_gcombDD", "cmb_equ", "rover and base observation dimensions are not equal!");
			return false;
		}
		t_glsqEquationMatrix equ_SD;

		for (int i = 0; i < base_obs_num ; i++)
		{
			vector<pair<int, double>> coef_SD;
			for (int j = 0; j < res_rove.B[i].size(); j++)
			{
				if (res_rove.B[i][j].first > npar_orig) continue;
				coef_SD.push_back(make_pair(res_rove.B[i][j].first, -res_rove.B[i][j].second));
			}
			for (int j = 0; j < res_base.B[i].size(); j++)
			{
				if (res_base.B[i][j].first > npar_orig) continue;
				coef_SD.push_back(make_pair(res_base.B[i][j].first, res_base.B[i][j].second));
			}
			double P_SD = 1 / (1.0 / res_rove.P[i] + 1.0 / res_base.P[i]);
			double l_SD = res_rove.l[i] - res_base.l[i];
			equ_SD.add_equ(coef_SD, P_SD, l_SD, res_rove.get_sitename(i), res_rove.get_satname(i), res_rove.get_obscombtype(i),false);
		}

		if (dynamic_cast<t_glsqEquationMatrix*>(&result))
		{
			t_glsqEquationMatrix *lsq_result = dynamic_cast<t_glsqEquationMatrix *>(&result);
			lsq_result->add_equ(equ_SD);
		}
		else
		{
			result = result + equ_SD;
		}

		for (int ipar = 0; ipar < params.parNumber(); ipar++)
			params[ipar] = params_temp[ipar];
		return true;
	}

	bool t_gcombDD::_temp_params(t_gallpar& params, t_gallpar& params_temp)
	{
		params_temp = params;
		t_gpar par_x_base; par_x_base.site = _site_base; par_x_base.parType = t_gpar::CRD_X; par_x_base.value(_crd_base[0]);
		t_gpar par_y_base; par_y_base.site = _site_base; par_y_base.parType = t_gpar::CRD_Y; par_y_base.value(_crd_base[1]);
		t_gpar par_z_base; par_z_base.site = _site_base; par_z_base.parType = t_gpar::CRD_Z; par_z_base.value(_crd_base[2]);
		t_gpar par_clk_rover; par_clk_rover.site = _site; par_clk_rover.parType = t_gpar::CLK; par_clk_rover.value(_clk_rover);
		t_gpar par_clk_base; par_clk_base.site = _site_base; par_clk_base.parType = t_gpar::CLK; par_clk_base.value(_clk_base);
		params_temp.addParam(par_x_base);
		params_temp.addParam(par_y_base);
		params_temp.addParam(par_z_base);
		params_temp.addParam(par_clk_rover);
		params_temp.addParam(par_clk_base);
		params_temp.reIndex();
		return true;
	}


	t_gcombMIX::t_gcombMIX(t_gsetbase *setting, t_glog* log, shared_ptr<t_gbiasmodel> bias_model, t_gallproc* data):
		t_gcombmodel(setting, log, bias_model, data),
		t_gcombIF(setting, log, bias_model, data)
	{
		_observ = dynamic_cast<t_gsetproc*>(setting)->obs_combin();
		_gdata_ionex = dynamic_cast<t_gionex*>((*data)[t_gdata::IONEX]);    //< IONEX Grid data, glfeng
	}

	t_gcombMIX::~t_gcombMIX()
	{

	}

	bool t_gcombMIX::cmb_equ(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, t_gbaseEquation& result)
	{
		obsdata.apply_bias(*_gallbias);
		if (_observ == IF_P1)
		{
			return this->_cmb_equ_IF_P1(epoch, params, obsdata, result);
		}
		return false;
	}


	bool t_gcombMIX::_cmb_equ_IF_P1(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, t_gbaseEquation& result)
	{
		GSYS sys = obsdata.gsys();
		set<GOBSBAND> bands = obsdata.band_avail();
		GOBSBAND b1 = _band_index[sys][FREQ_1];
		GOBSBAND b2 = _band_index[sys][FREQ_2];
		
		// P1 单双频模式，只处理FREQ1 range
		bool isIF = false, isP1 = false;

		if (bands.find(b1) != bands.end())	
		{
			t_glsqEquationMatrix equ_ALL;
			t_gobs obsP = obsdata.select_range(b1, true);

			t_gbaseEquation tempP;
			if (_bias_model->cmb_equ(epoch, params, obsdata, obsP, tempP))  isP1 = true;

			if (isP1)
			{
				// ion_delay
				if (_gdata_ionex)
				{
					t_giono_tecgrid iono_tec;
					double sion_L1, sion_rms;
					t_gtriple trs_rec_crd = _gprecisemodel->_trs_rec_crd;
					if (iono_tec.getIonoDelay(_gdata_ionex, obsdata, epoch, trs_rec_crd, sion_L1, sion_rms))
					{
						double f1 = obsdata.frequency(t_gsys::band_priority(obsdata.gsys(), FREQ_1));
						double fk = obsdata.frequency(b1);
						double alfa = (f1 * f1) / (fk * fk);
						double iono_delay = alfa * sion_L1;
						tempP.l.back() -= iono_delay;
					}
				}

				t_gobscombtype type(obsP, OBSCOMBIN::IF_P1);
				equ_ALL.add_equ(tempP.B[0], tempP.P[0], tempP.l[0], obsdata.site(), obsdata.sat(), type, false);

				// add result
				if (dynamic_cast<t_glsqEquationMatrix*>(&result))
				{
					t_glsqEquationMatrix *lsq_reuslt = dynamic_cast<t_glsqEquationMatrix *>(&result);
					lsq_reuslt->add_equ(equ_ALL);
				}
				else
				{
					result = result + equ_ALL;
				}
			}
		}

		// IF type 
		if (bands.find(b1) != bands.end() && bands.find(b2) != bands.end())
		{
			isIF = this->cmb_equ_IF(epoch, params, obsdata, b1, b2, result);
		}

		if (isP1 || isIF) return true;
		else      return false;
	}
}