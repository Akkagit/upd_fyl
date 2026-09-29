/**
*
* The format :
* @verbatim
  post-processing mode:
  <turboedit lite_mode="false" >
    <amb_output  valid="true"  />     <!-- 事后的静态与动态：都会储存到cycle_slip，判断是否会输出到log -->
    <ephemeris   valid="true"  />     <!-- 是否是使用广播星历 -->
    <check_pc   pc_limit="250" valid="true"  />
    <check_mw   mw_limit="4"   valid="true"  />
    <check_gf   gf_limit="1"   gf_rms_limit="2"   valid="true" />
    <check_sf   sf_limit="1"   valid="false" />
    <check_gap    gap_limit="20"    valid="true" />                                             <!-- unit:epoch, 可用于实时解，也可用于事后解 -->
    <check_short  short_limit="10"  valid="true" />                                             <!-- unit:epoch, 用于事后解 -->
    <check_statistics  min_percent="60"  min_mean_nprn="4"  max_mean_namb="3"  valid="true" />  <!-- 用于事后静态：网解 -->
  </turboedit>

  real-time mode or lite-mode:
  <turboedit lite_mode="true" >
	<check_mw   valid="true"  />
	<check_gf   valid="true"  />
	<smooth_win  value="25"   />  smooth windows length
	<check_gap  gap_limit="1"  valid="true" />           
  </turboedit>

* @endverbatim
*
* @verbatim
	History
	 -1.0 Guolong Feng	 2020-03-23  creat the file.
	 -2.0 Guolong Feng   2020-05-09  add real-time mode, separate
  @endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file			gsetturboedit.h
* @brief		turboedit set 
* @author       Guolong Feng, Wuhan University
* @version		1.0.0
* @date			2020-03-23
*
*/

#ifndef GSETTURBOEDIT_H
#define GSETTURBOEDIT_H

#define XMLKEY_TURBOEDIT             "turboedit"

#include <map>
#include <string>
#include <iostream>

#include "gset/gsetbase.h"
#include "gutils/ExportBase.h"

using namespace std;
using namespace gnut;

namespace great
{


	class BASE_LIBRARY_EXPORT t_gsetturboedit : public virtual t_gsetbase
	{
	public:
		/** @brief default constructor. */
		t_gsetturboedit();

		/** @brief default destructor. */
		virtual ~t_gsetturboedit() {};

		bool liteMode();

		bool isAmbOutput();

		bool isEphemeris();

		bool checkPC(double& pc_limit);

		bool checkMW(double& mw_limit);

		bool checkGF(double& gf_limit, double& gf_rms_limit);

		bool checkSingleFreq(double& sf_limit);

		bool checkGap(int& gap_limit);

		int smoothWindows();

		bool checkShort(int& short_limit);

		bool checkStatistics(double& minPercent, int& minMeanNprn, int& maxMeanNamb);

		/** @brief settings help. */
		void help();

	protected:
		
		double _defaultPCLimit;
		double _defaultMWLimit;
		double _defaultGFLimit;
		double _defaultGFRmsLimit;
		double _defaultSingleFreqLimit;
		int    _defaultGapArcLimit;
		int    _defaultShortArcLimit;
		double _defaultMinPercent;
		int    _defaultMinMeanNprn;
		int    _defaultMaxMeanNamb;

	};
}
#endif // !SETTURBOEDIT_H
