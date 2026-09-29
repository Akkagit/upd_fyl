/**
*
* @verbatim
	The format of this block:
	<ambiguity>
		<dd_mode>  IF_CB_WN, RAW_CB_WN, RAW_SG  </dd_mode>
		<!-- <upd> WL </upd> -->                    <!-- EWL WL NL / NONE -->
		<is_ppprtk> YES </is_ppprtk>				<!-- YES / NO -->
		<fix_mode>  SEARCH </fix_mode>				<!-- ROUND/SEARCH/NO -->
		<part_fix>  YES  </part_fix>			    <!--  YES/NO -->
		<boot>  2.0  </boot>						<!-- Threshold of bootstrapping in amb fix -->
		<ratio>  2.0  </ratio>						<!-- Threshold in LAMBDA method -->
		<add_leo>  NO   </add_leo>
		<all_baselines>  NO  </all_baselines>
		<min_common_time>  900  </min_common_time>                 <!-- Minimum common time/seconds -->
		<baseline_length_limit>  3500  </baseline_length_limit>   <!-- Baseline length limit/km-->
		<widelane_interval> 30 </widelane_interval>               <!-- seconds -->
		<widelane_decision     maxdev = "0.275" maxsig = "0.10" alpha = "1000" />   <!-- deriation, sigma in WL-cycle -->
		<narrowlane_decision   maxdev = "0.375" maxsig = "0.10" alpha = "1000" />   <!-- deriation, sigma in NL-cycle -->
	</ambiguity>
	@endverbatim

* @verbatim
	History
		-1.0 glFeng   2019-02-27 creat the file.
		-1.1 BoWong   2019-04-10 Adding Doxygen Style Code Remarks
		-1.2 glFeng   2019-04-11 Add function part_ambfix
	@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		   gsetamb.h
* @brief	   control set from XML
*
* @author      glFeng, Wuhan University
* @version	   1.0.0
* @date		   2019-04-10
*
*/

#ifndef GSETAMB_H
#define GSETAMB_H

#define XMLKEY_AMBIGUITY "ambiguity"

#include "gdata/gupd.h"
#include "gset/gsetbase.h"
#include "gutils/ExportAMB.h"

using namespace std;
using namespace gnut;
namespace great
{
	enum class DD_MODEL
	{
		IF_CB_WN = 0, // 通过MW+NL固定无电离层
		RAW_CB_WN= 1, // 通过MW+NL固定非组合
		RAW_CB   = 2, // 通过WL+NL固定非组合
		NONE     = 3
	};
	/**
	*@brief	   Class for set ambiguity fixed xml
	*/
	class AMB_LIBRARY_EXPORT t_gsetamb : public virtual t_gsetbase
	{
	public:
		/** @brief enum of FIX mode. */
		enum class FIX_MODE {
			NO,
			ROUND,
			SEARCH,
			HOLD
		};

		/** @brief enum of UPD mode. */
		enum class UPD_MODE {
			UPD,  // wl upd + nl upd
			IRC,  // wl upd
			OSB   // no upd
		};

		/** @brief default constructor. */
		t_gsetamb();

		/** @brief default destructor. */
		virtual ~t_gsetamb();

		/**
		* @brief settings check.
		* @return void
		*/
		void check();

		/**
		* @brief settings help.
		* @return void
		*/
		void help();

		/**
		* @brief  get upd type.
		* @return upd type
		*/
		UPDTYPE get_updtype();

		/**
		* @brief  get ambiguity fixing mode .
		* @return ambiguity fixing mode
		*/
		FIX_MODE fix_mode();

		/**
		* @brief  get upd mode .
		* @return upd mode
		*/
		UPD_MODE upd_mode();

		DD_MODEL dd_mode();

		/**
		* @brief  get lambda ratio .
		* @return lambda ratio
		*/
		double lambda_ratio();

		/**
		* @brief  get bootstrapping rate.
		* @return bootstrapping rate
		*/
		double bootstrapping();

		/**
		* @brief  check whether add leo.
		* @return true add leo false not add leo
		*/
		bool addleo();

		/**
		* @brief  check whether the baseline is independ.
		* @return true independ baseline
		*/
		bool independ_baseline();

		/**
		* @brief  get minimum common time of two observation arc.
		* @return minimum common time of two observation arc
		*/
		double min_common_time();

		/**
		* @brief  get the limit of baseline length.
		* @return the limit of baseline length.
		*/
		double max_baseline_length();

		/**
		* @brief   widelane ambiguity interval.
		* @return  widelane ambiguity interval.
		*/
		double wl_interval();

		/**
		* @brief   get ambiguity decision.
		* @return  ambiguity decision.
		*/
		map<string, double> get_amb_decision(string str);

		/**
		* @brief   whether take partial ambiguity fixed mode.
		* @return  true or false.
		*/
		bool part_ambfix();

		/**
		* @brief   value's size which take partial ambiguity fixed mode.
		* @return  true or false.
		*/
		int part_ambfix_num();

		/**
        * @brief   whether apply IRC products
        * @return  true or false.
        */
		// bool apply_irc();

		/**
		* @brief   whether write carrier range
		* @return  true or false.
		*/
		bool carrier_range_out();

		/**
		* @brief   whether use carrier range
		* @return  true or false.
		*/
		bool apply_carrier_range();

		/**
		* @brief change string to FIX_MODE.
		* @param[in]  str   fix mode string
		* @return	  FIX_MODE
        */
		FIX_MODE str2fixmode(string str);

		/**
		* @brief change string to UPD_MODE.
		* @param[in]  str   upd mode string
		* @return	  UPD_MODE
		*/
		UPD_MODE str2upd_mode(string str);

		/**
		* @brief change string to UPDTYPE.
		* @param[in]  str   upd mode string
		* @return	  UPDTYPE
		*/
		UPDTYPE str2updmode(string str);

		bool isSetRefSat();

	protected:

	private:

	};

	/** @brief change string to FIXMODE. */
	AMB_LIBRARY_EXPORT t_gsetamb::FIX_MODE str2fixmode(string str);
}
#endif
