/**
*
* @verbatim
	History
	 -1.0 Jiande Huang  2019-07-08  created
	 -2.0 Jiande huang  2020-04-12  add wn wl sg fix mode for pod all
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gpodlsqALL.h
* @brief		main control class for pod estimate
*
* @author       Jiande Huang, Wuhan University
* @version		1.0.0
* @date		    2019-04-18
*
*/


#ifndef GPODLSQALL_H
#define GPODLSQALL_H

#include "gmodels/glsqprocALL.h"
#include "gutils/ExportLSQ.h"
#include "gdata/gambcon.h"
#include "gset/gsetamb.h"

namespace great
{

	/**
	* @brief main control class for PPP
	*/
	class LSQ_LIBRARY_EXPORT t_gpodlsqALL :public t_glsqprocALL
	{
	public:
		/**
		* @brief constructor .
		*/
		t_gpodlsqALL();
		
		/**
		* @brief constructor.
		* @param[in] set   cfg settings.
		* @param[in] data  process data used in process.
		* @param[in] log   log file pointer.
		*/
		t_gpodlsqALL(t_gsetbase* set, t_gallproc* data, t_glog* log);

		/**
		* @brief destructor .
		*/
		virtual ~t_gpodlsqALL();

		/**
		* @brief batch processing .
		* @param[in] data   data used in proc.
		* @param[in] beg    beg  time.
		* @param[in] end    end  time.
		* @return	  true - successfully;
		*            false - failure
		*/
		bool ProcessBatch(t_gallproc* data, const t_gtime&beg, const t_gtime& end) override;



		/**
		* @brief  output the product file.
		* @return	  true - successfully;
		*            false - failure
		*/
		bool GenerateProduct() override;
	protected:
		/**
		* @brief init the product data class.
		* @param[in] data   the container for product data class.
		*/
		bool _initLsqProdData(t_gallprod* data) override;

		/**
		* @brief init the pars for processing.
		* @param[in] lsq   the class for lsq.
		*/
		bool _initLsqProcPars(t_glsq* lsq)      override;

		/**
		* @brief init the product data class.
		* @param[in] data   the container for process data class.
		*/
		bool _initLsqProcData(t_gallproc* data) override;

		/**
		* @brief update the ics value with lsq result.
		*/
		void _update_icsfile();
		void _readOfiles(int step, t_gallobs*& obs);
	protected:
		/**
		* @brief check the data or settings, and return checking result.
		* @return	  true - successfully;
		*            false - failure
		*/
		bool _check();
	};




}

#endif /*  GPPPLSQ_H */