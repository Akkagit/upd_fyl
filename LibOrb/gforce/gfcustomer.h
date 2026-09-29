/**
 *
 * @verbatim
     History
      -1.0 Yujie Qin     2019-01-06 Creat the file.
      -1.1 Yujie Qin     2019-03-05 edit for final version.
      -1.2 Yujie Qin     2019-03-15 Code comments are added in Doxygen format.
	  -1.3 Jiande Huang  2019-03-28 #ifndef name was wrong, from GCUSTOMER_H to GFCUSTOMER_H
	  -1.4 Yujie Qin     2019-07-03 Adding Doxygen Style Code Remarks
	  -1.5 Jiande Huang  2020-10-01 Standard code style and Optimize code logic
   @endverbatim
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		     gfcustomer.h
 * @brief		     The class for customer model, get the acceleration caused by customer.
 *
 * @author           Yujie Qin, Wuhan University
 * @version		     1.0.0
 * @date		     2019-01-06
 *
 */
#ifndef GFCUSTOMER_H
#define GFCUSTOMER_H

#include "gforce/gforce.h"
#include "gforce/gfarg_init.h"
#include "gio/glog.h"
#include "gutils/gconst_oi.h"
#include "newmat/newmat.h"
using namespace std;

namespace great
{
	/**
	 * @brief This is the base class for customer model.
	 */
	class OI_LIBRARY_EXPORT t_gfcustomer :public t_gforce
    {
    public:
	   //todo ×¢ÊÍÓÐÎó
	   /** @brief constructor. */
	   t_gfcustomer();

	   /** @brief constructor. 
	   * @param[in] arg   the init argument for init one force model.
	   * @param[in] log   the log file for write log.
	   */
	   t_gfcustomer(t_init_arg* arg,t_glog* log);

	   /** @brief destructor. */
	   virtual ~t_gfcustomer();

	   /** 
	   * @brief calculate the accleration and the corresponding partial matrix.
	   * @param[in] argument Parameters to be used in calculation
	   * @return bool
			@retval true     calculate normally
			@retval false    have some error
	   */
	   bool calculate(t_gforce_calarg* argument);
	  
	  
   protected:
	   //todo ÃüÃû
	  /** 
	  * @brief set mode for the customer model: constant, linear, one cpr and two cpr.
	  * @param[in] period the period of sat.
	  * @param[in] dt     the second of dynamic time.
	  */
       void set_mode(double period, double dt);

	   double               _mass = 0.0;          ///< the mass of sat.
	   double               _period_time = 0.0;
	   bool                 _same_sat  = false;     ///< the same sat or not
	   bool                 _isPartial = false;     ///< calculate the partial or not.
       double               _mode[6];       ///< the customer parameter mode(constant, linear, one cpr and two cpr)
	   map<ICS_PAR, double> _customer_pars; ///< the customer par from ics.
	   vector<ICS_PAR>      _cus_par;        ///< 18 cusstomer parameters
	   map<int, int>        _par_index;      ///< Modified by zhangwei
	 };
} // namespace

#endif