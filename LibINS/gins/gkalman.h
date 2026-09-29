#ifndef GKALMAN_H
#define GKALMAN_H

/**
* @file gkalman.h
* @details 
* @note t_gkalman Class for Kalman Filter Algorithm

* @verbatim
		History
		-0.1	ShenZhiheng	       2018-10-15 creat the file.
		-0.3	ShenZhiheng	       2019-03-20 add some base functions
		-0.5	ShenZhiheng	       2019-06-10 rename class name
		-0.7	ShenZhiheng	       2019-12-02 restruct the file
  @endverbatim

* @author ShenZhiheng
* @version 0.7
* @date 2018-10-15
* @license
*/

#include "gutility.h"
#include "gutils/ExportINS.h"

using namespace std;

namespace gsins
{

	/**
	* @class t_gkalman
	* @brief t_gkalman Class for Kalman Filter Algorithm
	*
	* Adaptive,Sequential,Variance constraint
	*/
	class INS_LIBRARY_EXPORT t_gkalman
	{

	public:

		/** @brief default constructor. */
		t_gkalman();

		/**
		* @brief Constructor
		* @note set parameter value by gset
		*/
		t_gkalman(gnut::t_gsetbase* gset);

		/**
		* @brief Constructor
		* @note set parameter value by dimension
		*/
		t_gkalman(int _nq, int _nr);


		/**
		* @brief Initial Function
		*
		* purely virual function,and should be inherited
		*/
		virtual void init(void) = 0;

		/**
		* @brief set process matrix Ft related to time,need to be Discreted
		* @note pure virtual function
		* purely virual function,and should be inherited
		*/
		virtual void set_Ft(void) = 0;
		virtual void set_Ft1(void) = 0;

		/**
		* @brief set measurement matrix Hk
		*
		* purely virual function,and should be inherited
		*
		*/
		virtual void set_Hk(void) = 0;

		/**
		* @brief set measurement Zk
		* @note pure virtual function
		* purely virual function,and should be inherited
		*/
		virtual void set_Meas(void) = 0;

		/**
		* @brief feed back
		*
		* @note parameter is feedbacking time
		* @param[in] fbts  interval of feedback
		*/
		virtual void feedback(double fbts);


		/**
		* @brief detect outlier by post-residual
		*
		* @param[in] threshold			post-residual
		* @return 
			outlier's index
		*/
		virtual int outlier_detect(double threshold = 0);

		virtual int _IGG_detect();

		/**
		* @brief Kalman TimeUpdate
		*
		* @note include:State transition,Variance update,(Feedback)
		* @param[in] kfts
		* @param[in] fback
		* @return 
		*/
		void time_update(double kfts, int fback = 1);

		/**
		* @brief Kalman measurement update
		*
		* @note include:calculate filter gain,update variance,(Feedback)
		* @param[in] fading   factor of fading
		* @return
		*/
		virtual void meas_update(double fading = 1.0);

		/**
		* @brief measurement noise adaptive
		*
		* @param[in] i		index of measurement
		* @param[in] r		residual
		* @param[in] Pr		variance of predicting
		*/
		void R_adaptive(int i, double r, double Pr);
		

		/** @brief state variance constraint */
		void Pk_constrain(void);

		/** @brief set measurement flag */
		void set_meas_flag(int flag);


	protected:
		Eigen::VectorXd _v_post;							/// post residual

	public:
		double kftk;										/// now
		int nq, nr;											/// dimension of state and meas
		int measflag;										/// flag of measurement index
		Eigen::MatrixXd Ft, Pk, Hk, Rk, Pk_Sav, Phik;				/// Kalman Matrix
		Eigen::VectorXd Xk, Zk, Qt, Rt, rts, Pmax, Pmin,
			Rmax, Rmin, Rb, Rbeta,
			FBTau, FBMax, FBXk, FBTotal;					/// Kalman Vector

	};

}


#endif