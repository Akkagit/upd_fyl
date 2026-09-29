#ifndef GSETINS_H
#define GSETINS_H

//todo  more description
/**
* @file gsetins.h
* @brief implements settings for ins
* @mainpage
* @author ShenZhiheng
* @email
* @version 1.0
* @date 2018-12-15
* @license
*/

#define XMLKEY_INS "ins"
#include "gins/gutility.h"
#include "gset/gsetbase.h"
using namespace gsins;
using namespace gnut;

namespace gsins
{
	//todo ×¢ÊÍºÍÃüÃû
	/**
	* @class t_gsetins
	* @brief ins configure
	* @author ShenZhiheng
	*/
	class INS_LIBRARY_EXPORT t_gsetins : public virtual t_gsetbase
	{
	public:
		t_gsetins();
		~t_gsetins();

		void check();
		void help();

		double ts();
		int freq();
		gsins::UNIT_TYPE GyroUnit();
		gsins::UNIT_TYPE AcceUnit();
		gsins::UNIT_TYPE MagUnit();
		string order();
		ODR_TYPE out_order();
		Eigen::Vector3d pos();
		Eigen::Vector3d vel();
		Eigen::Vector3d att();
		bool align();

		bool out_intsec();
		int out_freq();

		int cps();
		int subsample();

		double start();
		double end();
		double align_time();
		ALIGN_TYPE align_type();

		vector<string> motion();
		Eigen::Vector3d acce_bias(); /**< constant accelerator bias. for simulation use. */
		Eigen::Vector3d gyro_bias(); /**< constant gyroscope bias. for simulation use. */
		

	protected:
		double _ts;						/**< sample interval. */
		int _freq ;						/**< sample freq. */
		string _cps;						/**< compensation mode. */
		int _subsample;					/**< subsample number. */
		string _order;					/**< imu data order. */
		string _out_order;				/**< output data format.*/
		int _out_freq;				/**< output data format.*/
		bool _int_sec;				/**< output data format.*/
		string _GyroUnit;     /**< Gyro data unit. */
		string _AcceUnit;		/**< Acce data unit. */
		Eigen::Vector3d _pos, _vel, _att;  /**< pos is XYZ,vel is Vxyz,att is P,R,Y(deg). */
		double _start, _end;				/**< start time and end time. */
		double _align_time;			  	/**< coarse align time. */
		string _align_type;			/**< coarse align type. */
		Eigen::Vector3d _acce_bias; /**< constant accelerator bias. for simulation use. */
		Eigen::Vector3d _gyro_bias; /**< constant gyroscope bias. for simulation use. */
	};
}
#endif