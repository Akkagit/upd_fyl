#ifndef GPOSDATA_H
#define GPOSDATA_H

//todo  more description
/**
* @file gposdata.h
* @brief lci_pos data structure for storing pos data
* @details Store data by vector

* @verbatim
	History
	-1.0	ShenZhiheng	      2018-12-30 creat the file.
  @endverbatim

* @author ShenZhiheng
* @version 1.0
* @date 2018-12-30
* @license
*/


#include "gdata/gdata.h"
#include "gins/gutility.h"
#include "gutils/ExportINS.h"
using namespace gnut;

namespace gsins
{
	/**
	*@brief t_gposdata Class for storing the pos data
	*
	* t_gposdata is used for loosely coupled integration.
	* The data contains information about position,velocity,their std and so on.
	* The gcoder t_gposdata corresponding to the gcoder t_posfile.
	*/
	class INS_LIBRARY_EXPORT t_gposdata : public t_gdata
	{
	public:

		/** @brief default constructor. */
		t_gposdata();

		/** @brief default destructor. */
		~t_gposdata();


		/**
		* @struct data_pos
		* @brief describe pos data information
		*/
		struct data_pos
		{
			double t;
			Eigen::Vector3d pos, vn;
			Eigen::Vector3d Rpos, Rvn;
			double PDOP;
			int nSat;
			bool amb_state;
			bool operator<(const data_pos& dp)const
			{
				return (t < dp.t);
			}
			bool operator<(const double& d)const
			{
				return (t < d);
			}
			bool operator>(const double& d)const
			{
				return (t > d);
			}
			void operator=(const double& d)
			{
				t = d; pos = Eigen::Vector3d::Zero(); vn = Eigen::Vector3d::Zero();
			}
		};


		/**
		* @brief add one data into vector
		* @note data includes time and position and velocity
		* @param[in] t			time
		* @param[in] pos		position
		* @param[in] v			velocity
		* @return int			value reprents success or failure of add
		*/
		int add_pos(const double& t, const Eigen::Vector3d& pos, const Eigen::Vector3d& v);


		/**
		* @brief add one data into vector
		* @note data includes time and position and velocity
		* @param[in] t
		* @param[in] pos
		* @param[in] v
		* @return int
		*/
		int add_pos(const double& t, const Eigen::Vector3d& pos, const Eigen::Vector3d& v,
			const Eigen::Vector3d& R_pos, const Eigen::Vector3d& R_v, const double& pdop,const int& n);

		/**
		* @brief add one data into vector
		* @note data includes time and position and velocity
		* @param[in] t
		* @param[in] pos
		* @param[in] v
		* @return int
		*/
		int add_pos(const double& t, const Eigen::Vector3d& pos, const Eigen::Vector3d& v,
			const Eigen::Vector3d& R_pos, const Eigen::Vector3d& R_v, const double& pdop, const int& n,const bool& state);

		/**
		* @brief load data into pos and v and t accroding to direction
		* @note direction stands for loading direction
		* @param[in] direction    direction = forward(true)/backward(false)
		* @param[out] pos
		* @param[out] v
		* @param[out] t
		*/
		void load(double& t, Eigen::Vector3d& pos, Eigen::Vector3d& v, bool direction);  
		

		/**
		* @brief load data into pos and v and t accroding to direction
		* @param[in] t1    find the data nearest to t1 
		* @param[out] pos
		* @param[out] v
		* @param[out] t
		*/
		void load(double& t, Eigen::Vector3d& pos, Eigen::Vector3d& v, bool direction, double t1);
		
		/**
		* @brief load data into pos and v and t accroding to direction
		* @param[in] t1    find the data nearest to t1
		* @param[out] pos		position
		* @param[out] v			velocity
		* @param[out] t
		* @param[out] Qpos		position variance
		* @param[out] Qv		velocity variance
		* @param[out] PDOP		PDOP
		* @param[out] nSat		number of satellite
		*/
		void load(double& t, Eigen::Vector3d& pos, Eigen::Vector3d& v, Eigen::Vector3d& Qpos, Eigen::Vector3d& Qv, 
			double& PDOP, int &nSat, bool direction, double t1);
		
		/**
		* @brief load data into pos and v and t accroding to direction
		* @param[in] t1    find the data nearest to t1
		* @param[out] pos		position
		* @param[out] v			velocity
		* @param[out] t
		* @param[out] Qpos		position variance
		* @param[out] Qv		velocity variance
		* @param[out] PDOP		PDOP
		* @param[out] nSat		number of satellite
		*/
		void load(data_pos& data, double t1);


		/**
		* @brief judge the info valid
		* @note data includes time and position and velocity
		* @param[in] t
		* @param[in] pos
		* @param[in] tot
		* @param[in] threshold		check the data is valid or not by threshold
		* @return true or false is valid or invalid
		*/
		static bool IsValid(const double& t, const Eigen::Vector3d& pos, const double& tot, const double& threshold);
		
		virtual double beg_obs();  // get first epoch
		virtual double end_obs();  // get end epoch


	private:
		vector<data_pos> _vecpos;   /// latitude and longitude are saved by radian
		set<data_pos> _setpos;		/// latitude and longitude are saved by radian
		int _ptr;
	};


}



#endif