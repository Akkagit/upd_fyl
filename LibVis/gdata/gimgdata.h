
#ifndef GIMGDATA_H
#define GIMGDATA_H

/**
* @file gimgdata.h
* @brief img data structure for storing img data
* @details Store data by vector

* @verbatim
	History
	-1.0	LiShengyu	      2020-05-20 create the file.
  @endverbatim

* @author LiShengyu
* @version 1.0
* @date 2020-05-20
* @license
*/

#include "gdata/gdata.h"
#include "gvis/gutility.h"
#include "gutils/ExportVIS.h"


using namespace gnut;

namespace gvis
{
	
	/**
	*@brief t_gimgdata Class for storing the IMG data,just the timestamp and imgpath
	*
	* imgfile is used for Front end tracking.
	* The data contains information about  timestamp and imgpath.
	* The gcoder t_gimgdata corresponding to the gcoder t_imgfile.
	*/
	class VIS_LIBRARY_EXPORT t_gimgdata : public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gimgdata();

		/** @brief default destructor. */
		virtual ~t_gimgdata() {}

		/**
		* @brief add one data into vector
		* @note data includes time and img0_path and img1_path
		* @param[in] 't' is time
		* @param[in] 'img0_path' is left image path
		* @param[in] 'img1_path' is right image path
		* @return value reprents success or failure of add
		*/
		int add_IMG(const double& t, const string& img0_path, const string& img1_path);
		
		/**
		* @brief load data into img_path according to imu_t
		* @param[in] imu_t
		* @param[out] img_t
		* @param[out] img_path
		* @return void
		*/
		void load(double& imu_t, double& img_t, imgPath& img_path);

		/** 
		* @brief return imgpath data vector size
		* @note not use
		* @return int
		*/
		int size();

	private:
		vector<imgPath> _vecimg;      //imgpath data
		vector<dataIMG> _rtvectimg;   //imgdata for real time
	
	};

}



#endif