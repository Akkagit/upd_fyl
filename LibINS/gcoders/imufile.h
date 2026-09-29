#ifndef IMUFILE_H
#define IMUFILE_H


/**
* @file imufile.h
* @brief read and decode imu file 
* @verbatim
	imufile is the INS observation file,including Gyro and Acce observation.
	The example of th imu file( No File Head ):
	===================================================================================
		Sow		Gyro_X		Gyro_Y		Gyro_Z		Acce_X		Acce_Y		Acce_Z
	 11111.01   0.2500		0.2500		0.2500		0.2500		0.2500		10.000
	 11111.02   0.2500		0.2500		0.2500		0.2500		0.2500		10.000
	 ...
	 11200.90   0.2500		0.2500		0.2500		0.2500		0.2500		10.000
	===================================================================================
  @endverbatim
* @detials
	The order of data record is different in different file, so provide some orders for different file,
	including a->g,g->a; rfu,flu,frd and so on.You can add extra order if need.

* @verbatim
	History
	-1.0	ShenZhiheng	      2018-12-15 creat the file.
  @endverbatim

* @author ShenZhiheng
* @version 1.0
* @date 2018-12-15
* @license 
*/


#include "gins/gutility.h"
#include "gcoders/gcoder.h"
#include "gutils/gstring.h"
#include <Eigen/Eigen>
#include "gutils/ExportINS.h"
using namespace std;
using namespace gnut;
using namespace great;

namespace gsins {

	/**
	* @class t_imufile
	* @brief t_imufile Class for decoding and encoding the IMU file
	*
	* imufile is used for sins mechanical arrangement.
	* The document contains information about gyro and acce observation.
	* The gcoder t_imufile corresponding to the gdata t_gimudata.
	*/
	class INS_LIBRARY_EXPORT t_imufile : public t_gcoder
	{
	public:

		/**
		* @brief constructor.
		*
		* @param[in]  s        setbase control
		* @param[in]  version  version of the gcoder
		* @param[in]  sz       size of the buffer
		*/
		t_imufile(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief destructor. */
		~t_imufile() {}


		/**
		* @brief decode the header of the IMU data file.
		* imu file doesn't include head block
		*
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return
			@retval >=0 consume size of header decoding
			@retval <0  finish reading
		*/
		virtual int decode_head(char* buff, int sz, vector<string>& errmsg);

		/**
		* @brief decode the data body of the IMU data file.
		*
		* decode data body of IMU file, all the data read will store in the t_gimudata
		*
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return
			@retval >=0 consume size of header decoding
			@retval <0  finish reading
		*/
		virtual int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

		virtual int decode_unicore(const string& line, double& t, vector<double>& v);

		virtual int decode_imr_fsas(const string& line, double& t, vector<double>& v);
		
		virtual int decode_starneto(const char* block, int sz, vector<pair<double, vector<double>>>& v);

		virtual bool available(const t_gtime& now);

		virtual  int encode_head(char* buff, int sz, vector<string>& errmsg);

		virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);


	private:
		double _tt;
		bool _complete;		 /// reading is complete 
		string _order;		 /// the imu record order
		double _ts;			 /// data interval
		double _freq;		 /// data interval
		UNIT_TYPE _GyroUnit; /// Gyro data Unit
		UNIT_TYPE _AcceUnit; /// Acce data Unit
		UNIT_TYPE _MagUnit;  /// Mag data Unit
	};
}






#endif