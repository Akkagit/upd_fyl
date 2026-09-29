#ifndef LIDARFILE_H
#define LIDARFILE_H


/**
* @file lidarfile.h
* @brief read and decode lidar file
* @verbatim
	lidarfile is the lidar observation file,including TimeStamp and corresponding LIDARPATH
	The example of th lidar file( No File Head ):
	===================================================================================
		TimeStamp		Name
	    15.15           000001.txt/000001
	    15.25           000002.txt/000002

	    ...
	    120.55          000100.txt/000100
	===================================================================================
  @endverbatim
* @detials
* @verbatim
	History
	-1.0	WangHuidan	      2020-9-29 create the file.
  @endverbatim

* @author WangHuidan
* @version 1.0
* @date 2020-9-29
* @license
*/


#include "gins/gutility.h"
#include "gcoders/gcoder.h"
#include "gutils/gstring.h"
#include <Eigen/Eigen>
#include "gutils/ExportLIDAR.h"
using namespace std;
using namespace gnut;


namespace glidar {

	/**
	* @class t_lidarfile
	* @brief t_lidarfile Class for decoding and encoding the LIDAR file
	*
	* lidarfile is used for Front end tracking.
	* The document contains information about timestamp and lidarpath.
	* The gcoder t_lidarfile corresponding to the gdata t_glidardata.
	*/
	class LIDAR_LIBRARY_EXPORT t_lidarfile : public t_gcoder
	{
	public:

		/**
		* @brief constructor.
		*
		* @param[in]  s        setbase control
		* @param[in]  version  version of the gcoder
		* @param[in]  sz       size of the buffer
		*/
		t_lidarfile(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief destructor. */
		~t_lidarfile() {}


		/**
		* @brief decode the header of the LIDAR data file.
		* lidar file doesn't include head block
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
		* @brief decode the data body of the LIDAR data file.
		*
		* decode data body of LIDAR file, all the data read will store in the t_glidardata
		*
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return
			@retval >=0 consume size of header decoding
			@retval <0  finish reading
		*/
		virtual int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

	private:
		double _ts;			 // data interval
	};
}






#endif
