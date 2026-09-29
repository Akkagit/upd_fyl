#ifndef IMGFILE_H
#define IMGFILE_H


/**
* @file imgfile.h
* @brief read and decode img file
* @verbatim
	imgfile is the IMG observation file,including TimeStamp and corresponding IMGPATH
	The example of th img file( No File Head ):
	===================================================================================
		TimeStamp		Name
	    15.15           000001.jpg/000001
	    15.25           000002.jpg/000002

	    ...
	    120.55          000100.jpg/000100
	===================================================================================
  @endverbatim
* @detials
* @verbatim
	History
	-1.0	LiShengyu	      2020-5-20 create the file.
  @endverbatim

* @author LiShengyu
* @version 1.0
* @date 2020-5-20
* @license
*/


#include "gins/gutility.h"
#include "gcoders/gcoder.h"
#include "gutils/gstring.h"
#include <Eigen/Eigen>
#include "gutils/ExportVIS.h"
using namespace std;
using namespace gnut;


namespace gvis {

	/**
	* @class t_imgfile
	* @brief t_imgfile Class for decoding and encoding the IMG file
	*
	* imgfile is used for Front end tracking.
	* The document contains information about timestamp and imgpath.
	* The gcoder t_imgfile corresponding to the gdata t_gimgdata.
	*/
	class VIS_LIBRARY_EXPORT t_imgfile : public t_gcoder
	{
	public:

		/**
		* @brief constructor.
		*
		* @param[in]  s        setbase control
		* @param[in]  version  version of the gcoder
		* @param[in]  sz       size of the buffer
		*/
		t_imgfile(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief destructor. */
		~t_imgfile() {}


		/**
		* @brief decode the header of the IMG data file.
		* img file doesn't include head block
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
		* @brief decode the data body of the IMG data file.
		*
		* decode data body of IMG file, all the data read will store in the t_gimgdata
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
