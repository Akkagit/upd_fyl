#ifndef POSFILE_H
#define POSFILE_H


/**
* @file posfile.h
* @brief read and decode lci_pos file
* @verbatim
	pos is the GNSS output file,including Position and Velocity.
	The example of th pos file( No File Head ):
	==========================================================================================================================================================
	  Sow         X       Y       Z       VX       VY       VZ      X_RMS      Y_RMS      Z_RMS      VX_RMS      VY_RMS      VZ_RMS    nSat	     PDOP   
	11111.00   1234567 1234567 1234567	  5.0	   5.0	    5.0      0.1        0.1        0.1         0.1         0.1         0.1     	15	      1.5	
	11112.00   1234567 1234567 1234567	  5.0	   5.0	    5.0      0.1        0.1        0.1         0.1         0.1         0.1     	15	      1.5	
	11113.00   1234567 1234567 1234567	  5.0	   5.0	    5.0      0.1        0.1        0.1         0.1         0.1         0.1     	15	      1.5	
	...
	11120.00   1234567 1234567 1234567	  5.0	   5.0	    5.0      0.1        0.1        0.1         0.1         0.1         0.1     	15	      1.5	
	==========================================================================================================================================================
  @endverbatim

* @detials
	The order of data record is different in different file, so provide some orders for different file,
	including XYZXYZ,XYZ,BLHENU,LBHENU,LBHNEU,LBHRPYENU,XYZ_VXYZ_RXYZ_RVXYZ_PDOP_NSAT and so on.
	You can add extra order if need.

* @verbatim
	History
	-1.0	ShenZhiheng	      2018-12-20 creat the file.
  @endverbatim

* @author ShenZhiheng
* @version 1.0
* @date 2018-12-20
* @license
*/

#include "gins/gutility.h"
#include "gcoders/gcoder.h"
#include "gdata/gposdata.h"
#include <Eigen/Eigen>
#include "gutils/ExportINS.h"
using namespace std;
using namespace gnut;

namespace gsins {

	/**
	*@brief t_posfile Class for decoding and encoding the pos file
	*
	* imufile is used for loosely coupled integration.
	* The record contains information about position,velocity,their std and so on.
	* The gcoder t_posfile corresponding to the gdata t_gposdata.
	*/
	class INS_LIBRARY_EXPORT t_posfile : public t_gcoder
	{
	public:

		/**
		* @brief constructor.
		*
		* @param[in]  s        setbase control
		* @param[in]  version  version of the gcoder
		* @param[in]  sz       size of the buffer
		*/
		t_posfile(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);
		
		
		/** @brief destructor. */
		~t_posfile();

		/**
		* @brief decode the header of the pos data file.
		* pos file doesn't include head block
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
		* @brief decode the data body of the pos data file.
		*
		* decode data body of pos file, all the data read will store in the t_gposdata
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
		double TS;		/// pos data interval
		string _order;  /// the pos record order


	};

}






#endif