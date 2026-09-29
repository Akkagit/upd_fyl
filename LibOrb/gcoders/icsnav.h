/**
 * The ics file of navigation satellite is used for the orbit intergration.\n
 * Following is an example for the header of ics file.
 *
 * @verbatim
	 ## Header Of Initial Control State Of Satellite		! Start of header
	 %Satellite  = GNS L G 18                               ! GNS : type of satellite; G : means GPS; 18 : the number of satellite.
	 %PRN   = G01 G02 G05 G06 G07 G08 G09 G10 G11 ++		! Satellite list, end with ++
	 %PRN   = G12 G13 G14 G15 G16 G17 G18 G19 G20 ++		! Satellite list, end with ++
	 %Time  = GPS  57770  0.000  57771  0.000				! The period of ics file
	 ## End Of Header										! End of header
   @endverbatim
 *
 * An example for the record of the ics file
 * @verbatim
   =============================================================================================
   ! PRN           system type            satellite type
   =============================================================================================
	 G01           GPS                    GNS														! PRN and type of satellite
   =============================================================================================
   ! Par type      value                  sigma           MJD beg   SOD beg    MJD end   SOD end
   =============================================================================================
	 PXSAT         1.52838439706610e+04   2.672656e+01    57770     0.00000    57771     0.00000	! X in CRS
	 PYSAT         1.24351939819255e+03   2.672951e+01    57770     0.00000    57771     0.00000    ! Y in CRS
	 PZSAT        -2.18097748020348e+04   2.672115e+01    57770     0.00000    57771     0.00000	! Z in CRS
	 VXSAT        -8.54515649432574e-01   2.160778e-01    57770     0.00000    57771     0.00000    ! Vx in CRS
	 VYSAT         3.74652943658759e+00   6.971100e-02    57770     0.00000    57771     0.00000	! VY in CRS
	 VZSAT        -3.63836156422601e-01   1.661940e-01    57770     0.00000    57771     0.00000	! Vz in CRS
	 Kd_BERN      -2.67490233483971e+00   2.645756e+02    57770     0.00000    57771     0.00000	! Kd of ECOM 5
	 Ky_BERN       1.59007172727675e-01   2.672909e+02    57770     0.00000    57771     0.00000	! Ky of ECOM 5
	 Kb_BERN       3.20937856176232e-01   2.672860e+02    57770     0.00000    57771     0.00000	! Kb of ECOM 5
	 Kbc1_BERN    -4.88735260638650e-01   2.671676e+02    57770     0.00000    57771     0.00000	! Kbc1 of ECOM 5
	 Kbs1_BERN     7.95779696295994e-02   2.664011e+02    57770     0.00000    57771     0.00000	! Kbs1 of ECOM 5
	 END of SAT                                                                                     ! Mark for end of satellite
   @endverbatim
 *
 * @verbatim
	 History
	  -1.0	Jie Li	       2019-01-06 creat the file.
	  -1.1  Jiande Huang   2019-02-27 change the way to read the file.
	  -1.2  Hongjie Zheng  2019-03-22 fix the problem of limited size of buff.
	  -1.3  Jie Li	       2019-03-30 adding Doxygen Style Code Remarks
	  -1.4  Jiande Huang   2020-09-12 Remove the sat without data in ics head
	  -1.5  Jiande Huang   2020-09-13 Standard code style and Optimize code logic
   @endverbatim
 *
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		icsnav.h
 * @brief		The base class used to decode navigation satellite ics file information.
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-01-06
 */

#ifndef ICSNAV_H
#define ICSNAV_H

#include "gcoders/gcoder.h"
#include "gmodels/gpar.h"
#include "gutils/gtime.h"
#include "gdata/gics_nav.h"

#include "gutils/ExportOI.h"

namespace great
{
	/**
	*@brief	   Class for decoding and encoding the GNS ICS data
	* 
	* ics for nav sat is used for the oi of navgation satellite.\n 
	* The document contains information about the initial values of satellite position, 
	*    velocity and dynamic parameters.\n
	* The gcoder t_icsnav corresponding to the gdata icsnav and gallics_nav
	*/
	class OI_LIBRARY_EXPORT t_ics_nav : public t_gcoder
	{

	public:
		/** @brief default constructor. */
		t_ics_nav();

		/** 
		 * @brief constructor. 
		 *
		 * @param[in]  s        setbase control
		 * @param[in]  version  version of the gcoder
		 * @param[in]  sz       size of the buffer
		 */
		t_ics_nav(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_ics_nav();

		/** 
		 * @brief decode the header of the ICS data file. 
		 * 
		 * The function is used for decoding the head of ics file.\n
		 * pay attention to the buff and the size of buff which may cause some trouble when 
			using the wrong value in decoding.
		 *
	     * @brief decode header of ics file
	     * @param[in]  buff        buffer of the data
	     * @param[in]  sz          buffer size of the data
	     * @param[in]  errmsg      error message of the data decoding
	     * @return 
			@retval >=0 consume size of header decoding
			@retval <0  finish reading
	     */
		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg);

		/** 
		 * @brief decode the data body of the ICS data file. 
		 *
		 * decode data body of ics file, all the data read will store in the gallics_nav
		 *
	     * @param[in]  buff        buffer of the data
	     * @param[in]  sz          buffer size of the data
	     * @param[in]  errmsg      error message of the data decoding
	     * @return 
		    @retval >=0 consume size of body decoding
			@retval <0  finish reading
		 */
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

		/** 
		 * @brief encode the header of the ICS data file. 
	     *
	     * @brief encode header of ics file
	     * @param[in]  buff        buffer of the data
	     * @param[in]  sz          buffer size of the data
	     * @param[in]  errmsg      error message of the data decoding
	     * @return  size of header encoding
		 */
		virtual  int encode_head(char* buff, int sz, vector<string>& errmsg);

		/** 
		 * @brief decode the data body of the ICS data file. 
		 *
	     * @param[in]  buff        buffer of the data
	     * @param[in]  sz          buffer size of the data
	     * @param[in]  errmsg      error message of the data decoding
		 * @return
			@retval >=0 consume size of data encoding
			@retval <0  finish encodeing
		 */
		virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

	protected:
		double          _period;		///< last time(end time minus start time, unit: s)
		t_gtime         _beg_epo;		///< reference time(epoch)
		t_gtime         _end_epo;		///< end time(epoch)

		int             _sat_num;		///< number of satellites
		string          _sys_type;		///< satellite system type(GNS/LEO)
		set<string>		_sat_list;		///< PRN of satellites
		set<string>     _sys_list;		///< abbreviation of systems
		
		vector<string>	_par_list;		///< parameters of satellite initial control state
		bool            _isMixed;		///< the satellite system is mixed or single(M: true, L: false)
		string          _time_type;		///< type of time
		
		int             _tab;			///< mark for lines in head
		int             _tab_data;		///< mark for lines in data

		t_gics_nav*     _gics_nav;			///< point to ics data
		t_ics_record    _ics_record;		///< one record data of a satellite
	};
} //namespace

#endif