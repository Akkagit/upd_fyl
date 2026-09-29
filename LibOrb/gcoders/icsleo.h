/**
 * The ics file of LEO is used for the orbit intergration.\n
 * Following is an example for the header of ics file.
 *
  * @verbatim
	 ## Header Of Initial Control State Of Satellite		      ! Start of header
	 %Satellite       = LEO  3                                ! GNS : type of satellite; G : means GPS; 18 : the number of satellite.
	 %PRN			  = SWARM-A  SWARM-B  SWARM-C  		      ! Satellite list
	 %Time            = GPS  58178  300.000  58178  21900.000 ! The period of ics file
	 ## End Of Header										  ! End of header
    @endverbatim
 *
 * An example for the record of the ics file
 * @verbatim
   =============================================================================================
   ! leo name          system type   satellite type
   =============================================================================================
	 SWARM-A            SWARM-A      LEO	
   =============================================================================================
   ! Par type     value                  sigma           MJD beg   SOD beg    MJD end   SOD end
   =============================================================================================
	 PXSAT        1.56392967436266e+03   9.035848e-04    58178     300        58178    21900  ! X in CRS
	 PYSAT        6.63362110499519e+03   2.867381e-04    58178     300        58178    21900  ! Y in CRS
	 PZSAT        2.02912686497059e+02   6.402869e-04    58178     300        58178    21900  ! Z in CRS
	 VXSAT       -3.87483145866032e-01   1.731613e-07    58178     300        58178    21900  ! Vx in CRS
	 VYSAT       -1.53905295703186e-01   6.767191e-07    58178     300        58178    21900  ! VY in CRS
	 VZSAT        7.63650164888467e+00   2.437530e-07    58178     300        58178    21900  ! Vz in CRS
	 DRAG_c:360   8.88250378072887e-01   8.234903e-01    58178     300        58178    21900  ! parameter about atmosphere drag
	 EMP_Sc1:90  -1.59978351605521e-01   1.133621e-01    58178     300        58178    5700   ! parameter about customer
     EMP_Sc1:90   2.13644500672205e-01   2.842921e-01    58178     5700       58178    11100
	 EMP_Sc1:90   8.17284056593180e-01   0.000000e+00    58178     11100      58178    16500
     EMP_Sc1:90   5.55791640406790e-01   0.000000e+00    58178     16500      58178    21900
     EMP_Cc1:90   8.28868414978178e-01   6.377596e-03    58178     300        58178    5700   ! parameter about customer 
	 EMP_Cc1:90   6.05235030222337e-01   3.627868e-02    58178     5700       58178    11100
	 EMP_Cc1:90   6.25753159023848e-01   9.856595e-01    58178     11100      58178    16500
	 EMP_Cc1:90   9.10137748903950e-02   0.000000e+00    58178     16500      58178    21900
	 EMP_Sa1:90  -9.42208940506553e-02   2.228026e-03    58178     300        58178    5700   ! parameter about customer 
     EMP_Sa1:90  -3.47781468391525e-01   8.668912e-01    58178     5700       58178    11100
     EMP_Sa1:90   1.26707088967123e-02   7.050789e-02    58178     11100      58178    16500
     EMP_Sa1:90  -1.90549644713089e-02   1.124509e+00    58178     16500      58178    21900
	 EMP_Ca1:90   1.28244206784028e-01   1.858908e-02    58178     300        58178    5700   ! parameter about customer 
	 EMP_Ca1:90   1.97621829643196e-01   1.525137e+00    58178     5700       58178    11100
	 EMP_Ca1:90  -4.99437096850430e-01   1.042403e+00    58178     11100      58178    16500
	 EMP_Ca1:90   9.65433268035083e-02   4.037868e-02    58178     16500      58178    21900
	 EMP_Sr1:90  -1.61738489058531e+00   8.569292e-03    58178     300        58178    5700   ! parameter about customer 
	 EMP_Sr1:90  -1.39881711866727e+00   1.401725e+00    58178     5700       58178    11100
	 EMP_Sr1:90  -8.03502773559564e-01   1.158722e+00    58178     11100      58178    16500
	 EMP_Sr1:90   2.36758417102378e-01   1.239587e+00    58178     16500      58178    21900
	 EMP_Cr1:90   6.19517924135953e-02   1.128111e-02    58178     300        58178    5700  ! parameter about customer
	 EMP_Cr1:90  -1.28249905443894e-01   1.466815e+00    58178     5700       58178    11100
	 EMP_Cr1:90   1.57972731162602e+00   1.384547e+00    58178     11100      58178    16500
	 EMP_Cr1:90   6.56655599274101e-01   1.612287e+00    58178     16500      58178    21900
	 END of SAT                                                                                     ! Mark for end of satellite
   @endverbatim
 *
 * @verbatim
    History
     -1.0 Hongmin Zhang 2019-01-06 creat the file.
     -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
	 -1.2 Yujie Qin     2019-04-10 Change the format of file,add encode head and data
	 -1.3 Yujie Qin     2019-07-01 Adding Doxygen Style Code Remarks
	 -1.4 Jiande Huang  2020-09-12 Standard code style and Optimize code logic
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file			icsleo.h
 * @brief		Decode and encode the LEO ICS data
 *
 * @author		Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date			2019-03-27
 *
 */


#ifndef ICSLEO_H
#define ICSLEO_H

#include "gcoders/gcoder.h"
#include "gdata/gics_leo.h"
#include "gall/gallics_leo.h"

using namespace std;

namespace great
{
	/**
	*@brief	   Class for decoding and encoding the LEO ICS data
	*
	* ics for LEO is used for the oi of navgation satellite.\n
	* The document contains information about the initial values of satellite position,
	*    velocity and dynamic parameters.\n
	* The gcoder t_icsnav corresponding to the gdata icsleo and gallics_leo
	*/
	class OI_LIBRARY_EXPORT t_ics_leo : public t_gcoder
	{
	public:
		/** @brief default constructor. */
		t_ics_leo();

		/**
		 * @brief default constructor.
		 *
		 * @param[in]  s        setbase control
		 * @param[in]  version  version of the gcoder
		 * @param[in]  sz       size of the buffer
		 */
		t_ics_leo(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_ics_leo();

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
			@retval >=0 consume size of header encoding
			@retval <0  finish encodeing
		 */
		virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

	protected:
		double          _period;		    ///< last time(end time minus start time, unit: s)

		t_gtime         _beg_epo;		///< reference time(epoch)
		t_gtime         _end_epo;		///< end time(epoch)
		string          _sys_type;		///< satellite system type(LEO)
		string          _leo_type;		///< satellite type(SWARM-A)
		set<string>		_leo_list;		///< list of satellite name
		int             _leo_num;		///< number of LEO satellites
		string          _time_type;		///< type of time
		int             _tab;			///< mark for lines in head
		int             _tab_data;		///< mark for lines in data

		//By zhangwei
		t_gics_leo*      _gicsleo;      ///< storage ICS data for one LEO satellite
		t_icsleo_record  _record;       ///< one record of ICS data

	};
} //namespace
#endif