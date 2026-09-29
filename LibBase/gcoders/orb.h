/**
 * The orb file is the result of the orbit intergration, include position 
 *		and state trans matrix for satellites in every epoch, And orb file 
 *		will be used in pod process.
 * Following is an example for the header of orb file.
 *
 * @verbatim
		Sat Group: GPS											! Sat group (single sat type: sat type,two or more GNS sat:GNS,LEO:LEO)
		Sat Type: G												! Abbreviation for sat(GREAC)
		Sat List: 30											! Satellites' number
		  G01  G02  G05  G06  G07  G08  G09  G10  G11  G12		! Satellites' prn list
		  G13  G14  G15  G16  G17  G18  G19  G20  G21  G22
		  G23  G24  G25  G26  G27  G28  G29  G30  G31  G32
		Intege Control: true  11  60  300						! Integration settings (is extend time,adams order, adams input step, output step)
		Nequ: 36												! Number of equation
		Nepo: 297												! Number of epoch
		Beg Time: 58217 84900									! Begin time of oi
		End Time: 58219 1200									! End time of oi
		Ref Time: 58218 0										! reference time of oi
		END OF HEAD												! Mark for end of orb header
   @endverbatim
 *
 * An example for the record of the orb file
 * @verbatim
		*-5																					! Epoch (time = epoch * output_step + ref time)
		G01																					! Satellite prn
		       -1.48836816083279e+04     7.72072829515611e+02     2.18403055717414e+04		! Position of satellite
		        9.98811174091843e-01     4.75800336716186e-03    -3.34129392547388e-02		! Partial: dx/dx0,dy/dx0,dz/dx0
		        4.80130123042495e-03     9.77104183927932e-01    -6.99959626999930e-03		! Partial: dx/dy0,dy/dy0,dz/dy0
		       -3.34114933407117e-02    -6.93611992711632e-03     1.02466235202539e+00		! Partial: dx/dz0,dy/dz0,dz/dz0
		       -1.49943105660455e+03    -1.64525198153851e+00     1.68097076215478e+01		! Partial: dx/dvx0,dy/dvx0,dz/dvx0
		       -1.65610581393415e+00    -1.48825261229045e+03     2.41725371335385e+00		! Partial: dx/dvy0,dy/dvy0,dz/dvy0
		        1.68093452160285e+01     2.40134158453282e+00    -1.51249044165628e+03		! Partial: dx/dvz0,dy/dvz0,dz/dvz0
		        1.05364977493701e-03     3.47682589316881e-04     1.45501226316955e-04		! Partial: dx/dp0,dy/dp0,dz/dp0
		        3.39403603194023e-04    -1.05894211566652e-03     7.89522775210376e-05
		       -1.74349475482400e-04     2.88174375682843e-05     1.11167888291895e-03
		        2.19060025836463e-05    -7.03852982513234e-07    -1.45062680905524e-04
		       -1.72720583109935e-04     2.88888868111710e-05     1.10066974636375e-03
   @endverbatim
 *
 * @verbatim
	 History
	  -1.0	Jie Li	     2019-01-06 creat the file.
	  -1.1  Jie Li	     2019-03-30 adding doxygen style code remarks
	  -1.2  Huang Jiande 2020-09-12 Remove the sat without data in ics head
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		orb.h
 * @brief		The base class used to decode orb file information.
 *
 * @author      Jie Li, Wuhan University
 * @version		1.0.0
 * @date		2019-03-30
 */

#ifndef ORB_H
#define ORB_H

#include "gutils/ExportBase.h"

#include "gcoders/gcoder.h"

using namespace gnut;

namespace great
{
	/**
	*@brief	   Class for storaging head data of orb file
	*
	* The orb header contains the satellites list that attended OI 
	*	and the strategy taken in OI process.
	*/
	class BASE_LIBRARY_EXPORT t_orb_hd
	{
	public:
		/** @brief default constructor. */
		t_orb_hd();

		/** @brief default destructor. */
		virtual ~t_orb_hd();
		
		string			 sat_group;		///< satellite group single sys(GPS...) or Mix(GNS)
		set<char>		 sat_type;		///< satellite systems G R E C ...
		int				 sat_num;	    ///< number of satellites
		set<string>		 sat_prn;		///< satellites prn

		bool	         partial;		///< integrate variation eqution
		int		         admord;		///< order of adams numerical method
		int			     adminp;		///< integrate footstep (seconds)
		int			     admout;		///< output footstep (seconds)
		int	 	         nequ;			///< number of equation
		int		         nepo;			///< number of epoch

		t_gtime          beg;			///< begin time of data
		t_gtime          end;			///< end time of data
		t_gtime          t_ref;			///< reference time of ICS ( begin time of session)
	};

	/**
	*@brief	   Class for decoding and encoding the orb file
	*
	* Orb file is the result of orbit integration, mainly include the sat position 
	*		and state trans matrix.
	* The gcoder t_orb corresponding to the gdata gorb and gallorb
	*/
	class  BASE_LIBRARY_EXPORT t_orb : public t_gcoder
	{
	public:
		/**
		 * @brief default constructor.
		 *
		 * @param[in]  s        setbase control
		 * @param[in]  version  version of the gcoder
		 * @param[in]  sz       size of the buffer
		 */
		t_orb(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_orb();

		typedef vector< t_gtriple >  t_sat_funct;

		/**
		* @brief decode header of orb file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return 
			@retval >=0 consume size of header decoding
			@retval <0  finish reading
		*/
		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg) override;

		/**
		* @brief decode data body of orb file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return 
			@retval >=0 consume size of body decoding
			@retval <0  finish reading
		*/
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg) override;

		/**
		* @brief encode header of orb file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return  size of header encoding
		*/
		virtual  int encode_head(char* buff, int sz, vector<string>& errmsg);

		/**
		* @brief encode data body of orb file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return
			@retval >=0 consume size of body encoding
			@retval <0  finish encodeing
		*/
		virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

	protected:
	private:
		int       _epo;		 ///< current epoch
		t_orb_hd  _orb_hd;	 ///< orb head info
		string	  _sat;		 ///< satellite prn
		map<int, t_gtriple>	_funct_part;		///< part of funct (pos + dpos/dpos0 + dpos/dv0)
		vector<t_gtriple>	_partial_p;				///< partiial dpos/dp0
	};
}
#endif // !EGM2008_H
