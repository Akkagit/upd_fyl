/**
*
* @verbatim
History
-1.0 GREAT   		 2019-11-09 creat the file.
-1.1 Jiande Huang    2020-09-12 Standard code style and Optimize code logic
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gdeadata.h
* @brief		Class for storaging ocean_pole_tide data
*
* @author       GREAT, Wuhan University
* @version		1.0.0
* @date		    2019-11-09
*
*/

#ifndef  GDEADATA_H
#define  GDEADATA_H

#include <vector>

#include "gdata/gdata.h"
#include "gutils/ExportOI.h"

using namespace std;
using namespace gnut;

namespace great
{


	/**
	*@brief	   Class for storaging ocean_pole_tide data
	*
	* The class contains Doodson index and coeffcients of
	*		oceantide in different order and provide the
	*		interface to get data.
	*/
	class OI_LIBRARY_EXPORT t_gdeadata : public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gdeadata();

		/** @brief default destructor. */
		virtual ~t_gdeadata() {};

		/**
		* @brief add one oceanide record to _oct.
		* @param[in]   m		index m
		* @param[in]   n		index n
		* @param[in]   A B AI BI	coefficient
		*/
		void set_rec(const int& n, const int& m, double A,double B,double AI,double BI);

		/**
		* @brief add a load coefficient value to vector _load.
		* @param[in]   load		load coefficient value
		*/
		void set_load(double load);


		/** @brief whether the ocean tide data is empty.
		* @return  bool
		*	@retval   true   ocean tide data is empty
		*   @retval   false  ocean tide data is existent
		*/
		bool is_empty();

		/**
		* @brief get all oceantide coefficient.
		* @return   all oceantide coefficient
		*/
		map <int, map< int, double >>& Anm();
		map <int, map< int, double >>& Bnm();
		map <int, map< int, double >>& AInm();
		map <int, map< int, double >>& BInm();

		/**
		* @brief get all load coefficient.
		* @return   all load coefficient
		*/
		vector<double>& load();


	protected:
		map <int, map< int, double >> _Anm;
		map <int, map< int, double >> _Bnm;
		map <int, map< int, double >> _AInm;
		map <int, map< int, double >> _BInm;
		vector<double> _load;
	};
} //namespace

#endif