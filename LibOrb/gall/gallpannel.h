/**
*
* @verbatim
    History
     -1.0 Yujie Qin     2019-01-06 creat the file.
     -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
	 -1.2 Yujie Qin     2019-07-03 Adding Doxygen Style Code Remarks
	 -1.3 Huang Jiande  2020-09-12 Standard code style
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		      gallpannel.h
* @brief	      Storage the LEO pannel data(more than one satellite)
*
* @author         Yujie Qin created, Wuhan University
* @version	      1.0.0
* @date		      2019-03-27
*
*/


#ifndef GALLPANNEL_H
#define GALLPANNEL_H

#include "gdata/gpannel.h"
#include "gdata/gdata.h"

using namespace std;
using namespace gnut;

namespace great
{
	/**
	*@brief	   Class for all LEO pannel data storaging
	*/
	class OI_LIBRARY_EXPORT t_gallpannel : public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gallpannel();

		/** @brief default destructor. */
		virtual ~t_gallpannel();

		/**
	    * @brief push pannel data for one LEO into the data map.
	    * @param[in]   sat         LEO satellite name.
	    * @param[in]  pannels      pannel data for one LEO.
	    */
		void add_sat_pannel(string sat, t_gpannel pannel); 

		/** 
		* @brief whether the pannel data is empty.
	    * @return  bool
	    *	@retval   true    pannel data is empty
	    *   @retval   false   pannel data is existent
	    */
		bool isEmpty();

		/**
	    * @brief return all LEO sat's pannels.
	    * @return   all LEO sat's pannel data
	    */
		std::map<std::string, t_gpannel>& get_sat_pannels();

		/**
	    * @brief return one sat pannels.
	    * @param[in]   sat         LEO satellite name.
	    * @return      pannel data for one LEO.
	    */
		t_gpannel* get_sat_pannel(std::string sat);	

	protected:
		map<string, t_gpannel> _sat_pannels; ///< pannels for all satellties
	};
} //namespace

#endif