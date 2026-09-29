/**
 *
 * @verbatim
     History
     -1.0 Yujie Qin     created 2019-01-06 creat the file.
     -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
	 -1.2 Yujie Qin     2019-07-03 Adding Doxygen Style Code Remarks
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		    gallattitude.h
 * @brief		    Storage the LEO attitude data(more than one satellite)
*
 * @author          Yujie Qin created, Wuhan University
 * @version		    1.0.0
 * @date		    2019-03-27
 *
 */


#ifndef GALLATTITUDE_H
#define GALLATTITUDE_H

#include "gutils/ExportBase.h"
#include "gdata/gattitude.h"

namespace great
{
	/**
	*@brief	   Class for all LEO attitude data storaging
	*/
	class BASE_LIBRARY_EXPORT t_gallatt : public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gallatt();

		/** @brief default destructor. */
		virtual ~t_gallatt();

		/**
	    * @brief push attitude data for one LEO into the allattitude data map.
	    * @param [in]  sat         LEO satellite name.
	    * @param [in]  attitude    attitude data for one LEO.
	    */
		void add_attitude(string sat, t_gattitude& attitude);         

		/** 
		* @brief whether the attitude data is empty.
	    * @return  bool
	    *	@retval   true    attitude data is empty
	    *   @retval   false   attitude data is existent
	    */
		bool isEmpty();

		/**
	    * @brief return attitude data map for all LEO satellites.
	    * @return    attitude data map for all LEO satellites.
	    */
		map<string, t_gattitude>*	attitudes();	
		
		//todo 返回值指针引用统一
		/**
	    * @brief return attitude data for one LEO satellites.
	    * @return    attitude data for one LEO satellite.
	    */
		t_gattitude*				attitude(string sat);			

	protected:
		map<string, t_gattitude>    _attitudes;  ///< map of satellite name and attitude data for all LEO satellites

	};
} // namespace

#endif