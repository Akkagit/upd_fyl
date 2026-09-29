/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2019-04-02  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    glsqproc.h
* @brief		base class for lsq process
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-18
*
*/

#ifndef GLSQPROCALL_H
#define GLSQPROCALL_H

#include "gmodels/glsqproc.h"
#include "gdata/gsatdata.h"
#include "gutils/ExportLSQ.h"

namespace great
{

	/**
	* @brief base class for lsq process
	*/
	class LSQ_LIBRARY_EXPORT t_glsqprocALL : public t_glsqproc
	{
	public:
		t_glsqprocALL();
		t_glsqprocALL(t_gsetbase* set, t_gallproc* data, t_glog* log);
		virtual ~t_glsqprocALL();
	
	protected:


	protected:
		IONMODEL _ion_model;
	};
}


#endif