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

#ifndef GLSQPROC_H
#define GLSQPROC_H

#include "gdata/gdata.h"
#include "gall/gallproc.h"
#include "gall/gallprod.h"
#include "gproc/glsq.h"
#include "gmodels/gmodel.h"
#include "gmodels/gprecisemodel.h"
#include "gmodels/glsqproc.h"
#include "gproc/glsqmatrix.h"
#include "gdata/gsatdata.h"
#include "gall/gallobs.h"

#include "gutils/ExportLSQ.h"

namespace great
{

	/**
	* @brief base class for lsq process
	*/
	class LSQ_LIBRARY_EXPORT t_glsqprocIF : public t_glsqproc
	{
	public:
		t_glsqprocIF();
		t_glsqprocIF(t_gsetbase* set,t_gallproc* data, t_glog* log);
		virtual ~t_glsqprocIF();
		
	protected:

	};
}


#endif