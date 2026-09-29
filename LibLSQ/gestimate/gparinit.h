///**
//*
//* @verbatim
//	History
//	 -1.0 ZhengHJ  2018-12-02  created
//  @endverbatim
//* Copyright (c) 2018, Wuhan University. All rights reserved.
//*
//* @file		    gallprocdata.h
//* @brief		header files of least squares estimator
//*
//* @author       ZhengHJ, Wuhan University
//* @version		1.0.0
//* @date		    2019-04-17
//*
//*/
//
//#ifndef GPARINIT_H
//#define GPARINIT_H
//
//
//#include "gset/gsetbase.h"
//#include "gestimate/glsq.h"
//#include "gall/gallics_nav.h"
//#include "gdata/gdata.h"
//#include "gall/gallobj.h"
//
//namespace gnut
//{
//	/**
//	* @brief init par for lsq estimator
//	*/
//	class t_glsqparInit
//	{
//	public:
//		t_glsqparInit(t_glsq* lsq,t_gsetbase* set);
//		~t_glsqparInit();
//
//		// SET
//
//		// for single
//
//		// for whole 
//		bool initRecCrdPars(set<string>& sites, t_gallobj* allobj, double std_threshold=-1);
//		bool initRecTroPars(set<string>& sites);
//		bool initRecClkPars(set<string>& sites);
//		bool initSysBiaPars(set<string>& sites, set<string>& sys);
//		bool initSatClkPars(set<string>& sats);
//		bool initGnsIcsPars(set<string>& sats, t_gdata* icsdata);
//		bool initCrsEopPars(t_gdata* poledata);
//
//	private:
//		t_glsq* _lsq;
//		t_gsetbase* _set;
//		t_gtime _beg_time;
//		t_gtime _end_time;
//	};
//}
//
//#endif // !GPARINIT_H
