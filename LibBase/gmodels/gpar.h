// jdhuang 2019_07_09 remove  CLK_SAT

#ifndef GPAR_H
#define GPAR_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: parametres class
  Version: $ Rev: $

  2011-04-18 /PV: created

-*/

#include <string>

#include "gdata/gsatdata.h"
#include "gutils/gtriple.h"
#include "gset/gsetproc.h"



using namespace std;

namespace gnut
{


	class BASE_LIBRARY_EXPORT t_gparhead;

	class BASE_LIBRARY_EXPORT t_gtimearc
	{
	public:

		t_gtimearc(const t_gtime& beg,const t_gtime& end);
		~t_gtimearc();

		bool operator != (const t_gtimearc& Other) const;
		bool operator == (const t_gtimearc& Other) const;
		bool operator < (const t_gtimearc& Other) const ;
		bool operator <= (const t_gtimearc& Other) const;
		bool operator > (const t_gtimearc& Other) const;
		bool operator >= (const t_gtimearc& Other) const;

		bool inside(const t_gtimearc& Other) const;

		t_gtime begin;
		t_gtime end;
	};
	class BASE_LIBRARY_EXPORT t_gpar
	{
	public:

		const static string PAR_STR_SEP;
		const static vector<string> PAR_STR;

		enum t_type
		{
			CRD_X, CRD_Y, CRD_Z,                                        // coordinates
			TRP, GRD_N, GRD_E, SION, VION,                              // atmospheric parameters
			CLK, CLK_SAT, G15CLK, E15CLK, C15CLK, J15CLK,               // clocks
			G15CLK_SAT, E15CLK_SAT, C15CLK_SAT, J15CLK_SAT,             // satellite clocks for 13 combination (add by xiongyun)
			CLK_G, CLK_E, CLK_C, CLK_R, CLK_J,                          // receiver clocks for different systems (add by xiongyun)
			IFB_C3, IFB_C4, IFB_C5, IFCB_F3, IFCB_F4, IFCB_F5,          // inter-freq. code biases for FREQ_3, FREQ_4, FREQ_5, inter-freq. clock bias for GPS FREQ_3
			SAT_IFB_C3, SAT_IFB_C4, SAT_IFB_C5,
			GPS_REC_IFB_C3, GPS_REC_IFB_C4, GPS_REC_IFB_C5,
			GAL_REC_IFB_C3, GAL_REC_IFB_C4, GAL_REC_IFB_C5,
			BDS_REC_IFB_C3, BDS_REC_IFB_C4, BDS_REC_IFB_C5,
			QZS_REC_IFB_C3, QZS_REC_IFB_C4, QZS_REC_IFB_C5,
			GPS_IFB_C3, GPS_IFB_C4, GPS_IFB_C5,
			GAL_IFB_C3, GAL_IFB_C4, GAL_IFB_C5,
			BDS_IFB_C3, BDS_IFB_C4, BDS_IFB_C5,
			QZS_IFB_C3, QZS_IFB_C4, QZS_IFB_C5,			                                    // inter frequency 3 bias for gnss
			CLK_ICB, CLUSTERB,                                          // initial clock bias, cluster-dependent bias
			AMB_IF, AMB13_IF, AMB_WL, AMB_L1, AMB_L2, AMB_L3, AMB_L4, AMB_L5,             // ambiguities for indiv. freq. (number indicates freq not band)
			FCBS_IF, FCBS_L1, FCBS_L2, FCBS_L3, FCBS_L4, FCBS_L5,        // satellite fractional cycle biases for indiv. freq.
			FCBR_IF, FCBR_L1, FCBR_L2, FCBR_L3, FCBR_L4, FCBR_L5,	     // receiver  fractional cycle biases for indiv. freq.
			GLO_ISB, GLO_IFCB, GLO_IFPB, GAL_ISB, BDS_ISB, QZS_ISB,      // multi-GNSS
			IFB_GPS, IFB_BDS,  IFB_QZS, IFB_GAL, IFB_GAL_2, IFB_GAL_3,      // inter frequency bias for multi-frquency
			IFB_BDS_2, IFB_BDS_3,
			ATT_X, ATT_Y, ATT_Z,
			eb_X, eb_Y, eb_Z,
			db_X, db_Y, db_Z,
			CAM_ATT_X, CAM_ATT_Y, CAM_ATT_Z,
			CAM_CRD_X, CAM_CRD_Y, CAM_CRD_Z,
			EXTRINSIC_ATT_X, EXTRINSIC_ATT_Y, EXTRINSIC_ATT_Z,
			EXTRINSIC_CRD_X, EXTRINSIC_CRD_Y, EXTRINSIC_CRD_Z,
			LIDAR_ATT_X, LIDAR_ATT_Y, LIDAR_ATT_Z,
			LIDAR_CRD_X, LIDAR_CRD_Y, LIDAR_CRD_Z,
			GLO_IFB,													 // inter frequency bias for glonass
			VTEC, DCB_REC, DCB_SAT,										// for DCB estimator
			P1P2G_REC, P1P2E_REC, P1P2R_REC, P1P2C_REC,                 // GNSS-specific receiver code DCB P1-P2
			VEL_X, VEL_Y, VEL_Z,										 // velocity
	        ACC_X, ACC_Y, ACC_Z,										 // acceleration
			CLK_RAT,													 // satellite clock speed
			SCALE,GEOCX,GEOCY,GEOCZ,									 // ERP parameter add by lijie
			XPOLE, YPOLE, DXPOLE, DYPOLE, UT1, DUT1,						 // EOP parameter
			PXSAT, PYSAT, PZSAT, VXSAT, VYSAT, VZSAT,	                     // satellite parameter
			SR_scale, DRAG_c, 
			EMP_Sc1, EMP_Cc1, EMP_Sa1, EMP_Ca1, EMP_Sr1, EMP_Cr1,
			EMP_Sc2, EMP_Cc2, EMP_Sa2, EMP_Ca2, EMP_Sr2, EMP_Cr2,
			EMP_Ac, EMP_Bc, EMP_Aa, EMP_Ba, EMP_Ar, EMP_Br,
			ECOM_D0, ECOM_Dc, ECOM_Ds, ECOM_Y0, ECOM_Yc, ECOM_Ys, ECOM_B0, ECOM_Bc, ECOM_Bs, // ECOM Force parameter
			Xs,
			ECOM2_D0,ECOM2_D2c,ECOM2_D2s,ECOM2_D4c,ECOM2_D4s,			     // ECOM2 Force parameters
			ECOM2_Y0,ECOM2_X0,ECOM2_X1c,ECOM2_X1s,							 // ECOM2
			IFBR_RAT1, IFBR_RAT2,
			RB,                                                            //range biase for slr
			SION_3,
			NO_DEF
		};

		t_gpar(string site, t_type t, int i, string p, bool remove = true);
		t_gpar();
		virtual ~t_gpar();

		bool  operator==(const t_gpar&) const;
		bool  operator <(const t_gpar&) const;
		bool  operator >(const t_gpar&) const;
		t_gpar operator-(const t_gpar&) const;
		t_gpar operator+(const t_gpar&) const;

		void setTime(const t_gtime&, const t_gtime&);

		double partial(t_gsatdata&, t_gtime&, t_gtriple, t_gobs& gobs);
		double partial_doppler(t_gsatdata& satData, t_gtriple& groundXYZ, t_gtriple& groundVEL); // set conf for doppler (by zhshen)
		double partial_eop(t_gsatdata&);
		double partial_orb(t_gsatdata&);
		double partial_orb_leo(t_gsatdata&);


		void   value(double val) { _value = val; }
		double value()const { return _value; }

		void   apriori(double apr) { _apriori = apr; }
		double apriori()const { return _apriori; }

		t_type parType;
		int index;
		string prn;
		string site;
		t_gtime beg;
		t_gtime end;
		t_gtime stime;
		double  aprval;
		double pred;
		double  zhd; // for ztd par
		double  zwd;
		//for amb
		bool amb_ini = false;    //add for amb to be initialized
		// for LEO
		bool    lremove;  // add for LEO satellite(default: false)
		bool    lPWC;     // add for LEO PWC parameters(default: false)
		int     nPWC;     // number of parameters for this kind of PWC parameter
	    int     piecewise_time;// time 
		//
		string str_type() const;
		t_gparhead get_head() const;
		t_gtimearc get_timearc() const;

		void setMF(ZTDMPFUNC MF);
		void setMF(GRDMPFUNC MF);

	protected:
		double    _value;   // value
		ZTDMPFUNC _mf_ztd;      // mapping function for ZTD
		GRDMPFUNC _mf_grd;      // mapping function for GRD 
		double    _apriori;
		void      _getmf(t_gsatdata& satData, const t_gtriple& crd, const t_gtime& epoch, double& mfw, double& mfh, double& dmfw, double& dmfh);
	};

	class BASE_LIBRARY_EXPORT t_gparhead
	{
	public:
		t_gparhead(t_gpar::t_type type, string site, string sat);
		t_gparhead(const t_gparhead& Other);
		~t_gparhead();

		bool operator==(const t_gparhead& Other) const;
		bool operator<(const t_gparhead& Other) const;
		bool operator<=(const t_gparhead& Other) const;
		bool operator>(const t_gparhead& Other) const;
		bool operator>=(const t_gparhead& Other) const;

		t_gpar::t_type type;
		string	site;
		string	sat;
	};

	string BASE_LIBRARY_EXPORT gpar2str(const t_gpar& par);
	string BASE_LIBRARY_EXPORT ptype2str(const t_gpar::t_type& par);
	t_gpar BASE_LIBRARY_EXPORT str2gpar(const string& str_par);

} // namespace

#endif
