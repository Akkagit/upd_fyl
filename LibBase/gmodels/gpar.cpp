
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 
-*/

#include <iostream>
#include <cmath>

#include "newmat/newmat.h"

#include "gmodels/gpar.h"
#include "gmodels/ggmf.h"

using namespace std;

namespace gnut {   

	const string t_gpar::PAR_STR_SEP = "_";
	
	const vector<string> t_gpar::PAR_STR =
	{
			"CRD_X", "CRD_Y", "CRD_Z",                                        // coordinates
			"TRP", "GRD_N", "GRD_E", "SION", "VION",                              // atmospheric parameters
			"CLK", "CLK_SAT",  "G15CLK", "E15CLK", "C15CLK", "J15CLK",                                                     // clocks
			"G15CLK_SAT", "E15CLK_SAT", "C15CLK_SAT", "J15CLK_SAT",             // satellite clocks for 13 combination (add by xiongyun)
			"CLK_G", "CLK_E", "CLK_C", "CLK_R", "CLK_J",                          // receiver clocks for different systems (add by xiongyun)
			"IFB_C3", "IFB_C4", "IFB_C5", "IFCB_F3", "IFCB_F4", "IFCB_F5",          // inter-freq. code biases for FREQ_3, FREQ_4, FREQ_5, inter-freq. clock bias for GPS FREQ_3
			"SAT_IFB_C3", "SAT_IFB_C4", "SAT_IFB_C5",
			"GPS_REC_IFB_C3", "GPS_REC_IFB_C4", "GPS_REC_IFB_C5",
			"GAL_REC_IFB_C3", "GAL_REC_IFB_C4", "GAL_REC_IFB_C5",
			"BDS_REC_IFB_C3", "BDS_REC_IFB_C4", "BDS_REC_IFB_C5",
			"QZS_REC_IFB_C3", "QZS_REC_IFB_C4", "QZS_REC_IFB_C5",
			"GPS_IFB_C3", "GPS_IFB_C4", "GPS_IFB_C5",
			"GAL_IFB_C3", "GAL_IFB_C4", "GAL_IFB_C5",
			"BDS_IFB_C3", "BDS_IFB_C4", "BDS_IFB_C5",
			"QZS_IFB_C3", "QZS_IFB_C4", "QZS_IFB_C5",
			"CLK_ICB", "CLUSTERB",                                          // initial clock bias, cluster-dependent bias
			"AMB_IF", "AMB13_IF", "AMB_WL","AMB_L1", "AMB_L2", "AMB_L3", "AMB_L4", "AMB_L5",             // ambiguities for indiv. freq. (number indicates freq not band)
			"FCBS_IF", "FCBS_L1", "FCBS_L2", "FCBS_L3", "FCBS_L4", "FCBS_L5",       // satellite fractional cycle biases for indiv. freq.
			"FCBR_IF", "FCBR_L1", "FCBR_L2", "FCBR_L3", "FCBR_L4", "FCBR_L5",	     // receiver fractional cycle biases for indiv. freq.
			"GLO_ISB", "GLO_IFCB", "GLO_IFPB", "GAL_ISB", "BDS_ISB", "QZS_ISB",     // multi-GNSS
			"IFB_GPS","IFB_BDS","IFB_QZS","IFB_GAL","IFB_GAL_2", "IFB_GAL_3",
			"IFB_BDS_2","IFB_BDS_3",
			"ATT_X", "ATT_Y", "ATT_Z",
			"eb_X", "eb_Y", "eb_Z",
			"db_X", "db_Y", "db_Z",
			"CAM_ATT_X", "CAM_ATT_Y", "CAM_ATT_Z",
			"CAM_CRD_X", "CAM_CRD_Y", "CAM_CRD_Z",
			"EXTRINSIC_ATT_X", "EXTRINSIC_ATT_Y", "EXTRINSIC_ATT_Z",
			"EXTRINSIC_CRD_X", "EXTRINSIC_CRD_Y", "EXTRINSIC_CRD_Z",
			"GLO_IFB",													        // inter frequency bias for glonass
			"VTEC", "DCB_REC", "DCB_SAT",										 // for DCB estimator
			"P1P2G_REC", "P1P2E_REC", "P1P2R_REC", "P1P2C_REC",                  // GNSS-specific receiver code DCB P1-P2
			"VEL_X", "VEL_Y", "VEL_Z",										     // velocity
			"ACC_X", "ACC_Y", "ACC_Z",										     // acceleration
			"CLK_RAT",													         // satellite clock speed
			"SCALE","GEOCX","GEOCY","GEOCZ",									 // ERP parameter add by lijie
			"XPOLE", "YPOLE", "DXPOLE", "DYPOLE", "UT1", "DUT1",						 // EOP parameter
			"PXSAT", "PYSAT", "PZSAT", "VXSAT", "VYSAT", "VZSAT",	                     // satellite parameter
			"SR_scale", "DRAG_c", 
			"EMP_Sc1", "EMP_Cc1", "EMP_Sa1", "EMP_Ca1", "EMP_Sr1", "EMP_Cr1",
			"EMP_Sc2", "EMP_Cc2", "EMP_Sa2", "EMP_Ca2", "EMP_Sr2", "EMP_Cr2",
			"EMP_Ac", "EMP_Bc", "EMP_Aa", "EMP_Ba", "EMP_Ar", "EMP_Br",
			"D0", "Dc", "Ds", "Y0", "Yc", "Ys", "B0", "Bc", "Bs",							 // ECOM Force parameter
			"Xs",	
			"ECOM2_D0","ECOM2_D2c","ECOM2_D2s","ECOM2_D4c","ECOM2_D4s",						 // ECOM2
			"ECOM2_Y0","ECOM2_X0","ECOM2_X1c","ECOM2_X1s",									 // ECOM2
			"IFBR_RAT1" , "IFBR_RAT2",
			"RB",
			"SION_3",
	};

	t_gparhead::t_gparhead(t_gpar::t_type type, string site, string sat) :
		type(type),
		site(site),
		sat(sat)
	{
	}
	t_gparhead::t_gparhead(const t_gparhead & Other):
		type(Other.type),
		sat(Other.sat),
		site(Other.site)
	{
	}
	t_gparhead::~t_gparhead()
	{
	}
	bool t_gparhead::operator==(const t_gparhead & Other) const
	{
		if (this->type == Other.type &&
			this->sat == Other.sat &&
			this->site == Other.site)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	bool t_gparhead::operator<(const t_gparhead & Other) const
	{
		if (this->type < Other.type ||
			this->type == Other.type && this->site < Other.site ||
			this->type == Other.type && this->site == Other.site && this->sat < Other.sat)
		{
			return true;
		}
		else {
			return false;
		}
	}
	bool t_gparhead::operator<=(const t_gparhead & Other) const
	{
		if (*this == Other || *this < Other) {
			return true;
		}
		else {
			return false;
		}
	}
	bool t_gparhead::operator>(const t_gparhead & Other) const
	{
		if (*this <= Other) {
			return false;
		}
		else {
			return true;
		}
	}
	bool t_gparhead::operator>=(const t_gparhead & Other) const
	{
		if (*this < Other) {
			return false;
		}
		else {
			return true;
		}
	}
	t_gtimearc::t_gtimearc(const t_gtime& beg, const t_gtime& end):
		begin(beg),
		end(end)
	{
	}
	t_gtimearc::~t_gtimearc()
	{
	}
	bool t_gtimearc::operator!=(const t_gtimearc & Other) const
	{
		return !(*this == Other);
	}
	bool t_gtimearc::operator==(const t_gtimearc & Other) const
	{
		if (this->begin == Other.begin && this->end == Other.end) {
			return true;
		}
		else {
			return false;
		}
	}
	bool t_gtimearc::operator<(const t_gtimearc & Other) const
	{
		if (this->begin < Other.begin ||
			this->begin == Other.begin && this->end < Other.end)
		{
			return true;
		}
		else {
			return false;
		}
	}
	bool t_gtimearc::operator<=(const t_gtimearc & Other) const
	{
		if (*this == Other || *this < Other) {
			return true;
		}
		else {
			return false;
		}
	}
	bool t_gtimearc::operator>(const t_gtimearc & Other) const
	{
		if (*this <= Other) {
			return false;
		}
		else {
			return true;
		}

	}
	bool t_gtimearc::operator>=(const t_gtimearc & Other) const
	{
		if (*this < Other) {
			return false;
		}
		else {
			return true;
		}
	}
	bool t_gtimearc::inside(const t_gtimearc & Other) const
	{
		if (this->begin >= Other.begin && this->end <= Other.end)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

// constructor
// --------------------------------------------------------
t_gpar::t_gpar(string site, t_type t, int i, string p, bool remove)
{
   beg = FIRST_TIME;
   end = LAST_TIME;

   this->site = site;
   lremove = remove;  //added for LEO
   parType = t;
   index = i;
   prn = p;
   lPWC = false;
   amb_ini = false;
   //indPWC = 0;
   nPWC = 1;
   value( 0.0 );
   apriori( 0.0 );
}

t_gpar::t_gpar()
{
   beg = FIRST_TIME;
   end = LAST_TIME;   
}

// t_gpar destructor
// ---------------------------------------------------
t_gpar::~t_gpar()
{
}

// setting mapping function for ZTD
// ---------------------------------------------------
void t_gpar::setMF(ZTDMPFUNC MF)
{
  _mf_ztd = MF;   
}

// setting mapping function for GRD
// ---------------------------------------------------
void t_gpar::setMF(GRDMPFUNC MF)
{
  _mf_grd = MF;
}
   

// Partial derivatives 
// ----------------------------------------------------
double t_gpar::partial(t_gsatdata& satData, t_gtime& epoch, t_gtriple ground, t_gobs& gobs)
{
  double mfw, dmfw, mfh, dmfh;
  mfw = dmfw = mfh = dmfh = 0.0;   
  // modified by glfeng : satData.rotmat() is the rotmat Matrix of TR_epoch
  // change to partial_panda because earth rotation matirx
 // t_gtriple unit = satData.dloudx();
  //ColumnVector partial_panda = -(satData.rotmat().t() * unit.crd_cvect());
  switch( parType ){
    case CRD_X: if (satData.site() == this->site)
    {
		//cout << this->str_type() << " " << fixed << setprecision(10) << setw(25) << (value() - satData.satcrd().crd(0)) / satData.rho() << endl;
		return (value() - satData.satcrd().crd(0)) / satData.rho();
	   // cout << this->str_type() << " " << fixed << setprecision(10) << setw(25) << partial_panda(1) << endl;
	    //return partial_panda(1);
    } else return 0.0;
    case CRD_Y: if (satData.site() == this->site)
    {
		//cout << this->str_type() << " " << fixed << setprecision(10) << setw(25) << (value() - satData.satcrd().crd(1)) / satData.rho() << endl;
		return (value() - satData.satcrd().crd(1)) / satData.rho();
	  /*  cout << this->str_type() << " " << fixed << setprecision(10) << setw(25) << partial_panda(2) << endl;
	    return partial_panda(2);*/
    }else return 0.0;
    case CRD_Z: if (satData.site() == this->site)
    {
		//cout << this->str_type() << " " << fixed << setprecision(10) << setw(25) << (value() - satData.satcrd().crd(2)) / satData.rho() << endl;
		return (value() - satData.satcrd().crd(2)) / satData.rho();
	    //cout << this->str_type() << " " << fixed << setprecision(10) << setw(25) << partial_panda(3) << endl;
	    //return partial_panda(3);
    }else return 0.0;
  case CLK:
  {
	  if (satData.site() == this->site) return 1.0;
	  else return 0.0;
  }
  case CLK_SAT:
  {
	  if (satData.sat() == this->prn) return -1.0; 
	  else return 0.0;
  }
    case IFCB_F3 : if (t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_3 
                       && prn == satData.sat()) return -1.0; else return 0.0;
    case IFCB_F4 : if (t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_4
                       && prn == satData.sat()) return -1.0; else return 0.0;
    case IFCB_F5 : if (t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_5 
                       && prn == satData.sat()) return -1.0; else return 0.0;    
    case IFB_C3 : if (gobs.is_code() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_3) return 1.0; else return 0.0;
    case IFB_C4 : if (gobs.is_code() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_4) return 1.0; else return 0.0;
    case IFB_C5 : if (gobs.is_code() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_5) return 1.0; else return 0.0;    
    case TRP   :
		if (satData.site() == this->site) {
			_getmf(satData, ground, epoch, mfw, mfh, dmfw, dmfh);
#ifdef DEBUG
			cout << this->str_type() << " " << scientific << setprecision(10) << setw(25) << mfw << endl;

#endif // DEBUG

			return mfw;
		}
		else return 0.0;
    case SION  :
      {
		if (satData.site() == this->site) 
		{
			// jdhuang : change from G01_F
			double f1 = satData.frequency(t_gsys::band_priority(satData.gsys(), FREQ_1));
			double fk = satData.frequency(gobs.band());
			double alfa = 0.0;
			if (gobs.is_phase() && prn == satData.sat()) { alfa = -(f1*f1) / (fk*fk); }
			if (gobs.is_code() && prn == satData.sat()) { alfa = (f1*f1) / (fk*fk); }
			return alfa;
		}
		else return 0.0;
      }
    case VION  :
      {
		if (satData.site() == this->site) 
		{
			double f1 = satData.frequency(t_gsys::band_priority(satData.gsys(), FREQ_1));
			double fk = satData.frequency(gobs.band());
			double mf = 1.0 / sqrt(1.0 - pow(R_SPHERE / (R_SPHERE + 450000.0) * sin(G_PI / 2.0 - satData.ele()), 2));
			double alfa = 0.0;
			if (gobs.is_phase() && prn == satData.sat()) { alfa = -(f1*f1) / (fk*fk); }
			if (gobs.is_code() && prn == satData.sat()) { alfa = (f1*f1) / (fk*fk); }
			//cout << "Partial " << satData.sat() << " " << prn << "  " << gobs2str(gobs.gobs()) << "  " << gobs.band() << "  " << alfa << "  " << f1 << "  " << fk << endl;
			return alfa * mf;
		}
		else return 0.0;
      }
    case P1P2G_REC :
      {
		if (satData.site() == this->site) {
			double f1 = G01_F;
			double fk = satData.frequency(gobs.band());
			double alfa = (f1*f1) / (fk*fk);
			double beta = (G02_F*G02_F) / (G01_F*G01_F - G02_F * G02_F);
			FREQ_SEQ freq = t_gsys::band2freq(satData.gsys(), gobs.band());
			if (satData.gsys() == GPS && gobs.is_code() && (freq == FREQ_1 || freq == FREQ_2)) { return -alfa * beta; }
			else return 0.0;
		}
		else return 0.0;
      }
    case P1P2E_REC :
      {
		if (satData.site() == this->site) {
			double f1 = E01_F;
			double fk = satData.frequency(gobs.band());
			double alfa = (f1*f1) / (fk*fk);
			double beta = (E05_F*E05_F) / (E01_F*E01_F - E05_F * E05_F);
			FREQ_SEQ freq = t_gsys::band2freq(satData.gsys(), gobs.band());
			if (satData.gsys() == GAL && gobs.is_code() && (freq == FREQ_1 || freq == FREQ_2)) { return -alfa * beta; }
			else return 0.0;
		}
		else return 0.0;
      }    
    case GRD_N :
		if (satData.site() == this->site) {
			_getmf(satData, ground, epoch, mfw, mfh, dmfw, dmfh);
			if (_mf_grd == CHEN_HERRING) {
				double sinel = sin(satData.ele());
				double tanel = tan(satData.ele());
				double cosaz = cos(satData.azi());
				return (1.0 / (sinel*tanel + 0.0032)) * cosaz;
			}
			else if (_mf_grd == TILTING) {
				double cosaz = cos(satData.azi());
				return dmfw * cosaz;
			}
			else if (_mf_grd == BAR_SEVER) {
				double tanel = tan(satData.ele());
				double cosaz = cos(satData.azi());
				return mfw * (1.0 / tanel) * cosaz;
			}
			else cerr << "Grad N mapping function is not set up correctly!!!" << endl;
		}
		else return 0.0;
    case GRD_E :
		if (satData.site() == this->site) {
			_getmf(satData, ground, epoch, mfw, mfh, dmfw, dmfh);
			if (_mf_grd == CHEN_HERRING) {
				double sinel = sin(satData.ele());
				double tanel = tan(satData.ele());
				double sinaz = sin(satData.azi());
				return (1.0 / (sinel*tanel + 0.0032)) * sinaz;
			}
			else if (_mf_grd == TILTING) {
				double sinaz = sin(satData.azi());
				return dmfw * sinaz;
			}
			else if (_mf_grd == BAR_SEVER) {
				double tanel = tan(satData.ele());
				double sinaz = sin(satData.azi());
				return mfw * (1 / tanel) * sinaz;
			}
			else cerr << "Grad E mapping function is not set up correctly!!!" << endl;
		}
		else return 0.0;
     
    case AMB_IF  :
          if (satData.site() == this->site && gobs.is_phase() && prn == satData.sat()) return 1.0;
          else return 0.0;
    case AMB_L1  :
          if (satData.site() == this->site && gobs.is_phase() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_1 && prn == satData.sat()) return 1.0;
          else return 0.0;
    case AMB_L2  :
          if (satData.site() == this->site && gobs.is_phase() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_2 && prn == satData.sat()) return 1.0;
          else return 0.0;  
    case AMB_L3  :
          if (satData.site() == this->site && gobs.is_phase() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_3 && prn == satData.sat()) return 1.0;
          else return 0.0;     
    case AMB_L4  :
          if (satData.site() == this->site && gobs.is_phase() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_4 && prn == satData.sat()) return 1.0;
          else return 0.0;
    case AMB_L5  :
          if (satData.site() == this->site && gobs.is_phase() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_5 && prn == satData.sat()) return 1.0;
          else return 0.0;        
    case GLO_ISB :
          if (satData.site() == this->site && satData.gsys() == GLO) return 1.0;
          else return 0.0;
   case GLO_IFCB :
          if (satData.site() == this->site && !gobs.is_phase() && satData.gsys() == GLO && prn == satData.sat()) return 1.0;
          else return 0.0;
   case GLO_IFPB :
          if (satData.site() == this->site && gobs.is_phase() && satData.gsys() == GLO && prn == satData.sat()) return 1.0;
          else return 0.0;     
   case GLO_IFB:
	   if (satData.site() == this->site && satData.gsys() == GLO && prn == satData.sat()) return 1.0;
	   else return 0.0;
    case GAL_ISB :
          if (satData.site() == this->site && satData.gsys() == GAL) return 1.0;
          else return 0.0;
    case BDS_ISB :
          if (satData.site() == this->site && satData.gsys() == BDS) return 1.0;
          else return 0.0;
    case QZS_ISB :
          if (satData.site() == this->site && satData.gsys() == QZS) return 1.0;
          else return 0.0;
  }

  if (parType >= XPOLE && parType <= DUT1) {
	  cout << this->str_type() << " " << scientific << setprecision(15) << setw(30)<<  partial_eop(satData) << endl;
	  return partial_eop(satData);
  }
  else if (parType >= PXSAT && parType <= ECOM2_X1s) {
	  if (satData.sat() == prn) {
		  cout << this->str_type() << " " << scientific << setprecision(15) << setw(30) << partial_orb(satData) << endl;
		  return partial_orb(satData);
	  }
	  else if(satData.site() == site) {
		//cout << this->str_type() << " " << scientific << setprecision(15) << setw(30) << partial_orb_leo(satData) << endl;
		  return partial_orb_leo(satData);
	  }
	  else 
		  return 0.0;
  }

  return 0.0;
}

double t_gpar::partial_doppler(t_gsatdata& satData, t_gtriple& groundXYZ, t_gtriple& groundVEL)
{
	t_gtriple satcrd = satData.satcrd();
	t_gtriple satvel = satData.satvel();
	ColumnVector cSat = satcrd.crd_cvect();
	ColumnVector vSat = satvel.crd_cvect();
	ColumnVector cRec = groundXYZ.crd_cvect();
	ColumnVector vRec = groundVEL.crd_cvect();
	ColumnVector conf_crd = dotproduct(cSat - cRec, vSat - vRec)*(cSat - cRec) / pow(satData.rho(), 3);
	conf_crd -= (vSat - vRec) / satData.rho();
	ColumnVector e = (cSat - cRec) / satData.rho();
	satData._conf_crd.set(conf_crd);
	satData._e.set(e);
	switch (parType)
	{
	case CRD_X: return conf_crd(1);
	case CRD_Y: return conf_crd(2);
	case CRD_Z: return conf_crd(3);
	case VEL_X: return -e(1);
	case VEL_Y: return -e(2);
	case VEL_Z: return -e(3);
	case CLK_RAT: return 1.0;
	}
	return 0.0;
}

double t_gpar::partial_eop(t_gsatdata& satdata)
{

	// get REC crd in trs
	ColumnVector xyz_rec_trs(3);
	xyz_rec_trs = satdata.reccrd().crd_cvect();
	// get unit_vector in TRS
	ColumnVector unit_rec2sat(3);
	// get
	unit_rec2sat = (satdata.satcrd() - satdata.reccrd()).crd_cvect();
	// change TRS to CRS
	unit_rec2sat = (satdata.rotmat())*unit_rec2sat;
	unit_rec2sat = unit_rec2sat / unit_rec2sat.NormFrobenius();


	// testZHJ
	//cout << setw(20) << setprecision(13) << fixed << xyz_rec_trs;

	// get drotdxpole
	Matrix drdxpole = satdata.drdxpole();
	//drdxpole = drdxpole * xyz_rec_trs;
	drdxpole = drdxpole / RAD2SEC;
	// get drotdypole
	Matrix drdypole = satdata.drdypole();
	//drdypole = drdypole * xyz_rec_trs;
	drdypole = drdypole / RAD2SEC;
	// get drotdut1
	Matrix drdut1 = satdata.drdut1();
	//drdut1 = drdut1 * xyz_rec_trs;
	drdut1 = drdut1 / RAD2TSEC;
	double diff = (satdata.recTime() - beg) / 86400.0;


	// debug ZHJ
	//cout << setw(20) << setprecision(13) << fixed << drdxpole;
	//cout << setw(20) << setprecision(13) << fixed << drdypole;
	//cout << setw(20) << setprecision(13) << fixed << drdut1;
	//cout << setw(20) << setprecision(13) << fixed << unit_rec2sat;


	// dxdxpole dxdypole dxdut1
	Matrix dxdpoleut1(3, 6);
	// dx/xpole
	dxdpoleut1.Column(1) = drdxpole;
	// dx/dxpole
	dxdpoleut1.Column(2) = drdxpole * diff;
	// dx/ypole
	dxdpoleut1.Column(3) = drdypole;
	// dx/dypole
	dxdpoleut1.Column(4) = drdypole * diff;
	// dx/ut1
	dxdpoleut1.Column(5) = drdut1;
	// dx/dut1
	dxdpoleut1.Column(6) = drdut1 * diff;

	Matrix dr_derp(1, 6);
	dr_derp = unit_rec2sat.t() * dxdpoleut1;


	//bool debug = false;
	//if (debug) {
	//	for (int i = 0; i < 6; i++) {
	//		dr_derp(1,i + 1) = (*satdata.panda_info)[satdata.site()][satdata.sat()].epo_amat[i];
	//		if (i < 4) {
	//			dr_derp(1, i + 1) = -dr_derp(1, i + 1);
	//		}
	//	}
	//}

	switch (parType)
	{
	case XPOLE:
		return -dr_derp(1, 1);
		break;
	case YPOLE:
		return -dr_derp(1, 3);
		break;
	case DXPOLE:
		return - dr_derp(1, 2);
		break;
	case DYPOLE:
		return -dr_derp(1, 4);
		break;
	case UT1:
		return dr_derp(1, 5);
		break;
	case DUT1:
		//cout << setw(20) << setprecision(13) << fixed << satdata.drdut1();
		//cout << setw(20) << setprecision(13) << fixed << drdut1;
		//cout << setw(20) << setprecision(13) << fixed << dxdpoleut1;
		//cout << setw(20) << setprecision(13) << dr_derp << endl;
		return dr_derp(1, 6);
		break;
	default:
		return 0.0;
		break;
	}
}

double t_gpar::partial_orb(t_gsatdata& satdata)
{
	Matrix funct = satdata.orbfunct();

	//cout << setw(20) << setprecision(13) << funct << endl;

	int index_begin = satdata.satindex();

	// get unit_vector in TRS
	ColumnVector unit_rec2sat(3);
	// get
	unit_rec2sat = (satdata.satcrd() - satdata.reccrd()).crd_cvect();
	// change TRS to CRS
	unit_rec2sat = (satdata.rotmat())*unit_rec2sat;
	unit_rec2sat = unit_rec2sat / unit_rec2sat.NormFrobenius();

	// order: P X0 Y0 Z0 VX0..VZ0 solarPar

	//cout << "Before " << endl << setw(20) << setprecision(13) << scientific << funct << endl;
	funct = funct * unit_rec2sat;
	funct = funct * 1E3;

	//cout << "After " << endl << setw(20) << setprecision(13)<< scientific << funct << endl;
	//cout << "After " << endl << setw(20) << setprecision(13) << scientific << unit_rec2sat << endl;

	//bool debug = false;
	//if (debug) {
	//	for (int i = 0; i < 11; i++) {
	//		funct(i + 1, 1) = (*satdata.panda_info)[satdata.site()][satdata.sat()].orb_amat[i];
	//	}
	//}

	return funct(index - index_begin + 1, 1);

}


double t_gpar::partial_orb_leo(t_gsatdata &satdata)
{
	Matrix funct = satdata.orbfunct();
	int index_begin = satdata.satindex();

	// get unit_vector in crs
	ColumnVector unit_rec2sat(3);
	// get
	unit_rec2sat= -satdata.sat2reccrs().crd_cvect();
	
	unit_rec2sat = unit_rec2sat / unit_rec2sat.NormFrobenius();
	
	/*cout << endl;
	cout << satdata.epoch().mjd()<<"  "<< satdata.epoch().sod() <<"   u2s  "<<setw(30)<<setprecision(15)<<scientific<< unit_rec2sat << endl;*/
	// order: P X0 Y0 Z0 VX0..VZ0 solarPar
	funct = funct * unit_rec2sat;
	funct = funct * 1E3;

	//cout << setw(20) << setprecision(13)<< scientific << funct << endl;
	return funct(index - index_begin + 1, 1);
	/*switch (parType)
	{
	case PXSAT:case PYSAT:case PZSAT:
	case VXSAT:case VYSAT:case VZSAT:
		return funct(index - index_begin + 1, 1);
		break;
	case DRAG_c:
		return funct(7, 1); 
		break;
	case EMP_Sc1:
		return funct(8, 1); 
		break;
	case EMP_Cc1:
		return funct(9, 1); 
		break;
	case EMP_Sa1:
		return funct(10, 1); 
		break;
	case EMP_Ca1:
		return funct(11, 1); 
		break;
	case EMP_Sr1:
		return funct(12, 1); 
		break;
	case EMP_Cr1:
		return funct(13, 1); 
		break;
	default:
		return 0;

	}*/
		

	
}


// Operators for t_gpar
// -------------------------------------------------
bool t_gpar::operator==(const t_gpar& par) const
{
   if ( parType == par.parType      &&
      //  prn.compare(par.prn)   == 0 &&
      //  site.compare(par.site) == 0 &&
   beg == par.beg &&
   end == par.end ){
      return true;
   } else {
      return false;
   }
   
}

// Doplnil gabo aby par mohla byt pouzita ako klic v mape
bool t_gpar::operator <(const t_gpar& par) const
{
    if (( parType <  par.parType )                                        ||
        ( parType == par.parType && prn <  par.prn )                      ||
        ( parType == par.parType && prn == par.prn && site <  par.site )  ||
        ( parType == par.parType && prn == par.prn && site == par.site && beg < par.beg ) ||
        ( parType == par.parType && prn == par.prn && site == par.site && beg == par.beg && end < par.end )
        ){
       return true;
    }
    return false;
}

bool t_gpar::operator >(const t_gpar& par) const
{
	if ((parType > par.parType) ||
		(parType == par.parType && prn < par.prn) ||
		(parType == par.parType && prn == par.prn && site < par.site) ||
		(parType == par.parType && prn == par.prn && site == par.site && beg < par.beg) ||
		(parType == par.parType && prn == par.prn && site == par.site && beg == par.beg && end < par.end)
		) {
		return true;
	}
	return false;
}

t_gpar t_gpar::operator-(const t_gpar& p) const
{
   t_gpar par = (*this);
   par.value( value() - p.value() );
   return par;
}

t_gpar t_gpar::operator+(const t_gpar& p) const
{
   t_gpar par = (*this);
   par.value( value() + p.value() );
   return par;
}

// Setting begin and end time for validity object
// -------------------------------------------------
void t_gpar::setTime(const t_gtime& t1, const t_gtime& t2)
{
   this->beg = t1;
   this->end = t2;
}

// get data type
// ----------
string t_gpar::str_type() const
{
  // jdhuang : change to pretty recover file
  string type;
  switch( parType ){
   case  CRD_X    :  type = "CRD_X";     break;
   case  CRD_Y    :  type = "CRD_Y";     break;
   case  CRD_Z    :  type = "CRD_Z";     break;
   case  TRP      :  type = "TRP"  ;     break;
   case  VEL_X:  type = "VEL_X";     break;
   case  VEL_Y:  type = "VEL_Y";     break;
   case  VEL_Z:  type = "VEL_Z";     break;
   case  CLK_RAT:  type = "CLK_RAT";     break;
   case  SION     :  type = "SION_" + prn ;     break;
   case  VION     :  type = "VION_" + prn ;     break;  
   case  CLK      :  type = "CLK";              break;
   case  CLK_G    :  type = "CLK_G";            break;
   case  CLK_C    :  type = "CLK_C";            break;
   case  CLK_E    :  type = "CLK_E";            break;
   case  CLK_R    :  type = "CLK_R";            break;
   case  CLK_J    :  type = "CLK_J";            break;
   case  G15CLK   :  type = "G15CLK";           break;
   case  E15CLK   :  type = "E15CLK";           break;
   case  C15CLK   :  type = "C15CLK";           break;
   case  J15CLK   :  type = "J15CLK";           break;
   case  IFCB_F3  :  type = "IFCB_F3_" + prn;     break;
   case  IFCB_F4  :  type = "IFCB_F4_" + prn;     break;
   case  IFCB_F5  :  type = "IFCB_F5_" + prn;     break;    
   case  IFB_C3   :  type = "IFB_F3_" + prn;     break;
   case  IFB_C4   :  type = "IFB_F4_" + prn;     break;
   case  IFB_C5   :  type = "IFB_F5_" + prn;     break;    
   case  GPS_REC_IFB_C3:  type = "GPS_REC_IFB_C3";     break;
   case  GPS_REC_IFB_C4:  type = "GPS_REC_IFB_C4";     break;
   case  GPS_REC_IFB_C5:  type = "GPS_REC_IFB_C5";     break;
   case  GAL_REC_IFB_C3:  type = "GAL_REC_IFB_C3";     break;
   case  GAL_REC_IFB_C4:  type = "GAL_REC_IFB_C4";     break;
   case  GAL_REC_IFB_C5:  type = "GAL_REC_IFB_C5";     break;
   case  BDS_REC_IFB_C3:  type = "BDS_REC_IFB_C3";     break;
   case  BDS_REC_IFB_C4:  type = "BDS_REC_IFB_C4";     break;
   case  BDS_REC_IFB_C5:  type = "BDS_REC_IFB_C5";     break;
   case  QZS_REC_IFB_C3:  type = "QZS_REC_IFB_C3";     break;
   case  QZS_REC_IFB_C4:  type = "QZS_REC_IFB_C4";     break;
   case  QZS_REC_IFB_C5:  type = "QZS_REC_IFB_C5";     break;
   case  SAT_IFB_C3   :  type = "GPS_SAT_IFB_C3" + prn;     break;    
   case  SAT_IFB_C4   :  type = "GPS_SAT_IFB_C4" + prn;     break;    
   case  SAT_IFB_C5   :  type = "GPS_SAT_IFB_C5" + prn;     break;  
   case  CLK_SAT  :  type = "CLK_SAT_" + prn;   break;
   case  G15CLK_SAT: type = "G15CLK_SAT_" + prn;   break;
   case  E15CLK_SAT: type = "E15CLK_SAT_" + prn;   break;
   case  C15CLK_SAT: type = "C15CLK_SAT_" + prn;   break;
   case  J15CLK_SAT: type = "J15CLK_SAT_" + prn;   break;
   case  AMB_IF   :  type = "AMB_IF_" + prn ;   break;
   case  AMB13_IF:   type = "AMB13_IF_" + prn;   break;
   case  AMB_WL   :  type = "AMB_WL_" + prn;    break;
   case  AMB_L1   :  type = "AMB_L1_" + prn;    break;
   case  AMB_L2   :  type = "AMB_L2_" + prn;    break;
   case  AMB_L3   :  type = "AMB_L3_" + prn;    break;
   case  AMB_L4   :  type = "AMB_L4_" + prn;    break;  
   case  AMB_L5   :  type = "AMB_L5_" + prn;    break;
   case  FCBS_IF  :  type = "FCBS_IF";   break;
   case  FCBS_L1  :  type = "FCBS_L1";   break;
   case  FCBS_L2  :  type = "FCBS_L2";   break;
   case  FCBS_L3  :  type = "FCBS_L3";   break;
   case  FCBS_L4  :  type = "FCBS_L4";   break;
   case  FCBS_L5  :  type = "FCBS_L5";   break;
   case  FCBR_IF  :  type = "FCBR_IF";   break;
   case  FCBR_L1  :  type = "FCBR_L1";   break;
   case  FCBR_L2  :  type = "FCBR_L2";   break;
   case  FCBR_L3  :  type = "FCBR_L3";   break;
   case  FCBR_L4  :  type = "FCBR_L4";   break;
   case  FCBR_L5  :  type = "FCBR_L5";   break;
   case  CLK_ICB  :  type = "CLK_ICB";   break;
   case  CLUSTERB :  type = "CLUSTERB";  break;
   case  GRD_N    :  type = "GRD_N";     break;
   case  GRD_E    :  type = "GRD_E";     break;
   case  GLO_ISB  :  type = "GLO_ISB";   break;
   case  GLO_IFCB :  type = "GLO_IFCB";  break;
   case  GLO_IFPB :  type = "GLO_IFPB";  break;
   case  GLO_IFB  :  type = "GLO_IFB_"+prn; break;
   case  GPS_IFB_C3: type = "GPS_IFB_C3"; break;
   case  GAL_IFB_C3: type = "GAL_IFB_C3"; break;
   case  BDS_IFB_C3: type = "BDS_IFB_C3"; break;
   case  QZS_IFB_C3: type = "QZS_IFB_C3"; break;
   case  GAL_ISB  :  type = "GAL_ISB";   break;
   case  BDS_ISB  :  type = "BDS_ISB";   break;
   case  QZS_ISB  :  type = "QZS_ISB";   break;
   case  IFB_GPS  :  type = "IFB_GPS";   break;
   case  IFB_BDS:  type = "IFB_BDS";   break;
   case  IFB_GAL:  type = "IFB_GAL";   break;
   case  IFB_QZS:  type = "IFB_QZS";   break;
   case	 IFB_GAL_2: type = "IFB_GAL_2"; break;
   case  IFB_GAL_3: type = "IFB_GAL_3"; break;
   case	 IFB_BDS_2: type = "IFB_BDS_2"; break;
   case  IFB_BDS_3: type = "IFB_BDS_3"; break;
   case  P1P2G_REC :  type = "P1P2G_REC";  break;
   case  P1P2E_REC :  type = "P1P2E_REC";  break;  
   case ATT_X:  type = "ATT_X"; break;
   case ATT_Y:  type = "ATT_Y"; break;
   case ATT_Z:  type = "ATT_Z"; break;
   case eb_X:  type = "eb_X"; break;
   case eb_Y:  type = "eb_Y"; break;
   case eb_Z:  type = "eb_Z"; break;
   case db_X:  type = "db_X"; break;
   case db_Y:  type = "db_Y"; break;
   case db_Z:  type = "db_Z"; break;
   case EXTRINSIC_ATT_X: type = "EXTRINSIC_ATT_X"; break;
   case EXTRINSIC_ATT_Y: type = "EXTRINSIC_ATT_Y"; break;
   case EXTRINSIC_ATT_Z: type = "EXTRINSIC_ATT_Z"; break;
   case EXTRINSIC_CRD_X: type = "EXTRINSIC_CRD_X"; break;
   case EXTRINSIC_CRD_Y: type = "EXTRINSIC_CRD_Y"; break;
   case EXTRINSIC_CRD_Z: type = "EXTRINSIC_CRD_Z"; break;
   case CAM_ATT_X: type = "CAM_ATT_X"; break;
   case CAM_ATT_Y: type = "CAM_ATT_Y"; break;
   case CAM_ATT_Z: type = "CAM_ATT_Z"; break;
   case CAM_CRD_X: type = "CAM_CRD_X"; break;
   case CAM_CRD_Y: type = "CAM_CRD_Y"; break;
   case CAM_CRD_Z: type = "CAM_CRD_Z"; break;


   case  SCALE: type = "SCALE"; break;
   case  GEOCX: type = "GEOCX"; break;
   case  GEOCY: type = "GEOCY"; break;
   case  GEOCZ: type = "GEOCZ"; break;
   case  XPOLE: type = "XPOLE"; break;
   case  YPOLE: type = "YPOLE"; break;
   case  DXPOLE: type = "DXPOLE"; break;
   case  DYPOLE: type = "DYPOLE"; break;
   case  UT1: type = "UT1"; break;
   case  DUT1: type = "DUT1"; break;
   case  PXSAT: type = "PXSAT_"+prn; break;
   case  PYSAT: type = "PYSAT_"+prn; break;
   case  PZSAT: type = "PZSAT_"+prn; break;
   case  VXSAT: type = "VXSAT_"+prn; break;
   case  VYSAT: type = "VYSAT_"+prn; break;
   case  VZSAT: type = "VZSAT_"+prn; break;
   case  ECOM_D0: type = "D0_"+prn; break;
   case  ECOM_Dc: type = "Dc_"+prn; break;
   case  ECOM_Ds: type = "Ds_"+prn; break;
   case  ECOM_Y0: type = "Y0_"+prn; break;
   case  ECOM_Yc: type = "Yc_"+prn; break;
   case  ECOM_Ys: type = "Ys_"+prn; break;
   case  ECOM_B0: type = "X0_"+prn; break;
   case  ECOM_Bc: type = "Xc_"+prn; break;
   case  ECOM_Bs: type = "Xs_"+prn; break;
   case  SR_scale: type = "SR_scale_" + prn; break;
   case  DRAG_c:   type = "DRAG_c_" + prn; break;
   case  EMP_Sc1: type = "EMP_Sc1_" + prn; break;
   case  EMP_Cc1: type = "EMP_Cc1_" + prn; break;
   case  EMP_Sa1: type = "EMP_Sa1_" + prn; break;
   case  EMP_Ca1: type = "EMP_Ca1_" + prn; break;
   case  EMP_Sr1: type = "EMP_Sr1_" + prn; break;
   case  EMP_Cr1: type = "EMP_Cr1_" + prn; break;

   case  EMP_Sc2: type = "EMP_Sc2_" + prn; break;
   case  EMP_Cc2: type = "EMP_Cc2_" + prn; break;
   case  EMP_Sa2: type = "EMP_Sa2_" + prn; break;
   case  EMP_Ca2: type = "EMP_Ca2_" + prn; break;
   case  EMP_Sr2: type = "EMP_Sr2_" + prn; break;
   case  EMP_Cr2: type = "EMP_Cr2_" + prn; break;

   case  EMP_Ac: type = "EMP_Ac_" + prn; break;
   case  EMP_Bc: type = "EMP_Bc_" + prn; break;
   case  EMP_Aa: type = "EMP_Aa_" + prn; break;
   case  EMP_Ba: type = "EMP_Ba_" + prn; break;
   case  EMP_Ar: type = "EMP_Ar_" + prn; break;
   case  EMP_Br: type = "EMP_Br_" + prn; break;
   case  ECOM2_D0: type  =  "ECOM2_D0_" + prn; break;
   case  ECOM2_D2c: type = "ECOM2_D2c_" + prn; break;
   case  ECOM2_D2s: type = "ECOM2_D2s_" + prn; break;
   case  ECOM2_D4c: type = "ECOM2_D4c_" + prn; break;
   case  ECOM2_D4s: type = "ECOM2_D4s_" + prn; break;
   case  ECOM2_Y0: type = "ECOM2_Y0_" + prn; break;
   case  ECOM2_X0: type = "ECOM2_X0_" + prn; break;
   case  ECOM2_X1c: type = "ECOM2_X1c_" + prn; break;
   case  ECOM2_X1s: type = "ECOM2_X1s_" + prn; break;
   case RB:  type = "RB"; break;
   default        :  type = "UNDEF";
  }
   
return type;
}

t_gparhead t_gpar::get_head() const
{
	return t_gparhead(this->parType,this->site,this->prn);
}

t_gtimearc t_gpar::get_timearc() const
{
	return t_gtimearc(this->beg,this->end);
}

// get ZTD mf according to settings   
void t_gpar::_getmf(t_gsatdata& satData, const t_gtriple& crd, const t_gtime& epoch, double& mfw, double& mfh, double& dmfw, double& dmfh)
{
   if(parType != TRP && parType != GRD_N && parType != GRD_E) return;

   double ele = satData.ele();
     
   if (_mf_ztd == COSZ) {
      mfw = mfh = 1.0/sin(ele);
   }else if (_mf_ztd == GMF)  {
      t_gmf mf;
      mf.gmf( epoch.mjd(), crd[0], crd[1], crd[2], G_PI/2.0-ele,
         mfh, mfw, dmfh, dmfw );
   } else cerr << "ZTD mapping function is not set up correctly!!!" << endl;
}   


string gpar2str(const t_gpar& par)
{
	string partype = par.str_type();
	partype = par.site + t_gpar::PAR_STR_SEP + partype;
	if (par.prn == "") {
		partype += t_gpar::PAR_STR_SEP;
	}
	return partype;
}

string BASE_LIBRARY_EXPORT ptype2str(const t_gpar::t_type& parType)
{
	// get data type
	// ----------
		string type;
		switch (parType) 
		{
		case  t_gpar::t_type::CRD_X:  type = "CRD_X";     break;
		case  t_gpar::t_type::CRD_Y:  type = "CRD_Y";     break;
		case  t_gpar::t_type::CRD_Z:  type = "CRD_Z";     break;
		case  t_gpar::t_type::TRP:  type = "TRP";     break;
		case  t_gpar::t_type::SION:  type = "SION_";    break;
		case  t_gpar::t_type::VION:  type = "VION_";    break;
		case  t_gpar::t_type::CLK:  type = "CLK";              break;
		case  t_gpar::t_type::G15CLK:  type = "G15CLK";           break;
		case  t_gpar::t_type::E15CLK:  type = "E15CLK";           break;
		case  t_gpar::t_type::C15CLK:  type = "C15CLK";           break;
		case  t_gpar::t_type::J15CLK:  type = "J15CLK";           break;
		case  t_gpar::t_type::CLK_G:  type = "CLK_G";              break;
		case  t_gpar::t_type::CLK_E:  type = "CLK_E";              break;
		case  t_gpar::t_type::CLK_C:  type = "CLK_C";              break;
		case  t_gpar::t_type::CLK_R:  type = "CLK_R";              break;
		case  t_gpar::t_type::CLK_J:  type = "CLK_J";              break;
		case  t_gpar::t_type::IFCB_F3:  type = "IFCB_F3_";     break;
		case  t_gpar::t_type::IFCB_F4:  type = "IFCB_F4_";     break;
		case  t_gpar::t_type::IFCB_F5:  type = "IFCB_F5_";     break;
		case  t_gpar::t_type::IFB_C3:  type = "IFB_F3_";     break;
		case  t_gpar::t_type::IFB_C4:  type = "IFB_F4_";     break;
		case  t_gpar::t_type::IFB_C5:  type = "IFB_F5_";     break;
		case  t_gpar::t_type::SAT_IFB_C3:  type = "IFB_SAT_F3";    break;
		case  t_gpar::t_type::SAT_IFB_C4:  type = "IFB_SAT_F4";    break;
		case  t_gpar::t_type::SAT_IFB_C5:  type = "IFB_SAT_F5";    break;
		case  t_gpar::t_type::GPS_REC_IFB_C3:  type = "IFB_REC_F3";     break;
		case  t_gpar::t_type::GPS_REC_IFB_C4:  type = "IFB_REC_F4";     break;
		case  t_gpar::t_type::GPS_REC_IFB_C5:  type = "IFB_REC_F5";     break;
		case  t_gpar::t_type::GAL_REC_IFB_C3:  type = "GAL_REC_IFB_C3";     break;
		case  t_gpar::t_type::GAL_REC_IFB_C4:  type = "GAL_REC_IFB_C4";     break;
		case  t_gpar::t_type::GAL_REC_IFB_C5:  type = "GAL_REC_IFB_C5";     break;
		case  t_gpar::t_type::BDS_REC_IFB_C3:  type = "BDS_REC_IFB_C3";     break;
		case  t_gpar::t_type::BDS_REC_IFB_C4:  type = "BDS_REC_IFB_C4";     break;
		case  t_gpar::t_type::BDS_REC_IFB_C5:  type = "BDS_REC_IFB_C5";     break;
		case  t_gpar::t_type::QZS_REC_IFB_C3:  type = "QZS_REC_IFB_C3";     break;
		case  t_gpar::t_type::QZS_REC_IFB_C4:  type = "QZS_REC_IFB_C4";     break;
		case  t_gpar::t_type::QZS_REC_IFB_C5:  type = "QZS_REC_IFB_C5";     break;
		case  t_gpar::t_type::CLK_SAT:  type = "CLK_SAT_";   break;
		case  t_gpar::t_type::G15CLK_SAT:  type = "G15CLK_SAT_";   break;//add by xiongyun
		case  t_gpar::t_type::E15CLK_SAT:  type = "E15CLK_SAT_";   break;
		case  t_gpar::t_type::C15CLK_SAT:  type = "C15CLK_SAT_";   break;
		case  t_gpar::t_type::J15CLK_SAT:  type = "J15CLK_SAT_";   break;
		case  t_gpar::t_type::AMB_IF:  type = "AMB_IF_";   break;
		case  t_gpar::t_type::AMB13_IF:   type = "AMB13_IF_";   break;
		case  t_gpar::t_type::AMB_WL:  type = "AMB_WL_";    break;
		case  t_gpar::t_type::AMB_L1:  type = "AMB_L1_";    break;
		case  t_gpar::t_type::AMB_L2:  type = "AMB_L2_";    break;
		case  t_gpar::t_type::AMB_L3:  type = "AMB_L3_";    break;
		case  t_gpar::t_type::AMB_L4:  type = "AMB_L4_";    break;
		case  t_gpar::t_type::AMB_L5:  type = "AMB_L5_";    break;
		case  t_gpar::t_type::FCBS_IF:  type = "FCBS_IF";   break;
		case  t_gpar::t_type::FCBS_L1:  type = "FCBS_L1";   break;
		case  t_gpar::t_type::FCBS_L2:  type = "FCBS_L2";   break;
		case  t_gpar::t_type::FCBS_L3:  type = "FCBS_L3";   break;
		case  t_gpar::t_type::FCBS_L4:  type = "FCBS_L4";   break;
		case  t_gpar::t_type::FCBS_L5:  type = "FCBS_L5";   break;
		case  t_gpar::t_type::FCBR_IF:  type = "FCBR_IF";   break;
		case  t_gpar::t_type::FCBR_L1:  type = "FCBR_L1";   break;
		case  t_gpar::t_type::FCBR_L2:  type = "FCBR_L2";   break;
		case  t_gpar::t_type::FCBR_L3:  type = "FCBR_L3";   break;
		case  t_gpar::t_type::FCBR_L4:  type = "FCBR_L4";   break;
		case  t_gpar::t_type::FCBR_L5:  type = "FCBR_L5";   break;
		case  t_gpar::t_type::CLK_ICB:  type = "CLK_ICB";   break;
		case  t_gpar::t_type::CLUSTERB:  type = "CLUSTERB";  break;
		case  t_gpar::t_type::GRD_N:  type = "GRD_N";     break;
		case  t_gpar::t_type::GRD_E:  type = "GRD_E";     break;
		case  t_gpar::t_type::GLO_ISB:  type = "GLO_ISB";   break;
		case  t_gpar::t_type::GLO_IFCB:  type = "GLO_IFCB";  break;
		case  t_gpar::t_type::GLO_IFPB:  type = "GLO_IFPB";  break;
		case  t_gpar::t_type::GLO_IFB:  type = "GLO_IFB_"; break;
		case  t_gpar::t_type::GAL_ISB:  type = "GAL_ISB";   break;
		case  t_gpar::t_type::BDS_ISB:  type = "BDS_ISB";   break;
		case  t_gpar::t_type::QZS_ISB:  type = "QZS_ISB";   break;
		case  t_gpar::t_type::IFB_GPS:  type = "IFB_GPS";   break;
		case  t_gpar::t_type::IFB_BDS:  type = "IFB_BDS";   break;
		case  t_gpar::t_type::IFB_GAL:  type = "IFB_GAL";   break;
		case  t_gpar::t_type::IFB_GAL_2:  type = "IFB_GAL_2";   break;
		case  t_gpar::t_type::IFB_GAL_3:  type = "IFB_GAL_3";   break;
		case  t_gpar::t_type::IFB_BDS_2:  type = "IFB_BDS_2";   break;
		case  t_gpar::t_type::IFB_BDS_3:  type = "IFB_BDS_3";   break;
		case  t_gpar::t_type::IFB_QZS:  type = "IFB_QZS";   break;
		case  t_gpar::t_type::P1P2G_REC:  type = "P1P2G_REC";  break;
		case  t_gpar::t_type::P1P2E_REC:  type = "P1P2E_REC";  break;

		case  t_gpar::t_type::SCALE: type = "SCALE"; break;
		case  t_gpar::t_type::GEOCX: type = "GEOCX"; break;
		case  t_gpar::t_type::GEOCY: type = "GEOCY"; break;
		case  t_gpar::t_type::GEOCZ: type = "GEOCZ"; break;
		case  t_gpar::t_type::XPOLE: type = "XPOLE"; break;
		case  t_gpar::t_type::YPOLE: type = "YPOLE"; break;
		case  t_gpar::t_type::DXPOLE: type = "DXPOLE"; break;
		case  t_gpar::t_type::DYPOLE: type = "DYPOLE"; break;
		case  t_gpar::t_type::UT1: type = "UT1"; break;
		case  t_gpar::t_type::DUT1: type = "DUT1"; break;
		case  t_gpar::t_type::PXSAT: type = "PXSAT"; break;
		case  t_gpar::t_type::PYSAT: type = "PYSAT"; break;
		case  t_gpar::t_type::PZSAT: type = "PZSAT"; break;
		case  t_gpar::t_type::VXSAT: type = "VXSAT"; break;
		case  t_gpar::t_type::VYSAT: type = "VYSAT"; break;
		case  t_gpar::t_type::VZSAT: type = "VZSAT"; break;
		case  t_gpar::t_type::ECOM_D0: type = "D0"; break;
		case  t_gpar::t_type::ECOM_Dc: type = "Dc"; break;
		case  t_gpar::t_type::ECOM_Ds: type = "Ds"; break;
		case  t_gpar::t_type::ECOM_Y0: type = "Y0"; break;
		case  t_gpar::t_type::ECOM_Yc: type = "Yc"; break;
		case  t_gpar::t_type::ECOM_Ys: type = "Ys"; break;
		case  t_gpar::t_type::ECOM_B0: type = "X0"; break;
		case  t_gpar::t_type::ECOM_Bc: type = "Xc"; break;
		case  t_gpar::t_type::ECOM_Bs: type = "Xs"; break;
		case  t_gpar::t_type::SR_scale: type = "SR_scale"; break;
		case  t_gpar::t_type::DRAG_c:   type = "DRAG_c"  ; break;
		case  t_gpar::t_type::EMP_Sc1: type = "EMP_Sc1"; break;
		case  t_gpar::t_type::EMP_Cc1: type = "EMP_Cc1"; break;
		case  t_gpar::t_type::EMP_Sa1: type = "EMP_Sa1"; break;
		case  t_gpar::t_type::EMP_Ca1: type = "EMP_Ca1"; break;
		case  t_gpar::t_type::EMP_Sr1: type = "EMP_Sr1"; break;
		case  t_gpar::t_type::EMP_Cr1: type = "EMP_Cr1"; break;

		case  t_gpar::t_type::EMP_Sc2: type = "EMP_Sc2"; break;
		case  t_gpar::t_type::EMP_Cc2: type = "EMP_Cc2"; break;
		case  t_gpar::t_type::EMP_Sa2: type = "EMP_Sa2"; break;
		case  t_gpar::t_type::EMP_Ca2: type = "EMP_Ca2"; break;
		case  t_gpar::t_type::EMP_Sr2: type = "EMP_Sr2"; break;
		case  t_gpar::t_type::EMP_Cr2: type = "EMP_Cr2"; break;

		case  t_gpar::t_type::EMP_Ac: type = "EMP_Ac"; break;
		case  t_gpar::t_type::EMP_Bc: type = "EMP_Bc"; break;
		case  t_gpar::t_type::EMP_Aa: type = "EMP_Aa"; break;
		case  t_gpar::t_type::EMP_Ba: type = "EMP_Ba"; break;
		case  t_gpar::t_type::EMP_Ar: type = "EMP_Ar"; break;
		case  t_gpar::t_type::EMP_Br: type = "EMP_Br"; break;
		case  t_gpar::t_type::ECOM2_D0: type = "ECOM2_D0"   ; break;
		case  t_gpar::t_type::ECOM2_D2c: type = "ECOM2_D2c" ; break;
		case  t_gpar::t_type::ECOM2_D2s: type = "ECOM2_D2s" ; break;
		case  t_gpar::t_type::ECOM2_D4c: type = "ECOM2_D4c" ; break;
		case  t_gpar::t_type::ECOM2_D4s: type = "ECOM2_D4s" ; break;
		case  t_gpar::t_type::ECOM2_Y0: type = "ECOM2_Y0"; break;
		case  t_gpar::t_type::ECOM2_X0: type = "ECOM2_X0"; break;
		case  t_gpar::t_type::ECOM2_X1c: type = "ECOM2_X1c"; break;
		case  t_gpar::t_type::ECOM2_X1s: type = "ECOM2_X1s"; break;
		case  t_gpar::t_type::RB: type = "RB"; break;
		default:  type = "UNDEF";
		}
		return type;
}


t_gpar str2gpar(const string& str_par)
{
	int sep_first,sep_last = 0;
	string site, sat, str_type;
	sep_first = str_par.find(t_gpar::PAR_STR_SEP);
	site = str_par.substr(0, sep_first);

	sep_last = str_par.rfind(t_gpar::PAR_STR_SEP);
	sat = str_par.substr(sep_last + 1);

	str_type = str_par.substr(sep_first + 1, sep_last - sep_first - 1);

	t_gpar::t_type par_type(t_gpar::NO_DEF);

	for (int i=0;i <t_gpar::PAR_STR.size();i++)
	{
		if (t_gpar::PAR_STR[i] == str_type) {
			par_type = t_gpar::t_type(i);
			break;
		}
	}
	
	return t_gpar(site, par_type, 0, sat);
}

   
} // namespace
