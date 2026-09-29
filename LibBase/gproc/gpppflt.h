
#ifndef GPPPFLT_H
#define GPPPFLT_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
	  Research Institute of Geodesy, Topography and Cartography
	  Ondrejov 244, 251 65, Czech Republic

  Purpose: implements ppp client
  Version: $ Rev: $

  2011-06-06 /PV: created
  2011-06-06 /JD: revised for multi-threading

-*/

#include <fstream>

#include "gproc/gppp.h"
#include "gproc/gsppflt.h"
#include "gall/gallfltmat.h"
#include "gio/gxml.h"
#include "gio/gtcp.h"

namespace gnut
{

	class BASE_LIBRARY_EXPORT t_gpppflt : public t_gppp,
		public t_gsppflt,
		public t_gxml
	{
	public:
		t_gpppflt();
		t_gpppflt(string mark, t_gsetbase* set);
		virtual ~t_gpppflt();
		void add_coder(const vector<t_gcoder*>& coder);

		virtual int processBatch(const t_gtime& beg, const t_gtime& end, bool prtOut);
		virtual int processBatchFB(const t_gtime& beg, const t_gtime& end);
		void debugPrint();

	protected:
		// Process one epoch
		virtual int    _processEpoch(const t_gtime& runEpoch);

		// add one satellite data to A, l, P - carrier phase
		virtual int    _addObsL(t_gsatdata& satdata, unsigned int& iobs, t_gtriple& ell, Matrix& A, ColumnVector& l, DiagonalMatrix& P);

		// Predict state vector and covariance matrix
		virtual void   _predict();

		// Satelite position
		virtual int    _satPos(t_gtime&, t_gsatdata&);

		//   virtual int    _prepareData();

		   // Applying tides
		virtual int    _apply_tides(t_gtime& _epoch, t_gtriple& xRec);

		// Synchronize ambiguity with respect to observation data.
		void           _syncAmb();

		// Set output & return if smooth
		virtual void   _setOut();

		// Set output & return if smooth
		void           _prtOut(t_gtime& epo, t_gallpar& X, const SymmetricMatrix& Q, vector<t_gsatdata>& data, ostringstream& os, xml_node& node, bool saveProd = true);

		// Save apriory coordinates to products
		void           _saveApr(t_gtime& epo, t_gallpar& X, const SymmetricMatrix& Q);

		// Remove apriory coordinates from products
		void           _removeApr(t_gtime& epo);


		xml_node       line;
		bool          _read;
		t_giof*       _flt;
		t_giof*       _enufile;
		map<OFMT, t_gtcp*>  _maptcp;  // zhshen
		t_giof*       _smt;
		string        _kml_name;
		bool          _kml;
		bool          _restarted;
		bool          _beg_end;

		// Models   
		t_randomwalk* _grdStoModel;
		t_randomwalk* _ambStoModel;

		// Update time for RDW processes
		virtual void  _timeUpdate(const t_gtime& epo);

		// Data for smoothing (Rauch-Tung-Striebel algorithm)   
		t_gfltmat     _fltdata;
		t_gallfltmat  _vfltdata;
		vector<t_gcoder*>         _gcoder;  // zhShen

		virtual void  _backflt(bool prtOut);

		// Intervals for reseting parameters
		int           _reset_amb;
		void          _reset_ambig();
		void          _check_ambig();
		int           _reset_par;
		void          _reset_param();
		t_gtriple     xyz_standard;	    ///< enu coord reference
		map<string, int> _glofrq_num;   ///< glofrq number, set GLO frequency number
		pair<int, t_gtriple>   neu_sum_nepo;


		//   ofstream      _iono_ofstr;
		//   ofstream      _ifcb_ofstr;

		// satellite pco correction
		//   void          _satPCO_corr(t_gtime& epo, t_gsatdata& gsatdata);   
	};

} // namespace

#endif // GPPPFLT_H
