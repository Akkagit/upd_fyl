
#ifndef GSPPLSQ_H
#define GSPPLSQ_H

/* ----------------------------------------------------------------------
-*/

#include "gproc/gspp.h"
#include "gmodels/gdop.h"
#include "gproc/glsq.h"
#include "gall/gallpar.h"
#include "gall/gallbias.h"
#include "gall/gallprod.h"
#include "gutils/gsysconv.h"
#include "gmodels/gstochasticmodel.h"
#include "gmodels/gsppmodel.h"
#include "gset/gsetlsq.h"
#include "gutils/ExportLSQ.h"

#define SPP_MINSAT  (6)    // minimum number of satellites
namespace great
{
	class LSQ_LIBRARY_EXPORT t_gspplsq : public virtual t_gspp
	{
	public:
		// Constructor
		t_gspplsq(string mark, t_gsetbase* set, t_gallobj* obj, bool dynamic = false);

		// Destructor
		virtual ~t_gspplsq();

		// Main function of processing
		virtual int processBatch(const t_gtime& beg, const t_gtime& end);

		int ProcessOneEpoch(const t_gtime& now, vector<t_gsatdata>* vec_obs = NULL);


		double getRecClk(const t_gtime& time);
		t_gtriple getvel(const t_gtime& time);
		t_gtriple getCrd(const t_gtime& time);
		t_gtriple getAvergeCrd(); // for static site
		set<string> outcomsat(map<string, t_gsatdata>& satdata);
		map<t_gtime, int> getobsnum();

	protected:
		// result 
		map<t_gtime, t_gtriple> _map_crd;
		// Set Estimation objects   
		t_glsq*        _lsq = nullptr;
		// get settings
		int			  _get_own_settings();

		virtual void	  _setOut() override;
		// Prepare for data
		virtual int    _prepareData();

		// Processing functions for each epoch
		virtual int    _processEpoch(const t_gtime& runEpoch);
		// Adding pseudorange observations to construct B, OMC and P Matrix
		virtual int    _addObsP(t_gsatdata& satdata, unsigned int& iobs, t_gtriple& ell, Matrix& B, ColumnVector& OMC, DiagonalMatrix& P);
		// Adding phase observations to 
		virtual int	  _addObsL(t_gsatdata& satdata, unsigned int& iobs, t_gtriple& ell, Matrix& B, ColumnVector& OMC, DiagonalMatrix& P);
		/**
		 * @brief add doppler measurement to model
		 * @param[in] satdata
		 * @param[in] iobs
		 * @param[in] XYZ
		 * @param[out] A
		 * @param[out] l
		 * @param[out] P
		 * @return int
		 */
		int _addObsD(t_gsatdata& satdata, unsigned int& iobs, t_gtriple& XYZ, Matrix& A, ColumnVector& l, DiagonalMatrix& P);

		// Satelite position
		virtual int    _satPos(t_gtime& epo, t_gsatdata& gsatdata);

		virtual int	  _apply_tides(t_gtime& _epoch, t_gtriple& xRec);

		double		  _weightObs(t_gsatdata& satdata, OBSCOMBIN obstype, t_gobs* g1, t_gobs* g2 = 0);

		virtual void   _printout(const t_glsq& lsqestimator);
		// All the sat data
		vector<t_gsatdata> _data;

		// The minimum satellite number
		unsigned int       _minsat;
		// Sigma in unit weight
		double             _sig_unit;

		// Epoch time
		t_gtime           _epoch;
		// Begin time (in set)
		t_gtime			 _beg;
		// End  time (in set)
		t_gtime			 _end;
		// used gnss system
		set<string>	     _sys;
		// For caculate the dop value, including pdop,vdop,hdop and so on
		t_gdop            _dop;
		// output file
		t_giof*			 _result = nullptr;
		map<t_gtime, int> _obs_num;
		// estimate velocity
		bool                  _dynamics = false;
		bool                 _auto_band = false;
		CBIASCHAR            _cbiaschar;

		double               _sig_ion;
	};

} // namespace

#endif // GSPPFLT_H
