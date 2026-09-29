
#ifndef BNCCORR_H
#define BNCCORR_H

/* ----------------------------------------------------------------------

* @file bnccorr.h
* @brief Receive the BNC Data of correction from port,about the Orbit, Clock, Code bias, Phase bias, VTEC.
* @verbatim
	History
	-1.0	ShenZhiheng	      2019-11-10 creat the file.
	-1.1	ShenZhiheng	      2019-11-23 add the CORR_TYPE , intv , available .
	-1.2	ShenZhiheng	      2020-04-18 complete the decode_cd, decode_pb function .
  @endverbatim

* @author ShenZhiheng
* @version 1.0
* @date 2019-11-10
* @license
-*/

#include "gcoders/gcoder.h"
#include "gall/gallprec.h"
#include "gall/gallbias.h"

using namespace std;

namespace gnut {


	class BASE_LIBRARY_EXPORT t_bnccorr : public t_gcoder {

	public:

		enum CORR_TYPE
		{
			ORBIT,
			CLOCK,
			VTEC,
			CODE_BIAS,
			PHASE_BIAS,
			UPD_COR,
			AUG_COR,
			DEFAULT
		};


		map<int, int>intv{ {0,1}, {1,2}, {2,5},{3,10},{4,15},{ 5,30 },
						 { 6,60 },{ 7,120 },{ 8,240 },{ 9,300 },{10,600} ,
						 {11,900},{12,1800},{13,3600},{14,7200},{15,10800} };

		t_bnccorr(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);
		~t_bnccorr() {};

		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg);
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);


		virtual bool available(const t_gtime& now);
		virtual void available_false();
		virtual void available_true();

		string corr2str(const CORR_TYPE& type);
		CORR_TYPE str2corr(const string& str);

	protected:


		virtual int _decode_orbit(const t_gtime& epo, istringstream& ss);
		virtual int _decode_clock(const t_gtime& epo, istringstream& ss);
		virtual int _decode_cb(const t_gtime& epo, istringstream& ss);
		virtual int _decode_pb(const t_gtime& epo, istringstream& ss);
		virtual int _decode_vtec(const t_gtime& epo, istringstream& ss);
		virtual int _decode_upd(const t_gtime& epo, istringstream& ss);
		virtual int _decode_aug(const string& site, const t_gtime& epo, istringstream& ss);

	private:
		void    _check_recobj(string& site);

		bool    _validepo(const t_gtime& t);
		bool    _begepoch;
		// bool    _endepoch;
		t_gtime _tt;
		CORR_TYPE _type;
		unsigned int	  _udclkInt;	// SSR clk interval (added by zhShen)
		unsigned int	  _udorbInt;    // SSR orb interval (added by zhShen)
		bool _available;

		map< GSYS, map<FREQ_SEQ, GOBSBAND> > _band_index;


	};

} // namespace

#endif
