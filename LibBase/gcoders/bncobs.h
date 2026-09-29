
#ifndef BNCOBS_H
#define BNCOBS_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
	  Research Institute of Geodesy, Topography and Cartography
	  Ondrejov 244, 251 65, Czech Republic

  Purpose:
  Version: $Rev:$

  2011-04-20 /JD: created
  2012-09-27 /JD: L-obs in whole cycles instead of meters!
  2013-11-19 /JD: shared pointers
  2019-10-25 /zhShen: avaliable

-*/

#include "gcoders/gcoder.h"
#include "gdata/gobsgnss.h"
#include "gall/gallobs.h"

using namespace std;

namespace gnut {

	//const string begEpoch = "BEGEPOCH";
	//const string endEpoch = "ENDEPOCH";

	class BASE_LIBRARY_EXPORT t_bncobs : public t_gcoder {

	public:
		t_bncobs(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);
		~t_bncobs() {};

		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg);
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

		virtual  int encode_head(char* buff, int sz, vector<string>& errmsg);
		virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

		virtual bool available(const t_gtime& now);
		virtual void available_false();
		virtual void available_true();

	protected:

	private:
		void    _check_recobj(string& site);

		bool    _validepo(const t_gtime& t);
		bool    _begepoch;
		// bool    _endepoch;
		t_gtime _tt;

		bool _available;


	};

} // namespace

#endif
