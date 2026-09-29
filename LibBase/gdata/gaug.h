#ifndef GAUG_H
#define GAUG_H

#include "gdata/gdata.h"
#include "gutils/gnss.h"
#include "gutils/gtime.h"
#include "gutils/ExportBase.h"
#include "gutils/gtriple.h"

using namespace gnut;

namespace great
{
	enum AUGTYPE {
		TYPE_P = 1,
		TYPE_L = 2,
		TYPE_ION = 3,
		TYPE_TRP = 4,
		TYPE_NON = 999
	};

	typedef pair<AUGTYPE, GOBSBAND> t_augtype; ///< aug type
	typedef map<t_augtype, double>  t_type_aug;///< aug_type--value
	typedef map<string, t_type_aug> t_sat_aug; ///< satellite--aug_type
	typedef map<t_gtime, t_sat_aug> t_epo_aug; ///< epoch--satellite
	typedef map<string, t_epo_aug>  t_rec_aug; ///< site--epoch

	class BASE_LIBRARY_EXPORT t_gaug : public t_gdata
	{
	public:
		t_gaug();
		virtual ~t_gaug();
		t_gtime endTime_aug(string& site);

		void add_data(string site, const t_gtime& epo, string sat, t_augtype type, double value);
		t_rec_aug* getAug();
		t_epo_aug* getAug(string site);
		t_sat_aug* getAug(string site, const t_gtime& epo);
		t_type_aug* getAug(string site, const t_gtime& epo, string sat);
	protected:
		t_rec_aug _gaug;
	};

}
#endif