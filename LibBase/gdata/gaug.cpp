#include "gdata/gaug.h"

using namespace std;
namespace great
{
	t_gaug::t_gaug()
		:t_gdata()
	{
		gtrace("t_gaug::constructor");
		id_type(t_gdata::AUG);
	}

	t_gaug::~t_gaug()
	{
		gtrace("t_gaug::destructor");
	}
	t_gtime t_gaug::endTime_aug(string& site)
	{
		t_epo_aug* it_epo_aug = this->getAug(site);
		if (!it_epo_aug) return t_gtime();
		auto it_endepo = it_epo_aug->end();
		if (it_endepo != it_epo_aug->begin()) it_endepo--;
		else return t_gtime();

		return it_endepo->first;
	}
	void t_gaug::add_data(string site, const t_gtime& epo, string sat, t_augtype type, double value)
	{
		_gaug[site][epo][sat][type] = value;
	}
	t_rec_aug* t_gaug::getAug()
	{
		return &_gaug;
	}
	t_epo_aug* t_gaug::getAug(string site)
	{
		if (_gaug.find(site) != _gaug.end())return &_gaug[site];
		else return NULL;
	}
	t_sat_aug* t_gaug::getAug(string site, const t_gtime& epo)
	{
		if (_gaug.find(site) != _gaug.end() &&
			_gaug[site].find(epo) != _gaug[site].end())
			return &_gaug[site][epo];
		else return NULL;
	}
	t_type_aug* t_gaug::getAug(string site, const t_gtime& epo, string sat)
	{
		if (_gaug.find(site) != _gaug.end() &&
			_gaug[site].find(epo) != _gaug[site].end()&&
			_gaug[site][epo].find(sat)!= _gaug[site][epo].end())
			return &_gaug[site][epo][sat];
		else return NULL;
	}
}