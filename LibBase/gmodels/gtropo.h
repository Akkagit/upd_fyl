
#ifndef GTROPO_H
#define GTROPO_H 

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: implements troposphere model class
  Version: $ Rev: $

  2011-01-10 /JD: created

-*/

#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gprod/gprodcrd.h"
#include "gmodels/ggpt.h"
#include "gutils/ExportBase.h"

using namespace std;

namespace gnut
{

	class BASE_LIBRARY_EXPORT t_gtropo
	{
	public:
		//t_gtropo(string site);
		t_gtropo();
		virtual ~t_gtropo();

		virtual double getZHD(const t_gtriple& ell, const t_gtime& epo);  // ! Radians: Ell[0] and Ell[1]
		virtual double getZWD(const t_gtriple& ell, const t_gtime& epo);  // ! Radians: Ell[0] and Ell[1]

	protected:
		t_gpt        _gpt;
	};

	// ----------------------------- SITE MODELS ------------------------------------

	class BASE_LIBRARY_EXPORT t_saast : public t_gtropo {
	public:
		t_saast() {}
		~t_saast() {}

		virtual double getSTD(const double& ele, const double& hel);   // ! Radians: elevation
		virtual double getZHD(const t_gtriple& ell, const t_gtime& epo);  // ! Radians: Ell[0] and Ell[1]
		virtual double getZWD(const t_gtriple& ell, const t_gtime& epo);  // ! Radians: Ell[0] and Ell[1]
	};


	class BASE_LIBRARY_EXPORT t_davis : public t_gtropo {
	public:
		t_davis() {}
		~t_davis() {}

		virtual double getZHD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
		virtual double getZWD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]  
	};


	class BASE_LIBRARY_EXPORT t_hopf : public t_gtropo {
	public:
		t_hopf() {}
		~t_hopf() {}

		virtual double getZHD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
		virtual double getZWD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
	};


	class BASE_LIBRARY_EXPORT t_baby : public t_gtropo {
	public:
		t_baby() {}
		~t_baby() {}

		virtual double getZHD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
		virtual double getZWD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
	};


	class BASE_LIBRARY_EXPORT t_chao : public t_gtropo {
	public:
		t_chao() {}
		~t_chao() {}

		virtual double getZHD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
		virtual double getZWD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
	};


	class BASE_LIBRARY_EXPORT t_ifad : public t_gtropo {
	public:
		t_ifad() {}
		~t_ifad() {}

		virtual double getZHD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
		virtual double getZWD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
	};

} // namespace

#endif
