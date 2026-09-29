
#ifndef GREC_LEO_H
#define GREC_LEO_H

/* ----------------------------------------------------------------------

  2011-01-10 /HJZheng: created

-*/



#include "gdata/grec.h"
#include "gutils/ExportBase.h"


using namespace std;

namespace gnut {

	class BASE_LIBRARY_EXPORT t_grec_leo : public  t_grec{

	public:
		t_grec_leo();
		//   t_grec(const t_grec& obj);
		virtual ~t_grec_leo();


		//virtual void sync_orb(t_gallnav* orb) override;

	protected:

	private:
	};

} // namespace

#endif
