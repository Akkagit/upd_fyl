

#ifndef GSETORBFITLEO_H
#define GSETORBFITLEO_H

#define XMLKEY_ORBDIF		  "orbdifleo"
#define XMLKEY_ORBDIF_LEO     "leo"
#define XMLKEY_ORBDIF_FIT     "fit"
#define XMLKEY_ORBDIF_TRANS   "trans"
#define XMLKEY_ORBDIF_UNIT    "unit"
#define XMLKEY_ORBDIF_INTER   "inter"

#include "gset/gsetbase.h"
#include "gutils/ExportOI.h"
using namespace gnut;

namespace great
{
	class OI_LIBRARY_EXPORT t_gsetorbfitleo : public virtual t_gsetbase
	{
	public:
		/** @brief constructor. */
		t_gsetorbfitleo();

		/** @brief destructor. */
		virtual ~t_gsetorbfitleo();

		/**
		* @brief get the trans parameteres of orbit differ.
		* @return the trans parameteres(STRD)
		*/
		string trans();

		

		/**
		* @brief get the leos that need to differ.
		* @return the exclude sat
		*/
		set<string> leo();

		double unit();
		bool isfit();
		bool isinter();

		/** @brief settings check. */
		void check();

		/** @brief settings help. */
		void help();


	private:
		
	};
}
#endif
