#ifndef AUG_H
#define AUG_H

#include <string> 
#include <vector> 
#include <memory>

#include "gcoders/gcoder.h"
#include "gdata/gaug.h"

using namespace std;
using namespace gnut;

namespace great
{
	class  BASE_LIBRARY_EXPORT t_aug : public t_gcoder
	{
	public:

		/** @brief default constructor. */
		t_aug(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_aug();

		/** @brief decode the header of the aug data file. */
		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg)override;

		/** @brief decode the data body of the aug data file. */
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)override;

	protected:

		map<GSYS, vector<t_augtype>> _augtype;
		string _site;
		t_gtime _epoch;
	};
}

#endif
