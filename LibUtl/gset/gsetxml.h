/**
*

* @verbatim
	<!-->================================================================
	app_type    : PPP_IF,PPP_ALL,PCE_IF,PCE_ALL,POD_IF,POD_ALL
	if exist    : +YYYY+\+DOY+\ in the path will be replaced as year and doy
	rinexc_root : root dir, example rinexc_root/2019/001
	rinexo_root : root dir, example rinexo_root/2019/001
	project_root: root dir, example project_root/2019/001
	sysfile_root: root dir, example sysfile_root/oceanload
	================================================================<!-->
	<xml>
		<app_type>      </app_type>  <!--> <!-->
		<sp3_type>      </sp3_type>

		<time_interval>   </time_interval>
		<scan_rinexo>   </scan_rinexo>
		<sacn_sysfile>  </sacn_sysfile>
		<rinexc_type>   </rinexc_type>
		<sp3_root>		</sp3_root>
		<brd_root>		</brd_root>
		<dcb_root>      </dcb_root>
		<sinex_root>    </sinex_root>
		<rinexc_root>   </rinexc_root>
		<rinexo_root>   </rinexo_root>
		<sysfile_root>  </sysfile_root>
		<project_root>  </project_root>

	</xml>
 @endverbatim

* @verbatim
	History
	 -1.0 jdhuang	 2019-07-17 creat the xml file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gsetpar.h
* @brief		the set for creat XML files
*
* @author       jdhuang, Wuhan University
* @version		1.0.0
* @date			2019-04-07
*
*/

#ifndef GSETXML_H
#define GSETXML_H

#define XMLKEY_XML        "xml"

#include "gset/gsetbase.h"
#include "gutils/ExportUtl.h"

using namespace gnut;

namespace great
{
	class UTL_LIBRARY_EXPORT t_gsetxml : public virtual t_gsetbase
	{
	public:
		t_gsetxml();
		virtual ~t_gsetxml();

		void help();
		void check() {};

	public:
		string xml_eg();
		string app_type();
		string sp3_type();
		string rinexc_type();
		string sp3_root();
		string brd_root();
		string dcb_root();
		string sinex_root();
	    string rinexc_root();
		string rinexo_root();
		string sysfile_root();
		string project_root();

		double time_interval();
		bool scan_sysfile();
		bool scan_rinexo();
	};
}

#endif