/**
*
* @verbatim
	History
	 -1.0 Jiande Huang  2019-07-18  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gxmlfactory.cpp
* @brief		main control class for get xml
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-18
*
*/

#include "gproc/gxmlfactory.h"
#include "gset/gsetgen.h"
#include "gset/gsetgnss.h"
#include "gset/gsetxml.h"
#include "gset/gsetpar.h"
#include "gset/gsetrec.h"
#include "gio/gfile.h"

#include "gcoders/sinex.h"
#include "gcoders/gcoder.h"

#include "gutils/gtriple.h"
#include "gutils/gstring.h"


#include "gutils/ginfolog.h"

#ifdef WIN32
#include <io.h>
#include <direct.h> 
#else
#include <unistd.h>
#include <sys/stat.h>
#endif
#include <stdint.h>
#include <string>

#define MAX_PATH_LEN 256

#ifdef WIN32
#define ACCESS(fileName,accessMode) _access(fileName,accessMode)
#define MKDIR(path) _mkdir(path)
#else
#define ACCESS(fileName,accessMode) access(fileName,accessMode)
#define MKDIR(path) mkdir(path,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif
using namespace gnut;
namespace great
{
	t_gxmlfactory::t_gxmlfactory()
	{
	}

	t_gxmlfactory::t_gxmlfactory(t_gsetbase * set, t_glog * log) :
		_set(set),
		_log(log)
	{
		_beg_time = dynamic_cast<t_gsetgen*>(set)->beg();
		_end_time = dynamic_cast<t_gsetgen*>(set)->end();
		_rec_list = dynamic_cast<t_gsetgen*>(set)->recs();
		_sat_list = dynamic_cast<t_gsetgnss*>(set)->sat();
		_sys_list = dynamic_cast<t_gsetgen*>(set)->sys();

		_xml_eg       = dynamic_cast<t_gsetxml*>(set)->xml_eg();
		_app_type     = dynamic_cast<t_gsetxml*>(set)->app_type();
		_sp3_type     = dynamic_cast<t_gsetxml*>(set)->sp3_type();
		_rinexc_type  = dynamic_cast<t_gsetxml*>(set)->rinexc_type();
		_dcb_root     = dynamic_cast<t_gsetxml*>(set)->dcb_root();
		_sp3_root     = dynamic_cast<t_gsetxml*>(set)->sp3_root();
		_brd_root     = dynamic_cast<t_gsetxml*>(set)->brd_root();
		_sinex_root   = dynamic_cast<t_gsetxml*>(set)->sinex_root();
		_rinexc_root  = dynamic_cast<t_gsetxml*>(set)->rinexc_root();
		_rinexo_root  = dynamic_cast<t_gsetxml*>(set)->rinexo_root();
		_sysfile_root = dynamic_cast<t_gsetxml*>(set)->sysfile_root();
		_project_root = dynamic_cast<t_gsetxml*>(set)->project_root();
		_scan_sysfile = dynamic_cast<t_gsetxml*>(set)->scan_sysfile();
		_scan_rinexo  = dynamic_cast<t_gsetxml*>(set)->scan_rinexo();
		_time_interval  = dynamic_cast<t_gsetxml*>(set)->time_interval();
		_inp_list     = dynamic_cast<t_gsetinp*>(set)->inputs_all();

#ifdef WIN32
		_isWIN = true;
#else
		_isLINUX = true;
#endif

	}

	t_gxmlfactory::~t_gxmlfactory()
	{

	}

	bool t_gxmlfactory::ProcessBatch()
	{
		bool valid = false;
		if (_beg_time > _end_time)
		{
			write_log_info(_log, 0, "ERROR", "The beg time > end time !");
			write_log_info(_log, 0, "NOTE ", "The beg time is : " + _beg_time.str_ymdhms());
			write_log_info(_log, 0, "NOTE ", "The end time is : " + _end_time.str_ymdhms());
			return false;
		}

		if (_rinexc_type.empty() ||
			_rinexc_root.empty() ||
			_rinexo_root.empty() ||
			_app_type.empty() ||
			_sp3_type.empty() ||
			_rinexc_type.empty() ||
			_sysfile_root.empty() ||
			_project_root.empty())
		{
			write_log_info(_log, 0, "ERROR", "File path incompleteness !");
			return false;
		}

		valid = _init_proj_dirs(_project_root, _beg_time, _end_time);
		if (!valid)
		{
			write_log_info(_log, 0, "ERROR", "File path creation error, make sure the path is correct and privileged !");
			return false;
		}

		t_gtime crt_time = _beg_time;
		while (crt_time < _end_time)
		{
			string yyyydoy = crt_time.str_yyyydoy();
			_init_xml_file(_xml_list[yyyydoy], crt_time);

			crt_time = crt_time + 86400.0;
		}
		return true;
	}

	bool t_gxmlfactory::_init_proj_dirs(string root, const t_gtime & beg_time, const t_gtime & end_time)
	{
		write_log_info(_log, 0, "NOTE ", "Start checking and preparing file paths !");
		t_gtime crt_time = _beg_time;

		while (crt_time < _end_time)
		{
			write_log_info(_log, 0, "NOTE ", "The current processing date is : " + crt_time.str_ymd());
			string proj_path;
			string xml_path;
			string yyydoy = crt_time.str_yyyydoy();
			string yyyy = crt_time.str_yyyy();
			string doy = crt_time.str_doy();

			string crt_proj_path = _project_root;
			if (crt_proj_path.find("+YYYY+") != string::npos)
			{
				crt_proj_path.replace(crt_proj_path.find_first_of("+"), 6, yyyy);
			}

			if (crt_proj_path.find("+DOY+") != string::npos)
			{
				crt_proj_path.replace(crt_proj_path.find_first_of("+"), 5, doy);
			}

			if (_isWIN)
			{
				proj_path = crt_proj_path + "\\";
			}
			else if (_isLINUX)
			{
				proj_path = crt_proj_path + "/" + yyydoy + "/";
			}
			else
			{
				write_log_info(_log, 0, "ERROR", "Unrecognizable System !");
				return false;
			}
			
			xml_path = proj_path + _app_type + "_" + yyydoy + ".xml";
			createDirectory(proj_path);   // 返回 0 表示创建成功，-1 表示失败
			
			_proj_list[yyydoy] = proj_path;
			_xml_list[yyydoy] = xml_path;

			crt_time = crt_time + 86400.0;
		}

		return true;
	}

	bool t_gxmlfactory::_init_xml_file(string path, const t_gtime & crt_time)
	{
		xml_document xdoc;
		xdoc.append_child(pugi::node_comment).set_value("The xml file is the settings for GREAT software.");
		
		pugi::xml_node pre = xdoc.prepend_child(pugi::node_declaration);
		pre.append_attribute("version").set_value("1.0");
		pre.append_attribute("encoding").set_value("utf-8");
		

		if (xdoc.child("config"))
		{
			xdoc.remove_child("config");
		}
		xdoc.append_child("config");
		pugi::xml_node config = xdoc.child("config");
		config.append_child(pugi::node_comment).set_value("The main node for GREAT xml file.");

		if (_app_type == "POD_IF")
		{
			bool xml_valid = _init_podIF_xml_file(crt_time, config);
			if (!xml_valid) return false;
		}
		else if (_app_type == "PCE_IF")
		{
			bool xml_valid = _init_pceIF_xml_file(crt_time, config);
			if (!xml_valid) return false;
		}

		//xdoc.print(std::cout);
		cerr << " save the xml file as : " + path << endl;
		xdoc.save_file(path.c_str());

		return true;
	}

	bool t_gxmlfactory::_init_podIF_xml_file(const t_gtime& crt_time, pugi::xml_node & parent)
	{
		bool valid = true;
		valid = _insert_gen(crt_time, parent);			if (!valid) return false;
		valid = _insert_gnss(parent);			        if (!valid) return false;
		valid = _insert_inp_POD_IF(crt_time, parent);	if (!valid) return false;
		valid = _insert_par(crt_time, parent);			if (!valid) return false;
		valid = _insert_rec(crt_time, parent);			if (!valid) return false;
		valid = _insert_process(parent);			    if (!valid) return false;
		return true;
	}

	bool t_gxmlfactory::_init_pceIF_xml_file(const t_gtime& crt_time, pugi::xml_node & parent)
	{
		bool valid = true;
		valid = _insert_gen(crt_time, parent);			if (!valid) return false;
		valid = _insert_gnss(parent);			        if (!valid) return false;
		valid = _insert_inp_PCE_IF(crt_time, parent);	if (!valid) return false;
		valid = _insert_par_PCE(crt_time, parent);		if (!valid) return false;
		valid = _insert_rec(crt_time, parent);			if (!valid) return false;
		valid = _insert_process(parent);			if (!valid) return false;
		
		return true;
	}

	bool t_gxmlfactory::_insert_gen(const t_gtime& crt_time, pugi::xml_node& parent)
	{
		if (parent.child("gen")) parent.remove_child("gen");
		parent.append_child("gen");
		xml_node gen = parent.child("gen");
		gen.append_child(pugi::node_comment).set_value("The part for general settings.");

		gen.append_child(pugi::node_comment).set_value("beg : begin time for process period.");
		gen.append_child(pugi::node_comment).set_value("end : end   time for process period.");
		gen.append_child(pugi::node_comment).set_value("sys : satellite used in process period.");
		gen.append_child(pugi::node_comment).set_value("rec : receiver  used in process period.");
		gen.append_child(pugi::node_comment).set_value("int : sampling  used in process period.");


		if(!gen.child("beg")) gen.append_child("beg");
		if(!gen.child("end")) gen.append_child("end");
		if(!gen.child("sys")) gen.append_child("sys");
		if(!gen.child("rec")) gen.append_child("rec");
		if(!gen.child("int")) gen.append_child("int");

		xml_node gen_beg = gen.child("beg");
		xml_node gen_end = gen.child("end");
		xml_node gen_sys = gen.child("sys");
		xml_node gen_rec = gen.child("rec");
		xml_node gen_int = gen.child("int");

		// beg and end
		gen_beg.text() = crt_time.str_ymdhms().c_str();
		gen_end.text() = (crt_time + _time_interval).str_ymdhms().c_str();
		
		// sys
		std::set<string> sys = dynamic_cast<t_gsetgen*>(_set)->sys();
		gen_sys.append_child(pugi::node_pcdata).set_value(set2str(sys).c_str());
		
		// rec
		std::set<string> rec = dynamic_cast<t_gsetgen*>(_set)->recs();
		gen_rec.append_child(pugi::node_pcdata).set_value(set2str(rec).c_str());

		// intv
		string samp = dbl2str(dynamic_cast<t_gsetgen*>(_set)->sampling());
		gen_int.append_child(pugi::node_pcdata).set_value(samp.c_str());
		return true;
	}

	bool t_gxmlfactory::_insert_gnss(pugi::xml_node& parent)
	{
		if (_sys_list.count("GPS") > 0)
		{
			if (parent.child("gps")) parent.remove_child("gps");
			parent.append_copy(_set->xml_config_node("gps"));

			xml_node gps = parent.child("gps");
			gps.append_child(pugi::node_comment).set_value("GPS  : any GPS satellites ");
			gps.append_child(pugi::node_comment).set_value("sat  : list of GPS satellites: G01 G02 .. or empty(ALL) ");
			gps.append_child(pugi::node_comment).set_value("type : list of GPS  obs types: C L D S P  or empty(ALL) ");
			gps.append_child(pugi::node_comment).set_value("band : list of GPS  obs bands: 1 2 5      or empty(ALL)");
			gps.append_child(pugi::node_comment).set_value("attr : list of PGS attributes: A B C D I L M N P Q S W X Y Z or empty(ALL) ");
		}
		
		if (_sys_list.count("GAL") > 0)
		{
			if (parent.child("gal")) parent.remove_child("gal");
			parent.append_copy(_set->xml_config_node("gal"));
		}

		if (_sys_list.count("BDS") > 0)
		{
			if (parent.child("bds")) parent.remove_child("bds");
			parent.append_copy(_set->xml_config_node("bds"));
		}

		if (_sys_list.count("GLO") > 0)
		{
			if (parent.child("glo")) parent.remove_child("glo");
			parent.append_copy(_set->xml_config_node("glo"));

			xml_node glo = parent.child("glo");
			glo.append_child(pugi::node_comment).set_value("GLO  : any GLONASS satellites");
			glo.append_child(pugi::node_comment).set_value("sat  : list of GPS satellites: R01 R02 .. or empty(ALL)");
			glo.append_child(pugi::node_comment).set_value("type : list of GPS  obs types: C L D S P or empty(ALL)");
			glo.append_child(pugi::node_comment).set_value("band : list of GPS  obs bands: 1 2 3 or empty(ALL)");
			glo.append_child(pugi::node_comment).set_value("attr : list of PGS attributes: A B C D I L M N P Q S W X Y Z or empty(ALL)");

		}
		
		
		return true;
	}

	bool t_gxmlfactory::_insert_process(pugi::xml_node& parent)
	{
		parent.append_copy(_set->xml_config_node("process"));

		xml_node proc = parent.child("process");
		proc.append_child(pugi::node_comment).set_value(" process description :                                                                                     ");
		proc.append_child(pugi::node_comment).set_value(" phase[0, 1]   ..use carrier phase data																	");
		proc.append_child(pugi::node_comment).set_value(" tropo[0, 1]   ..estimate troposphere																		");
		proc.append_child(pugi::node_comment).set_value(" tropo_grad    ..tropospheric horizontal gradient models													");
		proc.append_child(pugi::node_comment).set_value(" tropo_model   ..tropospheric model(SAAS, GPT, ...)														");
		proc.append_child(pugi::node_comment).set_value(" tropo_slant   ..tropospheric slant delays produced														");
		proc.append_child(pugi::node_comment).set_value(" tropo_mf      ..tropospheric mapping function(COSZ, GMF, ...)												");
		proc.append_child(pugi::node_comment).set_value(" grad_mf       ..tropo gradient mapping function(TILTING, CHEN_HERRING, BAR_SEVER ...)						");
		proc.append_child(pugi::node_comment).set_value(" obs_weight    ..observation elevation dependant weighting(EQUAL, SINEL, SINEL2, SINEL4, CODPHA, MLTPTH)	");
		proc.append_child(pugi::node_comment).set_value(" residuals     ..type of residuals(RES_ORIG, RES_NORM, RES_ALL)											");
		proc.append_child(pugi::node_comment).set_value(" sig_init_crd  ..accuracy of initial coordinates[m]														");
		proc.append_child(pugi::node_comment).set_value(" sig_init_ztd  ..accuracy of initial zenith path delay[m]													");
		proc.append_child(pugi::node_comment).set_value(" sig_init_vion ..accuracy of initial vertical iono path delay[m]											");
		proc.append_child(pugi::node_comment).set_value(" minimum_elev  ..elevation angle cut - off[degree]															");
		proc.append_child(pugi::node_comment).set_value(" max_res_norm  ..maximal normalized residuals																");
		
		if (parent.child("process")) return true;
		else return false;
	}

	bool t_gxmlfactory::_insert_inp_POD_IF(const t_gtime& crt_time,pugi::xml_node& parent)
	{
		if (parent.child("inputs")) parent.remove_child("inputs");
		parent.append_child("inputs");
		xml_node inputs = parent.child("inputs");

		inputs.append_child(pugi::node_comment).set_value("The part for input files.       ");
		inputs.append_child(pugi::node_comment).set_value("rinexo  : observation file.     ");
		inputs.append_child(pugi::node_comment).set_value("rinexc  : clock file.           ");
		inputs.append_child(pugi::node_comment).set_value("ics     : orbit intergrate file.");
		inputs.append_child(pugi::node_comment).set_value("blq     : ocean tide file.      ");
		inputs.append_child(pugi::node_comment).set_value("atx     : antannal file.        ");
		inputs.append_child(pugi::node_comment).set_value("orb     : orbit file.           ");
		inputs.append_child(pugi::node_comment).set_value("poleut1 : ERP parameters file.  ");
		inputs.append_child(pugi::node_comment).set_value("sinex   : sinex file.           ");
		inputs.append_child(pugi::node_comment).set_value("biabern : DCB file.             ");
		inputs.append_child(pugi::node_comment).set_value("de      : planet orbit file.    ");
		inputs.append_child(pugi::node_comment).set_value("satpars : satellite parameters. ");

		_insert_inp_rnxo(crt_time, inputs);
		_insert_inp_rnxc(crt_time, inputs);
		_insert_inp_orb(crt_time, inputs);
		_insert_inp_ics(crt_time, inputs);
		_insert_inp_snx(crt_time, inputs);
		_insert_inp_dcb(crt_time, inputs);
		_insert_inp_brd(crt_time, inputs);
		_insert_inp_sp3(crt_time, inputs);
		_insert_inp_sys(inputs);
		return true;
	}

	bool t_gxmlfactory::_insert_inp_PCE_IF(const t_gtime& crt_time, pugi::xml_node& parent)
	{
		if (parent.child("inputs")) parent.remove_child("inputs");
		parent.append_child("inputs");
		xml_node inputs = parent.child("inputs");

		inputs.append_child(pugi::node_comment).set_value("The part for input files.       ");
		inputs.append_child(pugi::node_comment).set_value("rinexo  : observation file.     ");
		inputs.append_child(pugi::node_comment).set_value("rinexc  : sp3 clock file.       ");
		inputs.append_child(pugi::node_comment).set_value("rinexn  : navigation file.      ");
		inputs.append_child(pugi::node_comment).set_value("ics     : orbit intergrate file.");
		inputs.append_child(pugi::node_comment).set_value("blq     : ocean tide file.      ");
		inputs.append_child(pugi::node_comment).set_value("atx     : antannal file.        ");
		inputs.append_child(pugi::node_comment).set_value("orb     : orbit file.           ");
		inputs.append_child(pugi::node_comment).set_value("poleut1 : ERP parameters file.  ");
		inputs.append_child(pugi::node_comment).set_value("sinex   : sinex file.           ");
		inputs.append_child(pugi::node_comment).set_value("biabern : DCB file.             ");
		inputs.append_child(pugi::node_comment).set_value("de      : planet orbit file.    ");
		inputs.append_child(pugi::node_comment).set_value("satpars : satellite parameters. ");

		_insert_inp_rnxo(crt_time, inputs);
		_insert_inp_snx(crt_time, inputs);
		_insert_inp_dcb(crt_time, inputs);
		_insert_inp_brd(crt_time, inputs);
		_insert_inp_sp3(crt_time, inputs);
		_insert_inp_sys(inputs);
		return true;
	}

	bool t_gxmlfactory::_insert_inp_rnxo(const t_gtime& crt_time, pugi::xml_node& parent)
	{
		// crt rinexo path 
		string crt_path = _rinexo_root;
		crt_path = _replace_YYYY(crt_time, crt_path);
		crt_path = _replace_DOY(crt_time, crt_path);

		// rinexo
		parent.remove_child("rinexo");
		parent.append_child("rinexo");

		xml_node inputs_rinexo = parent.child("rinexo");
		bool isChangeLine = false;
		std::set<std::string>::const_iterator beg = _rec_list.begin();
		std::set<std::string>::const_iterator end = _rec_list.end();
		std::set<std::string>::const_iterator crt = beg;
		for (crt; crt != end; crt++)
		{
			std::string tmp;
			string yyyydoy = crt_time.str_yyyydoy();
			string doy = yyyydoy.substr(4, 3);
			string yy = yyyydoy.substr(2, 2);
			string file_name = str2lower(*crt) + doy + "0." + yy + "o";
			if (_isWIN)   tmp = tmp + crt_path + file_name + " ";
			if (_isLINUX) tmp = tmp + crt_path + file_name + " ";
			inputs_rinexo.append_child(pugi::node_pcdata).set_value(tmp.c_str());
		}
		return true;
	}

	bool t_gxmlfactory::_insert_inp_rnxc(const t_gtime& crt_time, pugi::xml_node& parent)
	{
		// crt rinexc path 
		string crt_path = _rinexc_root;
		crt_path = _replace_YYYY(crt_time, crt_path);
		crt_path = _replace_DOY(crt_time, crt_path);

		// rinexc	
		parent.remove_child("rinexc");
		parent.append_child("rinexc");
		xml_node inputs_rinexc = parent.child("rinexc");
		string rinexc1 = crt_path + str2lower(_rinexc_type) + (crt_time - 86400).str_gwkd() + ".clk_30";
		string rinexc2 = crt_path + str2lower(_rinexc_type) + crt_time.str_gwkd() + ".clk_30";
		string rinexc3 = crt_path + str2lower(_rinexc_type) + (crt_time + 86400).str_gwkd() + ".clk_30";

		inputs_rinexc.append_child(pugi::node_pcdata).set_value(rinexc1.c_str());
		inputs_rinexc.append_child(pugi::node_pcdata).set_value(rinexc2.c_str());
		inputs_rinexc.append_child(pugi::node_pcdata).set_value(rinexc3.c_str());
		return true;
	}

	bool t_gxmlfactory::_insert_inp_ics(const t_gtime& crt_time, pugi::xml_node& parent)
	{
		// ics
		parent.remove_child("ics");
		parent.append_child("ics");
		xml_node inputs_ics = parent.child("ics");
		inputs_ics.append_child(pugi::node_pcdata).set_value(sgg_file(crt_time, "ics").c_str());

		return true;
	}

	bool t_gxmlfactory::_insert_inp_orb(const t_gtime& crt_time, pugi::xml_node& parent)
	{
		// orb
		parent.remove_child("orb");
		parent.append_child("orb");
		xml_node inputs_orb = parent.child("orb");
		inputs_orb.append_child(pugi::node_pcdata).set_value(sgg_file(crt_time, "orb").c_str());

		return true;
	}

	bool t_gxmlfactory::_insert_inp_snx(const t_gtime& crt_time, pugi::xml_node& parent)
	{
		// crt proj path 
		string crt_path = _sinex_root;
		crt_path = _replace_YYYY(crt_time, crt_path);
		crt_path = _replace_DOY(crt_time, crt_path);

		parent.remove_child("sinex");
		parent.append_child("sinex");
		xml_node inputs_sinex = parent.child("sinex");
		string name = crt_path + "igs" + crt_time.str_yr() + "P" + crt_time.str_gwk() + ".snx";
		inputs_sinex.append_child(pugi::node_pcdata).set_value(name.c_str());

		return true;

	}

	bool t_gxmlfactory::_insert_inp_dcb(const t_gtime& crt_time, pugi::xml_node& parent)
	{
		// crt proj path 
		string crt_path = _dcb_root;
		crt_path = _replace_YYYY(crt_time, crt_path);
		crt_path = _replace_DOY(crt_time, crt_path);

		parent.remove_child("biabern");
		parent.append_child("biabern");
		xml_node inputs_dcb = parent.child("biabern");
		string name;
		//if (_inp_list.count(BIABERN_INP) == 0)
		//{
		string name1 = crt_path + "P1C1" + crt_time.str_yr() + crt_time.str_mon() + ".DCB";
		string name2 = crt_path + "P1P2" + crt_time.str_yr() + crt_time.str_mon() + ".DCB";
		string name3 = crt_path + "P2C2" + crt_time.str_yr() + crt_time.str_mon() + ".DCB";
		inputs_dcb.append_child(pugi::node_pcdata).set_value(name1.c_str());
		inputs_dcb.append_child(pugi::node_pcdata).set_value(name2.c_str());
		inputs_dcb.append_child(pugi::node_pcdata).set_value(name3.c_str());
		//}
		//else
		//{
		//	//方式2
		//	multimap<IFMT, string>::const_iterator beg;
		//	multimap<IFMT, string>::const_iterator end;
		//	multimap<IFMT, string>::const_iterator iter;
		//	beg = _inp_list.lower_bound(BIABERN_INP);
		//	end = _inp_list.upper_bound(BIABERN_INP);
		//	for (iter = beg; iter != end; iter++)
		//	{
		//		inputs_sinex.append_child(pugi::node_pcdata).set_value((*iter).second.substr(7).c_str());
		//	}
		//}

		return true;
	}

	bool t_gxmlfactory::_insert_inp_sp3(const t_gtime& crt_time, pugi::xml_node& parent)
	{
		// crt rinexc path 
		string crt_path = _sp3_root;
		crt_path = _replace_YYYY(crt_time, crt_path);
		crt_path = _replace_DOY(crt_time, crt_path);

		// rinexc	
		parent.remove_child("sp3");
		parent.append_child("sp3");
		xml_node inputs_sp3 = parent.child("sp3");
		string sp3_1 = crt_path + str2lower(_sp3_type) + (crt_time - 86400).str_gwkd() + ".sp3";
		string sp3_2 = crt_path + str2lower(_sp3_type) + crt_time.str_gwkd() + ".sp3";
		string sp3_3 = crt_path + str2lower(_sp3_type) + (crt_time + 86400).str_gwkd() + ".sp3";

		inputs_sp3.append_child(pugi::node_pcdata).set_value(sp3_1.c_str());
		inputs_sp3.append_child(pugi::node_pcdata).set_value(sp3_2.c_str());
		inputs_sp3.append_child(pugi::node_pcdata).set_value(sp3_3.c_str());
		return true;
	}

	bool t_gxmlfactory::_insert_inp_brd(const t_gtime& crt_time, pugi::xml_node& parent)
	{
		// crt proj path 
		string crt_path = _brd_root;
		crt_path = _replace_YYYY(crt_time, crt_path);
		crt_path = _replace_DOY(crt_time, crt_path);

		parent.remove_child("rinexn");
		parent.append_child("rinexn");
		xml_node inputs_brd = parent.child("rinexn");
		string name = crt_path + "brdm" + crt_time.str_doy() + "0." + crt_time.str_yr() + "p";
		inputs_brd.append_child(pugi::node_pcdata).set_value(name.c_str());

		return true;
	}

	bool t_gxmlfactory::_insert_inp_sys(pugi::xml_node& parent)
	{
		parent.remove_child("blq");
		parent.remove_child("poleut1");
		parent.remove_child("satpars");
		parent.remove_child("de");
		parent.remove_child("atx");

		parent.append_child("blq");
		parent.append_child("poleut1");
		parent.append_child("satpars");
		parent.append_child("de");
		parent.append_child("atx");
		xml_node inputs_blq     = parent.child("blq");
		xml_node inputs_poleut1 = parent.child("poleut1");
		xml_node inputs_satpars = parent.child("satpars");
		xml_node inputs_de      = parent.child("de");
		xml_node inputs_atx     = parent.child("atx");

		multimap<IFMT, string>::const_iterator itINP = _inp_list.begin();
		for (size_t i = 0; i < _inp_list.size() && itINP != _inp_list.end(); ++i, ++itINP)
		{
			// Get the file format/path, which will be used in decoder
			IFMT   ifmt(itINP->first);
			string path(itINP->second);
			switch (ifmt)
			{
			case ATX_INP :
				inputs_atx.append_child(pugi::node_pcdata).set_value((_sysfile_root + path.substr(7)).c_str());
				break;
			case BLQ_INP:
				inputs_blq.append_child(pugi::node_pcdata).set_value((_sysfile_root + path.substr(7)).c_str());
				break;
			case POLEUT1_INP:
				inputs_poleut1.append_child(pugi::node_pcdata).set_value((_sysfile_root + path.substr(7)).c_str());
				break;
			case SATPARS_INP:
				inputs_satpars.append_child(pugi::node_pcdata).set_value((_sysfile_root + path.substr(7)).c_str());
				break;
			case DE_INP:
				inputs_de.append_child(pugi::node_pcdata).set_value((_sysfile_root + path.substr(7)).c_str());
				break;
			case SINEX_INP:
			{
				parent.remove_child("sinex");
				parent.append_child("sinex");
				xml_node inputs_sinex = parent.child("sinex");
				inputs_sinex.append_child(pugi::node_pcdata).set_value((_sysfile_root + path.substr(7)).c_str());
				break;
			}
			default:
				break;
			}
		}

		if (_inp_list.count(DE_INP) == 0)
		{
			if (_isWIN) inputs_de.append_child(pugi::node_pcdata).set_value((_sysfile_root + _list_de_win).c_str());
			else inputs_de.append_child(pugi::node_pcdata).set_value((_sysfile_root + _list_de_lux).c_str());
		}

		if (_inp_list.count(ATX_INP) == 0)
		{
			inputs_atx.append_child(pugi::node_pcdata).set_value((_sysfile_root + _list_atx).c_str());
		}

		if (_inp_list.count(BLQ_INP) == 0)
		{
			inputs_blq.append_child(pugi::node_pcdata).set_value((_sysfile_root + _list_blq).c_str());
		}

		if (_inp_list.count(POLEUT1_INP) == 0)
		{
			inputs_poleut1.append_child(pugi::node_pcdata).set_value((_sysfile_root + _list_poleut1).c_str());
		}

		if (_inp_list.count(SATPARS_INP) == 0)
		{
			inputs_satpars.append_child(pugi::node_pcdata).set_value((_sysfile_root + _list_satpars).c_str());
		}

		return true;
	}


	bool t_gxmlfactory::_insert_par(const t_gtime& crt_time, pugi::xml_node& parent)
	{
		parent.remove_child("parameters");
		parent.append_child("parameters");
		xml_node pars = parent.child("parameters");

		pars.append_child(pugi::node_comment).set_value("The part for par settings.                        ");
		pars.append_child(pugi::node_comment).set_value("ID        : for index, PRN or Rec name.           ");
		pars.append_child(pugi::node_comment).set_value("sigCLK    : sigma for clk pars.                   ");
		pars.append_child(pugi::node_comment).set_value("sigPOS    : sigma for xyz pars.                   ");
		pars.append_child(pugi::node_comment).set_value("sigZTD    : sigma for ztd pars.                   ");
		pars.append_child(pugi::node_comment).set_value("sigVion   : sigma for vertical direction ion pars.");
		pars.append_child(pugi::node_comment).set_value("sigSion   : sigma for inclined direction ion pars.");
		pars.append_child(pugi::node_comment).set_value("sigIonoPd : sigma for ion par period.             ");
		pars.append_child(pugi::node_comment).set_value("sigTropPd : sigma for ztd par period.             ");
		pars.append_child(pugi::node_comment).set_value("<STA ID = \"XXXX\" sigCLK = \"9000\" sigPOS = \"0.100_0.100_0.100\" sigZTD = \"0.201\" sigTropPd = \"0.015\" sigVion = \"1\" sigSion = \"0\" sigIonoPd = \"0.015\" / >");
		pars.append_child(pugi::node_comment).set_value("<SAT ID = \"XXX\" sigCLK = \"5000\" sigPOS = \"10.000_10.000_10.000\" sigVEL = \"0.100_0.100_0.100\" sigECOM = \"10.000_10.000_10.000_10.000_10.000_10.000_10.000_10.000_10.000\" / >");

		t_gsetpar* set_par = dynamic_cast<t_gsetpar*>(_set);
		if (set_par == nullptr) return false;

		// for ERO =========================================>
		xml_node par_erp = pars.append_child("ERP");
		pugi::xml_attribute tmp_attri;
		tmp_attri = par_erp.append_attribute("sigXPOLE");
		double sigXPOLE = set_par->sigXpole();
		tmp_attri.set_value(sigXPOLE);
		tmp_attri = par_erp.append_attribute("sigYPOLE");
		double sigYPOLE = set_par->sigYpole();
		tmp_attri.set_value(sigYPOLE);
		tmp_attri = par_erp.append_attribute("sigUT1");
		double sigUT1 = set_par->sigUt1();
		tmp_attri.set_value(sigUT1);
		tmp_attri = par_erp.append_attribute("sigDXPOLE");
		double sigDXPOLE = set_par->sigDxpole();
		tmp_attri.set_value(sigDXPOLE);
		tmp_attri = par_erp.append_attribute("sigDYPOLE");
		double sigDYPOLE = set_par->sigDypole();
		tmp_attri.set_value(sigDYPOLE);
		tmp_attri = par_erp.append_attribute("sigDUT1");
		double sigDUT1 = set_par->sigDut1();
		tmp_attri.set_value(sigDUT1);

		set<string>::const_iterator beg;
		set<string>::const_iterator end;
		set<string>::const_iterator tmp;

		// for STA ===============================================>
		beg = _rec_list.begin();
		end = _rec_list.end();
		tmp = beg;
		for (tmp; tmp != end; tmp++)
		{
			string name = str2upper(*tmp);
			xml_node par_sta = pars.append_child("STA");
			pugi::xml_attribute tmp_attri;

			tmp_attri = par_sta.append_attribute("ID");
			tmp_attri.set_value(name.c_str());

			tmp_attri = par_sta.append_attribute("sigCLK");
			double sigRecCLK = set_par->sigRecCLK("XXXX");
			tmp_attri.set_value(sigRecCLK);

			tmp_attri = par_sta.append_attribute("sigPOS");
			map<string, double> sigPOS = set_par->sigRecPos("XXXX");
			string sigXYZ = trim(dbl2str(sigPOS["sigPX"])) + "_" +
				trim(dbl2str(sigPOS["sigPY"])) + "_" +
				trim(dbl2str(sigPOS["sigPZ"]));
			tmp_attri.set_value(sigXYZ.c_str());

			tmp_attri = par_sta.append_attribute("sigZTD");
			double sigZTD = set_par->sigZtd("XXXX");
			tmp_attri.set_value(sigZTD);

			tmp_attri = par_sta.append_attribute("sigTropPd");
			double sigTropPd = set_par->sigTropPd("XXXX");
			tmp_attri.set_value(sigTropPd);

			tmp_attri = par_sta.append_attribute("sigVion");
			double sigVion = set_par->sigVion("XXXX");
			tmp_attri.set_value(sigVion);

			tmp_attri = par_sta.append_attribute("sigSion");
			double sigSion = set_par->sigSion("XXXX");
			tmp_attri.set_value(sigSion);

			tmp_attri = par_sta.append_attribute("sigIonoPd");
			double sigIonoPd = set_par->sigIonoPd("XXXX");
			tmp_attri.set_value(sigIonoPd);
		}

		// for SAT ===============================================>
		beg = _sat_list.begin();
		end = _sat_list.end();
		tmp = beg;
		for (tmp; tmp != end; tmp++)
		{
			string name = str2upper(*tmp);
			xml_node par_sat = pars.append_child("SAT");
			pugi::xml_attribute tmp_attri;

			tmp_attri = par_sat.append_attribute("ID");
			tmp_attri.set_value(name.c_str());

			tmp_attri = par_sat.append_attribute("sigCLK");
			double sigSatCLK = set_par->sigSatCLK("XXX");
			tmp_attri.set_value(sigSatCLK);

			tmp_attri = par_sat.append_attribute("sigPOS");
			map<string, double> sigPOS = set_par->sigSatPos("XXX");
			string sigXYZ = trim(dbl2str(sigPOS["sigPX"])) + "_" +
				trim(dbl2str(sigPOS["sigPY"])) + "_" +
				trim(dbl2str(sigPOS["sigPZ"]));
			tmp_attri.set_value(sigXYZ.c_str());

			tmp_attri = par_sat.append_attribute("sigVEL");
			map<string, double> sigVEL = set_par->sigSatVel("XXX");
			string sigVXYZ = trim(dbl2str(sigVEL["sigVX"])) + "_" +
				trim(dbl2str(sigVEL["sigVY"])) + "_" +
				trim(dbl2str(sigVEL["sigVZ"]));
			tmp_attri.set_value(sigVXYZ.c_str());

			tmp_attri = par_sat.append_attribute("sigECOM");
			map<string, double> sigECOM = set_par->sigSatEcom("XXX");
			string sigSolar = trim(dbl2str(sigECOM["sigD0"])) + "_" +
				trim(dbl2str(sigECOM["sigDc"])) + "_" +
				trim(dbl2str(sigECOM["sigDs"])) + "_" +
				trim(dbl2str(sigECOM["sigY0"])) + "_" +
				trim(dbl2str(sigECOM["sigYc"])) + "_" +
				trim(dbl2str(sigECOM["sigYs"])) + "_" +
				trim(dbl2str(sigECOM["sigB0"])) + "_" +
				trim(dbl2str(sigECOM["sigBc"])) + "_" +
				trim(dbl2str(sigECOM["sigBs"]));
			tmp_attri.set_value(sigSolar.c_str());
		}

		return true;
	}

	bool t_gxmlfactory::_insert_par_PCE(const t_gtime& crt_time, pugi::xml_node& parent)
	{
		parent.remove_child("parameters");
		parent.append_child("parameters");
		xml_node pars = parent.child("parameters");

		pars.append_child(pugi::node_comment).set_value("The part for par settings.                        ");
		pars.append_child(pugi::node_comment).set_value("ID        : for index, PRN or Rec name.           ");
		pars.append_child(pugi::node_comment).set_value("sigCLK    : sigma for clk pars.                   ");
		pars.append_child(pugi::node_comment).set_value("sigPOS    : sigma for xyz pars.                   ");
		pars.append_child(pugi::node_comment).set_value("sigZTD    : sigma for ztd pars.                   ");
		pars.append_child(pugi::node_comment).set_value("sigVion   : sigma for vertical direction ion pars.");
		pars.append_child(pugi::node_comment).set_value("sigSion   : sigma for inclined direction ion pars.");
		pars.append_child(pugi::node_comment).set_value("sigIonoPd : sigma for ion par period.             ");
		pars.append_child(pugi::node_comment).set_value("sigTropPd : sigma for ztd par period.             ");
		pars.append_child(pugi::node_comment).set_value("<STA ID = \"XXXX\" sigCLK = \"9000\" sigPOS = \"0.100_0.100_0.100\" sigZTD = \"0.201\" sigTropPd = \"0.015\" sigVion = \"1\" sigSion = \"0\" sigIonoPd = \"0.015\" / >");
		pars.append_child(pugi::node_comment).set_value("<SAT ID = \"XXX\" sigCLK = \"5000\" sigPOS = \"10.000_10.000_10.000\" sigVEL = \"0.100_0.100_0.100\" sigECOM = \"10.000_10.000_10.000_10.000_10.000_10.000_10.000_10.000_10.000\" / >");

		t_gsetpar* set_par = dynamic_cast<t_gsetpar*>(_set);
		if (set_par == nullptr) return false;

		set<string>::const_iterator beg;
		set<string>::const_iterator end;
		set<string>::const_iterator tmp;

		// for STA ===============================================>
		beg = _rec_list.begin();
		end = _rec_list.end();
		tmp = beg;
		for (tmp; tmp != end; tmp++)
		{
			string name = str2upper(*tmp);
			xml_node par_sta = pars.append_child("STA");
			pugi::xml_attribute tmp_attri;

			tmp_attri = par_sta.append_attribute("ID");
			tmp_attri.set_value(name.c_str());

			tmp_attri = par_sta.append_attribute("sigCLK");
			double sigRecCLK = set_par->sigRecCLK("XXXX");
			tmp_attri.set_value(sigRecCLK);

			tmp_attri = par_sta.append_attribute("sigPOS");
			map<string, double> sigPOS = set_par->sigRecPos("XXXX");
			string sigXYZ = trim(dbl2str(sigPOS["sigPX"])) + "_" +
				trim(dbl2str(sigPOS["sigPY"])) + "_" +
				trim(dbl2str(sigPOS["sigPZ"]));
			tmp_attri.set_value(sigXYZ.c_str());

			tmp_attri = par_sta.append_attribute("sigZTD");
			double sigZTD = set_par->sigZtd("XXXX");
			tmp_attri.set_value(sigZTD);

			tmp_attri = par_sta.append_attribute("sigTropPd");
			double sigTropPd = set_par->sigTropPd("XXXX");
			tmp_attri.set_value(sigTropPd);

			tmp_attri = par_sta.append_attribute("sigVion");
			double sigVion = set_par->sigVion("XXXX");
			tmp_attri.set_value(sigVion);

			tmp_attri = par_sta.append_attribute("sigSion");
			double sigSion = set_par->sigSion("XXXX");
			tmp_attri.set_value(sigSion);

			tmp_attri = par_sta.append_attribute("sigIonoPd");
			double sigIonoPd = set_par->sigIonoPd("XXXX");
			tmp_attri.set_value(sigIonoPd);
		}

		// for SAT ===============================================>
		beg = _sat_list.begin();
		end = _sat_list.end();
		tmp = beg;
		for (tmp; tmp != end; tmp++)
		{
			string name = str2upper(*tmp);
			xml_node par_sat = pars.append_child("SAT");
			pugi::xml_attribute tmp_attri;

			tmp_attri = par_sat.append_attribute("ID");
			tmp_attri.set_value(name.c_str());

			tmp_attri = par_sat.append_attribute("sigCLK");
			double sigSatCLK = set_par->sigSatCLK("XXX");
			tmp_attri.set_value(sigSatCLK);

			tmp_attri = par_sat.append_attribute("sigPOS");
			map<string, double> sigPOS = set_par->sigSatPos("XXX");
			string sigXYZ = trim(dbl2str(sigPOS["sigPX"])) + "_" +
				trim(dbl2str(sigPOS["sigPY"])) + "_" +
				trim(dbl2str(sigPOS["sigPZ"]));
			tmp_attri.set_value(sigXYZ.c_str());

			tmp_attri = par_sat.append_attribute("sigVEL");
			map<string, double> sigVEL = set_par->sigSatVel("XXX");
			string sigVXYZ = trim(dbl2str(sigVEL["sigVX"])) + "_" +
				trim(dbl2str(sigVEL["sigVY"])) + "_" +
				trim(dbl2str(sigVEL["sigVZ"]));
			tmp_attri.set_value(sigVXYZ.c_str());

			tmp_attri = par_sat.append_attribute("sigECOM");
			map<string, double> sigECOM = set_par->sigSatEcom("XXX");
			string sigSolar = trim(dbl2str(sigECOM["sigD0"])) + "_" +
				trim(dbl2str(sigECOM["sigDc"])) + "_" +
				trim(dbl2str(sigECOM["sigDs"])) + "_" +
				trim(dbl2str(sigECOM["sigY0"])) + "_" +
				trim(dbl2str(sigECOM["sigYc"])) + "_" +
				trim(dbl2str(sigECOM["sigYs"])) + "_" +
				trim(dbl2str(sigECOM["sigB0"])) + "_" +
				trim(dbl2str(sigECOM["sigBc"])) + "_" +
				trim(dbl2str(sigECOM["sigBs"]));
			tmp_attri.set_value(sigSolar.c_str());
		}

		return true;
	}

	bool t_gxmlfactory::_insert_rec(const t_gtime& crt_time, pugi::xml_node& parent)
	{
		// remove for case
		if (parent.child("receiver")) parent.remove_child("receiver");

		parent.append_child("receiver");
		xml_node pars = parent.child("receiver");

		pars.append_child(pugi::node_comment).set_value("The part for receiver settings. ");
		pars.append_child(pugi::node_comment).set_value("id   : for index, Rec name.     ");
		pars.append_child(pugi::node_comment).set_value("X    : X pos.                   ");
		pars.append_child(pugi::node_comment).set_value("Y    : Y pos.                   ");
		pars.append_child(pugi::node_comment).set_value("Z    : Z pos.                   ");
	

		// Read sinex file =======================================================>
		// crt snx path 
		string crt_path = _sinex_root;
		crt_path = _replace_YYYY(crt_time, crt_path);
		crt_path = _replace_DOY(crt_time, crt_path);
		
		multimap<IFMT, string>::const_iterator itSNX = _inp_list.find(SINEX_INP);
		string name;
		t_gtime beg_time;
		if (itSNX != _inp_list.end())
		{
			name = crt_path + itSNX->second.substr(7);
			int gwk = str2int(itSNX->second.substr(13, 4));
			t_gtime snxtime(gwk, 0, 0, 0);
			beg_time = snxtime;
		}
		else
		{
			name = crt_path + "igs" + crt_time.str_yr() + "P" + crt_time.str_gwk() + ".snx";
		}
		t_gdata*  gdata = 0;
		t_gcoder* gcoder = new t_sinex(_set, "", 20480);
		t_gallobj* gobj  = new t_gallobj(nullptr, nullptr); gobj->glog(_log);
		if (gcoder) 
		{

			// prepare gio for the file
			t_gio*  gio = new t_gfile;
			gio->glog(_log);
			gio->path("file://"+ name);

			// Put the file into gcoder
			gcoder->clear();
			gcoder->path("file://" + name);
			gcoder->glog(_log);
			// Put the data container into gcoder
			gcoder->add_data("ID1", gdata);
			gcoder->add_data("OBJ", gobj);
			// Put the gcoder into the gio
			// Note, gcoder contain the gdata and gio contain the gcoder
			gio->coder(gcoder);

			// Read the data from file here
			gio->run_read();
			if(gio) delete gio; gio = nullptr;
			
		}

		
		set<string>::const_iterator beg;
		set<string>::const_iterator end;
		set<string>::const_iterator tmp;

		// for STA ===============================================>
		set<string> rec_tmp = gobj->objects_set(t_gdata::REC);
		beg = _rec_list.begin();
		end = _rec_list.end();
		tmp = beg;
		int valid_sta_num = 0;
		for (tmp; tmp != end; tmp++)
		{

			string rec = str2upper(*tmp);
			shared_ptr<t_gobj> grec = gobj->obj(rec);

			if (!grec) continue;

			t_gtriple rec_crd;
			if (itSNX != _inp_list.end())
			{
				 rec_crd = grec->crd(beg_time);
			}
			else
			{
				 rec_crd = grec->crd(_beg_time);
			}

			if (double_eq(rec_crd.crd(0), 0.0) &&
				double_eq(rec_crd.crd(1), 0.0) &&
				double_eq(rec_crd.crd(2), 0.0))
			{
				continue;
			}

			// one station valid
			valid_sta_num++;

			xml_node par_sta = pars.append_child("rec");
			pugi::xml_attribute tmp_attri;
			tmp_attri = par_sta.append_attribute("id");
			tmp_attri.set_value(rec.c_str());

			string crd_x = format("%16.6lf", rec_crd[0]);
			tmp_attri = par_sta.append_attribute("X");
			tmp_attri.set_value(crd_x.c_str());

			string crd_y = format("%16.6lf", rec_crd[1]);
			tmp_attri = par_sta.append_attribute("Y");
			tmp_attri.set_value(crd_y.c_str());

			string crd_z = format("%16.6lf", rec_crd[2]);
			tmp_attri = par_sta.append_attribute("Z");
			tmp_attri.set_value(crd_z.c_str());

		}

		if (gcoder) delete gcoder; gcoder = nullptr;
		if (gdata) delete gdata; gdata = nullptr;
		if (gobj) delete gobj; gobj = nullptr;

		// no station valid
		if (valid_sta_num == 0) return false;

		return true;
	}

	bool t_gxmlfactory::_insert_root()
	{
		return false;
	}

	int32_t t_gxmlfactory::createDirectory(const std::string & directoryPath)
	{
		uint32_t dirPathLen = directoryPath.length();
		if (dirPathLen > MAX_PATH_LEN)
		{
			return -1;
		}
		
		for (uint32_t i = 0; i < dirPathLen; ++i)
		{
			if (directoryPath[i] == '\\' || directoryPath[i] == '/')
			{
				string tmp = directoryPath.substr(0, i + 1);
				if (ACCESS(tmp.c_str(), 0) != 0)
				{
					int32_t ret = MKDIR(tmp.c_str());
					if (ret != 0)
					{
						return ret;
					}
				}
			}
		}
		return 0;
	}

	//string t_gxmlfactory::recs2rnxos(const t_gtime & crt_time, const set<string>& recs)
	//{
	//	bool isChangeLine = false;
	//	std::string tmp;
	//	std::set<std::string>::const_iterator beg = recs.begin();
	//	std::set<std::string>::const_iterator end = recs.end();
	//	std::set<std::string>::const_iterator crt = beg;

	//	if (recs.size() > 1)  isChangeLine = true;
	//	int num = 0;
	//	for (crt; crt != end; crt++)
	//	{
	//		if (isChangeLine) tmp = tmp + "\n\t\t\t";
	//		string yyyydoy = crt_time.str_yyyydoy();
	//		string doy = yyyydoy.substr(4, 3);
	//		string yy  = yyyydoy.substr(2, 2);
	//		string file_name = str2lower(*crt) + doy + "0." + yy + "o";
	//		num++;
	//		if(_isWIN)   tmp = tmp + _rinexo_root + file_name + " ";
	//		if(_isLINUX) tmp = tmp + _rinexo_root + file_name + " ";
	//	}
	//	tmp = tmp + "\n\t\t";
	//	return tmp;
	//}

	//string t_gxmlfactory::three_day_rnxcs(const t_gtime & crt_time)
	//{
	//	string rinexc;
	//	rinexc = rinexc + "\n\t\t\t" + _rinexc_root + str2lower(_rinexc_type) + (crt_time - 86400).str_gwkd() + ".clk_30";
	//	rinexc = rinexc + "\n\t\t\t" + _rinexc_root + str2lower(_rinexc_type) + crt_time.str_gwkd() + ".clk_30";
	//	rinexc = rinexc + "\n\t\t\t" + _rinexc_root + str2lower(_rinexc_type) + (crt_time + 86400).str_gwkd() + ".clk_30";
	//	rinexc = rinexc + "\n\t\t";
	//	return rinexc;
	//}

	string t_gxmlfactory::sgg_file(const t_gtime& crt_time, const string & type)
	{
		// crt proj path 
		string crt_path = _project_root;
		crt_path = _replace_YYYY(crt_time, crt_path);
		crt_path = _replace_DOY(crt_time, crt_path);

		string tmp;
		tmp = tmp + crt_path + "sgg_" + crt_time.str_yyyydoy() + "." + str2lower(type);
		return tmp;

	}

	string t_gxmlfactory::_replace_YYYY(const t_gtime& crt_time, const string& path)
	{
		// crt proj path 
		string crt_path = path;
		string yyyy = crt_time.str_yyyy();
		if (crt_path.find("+YYYY+") != string::npos)
		{
			crt_path.replace(crt_path.find_first_of("+"), 6, yyyy);
		}

		return crt_path;
	}

	string t_gxmlfactory::_replace_DOY(const t_gtime& crt_time, const string& path)
	{
		string crt_path = path;
		string doy = crt_time.str_doy();
		if (crt_path.find("+DOY+") != string::npos)
		{
			crt_path.replace(crt_path.find_first_of("+"), 5, doy);
		}
		return crt_path;

	}
}


