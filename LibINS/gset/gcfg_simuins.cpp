#include "gcfg_simuins.h"

gsins::t_gcfg_simuins::t_gcfg_simuins() :
	t_gcfg_ppp(),
	t_gsetins(),
	t_gsetign()
{
	_IFMT_supported.clear();
	_IFMT_supported.insert(RINEXC_INP);
	_IFMT_supported.insert(RINEXN_INP);
	_IFMT_supported.insert(ATX_INP);
	_IFMT_supported.insert(BLQ_INP);
	_IFMT_supported.insert(SP3_INP);
	_IFMT_supported.insert(BIABERN_INP);
	_IFMT_supported.insert(BIASINEX_INP);
	_IFMT_supported.insert(POLEUT1_INP);
	_IFMT_supported.insert(LEAPSECOND_INP);
	_IFMT_supported.insert(DE_INP);
	_IFMT_supported.insert(BIASINEX_INP);

	_IFMT_supported.insert(IMU_POS_INP);

	_OFMT_supported.insert(LOG_OUT);
	_OFMT_supported.insert(PPP_OUT);
	_OFMT_supported.insert(RINEXO_OUT);
	_OFMT_supported.insert(IMU_OUT);
	_OFMT_supported.insert(INS_OUT);

}

gsins::t_gcfg_simuins::~t_gcfg_simuins()
{
}

void gsins::t_gcfg_simuins::check()
{
	t_gcfg_ppp::check();
	t_gsetins::check();
	t_gsetign::check();
}

void gsins::t_gcfg_simuins::help()
{
	t_gcfg_ppp::help();
	t_gsetins::help();
	t_gsetign::help();
}

int gsins::t_gcfg_simuins::generateXML(stringstream & ss)
{
	char buffer[20480]="";

	t_gtime beg = dynamic_cast<t_gsetgen*>(this)->beg();
	t_gtime end = dynamic_cast<t_gsetgen*>(this)->end();

	ss << buffer; sprintf(buffer,"<?xml version = \"1.0\" encoding = \"UTF-8\" standalone = \"yes\" ?>\n");
	ss << buffer; sprintf(buffer,"	<!DOCTYPE config>																									\n");
	ss << buffer; sprintf(buffer,"	<config>																											\n");
	ss << buffer; sprintf(buffer,"																														\n");
	ss << buffer; sprintf(buffer,"	<gen>																												\n");
	ss << buffer; sprintf(buffer,"	<beg> \"%4d-%02d-%02d %02d:%02d:%02d\" </beg>																		\n", (int)beg.year(), (int)beg.mon(), (int)beg.day(), (int)beg.hour(), (int)beg.mins(), (int)beg.secs());
	ss << buffer; sprintf(buffer,"	<end> \"%4d-%02d-%02d %02d:%02d:%02d\" </end>																		\n", (int)end.year(), (int)end.mon(), (int)end.day(), (int)end.hour(), (int)end.mins(), (int)end.secs());
	ss << buffer; sprintf(buffer,"	<sys> GPS GAL BDS</sys>																								\n");
	ss << buffer; sprintf(buffer,"	<rec> %s </rec>																										\n", dynamic_cast<t_gsetgen*>(this)->recs().begin()->c_str());
	ss << buffer; sprintf(buffer,"	<int>  %.2lf </int>																									\n", dynamic_cast<t_gsetgen*>(this)->sampling());
	ss << buffer; sprintf(buffer,"	</gen>																												\n");
	ss << buffer; sprintf(buffer,"																														\n");
	ss << buffer; sprintf(buffer,"	<inputs>																											\n");
	if (this->outputs("rinexo").size() > 0)
	{
		ss << buffer; sprintf(buffer, "	<rinexo> %s </rinexo>																							\n", this->outputs("rinexo").c_str());
	}
	if (this->outputs("imu").size() > 0)
	{
		ss << buffer; sprintf(buffer, "	<imu> %s </imu>																									\n", this->outputs("imu").c_str());
	}
	for (auto ifmt : _IFMT_supported)
	{
		string tmp="";
		vector<string> inp_vec=this->inputs(ifmt);
		for (auto inp_file : inp_vec) {tmp += inp_file + " ";}
		string ifmt_str = ifmt2str(ifmt);
		transform(ifmt_str.begin(), ifmt_str.end(), ifmt_str.begin(), ::tolower);
		if (tmp != "")
		{
			ss << buffer; sprintf(buffer, "	<%s> %s </%s>																								\n", ifmt_str.c_str(), tmp.c_str(), ifmt_str.c_str());
		}
	}
	ss << buffer; sprintf(buffer,"	</inputs>																											\n");
	ss << buffer; sprintf(buffer,"																														\n");
	ss << buffer; sprintf(buffer,"	<outputs append = \"false\" verb = \"2\">																			\n");
	ss << buffer; sprintf(buffer,"	<log> LOGRT.log </log>																								\n");
	ss << buffer; sprintf(buffer,"	<ppp> $(rec)-TCI.ppp </ppp>																							\n");
	ss << buffer; sprintf(buffer,"	<flt> $(rec)-TCI.flt </flt>																							\n");
	ss << buffer; sprintf(buffer,"	<ins> great-TCI.ins  </ins>																							\n");
	ss << buffer; sprintf(buffer,"	</outputs>																											\n");
	ss << buffer; sprintf(buffer,"																														\n");
	ss << buffer; sprintf(buffer,"	<process																											\n");
	ss << buffer; sprintf(buffer,"	phase = \"true\"																									\n");
	ss << buffer; sprintf(buffer,"	tropo = \"true\"																									\n");
	ss << buffer; sprintf(buffer,"	iono = \"true\"																										\n");
	ss << buffer; sprintf(buffer,"	tropo_model = \"saastamoinen\"																						\n");
	ss << buffer; sprintf(buffer,"	tropo_mf = \"gmf\"																									\n");
	ss << buffer; sprintf(buffer,"	sig_init_crd = \"100\"																								\n");
	ss << buffer; sprintf(buffer,"	sig_init_vel = \"100\"																								\n");
	ss << buffer; sprintf(buffer,"	sig_init_ztd = \"100\"																								\n");
	ss << buffer; sprintf(buffer,"	sig_init_amb = \"30\"																								\n");
	ss << buffer; sprintf(buffer,"	minimum_elev = \"10\"																								\n");
	ss << buffer; sprintf(buffer,"	obs_combination = \"IONO_FREE\"																						\n");
	ss << buffer; sprintf(buffer,"	residuals = \"RES_NORM\"																							\n");
	ss << buffer; sprintf(buffer,"	max_res_norm = \"3\"																								\n");
	ss << buffer; sprintf(buffer,"	auto_band = \"false\"																								\n");
	ss << buffer; sprintf(buffer,"	pos_kin = \"true\"																									\n");
	ss << buffer; sprintf(buffer,"	realtime = \"false\"																								\n");
	ss << buffer; sprintf(buffer,"	obs_weight = \"SINEL\"																								\n");
	ss << buffer; sprintf(buffer,"	>																													\n");
	ss << buffer; sprintf(buffer,"	<freq_model> double </freq_model>																					\n");
	ss << buffer; sprintf(buffer,"	<ifb_model> none </ifb_model>																						\n");
	ss << buffer; sprintf(buffer,"	</process>																											\n");
	ss << buffer; sprintf(buffer,"																														\n");
	ss << buffer; sprintf(buffer,"	<filter																												\n");
	ss << buffer; sprintf(buffer,"	method_flt = \"kalman\"																								\n");
	ss << buffer; sprintf(buffer,"	method_smt = \"rtssvd\"																								\n");
	ss << buffer; sprintf(buffer,"	noise_crd = \"100\"																									\n");
	ss << buffer; sprintf(buffer,"	noise_vel = \"5\"																									\n");
	ss << buffer; sprintf(buffer,"	noise_dclk = \"100\"																								\n");
	ss << buffer; sprintf(buffer,"	rndwk_ztd = \"60\"																								\n");
	ss << buffer; sprintf(buffer,"	noise_vion = \"0.1\"																								\n");
	ss << buffer; sprintf(buffer,"	smt_delay = \"0\"																									\n");
	ss << buffer; sprintf(buffer,"	/>																													\n");
	ss << buffer; sprintf(buffer,"																														\n");
	ss << buffer; sprintf(buffer,"	<ins																												\n");
	ss << buffer; sprintf(buffer,"	align = \"true\"																									\n");
	ss << buffer; sprintf(buffer,"	align_time = \"50\"																									\n");
	ss << buffer; sprintf(buffer,"	start = \"%.3lf\"																									\n", (double)beg.sow());
	ss << buffer; sprintf(buffer,"	end = \"%.3lf\"																										\n", (double)end.sow());
	ss << buffer; sprintf(buffer,"	freq = \"%d\"																										\n", dynamic_cast<t_gsetins*>(this)->freq());
	ss << buffer; sprintf(buffer,"	order = \"%s\"																										\n", dynamic_cast<t_gsetins*>(this)->order().c_str());
	ss << buffer; sprintf(buffer,"	out_freq = \"%d\"																									\n", dynamic_cast<t_gsetins*>(this)->out_freq());
	ss << buffer; sprintf(buffer,"	out_order = \"%s\"																									\n", ((string)_doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("out_order").value()).c_str());
	ss << buffer; sprintf(buffer,"	compensation = \"poly\"																								\n");
	ss << buffer; sprintf(buffer,"	subsample = \"1\"																									\n");
	ss << buffer; sprintf(buffer,"	gyro_unit = \"%s\"																									\n", ((string)_doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("gyro_unit").value()).c_str());
	ss << buffer; sprintf(buffer,"	acce_unit = \"%s\"																									\n", ((string)_doc.child(XMLKEY_ROOT).child(XMLKEY_INS).attribute("acce_unit").value()).c_str());
	ss << buffer; sprintf(buffer,"	/>																													\n");
	ss << buffer; sprintf(buffer,"																														\n");
	ss << buffer; sprintf(buffer,"	<integration																										\n");
	ss << buffer; sprintf(buffer,"	ign_type = \"TCI\"																									\n");
	ss << buffer; sprintf(buffer,"	lever = \"%lf, %lf, %lf\"																							\n",dynamic_cast<t_gsetign*>(this)->lever()(0), dynamic_cast<t_gsetign*>(this)->lever()(1), dynamic_cast<t_gsetign*>(this)->lever()(2));
	ss << buffer; sprintf(buffer,"	max_pdop = \"3\"																									\n");
	ss << buffer; sprintf(buffer,"	min_sat = \"5\"																										\n");
	ss << buffer; sprintf(buffer,"	delay_t = \"0.0001\"																								\n");
	ss << buffer; sprintf(buffer,"	initial_misalignment_std = \"1,1,8\"																				\n");
	ss << buffer; sprintf(buffer,"	initial_vel_std = \"0.1,0.1,0.1\"																					\n");
	ss << buffer; sprintf(buffer,"	initial_pos_std = \"3,3,3\"																							\n");
	ss << buffer; sprintf(buffer,"	initial_gyro_std = \"10\"																							\n");
	ss << buffer; sprintf(buffer,"	initial_acce_std = \"1\"																							\n");
	ss << buffer; sprintf(buffer, "																														\n");
	ss << buffer; sprintf(buffer,"	misalignment_psd = \"%lf\"																							\n", dynamic_cast<t_gsetign*>(this)->misalignment_psd()(0));
	ss << buffer; sprintf(buffer,"	vel_psd = \"%lf\"																									\n", dynamic_cast<t_gsetign*>(this)->vel_psd()(0));
	ss << buffer; sprintf(buffer,"	pos_psd = \"%lf\"																									\n", dynamic_cast<t_gsetign*>(this)->pos_psd()(0));
	ss << buffer; sprintf(buffer,"	gyro_psd = \"%lf\"																									\n", dynamic_cast<t_gsetign*>(this)->gyro_psd()(0));
	ss << buffer; sprintf(buffer,"	acce_psd = \"%lf\"																									\n", dynamic_cast<t_gsetign*>(this)->acce_psd()(0));
	ss << buffer; sprintf(buffer,"	/>																													\n");
	ss << buffer; sprintf(buffer,"																														\n");
	ss << buffer; sprintf(buffer,"	<ambiguity>																											\n");
	ss << buffer; sprintf(buffer,"	<fix_mode> NO </fix_mode> <!--ROUND / SEARCH / NO-->																\n");
	ss << buffer; sprintf(buffer,"	<part_fix>  YES  </part_fix> <!--YES / NO-->																		\n");
	ss << buffer; sprintf(buffer,"	<ratio>  3.0  </ratio><!--Threshold in LAMBDA method-->																\n");
	ss << buffer; sprintf(buffer,"	<add_leo>  NO </add_leo>																							\n");
	ss << buffer; sprintf(buffer,"	<all_baselines>  NO  </all_baselines>																				\n");
	ss << buffer; sprintf(buffer,"	<min_common_time>  0 </min_common_time> <!--Minimum common time / seconds-->										\n");
	ss << buffer; sprintf(buffer,"	<baseline_length_limit>  3500 </baseline_length_limit> <!--Baseline length limit / km-->							\n");
	ss << buffer; sprintf(buffer,"	<widelane_interval> 30 </widelane_interval> <!--seconds-->															\n");
	ss << buffer; sprintf(buffer,"	<widelane_decision  maxdev = \"0.275\" maxsig = \"0.10\" alpha = \"1000\" /><!--deriation, sigma in WL - cycle-->	\n");
	ss << buffer; sprintf(buffer,"	<narrowlane_decision maxdev = \"0.375\" maxsig = \"0.10\" alpha = \"1000\" /><!--deriation, sigma in NL - cycle-->	\n");
	ss << buffer; sprintf(buffer,"	</ambiguity>																										\n");
	ss << buffer; sprintf(buffer,"																														\n");
	map< GSYS, map<FREQ_SEQ, GOBSBAND> > band_index;
	map< GSYS, map<GOBSBAND, FREQ_SEQ> > freq_index;

	band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(this)->band_index(gnut::GPS);
	band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(this)->band_index(gnut::GAL);
	band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(this)->band_index(gnut::BDS);
	freq_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(this)->freq_index(gnut::GPS);
	freq_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(this)->freq_index(gnut::GAL);
	freq_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(this)->freq_index(gnut::BDS);
	string tmp_freq = "", tmp_band = "";

	tmp_freq = "", tmp_band = "";
	for (auto f_b : band_index[GSYS::GPS]) { tmp_freq += to_string(f_b.first) + " "; tmp_band += to_string(f_b.second) + " "; }
	ss << buffer; sprintf(buffer,"	<gps sigma_L = \"%lf\" sigma_C = \"%lf\">																			\n", dynamic_cast<t_gsetgnss*>(this)->sigma_L(GSYS::GPS), dynamic_cast<t_gsetgnss*>(this)->sigma_C(GSYS::GPS));
	ss << buffer; sprintf(buffer, "	<freq> %s </freq>																									\n", tmp_freq.c_str());
	ss << buffer; sprintf(buffer,"	<band> %s </band>																									\n", tmp_band.c_str());
	ss << buffer; sprintf(buffer,"	</gps>																												\n");
	
	tmp_freq = "", tmp_band = "";
	for (auto f_b : band_index[GSYS::GAL]) { tmp_freq += to_string(f_b.first) + " "; tmp_band += to_string(f_b.second) + " "; }
	ss << buffer; sprintf(buffer,"	<gal sigma_L = \"%lf\" sigma_C = \"%lf\">																			\n", dynamic_cast<t_gsetgnss*>(this)->sigma_L(GSYS::GAL), dynamic_cast<t_gsetgnss*>(this)->sigma_C(GSYS::GAL));
	ss << buffer; sprintf(buffer, "	<freq> %s </freq>																									\n", tmp_freq.c_str());
	ss << buffer; sprintf(buffer,"	<band> %s </band>																									\n", tmp_band.c_str());
	ss << buffer; sprintf(buffer,"	</gal>																												\n");
	
	tmp_freq = "", tmp_band = "";
	for (auto f_b : band_index[GSYS::BDS]) { tmp_freq += to_string(f_b.first) + " "; tmp_band += to_string(f_b.second) + " "; }
	ss << buffer; sprintf(buffer,"	<bds sigma_L = \"%lf\" sigma_C = \"%lf\">																			\n", dynamic_cast<t_gsetgnss*>(this)->sigma_L(GSYS::BDS), dynamic_cast<t_gsetgnss*>(this)->sigma_C(GSYS::BDS));
	ss << buffer; sprintf(buffer, "	<freq> %s </freq>																									\n", tmp_freq.c_str());
	ss << buffer; sprintf(buffer,"	<band> %s </band>																									\n", tmp_band.c_str());
	ss << buffer; sprintf(buffer,"	</bds>																												\n");
	ss << buffer; sprintf(buffer,"																														\n");
	ss << buffer; sprintf(buffer,"	</config>																											\n");
	
	ss << buffer;
	return true;

}
