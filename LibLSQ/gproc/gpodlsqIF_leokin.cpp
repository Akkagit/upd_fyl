//yjqin change sth about kinematic
#include "gproc/gpodlsqIF_leokin.h"
#include "gset/gsetpar.h"
#include "gutils/ginfolog.h"
#include "gutils/gleo.h"


namespace great
{
	t_gpodlsqIF_leokin::t_gpodlsqIF_leokin(string site, t_gsetbase * set, t_gallproc * data, t_glog * log) :
		t_gppplsq(site, set, data, log),
		//_allkin(nullptr),
		_result_kin(nullptr),
		_sampling(0.0)
	{
		// est of crd must be KIN in leo_geometric esitmate
		_crd_est = KIN;

		_sampling = dynamic_cast<t_gsetgen*>(set)->sampling();
		

		_obs_limit = dynamic_cast<t_gsetout*>(set)->sp3_obslimit();


		_setOut();
		
	}
	t_gpodlsqIF_leokin::~t_gpodlsqIF_leokin()
	{
		/*if (_result_kin) {
			_result_kin->flush();
			_result_kin->close();
			delete _result_kin;
			_result_kin = nullptr;
		}*/

	}
	bool t_gpodlsqIF_leokin::GenerateProduct()
	{
		//if (_result_kin) {
		//	_extract_kin();
		//	//return true;
		//}
		//else {
		//	write_log_info(_glog, 0, "Don't have kin out file,Use kin_out_temp! Please Check out the xml file!!");
		//	return false;
		//}
		if (_result_kin)
		{
			_newsp3 = new t_gallprec(); _newsp3->overwrite(true);
			string tmp = dynamic_cast<t_gsetout*>(_gset)->outputs("sp3");
			if (tmp.empty())
			{
				tmp = "sp3_file";
			}
			_init_sp3(_newsp3);
			if (!_add_ppp_sp3(_newsp3))
			{
				return false;
			}
			t_gio* gout_sp3 = new t_gfile;  t_gcoder* gsp3_coder = new t_sp3(_gset, "", 4096);
			string sp3_path = tmp;
			//string ics_path = (dynamic_cast<t_gsetinp*>(_set)->inputs("ics"))[0];
			gout_sp3->path(sp3_path); gout_sp3->glog(_glog);
			gsp3_coder->clear(); gsp3_coder->glog(_glog); gsp3_coder->path(sp3_path);
			gsp3_coder->add_data("ID1", _newsp3);
			gout_sp3->coder(gsp3_coder);
			gout_sp3->run_write();
			delete gout_sp3; gout_sp3 = nullptr;
			delete gsp3_coder; gsp3_coder = nullptr;
			delete _newsp3; _newsp3 = nullptr;
		}

		// lsq recover recclk
		return _extract_clkfile(*_gall_rcv, t_gallprec::AR);

	}

	void t_gpodlsqIF_leokin::_setOut()
	{
		//t_gppplsqIF::_setOut();

		string kin_tmp;
		kin_tmp = dynamic_cast<t_gsetout*>(_gset)->outputs("sp3");
		if (kin_tmp.empty()) {
			write_log_info(_glog, 0, "Don't have kin out file,Use kin_out_temp! Please Check out the xml file!!");
			kin_tmp = "kin_out_temp";
		}

		substitute(kin_tmp, "$(date)", int2str(_beg_time.year()) + int2str(_beg_time.doy()), false);
		if (_result_kin) {
			if (_result_kin->is_open()) {
				_result_kin->close();
			}
			delete _result_kin;
		}
		_result_kin = new t_giof;
		_result_kin->tsys(t_gtime::GPS);
		_result_kin->mask(kin_tmp);
		_result_kin->append(dynamic_cast<t_gsetout*>(_gset)->append());


	}




	bool t_gpodlsqIF_leokin::_init_rec_crd_pars(t_glsq * lsq, bool isPPP)
	{

		int ipar_num = _lsq->_x_solve.parNumber();

		// È¡½»¼¯
		// Here change for leo sat
		set<string> rec_tmp = _gall_obj->objects_set(t_gdata::REC_LEO);
		set<string> rec_final;
		set_intersection(rec_tmp.begin(), rec_tmp.end(),
			_rec_list.begin(), _rec_list.end(),
			insert_iterator<set<string>>(rec_final, rec_final.begin()));

		_rec_list = rec_final;
		if (_rec_list.empty())
		{
			write_log_info(_glog, 0, "ERROR : The rec list is empty !");
			return false;
		}

		auto it_beg = _rec_list.begin();
		auto it_end = _rec_list.end();
		for (auto it = it_beg; it != it_end; ++it)
		{
			string site = *it;
			shared_ptr<t_gobj> _grec = _gall_obj->obj(site);

			// get the init value
			// @TODO : consider the crd and std is zeros
			t_gtriple rec_crd = _grec->crd(_beg_time);
			t_gtriple rec_std = dynamic_cast<t_gsetpar*>(_gset)->sigRecPosTriple(site);

			// init xyz 
			t_gpar par_x(site, t_gpar::CRD_X, ++ipar_num, "");
			t_gpar par_y(site, t_gpar::CRD_Y, ++ipar_num, "");
			t_gpar par_z(site, t_gpar::CRD_Z, ++ipar_num, "");

			// set the parameter tim
			par_x.value(rec_crd[0]);
			par_y.value(rec_crd[1]);
			par_z.value(rec_crd[2]);
			// set the inital value
			t_gtime endT = _end_time;
			if (_crd_est == KIN) {
				endT = _beg_time;
			}
			par_x.setTime(_beg_time, endT);
			par_y.setTime(_beg_time, endT);
			par_z.setTime(_beg_time, endT);
			// set the apriori  [crd read from rinex]
			par_x.apriori(rec_std[0]);
			par_y.apriori(rec_std[1]);
			par_z.apriori(rec_std[2]);
			//  add parameters
			lsq->add_parameter(par_x);
			lsq->add_parameter(par_y);
			lsq->add_parameter(par_z);
		}
		return true;
	}

	int t_gpodlsqIF_leokin::ProcessOneEpoch(const t_gtime& now, vector<t_gsatdata>* obs)
	{
		_crt_time = now;
		vector<t_gsatdata> crt_obs;
		if (obs) crt_obs = *obs;
		else crt_obs = _gall_obs->obs(_rec_list, now);//maybe bug

		if (crt_obs.size() == 0) {
			if (_glog && _glog->verb() >= 1) {
				_glog->comment(1, "gppplsq", _crt_rec + now.str_ymdhms(" no observation found at epoch: "));
			}
			return -1;
		}

		/*set glonass frequency number and correct dcb*/
		for (auto& iter : crt_obs)
		{
			if (iter.sys() == "R" && iter.channel() >= DEF_CHANNEL) iter.channel(_glofrq_num[iter.sat()]);
			iter.apply_bias(*_gall_bias);
		}

		/*set kinematic coordinate*/
		if (_crd_est == KIN)
		{
			if (_spplsq)
			{
				//_spplsq->processBatch(now, now);
				_kin_xyz[_crt_rec][now] = _spplsq->getCrd(now);
			}
			if (double_eq(_kin_xyz[_crt_rec][now][0] * _kin_xyz[_crt_rec][now][1] * _kin_xyz[_crt_rec][now][2], 0.0))
			{
				write_log_info(_glog, 0, "ERROR : epoch process failed" + now.str_ymdhms("Processing epoch: "));
				return -1;
			}
			if (_kin_xyz[_crt_rec][now][0]>1E10|| _kin_xyz[_crt_rec][now][1] > 1E10|| _kin_xyz[_crt_rec][now][2] > 1E10)
			{
				write_log_info(_glog, 0, "ERROR : epoch process failed" + now.str_ymdhms("Processing epoch: "));
				return -1;
			}
			if (_kin_xyz[_crt_rec][now][0]<-1E10 || _kin_xyz[_crt_rec][now][1] < -1E10 || _kin_xyz[_crt_rec][now][2] < -1E10)
			{
				write_log_info(_glog, 0, "ERROR : epoch process failed" + now.str_ymdhms("Processing epoch: "));
				return -1;
			}
			_quality_control.processOneEpoch(now, _crt_rec, _kin_xyz[_crt_rec][now], crt_obs);
		}
		else
		{
			_quality_control.processOneEpoch(now, _crt_rec, xyz_standard, crt_obs);
		}

		/*process one epoch*/
		bool epoch_valid = _processOneEpoch(now, crt_obs);
		if (!epoch_valid)
		{
			write_log_info(_glog, 0, "ERROR : epoch process failed" + now.str_ymdhms("Processing epoch: "));
			return -1;
		}

		if (_lsq->mode() == gnut::EPO && now < _end_time)
		{
			bool solve_valid = t_gppplsq::_SolveEpoch();
			if (!solve_valid)
			{
				write_log_info(_glog, 0, "ERROR : can not solve epoch, " + now.str_ymdhms("Processing epoch: "));
				return -1;
			}
			else if (_fix_mode == t_gsetamb::FIX_MODE::NO || !_ambfix->amb_fixed())
			{
				write_log_info(_glog, 0, "PPP FLOAT" + now.str_ymdhms("Processing epoch: "));
				return 0;
			}
		}
		return 1;
	}


	bool t_gpodlsqIF_leokin::_processOneRec(const t_gtime & crt_epoch, const std::string & crt_rec, std::vector<t_gsatdata>& crt_obs)
	{
		if (_kin_xyz[crt_rec].find(crt_epoch) == _kin_xyz[crt_rec].end() || _kin_xyz[crt_rec][crt_epoch].zero()) {
			write_log_info(_glog, 0, "The Init Crd is Zero! Skip!");
			return false;
		}

		// care the kin file
		t_gtriple xyz_leo = _kin_xyz[crt_rec][crt_epoch];
		t_gtriple xyz_leo_kin;
		if (_gall_nav->lagrange_pos(strs2strl(crt_rec), crt_epoch, xyz_leo_kin) > 0)
		{
			if ((xyz_leo_kin - xyz_leo).norm() < 500) {
				xyz_leo = xyz_leo_kin;
			}
			else {
				write_log_info(_glog, 0, "SPP And KIN Crd is So Diff too much! Skip!!!");
				return false;
			}
		}

		// use update par value
		int _idx, _idy, _idz;

		//string _name_short = strl2strs((_allkin->get_sat_name()).substr(0, 15));
		_idx = _lsq->_x_solve.getParam(_crt_rec, t_gpar::CRD_X, "", crt_epoch, crt_epoch);
		_idy = _lsq->_x_solve.getParam(_crt_rec, t_gpar::CRD_Y, "", crt_epoch, crt_epoch);
		_idz = _lsq->_x_solve.getParam(_crt_rec, t_gpar::CRD_Z, "", crt_epoch, crt_epoch);


		_lsq->_x_solve[_idx].value(xyz_leo[0]);
		_lsq->_x_solve[_idx].apriori(100.0);

		_lsq->_x_solve[_idy].value(xyz_leo[1]);
		_lsq->_x_solve[_idy].apriori(100.0);

		_lsq->_x_solve[_idz].value(xyz_leo[2]);
		_lsq->_x_solve[_idz].apriori(100.0);


		return t_gppplsq::_processOneRec(crt_epoch, crt_rec, crt_obs);
	}

	void t_gpodlsqIF_leokin::_extract_kin()
	{


		/// storage for kin
		map<t_gtime, map<string, double> > kin_map_res;
		map<string, double> kin_map;

		string temp;
		string time;

		/// get relative information from recover data
		auto& map_time_par = _gall_rcv->get_map_time_par();
		for (auto iter = map_time_par.begin(); iter != map_time_par.end(); iter++) {
			t_gtime epoch = iter->first;
			for (auto recover_par : iter->second){
				string identity = "";
				if (recover_par->par.parType == t_gpar::CRD_X) {
					identity = "PX";
				}
				else if (recover_par->par.parType == t_gpar::CRD_Y) {
					identity = "PY";
				}
				else if (recover_par->par.parType == t_gpar::CRD_Z) {
					identity = "PZ";
				}
				else {
					continue;
				}
				kin_map[identity] = recover_par->par.value() + recover_par->correct_value;
				kin_map_res[recover_par->par.beg] = kin_map;
			}
		}
		/// write header
		_write_kin_header();

		/// write data
		for (auto iter_time = kin_map_res.begin(); iter_time != kin_map_res.end(); iter_time++) {
			_write_kin_data(iter_time->first, iter_time->second);
		}

		_write_kin_end();

	}

	bool t_gpodlsqIF_leokin::_init_sp3(t_gallprec* sp3)
	{		
		
		
		sp3->set_beg(_beg_time);
		sp3->set_end(_end_time);
		sp3->set_data_type("V");

		vector<string> sats_used;
		
		sats_used.push_back(strs2strl(_crt_rec));
		int nsize = sats_used.size();
				
		string sat3 = "L" + int2str(1, 2);
		_prn3.clear();
		_prn3.push_back(sat3);
		sp3->add_interval(int(_sampling));
		sp3->add_interval(sat3,int(_sampling));
				
		for (int i = nsize; i < 102; i++)
		{
			_prn3.push_back("  0");
		}

		sp3->set_sat(sats_used);
		sp3->set_sat3(_prn3);
		sp3->set_sat_type("LEO");
		sp3->set_agency("GNUT");
		return true;

	}

	void t_gpodlsqIF_leokin::_add_spp_sp3(t_gallprec* sp3)
	{
		
		map<t_gtime, map<string, double> > kin_map_res;
		map<string, double> kin_map;
		map<string, double> vel_map;
		string temp;
		string time;
		t_gtime now;
		t_gtime end = _end_time;
		/// get relative information from recover data
		for (now = _beg_time; now <= _end_time; now = now + _sampling)
		{
			kin_map["PX"] = _spplsq->getCrd(now)[0];
			kin_map["PY"] = _spplsq->getCrd(now)[1];
			kin_map["PZ"] = _spplsq->getCrd(now)[2];
			kin_map_res[now] = kin_map;
		}

		for (now = _beg_time; now < _end_time; now = now + _sampling)
		{
			vel_map["VX"] = (kin_map_res[now + _sampling]["PX"] - kin_map_res[now]["PX"]) / _sampling;
			vel_map["VY"] = (kin_map_res[now + _sampling]["PY"] - kin_map_res[now]["PY"]) / _sampling;
			vel_map["VZ"] = (kin_map_res[now + _sampling]["PZ"] - kin_map_res[now]["PZ"]) / _sampling;
			kin_map_res[now]["VX"] = vel_map["VX"];
			kin_map_res[now]["VY"] = vel_map["VY"];
			kin_map_res[now]["VZ"] = vel_map["VZ"];
		}
		kin_map_res[end] = kin_map_res[end - _sampling];
		/// write data
		for (auto iter_time = kin_map_res.begin(); iter_time != kin_map_res.end(); iter_time++) {
			t_gtriple pos;
			t_gtriple dxyz(0, 0, 0);
			pos[0] = iter_time->second["PX"] / 1000.0;
			pos[1] = iter_time->second["PY"] / 1000.0;
			pos[2] = iter_time->second["PZ"] / 1000.0;
			double vel[4];
			vel[0] = iter_time->second["VX"] / 1000.0;
			vel[1] = iter_time->second["VY"] / 1000.0;
			vel[2] = iter_time->second["VZ"] / 1000.0;
			vel[3] = 0;
			double var[4] = { 0,0,0,0 };
			sp3->add_pos_vel(_prn3[0], iter_time->first, pos, 0.0, dxyz, 0.0, vel, var, _obs_num_spp[iter_time->first]);
		}


	}


	bool t_gpodlsqIF_leokin::_add_ppp_sp3(t_gallprec* sp3)
	{
		map<t_gtime, map<string, double> > kin_map_res;
		map<t_gtime, int> obs_data;
		map<string, double> kin_map;
		map<string, double> vel_map;
		string temp;
		string time;
		t_gtime now;
		t_gtime end;
		/// get relative information from recover data
		auto& map_time_par = _gall_rcv->get_map_time_par();
		auto& map_time_equ = _gall_rcv->get_map_time_equ();
		if (map_time_par.size() == 0)
		{
			return false;
		}
		for (auto iter = map_time_par.begin(); iter != map_time_par.end(); iter++) {
			t_gtime epoch = iter->first;
			end = epoch;
			for (auto recover_par : iter->second) {
				string identity = "";
				if (recover_par->par.parType == t_gpar::CRD_X) {
					identity = "PX";
				}
				else if (recover_par->par.parType == t_gpar::CRD_Y) {
					identity = "PY";
				}
				else if (recover_par->par.parType == t_gpar::CRD_Z) {
					identity = "PZ";
				}
				else {
					continue;
				}
				kin_map[identity] = recover_par->par.value() + recover_par->correct_value;
				kin_map_res[recover_par->par.beg] = kin_map;
			}
		}
		for (auto iter = map_time_equ.begin(); iter != map_time_equ.end(); iter++) {
			t_gtime epoch = iter->first;
			end = epoch;
			obs_data[epoch] = iter->second.size()/2;
		}
		for (now = map_time_par.begin()->first; now < end; now = now + _sampling)
		{
			vel_map["VX"] = (kin_map_res[now + _sampling]["PX"] - kin_map_res[now]["PX"]) / _sampling;
			vel_map["VY"] = (kin_map_res[now + _sampling]["PY"] - kin_map_res[now]["PY"]) / _sampling;
			vel_map["VZ"] = (kin_map_res[now + _sampling]["PZ"] - kin_map_res[now]["PZ"]) / _sampling;
			kin_map_res[now]["VX"] = vel_map["VX"];
			kin_map_res[now]["VY"] = vel_map["VY"];
			kin_map_res[now]["VZ"] = vel_map["VZ"];
		}
		kin_map_res[end] = kin_map_res[end - _sampling];
		/// write data
		for (auto iter_time = kin_map_res.begin(); iter_time != kin_map_res.end(); iter_time++) {
			t_gtime epo = iter_time->first;
			int num = obs_data[epo];
			if (num < _obs_limit)
			{
				continue;
			}
			t_gtriple pos;
			t_gtriple dxyz(0, 0, 0);
			pos[0] = iter_time->second["PX"] / 1000.0;
			pos[1] = iter_time->second["PY"] / 1000.0;
			pos[2] = iter_time->second["PZ"] / 1000.0;
			double vel[4];
			vel[0] = iter_time->second["VX"] / 1000.0;
			vel[1] = iter_time->second["VY"] / 1000.0;
			vel[2] = iter_time->second["VZ"] / 1000.0;
			vel[3] = 0;
			double var[4] = { 0,0,0,0 };
			sp3->add_pos_vel(_prn3[0], epo, pos, 0.0, dxyz, 0.0, vel, var, num);
		}
		return true;

	}

	/**
	* write kin file.
	* @param [in] epo     the current time
	* @param [in] kin_map the information stored for kin file.
	*/
	void t_gpodlsqIF_leokin::_write_kin_data(t_gtime epo, map<string, double> kin_map)
	{
		ostringstream os;

		int y, m, d, h, min, s;
		epo.ymd(y, m, d);
		epo.hms(h, min, s);
		os << "*  " << y;
		os << setw(3) << right << m;
		os << setw(3) << right << d;
		os << setw(3) << right << h;
		os << setw(3) << right << min;
		os << setw(12) << setiosflags(ios::fixed) << setprecision(8) << right << (double)s << endl;

		os << "PK01";
		os << setw(18) << setiosflags(ios::fixed) << setprecision(6) << right << kin_map["PX"] / 1000.0;
		os << setw(18) << setiosflags(ios::fixed) << setprecision(6) << right << kin_map["PY"] / 1000.0;
		os << setw(18) << setiosflags(ios::fixed) << setprecision(6) << right << kin_map["PZ"] / 1000.0;
		os << setw(18) << setiosflags(ios::fixed) << setprecision(6) << right << kin_map["T"] << endl;

		/// Print results
		if (_result_kin) {
			_result_kin->write(os.str().c_str(), os.str().size());
			_result_kin->flush();
		}
		else {
			_glog->comment(2, "t_gclklsq::_printout", "have no output file!");
		}

		return;
	}


	void t_gpodlsqIF_leokin::_write_kin_end()
	{
		ostringstream os;
		os << "EOF";

		if (_result_kin) {
			_result_kin->write(os.str().c_str(), os.str().size());
			_result_kin->flush();
		}
		else {
			_glog->comment(2, "t_gclklsq::_write_clk_header", "have no output file!");
		}
		return;

	}

	/**
	* write kin file header.
	*/
	void t_gpodlsqIF_leokin::_write_kin_header()
	{
		ostringstream os;
		int  minute, hour;
		int _nepochs;
		int second;
		int tmp = 17;
		bool first_line = true;
		int leo_num = 1;

		string name;

		name = *_rec_list.begin();
		name = strs2strl(name);
		_beg_time.hms(hour, minute, second);
		_nepochs = (int)((((_end_time.dmjd() - _beg_time.dmjd())*86400.0) / _sampling) + 1);
		///write kin head
		if (!_result_kin->is_open())
		{
			os << "#aP" << _beg_time.year() << setw(3) << right << _beg_time.mon();
			os << setw(3) << right << _beg_time.day();
			os << setw(3) << right << hour;
			os << setw(3) << right << minute;
			os << setw(12) << setiosflags(ios::fixed) << setprecision(8) << right << ((double)second);
			os << setw(8) << right << _nepochs;
			os << " ORBIT" << " IGS14" << "  FIT" << " KIN" << " GPST" << endl;

			os << "## " << _beg_time.gwk();
			os << setw(16) << setiosflags(ios::fixed) << setprecision(8) << (double)_beg_time.sow();
			os << setw(14) << setiosflags(ios::fixed) << setprecision(8) << (double)_sampling;
			os << "  " << _beg_time.mjd() << " ";
			os << setw(15) << setiosflags(ios::fixed) << setprecision(13) << (_beg_time.sod() / 86400.0) << endl;
			os << "+    " << leo_num << "   " << "K01";
			for (int j = 0; j < 16; j++)
			{
				os << "  0";
			}
			os << endl;

			for (int k = 0; k < 4; k++)
			{
				os << "+        ";
				for (int i = 0; i < tmp; i++)
				{
					os << "  0";
				}
				os << endl;
				//os << "+        ";
			}

			for (int k = 0; k < 5; k++)
			{
				os << "++       ";
				for (int i = 0; i < tmp; i++)
				{
					if (first_line&&k == 0)
					{
						os << "  5";
						first_line = false;
						continue;
					}
					os << "  0";
				}
				os << endl;
				//os << "+        ";
			}

			os << "%c L " << name << endl;
			os << "/*generated by GNUT" << endl;
		}
		if (_result_kin) {
			_result_kin->write(os.str().c_str(), os.str().size());
			_result_kin->flush();
		}
		else {
			_glog->comment(2, "t_gclklsq::_write_clk_header", "have no output file!");
		}
		return;
	}



	bool t_gpodlsqIF_leokin::_initLsqProcData(t_gallproc * data)
	{
		gtrace("t_gppplsq::_initLsqProcData");

		// prepare obs data
		_gall_bias = dynamic_cast<t_gallbias*>((*data)[t_gdata::ALLBIAS]);
		_gall_obs = dynamic_cast<t_gallobs*>((*data)[t_gdata::ALLOBS]);
		_gall_nav = dynamic_cast<t_gallnav*>((*data)[t_gdata::GRP_EPHEM]);
		_gall_obj = dynamic_cast<t_gallobj*>((*data)[t_gdata::ALLOBJ]);
		_gall_prec    = dynamic_cast<t_gallprec*>((*data)[t_gdata::ALLPREC]);
		
		if (!_read_kin())
		{
			if (_crd_est == KIN)
			{
				_spplsq = new t_gspplsq(_crt_rec, _gset, _gall_obj);
				_spplsq->setDAT(_gall_obs, _gall_nav);
				_spplsq->processBatch(_beg_time, _end_time);
				_obs_num_spp = _spplsq->getobsnum();
				_generate_spp_kin();
			}
		}
		if (_fix_mode == t_gsetamb::FIX_MODE::SEARCH || _fix_mode == t_gsetamb::FIX_MODE::ROUND)
		{
			_gupd = dynamic_cast<t_gupd*>((*data)[t_gdata::UPD]);
			//if (_frequency == 3)
			//{
			//	_gupd_epoch = dynamic_cast<t_gupd_epoch*>((*data)[t_gdata::UPD_EPOCH]);
			//}
		}

		// Whether data is valid
		if (!_gall_bias || !_gall_obs || !_gall_obj || (!_gall_nav && _crd_est == KIN) ||
			(_isClient && _fix_mode != t_gsetamb::FIX_MODE::NO && !_gupd))
		{
			write_log_info(_glog, 0, "ERROR £ºcan not init lsq process data");
			return false;
		}

		// init SPP for kinematic station
		

		// set glofrq number 
		_glofrq_num = _gall_obs->glo_freq_num();
		if (_glofrq_num.empty() && _gall_nav)
		{
			_glofrq_num = _gall_nav->glo_freq_num();
		}

		//// special for simulation PPP
		//if (_simulation)
		//{
		//	raw_CODE_bias(_beg_time, _allbias, _glofrq_num, _simu_raw_dcb);
		//}

		write_log_info(_glog, 4, "INIT : process data in lsq is prepared.");
		return true;
	}

	bool t_gpodlsqIF_leokin::_generate_spp_kin()
	{
		if (_result_kin)
		{
			_sppsp3 = new t_gallprec(); _sppsp3->overwrite(true);
			string tmp  = "sp3_file";
			_init_sp3(_sppsp3);
			_add_spp_sp3(_sppsp3);
			t_gio* gout_sp3 = new t_gfile;  t_gcoder* gsp3_coder = new t_sp3(_gset, "", 4096);
			string sp3_path = tmp;
			//string ics_path = (dynamic_cast<t_gsetinp*>(_set)->inputs("ics"))[0];
			gout_sp3->path(sp3_path); gout_sp3->glog(_glog);
			gsp3_coder->clear(); gsp3_coder->glog(_glog); gsp3_coder->path(sp3_path);
			gsp3_coder->add_data("ID1", _sppsp3);
			gout_sp3->coder(gsp3_coder);
			gout_sp3->run_write();
			delete gout_sp3; gout_sp3 = nullptr;
			delete gsp3_coder; gsp3_coder = nullptr;
			delete _sppsp3; _sppsp3 = nullptr;
		}

		// lsq recover recclk
		// return _extract_clkfile(*_allrcv, t_gallprec::AR);
		return true;
	}
	bool t_gpodlsqIF_leokin::_read_kin()
	{
		set<string> sat_sp3 = _gall_prec->satellites();
		string leo = strs2strl(_crt_rec);
		if (sat_sp3.count(leo) == 0)
		{
			write_log_info(_glog, 0, "read_kinfile : cannot find kin of leo.");
			return false;
		}
		t_gtime beg, end, epo;
		/*double pos[3],vel[3];
		int obs_num = 0;
		beg = _gall_prec->beg_data(leo);
		end = _gall_prec->end_data(leo);
		for (epo = beg; epo <= end; epo = epo + _sampling)
		{
			if (_gall_prec->get_pos_vel(leo, epo, pos, vel, obs_num))
			{
				_kin_xyz[_crt_rec][epo] = pos;
			}
			
		}*/
		t_gtriple pos, vel;
		beg = _gall_prec->beg_data(leo);
		end = _gall_prec->end_data(leo);
		for (epo = beg; epo <= end; epo = epo + _sampling)
		{
			_gall_prec->lagrange_pos(leo, epo, pos);
			_kin_xyz[_crt_rec][epo] = pos;
		}
		return true;
	}
}