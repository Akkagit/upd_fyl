#include "imufile.h"
#include "gdata/gimudata.h"
#include "gset/gsetins.h"
using namespace Eigen;

gsins::t_imufile::t_imufile(t_gsetbase* s, string version, int sz) :t_gcoder(s, version, sz),
_ts(0.005),
_GyroUnit(DPS),
_AcceUnit(MPS2),
_complete(false),
_order("garfu")
{
	gtrace("t_imufile::constructor");
	_ts = dynamic_cast<t_gsetins*>(s)->ts();
	_freq = dynamic_cast<t_gsetins*>(s)->freq();
	_GyroUnit = dynamic_cast<t_gsetins*>(s)->GyroUnit();
	_AcceUnit = dynamic_cast<t_gsetins*>(s)->AcceUnit();
	_MagUnit = dynamic_cast<t_gsetins*>(s)->MagUnit(); //add wh
	_order = dynamic_cast<t_gsetins*>(s)->order();
}

int gsins::t_imufile::decode_head(char* buff, int sz, vector<string>& errmsg)
{
#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_mutex.lock();

	// no header expected, but fill the buffer
	t_gcoder::_add2buffer(buff, sz);
	_mutex.unlock(); return -1;
}

int gsins::t_imufile::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
{
	gtrace("t_imufile::decode_data");
	_mutex.lock();

	if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

	int tmpsize = 0;

	bool isbinary = false;
	if (_order == "starneto" || _order == "imr") isbinary = true;

	if (!isbinary)
	{
		string line;
		while (!_complete && (tmpsize = t_gcoder::_getline(line, 0)) >= 0)
		{
			if (_order == "unicore" && line[0] != '#')
			{
				t_gcoder::_consume(tmpsize);
				continue;
			}
			if ((line[0] == '%'&&_order != "imr_fsas") || (line[0] == '#' && _order != "unicore") || line[0] == 'T')
			{
				t_gcoder::_consume(tmpsize);
				continue;
			}
			for (int i = 0; i < line.size(); i++)
			{
				if (line[i] == ',' || line[i] == '*' || line[i] == ';')line[i] = ' ';
			}
			double t, gx, gy, gz, ax, ay, az, mx, my, mz;
			Vector3d g_tmp, a_tmp, m_tmp;
			Vector3d wtmp, vtmp, mtmp;
			stringstream ss(line);
			if (_order[0] == 'g' && _order[1] == 'a') // g->a
			{
				ss >> t >> gx >> gy >> gz >> ax >> ay >> az;
				g_tmp = Vector3d(gx, gy, gz); a_tmp = Vector3d(ax, ay, az);
				///added wh
				if (_order[2] == 'm') // for mag data
				{
					ss >> mx >> my >> mz;
					m_tmp = Vector3d(mx, my, mz);
				}
			}
			else if (_order[0] == 'a' && _order[1] == 'g') // a->g
			{
				ss >> t >> ax >> ay >> az >> gx >> gy >> gz;
				g_tmp = Vector3d(gx, gy, gz); a_tmp = Vector3d(ax, ay, az);
				///added wh
				if (_order[2] == 'm') // for mag data
				{
					ss >> mx >> my >> mz;
					m_tmp = Vector3d(mx, my, mz);
				}
			}
			else if (_order == "unicore")
			{
				vector<double> v;
				if (this->decode_unicore(line, t, v) > 0 && line.size() < 140)  //usually one line's size<140 --sbs 1018
				{
					az = v[0]; ay = -v[1]; ax = v[2];
					gz = v[3]; gy = -v[4]; gx = v[5];
					g_tmp = Vector3d(gx, gy, gz); a_tmp = Vector3d(ax, ay, az);
					wtmp = g_tmp; vtmp = a_tmp;  //--sbs 1007
				}
				else
				{
					if (_log)  _log->comment(0, "decode unicore data failed!");
					t_gcoder::_consume(tmpsize);
					continue;
				}
			}
			else if (_order == "imr_fsas")
			{
				vector<double> v;
				if (this->decode_imr_fsas(line, t, v) > 0 && line.size() < 120)
				{
					az = v[0]; ay = -v[1]; ax = v[2];
					gz = v[3]; gy = -v[4]; gx = v[5];
					g_tmp = Vector3d(gx, gy, gz); a_tmp = Vector3d(ax, ay, az);
					wtmp = g_tmp*glv.sec; vtmp = a_tmp;
				}
				else
				{
					if (_log)  _log->comment(0, "decode unicore data failed!");
					t_gcoder::_consume(tmpsize);
					continue;
				}
			}
			else
			{
				if (_log)  _log->comment(0, "The imu data format is not exist!");
			}

			///<-- modified by wh  substr(2,3) -> substr(_order.size()-3,3) ;   for read mag data
			if (_order.substr(_order.size() - 3, 3) == "rfu")  // right -> forward -> up
			{
				wtmp = g_tmp; vtmp = a_tmp;

				if (_order[2] == 'm')
				{
					mtmp = m_tmp;
				}
			}
			else if (_order.substr(_order.size() - 3, 3) == "flu") // forward -> left -> up
			{
				Matrix3d Rotation; Rotation << 0, -1, 0, 1, 0, 0, 0, 0, 1;
				wtmp = Rotation * g_tmp; vtmp = Rotation * a_tmp;
				if (_order[2] == 'm')
				{
					mtmp = Rotation * m_tmp;
				}
			}
			else if (_order.substr(_order.size() - 3, 3) == "frd") // forward -> right -> down
			{
				Matrix3d Rotation; Rotation << 0, 1, 0, 1, 0, 0, 0, 0, -1;
				wtmp = Rotation * g_tmp; vtmp = Rotation * a_tmp;
				if (_order[2] == 'm')
				{
					mtmp = Rotation * m_tmp;
				}
			}
			else if (_order.substr(_order.size() - 3, 3) == "rbd") // right -> behind -> down
			{
				Matrix3d Rotation; Rotation << 1, 0, 0, 0, -1, 0, 0, 0, -1;
				wtmp = Rotation * g_tmp; vtmp = Rotation * a_tmp;
				if (_order[2] == 'm')
				{
					mtmp = Rotation * m_tmp;
				}
			}

			switch (_GyroUnit)
			{
			case RAD:
				break;
			case DEG:
				wtmp *= glv.deg; break;
			case RPS:
				wtmp *= _ts; break;
			case DPS:
				wtmp *= glv.dps * _ts; break;
			case RPH:
				wtmp /= glv.hur * _ts; break;
			case DPH:
				wtmp *= glv.dph * _ts; break;
			default:
				break;
			}
			switch (_AcceUnit)
			{
			case MPS:
				break;
			case MPS2:
				vtmp *= _ts; break;
			}

			///add wh
			switch (_MagUnit)
			{
			case MPS:
				break;
			case MPS2:
				mtmp *= _ts; break;
			}


			if (t < dynamic_cast<t_gsetins*>(_set)->start()) //-1反向可输出最后一个历元
			{
				t_gcoder::_consume(tmpsize);
				continue;
			}
			if (t > dynamic_cast<t_gsetins*>(_set)->end()) //+1正向可输出最后一个历元
				_complete = true;

			_tt = t;

			map<string, t_gdata*>::iterator it = _data.begin();
			while (it != _data.end())
			{
				((t_gimudata*)it->second)->set_ts(_ts);
				if (it->second->id_type() == t_gdata::IMUDATA)
				{
					///modified by wh  read mag data
					if (_order[2] == 'm')
						((t_gimudata*)it->second)->add_IMU(t, wtmp, vtmp, mtmp);
					else
						((t_gimudata*)it->second)->add_IMU(t, wtmp, vtmp);
				}
				it++;
			}
			if (ss.fail()) {
				if (_log) _log->comment(1, "imufile", "warning: incorrect IMU data record: " + ss.str());
				else       cerr << "warning: incorrect IMU data record: " << ss.str() << endl;
				t_gcoder::_consume(tmpsize);
				_mutex.unlock(); return -1;
			}
			t_gcoder::_consume(tmpsize);
			cnt++;
		}
	}
	else
	{
		const char* block;

		// decode all the data in the buffer in one time (instead of using a loop)
		if (!_complete && (tmpsize = t_gcoder::_getbuffer(block)) >= 0)
		{
			vector<pair<double, vector<double>>> v; // imu data should be RAD, MPS
			int consume_size = -1;

			if (_order == "starneto") consume_size = decode_starneto(block, tmpsize, v);
			else if (_order == "imr") consume_size = -1; /* not implemented yet */

			if (consume_size > 0)
			{
				tmpsize = consume_size;
			}
			else
			{
				if (_log)  _log->comment(0, "decode starneto data failed!");
				t_gcoder::_consume(tmpsize);
				_mutex.unlock();
				return -1;
			};

			// add data records to t_gdata
			if (v.size() > 0 && v.back().first >= dynamic_cast<t_gsetins*>(_set)->start())
			{
				if (v.front().first > dynamic_cast<t_gsetins*>(_set)->end() + 5)
					_complete = true;

				_tt = v.back().first;

				map<string, t_gdata*>::iterator it = _data.begin();
				while (it != _data.end())
				{
					if (it->second->id_type() == t_gdata::IMUDATA)
					{
						((t_gimudata*)it->second)->set_ts(_ts);
						for (int i = 0; i < v.size(); i++)
						{
							Vector3d wtmp, vtmp;
							for (int j = 0; j < 3; j++) wtmp(j) = v[i].second[j];
							for (int j = 0; j < 3; j++) vtmp(j) = v[i].second[j + 3];
							((t_gimudata*)it->second)->add_IMU(v[i].first, wtmp, vtmp);
						}
					}
					it++;
				}
			}
			t_gcoder::_consume(tmpsize);
			cnt++;
		}
	}
	_mutex.unlock();
	return 0;
}

int gsins::t_imufile::decode_unicore(const string& line, double& t, vector<double>& v)
{
	gtrace("t_imufile::decode_unicore");
	try {
		vector<string> ret;
		split(line, " ", ret);
		if (ret.size() != 22)  return -1 ;  // incomplete  --sbs

		int week = stoi(ret[5]);
		t = stod(ret[6]);
		for (int i = 15; i <= 17; i++)
		{
			//cout << stod(ret[i]) << " " << ACC_SCALE << endl;
			v.push_back(stod(ret[i]) * UNICORE_ACC_SCALE);
		}
		for (int i = 18; i <= 20; i++)
			v.push_back(stod(ret[i]) * UNICORE_GYRO_SCALE);  //--sbs 1007
		return 1;
	}
	catch (...)
	{
		return -1;
	}
}

int gsins::t_imufile::decode_imr_fsas(const string & line, double & t, vector<double>& v)
{
	gtrace("t_imufile::decode_unicore");
	try {
		vector<string> ret;
		split(line, " ", ret);
		if (ret.size() != 13)  return -1;  // incomplete  --sbs

		int week = stoi(ret[3]);
		t = stod(ret[4]);
		for (int i = 6; i <= 8; i++)
		{
			//cout << stod(ret[i]) << " " << ACC_SCALE << endl;
			v.push_back(stod(ret[i]) * IMR_FSAS_ACCE_SCALE);
		}
		for (int i = 9; i <= 11; i++)
			v.push_back(stod(ret[i]) * IMR_FSAS_GYRO_SCALE);  //--sbs 1007
		return 1;
	}
	catch (...)
	{
		return -1;
	}
}

int gsins::t_imufile::decode_starneto(const char* block, int sz, vector<pair<double, vector<double>>>& v)
{
	try {
		const int MARGIN = 50;
		int consume_size = 0;

		for (int i = 0; i < sz - MARGIN; i++)
		{
			if (block[i] == (char)0xAA && block[i + 1] == (char)0x55 && block[i + 2] == (char)0x05)
			{
				i += 3;
				int byten = 0;
				uint16_t    week;
				uint32_t    sow;
				double      gx;
				double      gy;
				double      gz;
				double      ax;
				double      ay;
				double      az;

				byten = sizeof(week); memcpy(&week, block + i, byten);  i += byten;
				byten = sizeof(sow);  memcpy(&sow, block + i, byten);  i += byten;
				byten = sizeof(gx);   memcpy(&gx, block + i, byten);  i += byten;
				byten = sizeof(gy);   memcpy(&gy, block + i, byten);  i += byten;
				byten = sizeof(gz);   memcpy(&gz, block + i, byten);  i += byten;
				byten = sizeof(ax);   memcpy(&ax, block + i, byten);  i += byten;
				byten = sizeof(ay);   memcpy(&ay, block + i, byten);  i += byten;
				byten = sizeof(az);   memcpy(&az, block + i, byten);  i += byten;

				gx *= glv.deg * _ts;
				gy *= glv.deg * _ts;
				gz *= glv.deg * _ts;
				ax *= STARNETO_G * _ts;
				ay *= STARNETO_G * _ts;
				az *= STARNETO_G * _ts;

				v.push_back(make_pair(sow / 1000.0, vector<double>({ gx,gy,gz,ax,ay,az })));
			}
			else if (block[i] == (char)0xAA && block[i + 1] == (char)0x44 && block[i + 2] == (char)0x13)
			{
				int byten = 0;
				i += 6;

				uint16_t    week;
				uint32_t    sow;
				uint32_t    week2;
				double      sow2;
				int32_t     state;
				int32_t     az;
				int32_t     ay;
				int32_t     ax;
				int32_t     gz;
				int32_t     gy;
				int32_t     gx;
				byten = sizeof(week);  memcpy(&week, block + i, byten); i += byten;
				byten = sizeof(sow);   memcpy(&sow, block + i, byten); i += byten;
				byten = sizeof(week2); memcpy(&week2, block + i, byten); i += byten;
				byten = sizeof(sow2);  memcpy(&sow2, block + i, byten); i += byten;
				byten = sizeof(state); memcpy(&state, block + i, byten); i += byten;
				byten = sizeof(az);    memcpy(&az, block + i, byten); i += byten;
				byten = sizeof(ay);    memcpy(&ay, block + i, byten); i += byten;
				byten = sizeof(ax);    memcpy(&ax, block + i, byten); i += byten;
				byten = sizeof(gz);    memcpy(&gz, block + i, byten); i += byten;
				byten = sizeof(gy);    memcpy(&gy, block + i, byten); i += byten;
				byten = sizeof(gx);    memcpy(&gx, block + i, byten); i += byten;

				double gyro[3] = { gx,-gy,gz };
				double accel[3] = { ax,-ay,az };

				for (int i = 0; i < 3; i++)  gyro[i] *= STARNETO_GYRO_SCALE;
				for (int i = 0; i < 3; i++) accel[i] *= STARNETO_ACC_SCALE;
				v.push_back(make_pair(sow / 1000.0, vector<double>({ gyro[0],gyro[1],gyro[2],accel[0],accel[1],accel[2] })));

			}
			consume_size = i;
		}

		return consume_size;
	}
	catch (...)
	{
		return -1;
	}
}

bool gsins::t_imufile::available(const t_gtime& now)
{
	double t = now.sow() + now.dsec();
	if (_tt > t)return true;
	else
	{
		return false;
	}
}

int gsins::t_imufile::encode_head(char* buff, int sz, vector<string>& errmsg)
{
#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_mutex.lock();

	_mutex.unlock(); return -1;
}
int gsins::t_imufile::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
{
#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_mutex.lock();
	try
	{
		if (_ss_position == 0)
		{
			t_gimudata* imudata = NULL;
			for (auto it = _data.begin(); it != _data.end(); ++it)
			{
				if (it->second->id_type() == t_gdata::IMUDATA)
				{
					imudata = dynamic_cast<t_gimudata*>(it->second);
				}
			}
			if (imudata)
			{
				while (true)
				{
					std::vector<Eigen::Vector3d> wm;
					std::vector<Eigen::Vector3d> vm;
					double t = -1.0, ts = -1.0;
					imudata->load(wm, vm, t, ts, 1); // rfu
					if (t < 0)
						break;
					Eigen::Vector3d g_tmp = wm[0]; // RAD
					Eigen::Vector3d a_tmp = vm[0]; // MPS

					Eigen::Vector3d wtmp, vtmp;

					if (_order.substr(2, 3) == "rfu")  // right -> forward -> up
					{
						wtmp = g_tmp; vtmp = a_tmp;
					}
					else if (_order.substr(2, 3) == "flu") // forward -> left -> up
					{
						Matrix3d Rotation; Rotation << 0, -1, 0, 1, 0, 0, 0, 0, 1;
						wtmp = Rotation.transpose() * g_tmp; vtmp = Rotation.transpose() * a_tmp;
					}
					else if (_order.substr(2, 3) == "frd") // forward -> right -> down
					{
						Matrix3d Rotation; Rotation << 0, 1, 0, 1, 0, 0, 0, 0, -1;
						wtmp = Rotation.transpose() * g_tmp; vtmp = Rotation.transpose() * a_tmp;
					}
					else if (_order.substr(2, 3) == "rbd") // right -> behind -> down
					{
						Matrix3d Rotation; Rotation << 1, 0, 0, 0, -1, 0, 0, 0, -1;
						wtmp = Rotation.transpose() * g_tmp; vtmp = Rotation.transpose() * a_tmp;
					}

					switch (_GyroUnit)
					{
					case RAD:
						break;
					case DEG:
						wtmp /= glv.deg; break;
					case RPS:
						wtmp /= _ts; break;
					case DPS:
						wtmp /= glv.dps * _ts; break;
					case RPH:
						wtmp *= glv.hur * _ts; break;
					case DPH:
						wtmp /= glv.dph * _ts; break;
					default:
						break;
					}
					switch (_AcceUnit)
					{
					case MPS:
						break;
					case MPS2:
						vtmp /= _ts; break;
					}

					if (_order[0] == 'g' && _order[1] == 'a') // g->a
					{
						_ss << setw(15) << setiosflags(ios::fixed) << setprecision(4) << t << " " << resetiosflags(ios::fixed)
							<< setw(20) << setprecision(10) << wtmp(0) << " "
							<< setw(20) << setprecision(10) << wtmp(1) << " "
							<< setw(20) << setprecision(10) << wtmp(2) << " "
							<< setw(20) << setprecision(10) << vtmp(0) << " "
							<< setw(20) << setprecision(10) << vtmp(1) << " "
							<< setw(20) << setprecision(10) << vtmp(2) << " " << std::endl;
						std::cout << setw(15) << setiosflags(ios::fixed) << setprecision(4) << t << " " << resetiosflags(ios::fixed)
							<< setw(20) << setprecision(10) << wtmp(0) << " "
							<< setw(20) << setprecision(10) << wtmp(1) << " "
							<< setw(20) << setprecision(10) << wtmp(2) << " "
							<< setw(20) << setprecision(10) << vtmp(0) << " "
							<< setw(20) << setprecision(10) << vtmp(1) << " "
							<< setw(20) << setprecision(10) << vtmp(2) << " " << std::endl;
					}
					else if (_order[0] == 'a' && _order[1] == 'g') // a->g
					{
						_ss << setw(15) << setiosflags(ios::fixed) << setprecision(4) << t << " " << resetiosflags(ios::fixed)
							<< setw(20) << setprecision(10) << vtmp(0) << " "
							<< setw(20) << setprecision(10) << vtmp(1) << " "
							<< setw(20) << setprecision(10) << vtmp(2) << " "
							<< setw(20) << setprecision(10) << wtmp(0) << " "
							<< setw(20) << setprecision(10) << wtmp(1) << " "
							<< setw(20) << setprecision(10) << wtmp(2) << " " << std::endl;
					}
					else
					{
						if (_log)
						{
							_log->comment(0, "ERROR : t_imufile::encode_data order bot supported.");
						}
						else
						{
							cout << "ERROR : t_imufile::encode_data order bot supported." << endl;
						}
					}
				}

			}
		}
	}
	catch (...)
	{
		if (_log)
		{
			_log->comment(0, "ERROR : t_imufile::encode_data throw exception");
		}
		else
		{
			cout << "ERROR : t_imufile::encode_data throw exception" << endl;
		}
		_mutex.unlock(); return -1;
	}

	int size = _fill_buffer(buff, sz);
	_mutex.unlock();
	return size;
}



