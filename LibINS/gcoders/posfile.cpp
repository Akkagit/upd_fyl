#include "posfile.h"
#include "gins/gbase.h"
#include "gset/gsetign.h"
using namespace Eigen;
gsins::t_posfile::t_posfile(t_gsetbase * s, string version, int sz) : t_gcoder(s, version, sz)
{
	gtrace("t_posfile::constructor");
	if (!s) return;
	TS = dynamic_cast<t_gsetign*>(s)->TS();
	_order = dynamic_cast<t_gsetign*>(s)->order();
}

gsins::t_posfile::~t_posfile()
{
}

int gsins::t_posfile::decode_head(char * buff, int sz, vector<string>& errmsg)
{
#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_mutex.lock();

	// no header expected, but fill the buffer
	t_gcoder::_add2buffer(buff, sz);
	_mutex.unlock(); return -1;
}

int gsins::t_posfile::decode_data(char * buff, int sz, int & cnt, vector<string>& errmsg)
{
	gtrace("t_posfile::decode_data");
	_mutex.lock();

	if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

	int tmpsize = 0;
	string line;
	while ((tmpsize = t_gcoder::_getline(line, 0)) >= 0)
	{
		if (line[0] == '%' || line[0] == '#')
		{
			t_gcoder::_consume(tmpsize);
			continue;
		}
		for (int i = 0; i < line.size(); i++)
		{
			if (line[i] == ',')line[i] = ' ';
		}
		int nSat;
		double t, XYZ[3], Vxyz[3], Att[3], RXYZ[3], RVxyz[3], PDOP, ratio;
		bool state;
		Vector3d ell, v, vn;
		Vector3d R_ell, R_v, R_vn;
		stringstream ss(line);
		if (_order == "XYZXYZ")
		{
			ss >> t >> XYZ[0] >> XYZ[1] >> XYZ[2] >> Vxyz[0] >> Vxyz[1] >> Vxyz[2];
			//ss >> QXYZ[0] >> QXYZ[1] >> QXYZ[2] >> QVxyz[0] >> QVxyz[1] >> QVxyz[2];
			ell = Cart2Geod(Vector3d(XYZ), false); v = Vector3d(Vxyz);
			vn = gsins::t_gbase::Cen(ell).transpose()*v;
		}
		else if (_order == "XYZ")
		{
			ss >> t >> XYZ[0] >> XYZ[1] >> XYZ[2];
			ell = Cart2Geod(Vector3d(XYZ), false);
			vn = Eigen::Vector3d::Zero();
		}
		else if (_order=="BLHENU")
		{
			ss >> t >> XYZ[0] >> XYZ[1] >> XYZ[2] >> Vxyz[0] >> Vxyz[1] >> Vxyz[2];
			ell = Vector3d(XYZ[0] * glv.deg, XYZ[1] * glv.deg, XYZ[2]);
			vn = Vector3d(Vxyz);
		}
		else if (_order == "LBHENU")
		{
			ss >> t >> XYZ[1] >> XYZ[0] >> XYZ[2] >> Vxyz[0] >> Vxyz[1] >> Vxyz[2];
			ell = Vector3d(XYZ[0] * glv.deg, XYZ[1] * glv.deg, XYZ[2]);
			vn = Vector3d(Vxyz);
		}
		else if (_order == "LBHNEU")
		{
			ss >> t >> XYZ[1] >> XYZ[0] >> XYZ[2] >> Vxyz[1] >> Vxyz[0] >> Vxyz[2];
			ell = Vector3d(XYZ[0] * glv.deg, XYZ[1] * glv.deg, XYZ[2]);
			vn = Vector3d(Vxyz);
		}
		else if (_order == "LBHRPYENU")
		{
			ss >> t >> XYZ[1] >> XYZ[0] >> XYZ[2] >> Att[0] >> Att[1] >> Att[2] >> Vxyz[0] >> Vxyz[1] >> Vxyz[2];
			ell = Vector3d(XYZ[0] * glv.deg, XYZ[1] * glv.deg, XYZ[2]);
			vn = Vector3d(Vxyz);
		}
		else if (_order == "PVT_POS")
		{
			double tmp; string s;
			ss >> t >> XYZ[0] >> XYZ[1] >> XYZ[2] >> Vxyz[0] >> Vxyz[1] >> Vxyz[2];
			ss >> RXYZ[0] >> RXYZ[1] >> RXYZ[2] >> RVxyz[0] >> RVxyz[1] >> RVxyz[2];
			ss >> nSat >> PDOP >> tmp >> s;
			if (s == "Fixed")state = true;
			else state = false;

			ell = Cart2Geod(Vector3d(XYZ), false); v = Vector3d(Vxyz);
			vn = gsins::t_gbase::Cen(ell).transpose()*v;
			
			Eigen::Matrix3d R_XYZ = Eigen::Vector3d(RXYZ[0], RXYZ[1], RXYZ[2]).array().matrix().asDiagonal();
			Eigen::Matrix3d A = t_gbase::dGeod2Cart(gsins::t_gearth(), ell);
			Eigen::Matrix3d R_ELL = A.inverse()*R_XYZ*(A.inverse()).transpose();
			R_ell = Eigen::Vector3d(R_ELL(0, 0), R_ELL(1, 1), R_ELL(2, 2));
			
			Eigen::Matrix3d R_VXYZ = Eigen::Vector3d(RVxyz[0], RVxyz[1], RVxyz[2]).array().matrix().asDiagonal();
			Eigen::Matrix3d C = t_gbase::Cen(ell).transpose();
			Eigen::Matrix3d R = C*R_VXYZ*C.transpose();
			R_vn= Eigen::Vector3d(R(0, 0), R(1, 1), R(2, 2));

			ell = Vector3d(XYZ); vn = Vector3d(Vxyz);
			R_ell = Eigen::Vector3d(R_XYZ(0, 0), R_XYZ(1, 1), R_XYZ(2, 2));
			R_vn = Eigen::Vector3d(R_VXYZ(0, 0), R_VXYZ(1, 1), R_VXYZ(2, 2));
		}
		else if (_order == "RTK_POS")
		{
			double tmp; string s;
			ss >> tmp >> t >> XYZ[0] >> XYZ[1] >> XYZ[2] >> s >> nSat >> tmp >> RXYZ[0] >> RXYZ[1] >> RXYZ[2];
			ss >> ratio >> PDOP >> tmp >> Vxyz[0] >> Vxyz[1] >> Vxyz[2] >> RVxyz[0] >> RVxyz[1] >> RVxyz[2];
			
			if (s == "FIX")state = true;
			else state = false;
			ell = Cart2Geod(Vector3d(XYZ), false); v = Vector3d(Vxyz);
			vn = gsins::t_gbase::Cen(ell).transpose()*v;

			Eigen::Matrix3d R_XYZ = Eigen::Vector3d(RXYZ[0], RXYZ[1], RXYZ[2]).array().matrix().asDiagonal();
			Eigen::Matrix3d A = t_gbase::dGeod2Cart(gsins::t_gearth(), ell);
			Eigen::Matrix3d R_ELL = A.inverse()*R_XYZ*(A.inverse()).transpose();
			R_ell = Eigen::Vector3d(R_ELL(0, 0), R_ELL(1, 1), R_ELL(2, 2));

			Eigen::Matrix3d R_VXYZ = Eigen::Vector3d(RVxyz[0], RVxyz[1], RVxyz[2]).array().matrix().asDiagonal();
			Eigen::Matrix3d C = t_gbase::Cen(ell).transpose();
			Eigen::Matrix3d R = C*R_VXYZ*C.transpose();
			R_vn = Eigen::Vector3d(R(0, 0), R(1, 1), R(2, 2));

			//ell = Vector3d(XYZ); vn = Vector3d(Vxyz);
			//R_ell = Eigen::Vector3d(R_XYZ(0, 0), R_XYZ(1, 1), R_XYZ(2, 2));
			//R_vn = Eigen::Vector3d(R_VXYZ(0, 0), R_VXYZ(1, 1), R_VXYZ(2, 2));
		}
		else if (_order == "RTKLIB_POS")
		{
			string s;
			ss >> t >> XYZ[0] >> XYZ[1] >> XYZ[2] >> s >> nSat >> RXYZ[0] >> RXYZ[1] >> RXYZ[2];
			ss >> s >> s >> s >> s;
			ss >> ratio >> Vxyz[0] >> Vxyz[1] >> Vxyz[2] >> RVxyz[0] >> RVxyz[1] >> RVxyz[2];
			ss >> s >> s >> s >> PDOP;

			ell = Cart2Geod(Vector3d(XYZ), false); v = Vector3d(Vxyz);
			vn = gsins::t_gbase::Cen(ell).transpose()*v;

			Eigen::Matrix3d R_XYZ = Eigen::Vector3d(RXYZ[0], RXYZ[1], RXYZ[2]).array().matrix().asDiagonal();
			Eigen::Matrix3d A = t_gbase::dGeod2Cart(gsins::t_gearth(), ell);
			Eigen::Matrix3d R_ELL = A.inverse()*R_XYZ*(A.inverse()).transpose();
			R_ell = Eigen::Vector3d(R_ELL(0, 0), R_ELL(1, 1), R_ELL(2, 2));

			Eigen::Matrix3d R_VXYZ = Eigen::Vector3d(RVxyz[0], RVxyz[1], RVxyz[2]).array().matrix().asDiagonal();
			Eigen::Matrix3d C = t_gbase::Cen(ell).transpose();
			Eigen::Matrix3d R = C*R_VXYZ*C.transpose();
			R_vn = Eigen::Vector3d(R(0, 0), R(1, 1), R(2, 2));
		}
		map<string, t_gdata*>::iterator it = _data.begin();
		while (it != _data.end())
		{
			if (it->second->id_type() == t_gdata::LCI_POS)
				if (_order == "XYZ_VXYZ_RXYZ_RVXYZ_NSAT_PDOP" || _order == "RTKLIB_POS")
				{
					((t_gposdata*)it->second)->add_pos(t, ell, vn, R_ell, R_vn, PDOP, nSat);
				}
				else if (_order == "RTK_POS" || _order == "PVT_POS")
				{
					((t_gposdata*)it->second)->add_pos(t, ell, vn, R_ell, R_vn, PDOP, nSat, state);
				}
				else
					((t_gposdata*)it->second)->add_pos(t, ell, vn);
			
			it++;
		}
		if (ss.fail()) {
			if (_log) _log->comment(1, "posfile", "warning: incorrect POS data record: " + ss.str());
			else       cerr << "warning: incorrect POS data record: " << ss.str() << endl;
			t_gcoder::_consume(tmpsize);
			_mutex.unlock(); return -1;
		}
		t_gcoder::_consume(tmpsize);
		cnt++;
	}
	_mutex.unlock();
	return 0;
}
