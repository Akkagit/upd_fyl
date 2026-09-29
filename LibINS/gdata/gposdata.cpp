#include "gposdata.h"

gsins::t_gposdata::t_gposdata()
{
	gtrace("t_gposdata::t_gposdata");
	id_type(LCI_POS);
	_ptr = 0;
}

gsins::t_gposdata::~t_gposdata()
{
}

int gsins::t_gposdata::add_pos(const double & t, const Eigen::Vector3d & pos, const Eigen::Vector3d & v)
{
	gtrace("t_gposdata::add_pos");

#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	data_pos tmp = { t,pos,v };
	_vecpos.push_back(tmp);
	_setpos.insert(tmp);
	_gmutex.unlock();
	return 0;
}

int gsins::t_gposdata::add_pos(const double & t, const Eigen::Vector3d & pos, const Eigen::Vector3d & v, const Eigen::Vector3d & R_p, const Eigen::Vector3d & R_v, const double & pdop, const int & n)
{
	gtrace("t_gposdata::add_pos");

#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	data_pos tmp = { t,pos,v ,R_p,R_v,pdop,n };
	_setpos.insert(tmp);
	_gmutex.unlock();
	return 0;
}

int gsins::t_gposdata::add_pos(const double & t, const Eigen::Vector3d & pos, const Eigen::Vector3d & v, const Eigen::Vector3d & R_p, const Eigen::Vector3d & R_v, const double & pdop, const int & n, const bool & state)
{
	gtrace("t_gposdata::add_pos");

#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	data_pos tmp = { t,pos,v ,R_p,R_v,pdop,n,state };
	_setpos.insert(tmp);
	_gmutex.unlock();
	return 0;
}

void gsins::t_gposdata::load(double & t, Eigen::Vector3d & pos, Eigen::Vector3d & v, bool direction)
{
	if (direction)
	{
		t = _vecpos[_ptr].t; pos = _vecpos[_ptr].pos; v = _vecpos[_ptr].vn;
		_ptr++;
	}
	else
	{
		_ptr--;
		t = _vecpos[_ptr].t; pos = _vecpos[_ptr].pos; v = _vecpos[_ptr].vn;
	}
}

void gsins::t_gposdata::load(double & t, Eigen::Vector3d & pos, Eigen::Vector3d & v, bool direction, double t1)
{
	set<data_pos>::const_iterator it_index, it_start = _setpos.begin();
	data_pos tmp{ t1, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero() };
	it_index = _setpos.upper_bound(tmp);
	if (it_index == _setpos.end())
	{
		t = 0.0; pos = v = Eigen::Vector3d::Zero();
	}
	it_index--;
	t = it_index->t; pos = it_index->pos; v = it_index->vn;
	_ptr++;
	
}

void gsins::t_gposdata::load(double & t, Eigen::Vector3d & pos, Eigen::Vector3d & v, Eigen::Vector3d& Qpos, Eigen::Vector3d& Qv, 
	double & PDOP, int & nSat, bool direction, double t1)
{
	set<data_pos>::const_iterator it_index, it_start = _setpos.begin();
	data_pos tmp{ t1, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero() };
	it_index = _setpos.upper_bound(tmp);
	if (it_index == _setpos.end())
	{
		t = 0.0; pos = v = Eigen::Vector3d::Zero();
	}
	it_index--;
	t = it_index->t; pos = it_index->pos; v = it_index->vn; PDOP = it_index->PDOP; nSat = it_index->nSat;
	Qpos = it_index->Rpos; Qv = it_index->Rvn;
	_ptr++;
}

void gsins::t_gposdata::load(data_pos & data, double t1)
{
	int  nSat; double t, PDOP; Eigen::Vector3d pos, v, Qpos, Qv; bool state;
	set<data_pos>::const_iterator it_index, it_start = _setpos.begin();
	data_pos tmp{ t1, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero() };
	it_index = _setpos.upper_bound(tmp);
	if (it_index == _setpos.end())
	{
		data = { 0 }; return;
	}
	it_index--;
	t = it_index->t; pos = it_index->pos; v = it_index->vn; PDOP = it_index->PDOP; nSat = it_index->nSat;
	Qpos = it_index->Rpos; Qv = it_index->Rvn;
	state = it_index->amb_state;
	data = { t,pos,v,Qpos,Qv,PDOP,nSat ,state };
}

bool gsins::t_gposdata::IsValid(const double & t, const Eigen::Vector3d & pos, const double& tot, const double& threshold)
{
	if (double_eq(t, 0.0) || pos == Eigen::Vector3d::Zero())return false;
	if (fabs(tot - t) < threshold)return true;
	return false;
}

double gsins::t_gposdata::beg_obs()
{
	if (_setpos.size() == 0)return 0;
	return _setpos.begin()->t;
}

double gsins::t_gposdata::end_obs()
{
	if (_setpos.size() == 0)return 0;
	return _setpos.rbegin()->t;
}

