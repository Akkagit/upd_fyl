/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambBlRec.cpp
* @brief	function for ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#include "gambfix/gambBlRec.h"
#include "gambfix/gambBl.h"
#include "gall/gallobj.h"
#include "gset/gsetamb.h"
#include "gset/gsetrec.h"
#include "gutils/gstring.h"
#include "gutils/gleo.h"
#include <algorithm>

#include "newmat/newmatap.h"

using namespace std;

namespace great
{
	great::t_gamb_baseline_recs::t_gamb_baseline_recs()
	{
	}

	t_gamb_baseline_recs::t_gamb_baseline_recs(const t_gtime& cut_time, double cut_length)
	{
		_cut_length = cut_length;
		_cut_time   = cut_time;
	}


	great::t_gamb_baseline_recs::~t_gamb_baseline_recs()
	{

	}

	void t_gamb_baseline_recs::add_cut(const t_gtime& cut_time, double cut_length)
	{
		_cut_length = cut_length;
		_cut_time = cut_time;
	}

	void t_gamb_baseline_recs::add_gobj(t_gallobj* obj)
	{
		this->_gobj = obj;

	}

	void t_gamb_baseline_recs::add_gset(t_gsetbase* set)
	{
		this->_gset = set;
		_beg_time = dynamic_cast<t_gsetgen*>(set)->beg();
		_end_time = dynamic_cast<t_gsetgen*>(set)->end();
	}

	void t_gamb_baseline_recs::add_glog(t_glog* log)
	{
		this->_glog = log;
	}

	void t_gamb_baseline_recs::add_rec(const string& rec, const string& sat)
	{
		if (std::find(_recs.begin(),_recs.end(),rec) == std::end(_recs))
		{
			_recs.push_back(rec);
			std::sort(_recs.begin(), _recs.end());
		}
		
		if (_sys[rec].count(sat.substr(0, 1)) < 1) _sys[rec].insert(sat.substr(0, 1));
	}

	vector<pair<string, string>> t_gamb_baseline_recs::baselines()
	{
		cout << "get the baselines for recs " << endl;

		vector<pair<string, string>> tmp;
		
		if (_baselines.empty())
		{
			cout << "get rec xyz  " << endl;
			this->_get_rec_xyz();

			cout << "create_baseline  " << endl;
			this->_create_baseline();
		}

		for (const auto& iter : _baselines)
		{
			tmp.push_back(iter.bl_pair());
		}


		// cout the base line
		for (const auto& iter : tmp)
		{
			cout << "The baseline is " + iter.first + "   " + iter.second << endl;
		}

		return tmp;

	}


	bool great::t_gamb_baseline_recs::_check_amb_depend()
	{
		// rec map
		map<string, int> _obj_index;
		for (int i = 0; i < _recs.size(); i++)
		{
			_obj_index[_recs[i]] = i + 1;
		}

		// judge
		vector<t_gamb_baseline> Independent_baseline;
		Matrix	vector_space;
		for (const auto& iter : _baselines)
		{
			const auto& line = iter.bl_pair();
			const auto& tmp1 = _obj_index.count(line.first);
			const auto& tmp2 = _obj_index.count(line.second);
			if (tmp1 <= 0 || tmp2 <= 0) continue;

			int ncol = _obj_index.size();
			if (vector_space.is_zero())
			{
				vector_space.resize(1, ncol);
				vector_space = 0.0;
				int index1 = _obj_index[line.first];
				int index2 = _obj_index[line.second];
				vector_space(1, index1) = 1;
				vector_space(1, index2) = -1;
				
				Independent_baseline.push_back(iter);
				continue;
			}

			RowVector new_vector(0.0, ncol);
			int index1 = _obj_index[line.first];
			int index2 = _obj_index[line.second];
			new_vector(index1) = 1;
			new_vector(index2) = -1;

			// 数学原理：奇异矩阵经过SVD分解后，D矩阵对角线的非零元素个数等于A矩阵的矩阵秩
			Matrix A = vector_space & new_vector;
			DiagonalMatrix D;
			if (A.nrows() < A.ncols())
			{
				SVD(A.t(), D);
			}
			else
			{
				SVD(A, D);
			}
			//cout << D << endl;
			bool isDepend = (fabs(D(D.nrows(), D.ncols())) > 0.01);

			if (A.nrows() != (vector_space.nrows() + new_vector.nrows())) return false;
			if (isDepend)
			{
				vector_space = A;	
				Independent_baseline.push_back(iter);
			}
		}

		_baselines = Independent_baseline;

		return true;
	}

	bool great::t_gamb_baseline_recs::_get_rec_xyz()
	{
		if (this->_gobj == nullptr && this->_gset == nullptr) return false;

		// 优先级别为，先从sinex找，再从obs找，最后从xml找

		for (const auto& crt : this->_recs)
		{
			bool   isfind = false;
			const  string&  rec = str2upper(crt);

			// 找xml
			if (this->_gset)
			{
				t_gtriple tmp = dynamic_cast<t_gsetrec*>(_gset)->get_crd_xyz(rec);
				if (!tmp.zero())
				{
					cout << rec << "  find xyz in xml file  " << endl;
					this->_xyz[rec] = tmp;
					isfind = true;
				}
			}

			if (isfind == false && this->_gobj)
			{
				shared_ptr<t_gobj> grec = this->_gobj->obj(rec);
				if (grec)
				{
					t_gtime tmp = _cut_time;
					while (tmp <= _end_time)
					{
						//t_gtriple rec_crd = grec->crd(_cut_time);
						t_gtriple rec_crd;
						t_gtriple rec_std;
						grec->get_recent_crd(_cut_time, 100, rec_crd, rec_std);
						if (!rec_crd.zero())
						{
							cout << rec << "  find xyz in snx file  " << endl;
							this->_xyz[rec] = rec_crd;
							isfind = true;
						}

						if (isfind) break;
						tmp = tmp + 300;
					}
				}
			}
		}
		return true;
	}

	bool t_gamb_baseline_recs::_create_baseline()
	{
		if (_recs.empty())
		{
			cout << "  rec list is empty, return false  " << endl;
			return false;
		}

		const auto& loop_1st = _recs;
		const auto& loop_2nd = _recs;

		// 存储有效测站列表
		int index = 0;
		vector<string> tmp;
		for (const auto& iter : _recs)
		{
			if (!IS_SUPPORTED_LEO(strs2strl(iter)) &&this->_xyz.count(iter) == 0) continue;
			tmp.push_back(iter);
		}
		_recs = tmp;

#ifdef DEBUG_AMBFIX
		cout << endl;
		cout << " ============> The baselines of recs(before sort by distance and check dependence) ================" << endl;
		// 输出有效测站
		for (auto crt_obj = objs.begin(); crt_obj != objs.end(); crt_obj++)
		{
			cout << setw(8) << crt_obj->first
				<< setw(8) << crt_obj->second
				<< endl;
		}
		cout << endl;
#endif

		//构建所有可能的基线
		for (int i = 0; i < _recs.size(); i++)
		{
			const string& rec_1st = _recs[i];
			//测站没有对应观测卫星
			if (_sys[rec_1st].size() == 0)
			{
				cout << "  no system for " + rec_1st + " skiped" << endl;
				continue;
			}
			for (int j = i + 1; j < _recs.size(); j++)
			{
				// 简单判断
				const string& rec_2nd = _recs[j];
				if (_sys[rec_2nd].size() == 0) continue;

				// 判断距离是否超限
				// 单位转换为千米
				t_gtriple first  = this->_xyz[rec_1st];
				t_gtriple second = this->_xyz[rec_2nd];
				t_gtriple D_value = first - second;
				double distance = D_value.norm() / 1000.0;
				if (IS_SUPPORTED_LEO(strs2strl(rec_1st)) || IS_SUPPORTED_LEO(strs2strl(rec_2nd))) distance = 1;
				if (distance > fabs(_cut_length)) continue;

				//判断两个测站之间有没有共同的卫星系统
				int nsys = 0;
				for (const auto& iter : _sys[rec_1st])
				{
					if (_sys[rec_2nd].find(iter) != std::end(_sys[rec_2nd])) nsys++;
				}
				if (nsys == 0) continue;

				//构建基线
				t_gamb_baseline bl(rec_1st, rec_2nd, distance, nsys);
				_baselines.push_back(bl);
			}
		}

		// 如果有有效基线，我们就赋给私有变量
		if (_baselines.empty()) return false;

#ifdef DEBUG_AMBFIX
		int tmp1 = 0;
		cout << endl;
		cout << " ============> The baselines of recs(before sort by distance and check dependence) ================" << endl;
		for (auto iter : sort_bl)
		{
			tmp1++;
			cout << setw(8) << tmp1
				<< setw(4) << iter.nsys
				<< setw(5) << iter.obj1
				<< setw(5) << iter.obj2
				<< setw(16) << setprecision(8) << iter.distance
				<< endl;
		}
		cout << endl;
#endif

		//按距离值从小到大排序
		//按系统数从大到小排序
		std::sort(_baselines.begin(), _baselines.end(), compare);

		//选取独立基线
		_check_amb_depend();
#ifdef DEBUG_AMBFIX
		int tmp2 = 0;
		cout << endl;
		cout << " ============> The baselines of recs(after sort by distance) ================" << endl;
		for (auto iter : sort_bl)
		{
			tmp2++;
			cout << " Rec BaseLines : "
				<< setw(8) << tmp2
				<< setw(4) << iter.nsys
				<< setw(5) << iter.obj1
				<< setw(5) << iter.obj2
				<< setw(16) << setprecision(8) << iter.distance
				<< endl;
		}
		cout << endl;
#endif

#ifdef DEBUG_AMBFIX
		int tmp3 = 0;
		cout << endl;
		cout << " ============> The baselines of recs(after _check_amb_depend) ================" << endl;
		for (auto iter : this->_baselines)
		{
			tmp3++;
			cout <<  " Rec BaseLines : "
				<< setw(8) << tmp3
				<< setw(4) << iter.nsys
				<< setw(5) << iter.obj1
				<< setw(5) << iter.obj2
				<< setw(16) << setprecision(8) << iter.distance
				<< endl;
		}
		cout << endl;
#endif
		return true;
	}

	bool compare(const t_gamb_baseline& bl1, const t_gamb_baseline& bl2)
	{
		if (bl1.nsys > bl2.nsys)
		{
			return true;
		}

		if ((bl1.nsys == bl2.nsys) && (bl1.distance < bl2.distance))
		{
			return true;
		}

		return false;
	}

}

