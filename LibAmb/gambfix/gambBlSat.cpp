/**
*
* @verbatim
	History
	-1.0 jdhuang  2019-10-24  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambBlSat.cpp
* @brief	function for ambiguity.
*
* @author   Jiande Huang , Wuhan University
* @version	1.0.0
* @date		2019-10-24
*/

#include "gambfix/gambBlSat.h"
#include "gset/gsetamb.h"
#include "gset/gsetrec.h"

#include "newmat/newmatap.h"
namespace great
{
	great::t_gamb_baseline_sats::t_gamb_baseline_sats()
	{
	}

	great::t_gamb_baseline_sats::t_gamb_baseline_sats(set<string>& sats)
	{
		// 遍历所有基线

		for (auto iter : sats)
		{
			auto for_iter = ++sats.find(iter);
			while (for_iter != sats.end())
			{
				t_gamb_baseline tmp(iter, *for_iter, 1, 0.0);
				_baselines.push_back(tmp);
				for_iter++;
			}
		}

#ifdef DEBUG_AMBFIX
		int base_num = 0;
		cout << " ============> The baselines of sats ================" << endl;
		for (auto iter : _baselines)
		{
			base_num++;
			cout << setw(8) << base_num
				<< setw(4) << iter.nsys
				<< setw(5) << iter.obj1
				<< setw(5) << iter.obj2
				<< setw(16) << setprecision(8) << iter.distance 
				<< endl;
		}
#endif
	}

	great::t_gamb_baseline_sats::~t_gamb_baseline_sats()
	{

	}

	void t_gamb_baseline_sats::add_sat(const string& sat)
	{
		if (std::find(_sats.begin(), _sats.end(), sat) == std::end(_sats))
		{
			_sats.push_back(sat);
			std::sort(_sats.begin(), _sats.end());
		}
	}

	int great::t_gamb_baseline_sats::baseline_num()
	{
		return _baselines.size();
	}

	vector<pair<string,string>> great::t_gamb_baseline_sats::baselines()
	{
		vector<pair<string, string>> tmp;
		if (_baselines.empty())
		{
			if (_sats.empty()) return tmp;
			for (int i = 0; i < _sats.size(); i++)
			{
				for (int j = i + 1; j < _sats.size(); j++)
				{
					if (_sats[i].substr(0, 1) != _sats[j].substr(0, 1)) continue;
					t_gamb_baseline line(_sats[i], _sats[j], 0, 1);
					_baselines.push_back(line);
				}
			}
		}

		for (const auto& iter : _baselines)
		{
			tmp.push_back(iter.bl_pair());
		}

		// cout the base lines
		for (const auto& iter : tmp)
		{
			cout << "The baseline is " + iter.first + "   " + iter.second << endl;
		}
		return tmp;
	}

	vector<tuple<string, string, string, string>> t_gamb_baseline_sats::get_depend_Dd_baseline(vector<tuple<string, string, string, string>>& baselines)
	{
		
		vector<tuple<string, string, string, string>> result;

		// get all the objs
		map<string, int> index;
		int count = 0;
		for (auto& it : baselines)
		{
			string tmp_1st = get<0>(it);
			string tmp_2nd = get<1>(it);
			string tmp_3rd = get<2>(it);
			string tmp_4th = get<3>(it);

			if (index.count(tmp_1st) <= 0)
			{
				count++;
				index[tmp_1st] = count;
			}

			if (index.count(tmp_2nd) <= 0)
			{
				count++;
				index[tmp_2nd] = count;
			}
			
			if (index.count(tmp_3rd) <= 0)
			{
				count++;
				index[tmp_3rd] = count;
			}
			
			if (index.count(tmp_4th) <= 0)
			{
				count++;
				index[tmp_4th] = count;
			}

		}


		// check the depend
		Matrix vector_space;
		int ncol = index.size();
		for (auto& it : baselines)
		{
			string tmp_1st = get<0>(it);
			string tmp_2nd = get<1>(it);
			string tmp_3rd = get<2>(it);
			string tmp_4th = get<3>(it);

			auto tmp1 = index.count(tmp_1st);
			auto tmp2 = index.count(tmp_2nd);
			auto tmp3 = index.count(tmp_3rd);
			auto tmp4 = index.count(tmp_4th);
			int index1 = index[tmp_1st];
			int index2 = index[tmp_2nd];
			int index3 = index[tmp_3rd];
			int index4 = index[tmp_4th];

			if (tmp1 <= 0 || tmp2 <= 0 || tmp3 <= 0 || tmp4 <= 0) continue;

			if (vector_space.is_zero())
			{
				vector_space.resize(1, ncol);
				vector_space = 0.0;
				vector_space(1, index1) = 1;
				vector_space(1, index2) = -1;
				vector_space(1, index3) = -1;
				vector_space(1, index4) = 1;
				result.push_back(it);
			}

			RowVector new_vector(0.0, ncol);
			vector_space(1, index1) = 1;
			vector_space(1, index2) = -1;
			vector_space(1, index3) = -1;
			vector_space(1, index4) = 1;
			new_vector(index1) = 1;
			new_vector(index2) = -1;
			new_vector(index3) = -1;
			new_vector(index4) = 1;

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
			bool isDepend = (fabs(D(D.nrows(), D.ncols())) > 1.0e-2);

			if (A.nrows() != (vector_space.nrows() + new_vector.nrows())) continue;
			if (isDepend)
			{
				//cout << D << endl;
				vector_space = A;
				result.push_back(it);
			}

		}

		return result;

	}
}











