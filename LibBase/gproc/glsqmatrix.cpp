/**
*
* @verbatim
	History
	 -1.0 ZhengHJ  2018-12-02  created
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    glsqmatrix.cpp
* @brief		Matrix for NEQ,W,observ_equations
*
* @author       ZhengHJ, Wuhan University
* @version		1.0.0
* @date		    2019-04-17
*
*/

#include "gproc/glsqmatrix.h"
#include "Eigen/Dense"
#include <stdexcept>
#include <algorithm>
#include <emmintrin.h>
namespace great
{

	t_glsqEquationMatrix::t_glsqEquationMatrix()
	{

	}

	t_glsqEquationMatrix::~t_glsqEquationMatrix()
	{

	}



	void t_glsqEquationMatrix::add_equ(const vector<pair<int, double> >& B_value, const double& P_value, const double& l_value, const string& stie_name, const string& sat_name, const t_gobscombtype& obscombtype, const bool& is_newamb)
	{
		B.push_back(B_value);
		P.push_back(P_value);
		l.push_back(l_value);
		this->_site_sat_pairlist.push_back(make_pair(stie_name, sat_name));
		this->_obstypelist.push_back(obscombtype);
		this->_newamb_list.push_back(is_newamb);
	}

	void t_glsqEquationMatrix::add_equ(const t_glsqEquationMatrix & Other)
	{
		for (int i = 0; i < Other.num_equ(); i++)
		{
			this->add_equ(Other.B[i], Other.P[i], Other.l[i], Other._site_sat_pairlist[i].first, Other._site_sat_pairlist[i].second, Other._obstypelist[i], false);
		}
	}


	void t_glsqEquationMatrix::add_equ(const Matrix& B_value, const DiagonalMatrix& P_value, const ColumnVector& l_value, const vector<string>& site_name, const vector<string>&  sat_name, const vector<t_gobscombtype>& obstype)
	{

		if (B_value.Nrows() != P_value.Nrows() || B_value.Nrows() != l_value.Nrows() ||
			B_value.Nrows() != site_name.size() || B_value.Nrows() != sat_name.size() || B_value.Nrows() != obstype.size())
		{
			throw exception();
		}

		int num_equ = B_value.Nrows();
		//for B
		vector<pair<int, double> > temp;
		for (int row = 1; row <= B_value.Nrows(); row++)
		{
			temp.clear();
			for (int col = 1; col <= B_value.Ncols(); col++)
			{
				if (B_value(row, col) == 0.0)
				{
					continue;
				}
				temp.push_back(make_pair(col, B_value(row, col)));
			}

			B.push_back(temp);


			P.push_back(P_value(row));
			l.push_back(l_value(row));

			// Note: row from 1
			_site_sat_pairlist.push_back(make_pair(site_name[row - 1], sat_name[row - 1]));
			_obstypelist.push_back(obstype[row - 1]);
			_newamb_list.push_back(false);
		}
	}

	void t_glsqEquationMatrix::set_newamb(const int& idx, const bool& is_newamb)
	{
		if (idx >= this->num_equ())
		{
			throw exception(logic_error("input idx is more than number of equ"));
		}

		_newamb_list[idx] = is_newamb;
	}

	void t_glsqEquationMatrix::set_newamb(t_gturboedit* tb_slip)
	{
		if (!tb_slip) {
			return;
		}

		set<pair<string, string> > site_sat_list;
		for (int i = 0; i < _site_sat_pairlist.size(); i++) {
			if (!this->get_obscombtype(i).is_phase()) {
				continue;
			}
			string rec = this->get_sitename(i);
			string sat = this->get_satname(i);

			if (tb_slip->new_amb(rec, sat)) {
				_newamb_list[i] = true;
				site_sat_list.insert(make_pair(rec, sat));
			}
		}
		for (auto iter : site_sat_list) {
			tb_slip->set_new_amb(iter.first, iter.second, false);
		}
	}

	void t_glsqEquationMatrix::remove(const int& idx)
	{
		B.erase(B.begin() + idx - 1);
		P.erase(P.begin() + idx - 1);
		l.erase(l.begin() + idx - 1);

		_site_sat_pairlist.erase(_site_sat_pairlist.begin() + idx - 1);
		_obstypelist.erase(_obstypelist.begin() + idx - 1);
		_newamb_list.erase(_newamb_list.begin() + idx - 1);
	}

	void t_glsqEquationMatrix::remove_last_equation()
	{
		B.pop_back();
		P.pop_back();
		l.pop_back();
		_site_sat_pairlist.pop_back();
		_obstypelist.pop_back();
		_newamb_list.pop_back();
	}


	void t_glsqEquationMatrix::chageNewMat(Matrix& B_value, DiagonalMatrix& P_value, ColumnVector& l_value, const int& par_num)
	{
		B_value.ReSize(B.size(), par_num); B_value = 0.0;
		for (int row = 0; row < B.size(); row++)
		{
			for (int col = 0; col < B[row].size(); col++)
			{
				B_value(row + 1, B[row][col].first) = B[row][col].second;
			}
		}
		P_value.ReSize(P.size()); P_value = 0.0;
		for (int row = 0; row < P.size(); row++) {
			P_value(row + 1) = P[row];
		}
		l_value.ReSize(l.size()); l_value = 0.0;
		for (int row = 0; row < l.size(); row++) {
			l_value(row + 1) = l[row];
		}
	}

	void t_glsqEquationMatrix::chageNewMat(Matrix & B_value, SymmetricMatrix & P_value, ColumnVector & l_value, const int& par_num)
	{
		B_value.ReSize(B.size(), par_num); B_value = 0.0;
		for (int row = 0; row < B.size(); row++)
		{
			for (int col = 0; col < B[row].size(); col++)
			{
				B_value(row + 1, B[row][col].first) = B[row][col].second;
			}
		}
		P_value.ReSize(P.size()); P_value = 0.0;
		for (int row = 0; row < P.size(); row++) {
			P_value(row + 1, row + 1) = P[row];
		}
		l_value.ReSize(l.size()); l_value = 0.0;
		for (int row = 0; row < l.size(); row++) {
			l_value(row + 1) = l[row];
		}
	}

	void t_glsqEquationMatrix::print()
	{
		int count = 0;
		cout << "B is" << endl;
		for (int i = 1; i < B.size();) {
			count = 0;
			for (int j = 0; j < B[i].size(); j++) {
				count++;
				cout << scientific << setprecision(15) << setw(30) << B[i][j].second << "   ";
				if (count % 3 == 0)
					cout << endl;
			}
			cout << endl;
			i += 2;
		}
		cout << "P is" << endl;
		for (int i = 0; i < B.size(); i++) {
			cout << scientific << setprecision(15) << setw(30) << P[i] << endl;
		}
		cout << "l is" << endl;
		for (int i = 0; i < B.size(); ) {
			cout << scientific << setprecision(15) << setw(30) << l[i] << endl;
			i += 2;
		}
	}

	int t_glsqEquationMatrix::num_equ() const
	{
		return B.size();
	}

	double t_glsqEquationMatrix::res_equ() const
	{
		double ans = 0;
		for (int i = 0; i < num_equ(); i++) {
			ans += P[i] * l[i] * l[i];
		}
		return ans;
	}

	double t_glsqEquationMatrix::res_equ(bool phase) const
	{
		if (!phase)
		{
			return this->res_equ();
		}
		double ans = 0;
		for (int i = 0; i < num_equ(); i++) {
			if (get_obscombtype(i).is_code()) {
				continue;
			}
			ans += P[i] * l[i] * l[i];
		}
		return ans;
	}

	string t_glsqEquationMatrix::get_sitename(int equ_idx) const
	{
		return _site_sat_pairlist[equ_idx].first;
	}

	set<string> t_glsqEquationMatrix::get_sitelist(string rec) const
	{
		set<string> sat_list;
		int num_equ = this->num_equ();
		for (int i = 0; i < num_equ; i++) {
			if (_site_sat_pairlist[i].first == rec) {
				sat_list.insert(_site_sat_pairlist[i].second);
			}
		}

		return sat_list;
	}

	string t_glsqEquationMatrix::get_satname(int equ_idx) const
	{
		return _site_sat_pairlist[equ_idx].second;
	}

	string t_glsqEquationMatrix::get_obscombtype2str(int equ_idx) const
	{
		return _obstypelist[equ_idx].convert2str();
	}

	t_gobscombtype t_glsqEquationMatrix::get_obscombtype(int equ_idx) const
	{
		return _obstypelist[equ_idx];
	}

	vector<double> t_glsqEquationMatrix::get_codeomc() const
	{
		int equ_num = this->num_equ();
		vector<double> codeomc;
		for (int i = 0; i < equ_num; i++) 
		{
			if (_obstypelist[i].is_code()) 
			{
				codeomc.push_back(l[i]);
			}
		}
		return codeomc;
	}
	vector<double> t_glsqEquationMatrix::get_phaseomc() const
	{
		int equ_num = this->num_equ();
		vector<double> phaseomc;
		for (int i = 0; i < equ_num; i++)
		{
			if (_obstypelist[i].is_phase())
			{
				phaseomc.push_back(l[i]);
			}
		}
		return phaseomc;
	}

	

	bool t_glsqEquationMatrix::is_newamb(int equ_idx) const
	{
		return _newamb_list[equ_idx];
	}
	void print_equ_debinfo(const t_glsqEquationMatrix& equ)
	{
		for (int i = 0; i < equ.num_equ() / 2.0; i++)
		{
			cout << setw(5) << equ._site_sat_pairlist[2 * i].first
				<< setw(5) << equ._site_sat_pairlist[2 * i].second
				<< setiosflags(ios::fixed) << setprecision(4)
				<< setw(13) << equ.l[2 * i]
				<< setw(13) << equ.l[2 * i + 1]
				<< endl;
		}
	}

	vector<int> t_glsqEquationMatrix::find_equ(const string& site)
	{
		vector<int> ans;
		int num_equ = this->num_equ();
		for (int i = 0; i < num_equ; i++) {
			if (_site_sat_pairlist[i].first == site) {
				ans.push_back(i);
			}
		}
		return ans;
	}

	vector<int> t_glsqEquationMatrix::find_equ(const string& site, const string & sat)
	{
		vector<int> ans;
		int num_equ = this->num_equ();
		auto site_sat = make_pair(site, sat);
		for (int i = 0; i < num_equ; i++) {
			if (_site_sat_pairlist[i] == site_sat) {
				ans.push_back(i);
			}
		}
		return ans;
	}

	int t_glsqEquationMatrix::find_equ(const string& site, const string& sat, const t_gobscombtype& obscomtype)
	{
		int ans = -1;
		int num_equ = this->num_equ();
		auto site_sat = make_pair(site, sat);
		for (int i = 0; i < num_equ; i++) {
			if (_site_sat_pairlist[i] == site_sat && _obstypelist[i] == obscomtype) {
				ans = i;
				break;
			}
		}
		return ans;
	}

	void t_glsqEquationMatrix::clear_allequ()
	{
		B.clear();
		P.clear();
		l.clear();
		_site_sat_pairlist.clear();
		_obstypelist.clear();
		_newamb_list.clear();
	}

	L_SymmetricMatrix::L_SymmetricMatrix() :
		row_record(0),
		col_record(0)
	{

	}

	L_SymmetricMatrix::~L_SymmetricMatrix()
	{

	}

	int L_SymmetricMatrix::num() const
	{
		return _element.size();
	}

	double & L_SymmetricMatrix::num(int a, int b)
	{
		int row, col;
		if (a < b) {
			row = b;
			col = a;
		}
		else
		{
			row = a;
			col = b;
		}

		if (row >= row_record) {
			advance(row_it, row - row_record);
		}
		else {
			row_it = _element.begin();
			advance(row_it, row - 1);
		}


		if (row == row_record && col >= col_record) {
			advance(col_it, col - col_record);
		}
		else {
			col_it = row_it->begin();
			advance(col_it, col - 1);
		}
		row_record = row;
		col_record = col;

		return *col_it;
	}

	//double L_SymmericMatrix::num(int a, int b) const
	//{
	//	int row, col;
	//	if (a < b) {
	//		row = b;
	//		col = a;
	//	}
	//	else
	//	{
	//		row = a;
	//		col = b;
	//	}
	//	if (row >= row_record) {
	//		advance(row_iter, row - row_record);
	//	}
	//	else{
	//		advance(_element.begin(), row - 1);
	//	}
	//	row_record = row;

	//	if (col >= col_record) {
	//		advance(col_iter,col-b)
	//	}

	//	auto row_iter = _element.begin();
	//	advance(row_iter, row - 1);
	//	auto col_iter = row_iter->begin();
	//	advance(col_iter, col - 1);
	//	return *col_iter;
	//}

	void L_SymmetricMatrix::resize(int size)
	{
		for (int i = 1; i <= size; i++)
		{
			_element.push_back(list<double>(i, 0.0));
		}
		row_it = _element.begin();
		row_record = 1;
		col_it = row_it->begin();
		col_record = 1;
	}

	void L_SymmetricMatrix::add(const t_glsqEquationMatrix& equ)
	{
		// TODO
		for (int num = 0; num < equ.num_equ(); num++)
		{
			auto row_iter = _element.begin();
			advance(row_iter, equ.B[num][0].first - 1);
			for (int num_row = 0; num_row < equ.B[num].size(); num_row++)
			{
				if (num_row > 0) {
					advance(row_iter, equ.B[num][num_row].first - equ.B[num][num_row - 1].first);
				}
				auto col_iter = row_iter->begin();
				advance(col_iter, equ.B[num][0].first - 1);
				for (int num_col = 0; num_col < num_row + 1; num_col++)
				{
					if (num_col > 0) {
						advance(col_iter, equ.B[num][num_col].first - equ.B[num][num_col - 1].first);
					}
					*col_iter = *col_iter + equ.B[num][num_row].second * equ.P[num] * equ.B[num][num_col].second;
				}
			}
		}
		row_record = 1;
		row_it = _element.begin();
		col_record = 1;
		col_it = row_it->begin();
	}

	void L_SymmetricMatrix::addBackZero()
	{
		_element.push_back(list<double>(_element.size() + 1, 0.0));
		row_record = 1;
		row_it = _element.begin();
		col_record = 1;
		col_it = row_it->begin();
	}

	void L_SymmetricMatrix::remove(int idx)
	{
		if (idx < 1 || idx >_element.size()) {
			throw exception();
		}

		if (idx == 1)
		{
			_element.pop_front();
			for (auto iter = _element.begin(); iter != _element.end(); iter++)
			{
				iter->pop_front();
			}
		}
		else if (idx == _element.size())
		{
			_element.pop_back();
		}
		else {
			auto iterRow = _element.begin();
			advance(iterRow, idx - 1);
			iterRow = _element.erase(iterRow);
			while (iterRow != _element.end())
			{
				auto iterCol = iterRow->begin();
				advance(iterCol, idx - 1);
				iterRow->erase(iterCol);
				iterRow++;
			}
		}

		row_record = 1;
		row_it = _element.begin();
		col_record = 1;
		col_it = row_it->begin();
	}

	void L_SymmetricMatrix::print()
	{

		cout << setw(20) << setprecision(5);
		for (auto Row = _element.begin(); Row != _element.end(); Row++)
		{
			for (auto Col = Row->begin(); Col != Row->end(); Col++)
			{
				cout << setw(20) << *Col;
			}
			cout << endl;
		}

	}

	double L_SymmetricMatrix::center_value(int idx) const
	{
		if (idx<1 || idx > _element.size())
		{
			throw exception();
		}


		auto iter = _element.begin();
		advance(iter, idx - 1);
		return *(iter->rbegin());

	}

	SymmetricMatrix L_SymmetricMatrix::changeNewMat() const
	{
		SymmetricMatrix temp(_element.size());
		int row = 1, col = 1;
		for (auto Row = _element.begin(); Row != _element.end(); Row++)
		{
			col = 1;
			for (auto Col = Row->begin(); Col != Row->end(); Col++)
			{
				temp(row, col) = *Col;
				col++;
			}
			row++;
		}
		return temp;
	}


	L_ColumnVector::L_ColumnVector() :
		row_it(_element.begin()),
		row_record(1)
	{
	}

	L_ColumnVector::~L_ColumnVector()
	{

	}

	int L_ColumnVector::num() const
	{
		return _element.size();
	}

	double & L_ColumnVector::num(int a)
	{
		//auto row_iter = _element.begin();
		//advance(row_iter, a - 1);
		//return *row_iter;

		if (a >= row_record) {
			advance(row_it, a - row_record);
		}
		else {
			row_it = _element.begin();
			advance(row_it, a - 1);
		}
		row_record = a;

		return *row_it;
	}

	//double L_ColumnVector::num(int a) const
	//{
	//	auto row_iter = _element.begin();
	//	advance(row_iter, a - 1);
	//	return *row_iter;
	//}


	void L_ColumnVector::add(const t_glsqEquationMatrix& equ)
	{
		for (int num = 0; num < equ.num_equ(); num++)
		{
			auto iter_element = _element.begin();
			advance(iter_element, equ.B[num][0].first - 1);
			for (int num_par = 0; num_par < equ.B[num].size(); num_par++)
			{
				if (num_par > 0) {
					advance(iter_element, equ.B[num][num_par].first - equ.B[num][num_par - 1].first);
				}
				*iter_element = *iter_element + equ.B[num][num_par].second * equ.P[num] * equ.l[num];
			}
		}

		row_it = _element.begin();
		row_record = 1;
	}

	void L_ColumnVector::addBackZero()
	{
		_element.push_back(0.0);

		row_it = _element.begin();
		row_record = 1;
	}

	void L_ColumnVector::resize(int size)
	{
		for (int i = 1; i <= size; i++)
		{
			_element.push_back(0.0);
		}

		row_it = _element.begin();
		row_record = 1;
	}

	void L_ColumnVector::remove(int idx)
	{
		auto row = _element.begin();
		advance(row, idx - 1);
		_element.erase(row);

		row_it = _element.begin();
		row_record = 1;
	}


	void L_ColumnVector::print()
	{
		cout << setprecision(5);
		for (auto row = _element.begin(); row != _element.end(); row++)
		{
			cout << setw(20) << *row << endl;
		}
	}

	ColumnVector L_ColumnVector::changeNewMat()
	{
		ColumnVector temp(_element.size());
		int row = 1;
		for (auto Row = _element.begin(); Row != _element.end(); Row++)
		{
			temp(row) = *Row;
			row++;
		}
		return temp;
	}

	int V_SymmetricMatrix::num() const
	{
		return _element.size();
	}

	//double& V_SymmetricMatrix::num(int a, int b)
	//{
	//	int row = a, col = b;
	//	if (a < b) {
	//		row = b;
	//		col = a;
	//	}
	//	return _element[row - 1][col - 1];
	//}

	void V_SymmetricMatrix::resize(int num)
	{
		_element.clear();
		_element.reserve(num);
		for (int i = 0; i < num; i++)
		{
			_element.push_back(vector<double>(i + 1, 0.0));
		}
	}

	void V_SymmetricMatrix::add(const t_glsqEquationMatrix & equ)
	{
		// cycle equ
		for (int num = 0; num < equ.num_equ(); num++)
		{
			// cycle par
			auto B_temp = equ.B[num];
			sort(B_temp.begin(), B_temp.end());
			for (int ipar = 0; ipar < B_temp.size(); ipar++)
			{
				auto row = B_temp[ipar];
				for (int jpar = 0; jpar <= ipar; jpar++)
				{
					auto col = B_temp[jpar];
					_element[row.first - 1][col.first - 1] += row.second * equ.P[num] * col.second;
				}
			}
		}
	}

	void V_SymmetricMatrix::del(const t_glsqEquationMatrix& equ)
	{
		// cycle equ
		for (int num = 0; num < equ.num_equ(); num++)
		{
			// cycle par
			auto B_temp = equ.B[num];
			sort(B_temp.begin(), B_temp.end());
			for (int ipar = 0; ipar < B_temp.size(); ipar++)
			{
				auto row = B_temp[ipar];
				for (int jpar = 0; jpar <= ipar; jpar++)
				{
					auto col = B_temp[jpar];
					_element[row.first - 1][col.first - 1] -= row.second * equ.P[num] * col.second;
				}
			}
		}
	}

	void V_SymmetricMatrix::add(const t_glsqEquationMatrix& equ, bool phase)
	{
		if (!phase)
		{
			this->add(equ);
			return;
		}
		// cycle equ
		for (int num = 0; num < equ.num_equ(); num++)
		{
			if (equ.get_obscombtype(num).is_code()) {
				continue;
			}
			// cycle par 
			auto B_temp = equ.B[num];
			sort(B_temp.begin(), B_temp.end());
			for (int ipar = 0; ipar < B_temp.size(); ipar++)
			{
				auto row = B_temp[ipar];
				for (int jpar = 0; jpar <= ipar; jpar++)
				{
					auto col = B_temp[jpar];
					_element[row.first - 1][col.first - 1] += row.second * equ.P[num] * col.second;
				}
			}
		}
	}

	void V_SymmetricMatrix::addBackZero()
	{
		_element.push_back(vector<double>(_element.size() + 1, 0.0));
	}

	void V_SymmetricMatrix::remove(int idx)
	{
		for (int row = _element.size(); row > idx; row--) 
		{
			_element[row - 1].erase(_element[row - 1].begin() + idx - 1);
		}
		_element.erase(_element.begin() + idx - 1);
	}

	void V_SymmetricMatrix::print()
	{
		cout << setw(20) << setprecision(5);
		for (auto Row = _element.begin(); Row != _element.end(); Row++)
		{
			for (auto Col = Row->begin(); Col != Row->end(); Col++)
			{
				cout << setw(20) << *Col;
			}
			cout << endl;
		}
	}


	SymmetricMatrix V_SymmetricMatrix::changeNewMat() const
	{
		SymmetricMatrix temp(_element.size());
		int row = 1, col = 1;
		for (auto Row = _element.begin(); Row != _element.end(); Row++)
		{
			col = 1;
			for (auto Col = Row->begin(); Col != Row->end(); Col++)
			{
				temp(row, col) = *Col;
				col++;
			}
			row++;
		}
		return temp;
	}

	double V_SymmetricMatrix::center_value(int idx) const
	{
		return _element[idx - 1][idx - 1];
	}

	void V_SymmetricMatrix::swap(int a, int b)
	{
		double temp_a = num(a, a);
		double temp_b = num(b, b);
		double temp_ab = num(a, b);
		for (int i = 0; i < num(); i++)
		{
			double temp = num(a, i + 1);
			num(a, i + 1) = num(b, i + 1);
			num(b, i + 1) = temp;
		}
		num(a, a) = temp_b;
		num(b, b) = temp_a;
		num(a, b) = temp_ab;
	}

	int V_ColumnVector::num() const
	{
		return _element.size();
	}

	double & V_ColumnVector::num(int a)
	{
		return _element[a - 1];
	}

	void V_ColumnVector::resize(int num)
	{
		_element.resize(num, 0.0);
	}

	void V_ColumnVector::add(const t_glsqEquationMatrix & equ)
	{
		for (int num = 0; num < equ.num_equ(); num++)
		{
			for (int ipar = 0; ipar < equ.B[num].size(); ipar++)
			{
				_element[equ.B[num][ipar].first - 1] += equ.B[num][ipar].second * equ.P[num] * equ.l[num];
			}
		}
	}

	void V_ColumnVector::del(const t_glsqEquationMatrix& equ)
	{
		for (int num = 0; num < equ.num_equ(); num++)
		{
			for (int ipar = 0; ipar < equ.B[num].size(); ipar++)
			{
				_element[equ.B[num][ipar].first - 1] -= equ.B[num][ipar].second * equ.P[num] * equ.l[num];
			}
		}
	}

	void V_ColumnVector::add(const t_glsqEquationMatrix& equ, bool phase)
	{
		if (!phase)
		{
			this->add(equ);
			return;
		}
		for (int num = 0; num < equ.num_equ(); num++)
		{
			if (equ.get_obscombtype(num).is_code()) {
				continue;
			}
			for (int ipar = 0; ipar < equ.B[num].size(); ipar++)
			{
				_element[equ.B[num][ipar].first - 1] += equ.B[num][ipar].second * equ.P[num] * equ.l[num];
			}
		}
	}

	void V_ColumnVector::addBackZero()
	{
		_element.push_back(0.0);
	}

	void V_ColumnVector::remove(int idx)
	{
		_element.erase(_element.begin() + idx - 1);
	}

	void V_ColumnVector::print()
	{
		cout << setprecision(5);
		for (auto row = _element.begin(); row != _element.end(); row++)
		{
			cout << setw(20) << *row << endl;
		}
	}

	ColumnVector V_ColumnVector::changeNewMat()
	{
		ColumnVector temp(_element.size());
		int row = 1;
		for (auto Row = _element.begin(); Row != _element.end(); Row++)
		{
			temp(row) = *Row;
			row++;
		}
		return temp;
	}

	void V_ColumnVector::swap(int a, int b)
	{
		double temp;
		temp = _element[a - 1];
		_element[a - 1] = _element[b - 1];
		_element[b - 1] = temp;
	}

	bool remove_lsqmatrix(int idx, V_SymmetricMatrix& NEQ, V_ColumnVector& W)
	{
		if (NEQ.num() != W.num())
			return false;
		if (idx <1 || idx > NEQ.num())
			return false;

		double center_value = NEQ.center_value(idx);

		if (!double_eq(center_value, 0.0))
		{
			int NEQ_num = NEQ.num();

			// record need remove element
			vector<double>  remove_element;
			remove_element.reserve(NEQ_num);
			
			double w_remove = 0.0;
			vector<int> remove_rows;
			for (int col = 1; col <= NEQ_num; col++)
			{
				remove_element.push_back(NEQ.num(idx, col));
				if (col == idx || double_eq(NEQ.num(idx, col), 0.0)) 
				{
					continue;
				}
				remove_rows.push_back(col - 1);
			}

			w_remove = W._element[idx - 1];
			
			// remove the par in NEQ and W by for cycle
#ifdef USE_OPENMP
#pragma	omp parallel for schedule(dynamic)
#endif

			int remove_rows_size = remove_rows.size();
			for (int i = 0; i < remove_rows_size; ++i)
			{
				int row = remove_rows[i];
				double temp_coeff = -remove_element[row] / center_value;
				for (int col = 0; col <= row; col++)
				{
					NEQ._element[row][col] += remove_element[col] * temp_coeff;
				}

				W._element[row] += w_remove * temp_coeff;
			}

		}

		NEQ.remove(idx);
		W.remove(idx);

		return true;
	}
	bool BASE_LIBRARY_EXPORT rearrange_lsqmatrix(const vector<int>& remove_idx, V_SymmetricMatrix& NEQ, V_ColumnVector& W, Eigen::MatrixXd& N11, Eigen::MatrixXd& N21, Eigen::MatrixXd& N22, Eigen::VectorXd& W1, Eigen::VectorXd& W2, bool idx_from_zero)
	{
		int matrix_begin = idx_from_zero ? 0 : 1;
		if (remove_idx.size() == 0) {
			return true;
		}
		int rows = NEQ.num();
		int rows_new = rows - remove_idx.size();


		vector<bool> idx_map(rows, false);
		for (int i = 0; i < remove_idx.size(); i++)
		{
			idx_map[remove_idx[i] - 1] = true;
		}

		// rearrange
		int i_skip = 0;
		for (int i = 0; i < rows; i++)
		{
			int j_skip = 0;
			if (idx_map[i] == true)
			{
				i_skip++;
				W2(i_skip - 1 + matrix_begin) = W._element[i];
				for (int j = 0; j <= i; j++)
				{
					if (idx_map[j] == true)
					{
						j_skip++;
						N22(i_skip - 1 + matrix_begin, j_skip - 1 + matrix_begin) = NEQ._element[i][j];
						N22(j_skip - 1 + matrix_begin, i_skip - 1 + matrix_begin) = NEQ._element[i][j];
					}
					else
					{
						N21(i_skip - 1 + matrix_begin, j - j_skip + matrix_begin) = NEQ._element[i][j];
					}
				}
			}
			else
			{
				W1(i - i_skip + matrix_begin) = W._element[i];
				for (int j = 0; j <= i; j++)
				{
					if (idx_map[j] == true)
					{
						j_skip++;
						N21(j_skip - 1 + matrix_begin, i - i_skip + matrix_begin) = NEQ._element[i][j];
					}
					else
					{
						N11(i - i_skip + matrix_begin, j - j_skip + matrix_begin) = NEQ._element[i][j];
						N11(j - j_skip + matrix_begin, i - i_skip + matrix_begin) = NEQ._element[i][j];
					}
				}
			}
		}
		return true;
	}
	bool BASE_LIBRARY_EXPORT rearrange_lsqmatrix(vector<int>& remove_idx, V_SymmetricMatrix & NEQ, V_ColumnVector & W, t_gallpar& allpar)
	{
		if (remove_idx.size() == 0) 
		{
			return true;
		}

		int rows = NEQ.num();
		int rows_new = rows - remove_idx.size();


		vector<bool> idx_map(rows, false);
		int remove_size = remove_idx.size();
		for (int i = 0; i < remove_size; i++)
		{
			idx_map[remove_idx[i] - 1] = true;
		}

		// rearrange 
		t_gallpar allpar_orig = allpar; allpar.delAllParam();
		for (int i = 0; i < rows; i++) 
		{
			if (idx_map[i]) 
			{
				continue;
			}
			allpar.addParam(allpar_orig[i]);
		}
		for (int i = 0; i < remove_size; i++) 
		{
			allpar.addParam(allpar_orig[remove_idx[i]-1]);
		}

		// rearrange remove_idx
		remove_idx.clear();
		for (int i = 0; i < remove_size; i++) 
		{
			remove_idx.push_back(rows_new + 1 + i);
		}

		V_SymmetricMatrix NEQ_orig = NEQ;
		V_ColumnVector W_orig = W;
		
		// rearrange
		int i_skip = 0;
		for (int i = 0; i < rows; i++)
		{
			int j_skip = 0;
			if (idx_map[i] == true)
			{
				i_skip++;
				W._element[rows_new+i_skip - 1] = W_orig._element[i];
				
				for (int j = 0; j <= i; j++)
				{
					if (idx_map[j] == true)
					{
						j_skip++;
						NEQ.num(rows_new + i_skip , rows_new + j_skip) = NEQ_orig._element[i][j];
					}
					else
					{
						NEQ._element[ rows_new + i_skip - 1][j - j_skip ] = NEQ_orig._element[i][j];
					}
				}
			}
			else
			{
				W._element[i - i_skip] = W_orig._element[i];
				for (int j = 0; j <= i; j++)
				{
					if (idx_map[j] == true)
					{
						j_skip++;
						NEQ.num(rows_new + j_skip ,i - i_skip +1) = NEQ_orig._element[i][j];
					}
					else
					{
						NEQ._element[i - i_skip][j - j_skip] = NEQ_orig._element[i][j];
					}
				}
			}
		}
		return true;
	}
	
	bool BASE_LIBRARY_EXPORT remove_lsqmatrix(const vector<int>& remove_idx, V_SymmetricMatrix& NEQ, V_ColumnVector& W)
	{
		if (remove_idx.size() == 0) 
		{
			return true;
		}
		int rows = NEQ.num();
		int rows_new = rows - remove_idx.size();

		// re-arrange
		//Eigen::MatrixXd NEQ_re = Eigen::MatrixXd::Zero(rows, rows);
		//Eigen::VectorXd W_re = Eigen::VectorXd::Zero(rows);

		//vector<bool> idx_map(rows, false);
		//for (int i = 0; i < remove_idx.size(); i++)
		//{
		//	idx_map[remove_idx[i]-1] = true;
		//	//std::cout << "remove:" << remove_idx[i] << std::endl;
		//}

		//// ÖØÐÂÅÅÁÐ¾ØÕó
		//int i_skip = 0;
		//for (int i = 0; i < rows; i++)
		//{
		//	int j_skip = 0;
		//	if (idx_map[i] == true)
		//	{
		//		i_skip++;
		//		W_re(rows_new + i_skip - 1) = W.num(i + 1);
		//		for (int j = 0; j <= i; j++)
		//		{
		//			if (idx_map[j] == true)
		//			{
		//				j_skip++;
		//				NEQ_re(rows_new + i_skip - 1, rows_new + j_skip - 1) = NEQ.num(i + 1, j + 1);
		//				NEQ_re(rows_new + j_skip - 1, rows_new + i_skip - 1) = NEQ.num(i + 1, j + 1);

		//			}
		//			else
		//			{
		//				NEQ_re(rows_new + i_skip - 1, j - j_skip) = NEQ.num(i + 1, j + 1);
		//				NEQ_re(j - j_skip, rows_new + i_skip - 1) = NEQ.num(i + 1, j + 1);

		//			}
		//		}
		//	}
		//	else
		//	{
		//		W_re(i - i_skip) = W.num(i + 1);
		//		for (int j = 0; j <= i; j++)
		//		{
		//			if (idx_map[j] == true)
		//			{
		//				j_skip++;
		//				NEQ_re(i - i_skip, rows_new + j_skip - 1) = NEQ.num(i + 1, j + 1);
		//				NEQ_re(rows_new + j_skip - 1, i - i_skip) = NEQ.num(i + 1, j + 1);

		//			}
		//			else
		//			{
		//				NEQ_re(i - i_skip, j - j_skip) = NEQ.num(i + 1, j + 1);
		//				NEQ_re(j - j_skip, i - i_skip) = NEQ.num(i + 1, j + 1);

		//			}
		//		}
		//	}
		//}


		Eigen::MatrixXd NEQ_new;
		Eigen::MatrixXd W_new;

		Eigen::MatrixXd A = Eigen::MatrixXd::Zero(rows_new, rows_new);
		Eigen::MatrixXd B = Eigen::MatrixXd::Zero(rows - rows_new, rows_new);
		Eigen::MatrixXd C = Eigen::MatrixXd::Zero(rows - rows_new, rows - rows_new);

		Eigen::VectorXd W1 = Eigen::VectorXd::Zero(rows_new);
		Eigen::VectorXd W2 = Eigen::VectorXd::Zero(rows - rows_new);

		// rearrange
		rearrange_lsqmatrix(remove_idx, NEQ, W, A, B, C, W1, W2, true);


		// remove 
		Eigen::MatrixXd Ci = C.inverse();
		Eigen::MatrixXd BtCi = B.transpose()*Ci;

		NEQ_new = A - BtCi * B;
		W_new = W1 - BtCi * W2;

		// get new res_obs
		//if (res_obs)
		//{
		//	Eigen::MatrixXd W_remove = W_re.tail(rows - rows_new);
		//	*res_obs = *res_obs - (W_remove.transpose()*Ci*W_remove)(1,1);
		//}

		// assing ouptut matrix
		for (int i = 0; i < rows_new; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				NEQ.num(i + 1, j + 1) = NEQ_new(i, j);
				NEQ.num(j + 1, i + 1) = NEQ_new(i, j);
			}
			W.num(i + 1) = W_new(i);
		}

		// Resize
		for (int i = rows - 1; i >= rows_new; i--)
		{
			NEQ.remove(i + 1);
			W.remove(i + 1);
		}

		return true;
	}
}