/**
*
* @verbatim
    History
     -1.0 Wei Zhang     2019-01-20 creat the file.
     -1.1 Wei Zhang     2019-03-27 Adding Doxygen Style Code Remarks
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gscf2crsleo.cpp
* @brief    calculate the rotation matrix from SCF to CRS
*
* @author      Wei Zhang, Wuhan University
* @version		1.0.0
* @date		2019-03-28
*
*/


#include "gutils/gscf2crsleo.h"
#include "gutils/gtypeconv.h"
#include "gdata/gattitude.h"
#include "gmodels/ginterp.h"

using namespace std;
namespace great
{
	/** @brief default constructor. */
	t_gscf2crsleo::t_gscf2crsleo(string satname)
	{
		_sat_name = satname;
		_isAttValid = false;
		_rotMat.resize(3, 3);
		_rotMatScf2crs.resize(3, 3);
		_epoch = t_gtime();
	}

	/** @brief default destructor. */
	t_gscf2crsleo::~t_gscf2crsleo()
	{
	}

	/**
	* @brief calculate the nominal rotation matrix from SCF to CRS(without attitude data)
	* @param[in]  pos  the position of satellite in CRS
	* @param[in]  vel  the velocity of satellite in CRS
	* @return  calculate the rotation matrix correctly or not
	*/
	int t_gscf2crsleo::_calcNomMat(ColumnVector pos, ColumnVector vel)
	{
		try
		{
			if (pos.SumSquare() == 0.0 || vel.SumSquare() == 0.0)
			{
				cerr << "Wrong input: satellite positon/velocity" << endl;
				return -1;
			}
			int mode = 1;
			if (_sat_name == "GRACE-B"|| _sat_name == "GRACE-C"||_sat_name =="SENTINEL-3A"|| _sat_name == "SENTINEL-3B" )
			{
				mode = -1;
			}
			pos = pos * -1;
			vel = vel * mode;
			ColumnVector unitpos, unitvel, cross, unitcross;
			unitpos.resize(3);
			unitvel.resize(3);
			unitpos = pos / pos.NormFrobenius();
			unitvel = vel / vel.NormFrobenius();
			cross = crossproduct(unitpos, unitvel);
			unitcross = cross / cross.NormFrobenius();
			unitpos = crossproduct(unitvel, unitcross);
			_rotMatScf2crs.column(1) = unitvel;
			_rotMatScf2crs.column(2) = unitcross;
			_rotMatScf2crs.column(3) = unitpos;
			return 1;
		}
		catch (...)
		{
			cout << "ERROR:  t_gscf2crsleo::_calNomMat  unknown mistake" << endl;
			return -1;
		}
	}

	/**
	* @brief calculate the true rotation matrix from SCF to CRS(with attitude data)
	* @param[in]  time     the current epoch
	* @param[in]  attdata  the attitude data for the LEO satellite
	* @return   calculate the rotation matrix correctly or not
	*/
	int  t_gscf2crsleo::_calcAttMat(t_gtime time, t_gattitude* attdata)
	{
		try
		{
			if (attdata->beg() > time || attdata->end() < time)
			{
				/*cerr << "Input time beyond the time table of attitude file: "
					<< "Request time: " << dbl2str(time.dmjd()) << "   "
					<< "Begin time: " << dbl2str(attdata->beg().dmjd()) << "   "
					<< "End time: " << dbl2str(attdata->end().dmjd()) << endl;*/
				return -1;
			}
			if (attdata == NULL)
			{
				return -1;
			}
			t_attitude_record t_q;
			map<t_gtime, double> t_q1, t_q2, t_q3, t_q4;
			map<t_gtime, t_attitude_record>* attitude = attdata->attitude();
			map<t_gtime, t_attitude_record>::iterator beg, end;
			beg = attitude->lower_bound(time);
			end = attitude->upper_bound(time);
			if (beg == end)
			{
				beg--;
			}
			//lfound = false;
			t_q1[beg->first] = beg->second.q1;
			t_q1[end->first] = end->second.q1;
			t_q2[beg->first] = beg->second.q2;
			t_q2[end->first] = end->second.q2;
			t_q3[beg->first] = beg->second.q3;
			t_q3[end->first] = end->second.q3;
			t_q4[beg->first] = beg->second.q4;
			t_q4[end->first] = end->second.q4;

			double q1, q2, q3, q4;

			double alpha = (time.dmjd() - beg->first.dmjd()) / (end->first.dmjd() - beg->first.dmjd());
			q1 = t_q1[beg->first] + alpha*(t_q1[end->first] - t_q1[beg->first]);
			q2 = t_q2[beg->first] + alpha*(t_q2[end->first] - t_q2[beg->first]);
			q3 = t_q3[beg->first] + alpha*(t_q3[end->first] - t_q3[beg->first]);
			q4 = t_q4[beg->first] + alpha*(t_q4[end->first] - t_q4[beg->first]);

			if (abs(q1) < 1e-10&&abs(q2) < 1e-10&&abs(q3) < 1e-10&&abs(q4) < 1e-10)
			{
				return 0;
			}
			double s = q1*q1 + q2*q2 + q3*q3 + q4*q4;
			double q11 = q1*q1 / s;
			double q22 = q2*q2 / s;
			double q33 = q3*q3 / s;
			double q44 = q4*q4 / s;
			double q12 = q1*q2 / s;
			double q13 = q1*q3 / s;
			double q14 = q1*q4 / s;
			double q23 = q2*q3 / s;
			double q24 = q2*q4 / s;
			double q34 = q3*q4 / s;
			_rotMat(1, 1) = q11 + q44 - q22 - q33;
			_rotMat(2, 2) = q22 + q44 - q11 - q33;
			_rotMat(3, 3) = q33 + q44 - q11 - q22;
			_rotMat(1, 2) = 2.0*(q12 - q34);
			_rotMat(1, 3) = 2.0*(q13 + q24);
			_rotMat(2, 1) = 2.0*(q12 + q34);
			_rotMat(2, 3) = 2.0*(q23 - q14);
			_rotMat(3, 1) = 2.0*(q13 - q24);
			_rotMat(3, 2) = 2.0*(q23 + q14);
			return 1;
		}
		catch (...)
		{
			cout << "ERROR:  t_gscf2crsleo::_calAttMat  unknown mistake" << endl;
			return -1;
		}
	}

	t_gtime t_gscf2crsleo::getCurtEpoch()
	{
		return _epoch;
	}

	string t_gscf2crsleo::getSatName()
	{
		return _sat_name;
	}

	/**
	* @brief main function for rotation calculation
	* @param[in]  time    the current epoch
	* @param[in]  attdata the attitude data for the LEO satellite
	* @param[in]  pos     the position of satellite in CRS
	* @param[in]  vel     the velocity of satellite in CRS
	* @param[in]  trs2crs the rotation matrix from TRS to CRS
	* @param[out] scf2crs the rotaiton matrix from SCF to CRS
	* @return   calculate the rotation matrix correctly or not
	*/
	int t_gscf2crsleo::calcRotMat(t_gtime time, t_gattitude* attdata,const ColumnVector& pos,const ColumnVector& vel, Matrix& trs2crs, Matrix& scf2crs)
	{
		try
		{
			if (attdata != NULL)
			{
				_isAttValid = true;
			}
			
			if (_isAttValid && _calcAttMat(time, attdata) > 0)
			{
				if (_sat_name.substr(0, 4) == "COSM")
				{
					if (_calcNomMat(pos, vel) > 0)
					{
						_rotMat *= _rotMatScf2crs;
					}
				}
				if (_sat_name.substr(0, 4) == "SWAR")
				{
					_rotMat = trs2crs * _rotMat.t();
				}
				//cout << "Epoch: " << dbl2str(time.dmjd()) << ":  Attitude from att file" << endl;
			}
			else
			{
				if (_calcNomMat(pos, vel) > 0)
				{
					_rotMat = _rotMatScf2crs;
					//cout << "Epoch: " << dbl2str(time.dmjd()) << ":  Nominal attitude" << endl;
				}
				else
				{
					cerr << "Cannot determine the attitude for satellite " << _sat_name
						<< ", at epoch: " << dbl2str(time.dmjd()) << endl;
					scf2crs = _rotMat;
					return -1;
				}
			}
			scf2crs = _rotMat;
			return 1;
		}
		catch (...)
		{
			cout << "ERROR:  t_gscf2crsleo::cal_scf2crs  unknown mistake" << endl;
			return -1;
		}
	}

} // namespace
