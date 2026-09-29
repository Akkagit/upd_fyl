
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  This file is part of the G-Nut C++ library.
 
-*/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "gmodels/gfunclin.h"

using namespace std;

namespace gnut {   

// function (linH0)
// ---------
int t_gfunc_linH0::value( const ColumnVector& dat, const ColumnVector& fit, ColumnVector& val )
{
  if( fit.size() != 2 ){
    cout << "# warning: not properly initialized function (linH0)\n"; return -1;
  }
   
  for( int i=1; i<=dat.size(); ++i ){
   val(i) = fit(1) - fit(2)*dat(i);
  }
  return 0;
}
// ---------
int t_gfunc_linH0::deriv( const ColumnVector& dat, const ColumnVector& fit, Matrix& val )
{
  if( fit.size() != 2 ){
    cout << "# warning: not properly initialized function (linH0)\n"; return -1;
  }
  for( int i=1; i<=dat.size(); ++i ){
    val(i,1) =   1;
    val(i,2) = - dat(i);
  }
  return 0;
}


// function (linHS)
// ---------
int t_gfunc_linHS::value( const ColumnVector& dat, const ColumnVector& fit, ColumnVector& val )
{
  if( _coeff.find("H0") == _coeff.end() || fit.size() != 2 ){
    cout << "# warning: not properly initialized function (linHS)\n"; return -1;
  }

  double H0 = _coeff["H0"];

  for( int i=1; i<=dat.size(); ++i ){
   val(i) = fit(2) - fit(1)*(dat(i) - H0);
  }
  return 0;
}
// ---------
int t_gfunc_linHS::deriv( const ColumnVector& dat, const ColumnVector& fit, Matrix& val )
{
  if( _coeff.find("H0") == _coeff.end() || fit.size() != 2 ){
    cout << "# warning: not properly initialized function (linHS)\n"; return -1;
  }

  for( int i=1; i<=dat.size(); ++i ){
    val(i,1) = - dat(i);
    val(i,2) =   1;
  }
  return 0;
}


// function (linH1)
// ---------
int t_gfunc_linH1::value( const ColumnVector& dat, const ColumnVector& fit, ColumnVector& val )
{
  if( _coeff.find("a0") == _coeff.end() ||
      _coeff.find("H0") == _coeff.end() || fit.size() != 1 )
  {
    cout << "# warning: not properly initialized function (linH1)\n"; return -1;
  }

  double a0 = _coeff["a0"];
  double H0 = _coeff["H0"];
   
  for( int i=1; i<=dat.size(); ++i ){
    val(i) = a0 - fit(1)*(dat(i)-H0);
  }
  return 0;
}
// ---------
int t_gfunc_linH1::deriv( const ColumnVector& dat, const ColumnVector& fit, Matrix& val )
{
  if( _coeff.find("H0") == _coeff.end() || fit.size() != 1 )
  {
    cout << "# warning: not properly initialized function (linH1)\n"; return -1;
  }

  double H0 = _coeff["H0"];
      
  for( int i=1; i<=dat.size(); ++i ){
    val(i,1) = - (dat(i)-H0);
  }
  return 0;
}

// function (linT)
// ---------
int t_gfunc_linT::value( const ColumnVector& dat, const ColumnVector& fit, ColumnVector& val )
{
//  if( _coeff.find("a0") == _coeff.end() ||
//      _coeff.find("H0") == _coeff.end() || fit.size() != 1 )
//  {
//    cout << "# warning: not properly initialized function (linH1)\n"; return -1;
//  }

//  double a0 = _coeff["a0"];
//  double H0 = _coeff["H0"];
   for( int i=1; i<=dat.size(); ++i )
     {
	val(i) = fit(1)+fit(2)*dat(i);
     }
   return 0;
}
// ---------
int t_gfunc_linT::deriv( const ColumnVector& dat, const ColumnVector& fit, Matrix& val )
{
//  if( _coeff.find("H0") == _coeff.end() || fit.size() != 1 )
//{
//    cout << "# warning: not properly initialized function (linH1)\n"; return -1;
//  }

//  double H0 = _coeff["H0"];
   for( int i=1; i<=dat.size(); ++i )
     {
	val(i,1) = 1;
	val(i,2) = dat(i);
     }    
  return 0;
}  
   
// function (linCS)
// ---------
int t_gfunc_linCS::value( const ColumnVector& dat, const ColumnVector& fit, ColumnVector& val )
{
  if( fit.size() != 2 )
  {
    cout << "# warning: not properly initialized function (linCS)\n"; return -1;
  }

  for( int i=1; i<=dat.size(); ++i ){
    val(i) = fit(1)*cos(dat(i)) + fit(2)*sin(dat(i));
  }

  return 0;
}
// ---------
int t_gfunc_linCS::deriv( const ColumnVector& dat, const ColumnVector& fit, Matrix& val )
{
  if( fit.size() != 2 )
  {
    cout << "# warning: not properly initialized function (linCS)\n"; return -1;
  }

  for( int i=1; i<=dat.size(); ++i ){
    val(i,1) = cos(dat(i));
    val(i,2) = sin(dat(i));
  }

  return 0;
}


// function (linCSH)
// ---------
int t_gfunc_linCSH::value( const ColumnVector& dat, const ColumnVector& fit, ColumnVector& val )
{
  if( fit.size() != 3 )
  {
    cout << "# warning: not properly initialized function (linCSH)\n"; return -1;
  }

  for( int i=1; i<=dat.size(); ++i ){
    val(i) = ( fit(1)*cos(dat(i)) + fit(2)*sin(dat(i)) ) / fit(3);
  }

  return 0;
}
// ---------
int t_gfunc_linCSH::deriv( const ColumnVector& dat, const ColumnVector& fit, Matrix& val )
{
  if( fit.size() != 3 )
  {
    cout << "# warning: not properly initialized function (linCSH)\n"; return -1;
  }

  for( int i=1; i<=dat.size(); ++i ){
    val(i,1) = cos( dat(i) )/fit(3); // fit(3) = H_tropo
    val(i,2) = sin( dat(i) )/fit(3); // fit(3) = H_tropo
    val(i,3) = -(fit(1)*cos(dat(i)) + fit(2)*sin(dat(i))) /fit(3) /fit(3); // fit(3) = H_tropo
  }

  return 0;
}

} // namespace
