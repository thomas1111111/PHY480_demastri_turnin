//  file: gsl_spline_test_class.cpp
// 
//  Test program for the gsl spline routines using the Spline class
//
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//
//  Revision history:
//      02/10/09 -- created from gsl_cubic_spline_test.cpp
//
//  Notes:
//   * uses the GSL interpolation functions (see online documentation) 
//
//*****************************************************************
// include files
#include <iostream>    // cout and cin
#include <iomanip>     // manipulators like setprecision
#include <fstream>
#include <cmath>
#include <string>     // C++ strings                                 
using namespace std;    
#include "GslSpline.h"  // Header file for the GSL Spline class

inline double sqr (double z) {return z*z;}  // inline function for z^2
double sigma_th(double sigma0, double E, double Er, double gamma);

int
main (void)
{
  const int NMAX = 300;   // maximum number of array points 
  double x_values[NMAX], y_values[NMAX];
	
  int npts = 9;
  double sigma0 = 63900;
  double Er = 78.0;
  double gamma = 55.0;
  double deltax = 25.0;


  double y_start[npts] = {10.6,16.0,45.0,85.5,52.8,19.9,10.8,8.25,4.7};

  for(int i = 0; i < npts; i++)
  {
	double curr_x = deltax*i;
	x_values[i] = curr_x;
	y_values[i] = y_start[i];
  }

  // Make the spline objects
  Spline my_cubic_spline (x_values, y_values, npts, "cubic");
  Spline my_linear_spline (x_values, y_values, npts, "linear");
  Spline my_poly_spline (x_values, y_values, npts, "polynomial");

  
  // Evaluate the spline and derivatives
  
  ofstream my_out;
  my_out.open("spline.dat");

  my_out << "    x     sig_the   sig_cubi   sig_line   sig_poly   " << endl;

  for(double x = 0.; x < 201.; x +=5.)
  {
  double y_cub = my_cubic_spline.y (x);
  double y_lin = my_linear_spline.y (x);
  double y_pol = my_poly_spline.y (x);

  double sig_the = sigma_th(sigma0, x, Er, gamma);
  
  my_out << fixed << setprecision(6) 
       << x << "   " << sig_the << "  "<< y_cub << "  " << y_lin << "  "<< y_pol << endl;
  
  }
  
  my_out.close();

  return (0);      // successful completion 
}

double
sigma_th(double sigma0, double E, double Er, double gamma)
{
  return sigma0/(sqr(E-Er)+sqr(gamma/2));
}

