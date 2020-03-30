//  Name: Thomas DeMastri     email: demastri@msu.edu
//
//  filename: integ.cpp
//
//  This puprose of this program is to test the accuracy
//  of the integration routines we wrote.
// 
//  This program is based off of integ_test.cpp from session_03
//
//  Revision History:
//  Created: 2/27/20 - Added integrand and function of integral
//  		       Chose to use cos(x) as integrand since
//                     sin(x) is the analytical definite integral
//                     and there are routines to calculate this with
//                     high precision.
//  
//  Version 2:
//
//  Got my Simpson's Rule and Milne's rule to work. Adding the GSL Routine now.
//  Updated my integrand to A*cos(k*x+b) since that is a more interesting problem.
//
//  Version 3:
//
//  I had a lot of trouble getting the gsl integration routine to work,
//  since it requires a certain format for the integrand function to be in
//  that didn't play well with what I had already written.
//  
//  After many attempts to try to get one integrand function that works for
//  my rules and for the GSL library, I've decided it's a task for a more able coder.
//  I've written a new integrand function just for passing into the GSL routine.
//  I'm closely following an example published here: http://www.bnikolic.co.uk/nqm/1dinteg/gslgk.html
//  worked out by Bojan Nikolic, email: bojan@bnikolic.co.uk

//We need to include some packages:
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;

//We'll include the routines we wrote as well:
#include "integ_methods.h"

//And declare our functions:
double test_integrand(double x, double A, double k, double b);
double analytic_integral(double x, double A, double k, double b);

int
main()
{
  //We'll start by putting some bounds on our problem:
  double xmin = 0.;
  double xmax = M_PI*2.;

  //And set up the parameters for our integrand:
  double A = 1.5; //amplitude
  double k = 2.0; //wave number
  double b = M_PI/6.; //phase 

  double exact_answer = analytic_integral(xmax,A,k,b)-analytic_integral(xmin,A,k,b); //calculating the analytical answer for error purposes

  //Open up a data file to write to:
  ofstream my_out ("integ.dat");
 
  my_out << setw(16) << "N" << setw(16) << "Simpson's" << setw(16) << "Milne's" << setw(16) << "GSL" << setw(16) << "Exact" << endl;
  cout << exact_answer << endl;
  //

}

	
double 
test_integrand(double x, double A, double k, double b)
{
   //declaring our function to integrate over	
   
   return(A*cos(k*x+b));


}

double 
analytic_integral(double x, double A, double k, double b)
{   
   //this is the exact analytical answer to our test integrand,
   //included here for error analysis

   return((A/k)*sin(k*x+b));


}


//GETTING THE GSL ROUTINE TO WORK:

//we're using the integ_params structure defined in the header file here, and defining the integrand
//declared in the header as well.

double
integrand_gsl(double x, void *p)
{
//this is the format that gsl integration likes its integrands in
//so we're stuck with it as far as I can tell
   integ_params &params = *reinterpret_cast<integ_params *>(p); //this is recasting the void parameters we passed in as parameters in the structure
   								//described above
   return (params.A*cos(params.k*x+params.b)); //calculating our integrand

}











	
