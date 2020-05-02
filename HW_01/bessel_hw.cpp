//filename: bessel_hw.cpp
//
//Author: Thomas DeMastri, email: demastri@msu.edu
//
//I edited the program made by Dick Furnstahl, email: furnstahl.1@osu.edu
//
//Revision History:
//
//Made: 2/7/2020:
//
//  Added a calculation of the relative difference and saved it to the file.
//
//  Analysis of my graph:
//  	The logscale plots are much more informative than my 
//  	In the region from x= 0 to 1, the upward and downard recursion heavily disagree.
//  	The error then falls off like a power law in the region from x=1 to 10, and 
//  	in the region between x=10 and x=40, the only difference between the two algorithms is from
//  	machine precision. 
//  	The error then steeply climbs after x=40, meaing it is very difficult to numerically approximate
//  	bessel functions for such large x.
//	One should use upward recursion for small  x, and really shouldn't use any approximations for large x.
//************************************************************************

// include files
#include <iostream>		// note that .h is omitted
#include <iomanip>		// note that .h is omitted
#include <fstream>		// note that .h is omitted

#include <gsl/gsl_sf.h>		//included for 3. c.)

#include <cmath>
using namespace std;		// we need this when .h is omitted

// function prototypes 
double down_recursion (double x, int n, int m);	// downward algorithm 
double up_recursion (double x, int n);	        // upward algorithm 

// global constants  
const double xmax = 100.0;	// max of x  
const double xmin = 0.1;	// min of x >0  
const double step = 0.1;	// delta x  
const int order = 10;		// order of Bessel function 
const int start = 50;		// used for downward algorithm 

//********************************************************************
int
main ()
{
  double ans_down, ans_up, relative_diff;

  // open an output file stream
  ofstream my_out ("bessel_hw.dat");

  my_out << "# Spherical Bessel functions via up and down recursion" 
         << endl;

  my_out << "        x             ans_up             ans_down             rel_diff             log10(x)          log10(rel_diff)            j10(x)" << endl;

  // step through different x values
  for (double x = xmin; x <= xmax; x += step)
    {
      ans_down = down_recursion (x, order, start);
      ans_up = up_recursion (x, order);

      relative_diff = abs(ans_up-ans_down)/(abs(ans_up)+abs(ans_down));
      
      my_out << fixed << setprecision (13) << setw (8) << x << " "
	<< scientific << setprecision (13)
	<< setw (13) << ans_down << " "
	<< setw (13) << ans_up << " "
	<< setw (13) << relative_diff << " "
	<< setw (13) << log10(x) << " "
	<< setw (13) << log10(relative_diff) << " "
	<< setw (13) << gsl_sf_bessel_jl(order, x) << " " //gsl routine added for part 3
        << endl;
    }
  cout << "data stored in bessel_hw.dat." << endl;

  // close the output file
  my_out.close ();
  return (0);
}


//------------------------end of main program----------------------- 

// function using downward recursion  
double
down_recursion (double x, int n, int m)
{
  double j[start + 2];		// array to store Bessel functions 
  j[m + 1] = j[m] = 1.;		// start with "something" (choose 1 here) 
  for (int k = m; k > 0; k--)
    {
      j[k - 1] = ((2.* double(k) + 1.) / x) * j[k] - j[k + 1];  // recur. rel.
    }
  double scale = (sin (x) / x) / j[0];	// scale the result 
  return (j[n] * scale);
}


//------------------------------------------------------------------ 

// function using upward recursion  
double
up_recursion (double x, int n)
{
  double term_three = 0.;
  double term_one = (sin (x)) / x;	// start with lowest order 
  double term_two = (sin (x) - x * cos (x)) / (x * x);	// next order
  for (int k = 1; k < n; k += 1)	// loop for order of function     
    { // recurrence relation
      term_three = ((2.*double(k) + 1.) / x) * term_two - term_one;	       
      term_one = term_two;
      term_two = term_three;
    }
  return (term_three);
}
