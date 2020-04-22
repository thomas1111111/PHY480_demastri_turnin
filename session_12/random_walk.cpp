//  file: random_walk.cpp
// 
//  Program to illustrate random walks
//
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//
//  Revision history:
//      03/05/04  translated from random_walk.c (11/25/02)
//      02/19/05  added more comments and math.h
//      02/14/06  added output comment
//
//  Notes:
//   * implements method 2 from the list in section 6.10
//      of the Landau/Paez text.
//   * random numbers are generated uniformly from a to b
//   * uses the GSL random number functions
//   * both the gsl_rng.h and gsl_randist.h header files are needed
//   * the current version uses the gsl_rng_taus random number
//      generator.  There are many other choices (just change
//      the name in the gsl_rng_alloc statement).  See the GSL
//      manual for a list of generators and their properties.
//
//******************************************************************

// include files
#include <iostream>		// cout and cin
#include <iomanip>		// manipulators like setprecision
#include <fstream>		// file input and output
#include <cmath>
using namespace std;		// we need this when .h is omitted

#include <gsl/gsl_rng.h>	// GSL random number generators
#include <gsl/gsl_randist.h>	// GSL random distributions

// function prototypes
extern unsigned long int random_seed ();   // routine to generate a seed

//********************************************************************
int
main (void)
{ 
  unsigned long int seed;	// seed for random number generator

  double lower = -sqrt (2.);	// lower limit of uniform range 
  double upper = sqrt (2.);	// upper limit of uniform range 

  double delta_x = 0.;		// uniform random number from [lower,upper] 
  double delta_y = 0.;		// 2nd random number from [lower,upper] 

  double x = 0.;
  double y = 0.;

  ofstream out;
  out.open ("random_walk_final.dat");

  out << "log10(N)     log10(R_avg)"  << endl;

  for (int npts = 100; npts<2e5; npts=npts*2)
  {
	
	double R_avg = 0.;

  	for (int N = 0; N < npts/2; N++) 
  	{
  		x = 0.;		// current x 
  		y = 0.;		// current y 

		double R = 0.;

  		gsl_rng *rng_ptr;		// declare pointer to random number 
                                //   generator (rng) 

  		rng_ptr = gsl_rng_alloc (gsl_rng_taus);	// allocate the rng 
 		seed = random_seed ();
  		gsl_rng_set (rng_ptr, seed); 

  		// do the walk and output to a file     
  		for (int i = 0; i < npts; i++)
    		{
      			delta_x = gsl_ran_flat (rng_ptr, lower, upper);
      			delta_y = gsl_ran_flat (rng_ptr, lower, upper);
      			x += delta_x;
      			y += delta_y;
    		} 

		R = sqrt(x*x + y*y);
		R_avg += R;
			
  		gsl_rng_free (rng_ptr);	// free the random number generator
  
  	}
        cout << npts << endl;
	R_avg = R_avg/(npts/2);
        out << log10(npts) << "     " << log10(R_avg) << endl;

  }

  out.close ();			// close the output file

  return (0);
}
